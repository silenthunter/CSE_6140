#include <cuda.h>
#include <stdio.h>
#include <iostream>
#include <sys/time.h>

using namespace std;

const int BLOCK_WIDTH = 2;
const int BLOCK_HEIGHT = 2;
const int DEFAULT_ELE = 16;
extern __shared__ int shmem[];

typedef struct __align__(8) linkNode
{
	int edge;
	int next;
} linkNode;

//Create CSR edge storage
__global__ void convertEdges(int* edges, int numEdge, int numVert, int* newArrays)
{
	int* val = newArrays;
	int* col = &val[numEdge];// + numEdge;
	int* row = &col[numEdge];

	int lastRow = -1;
	int rowNum = 0;
	for(int i = 0; i < numEdge; i++)
	{
		int v1 = edges[i * 2];//Row
		int v2 = edges[i * 2 + 1];//Col

		val[i] = 1;//Weight
		col[i] = v2;

		if(lastRow != v1)
		{
			//Make sure to skip unconnected vertices
			int diff = v1 - lastRow;
			for(int j = 0; j < diff - 1; j++)
				row[rowNum++] = i;
			row[rowNum++] = i;
			lastRow = v1;
		}
	}

	//Fill in final row columns
	for(int i = rowNum; i <= numVert; i++)
		row[i] = numEdge;
}

__device__ int findNext(int* edges, int numEdge, int numVert, int v, int* destination)
{
	int* val = edges;
	int* col = &val[numEdge];
	int* row = &col[numEdge];

	int idx = row[v];
	int nextIdx = row[v + 1];

	int count = 0;
	for(int i = idx; i < nextIdx; i++)
	{
		destination[count++] = col[i];
	}

	return count;
}

__device__ void pushQueue(int element, int* queue, int queueSize, unsigned int* head, unsigned int* tail)
{
	int idx = atomicInc(tail, queueSize);
	queue[idx] = element;
}

__device__ int popQueue(int* queue, int queueSize, unsigned int* head, unsigned int* tail)
{
	int retn;
	if(*head == *tail) retn = -1;
	else
	{
		int idx = atomicInc(head, queueSize);
		retn = queue[idx];
	}
	
	return retn;
}

__device__ void pushStack(int element, int* stack, int* head)
{
	int idx = atomicAdd(head, 1);
	stack[idx] = element;
}

__device__ int popStack(int* stack, int* head)
{
	int retn;
	if(*head == 0) retn = -1;
	else
	{
		int idx = atomicSub(head, 1);
		retn = stack[idx - 1];
	}
	
	return retn;
}

__device__ int ELEMENTS = DEFAULT_ELE;
__device__ int S_SIZE;
__device__ int P_SIZE;
__device__ int PATH_SIZE;
__device__ int D_SIZE;
__device__ int Q_SIZE;

__device__ void doAlg(int numVert, int* edges, int numEdges, linkNode* pList, float* BC, int* glob, float* globDep)
{
	S_SIZE = numVert;
	P_SIZE = numEdges + numVert;
	D_SIZE = numVert;
	Q_SIZE = numVert;
	PATH_SIZE = numVert;
	ELEMENTS = numVert;

	int block_idx = gridDim.x * blockIdx.y + blockIdx.x;
	int localIdx = threadIdx.x + threadIdx.y * blockDim.x;

	unsigned int PTR_OFFSET = block_idx * (S_SIZE + D_SIZE + Q_SIZE + PATH_SIZE);
	
	int* S = &glob[PTR_OFFSET];
	int* S_head = &shmem[4];
	*S_head = 0;
	PTR_OFFSET += S_SIZE;
	
	linkNode* P = &pList[block_idx * P_SIZE];
	//Blank the previous items
	for(int i = 0; i < P_SIZE && localIdx == 0; i++)
	{
		P[i].edge = -1;
		P[i].next = -1;
	}

	int* pathCount = &glob[PTR_OFFSET];
	for(int i = 0; i < PATH_SIZE; i++)
	{
		pathCount[i] = 0;
	}
	pathCount[block_idx] = 1;
	PTR_OFFSET += PATH_SIZE;

	int* d = &glob[PTR_OFFSET];
	for(int i = 0; i < D_SIZE; i++)
	{
		d[i] = -1;
	}
	d[block_idx] = 0;
	PTR_OFFSET += D_SIZE;
	
	int* Q = &glob[PTR_OFFSET];
	unsigned int* Q_head = (unsigned int*)&shmem[2];
	unsigned int* Q_tail = (unsigned int*)&shmem[3];
	if(localIdx == 0)
	{
		*Q_head = 0;
		*Q_tail = 0;
	}
	PTR_OFFSET += Q_SIZE;
	
	if(localIdx == 0)
		pushQueue(block_idx, Q, Q_SIZE, Q_head, Q_tail);

	int* front = &shmem[0];
	int* nextFront = &shmem[1];
	if(localIdx == 0)
	{
		*front = 1;
		*nextFront = 0;
	}

	int* testCnt = &shmem[5];
	*testCnt = 0;

	while(*Q_head != *Q_tail || *front != 0 || *nextFront != 0)
	{
		__syncthreads();
		if(*front == 0 && localIdx == 0)
		{
			*front = *nextFront;
			*nextFront = 0;
		}
		__syncthreads();
		int v = -1;
		int atmSub = atomicSub(front, 1);
		if(atmSub > 0)
			v = popQueue(Q, Q_SIZE, Q_head, Q_tail);
		else
			atomicAdd(front, 1);//Don't let the counter go below 0
		if(v < 0) continue;

		pushStack(v, S, S_head);

		int w[1024];
		int edgeCount = findNext(edges, numEdges, numVert, v, w);

		for(int i = 0; i < edgeCount; i++)
		{
			int wNode = w[i];
			if(d[wNode] < 0 && atomicCAS(&d[wNode], -1, d[v] + 1) == -1)
			{
				pushQueue(wNode, Q, Q_SIZE, Q_head, Q_tail);
				atomicAdd(nextFront, 1);//Frontier expands
			}

			if(d[wNode] == d[v] + 1)
			{
				atomicAdd(&pathCount[wNode], pathCount[v]);
				//pathCount[wNode] = pathCount[wNode] + pathCount[v];
				
				//Append v to the PrevNode list

				//Find the next empty slot. Start after the initial lookup indices
				if(atomicCAS(&P[wNode].edge, -1, v) < 0);
				else
				{
					atomicAdd(testCnt, 1);
					int empty;
					int casToken = 0;
					for(empty = numVert; casToken != -1; empty++)
					{
						casToken = atomicCAS(&P[empty].edge, -1, v);
					}
					empty--;

					linkNode* j = &P[wNode];
					//while(j !=  NULL)
					while(atomicCAS(&(j->next), -1, empty) >= 0)
					{
						j = &P[j->next];
					}
				}
				
				
			}
		}

	}
	
	float* dep = &globDep[block_idx * numVert];

	__syncthreads();
	
	while(*S_head > 0 && localIdx == 0)
	{
		int w = popStack(S, S_head);
		if(w == block_idx) continue;
		
		//Loop through each v in P[w]
		linkNode* node = &P[w];
		while(node != NULL)
		{
			int v = node->edge;
			if(v < 0) continue;

			dep[v] = dep[v] + ((float)pathCount[v]/(float)pathCount[w]) * (1 + dep[w]);

			node = node->next < 0 ? NULL : &P[node->next];
		}
		
		if(w != block_idx)
		{
			atomicAdd(&BC[w], dep[w]);
		}
	}
	
}

__global__ void betweennessCentrality(int numVert, int numEdges, int *edges, linkNode* pList, float* BC, int* glob, float* dep)
{
	//sortEdges(edges, path);
	int block_idx = blockIdx.x + blockIdx.y * gridDim.x;

	if(block_idx >= numVert) return;

	shmem[0] = 0;

	__syncthreads();

	doAlg(numVert, edges, numEdges, pList, BC, glob, dep);

		
}

int edgeCompare(const void* a, const void* b)
{
	int* av1 = (int*)a;
	int* av2 = av1 + 1;
	int* bv1 = (int*)b;
	int* bv2 = bv1 + 1;

	if(*av1 < *bv1) return -1;
	else if(*av1 > *bv1) return 1;
	else if(*av2 < *bv2) return -1;
	else if(*av2 > *bv2) return 1;
	else return 0;
}

int main(int argc, char* argv[])
{
	int elements = DEFAULT_ELE;

	//cudaProfilerStart();
	int *d_mem;
	int *h_edge;
	int *d_edge;
	int *d_optEdge;
	linkNode* pList;
	float *d_bc;
	float *h_bc;
	int *d_glob;
	float *d_dep;

	int numVert = elements;
	int numEdge = elements - 1;

	if(argc < 2)
	{

		FILE *grFile;
		grFile = fopen("test.gr", "w");
		fprintf(grFile, "p %d %d d u 0\n", numVert, numEdge);
		h_edge = (int*)malloc(sizeof(int) * numEdge * 2);
		
		//Init edges
		for(int i = 0; i < numEdge; i++)
		{
			h_edge[i * 2] = i % numVert;
			h_edge[i * 2 + 1] = (i + 1) % numVert;
			fprintf(grFile, "%d %d\n", h_edge[i * 2], h_edge[i * 2 + 1]);
		}
		fclose(grFile);
	}
	else
	{
		FILE *grFile;
		char buff[1024];
		int edgeCnt = 0;
		grFile = fopen(argv[1], "r");
		while(!feof(grFile))
		{
			fgets(buff, 1024, grFile);
			if(feof(grFile))break;

			//This is the  "problem" line
			if(buff[0] == 'p')
			{
				char* token = strtok(buff, " ");

				token = strtok(NULL, " ");
				token = strtok(NULL, " ");
				numVert = atoi(token);
				
				token = strtok(NULL, " ");
				numEdge = atoi(token) * 2;

				h_edge = (int*)malloc(sizeof(int) * numEdge * 2);

			}
			else if(buff[0] == '#' || buff[0] == 'c')
				continue;
			else if(buff[0] == 'a')
			{
				char* token = strtok(buff, " ");

				//Skip 'a'
				token = strtok(NULL, " ");
				int e1 = atoi(token) - 1;
				h_edge[edgeCnt] = e1;
				h_edge[edgeCnt + 3] = e1;
				token = strtok(NULL, " ");
				int e2 = atoi(token) - 1;
				h_edge[edgeCnt + 1] = e2;
				h_edge[edgeCnt + 2] = e2;
				edgeCnt += 4;
			}
		}
	}

	//Sort the arrays for CSR
	qsort(h_edge, numEdge, sizeof(int) * 2, edgeCompare);
	
	long totalMem = 0;
	cudaMalloc((void**)&d_optEdge, sizeof(int) * (numVert + 1) * (numEdge * 2));
	totalMem += sizeof(int) * (numVert + 1) * (numEdge * 2);


	struct timeval totalStart, totalEnd;
	gettimeofday(&totalStart, NULL);

	//First convert edges
	cudaMalloc((void**)&d_edge, sizeof(int) * numEdge * 2);
	//totalMem += sizeof(int) * numEdge * 2;
	cudaMemcpy(d_edge, h_edge, sizeof(int) * numEdge * 2, cudaMemcpyHostToDevice);
	convertEdges<<<1,1>>>(d_edge, numEdge, numVert, d_optEdge);
	cudaDeviceSynchronize();
	cudaFree(d_edge);

	//Test code
	int* test = (int*)malloc(sizeof(int) * (numVert + 1) * (numEdge * 2));
	cudaMemcpy(test, d_optEdge, sizeof(int) * (numVert + 1) * (numEdge * 2), cudaMemcpyDeviceToHost);
	int* val = test;
	int* col = test + numEdge;
	int* row = col + numEdge;
	

	h_bc = (float*)malloc(sizeof(float) * numVert);
	cudaMalloc((void**)&d_mem, sizeof(int) * numVert);
	totalMem += sizeof(int) * numVert;
	
	cudaMalloc((void**)&d_bc, sizeof(float) * numVert);
	cudaMemset(d_bc, 0, sizeof(float) * numVert);
	totalMem += sizeof(float) * numVert;

	cudaMalloc((void**)&d_glob, sizeof(int) * numVert * (numVert * 5));
	totalMem += sizeof(int) * numVert * ((numVert * 5));

	cudaMalloc((void**)&pList, sizeof(linkNode) * numVert * (numVert + numEdge));
	totalMem += sizeof(linkNode) * numVert * (numVert + numEdge);

	cudaMalloc((void**)&d_dep, sizeof(float) * numVert * numVert);
	totalMem += sizeof(float) * numVert * numVert;
	
	dim3 block(BLOCK_WIDTH, BLOCK_HEIGHT);
	int gridSize = numVert;//ceil(numVert / (float)(BLOCK_WIDTH * BLOCK_HEIGHT));
	dim3 grid(gridSize);

	struct timeval start, end;

	gettimeofday(&start, NULL);
	betweennessCentrality<<<grid,block, 512>>>(numVert, numEdge, d_optEdge, pList, d_bc, d_glob, d_dep);
	cudaDeviceSynchronize();
	gettimeofday(&end, NULL);
	long elapsed = (end.tv_sec * 1000000 + end.tv_usec) - (start.tv_sec * 1000000 + start.tv_usec);

	cudaError_t error = cudaGetLastError();
	
	int* h_mem = (int*)malloc(sizeof(int) * numVert);
	cudaMemcpy(h_mem, d_mem, sizeof(int) * numVert, cudaMemcpyDeviceToHost);
	cudaMemcpy(h_bc, d_bc, sizeof(float) * numVert, cudaMemcpyDeviceToHost);

	gettimeofday(&totalEnd, NULL);
	long totalElapsed = (totalEnd.tv_sec * 1000000 + totalEnd.tv_usec) - 
		( totalStart.tv_sec * 1000000 + totalStart.tv_usec); 

	for(int i = 0; i < numVert; i++)
	{
		cout << h_bc[i] << endl;
	}
	//cout<<elements<<endl;
	
	//cudaProfilerStop();
	
	cudaDeviceReset();
	cout << cudaGetErrorString(error) << endl;
	cout << "Mem Used: " << totalMem << endl;
	cout << "Time(usec): " << elapsed << endl;
	cout << "Total Time(usec): " << totalElapsed << endl;
	
	return 0;
}
