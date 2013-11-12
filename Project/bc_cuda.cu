#include <cuda.h>
#include <stdio.h>
#include <iostream>

using namespace std;

__device__ const int MAX_DEGREE = 50;

const int BLOCK_WIDTH = 2;
const int BLOCK_HEIGHT = 2;
const int DEFAULT_ELE = 1024;

__device__ void sortEdges(int* edges, int* sorted)
{
	int x = blockDim.x * blockIdx.x + threadIdx.x;
	int y = blockDim.y * blockIdx.y + threadIdx.y;
	int i = x + y * gridDim.x * blockDim.x;
	int n1 = edges[i * 2];
	int n2 = edges[i * 2 + 1];
	
	int* arrStart = &sorted[n1 * MAX_DEGREE];
	int retnVal = 1;
	while(retnVal != 0)
		retnVal = atomicCAS(arrStart, 0, n2);
}


//HACK: This will be incredibly slow  on CUDA!
__device__ int findNext(int* edges, int numEdge, int v, int* destination)
{
	int count = 0;

	for(int i = 0; i < numEdge * 2; i+=2)
	{
		if(edges[i] == v)
			destination[count++] = edges[i + 1];
	}

	return count;
}

__device__ void pushQueue(int element, int* queue, int queueSize, int* head, int* tail)
{
	queue[*tail] = element;
	*tail = (*tail + 1) % queueSize;
}

__device__ int popQueue(int* queue, int queueSize, int* head, int* tail)
{
	int retn = queue[*head];
	*head = (*head + 1) % queueSize;
	
	return retn;
}

__device__ void pushStack(int element, int* stack, int* head)
{
	stack[*head] = element;
	*head = *head + 1;
}

__device__ int popStack(int* stack, int* head)
{
	*head = *head - 1;
	int retn = stack[*head];
	
	return retn;
}

__device__ int ELEMENTS = DEFAULT_ELE;
__device__ int S_SIZE;
__device__ int P_SIZE;
__device__ int PATH_SIZE;
__device__ int D_SIZE;
__device__ int Q_SIZE;

__device__ void doAlg(int numVert, int* edges, int numEdges, float* BC, int* glob, float* globDep)
{
	S_SIZE = numVert;
	P_SIZE = numVert;
	D_SIZE = numVert;
	Q_SIZE = numVert;
	PATH_SIZE = numVert;
	ELEMENTS = numVert;

	int x = blockDim.x * blockIdx.x + threadIdx.x;	
	int y = blockDim.y * blockIdx.y + threadIdx.y;	
	int idx = x + y * blockDim.x * gridDim.x;

	int PTR_OFFSET = idx * (S_SIZE + (P_SIZE * MAX_DEGREE) + D_SIZE + Q_SIZE + PATH_SIZE);
	
	int* S = &glob[PTR_OFFSET];
	int S_head = 0;
	PTR_OFFSET += S_SIZE;
	
	int* P = &glob[PTR_OFFSET];
	//Blank the previous items
	for(int i = 0; i < P_SIZE; i++)
		for(int j = 0; j < MAX_DEGREE; j++)
		{
			P[i + P_SIZE * j] = -1;
		}
	PTR_OFFSET += P_SIZE * MAX_DEGREE;

	int* pathCount = &glob[PTR_OFFSET];
	for(int i = 0; i < PATH_SIZE; i++)
	{
		pathCount[i] = 0;
	}
	pathCount[idx] = 1;
	PTR_OFFSET += PATH_SIZE;

	int* d = &glob[PTR_OFFSET];
	for(int i = 0; i < D_SIZE; i++)
	{
		d[i] = -1;
	}
	d[idx] = 0;
	PTR_OFFSET += D_SIZE;
	
	int* Q = &glob[PTR_OFFSET];
	int Q_head = 0;
	int Q_tail = 0;
	PTR_OFFSET += Q_SIZE;
	
	pushQueue(idx, Q, Q_SIZE, &Q_head, &Q_tail);

	while(Q_head != Q_tail)
	{
		int v = popQueue(Q, Q_SIZE, &Q_head, &Q_tail);
		pushStack(v, S, &S_head);

		int w[MAX_DEGREE];
		int edgeCount = findNext(edges, numEdges, v, w);

		for(int i = 0; i < edgeCount; i++)
		{
			int wNode = w[i];
			if(d[wNode] < 0)
			{
				pushQueue(wNode, Q, Q_SIZE, &Q_head, &Q_tail);
				d[wNode] = d[v] + 1;
			}
			
			if(d[wNode] == d[v] + 1)
			{
				pathCount[wNode] = pathCount[wNode] + pathCount[v];
				
				//Append v to the PrevNode list
				for(int j = 0; j < MAX_DEGREE; j++)
				{
					if(P[wNode + P_SIZE * j] < 0)
					{
						P[wNode + P_SIZE * j] = v;
						break;
					}
				}
			}
		}

	}
	
	float* dep = &globDep[idx * ELEMENTS];
	
	while(S_head > 0)
	{
		int w = popStack(S, &S_head);
		
		//Loop through each v in P[w]
		for(int i = 0; i < MAX_DEGREE; i++)
		{
			int v = P[w + P_SIZE * i];
			if(v < 0) continue;

			dep[v] = dep[v] + ((float)pathCount[v]/(float)pathCount[w]) * (1 + dep[w]);
		}
		
		if(w != idx)
		{
			atomicAdd(&BC[w], dep[w]);
		}
	}
	
}

__global__ void betweennessCentrality(int numVert, int numEdges, int *edges, float* BC, int* glob, float* dep)
{
	extern __shared__ int path[];
	
	//sortEdges(edges, path);
	int x = blockDim.x * blockIdx.x + threadIdx.x;	
	int y = blockDim.y * blockIdx.x + threadIdx.y;	
	int idx = x + y * blockDim.x * gridDim.x;

	if(idx >= numVert) return;

	BC[idx] = 0.0f;

	__syncthreads();

	doAlg(numVert, edges, numEdges, BC, glob, dep);

		
}

int main(int argc, char* argv[])
{
	int elements = DEFAULT_ELE;

	//cudaProfilerStart();
	int *d_mem;
	int *h_edge;
	int *d_edge;
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
	h_bc = (float*)malloc(sizeof(float) * numVert);
	cudaMalloc((void**)&d_mem, sizeof(int) * numVert);
	
	cudaMalloc((void**)&d_edge, sizeof(int) * numEdge * 2);

	cudaMalloc((void**)&d_bc, sizeof(float) * numVert);

	cudaMalloc((void**)&d_glob, sizeof(int) * numVert * numVert * 8);
	cudaMalloc((void**)&d_dep, sizeof(float) * numVert * numVert);
	cudaMemcpy(d_edge, h_edge, sizeof(int) * numEdge * 2, cudaMemcpyHostToDevice);
	
	dim3 block(BLOCK_WIDTH, BLOCK_HEIGHT);
	int gridSize = ceil(numVert / (float)(BLOCK_WIDTH * BLOCK_HEIGHT));
	dim3 grid(gridSize);
	//test<<<grid,block>>>(d_mem);
	betweennessCentrality<<<grid,block>>>(numVert, numEdge, d_edge, d_bc, d_glob, d_dep);
	cudaError_t error = cudaGetLastError();
	
	int* h_mem = (int*)malloc(sizeof(int) * numVert);
	cudaMemcpy(h_mem, d_mem, sizeof(int) * numVert, cudaMemcpyDeviceToHost);
	cudaMemcpy(h_bc, d_bc, sizeof(float) * numVert, cudaMemcpyDeviceToHost);
	

	for(int i = 0; i < numVert; i++)
	{
		cout << h_bc[i] << endl;
	}
	//cout<<elements<<endl;
	
	//cudaProfilerStop();
	
	cudaDeviceReset();
	cout << cudaGetErrorString(error) << endl;
	
	return 0;
}
