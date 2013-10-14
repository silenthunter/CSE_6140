#include <cuda.h>
#include <stdio.h>
#include <iostream>

using namespace std;

__device__ const int MAX_DEGREE = 4;

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

const int ELEMENTS = 64;
const int S_SIZE = ELEMENTS;
const int P_SIZE = ELEMENTS;
const int P_LIST_SIZE = ELEMENTS;
const int PATH_SIZE = ELEMENTS;
const int D_SIZE = ELEMENTS;
const int Q_SIZE = ELEMENTS;


__device__ void doAlg(int numVert, int* edges, int numEdges, int* BC, int* glob)
{
	int x = blockDim.x * blockIdx.x + threadIdx.x;	
	int y = blockDim.y * blockIdx.x + threadIdx.y;	
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
			P[i * P_SIZE + j] = -1;
		}
	PTR_OFFSET += P_SIZE * MAX_DEGREE;

	int* pathCount = &glob[PTR_OFFSET];
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

		int w[8];
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
					if(P[wNode * P_SIZE + j] < 0)
					{
						P[wNode * P_SIZE + j] = v;
						break;
					}
				}
			}
		}

	}
	
	float* dep = (float*)&glob[PTR_OFFSET];
	
	while(S_head > 0)
	{
		int w = popStack(S, &S_head);
		
		//Loop through each v in P[w]
		for(int i = 0; i < MAX_DEGREE; i++)
		{
			int v = P[0];//P[w * P_SIZE + i];
			//dep[v] = dep[v] + (pathCount[v]/pathCount[w]) * (1 + dep[w]);
		}
		
		if(w != idx)
		{
			//TODO: Atomic
			BC[w] = 1;//BC[w] + dep[w];
		}
	}
	
}

__global__ void betweennessCentrality(int numVert, int numEdges, int *edges, int* BC, int* glob)
{
	extern __shared__ int path[];
	
	//sortEdges(edges, path);
	int x = blockDim.x * blockIdx.x + threadIdx.x;	
	int y = blockDim.y * blockIdx.x + threadIdx.y;	
	int idx = x + y * blockDim.x * gridDim.x;

	doAlg(numVert, edges, numEdges, BC, glob);

		
}

int main()
{
	const int elements = 64;

	//cudaProfilerStart();
	int *d_mem;
	int *h_edge;
	int *d_edge;
	int *d_bc;
	int *h_bc;
	int *d_glob;
	
	cudaMalloc((void**)&d_mem, sizeof(int) * elements);
	
	h_edge = (int*)malloc(sizeof(int) * elements * 2);
	cudaMalloc((void**)&d_edge, sizeof(int) * elements * 2);

	h_bc = (int*)malloc(sizeof(int) * elements);
	cudaMalloc((void**)&d_bc, sizeof(int) * elements);

	cudaMalloc((void**)&d_glob, sizeof(int) * elements * elements * 20);
	
	//Init edges
	for(int i = 0; i < elements; i++)
	{
		h_edge[i * 2] = i % elements;
		h_edge[i * 2 + 1] = (i + 1) % elements;
	}
	cudaMemcpy(d_edge, h_edge, sizeof(int) * elements * 2, cudaMemcpyHostToDevice);
	
	dim3 block(8,8);
	dim3 grid(elements / 64);
	//test<<<grid,block>>>(d_mem);
	betweennessCentrality<<<grid,block>>>(elements, elements, d_edge, d_bc, d_glob);
	cudaError_t error = cudaGetLastError();
	
	int* h_mem = (int*)malloc(sizeof(int) * elements);
	cudaMemcpy(h_mem, d_mem, sizeof(int) * elements, cudaMemcpyDeviceToHost);
	cudaMemcpy(h_bc, d_bc, sizeof(int) * elements, cudaMemcpyDeviceToHost);
	

	for(int i = 0; i < elements; i++)
	{
		cout << h_bc[i] << endl;
	}
	//cout<<elements<<endl;
	
	//cudaProfilerStop();
	
	cudaDeviceReset();
	cout << cudaGetErrorString(error) << endl;
	
	return 0;
}
