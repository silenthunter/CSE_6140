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
	*tail = (*tail + 1) % queueSize;
	queue[tail] = element;
}

__device__ int popQueue(int* queue, int queueSize, int* head, int* tail)
{
	int retn = queue[*head];
	*head = (*head + 1) % queueSize;
	
	return retn;
}

__device__ void pushStack(int element, int* stack, int* head)
{
	*head = *head + 1;
	stack[tail] = element;
}

__device__ int popStack(int* stack, int* head)
{
	int retn = stack[*head];
	*head = *head - 1;
	
	return retn;
}


__device__ void doAlg(int numVert)
{
	int x = blockDim.x * blockIdx.x + threadIdx.x;	
	int y = blockDim.y * blockIdx.x + threadIdx.y;	
	int idx = x + y * blockDim.x * gridDim.x;
	
	int S[];
	int S_head = 0;
	
	int P[];
	pathCount[]; pathCount[idx] = 1;
	int d[]; d[idx] = 0;
	
	int Q_size;
	int Q[Q_size];
	int Q_head = 0;
	int Q_tail = 0;
	
	pushQueue(idx, Q, Q_size, &Q_head, &q_tail);
	
	
}

__global__ void betweennessCentrality(int numVert, int numEdges, int *edges, int* BC)
{
	extern __shared__ int path[];
	
	//sortEdges(edges, path);
	int x = blockDim.x * blockIdx.x + threadIdx.x;	
	int y = blockDim.y * blockIdx.x + threadIdx.y;	
	int idx = x + y * blockDim.x * gridDim.x;

	int arr[8];
	int count = findNext(edges, numEdges, idx, arr);
	if(count > 0)
		BC[idx] = arr[0];
	else
		BC[idx] = -1;
		
}

int main()
{
	const int elements = 1024;

	//cudaProfilerStart();
	int *d_mem;
	int *h_edge;
	int *d_edge;
	int *d_bc;
	int *h_bc;
	
	cudaMalloc((void**)&d_mem, sizeof(int) * elements);
	
	h_edge = (int*)malloc(sizeof(int) * elements * 2);
	cudaMalloc((void**)&d_edge, sizeof(int) * elements * 2);

	h_bc = (int*)malloc(sizeof(int) * elements);
	cudaMalloc((void**)&d_bc, sizeof(int) * elements);
	
	//Init edges
	for(int i = 0; i < elements; i++)
	{
		h_edge[i * 2] = i % elements;
		h_edge[i * 2 + 1] = (i + 1) % elements;
	}
	cudaMemcpy(d_edge, h_edge, sizeof(int) * elements * 2, cudaMemcpyHostToDevice);
	
	dim3 block(32,32);
	dim3 grid(elements / 1024);
	//test<<<grid,block>>>(d_mem);
	betweennessCentrality<<<grid,block, sizeof(int) * elements * MAX_DEGREE>>>(elements, elements, d_edge, d_bc);
	cudaError_t error = cudaGetLastError();
	cout << cudaGetErrorString(error) << endl;
	
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
	
	return 0;
}
