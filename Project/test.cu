#include <cuda.h>
#include <stdio.h>
#include <iostream>
#include <thrust/device_vector.h>

using namespace std;

__device__ const int MAX_DEGREE = 4;

__global__ void test(int* d_mem)
{
	int x = blockIdx.x * blockDim.x + threadIdx.x;
	int y = blockIdx.y * blockDim.y + threadIdx.y;
	int elementPitch = blockDim.x * gridDim.x;
	int i = y * elementPitch + x;
	for(int j = 0; j < 1024; j++)
		d_mem[i] = warpSize;
}

__device__ void sortEdges(int* edges, int* sorted)
{
	int x = blockDim.x * blockIdx.x + threadIdx.x;
	int y = blockDim.y * blockIdx.y + threadIdx.y;
	int i = x + y * gridDim.x * blockDim.x;
	int n1 = edges[i * 2];
	int n2 = edges[i * 2 + 1];
	
	int* arrStart = &sorted[n1 * MAX_DEGREE];
	int retnVal = 1;
	//while(retnVal != 0)
		retnVal = atomicCAS(arrStart, 0, n2);
}

__global__ void betweennessCentrality(int numVert, int numEdges, int *edges, int* BC)
{
	extern __shared__ int path[];
	
	sortEdges(edges, path);
	//int y = blockDim.x * threadIdx.y;
	//int x = threadIdx.x;
}

int main()
{
	const int elements = 1024;

	//cudaProfilerStart();
	int *d_mem;
	int *h_edge;
	int *d_edge;
	
	cudaMalloc((void**)&d_mem, sizeof(int) * elements);
	
	h_edge = (int*)malloc(sizeof(int) * elements * 2);
	cudaMalloc((void**)&d_edge, sizeof(int) * elements * 2);
	
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
	betweennessCentrality<<<grid,block, sizeof(int) * elements * MAX_DEGREE>>>(elements, elements, d_edge, 0x0);
	cudaError_t error = cudaGetLastError();
	cout << cudaGetErrorString(error) << endl;
	
	int* h_mem = (int*)malloc(sizeof(int) * elements);
	cudaMemcpy(h_mem, d_mem, sizeof(int) * elements, cudaMemcpyDeviceToHost);
	
	cout<<elements<<endl;
	
	//cudaProfilerStop();
	
	cudaDeviceReset();
	
	return 0;
}