#include "bc_cuda.h"
#include <string.h>
#include <stdio.h>
#include <sys/time.h>
#include <stdlib.h>
#include <mpi.h>

int main(int argc, char* argv[])
{

	int numprocs, rank;
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);


	FILE *grFile;
	char buff[1024];
	int edgeCnt = 0;
	int numVert;
	int numEdge;
	int* h_edge;
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

	int numSegment = 1;
	int seg = 0;
	float *bc, *bc_recv;

	bc = (float*)malloc(sizeof(float) * numVert);
	bc_recv = (float*)malloc(sizeof(float) * numVert);
	betweennessCentralityCU(numVert, numEdge, h_edge, bc, numprocs, rank);
	MPI_Reduce(bc, bc_recv, numVert, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);

	MPI_Finalize();

	int i;
	for(i = 0; i < numVert; i++)
	{
		printf("%f\n", bc_recv[i]);
	}
}
