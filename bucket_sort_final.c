#include "mpi.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

long INTERVAL;
long min;


int compareIntegers(const void* first, const void* second) {
    int x = *((int*)first), y =  *((int*)second);
    if (x == y) {
        return 0;
    }
    else if (x < y) {
        return -1;
    }
    else {
        return 1;
    }
}


int find_bucket(int num, int p_num){
	return ((-min + num)/INTERVAL)%p_num;
}

int main(int argc, char *argv[]){
	int myrank, P;
	int sub_count[4];

	if(argc != 2){
		printf("\nPlease include N, problem size\n");
		return 0;
	}


	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &P);
	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);	


	int N = atoi(argv[1]);
	int count[P];
	memset(count, 0, sizeof count);
	int list[N];
	int displs[P];
	int dist_list[N];
	int index[P];

	double t1 = MPI_Wtime();


	if(myrank == 0){

		int r;
		int max = 0, min = 0;
		for(int i = 0; i < N; i++){

			r = rand() % 10000;
			list[i] = r;
        
        	if(r > max) max = r;
        	if(r < min) min = r;

		}

		INTERVAL = ((max - min) / P) + 1;

		for (int i = 0; i < N; i++) {
			int bucket = find_bucket(list[i], P);
			count[bucket]++;
		}

	}

	double t2 = MPI_Wtime();


	MPI_Scatter(count, 1, MPI_INT, sub_count, 1, MPI_INT, 0, MPI_COMM_WORLD );
			

	int bucket_list[sub_count[0]];


	if(myrank == 0){

		displs[0] = 0;
		for(int j = 1; j < P; j++){
			displs[j] = count[j-1] + displs[j-1];
		}


		for(int j = 0; j < N; j++){
			int bucket = find_bucket(list[j], P);
			dist_list[displs[bucket] + index[bucket]] = list[j];
			index[bucket]++;
		}
	}

	double t3 = MPI_Wtime();



	MPI_Scatterv(dist_list, count, displs, MPI_INT, bucket_list, sub_count[0], MPI_INT, 0, MPI_COMM_WORLD);

	double t4 = MPI_Wtime();


	qsort(bucket_list, sub_count[0],sizeof(int), &compareIntegers);

	double t5 = MPI_Wtime();


	MPI_Gatherv(bucket_list,sub_count[0], MPI_INT, dist_list, count, displs, MPI_INT, 0, MPI_COMM_WORLD );

	double t6 = MPI_Wtime();

	if(myrank == 0){
		int sorted = 1;
		for(int k = 0; k < N - 2; k++){
			if(dist_list[k] > dist_list[k+1]){
				sorted = 0;
			}
		}

		if(sorted == 1) {
			printf("\nSORTING CORRECT\n");
		} 
		else {
			printf("\nSORTING NOT CORRECT\n");
		}
		
		
		printf("\nTime to Generate: %f\n", t2-t1);
		printf("\nTime to Bin: %f\n", t3-t2);
		printf("\nTime to Distribute: %f\n", t4-t3);
		printf("\nTime to Sort: %f\n", t5-t4);
		printf("\nTime to Gather: %f\n", t6-t5);
		printf("\nTotal Execution Time: %f\n", t6-t1);

	}

	MPI_Finalize();

}