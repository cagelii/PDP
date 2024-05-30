#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <mpi.h>

#define writeOutput 1

int numValues;
int nreturnValues;

int read_input(const char *file_name, int **values) {
	FILE *file;
	if (NULL == (file = fopen(file_name, "r"))) {
		perror("Couldn't open input file");
		return -1;
	}
	int num_values;
	if (EOF == fscanf(file, "%d", &num_values)) {
		perror("Couldn't read element count from input file");
		return -1;
	}
	if (NULL == (*values = malloc(num_values * sizeof(int)))) {
		perror("Couldn't allocate memory for input");
		return -1;
	}
	for (int i=0; i<num_values; i++) {
		if (EOF == fscanf(file, "%d", &((*values)[i]))) {
			perror("Couldn't read elements from input file");
			return -1;
		}
	}
	if (0 != fclose(file)){
		perror("Warning: couldn't close input file");
	}
	return num_values;
}

int write_output(char *file_name, const int *output, int num_values) {
	FILE *file;
	if (NULL == (file = fopen(file_name, "w"))) {
		perror("Couldn't open output file");
		return -1;
	}
	for (int i = 0; i < num_values; i++) {
		if (0 > fprintf(file, "%d ", output[i])) {
			perror("Couldn't write to output file");
		}
	}
	if (0 > fprintf(file, "\n")) {
		perror("Couldn't write to output file");
	}
	if (0 != fclose(file)) {
		perror("Warning: couldn't close output file");
	}
	return 0;
}

int cmpfunc (const void* a, const void* b) {
   return (*(int*)a - *(int*)b);
}

void mergeArrays(int *arr1, int *arr2, int n1, int n2, int *arr3)
{
    int i = 0, j = 0, k = 0;
    while (i<n1 && j <n2)
    {
        if (arr1[i] < arr2[j])
            arr3[k++] = arr1[i++];
        else
            arr3[k++] = arr2[j++];
    }

    while (i < n1)
        arr3[k++] = arr1[i++];

    while (j < n2)
        arr3[k++] = arr2[j++];
}

int *sort(int *list, int elements, int pivotStrategy, int commsize, int localrank, MPI_Comm communicator, int ogrank){
	if (commsize == 1){ // base case
		nreturnValues = elements;
		return list;
	}
	
	int pivotElement;
	int median = elements ? list[elements/2] : 0;
	
	switch (pivotStrategy){ // choose pivot for the group
	case 1:
	if (localrank==0) {pivotElement = median;}
		MPI_Bcast(&pivotElement, 1, MPI_INT, 0, communicator);
		break;
	case 2:
		int *medians = (int *) malloc(commsize*sizeof(int));
		MPI_Allgather(&median, 1, MPI_INT, medians, 1, MPI_INT, communicator);
		qsort(medians, commsize, sizeof(int), cmpfunc);
		pivotElement = medians[commsize/2];
		free(medians);
		break;
	case 3:
		MPI_Allreduce(&median, &pivotElement, 1, MPI_INT, MPI_SUM, communicator);
		pivotElement /= commsize;
		break;
	}
	
	int *smaller = (int *) malloc(sizeof(int)*elements); // elements smaller than the pivot
	int i;
	for (i = 0; i<elements; i++){
		if (list[i]>pivotElement) {break;}
		smaller[i] = list[i];
	}

	smaller = (int *) realloc(smaller, sizeof(int)*i); //reallocate so that it doesn't use unneccessary space
	int *bigger = (int *) malloc(sizeof(int)*(elements-i));
	for (int j = 0; j<(elements-i); j++){
		bigger[j] = list[j+i];
	}

	int color = localrank*2/commsize; // divide the communicator into two new ones
	// printf("rank: %d, color: %d\n", ogrank, color);
	int newrank, newsize;
	MPI_Comm newcomm;
	MPI_Comm_split(communicator, color, localrank, &newcomm);
	MPI_Comm_rank(newcomm, &newrank);
	MPI_Comm_size(newcomm, &newsize);

	int elementsToReceive;
	int *received, *result;
	
	if (color == 1){ // color 1 receives values larger than the pivot
		int target = localrank-(commsize/2);
		MPI_Send(&i, 1, MPI_INT, target, 0, communicator); // need to send numbers to receive
		MPI_Recv(&elementsToReceive, 1, MPI_INT, target, 0, communicator, MPI_STATUS_IGNORE);

		MPI_Send(smaller, i, MPI_INT, target, 0, communicator);

		received = (int *) malloc(sizeof(int)*elementsToReceive);
		list = (int *) realloc(list, sizeof(int)*(elementsToReceive+elements-i));
		MPI_Recv(received, elementsToReceive, MPI_INT, target, 0, communicator, MPI_STATUS_IGNORE);
		mergeArrays(received, bigger, elementsToReceive, elements-i, list);
		
		free(received);
				
		free(smaller); free(bigger); 
		result = sort(list, elementsToReceive+elements-i, pivotStrategy, newsize, newrank, newcomm, ogrank);
	}
	else {
		int elementsToSend = elements-i;
		int target = localrank+(commsize/2);
		MPI_Recv(&elementsToReceive, 1, MPI_INT, target, 0, communicator, MPI_STATUS_IGNORE);
		MPI_Send(&elementsToSend, 1, MPI_INT, target, 0, communicator);

		received = (int *) malloc(sizeof(int)*elementsToReceive);
		list = (int *) realloc(list, sizeof(int)*(elementsToReceive+i));
		MPI_Recv(received, elementsToReceive, MPI_INT, target, 0, communicator, MPI_STATUS_IGNORE);
		mergeArrays(received, smaller, elementsToReceive, i, list);
		free(received);

		MPI_Send(bigger, elementsToSend, MPI_INT, target, 0, communicator);

		free(smaller); free(bigger); 
		result = sort(list, elementsToReceive+i, pivotStrategy, newsize, newrank, newcomm, ogrank);
	}
	MPI_Comm_free(&newcomm);
	return result;
}

int main(int argc, char **argv){
    int rank, numProcesses;
    int* input;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &numProcesses);
	if (argc != 4){
        	printf("Usage: ./quicksort input_file output_file pivot_selection\n");
        	MPI_Abort(MPI_COMM_WORLD, 1);
    	}

    char *inputName, *outputName;
    if (rank == 0){
        inputName = argv[1];
        outputName = argv[2];
        if (0 > (numValues = read_input(inputName, &input))) {
			MPI_Abort(MPI_COMM_WORLD, 1);
		}
    }
	double start = MPI_Wtime();
	MPI_Bcast(&numValues, 1, MPI_INT, 0, MPI_COMM_WORLD);
	int rest = numValues%numProcesses;
	int chunksize = rest > rank ? numValues/numProcesses+1 : numValues/numProcesses; //extend block if necessary
	int *array = (int *) malloc(sizeof(int)*chunksize);
	int *sendcount = (int*) malloc(sizeof(int)*numProcesses);
	int *disps = (int*) malloc(sizeof(int)*numProcesses); // scattering process
	int temp = 0;
	for (int i = 0; i < numProcesses; i++){
		sendcount[i] = rest > i ? numValues/numProcesses+1 : numValues/numProcesses;
		disps[i] = temp;
		temp += sendcount[i];
	}
	MPI_Scatterv(input, sendcount, disps, MPI_INT, array, chunksize, MPI_INT, 0, MPI_COMM_WORLD);
	free(disps); free(sendcount); free(input);

	int pivotMethod = atoi(argv[3]);
	if (pivotMethod != 1 && pivotMethod != 2 && pivotMethod != 3){
		printf("Invalid pivot method\n");
		MPI_Abort(MPI_COMM_WORLD, 1);
	}
	
	MPI_Barrier(MPI_COMM_WORLD);

    qsort(array, chunksize, sizeof(int), cmpfunc);
 	int *result = sort(array, chunksize, pivotMethod, numProcesses, rank, MPI_COMM_WORLD, rank);
	double time = MPI_Wtime()-start;
	if (rank != 0){ // gather data, could've probably used gatherv
		MPI_Send(&nreturnValues, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
		MPI_Send(result, nreturnValues, MPI_INT, 0, 0, MPI_COMM_WORLD);
	}
	else{
		int j;
		input = (int *) malloc(sizeof(int)*numValues);
		if (input == NULL){printf("can't allocate output"); MPI_Abort(MPI_COMM_WORLD, 1);}
		for (j = 0; j < nreturnValues; j++){
			input[j] = result[j];
		}
		for (int i = 1; i < numProcesses; i++){
			int values;
			MPI_Recv(&values, 1, MPI_INT, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			MPI_Recv(&input[j], values, MPI_INT, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			j += values;
		}	
	}
	double longestTime;
	MPI_Reduce(&time, &longestTime, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
	if (rank == 0){
		#if writeOutput
		if (0 != write_output(outputName, input, numValues)) {
			return 2;
		}
		#endif
		free(input);
		printf("%f\n", time);
	}
	free(result);
	MPI_Finalize();
}
