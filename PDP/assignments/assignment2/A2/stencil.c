#include "stencil.h"

int main(int argc, char **argv) {	
	int rank, numProcesses;
	int chunksize, num_values;
	double maxTime, *input;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &numProcesses);
	MPI_Status stat;

	if (4 != argc) {
		if (rank==0){
			printf("Usage: stencil input_file output_file number_of_applications\n");
			MPI_Abort(MPI_COMM_WORLD,1);
		}
	}

	char *input_name, *output_name;

	int num_steps = atoi(argv[3]);
	if (rank==0){
		input_name = argv[1];
		output_name = argv[2];

		// Read input file
		if (0 > (num_values = read_input(input_name, &input))) {
			MPI_Abort(MPI_COMM_WORLD, 1);
		}
		if (num_values%numProcesses != 0){
			printf("Array length must be divisible by number of processes!\n");
			MPI_Abort(MPI_COMM_WORLD,1);
		}
	}
	MPI_Bcast(&num_values, 1, MPI_INT, 0, MPI_COMM_WORLD);
	chunksize = num_values/numProcesses;
	double *array = malloc(sizeof(double)*chunksize);
	MPI_Scatter(input, chunksize, MPI_DOUBLE, array, chunksize, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	
	// Stencil values
	double h = 2.0*PI/num_values;
	const int STENCIL_WIDTH = 5;
	const int EXTENT = STENCIL_WIDTH/2;
	const double STENCIL[] = {1.0/(12*h), -8.0/(12*h), 0.0, 8.0/(12*h), -1.0/(12*h)};

	if (STENCIL_WIDTH>chunksize){
		if (rank==0){
			printf("Too big stencil on each chunk!\n");
			MPI_Abort(MPI_COMM_WORLD, 1);
		}
	}

	// Start timer
	MPI_Barrier(MPI_COMM_WORLD); // ensure "fair" start, I get really bad results if I don't do this
	double start = MPI_Wtime();

	// Allocate data for result
	double *output;
	if (NULL == (output = malloc(chunksize * sizeof(double)))) {
		perror("Couldn't allocate memory for output");
		MPI_Abort(MPI_COMM_WORLD, 1);
	}
	// Repeatedly apply stencil
	double edgeValues[2*EXTENT];
	for (int s=0; s<num_steps; s++) {
		double timeXd = MPI_Wtime();
		int left = rank%numProcesses ? rank-1:numProcesses-1;
		int right = (rank+1)%numProcesses;
		
		MPI_Request sendrecv_request[4]; //use non blocking communication so that processes doesn't wait around
		MPI_Isend(&array[chunksize-EXTENT], EXTENT, MPI_DOUBLE, right, 0, MPI_COMM_WORLD, &sendrecv_request[0]);
		MPI_Irecv(edgeValues, EXTENT, MPI_DOUBLE, left, 0, MPI_COMM_WORLD, &sendrecv_request[1]);
		MPI_Isend(array, EXTENT, MPI_DOUBLE, left, 0, MPI_COMM_WORLD, &sendrecv_request[2]);
		MPI_Irecv(&edgeValues[EXTENT], EXTENT, MPI_DOUBLE, right, 0, MPI_COMM_WORLD, &sendrecv_request[3]);

		// the middle
		for (int i = EXTENT; i < chunksize - EXTENT; i++) {
			double result = 0;
			for (int j = 0; j < STENCIL_WIDTH; j++) {
				int index = i - EXTENT + j;
				result += STENCIL[j] * array[index];
			}
			output[i] = result;
		}
		// wait for non-blocking communication to complete
		MPI_Waitall(4, sendrecv_request, MPI_STATUSES_IGNORE);
		
		for (int i = 0; i<EXTENT; i++){
			double result = 0;
			for (int j = 0; j<STENCIL_WIDTH; j++){
				double prod = (j<EXTENT-i ? edgeValues[j+i] : array[i+j-EXTENT]);
				result += STENCIL[j]*prod;
			}
			output[i] = result;
		}
		
		for (int i = chunksize-EXTENT; i < chunksize; i++){
			double result = 0;
			for (int j = 0; j < STENCIL_WIDTH; j++){
				double prod = chunksize-i+EXTENT>j ? array[i-EXTENT+j] : edgeValues[j-chunksize+i];
				result += STENCIL[j]*prod;
			}
			output[i] = result;
		}
		
		// Swap input and output
		if (s < num_steps-1) {
			double *tmp = array;
			array = output;
			output = tmp;
		}
	}
	// Stop timer
	double my_execution_time = MPI_Wtime() - start;
	MPI_Reduce(&my_execution_time, &maxTime, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
	// Write result
	if (rank==0){
		printf("%f\n", maxTime); 
	}
	MPI_Gather(output, chunksize, MPI_DOUBLE, input, chunksize, MPI_DOUBLE, 0, MPI_COMM_WORLD);
#ifdef PRODUCE_OUTPUT_FILE
	if (rank == 0){
		if (0 != write_output(output_name, input, num_values)) {
			return 2;
		}
	}
#endif
	// Clean up
	free(array);
	free(output);
	if (rank == 0){
		free(input);
	}
	
	MPI_Finalize();
	return 0;
}


int read_input(const char *file_name, double **values) {
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
	if (NULL == (*values = malloc(num_values * sizeof(double)))) {
		perror("Couldn't allocate memory for input");
		return -1;
	}
	for (int i=0; i<num_values; i++) {
		if (EOF == fscanf(file, "%lf", &((*values)[i]))) {
			perror("Couldn't read elements from input file");
			return -1;
		}
	}
	if (0 != fclose(file)){
		perror("Warning: couldn't close input file");
	}
	return num_values;
}

int write_output(char *file_name, const double *output, int num_values) {
	FILE *file;
	if (NULL == (file = fopen(file_name, "w"))) {
		perror("Couldn't open output file");
		return -1;
	}
	for (int i = 0; i < num_values; i++) {
		if (0 > fprintf(file, "%.4f ", output[i])) {
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