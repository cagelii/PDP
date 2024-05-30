#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <mpi.h>

int main(int argc, char **argv){
    /*Parallel implementaion to solve the poisson equation in two
    dimensions using finite differences*/

    int id, p;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &id);
    MPI_Comm_size(MPI_COMM_WORLD, &p);
    if (argc != 2 && id == 0){
        printf("Usage: ./project n");
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
    int n = atoi(argv[1]);
    int base_chunk = n / p;
    int rest = n % p;
    int chunk = base_chunk + (id < rest ? 1 : 0);

    double* d = (double *) malloc(sizeof(double)*chunk*n); //divides so that each process gets "chunk" rows (each of length n)
    double* g = (double *) malloc(sizeof(double)*chunk*n);
    double* u = (double *) calloc(chunk*n, sizeof(double)); //initialize u as 0
    double* q = (double *) malloc(sizeof(double)*chunk*n);

    MPI_Barrier(MPI_COMM_WORLD);
    double start = MPI_Wtime();

    double h = 1.0/(n+1);
    double q0, q1;
    double qlocal = 0.0;
    double x, y, bx, by;
    int k = 0; int kn;
    double h2 = 2*h*h;

    for (int j = id * base_chunk + (id < rest ? id : rest); j < (id + 1) * base_chunk + (id < rest ? id + 1 : rest); j++){ // initialize g and d, "natural" ordering according to wiki
        y = (j+1)*h; //interior point
        by = y*(1-y);
        kn = k*n;
        for (int i = 0; i < n; i++){ 
            x = (i+1) * h; //interior point
            bx = x*(1-x);
            d[kn+i] = h2*(bx + by); // d_ij
            g[kn+i] = -d[kn+i];
            qlocal += g[kn+i]*g[kn+i];
        }
        k++;
    }
    
    MPI_Allreduce(&qlocal, &q0, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    double* uppervals = (double *) malloc(sizeof(double)*n); //overlapping values needed for stencil operation, technically first and last rank doesn't need both
    double* lowervals = (double *) malloc(sizeof(double)*n);
    double tau, beta, temp, enTillTemp;
    int iterations = 0;

    do { //conjugate gradient
        //4  
        MPI_Request req[4];
        if (id != 0){
            MPI_Isend(d, n, MPI_DOUBLE, id - 1, 1, MPI_COMM_WORLD, &req[0]);
            MPI_Irecv(lowervals, n, MPI_DOUBLE, id - 1, 0, MPI_COMM_WORLD, &req[1]);
        }
        if (id != p - 1){
            MPI_Irecv(uppervals, n, MPI_DOUBLE, id + 1, 1, MPI_COMM_WORLD, &req[2]);
            MPI_Isend(&d[n * (chunk - 1)], n, MPI_DOUBLE, id + 1, 0, MPI_COMM_WORLD, &req[3]);
        }

        for (int j = 1; j < chunk-1; j++){ //middle points can be calculated without overlapping values
            q[j*n] = 4*d[j*n] - d[j*n+1] - d[(j-1)*n] - d[(j+1)*n]; //first element has no element to the left
            for (int i = 1; i < n-1; i++){
                q[j*n+i] = 4*d[j*n+i] - d[j*n+i-1] - d[j*n+i+1] - d[(j-1)*n+i] - d[(j+1)*n+i]; //middle points
            }
            q[(j+1)*n-1] = 4*d[(j+1)*n-1] - d[(j+1)*n-2] - d[(j+2)*n-1] - d[j*n-1]; // last element has no element to the right
        }
        
        if (id != 0) { //ensure all processes has got the overlapping values
            MPI_Wait(&req[0], MPI_STATUS_IGNORE);
            MPI_Wait(&req[1], MPI_STATUS_IGNORE);
        }
        if (id != p - 1) {
            MPI_Wait(&req[2], MPI_STATUS_IGNORE);
            MPI_Wait(&req[3], MPI_STATUS_IGNORE);
        }
        
        for (int i = 0; i < n; i++){ //first row
            q[i] = 4*d[i]-(i ? d[i-1]:0) - (i != n-1 ? d[i+1]:0) - (id ? lowervals[i]:0) - d[i+n];
        }
        
        int upperIndex = 0;
        for (int i = (chunk-1)*n; i < chunk*n; i++){ // last row
            q[i] = 4*d[i] - (i == (chunk-1)*n ? 0:d[i-1]) - (i == chunk*n-1 ? 0:d[i+1]) - (id == p-1 ? 0 : uppervals[upperIndex++]) - d[i-n];
        }
        
        //5
        temp = 0.0;
        for (int i = 0; i < n*chunk; i++){
            temp += d[i]*q[i];
        }
        MPI_Allreduce(&temp, &enTillTemp, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        tau = q0/enTillTemp;
        
        //6, 7 and 8
        temp = 0;
        for (int i = 0; i < n*chunk; i++){//these steps can be clumped into one for loop because they do not depend on each other
            u[i] += tau*d[i];
            g[i] += tau*q[i];
            temp += g[i]*g[i];
        }
        MPI_Allreduce(&temp, &q1, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);    

        //9
        beta = q1/q0;

        //10
        for (int i = 0; i < chunk*n; i++){
            d[i] = -g[i] + beta*d[i];
        }

        //11
        q0 = q1;
        
        iterations++;
    } while (iterations<200); // did a while loop to experiment with a converging stopping condition
    double executiontime = MPI_Wtime()-start;
    free(uppervals); free(lowervals); free(g); free(q); free(d); free(u);
   
    double maxtime;
    MPI_Reduce(&executiontime, &maxtime, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    if (id == 0){
        // printf("Time: %f\n", maxtime);
        printf("%.3e\n", sqrt(q0));
    }

    MPI_Finalize();
}