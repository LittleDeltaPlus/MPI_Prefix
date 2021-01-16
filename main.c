#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <mpi.h>
#include <math.h>

struct metadata{
    int length;
    int N;
    int levels;
    int runningTotal;
    int loops;
};

void lin_Scan(const int* input_array, int* output_array,const int sz){
    output_array[0] = input_array[0];
    for (int i = 1; i < sz; ++i) {
        output_array[i]=output_array[i-1]+input_array[i];
    }
}

void gen_sequential(int* input_array, const int size, const int upBound){
    for (int i = 0; i < size; ++i) {
        input_array[i]=i+1;
    }
    for (int i = size; i < upBound; ++i) {
        input_array[i] = 0;
    }
}

void gen_rand(int* input_array, const int size, const int upBound){
    srand(time(NULL));
    for (int i = 0; i < size; ++i) {
        input_array[i] = rand();
    }
    for (int i = size; i < upBound; ++i) {
        input_array[i] = 0;
    }
}

int get_UI(){
    char val[10], *eptr; //Input string for scanf
    int N; //The size of the array to be evaluated

    printf("Enter Array Length:");
    fflush(stdout);
    scanf("%s", &val); //Get input from user
    N = strtol(val, &eptr, 10); //Sanitize input
//    N = 16;

    // Error out if array is too small
    if(N == 0){
        printf("\nplease choose a number larger than 0\n");
        fflush(stdout);
        return 0;
    }
}

int main(int argc, char **argv)
{
    MPI_Init(&argc, &argv);

    //Init Local Variables
    struct metadata meta; //MetaData Struct
    int rankElement, rankOutput, buf; //Local variables for scan
    int *comp_domain, *output; //input and output Arrays
    const int root = 0; //make root 0

    //Get MPI info
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    int id = rank+1;

    int user_tag = rank;

    if(rank==root){
        //Get Valid User Input
        meta.N = 0;
        while (meta.N == 0){
            meta.N = get_UI();
        }

        //pad array to nearest 2^x
        meta.levels = ceil(log10(meta.N) / log10(2));
        meta.length = ceil(pow(2, meta.levels));
        meta.runningTotal = 0;
        size = 16;
        meta.loops = ceil((double)meta.length/(double)size);

        if(meta.loops > 1){
            meta.levels = ceil(log10(size* meta.loops)/log10(2));
            meta.length = ceil(pow(2, meta.levels));
        }


        // Allocate Input and Output array sizes
        comp_domain = malloc(meta.length * sizeof(int));
        output = malloc(meta.length * sizeof(int));

        //Generate Input array
        gen_sequential(comp_domain, meta.N, meta.length);
        //Populate Output with input
        memcpy(output, comp_domain, meta.length * sizeof(int));
        //If only one thread, evaluate linearly
        if (size == 1){
            lin_Scan(comp_domain, output, meta.N);
            for (int i = 0; i < meta.length; ++i) {
                printf("%d, %d\n", comp_domain[i], output[i]);
                fflush(stdout);
            }
            MPI_Finalize();
            return 0;
        }
    }



    //giveData
    int *loopOutput;
    loopOutput = calloc(size, sizeof(int));

    MPI_Bcast(&meta, 5, MPI_INT, 0, MPI_COMM_WORLD);
    for (int i = 0; i < meta.loops; ++i) {

        if(rank == root){
            memcpy(loopOutput, output +i*size, size * sizeof(int));
            loopOutput[0] += meta.runningTotal;
        }

        MPI_Scatter(loopOutput, 1, MPI_INT, &rankElement, 1, MPI_INT, 0, MPI_COMM_WORLD);

        rankOutput = rankElement;

        //upPhase
        for (int j = 1; j <= meta.levels; j++) {
            rankElement = rankOutput;
            // Find the 'step' for this iteration, and the corresponding data
            int iterator = ceil(pow(2, j));
            int diff = ceil(pow(2, j - 1));

            //Get work Based off rank
            if (id % iterator == 0 && rank - diff >= 0) {
                //recv
                MPI_Recv(&buf, 1, MPI_INT, rank - diff, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                rankOutput = rankElement + buf;
            } else if (id % iterator == diff && rank + diff < size) {
                //send
                MPI_Send(&rankElement, 1, MPI_INT, rank + diff, 0, MPI_COMM_WORLD);
            } else {
                //wait
            }
            //Ensure all finished
            MPI_Barrier(MPI_COMM_WORLD);
        }

        //downPhase
        for (int j = meta.levels - 1; j >= 1; j--) {
            rankElement = rankOutput;
            // Find the 'step' for this iteration, and the corresponding data
            int iterator = ceil(pow(2, j));
            int diff = ceil(pow(2, j - 1));

            //Get work Based off rank
            if (id % iterator == 0 && rank + diff < size) {
                //recv
                MPI_Send(&rankElement, 1, MPI_INT, rank + diff, 0, MPI_COMM_WORLD);
            } else if (id % iterator == diff && rank - diff >= 0 && id > diff) {
                //send
                MPI_Recv(&buf, 1, MPI_INT, rank - diff, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                rankOutput = rankElement + buf;

            } else {
                //wait
            }
            //ensure all finished
            MPI_Barrier(MPI_COMM_WORLD);
        }

        //Collect all rankOutputs into loopOutput
        MPI_Gather(&rankOutput, 1, MPI_INT, loopOutput, 1, MPI_INT, 0, MPI_COMM_WORLD);

        if(rank == root){
            memcpy(output + i*size, loopOutput, size * sizeof(int));
            meta.runningTotal = loopOutput[size-1];
        }
    }

    if(rank == root){
        //generate linear eq. for testing
        int *linear_output = malloc(meta.length * sizeof(int));
        lin_Scan(comp_domain, linear_output, meta.length);

        //List input vs. MPI vs. Linear
        for (int i = 0; i < meta.N; ++i) {
            printf("%d, %d, %d\n", comp_domain[i], output[i] ,linear_output[i]);
            fflush(stdout);
        }
    }

    //Exit
    MPI_Finalize();
    return 0;
}



/**
 * 1. get size
 * 2. threads = size
 * 3. get upper bounds next 2^x
 * 4. gen data
 * 5. broadcast meta-data
 * 6. send each thread it's own data
 * 7. (par) foreach level in algotitm -> {
 *          am I receiving?
 *              yes - getData
 *              no - sendData or Wait
 *          on receive -> {
 *              perform update
 *          }
 * }
 *
 */
