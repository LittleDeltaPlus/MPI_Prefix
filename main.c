

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <mpi.h>
#include <math.h>

/**
 * Data About The Scan Array And MPI
 */
struct metadata{
    int N; ///The Chosen Array Length
    int length; ///The Length of the array with padding
    int levels; ///The Number of Steps Need to complete the Up Phase
    int loops;  ///The Number of Loops required to sum all values
    int size;   ///The Number of Active Threads
    int root;   ///The Root Thread
};
/**
 * Local Data For Each Thread
 */
struct thread {
    int rank; ///Thread's MPI Rank
    int runningTotal; ///Running Total Used by thread 0
    int element; ///Output Value for each thread
    int buffer; ///Space for a received element
};
/**
 * Handles CLI argument calls, allows sequential numbers to be used
 * @param [in] argc
 * @param [in] argv
 * @return AppMode 0: exit 1: Use Sequential Numbers 2: Use Random Numbers
 */
int handleArgs(int argc, char **argv);

/**
 * A linear implementation of the prefix scan
 * @param [in] input_array Values to be summed
 * @param [out] output_array The summed values
 * @param [in] sz How many values are present
 */
void LinScan(const int* input_array, int* output_array, int sz);
/**
 * Generates an array of sequential numbers starting at 1, Used to Generate Predictable Test Data
 * @param [out] output_array The List of sequential Numbers
 * @param [in] size The Size of the output array
 * @param [in] upBound The size + any padding needed
 */
void GenSequential(int* output_array, int size, int upBound);
/**
 * Generates a list of Random numbers
 * @param [out] output_array The List of generated Random Numbers
 * @param [in] size The Size of the output array
 * @param [in] upBound The size + any padding needed
 */
void GenRand(int* output_array, int size, int upBound);
/**
 * Use MPI metadata to calculate the number of levels & padding needed
 * @param [in,out] metaData
 */
void GetArraySize(struct metadata *metaData);
/**
 * An Alternative Implementation of Prefix Scan Sum
 * @param [in] metaData Data about the size of the scanned array & MPI Size
 * @param [in,out] threadData Local variables for each thread
 * @param [out] output The scanned array
 */
void PrefixScan(struct metadata *metaData, struct thread *threadData, int *output);
/**
 * The Up Phase of the Prefix scan
 * @param [in] metaData Data about the size of the scanned array & MPI Size
 * @param [in,out] threadData Local variables for each thread, including that thread's momentary total
 */
void UpPhase(struct metadata *metaData, struct thread *threadData);
/**
 * The Down Phase of the Prefix scan
 * @param [in] metaData Data about the size of the scanned array & MPI Size
 * @param [in,out] threadData Local variables for each thread, including that thread's output
 */
void DownPhase(struct metadata *metaData, struct thread *threadData);
/**
 * Gets a sanitised array length from the user
 * @return N, The Chosen Array Length
 */
int GetUI();



int main(int argc, char **argv)
{
    MPI_Init(&argc, &argv);

    //Handle Arguments, perform accordingly
    int appMode = handleArgs(argc, argv);
    //Broadcast in-case of Exit
    MPI_Bcast(&appMode, 1, MPI_INT, 0, MPI_COMM_WORLD);
    if(appMode == 0){
        MPI_Finalize();
        return 1;
    }


    //Init Local Variables
    struct metadata meta; //MetaData Struct
    struct thread threadData; //local variables
    int *comp_domain, *output; //input and output Arrays
    meta.root = 0;//Set root 0
    threadData.runningTotal = 0;//Start counting from 0

    //Get MPI info
    int size, rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    threadData.rank = rank;


    if(threadData.rank==meta.root){
        meta.size = size;
        //Calculate the size of the computation
        GetArraySize(&meta);

        // Allocate Input and Output array sizes
        comp_domain = malloc(meta.length * sizeof(int));
        output = malloc(meta.length * sizeof(int));

        //Generate Input array
        if(appMode == 1){
            GenSequential(comp_domain, meta.N, meta.length);
        } else {
            GenRand(comp_domain, meta.N, meta.length);
        }
        //Populate Output with input
        memcpy(output, comp_domain, meta.length * sizeof(int));
    }

    //If only one thread, evaluate linearly
    if (meta.size == 1 && threadData.rank == meta.root){

        //perform a single threaded prefix scan
        LinScan(comp_domain, output, meta.N);
        //print output to screen
        for (int i = 0; i < meta.N; ++i) {
            printf("%d, %d\n", comp_domain[i], output[i]);
            fflush(stdout);
        }
    }
    else {
        //Perform Multi-threaded Prefix scan
        PrefixScan(&meta, &threadData, output);

        //Output to screen
        if(threadData.rank == meta.root){
            //generate linear eq. for testing

            int *linear_output = malloc(meta.length * sizeof(int));
            LinScan(comp_domain, linear_output, meta.length);

            //List input vs. MPI vs. Linear
            for (int i = 0; i < meta.N; i++) {
                printf("%d, %d, %d\n", comp_domain[i], output[i]);
                fflush(stdout);

                if(output[i] != linear_output[i]){
                    printf("error @ i = %d: %d != %d\nstopping\n"
                           "please raise an issue at:"
                           " https://github.com/LittleDeltaPlus/MPI_Prefix/issues\n", i, output[i] ,linear_output[i]);
                    fflush(stdout);
                    i = meta.N;
                }
            }
        }
    }
    //Exit
    MPI_Finalize();
    return 0;
}

void GetArraySize(struct metadata *metaData){
    //Get Valid User Input
    metaData->N = 0;
    while (metaData->N == 0){
        metaData->N = GetUI();
    }

    //pad array to nearest 2^x
    metaData->levels = ceil(log10(metaData->N) / log10(2));
    metaData->length = ceil(pow(2, metaData->levels));
    metaData->loops = ceil((double)metaData->length/(double)metaData->size);

    //If fit N to a multiple of size
    if(metaData->loops > 1){
        metaData->levels = ceil(log10(metaData->size* metaData->loops)/log10(2));
        metaData->length = ceil(pow(2, metaData->levels));
    }
}

void PrefixScan(struct metadata *metaData, struct thread *threadData, int *output) {
    //giveData
    MPI_Bcast(metaData, 6, MPI_INT, 0, MPI_COMM_WORLD);
    int *loopOutput;
    for (int i = 0; i < metaData->loops; ++i) {

        if (threadData->rank == metaData->root) {
            loopOutput = malloc(metaData->size * sizeof(int));
            memcpy(loopOutput, output + i * metaData->size, metaData->size * sizeof(int));
            loopOutput[0] += threadData->runningTotal;
        }

        //Scatter data for this loop
        MPI_Scatter(loopOutput, 1, MPI_INT, &threadData->element, 1, MPI_INT, 0, MPI_COMM_WORLD);
        //Up Phase
        UpPhase(metaData, threadData);
        //Dow Phase
        DownPhase(metaData, threadData);
        //Collect all rankOutputs into loopOutput
        MPI_Gather(&threadData->element, 1, MPI_INT, loopOutput, 1, MPI_INT, 0, MPI_COMM_WORLD);

        //Copy evaluation to output, update running total
        if (threadData->rank == metaData->root) {
            memcpy(output + i * metaData->size, loopOutput, metaData->size * sizeof(int));
            threadData->runningTotal = loopOutput[metaData->size - 1];
        }
    }
}

void UpPhase(struct metadata *metaData, struct thread *threadData){
    int id = threadData->rank+1;
    for (int j = 1; j <= metaData->levels; j++) {
        // Find the 'step' for this iteration, and the corresponding data
        int iterator = ceil(pow(2, j));
        int diff = ceil(pow(2, j - 1));

        //Get work Based off rank
        if (id % iterator == 0 && threadData->rank - diff >= 0) {
            //recv
            MPI_Request recv_request;
            MPI_Irecv(&threadData->buffer, 1, MPI_INT, threadData->rank - diff, 0, MPI_COMM_WORLD, &recv_request);
            MPI_Wait(&recv_request, MPI_STATUS_IGNORE);
            threadData->element += threadData->buffer;
        } else if (id % iterator == diff && threadData->rank + diff < metaData->size) {
            //send
            MPI_Request send_request;
            MPI_Isend(&threadData->element, 1, MPI_INT, threadData->rank + diff, 0, MPI_COMM_WORLD, &send_request);
        }
    }
};

void DownPhase(struct metadata *metaData, struct thread *threadData){
    int id = threadData->rank+1;
    for (int j = metaData->levels - 1; j >= 1; j--) {
        // Find the 'step' for this iteration, and the corresponding data
        int iterator = ceil(pow(2, j));
        int diff = ceil(pow(2, j - 1));

        //Get work Based off rank
        if (id % iterator == 0 && threadData->rank + diff < metaData->size) {
            //send
            MPI_Request send_request;
            MPI_Isend(&threadData->element, 1, MPI_INT, threadData->rank + diff, 0, MPI_COMM_WORLD, &send_request);
        } else if (id % iterator == diff && threadData->rank - diff >= 0 && id > diff) {
            //recv
            MPI_Request recv_request;
            MPI_Irecv(&threadData->buffer, 1, MPI_INT, threadData->rank - diff, 0, MPI_COMM_WORLD, &recv_request);
            MPI_Wait(&recv_request, MPI_STATUS_IGNORE);
            threadData->element += threadData->buffer;
        }
    }
};


void LinScan(const int* input_array, int* output_array, const int sz){
    //start the count at the first value
    output_array[0] = input_array[0];
    for (int i = 1; i < sz; ++i) {
        //add the previous output value to the following input value
        output_array[i]=output_array[i-1]+input_array[i];
    }
}

void GenSequential(int* output_array, const int size, const int upBound){
    for (int i = 0; i < size; ++i) {
        //Set each number to a sequential value, starting at one
        output_array[i]= i + 1;
    }
    for (int i = size; i < upBound; ++i) {
        //pad if necessary
        output_array[i] = 0;
    }
}

void GenRand(int* output_array, const int size, const int upBound){
    //seed random number generator
    //TODO: Find Cryptographic compliant variant
    srand(time(NULL));
    for (int i = 0; i < size; ++i) {
        //set each 'useful' element to a random number
        output_array[i] = rand();
    }
    for (int i = size; i < upBound; ++i) {
        //pad if necessary
        output_array[i] = 0;
    }
}

int GetUI(){
    char val[10], *eptr; //Input string for scanf
    int N; //The size of the array to be evaluated

    printf("Enter Array Length:");
    fflush(stdout);
    scanf("%s", &val); //Get input from user
    N = strtol(val, &eptr, 10); //Sanitize input

    // Error out if array is too small
    if(N == 0){
        printf("\nplease choose a number larger than 0\n");
        fflush(stdout);
        return 0;
    }
    return N;
}

int handleArgs(int argc, char **argv){
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if (rank == 0){
        int genSequential = 2;
        size_t optind;
        char * eptr;
        for (optind = 1; optind < argc && argv[optind][0] == '-'; optind++) {
            switch (argv[optind][1]) {
                case 's': genSequential = 1; break;
                case 'r': genSequential = 2; break;
                case '?': genSequential = 0;
                    printf("Usage: %s [-rs?]\n"
                           "    -r : generate random sequence\n"
                           "    -s : generate sequential sequence\n"
                           "    -? : display this help message\n", argv[0]);
                    fflush(stdout); break;
                default:
                    fprintf(stderr, "Usage: %s [-rs?] (use -? for help)\n", argv[0]);
                    genSequential = 0;
            }
        }
        return genSequential;
    } else {
        return 2;
    }
}