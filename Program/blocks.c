/*
	CITS3402 Project 1 2016
	Name:			Ammar Abu Shamleh, Pradyumn Vij
	Student number: 21521274, 21469477
    Date:           September 2016
*/
#include "blocks.h"
//	Program Name
const char* programName;

int main(int argc, char** argv){
    
    //Initialize MPI
    MPI_Init(NULL, NULL);
    //Get comm rank and size
    int commRank;
    int commSize;
    MPI_Comm_rank(MPI_COMM_WORLD, &commRank);
    MPI_Comm_size(MPI_COMM_WORLD, &commSize);
    
    
	// Program name without /, cast to constant for file
	programName = (const char*) strrchr(argv[0], '/') + 1;

    //Time overall execution of program
    double programTimeStart = omp_get_wtime();
    

    /* READ IN MATRIX DATA AND KEYS (ALL PROCESSES DO THIS)*/

    //Allocate memory for matrix database
	double **dataMatrix = (double**)malloc(1*sizeof(double*));
    //Read in matrix, and record number of rows and columns in ROWS and COLS
    double startTimeForIO = omp_get_wtime();
    dataMatrix = readMatrix(DATA_FILE, dataMatrix);
    //Transpose data matrix
    double **transposedData = transposeMatrix(dataMatrix);
    //Free original matrix
    freeData(dataMatrix);

    // Allocate memory for keys
	long long *keyDatabase = malloc(ROWS*sizeof(long long));
	// Read in keys
	readKeys(KEY_FILE, keyDatabase);
    double timeForIO = omp_get_wtime() - startTimeForIO;

    //Set number of threads for parallel regions on each node
    omp_set_num_threads(NUM_THREADS);

    
    /* THESE VARIABLES ARE FOR TIMING EXECUTION AND RETURNING EXPERIMENTAL RESULTS */
	double timeForSequentialOptimisedBlockGeneration = 0;
    double timeForParallelOptimisedBlockGeneration = 0;
    double timeForSequentialOptimisedCollisionDetection = 0;
    double timeForParallelOptimisedCollisionDetection = 0;
    double timeForMPIOptimisedBlockGeneration = 0;
    double timeForMPIOptimisedCollisionDetection = 0;
    double startTime;
    double execTime;


    Block *blockDatabase;
    int numBlocks;
    Collision *collisions;
    int numCollisions;
    

    /* FIND ALL BLOCKS SEQUENTIALLY USING OPTIMISED CODE AND TIME EXECUTION */
    if(commRank == 0) {
        //First free previous block database and reset block counter
        printf("Finding all blocks using sequential optimised code\n");
        blockDatabase = malloc(sizeof(Block));
        numBlocks = 0;
        
        //Now, find all blocks using sequential optimised code
        startTime = omp_get_wtime();
        blockDatabase = findBlocksOptimised(blockDatabase, transposedData, keyDatabase, &numBlocks);
        execTime = omp_get_wtime() - startTime;
        timeForSequentialOptimisedBlockGeneration = execTime;
        printf("Sequential optimised block generation complete\n\n");
    }
    
    
    /* FIND ALL BLOCKS IN PARALLEL (OPENMP AND MPI) USING OPTIMISED CODE AND TIME EXECUTION */
    //BARRIER HERE, so slave threads don't start computation until master reaches this point (for fair benchmarking)
    MPI_Barrier(MPI_COMM_WORLD);
    if(commRank == 0) {
        freeBD(blockDatabase, numBlocks);
    }
    numBlocks = 0;
    
    //VARIABLES FOR SLAVE PROCESSES (will be used later in code, so they're not declared  inside slave region)
    Block *partialBlockDatabase;
    int numBlocksPartial;
    
    if(commRank == 0) {
        /* MASTER PROCESS*/
        printf("Finding all blocks using MPI optimised code\n");
        startTime = omp_get_wtime();
        blockDatabase = malloc(sizeof(Block));

        //Master process. All slave processes will perform their work and send their results
        for(int i=1; i<commSize; i++) {
            //The process will send back an array of blocks that it found between its columns. It will first send back the number of blocks it found
            int numBlocksPartial;
            MPI_Status status;
            //Read number of blocks slave process found
            MPI_Recv(&numBlocksPartial, 1, MPI_INT, i, 123, MPI_COMM_WORLD, &status);
            for(int j=0; j<numBlocksPartial; j++) {
                //Read current block being sent by slave process
                long long signature;
                int column;
                int rows[4];
                MPI_Recv(&signature, 1, MPI_LONG_LONG_INT, i, 123, MPI_COMM_WORLD, &status);
                MPI_Recv(&column, 1, MPI_INT, i, 123, MPI_COMM_WORLD, &status);
                MPI_Recv(rows, 4, MPI_INT, i, 123, MPI_COMM_WORLD, &status);
                //Store current block in block database
                blockDatabase = (Block *) realloc(blockDatabase, sizeof(Block) * (numBlocks+1));
                blockDatabase[numBlocks].signature = signature;
                blockDatabase[numBlocks].column = column;
                blockDatabase[numBlocks].rows[0] = rows[0];
                blockDatabase[numBlocks].rows[1] = rows[1];
                blockDatabase[numBlocks].rows[2] = rows[2];
                blockDatabase[numBlocks++].rows[3] = rows[3];
            }
        }
        
        //Record execution time
        execTime = omp_get_wtime() - startTime;
        timeForMPIOptimisedBlockGeneration = execTime;
        printf("MPI optimised block generation complete\n\n");
    }
    
    else {
        //Slave process. Compute start and end values for work to be performed
        int chunkSize = COLS / (commSize-1);
        int start = (commRank-1) * chunkSize;
        int end = start + chunkSize;
        if(end > COLS) end = COLS-1;
        
        //The number of blocks this process finds
        numBlocksPartial = 0;
        
        //Now, find all blocks between these columns
        
        partialBlockDatabase = (Block *) malloc(sizeof(Block));
        partialBlockDatabase = findAllBlocksOptimisedMPI(partialBlockDatabase, transposedData, keyDatabase, &numBlocksPartial, start, end);
        
        //Send all located blocks back to master. First send number of blocks found
        MPI_Send(&numBlocksPartial, 1, MPI_INT, 0, 123, MPI_COMM_WORLD);
        for(int j=0; j<numBlocksPartial; j++) {
            //Send the block's signature, column, and row array
            MPI_Send(&partialBlockDatabase[j].signature, 1, MPI_LONG_LONG_INT, 0, 123, MPI_COMM_WORLD);
            MPI_Send(&partialBlockDatabase[j].column, 1, MPI_INT, 0, 123, MPI_COMM_WORLD);
            MPI_Send(partialBlockDatabase[j].rows, 4, MPI_INT, 0, 123, MPI_COMM_WORLD);
        }
    }
    
    
    
    /* FIND ALL BLOCKS IN PARALLEL (openmp only on master node) USING OPTIMISED CODE AND TIME EXECUTION */
    
    //First, free previous block database and reset block counter
    if(commRank == 0) {
        printf("Finding all blocks using parallelised optimised code\n");
        freeBD(blockDatabase,numBlocks);
        blockDatabase = malloc(sizeof(Block));
        numBlocks = 0;
        
        //Now, find all blocks using parallel optimised code
        startTime = omp_get_wtime();
        blockDatabase = findBlocksOptimisedParallel(blockDatabase, transposedData, keyDatabase, &numBlocks);
        execTime = omp_get_wtime() - startTime;
        timeForParallelOptimisedBlockGeneration = execTime;
        printf("Parallelised optimised block-generation complete\n\n");
    }
    

    
    /* FIND ALL COLLISIONS SEQUENTIALLY USING OPTIMISED CODE AND TIME EXECUTION */

    //First, free previous collision database and reset collision counter
    if(commRank == 0) {
        printf("Finding all collisions using sequential optimised code\n");
        int numCollisions = 0;
        
        //Now, find all collisions using optimised sequential code
        startTime = omp_get_wtime();
        collisions = findCollisionsOptimised(blockDatabase, numBlocks, &numCollisions);
        execTime = omp_get_wtime() - startTime;
        timeForSequentialOptimisedCollisionDetection = execTime;
        printf("Sequential optimised collision detection complete\n\n");
    }

    
    

    /* FIND ALL COLLISIONS IN PARALLEL (Openmp only on master node) USING OPTIMISED CODE AND TIME EXECUTION */

    //First, free previous collision database and reset collision counter
    if(commRank == 0) {
        printf("Finding all collisions using parallelised optimised collision detection\n");
        freeCollisionDB(collisions, numCollisions);
        numCollisions = 0;
        
        //Now, find all collisions using optimised parallel code
        startTime = omp_get_wtime();
        collisions = findCollisionsOptimisedParallel(blockDatabase, numBlocks, &numCollisions);
        execTime = omp_get_wtime() - startTime;
        timeForParallelOptimisedCollisionDetection = execTime;
        printf("Parallelised optimised collision detection complete\n\n");
    }
    
    
    /* FIND ALL COLLISIONS IN PARALLEL USING MPI AND TIME EXECUTION */
    MPI_Barrier(MPI_COMM_WORLD);
    if(commRank == 0) {
        freeCollisionDB(collisions, numCollisions);
    }
    numCollisions = 0;
    
    if(commRank == 0) {
        /* MASTER PROCESS */
        startTime = omp_get_wtime();
        collisions = malloc(sizeof(Collision));
        numCollisions = 0;
        
        //All slave processes will compute their collisions, and send them back
        for(int i=1; i<commSize; i++) {
            //First, they will send the number of collisions they found
            int numCollisionsFoundBySlave;
            MPI_Status status;
            MPI_Recv(&numCollisionsFoundBySlave, 1, MPI_INT, i, 123, MPI_COMM_WORLD, &status);
            
            //Now, they will send all collisions
            for(int j=0; j<numCollisionsFoundBySlave; j++) {
                //Reallocate memory for collision database
                collisions = (Collision *) realloc(collisions, sizeof(Collision) * (numCollisions + 1));
                //The slave will start by sending the number of blocks found in the collision
                int numBlocksFound;
                MPI_Recv(&numBlocksFound, 1, MPI_INT, i, 123, MPI_COMM_WORLD, &status);
                collisions[numCollisions].numBlocksInCollision = numBlocksFound;
                collisions[numCollisions].columns = (int *) malloc(sizeof(int) * numBlocksFound);
                collisions[numCollisions].blocks = (Block *) malloc(sizeof(Block) * numBlocksFound);
                //For each block, the slave will send the block's details
                for(int k=0; k<numBlocksFound; k++) {
                    int curColumn;
                    long long curSig;
                    int curRows[4];
                    MPI_Recv(&curColumn, 1, MPI_INT, i, 123, MPI_COMM_WORLD, &status);
                    MPI_Recv(&curSig, 1, MPI_LONG_LONG_INT, i, 123, MPI_COMM_WORLD, &status);
                    MPI_Recv(curRows, 4, MPI_INT, i, 123, MPI_COMM_WORLD, &status);
                    
                    //Store the current column number in the collision
                    collisions[numCollisions].columns[k] = curColumn;
                    //Store the current block details
                    collisions[numCollisions].blocks[k].signature = curSig;
                    collisions[numCollisions].blocks[k].column = curColumn;
                    collisions[numCollisions].blocks[k].rows[0] = curRows[0];
                    collisions[numCollisions].blocks[k].rows[1] = curRows[1];
                    collisions[numCollisions].blocks[k].rows[2] = curRows[2];
                    collisions[numCollisions].blocks[k].rows[3] = curRows[3];
                }
                numCollisions++;
            }
        }
        //Record execution time
        execTime = omp_get_wtime() - startTime;
        timeForMPIOptimisedCollisionDetection = execTime;
        printf("MPI Optimised collision detection complete\n\n");
    }
    
    else {
        /* SLAVE PROCESS */
        //First, find all the collisions I'm responsible for
        int numPartialCollisions;
        Collision *partialCollisionDatabase = findCollisionsOptimisedParallel(partialBlockDatabase, numBlocksPartial, &numPartialCollisions);
        
        //Now, send to master how many collisions I found
        MPI_Send(&numPartialCollisions, 1, MPI_INT, 0, 123, MPI_COMM_WORLD);
        //Now, send each collision's details
        for(int j=0; j<numPartialCollisions; j++) {
            //Send the number of blocks found in this collision
            MPI_Send(&partialCollisionDatabase[j].numBlocksInCollision, 1, MPI_INT, 0, 123, MPI_COMM_WORLD);
            //For each block, send its details
            for(int k=0; k<partialCollisionDatabase[j].numBlocksInCollision; k++) {
                //Send column number first
                MPI_Send(&partialCollisionDatabase[j].columns[k], 1, MPI_INT, 0, 123, MPI_COMM_WORLD);
                //Send the block's signature and rows
                MPI_Send(&partialCollisionDatabase[j].blocks[k].signature, 1, MPI_LONG_LONG_INT, 0, 123, MPI_COMM_WORLD);
                MPI_Send(partialCollisionDatabase[j].blocks[k].rows, 4, MPI_INT, 0, 123, MPI_COMM_WORLD);
            }
        }
        //Now, free all memory for partial databases
        free(partialBlockDatabase);
        freeCollisionDB(partialCollisionDatabase, numPartialCollisions);
    }
    
    
    
    /* PRINT ALL RESULTS TO RESULTS FILE (ONLY MASTER DOES THIS) */
    if(commRank == 0) {
        FILE *resultsOutput = fopen("results.txt", "w");
        if(resultsOutput == NULL) {
            fprintf(stderr, "Error. Can't create file \"results.txt\" to print results to\n");
            exit(EXIT_FAILURE);
        }
        fprintf(resultsOutput,"I/O took                                                   %10lf seconds\n\n", timeForIO);
        fprintf(resultsOutput,"Sequential optimised block generation took                 %10lf seconds\n", timeForSequentialOptimisedBlockGeneration);
        fprintf(resultsOutput,"Openmp optimised block generation took                     %10lf seconds\n", timeForParallelOptimisedBlockGeneration);
        fprintf(resultsOutput,"MPI optimised block generation took                        %10lf seconds\n", timeForMPIOptimisedBlockGeneration);
        fprintf(resultsOutput,"Speed-up factor on optimised block generation with omp     %10lf\n", (double) timeForSequentialOptimisedBlockGeneration / timeForParallelOptimisedBlockGeneration);
        fprintf(resultsOutput,"Speed-up factor on optimised block generation with MPI     %10lf\n\n", (double) timeForParallelOptimisedBlockGeneration / timeForMPIOptimisedBlockGeneration);
        fprintf(resultsOutput,"Sequential optimised collision detection took              %10lf seconds\n", timeForSequentialOptimisedCollisionDetection);
        fprintf(resultsOutput,"Parallel optimised collision detection took                %10lf seconds\n", timeForParallelOptimisedCollisionDetection);
        fprintf(resultsOutput,"MPI optimised collision detection took                     %10lf seconds\n", timeForMPIOptimisedCollisionDetection);
        fprintf(resultsOutput,"Speed-up factor on openmp optimised collision detection is %10lf\n", (double) timeForSequentialOptimisedCollisionDetection / timeForParallelOptimisedCollisionDetection);
        fprintf(resultsOutput,"Speed-up factor on MPI optimisied collision detection is   %10lf\n\n", (double) timeForParallelOptimisedCollisionDetection / timeForMPIOptimisedCollisionDetection);
        fprintf(resultsOutput,"%d Blocks, %d Collisions\n", numBlocks, numCollisions);
        fclose(resultsOutput);
        
        /* PRINT ALL COLLISIONS FOUND TO OUTPUT FILE */
        FILE *dataOutput = fopen("output.txt", "w");
        if(dataOutput == NULL) {
            fprintf(stderr, "Error. Can't create file \"output.txt\" to print located collisions to\n");
            exit(EXIT_FAILURE);
        }
        //Loop through all collisions
        for(int i=0; i<numCollisions; i++) {
            fprintf(dataOutput, "Collision:         %d\n", i+1);
            fprintf(dataOutput, "Signature:         %lld\n", collisions[i].blocks[0].signature);
            fprintf(dataOutput, "Number of blocks:  %d\n", collisions[i].numBlocksInCollision);
            fprintf(dataOutput, "Columns:           %d", collisions[i].columns[0]);
            //Print all column numbers
            for(int j=1; j<collisions[i].numBlocksInCollision; j++) {
                fprintf(dataOutput, ", %d", collisions[i].columns[j]);
            }
            fprintf(dataOutput, "\n");
            fprintf(dataOutput, "Rows:              %d, %d, %d, %d\n", collisions[i].blocks[0].rows[0], collisions[i].blocks[0].rows[1], collisions[i].blocks[0].rows[2], collisions[i].blocks[0].rows[3]);
            //Print newlines
            fprintf(dataOutput, "\n\n\n\n");
        }
        fclose(dataOutput);
        
        /* Finish program timer and print completion message to screen */
        double programTimeEnd = omp_get_wtime();
        double programTime = programTimeEnd - programTimeStart;
        printf("Program completed execution. Execution took %lf seconds\n", programTime);
        printf("All located collisions printed to \"output.txt\"\n");
        printf("Detailed analysis of execution time and parallelisation speed-up factors printed to \"results.txt\"\n");
    }


    /* FINAL MEMORY CLEANUP */
    
    //Free all dynamically allocated memory for key and matrix databases
    freeTransposedData(transposedData);
    free(keyDatabase);
    
    if(commRank == 0) {
        //Free dynamically allocated memory for block database
        free(blockDatabase);
        //Free dynamically allocated memory for collision database
        freeCollisionDB(collisions, numCollisions);
        //Free dynamically allocated memory for merged collision database
        //    freeMergedDB(merged, numMerged);
    }


    MPI_Finalize();
    /* EXIT PROGRAM */
	return EXIT_SUCCESS;
}
