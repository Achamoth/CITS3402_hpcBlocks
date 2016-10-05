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
	// Program name without /, cast to constant for file
	programName = (const char*) strrchr(argv[0], '/') + 1;
    
    //Time overall execution of program
    double programTimeStart = omp_get_wtime();
    
    /* READ IN MATRIX DATA AND KEYS */

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



    /* THESE VARIABLES ARE FOR TIMING EXECUTION AND RETURNING EXPERIMENTAL RESULTS */
    double timeForSequentialBlockGeneration = 0;
    double timeForParallelBlockGeneration = 0;
    double timeForSequentialBruteForceCollisionDetection = 0;
	double timeForSequentialOptimisedBlockGeneration = 0;
    double timeForParallelOptimisedBlockGeneration = 0;
    double timeForParallelBruteForceCollisionDetection = 0;
    double timeForSequentialOptimisedCollisionDetection = 0;
    double timeForParallelOptimisedCollisionDetection = 0;
    double startTime;
    double execTime;


    //Set number of threads for parallel regions
    omp_set_num_threads(NUM_THREADS);

    /* FIND ALL BLOCKS SEQUENTIALLY USING BRUTE FORCE CODE AND TIME EXECUTION */

    //Create pool of blocks
    Block *blockDatabase = malloc(1*sizeof(Block));;

    //Find all blocks in matrix
    printf("Finding all blocks using sequential brute-force code\n");
    int numBlocks = 0;
    startTime = omp_get_wtime();
    blockDatabase = findBlocks(blockDatabase, transposedData, keyDatabase, &numBlocks);
    execTime = omp_get_wtime() - startTime;
    timeForSequentialBlockGeneration = execTime;
    printf("Sequential brute-force block generation completed\n\n");



    /* FIND ALL COLLISIONS SEQUENTIALLY USING BRUTE FORCE CODE AND TIME EXECUTION */

    //Find all collisions among blocks
    printf("Finding all collisions using sequential brute-force code\n");
    int numCollisions = 0;
    startTime = omp_get_wtime();
    Collision *collisions = findCollisions(blockDatabase, numBlocks, &numCollisions);
    execTime = omp_get_wtime() - startTime;
    timeForSequentialBruteForceCollisionDetection = execTime;
    printf("Sequential brute-force collision detection finished\n\n");



    /* FIND ALL BLOCKS IN PARALLEL USING BRUTE FORCE CODE AND TIME EXECUTION */

    //First free previous block database and reset block counter
    printf("Finding all blocks using parallelised brute-force code\n");
    freeBD(blockDatabase, numBlocks);
    blockDatabase = malloc(sizeof(Block));
    numBlocks = 0;

    //Now, find all blocks using parallel brute force code
    startTime = omp_get_wtime();
    blockDatabase = findBlocksParallel(blockDatabase, transposedData, keyDatabase, &numBlocks);
    execTime = omp_get_wtime() - startTime;
    timeForParallelBlockGeneration = execTime;
    printf("Parallelised brute-force block generation complete\n\n");

    
    
    /* FIND ALL BLOCKS SEQUENTIALLY USING OPTIMISED CODE AND TIME EXECUTION */

	//First free previous block database and reset block counter
    printf("Finding all blocks using sequential optimised code\n");
    freeBD(blockDatabase, numBlocks);
    blockDatabase = malloc(sizeof(Block));
    numBlocks = 0;

	//Now, find all blocks using sequential optimised code
    startTime = omp_get_wtime();
    blockDatabase = findBlocksOptimised(blockDatabase, transposedData, keyDatabase, &numBlocks);
    execTime = omp_get_wtime() - startTime;
    timeForSequentialOptimisedBlockGeneration = execTime;
    printf("Sequential optimised block generation complete\n\n");
    

    
    /* FIND ALL BLOCKS IN PARALLEL USING OPTIMISED CODE AND TIME EXECUTION */
    
    //First, free previous block database and reset block counter
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
    

    /* FIND ALL COLLISIONS IN PARALLEL USING BRUTE FORCE CODE AND TIME EXECUTION */

    //First, free previous collision database and reset collision counter
    printf("Finding all collisions using parallelised brute-force code\n");
    freeCollisionDB(collisions, numCollisions);
    numCollisions = 0;

    //Now, find all collisions using parallel brute force code
    startTime = omp_get_wtime();
    collisions = findCollisionsParallel(blockDatabase, numBlocks, &numCollisions);
    execTime = omp_get_wtime() - startTime;
    timeForParallelBruteForceCollisionDetection = execTime;
    printf("Parallelised brute-force collision detection complete\n\n");



    /* FIND ALL COLLISIONS SEQUENTIALLY USING OPTIMISED CODE AND TIME EXECUTION */

    //First, free previous collision database and reset collision counter
    printf("Finding all collisions using sequential optimised code\n");
    freeCollisionDB(collisions, numCollisions);
    numCollisions = 0;

    //Now, find all collisions using optimised sequential code
    startTime = omp_get_wtime();
    collisions = findCollisionsOptimised(blockDatabase, numBlocks, &numCollisions);
    execTime = omp_get_wtime() - startTime;
    timeForSequentialOptimisedCollisionDetection = execTime;
    printf("Sequential optimised collision detection complete\n\n");



    /* FIND ALL COLLISIONS IN PARALLEL USING OPTIMISED CODE AND TIME EXECUTION */

    //First, free previous collision database and reset collision counter
    printf("Finding all collisions using parallelised optimised collision detection\n");
    freeCollisionDB(collisions, numCollisions);
    numCollisions = 0;

    //Now, find all collisions using optimised parallel code
    startTime = omp_get_wtime();
    collisions = findCollisionsOptimisedParallel(blockDatabase, numBlocks, &numCollisions);
    execTime = omp_get_wtime() - startTime;
    timeForParallelOptimisedCollisionDetection = execTime;
    printf("Parallelised optimised collision detection complete\n\n");
    
    
    /* FIND ALL MERGED COLLISIONS */
    
    //Find all merged collisions
//    int numMerged = 0;
//    MergedCollision *merged = mergeCollisions(collisions, numCollisions, &numMerged);
    
    

    
    /* PRINT ALL RESULTS TO RESULTS FILE */
    FILE *resultsOutput = fopen("results.txt", "w");
    if(resultsOutput == NULL) {
        fprintf(stderr, "Error. Can't create file \"results.txt\" to print results to\n");
        exit(EXIT_FAILURE);
    }
    fprintf(resultsOutput,"I/O took                                                %10lf seconds\n\n", timeForIO);
    fprintf(resultsOutput,"Sequential brute-force block generation took            %10lf seconds\n", timeForSequentialBlockGeneration);
    fprintf(resultsOutput,"Parallel brute-force block generation took              %10lf seconds\n", timeForParallelBlockGeneration);
    fprintf(resultsOutput,"Speed-up factor on brute-force block generation is      %10lf\n\n", (double) timeForSequentialBlockGeneration/timeForParallelBlockGeneration);
    fprintf(resultsOutput,"Sequential brute-force collision detection took         %10lf seconds\n", timeForSequentialBruteForceCollisionDetection);
    fprintf(resultsOutput,"Parallel brute-force collision detection took           %10lf seconds\n", timeForParallelBruteForceCollisionDetection);
    fprintf(resultsOutput,"Speed-up factor on brute-force collision detection is   %10lf\n\n", (double) timeForSequentialBruteForceCollisionDetection/timeForParallelBruteForceCollisionDetection);
	fprintf(resultsOutput,"Sequential optimised block generation took              %10lf seconds\n", timeForSequentialOptimisedBlockGeneration);
    fprintf(resultsOutput,"Parallel optimised block generation took                %10lf seconds\n", timeForParallelOptimisedBlockGeneration);
    fprintf(resultsOutput,"Speed-up factor on optimised block generation           %10lf\n\n", (double) timeForSequentialOptimisedBlockGeneration / timeForParallelOptimisedBlockGeneration);
    fprintf(resultsOutput,"Sequential optimised collision detection took           %10lf seconds\n", timeForSequentialOptimisedCollisionDetection);
    fprintf(resultsOutput,"Parallel optimised collision detection took             %10lf seconds\n", timeForParallelOptimisedCollisionDetection);
    fprintf(resultsOutput,"Speed-up factor on optimised collision detection is     %10lf\n\n", (double) timeForSequentialOptimisedCollisionDetection / timeForParallelOptimisedCollisionDetection);
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
    
    
    /* FINAL MEMORY CLEANUP */

    //Free all dynamically allocated memory for key and matrix databases
    freeTransposedData(transposedData);
    free(keyDatabase);
    //Free dynamically allocated memory for block database
    free(blockDatabase);
    //Free dynamically allocated memory for collision database
    freeCollisionDB(collisions, numCollisions);
    //Free dynamically allocated memory for merged collision database
//    freeMergedDB(merged, numMerged);



    /* EXIT PROGRAM */
	return EXIT_SUCCESS;
}
