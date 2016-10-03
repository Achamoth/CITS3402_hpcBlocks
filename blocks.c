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



    /* READ IN MATRIX DATA AND KEYS */

    //Allocate memory for matrix database
	double **dataMatrix = (double**)malloc(1*sizeof(double*));
    //Read in matrix, and record number of rows and columns in ROWS and COLS
    dataMatrix = readMatrix(DATA_FILE, dataMatrix);
    //Transpose data matrix
    double **transposedData = transposeMatrix(dataMatrix);
    //Free original matrix
    freeData(dataMatrix);

    // Allocate memory for keys
	long long *keyDatabase = malloc(ROWS*sizeof(long long));
	// Read in keys
	readKeys(KEY_FILE, keyDatabase);



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



    /* FIND ALL BLOCKS SEQUENTIALLY USING BRUTE FORCE CODE AND TIME EXECUTION */

    //Create pool of blocks
    Block *blockDatabase = malloc(1*sizeof(Block));;

    //Find all blocks in matrix
    int numBlocks = 0;
    startTime = omp_get_wtime();
    blockDatabase = findBlocks(blockDatabase, transposedData, keyDatabase, &numBlocks);
    execTime = omp_get_wtime() - startTime;
    timeForSequentialBlockGeneration = execTime;



    /* FIND ALL COLLISIONS SEQUENTIALLY USING BRUTE FORCE CODE AND TIME EXECUTION */

    //Find all collisions among blocks
    int numCollisions = 0;
    startTime = omp_get_wtime();
    Collision *collisions = findCollisions(blockDatabase, numBlocks, &numCollisions);
    execTime = omp_get_wtime() - startTime;
    timeForSequentialBruteForceCollisionDetection = execTime;



    /* FIND ALL BLOCKS IN PARALLEL USING BRUTE FORCE CODE AND TIME EXECUTION */

    //First free previous block database and reset block counter
    freeBD(blockDatabase, numBlocks);
    blockDatabase = malloc(sizeof(Block));
    numBlocks = 0;

    //Now, find all blocks using parallel brute force code
    startTime = omp_get_wtime();
    blockDatabase = findBlocksParallel(blockDatabase, transposedData, keyDatabase, &numBlocks);
    execTime = omp_get_wtime() - startTime;
    timeForParallelBlockGeneration = execTime;

    
    
    /* FIND ALL BLOCKS SEQUENTIALLY USING OPTIMISED CODE AND TIME EXECUTION */

	//First free previous block database and reset block counter
    freeBD(blockDatabase, numBlocks);
    blockDatabase = malloc(sizeof(Block));
    numBlocks = 0;

	//Now, find all blocks using sequential optimised code
    startTime = omp_get_wtime();
    blockDatabase = findBlocksOptimised(blockDatabase, transposedData, keyDatabase, &numBlocks);
    execTime = omp_get_wtime() - startTime;
    timeForSequentialOptimisedBlockGeneration = execTime;
    

    
    /* FIND ALL BLOCKS IN PARALLEL USING OPTIMISED CODE AND TIME EXECUTION */
    
    //First, free previous block database and reset block counter
    freeBD(blockDatabase,numBlocks);
    blockDatabase = malloc(sizeof(Block));
    numBlocks = 0;
    
    //Now, find all blocks using parallel optimised code
    startTime = omp_get_wtime();
    blockDatabase = findBlocksOptimisedParallel(blockDatabase, transposedData, keyDatabase, &numBlocks);
    execTime = omp_get_wtime() - startTime;
    timeForParallelOptimisedBlockGeneration = execTime;
    

    /* FIND ALL COLLISIONS IN PARALLEL USING BRUTE FORCE CODE AND TIME EXECUTION */

    //First, free previous collision database and reset collision counter
    freeCollisionDB(collisions, numCollisions);
    numCollisions = 0;

    //Now, find all collisions using parallel brute force code
    startTime = omp_get_wtime();
    collisions = findCollisionsParallel(blockDatabase, numBlocks, &numCollisions);
    execTime = omp_get_wtime() - startTime;
    timeForParallelBruteForceCollisionDetection = execTime;



    /* FIND ALL COLLISIONS SEQUENTIALLY USING OPTIMISED CODE AND TIME EXECUTION */

    //First, free previous collision database and reset collision counter
    freeCollisionDB(collisions, numCollisions);
    numCollisions = 0;

    //Now, find all collisions using optimised sequential code
    startTime = omp_get_wtime();
    collisions = findCollisionsOptimised(blockDatabase, numBlocks, &numCollisions);
    execTime = omp_get_wtime() - startTime;
    timeForSequentialOptimisedCollisionDetection = execTime;



    /* FIND ALL COLLISIONS IN PARALLEL USING OPTIMISED CODE AND TIME EXECUTION */

    //First, free previous collision database and reset collision counter
    freeCollisionDB(collisions, numCollisions);
    numCollisions = 0;

    //Now, find all collisions using optimised parallel code
    startTime = omp_get_wtime();
    collisions = findCollisionsOptimisedParallel(blockDatabase, numBlocks, &numCollisions);
    execTime = omp_get_wtime() - startTime;
    timeForParallelOptimisedCollisionDetection = execTime;



    /* PRINT ALL RESULTS */
    printf("\n\n\n\n");
    printf("Sequential brute-force block generation took         %10lf seconds\n", timeForSequentialBlockGeneration);
    printf("Sequential brute-force collision detection took      %10lf seconds\n", timeForSequentialBruteForceCollisionDetection);
    printf("Parallel brute-force block generation took           %10lf seconds\n", timeForParallelBlockGeneration);
	printf("Sequential optimised block generation took           %10lf seconds\n", timeForSequentialOptimisedBlockGeneration);
    printf("Parallel optimised block generation took             %10lf seconds\n", timeForParallelOptimisedBlockGeneration);
    printf("Parallel brute-force collision detection took        %10lf seconds\n", timeForParallelBruteForceCollisionDetection);
    printf("Sequential optimised collision detection took        %10lf seconds\n", timeForSequentialOptimisedCollisionDetection);
    printf("Parallel optimised collision detection took          %10lf seconds\n", timeForParallelOptimisedCollisionDetection);
    printf("\n\n\n\n");

//    //QUICK TEST
//    printBlock(collisions[0].blocks[0]);
//    printBlock(collisions[0].blocks[1]);
    
    /* FINAL MEMORY CLEANUP */

    //Free all dynamically allocated memory for key and matrix databases
    freeTransposedData(transposedData);
    free(keyDatabase);
    //Free dynamically allocated memory for block database
    free(blockDatabase);
    //Free dynamically allocated memory for collision database
    freeCollisionDB(collisions, numCollisions);



    /* EXIT PROGRAM */
	return EXIT_SUCCESS;
}
