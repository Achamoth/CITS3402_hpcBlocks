/*
 CITS3402 Project 1 2016
 Name:           Ammar Abu Shamleh, Pradyumn Vij
 Student number: 21521274, 21469477
 Date:           September 2016
 */
#include "blocks.h"

/*
 findSig
 
 input key database, and four row numbers
 Finds sum of signatures at specified row numbers
 */
long long findSig(int r1, int r2, int r3, int r4, long long *kd) {
    long long sum = kd[r1]+kd[r2]+kd[r3]+kd[r4];
    return sum;
}

/*
 mergeBlockDatabases
 
 input partial block database and complete block database
 Copies all data (actually copies pointers) from partial block database to complete database
 */
Block *mergeBlockDatabases(Block *completeDB, Block *partialDB, int numBlockInPartial, int *numCopiedToComplete) {
    //Reallocate more memory for complete database
    completeDB = (Block *) realloc(completeDB, sizeof(Block) * (*numCopiedToComplete+numBlockInPartial));
    for(int i=0; i<numBlockInPartial; i++) {
        completeDB[i + *numCopiedToComplete] = partialDB[i];
    }
    *numCopiedToComplete += numBlockInPartial;
    return completeDB;
}

/*
 findBlocks
 
 input blockDatabase and matrixDatabase
 Finds all blocks in matrixDatabase and stores them in blockDatabase
 */
Block *findBlocks(Block *blockDB, double **mat, long long *kd, int *numBlocks) {
    int nextBlock = 0;
    //Loop through matrix columns (excluding last column)
    for(int col=0; col<COLS-1; col++) {
        //Loop through matrix rows in parallel
        omp_set_num_threads(NUM_THREADS);
        #pragma omp parallel
        {
            Block *partialBlockDB = (Block *) malloc(1 * sizeof(Block));
            int localNextBlock = 0; //Thread private counter for next index into block database
            int ID = omp_get_thread_num();
            int numThreads = omp_get_num_threads();
            for(int row1=ID; row1<ROWS; row1+=numThreads) {
                for(int row2=row1+1; row2<ROWS; row2++) {
                    //Check if they're in the same neighbourhood
                    if(fabs(mat[col][row1] - mat[col][row2])>DIA) continue;
                    for(int row3=row2+1; row3<ROWS; row3++) {
                        //Check they're in the same neighbourhood
                        if(fabs(mat[col][row3]-mat[col][row2])>DIA || fabs(mat[col][row3]-mat[col][row1])>DIA) continue;
                        for(int row4=row3+1; row4<ROWS; row4++) {
                            //Check they're in the same neighbourhood
                            if(fabs(mat[col][row4]-mat[col][row1])>DIA || fabs(mat[col][row4]-mat[col][row2])>DIA || fabs(mat[col][row4]-mat[col][row3])>DIA) continue;
                            //We have found a block, and must store it in the block database
                            partialBlockDB = (Block *) realloc(partialBlockDB, (localNextBlock+1)*sizeof(Block));
                            //partialBlockDB[localNextBlock] = (Block *) malloc(sizeof(Block));
                            partialBlockDB[localNextBlock].signature = findSig(row1, row2, row3, row4, kd);
                            partialBlockDB[localNextBlock].column = col;
                            localNextBlock++;
                            //TEST OUTPUT
                            printf("Thread %d: Found block at column %d on rows %d, %d, %d, %d\n", ID, col, row1, row2, row3, row4);
                        }
                    }
                }
            }
            //Merge partial block database with complete block database
            #pragma omp critical
            {
                blockDB = mergeBlockDatabases(blockDB, partialBlockDB, localNextBlock, &nextBlock);
            }
        }
    }
    
    //Find all blocks in last column using parallel tasks
    /*int col = 499;
    #pragma omp parallel
    {
        //Have a thread generate tasks to split block generation of last column up
        #pragma omp single
        {
            //int numThreads = omp_get_num_threads();
            //int increment = ceil((ROWS/numThreads));
            int increment = 300;
            int init = 0;
            int final = init+increment;
            while(final < ROWS) {
                //Generate a task for each chunk of the column that finds all the blocks in that chunk
                #pragma omp task firstprivate(init, final)
                {
                    //Private thread variables
                    Block **partialBlockDB = (Block **) malloc(1 * sizeof(Block *));
                    int localNextBlock = 0;
                    //Thread number
                    int ID = omp_get_thread_num();
                    for(int r1=init; r1<final; r1++) {
                        for(int r2=r1+1; r2<ROWS; r2++) {
                            //Check if they're in the same neighbourhood
                            if(fabs(mat[col][r1] - mat[col][r2]) > DIA) continue;
                            for(int r3 = r2+1; r3<ROWS; r3++) {
                                //Check if they're in the same neighbourhood
                                if(fabs(mat[col][r1]-mat[col][r3])>DIA || fabs(mat[col][r2]-mat[col][r3])>DIA) continue;
                                for(int r4=r3+1; r4<ROWS; r4++) {
                                    //Check they're in the same neighbourhood
                                    if(fabs(mat[col][r4]-mat[col][r1])>DIA || fabs(mat[col][r4]-mat[col][r2])>DIA || fabs(mat[col][r4]-mat[col][r4])>DIA) continue;
                                    //We have found a block, and must store it in the block database
                                    partialBlockDB = (Block **) realloc(partialBlockDB, (localNextBlock+1) * sizeof(Block *));
                                    partialBlockDB[localNextBlock] = (Block *) malloc(sizeof(Block));
                                    partialBlockDB[localNextBlock]->signature = findSig(r1, r2, r3, r4, kd);
                                    partialBlockDB[localNextBlock]->column = col;
                                    localNextBlock++;
                                    
                                    //TEST OUTPUT
                                    printf("Thread %d: Found block at column %d on rows %d, %d, %d, %d\n", ID, col, r1, r2, r3, r4);
                                }
                            }
                        }
                    }
                    //Add partial block database to complete block database
                    #pragma omp critical
                    {
                        blockDB = mergeBlockDatabases(blockDB, partialBlockDB, localNextBlock, &nextBlock);
                    }
                }
                //Calculate next chunk size and starting and ending index values
                final += increment;
                init += increment;
                if(final > ROWS) final = ROWS;
            }
        }
    }*/
    
    *numBlocks = nextBlock;
    printf("%d\n", *numBlocks);
    return blockDB;
}

/*
 findCollisions
 
 input blockDatabase
 Finds all collisions between generated blocks and return collision database; also store number of collisions found
 */
Collision *findCollisions(Block *blockDB, int numBlocks, int *numberCollisionsFound) {
    //Set up collision database
    int numCollisions = 0;
    Collision *collisions = (Collision *) malloc((numCollisions+1) * sizeof(Collision));
    //Allocate memory for first entry of collision database
    //collisions[0] = (Collision *) malloc(sizeof(Collision));
    //Record whether or not block has already been detected in a collision
    bool *collided = malloc(sizeof(bool) * numBlocks);
    for(int i=0; i<numBlocks; i++) collided[i] = false;
    
    //Loop over all blocks
    for(int i=0; i<numBlocks; i++) {
        Block curBlock = blockDB[i];
        long long curSig = curBlock.signature;
        int curCollisions = 0;
        //If current block has already been detected in a collision, skip it
        if(collided[i]) continue;
        //Start inner loop to compare to all other blocks
        for(int j=i+1; j<numBlocks; j++) {
            Block compBlock = blockDB[j];
            //Ensure that blocks are not in same column
            if(curBlock.column == compBlock.column) {
                continue;
            }
            //Check for collision
            if(compBlock.signature == curSig) {
                //Add block to current collision
                if(curCollisions == 0) {
                    collided[j] = true;
                    //Allocate memory for collision's column database, and store first column
                    collisions[numCollisions].columns = (int *) malloc(5 * sizeof(int));
                    collisions[numCollisions].numBlocksInCollision = 1;
                    collisions[numCollisions].columns[curCollisions] = curBlock.column;
                    curCollisions++;
                    //Allocate more memory for collision database
                    numCollisions++;
                    collisions = (Collision *) realloc(collisions, (numCollisions+1)*sizeof(Collision));
                }
                printf("%d %d: ", numCollisions-1, curCollisions);
                collisions[numCollisions-1].numBlocksInCollision += 1;
                collisions[numCollisions-1].columns[curCollisions] = compBlock.column;
                curCollisions++;
                //Allocate more memory for current collision column database
                collisions[numCollisions-1].columns = (int *) realloc(collisions[numCollisions-1].columns, ((curCollisions+1)*sizeof(int)));
                
                //TEST OUTPUT
                printf("Found collision at blocks %d and %d. Cols are: %d and %d. Sigs are: %lld and %lld \n", i, j, curBlock.column, collisions[numCollisions-1].columns[curCollisions-1], curSig, compBlock.signature);
                
            }
        }
    }
    free(collided);
    *numberCollisionsFound = numCollisions;
    return collisions;
}
