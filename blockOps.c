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

long long findSig(int* r1, int* r2, int* r3, int* r4, long long* kd) {
    long long sum = kd[*r1]+kd[*r2]+kd[*r3]+kd[*r4];
    return sum;
}

void addBlock(Block** blockDB, int* blockIndex, long long* kd, int* col, int* row1, int* row2, int* row3, int* row4){
    // Add memory for block struct
    blockDB[*blockIndex] = (Block *) malloc(sizeof(Block));
    blockDB[*blockIndex]->signature = findSig(row1, row2, row3, row4, kd);
    blockDB[*blockIndex]->column = *col;
}


/*

    compareDoubles

    Comparator function for qSort, to sort doubles in ascending order

*/
int compareDoubles(const void* a, const void* b){
    return *(double*)a - *(double*)b;
}

/*
    generateBlocksSlide

    array - all values in a column from the  data matrix
    blockDB - 2d array of blocks being created filled by blocks
    kd - Key database that represent rows
    blockIndex - Current position in blockDB to fill
*/
void generateBlocksSlide(double* array, Block** blockDB, long long* kd, int* col, int* blockIndex){
    // lower bound in array
    int lower = 0;
    // sort the entire array in ascending order O(nlgn)
    qsort(array, ROWS, sizeof(double), compareDoubles);
    bool lowerChanged = false;
    // Upper bound incrementing is growing the size of the window, default start at 1
    for(int upper = 1; upper <= ROWS; ++upper){
        // Check if window contains elements in the same neighbourhood.
        // If not move the lower bound up till it does or until upper == lower
        while(array[upper-1] - array[lower] > DIA){
            ++lower;
            // Take note that the lower bound i.e. new window instance
            lowerChanged = true;
        }
        // If the lower bound has jut been changed
        // OR
        // We have reached the end of the array
        if(lowerChanged || upper == ROWS){
            // Get all combinations within the window of size 4
            for(int i = lower+1; i < upper; ++i){
                for(int j = i+1; j < upper; ++j){
                    for(int k = j+1; k < upper; ++k){
                        //  Increase memory for the new block pointer to the database
                        blockDB = (Block **) realloc(blockDB, (*blockIndex+1)*sizeof(Block*));
                        addBlock(blockDB, blockIndex, kd, col, &lower, &i, &j, &k);
                        //  Increment number of blocks
                        blockIndex++;
                    }
                }
            }
            // Reset flag for lower bounds change
            lowerChanged = false;
        }
    }    
}



/*

    generateBlocksBrute

    array - all values in a column from the  data matrix
    blockDB - 2d array of blocks being created filled by blocks
    kd - Key database that represent rows
    blockIndex - Current position in blockDB to fill
*/
void generateBlocksBrute(double* array, Block** blockDB, long long* kd, int* col, int* blockIndex){
    //Loop through matrix rows
    for(int row1=0; row1<ROWS; row1++) {
        for(int row2=row1+1; row2<ROWS; row2++) {
            //Ensure row1 and row2 are unique
            if(row2 == row1) continue;
            //Check if they're in the same neighbourhood
            if(fabs(array[row1] - array[row2])>DIA) continue;
            for(int row3=row2+1; row3<ROWS; row3++) {
                //Ensure rows 1, 2 and 3 are unique
                if(row3 == row2 || row3 == row1) continue;
                //Check they're in the same neighbourhood
                if(fabs(array[row3]-array[row2])>DIA || fabs(array[row3]-array[row1])>DIA) continue;
                for(int row4=row3+1; row4<ROWS; row4++) {
                    //Ensure all rows are unique
                    if(row4==row1 || row4==row2 || row4==row3) continue;
                    //Check they're in the same neighbourhood
                    if(fabs(array[row4]-array[row1])>DIA || fabs(array[row4] - array[row2])>DIA || fabs(array[row4]-array[row3])>DIA) continue;
                    //  We have found a block, and must store it in the block database
                    //  Increase memory for the new block pointer to the database
                    blockDB = (Block **) realloc(blockDB, (*blockIndex+1)*sizeof(Block*));
                    addBlock(blockDB, blockIndex, kd, col, &row1, &row2, &row3, &row4);
                    // Increment blockIndex;
                    blockIndex++;
                    //TEST OUTPUT
                    printf("Found block at column %d on rows %d, %d, %d, %d\n", *col, row1, row2, row3, row4);
                }
            }
        }
    }
}

/*
    findBlocks

    input blockDatabase and matrixDatabase
    Finds all blocks in matrixDatabase and stores them in blockDatabase
*/

void findBlocks(Block **blockDB, double **mat, long long *kd, int *numBlocks) {
    int nextBlock = 0;
    //Loop through matrix columns
    double* tempContainer = malloc(ROWS * sizeof(double));
    for(int col=0; col<COLS; col++) {
        // Create an array of doubles
        for(int row = 0; row < ROWS; ++row){
            tempContainer[row] = mat[row][col];
        }
        // Use the double Array to generate blocks and fill the BlockDB
        generateBlocksBrute(tempContainer, blockDB, kd, &col, numBlocks);
    }
    // Free the utility container
    free(tempContainer);
    *numBlocks = nextBlock;
}



/*
    findCollisions

    input blockDatabase
    Finds all collisions between generated blocks and return collision database; also store number of collisions found
*/

Collision **findCollisions(Block **blockDB, int numBlocks, int *numberCollisionsFound) {
    //Set up collision database
    int numCollisions = 0;
    Collision **collisions = (Collision **) malloc((numCollisions+1) * sizeof(Collision));
    //Allocate memory for first entry of collision database
    collisions[0] = (Collision *) malloc(sizeof(Collision));
    //Record whether or not block has already been detected in a collision
    bool *collided = malloc(sizeof(bool) * numBlocks);
    for(int i=0; i<numBlocks; i++) collided[i] = false;
    //Loop over all blocks
    for(int i=0; i<numBlocks; i++) {
        Block *curBlock = blockDB[i];
        long long curSig = curBlock->signature;
        int curCollisions = 0;
        //If current block has already been detected in a collision, skip it
        if(collided[i]) continue;
        //Start inner loop to compare to all other blocks
        for(int j=i+1; j<numBlocks; j++) {
            Block *compBlock = blockDB[j];
            //Ensure that blocks are not in same column
            if(curBlock->column == compBlock->column) {
                continue;
            }
            //Check for collision
            if(compBlock->signature == curSig) {
                //Add block to current collision
                if(curCollisions == 0) {
                    collided[j] = true;
                    collisions[numCollisions] = (Collision *) malloc(sizeof(Collision));
                    collisions[numCollisions]->collidingBlocks = (Block **) malloc(sizeof(Block));
                    collisions[numCollisions]->numBlocksInCollision = 0;
                    //Allocate more memory for collision database
                    numCollisions++;
                    collisions = (Collision **) realloc(collisions, (numCollisions+1)*sizeof(Collision *));
                }
                printf("%d %d: ", numCollisions-1, curCollisions);
                collisions[numCollisions-1]->collidingBlocks[curCollisions] = compBlock;
                collisions[numCollisions-1]->numBlocksInCollision += 1;
                curCollisions++;
                //Allocate more memory for current collision blocks
                collisions[numCollisions-1]->collidingBlocks = (Block **) realloc(collisions[numCollisions-1]->collidingBlocks, ((curCollisions+1)*sizeof(Block *)));           
                //TEST OUTPUT
                printf("Found collision at blocks %d and %d. Sigs are: %lld and %lld \n", i, j, curSig, collisions[numCollisions-1]->collidingBlocks[curCollisions-1]->signature);
            }
        }
    }
    free(collided);
    *numberCollisionsFound = numCollisions;
    return collisions;
}