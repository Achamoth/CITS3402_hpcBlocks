/*
    CITS3402 Project 1 2016
    Name:           Ammar Abu Shamleh, Pradyumn Vij
    Student number: 21521274, 21469477
    Date:           September 2016
*/

#include "blocks.h"

/*
    compareDoubles

    Comparator function for qSort, to sort doubles in ascending order
*/
int compareDoubles(const void* a, const void* b){
    pair resA = *(pair*)a;
    pair resB = *(pair*)b;
    if(resA.value < resB.value) return -1;
    if(resA.value > resB.value) return 1;
    return 0;
}

/*
    findBlocks

    input blockDatabase and matrixDatabase
    Finds all blocks in matrixDatabase and stores them in blockDatabase
*/
Block *findBlocks(Block *blockDB, double **mat, long long *kd, int *numBlocks) {
    int nextBlock = 0;
    //Loop through matrix columns
    // double* tempContainer = malloc(ROWS * sizeof(double));
    pair* pairContainer = malloc(ROWS * sizeof(pair));
    for(int col=0; col< COLS; col++) {
        // Create an array of doubles
        for(int row = 0; row < ROWS; ++row){
            //tempContainer[row] = mat[row][col];
            pairContainer[row].value = mat[col][row];
            pairContainer[row].key = kd[row];
        }

        // Use sliding technique to fill BlockDB
        // lower bound in array
        int lower = 0;
        // sort the entire array in ascending order O(nlgn)
        qsort(pairContainer, ROWS, sizeof(pair), compareDoubles);
        // Upper bound incrementing is growing the size of the window, default start at 1
        // Go out of bounds with inclusion of boundary ROWS
        for(int upper = 1; upper <= ROWS; ++upper){
            // Check if window contains elements in the same neighbourhood.
            // If not move the lower bound up till it does or until upper == lower
            while(pairContainer[upper-1].value - pairContainer[lower].value > DIA){
                ++lower;
                // Take note that the lower bound i.e. new window instance
            }
            // Get all combinations within the window of size 4
            for(int i = lower; i < upper-1; ++i){
                for(int j = i+1; j < upper-1; ++j){
                    for(int k = j+1; k < upper-1; ++k){
                        //  Increase memory for the new block pointer to the database
                        blockDB = (Block *) realloc(blockDB, (nextBlock+1)*sizeof(Block));
                        blockDB[nextBlock].signature = pairContainer[i].key + pairContainer[j].key + pairContainer[k].key + pairContainer[upper-1].key;
                        blockDB[nextBlock].column = col;
                        //  Increment number of blocks
                        nextBlock++;
                        // Uncomment to print all rows / indexes being found
                        //printf("Found block at column %d on rows %lld, %lld, %lld, %lld\n", col, pairContainer[i].key, pairContainer[j].key, pairContainer[k].key, pairContainer[upper-1].key);
                    }
                }
            }
        }
    }
    // Free the utility container
    *numBlocks = nextBlock;
    free(pairContainer);
    printf("%d\n", *numBlocks);
    return blockDB;
}




/*
 cmpfunc

 Input; two variables of same datatype
 Performs comparison for sorting function. Returns positive if first value is larger; negative if second value is larger
 */
int cmpfunc(const void *a, const void *b) {
    //Find value of block a's signature
    long sigA = 0;
    Block ba = *(Block *) a;
    sigA = ba.signature;

    //Find value of block b's signature
    long sigB = 0;
    Block bb = *(Block *) b;
    sigB = bb.signature;

    //Return comparison
    return sigA - sigB;
}

/*
 findCollisionsOptimised

 input blockDatabase
 Finds all collisions between generated blocks and return collision database; also store number of collisions found. Using sorting method instead of brute force
 */
Collision *findCollisions(Block *blockDB, int numBlocks, int *numberCollisionsFound) {
    //Sort block database
    qsort(blockDB, numBlocks, sizeof(Block), cmpfunc);

    //Set up collision database
    Collision *collisions = (Collision *) malloc(1 * sizeof(Collision));

    //Linearly loop through block database, storing collisions as they're found
    int numCollisions = 0;
    long previousSig = blockDB[0].signature;
    int curBlocksInCollision = 0;
    for(int i=1; i<numBlocks; i++) {
        if(blockDB[i].signature == previousSig) {
            //Collision detected
            if(curBlocksInCollision == 0) {
                //New collision. Reallocate more memory for collision database
                collisions = (Collision *) realloc(collisions, sizeof(Collision) * (numCollisions + 1));
                //Increment counter
                numCollisions++;
                //Store first two blocks of current collision
                curBlocksInCollision = 2;
                collisions[numCollisions-1].numBlocksInCollision = 2;
                collisions[numCollisions-1].columns = (int *) malloc(sizeof(int) * 5);
                collisions[numCollisions-1].columns[0] = blockDB[i-1].column;
                collisions[numCollisions-1].columns[1] = blockDB[i].column;
                //printf("%d: Found collision on signature %ld with %d blocks in it\n", numCollisions-1, previousSig, curBlocksInCollision);
            }
            else {
                //Increment counter
                curBlocksInCollision++;
                collisions[numCollisions-1].numBlocksInCollision = curBlocksInCollision;
                //Store current block's column number, and rellocate collisions column array
                collisions[numCollisions-1].columns = (int *) realloc(collisions[numCollisions-1].columns, sizeof(int) * curBlocksInCollision);
                collisions[numCollisions-1].columns[curBlocksInCollision-1] = blockDB[i].column;
                //printf("%d: Found collision on signature %ld with %d blocks in it\n", numCollisions-1, previousSig, curBlocksInCollision);
            }
        }
        else {
            //No collision between previous block and current block. Record new signature
            previousSig = blockDB[i].signature;
            curBlocksInCollision = 0;
        }
    }
    //Return collision database
    return collisions;
}
