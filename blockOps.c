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

// /*
//     generateBlocksSlide

//     array - all values in a column from the  data matrix
//     blockDB - 2d array of blocks being created filled by blocks
//     kd - Key database that represent rows
//     blockIndex - Current position in blockDB to fill
// */
// void generateBlocksSlide(pair* array, Block** blockDB, long long* kd, int* col, int* blockIndex){
//     // lower bound in array
//     int lower = 0;
//     // sort the entire array in ascending order O(nlgn)
//     qsort(array, ROWS, sizeof(pair), compareDoubles);
//     // Upper bound incrementing is growing the size of the window, default start at 1
//     // Go out of bounds with inclusion of boundary ROWS
//     for(int upper = 1; upper <= ROWS; ++upper){
//         // Check if window contains elements in the same neighbourhood.
//         // If not move the lower bound up till it does or until upper == lower
//         while(array[upper-1].value - array[lower].value > DIA){
//             ++lower;
//             // Take note that the lower bound i.e. new window instance
//         }
//         // Get all combinations within the window of size 4
//         for(int i = lower; i < upper-1; ++i){
//             for(int j = i+1; j < upper-1; ++j){
//                 for(int k = j+1; k < upper-1; ++k){
//                     //  Increase memory for the new block pointer to the database
//                     blockDB = (Block *) realloc(blockDB, (*blockIndex+1)*sizeof(Block));
//                     addBlock(blockDB, blockIndex, kd, col, &array[i].index, &array[j].index, &array[k].index, &array[upper-1].index);

//                     // Uncomment to print all rows / indexes being found
//                     printf("Found block at column %d on rows %d, %d, %d, %d\n", *col, array[i].index, array[j].index, array[k].index, array[upper-1].index);
//                     //  Increment number of blocks
//                     (*blockIndex)++;
//                 }
//             }
//         }
//     }
// }




//     generateBlocksBrute

//     array - all values in a column from the  data matrix
//     blockDB - 2d array of blocks being created filled by blocks
//     kd - Key database that represent rows
//     blockIndex - Current position in blockDB to fill

// void generateBlocksBrute(double* array, Block** blockDB, long long* kd, int* col, int* blockIndex){
//     //Loop through matrix rows
//     for(int row1=0; row1<ROWS; row1++) {
//         for(int row2=row1+1; row2<ROWS; row2++) {
//             //Ensure row1 and row2 are unique
//             if(row2 == row1) continue;
//             //Check if they're in the same neighbourhood
//             if(fabs(array[row1] - array[row2])>DIA) continue;
//             for(int row3=row2+1; row3<ROWS; row3++) {
//                 //Ensure rows 1, 2 and 3 are unique
//                 if(row3 == row2 || row3 == row1) continue;
//                 //Check they're in the same neighbourhood
//                 if(fabs(array[row3]-array[row2])>DIA || fabs(array[row3]-array[row1])>DIA) continue;
//                 for(int row4=row3+1; row4<ROWS; row4++) {
//                     //Ensure all rows are unique
//                     if(row4==row1 || row4==row2 || row4==row3) continue;
//                     //Check they're in the same neighbourhood
//                     if(fabs(array[row4]-array[row1])>DIA || fabs(array[row4] - array[row2])>DIA || fabs(array[row4]-array[row3])>DIA) continue;
//                     //  We have found a block, and must store it in the block database
//                     //  Increase memory for the new block pointer to the database
//                     blockDB = (Block **) realloc(blockDB, (*blockIndex+1)*sizeof(Block*));
//                     addBlock(blockDB, blockIndex, kd, col, &row1, &row2, &row3, &row4);
//                     // Increment blockIndex;
//                     (*blockIndex)++;
//                     //TEST OUTPUT
//                     printf("Found block at column %d on rows %d, %d, %d, %d\n", *col, row1, row2, row3, row4);
//                 }
//             }
//         }
//     }
// }

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
    for(int col=0; col< COLS - 1; col++) {
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
                        nextBlock++;
                        // Uncomment to print all rows / indexes being found
                        //printf("Found block at column %d on rows %lld, %lld, %lld, %lld\n", col, pairContainer[i].key, pairContainer[j].key, pairContainer[k].key, pairContainer[upper-1].key);
                        //  Increment number of blocks
                    }
                }
            }
        }











        // Use the double Array to generate blocks and fill the BlockDB
        //generateBlocksSlide(pairContainer, blockDB, kd, &col, numBlocks);
        // generateBlocksBrute(tempContainer, blockDB, kd, &col, numBlocks);
    }
    // Free the utility container
    // free(tempContainer);
    *numBlocks = nextBlock;
    free(pairContainer);

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