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
    findSum
 
    inpput matrixDatabase and four row numbers, and a column number
    Finds sum of four elements at specified row and column number
*/
double findSum(int r1, int r2, int r3, int r4, int col, double **mat) {
    double sum = mat[r1][col] + mat[r2][col] + mat[r3][col] + mat[r4][col];
    return sum;
}

/*
    findBlocks
 
    input blockDatabase and matrixDatabase
    Finds all blocks in matrixDatabase and stores them in blockDatabase
*/
Block **findBlocks(Block **blockDB, double **mat, long long *kd, int *numBlocks) {
    int nextBlock = 0;
    //Loop through matrix columns
    for(int col=0; col<COLS; col++) {
        //Loop through matrix rows
        for(int row1=0; row1<ROWS; row1++) {
            for(int row2=row1+1; row2<ROWS; row2++) {
                //Ensure row1 and row2 are unique
                if(row2 == row1) continue;
                //Check if they're in the same neighbourhood
                if(fabs(mat[row1][col] - mat[row2][col])>DIA) continue;
                for(int row3=row2+1; row3<ROWS; row3++) {
                    //Ensure rows 1, 2 and 3 are unique
                    if(row3 == row2 || row3 == row1) continue;
                    //Check they're in the same neighbourhood
                    if(fabs(mat[row3][col]-mat[row2][col])>DIA || fabs(mat[row3][col]-mat[row1][col])>DIA) continue;
                    for(int row4=row3+1; row4<ROWS; row4++) {
                        //Ensure all rows are unique
                        if(row4==row1 || row4==row2 || row4==row3) continue;
                        //Check they're in the same neighbourhood
                        if(fabs(mat[row4][col]-mat[row1][col])>DIA || fabs(mat[row4][col]-mat[row2][col])>DIA || fabs(mat[row4][col]-mat[row3][col])>DIA) continue;
                        //We have found a block, and must store it in the block database
                        blockDB = (Block **) realloc(blockDB, (nextBlock+1)*sizeof(Block*));
                        blockDB[nextBlock] = (Block *) malloc(sizeof(Block));
                        blockDB[nextBlock]->signature = findSig(row1, row2, row3, row4, kd);
                        blockDB[nextBlock]->sumOfElements = findSum(row1, row2, row3, row4, col, mat);
                        blockDB[nextBlock]->column = col;
                        nextBlock++;
                        //TEST OUTPUT
                        printf("Found block at column %d on rows %d, %d, %d, %d\n", col, row1, row2, row3, row4);
                    }
                }
            }
        }
    }
    *numBlocks = nextBlock;
    return blockDB;
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
