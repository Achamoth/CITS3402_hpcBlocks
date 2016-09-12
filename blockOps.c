/*
    CITS3402 Project 1 2016
    Name:           Ammar Abu Shamleh, Pradyumn Vij 
    Student number: 21521274, 21469477
    Date:           September 2016
*/
#include "blocks.h"

/*
    freeBD
 
    input blockDatabase: dynamically allocated database of Blocks
    Frees memory allocated to block database
 
*/
void freeBD(Block **bd, int numBlocks) {
    //Loop through all blocks
    for(int i=0; i<numBlocks; i++) {
        free(bd[i]);
    }
    free(bd);
}

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
    for(int col=0; col < allCols; col++) {
        //Loop through matrix allRows
        for(int row1=0; row1<allRows; row1++) {
            for(int row2=row1+1; row2<allRows; row2++) {
                //Ensure row1 and row2 are unique
                if(row2 == row1) continue;
                //Check if they're in the same neighbourhood
                if(fabs(mat[row1][col] - mat[row2][col])>DIA) continue;
                for(int row3=row2+1; row3<allRows; row3++) {
                    //Ensure allRows 1, 2 and 3 are unique
                    if(row3 == row2 || row3 == row1) continue;
                    //Check they're in the same neighbourhood
                    if(fabs(mat[row3][col]-mat[row2][col])>DIA || fabs(mat[row3][col]-mat[row1][col])>DIA) continue;
                    for(int row4=row3+1; row4<allRows; row4++) {
                        //Ensure all allRows are unique
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
                        printf("Found block at column %d on allRows %d, %d, %d, %d\n", col, row1, row2, row3, row4);
                    }
                }
            }
        }
    }
    *numBlocks = nextBlock;
    return blockDB;
}