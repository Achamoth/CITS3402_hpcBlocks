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
    
    //Create pool of blocks
    Block **blockDatabase = malloc(1*sizeof(Block *));;
    
    //Find all blocks in matrix
    int numBlocks = 0;
    blockDatabase = findBlocks(blockDatabase, transposedData, keyDatabase, &numBlocks);
    
    //Find all collisions among blocks
    int numCollisions = 0;
    Collision **collisions = findCollisions(blockDatabase, numBlocks, &numCollisions);
    
    //Free all dynamically allocated memory for key and matrix databases
    freeTransposedData(transposedData);
    free(keyDatabase);
    //Free dynamically allocated memory for block database
    freeBD(blockDatabase, numBlocks);
    //Free dynamically allocated memory for collision database
    freeCollisionDB(collisions, numCollisions);
    //Exit program
	return EXIT_SUCCESS;
}
