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
    //Read in matrix, and record number of allRows and columns in allRows and allCols
    dataMatrix = readMatrix(DATA_FILE, dataMatrix);

	// Allocate memory for keys
	long long *keyDatabase = malloc(allRows*sizeof(long long));
	// Read in keys
	readKeys(KEY_FILE, keyDatabase);
    
    
    //Create pool of blocks
    Block **blockDatabase = malloc(1*sizeof(Block *));;
    
    //Find all blocks in matrix
    int numBlocks = 0;
    blockDatabase = findBlocks(blockDatabase, dataMatrix, keyDatabase, &numBlocks);
    
    //Free all dynamically allocated memory for key and matrix databases
    freeData(dataMatrix, keyDatabase);
    //Free dynamically allocated memory for block database
    freeBD(blockDatabase, numBlocks);
    //Exit program
	return EXIT_SUCCESS;
}
