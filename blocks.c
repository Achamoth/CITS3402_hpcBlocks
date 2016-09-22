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

	// Open KEY_FILE
	FILE* keyFile= openFile(KEY_FILE);
	// Count keys
	countKeys(keyFile);
	// Allocate memory for key array
	long long* keyDatabase = malloc(ROWS * sizeof(long long));
	// Read in keys
	readKeys(keyFile, keyDatabase);
	// Close file
	closeFile(keyFile);

	// Open DATA_FILE
	FILE* dataFile = openFile(DATA_FILE);
	// Count rows and columns in data file
	countRowsCols(dataFile);
	// Allocate contiguous memory for 2d array
	double** dataMatrix = malloc(ROWS*sizeof(double*));
	double* temp = malloc(ROWS * COLS * sizeof(double));
	for(int i = 0; i < ROWS; ++i){
		dataMatrix[i] = temp + (i * COLS);
	}
	// Read filestream into 2d array
	readData(dataFile, dataMatrix);
	// Close DATA filestream
	closeFile(dataFile);

    //Create pool of blocks
    Block **blockDatabase = malloc(1*sizeof(Block *));;

    //Find all blocks in matrix
    int numBlocks = 0;
    findBlocks(blockDatabase, dataMatrix, keyDatabase, &numBlocks);

    //Find all collisions among blocks
    int numCollisions = 0;
    Collision **collisions = findCollisions(blockDatabase, numBlocks, &numCollisions);

    //Free all dynamically allocated memory for key and matrix databases
	free(keyDatabase);
	free(dataMatrix);
	free(temp);
    //freeData(dataMatrix, keyDatabase);
    //Free dynamically allocated memory for block database
    freeBD(blockDatabase, numBlocks);
    //Free dynamically allocated memory for collision database
    freeCollisionDB(collisions, numCollisions);
    //Exit program
	return EXIT_SUCCESS;
}
