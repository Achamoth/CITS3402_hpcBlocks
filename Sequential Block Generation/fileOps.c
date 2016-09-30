/*
	CITS3402 Project 1 2016
	Name:			Ammar Abu Shamleh, Pradyumn Vij 
	Student number: 21521274, 21469477
    Date:           September 2016
*/
#include "blocks.h"

int ROWS;
int COLS;
    
/*
	readData

	input char pointer to file containing comma separated floating point numbers
	Reads data from data.txt into 2d array of doubles ROWS x COLS in size
*/
void readData(char *filename, double **dataMatrix){
	FILE *dataFile = fopen(filename, "r");

	//	Null file pointer check
	if(dataFile == NULL){
		fprintf(stderr, "%s - The file \"%s\" could not be opened\n", 
			programName, filename);

		exit(EXIT_FAILURE);
	}

	for(int r = 0; r < ROWS; ++r){
		for(int c = 0; c < COLS; ++c){
			//  Use fscanf to skip the comma characters
			fscanf(dataFile, "%lf%*c", &dataMatrix[r][c]);
            //TEST OUTPUT
			//printf("%lf found at index %d : %d\n", dataMatrix[r][c], r, c);
		}
	}

	//	Close file pointer if open
	if(dataFile != NULL){
		fclose(dataFile);
	}
}

/*
    readMatrix
        
    input char pointer to file containing comma seperated floating point numbers
    input dataMatrix to read data into
    Reads data from data.txt into 2d array of doubles
    Does the same as readData, but works with matrix of any dimensions
*/
double **readMatrix(char *filename, double **mat) {
    //Open data file
    FILE *fp = fopen(filename, "r");
    
    //Ensure file could be opened successfully
    if(fp == NULL) {
        fprintf(stderr, "Failed to open data file\n");
        exit(EXIT_FAILURE);
    }
    
    //Count number of rows and columns in matrix
    int curRow = 0;
    int curCol = 0;
    
    //Read file row by row
    char line[10000];
    while(fgets(line, 10000, fp)) {
        //Counts number of columns read on current line
        curCol = 0;
        
        //Reallocate memory (provide additional memory for another row)
        mat = (double **) realloc(mat, (curRow+1) * sizeof(double *));
        
        //Allocate memory for current row, starting with 1 column
        mat[curRow] = (double *) malloc(sizeof(double) * (curCol +1));
        
        //Store string in pointer, rather than array
        char *data = line;
        //Offset into string (for scanf)
        int offset;
        
        //Read all entries of current row into array
        while(sscanf(data, "%lf%*c%n ", &mat[curRow][curCol], &offset) == 1) {
            //Reallocate larger memory block for current row (add another column space)
            //TEST OUTPUT
            //printf("Read %f at index %d:%d\n", mat[curRow][curCol], curRow, curCol);
            curCol++;
            mat[curRow] = realloc(mat[curRow], sizeof(double) * (curCol+1));
            data += offset;
        }
        //Increment row counter
        curRow++;
    }
    
    //Calculate number of rows and columns
    ROWS = curRow;
    COLS = curCol;
    
    //Return matrix
    return mat;
}

/*
	readKeys

	input char pointer to file containing space separated long longs
	Reads data from KEY_FILE into array of long long ROW in size
*/
void readKeys(char *filename, long long *keyDatabase){
	FILE *dataFile = fopen(filename, "r");

	//	Null file pointer check
	if(dataFile == NULL){
		fprintf(stderr, "%s - The file \"%s\" could not be opened\n", 
			programName, filename);

		exit(EXIT_FAILURE);
	}

	//	Read keys, separeted by whitespace
	for(int k = 0; k < ROWS; ++k){
		fscanf(dataFile, "%lld", &keyDatabase[k]);
	}

	//	Close file pointer if open
	if(dataFile != NULL){
		fclose(dataFile);
	}
}

/*
    freeData
 
    input dynamically allocated matrix and keyDatabase
    Frees dynamically allocated memory
*/
void freeData(double **mat, long long *keys) {
    //Free all rows of matrix
    for(int i=0; i<ROWS; i++) {
        free(mat[i]);
    }
    //Free matrix
    free(mat);
    
    //Free key database
    free(keys);
}

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
 freeCollisionDB
 
 input collision database and number of collisions
 Frees all memory allocated for collision database, and individual collisions
 */
void freeCollisionDB(Collision **cdb, int numCollisions) {
    //Free each collision
    for(int i=0; i<numCollisions; i++) {
        //First, free each collision's block database (don't need to free individual blocks, as they're already freed in freeBD())
        free(cdb[i]->collidingBlocks);
        //Now, free collision
        free(cdb[i]);
    }
    //Now, free database memory
    free(cdb);
}