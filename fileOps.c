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
    openFile

    Opens file into a filestream with error checking
*/
FILE *openFile(char *filename){
    FILE *dataFile = fopen(filename, "r");
    if(dataFile == NULL){
		fprintf(stderr, "%s - The file \"%s\" could not be opened\n",
			programName, filename);

		exit(EXIT_FAILURE);
	}

    return dataFile;
}

/*
    closeFile

    Closes file streams, performing error check to see if open
*/
void closeFile(FILE *filePointer){
    if(filePointer == NULL){
        fclose(filePointer);
    }
}

/*
    countKeys

    Counts the number of large integers in the KEY_FILE
*/
void countKeys(FILE *dataFile){
    // Count rows
    long long temp;
    ROWS = 0;
    while(fscanf(dataFile, "%lld", &temp) > 0){
        ROWS++;
    }
}

/*
	readKeys

	input char pointer to file containing space separated long longs
	Reads data from KEY_FILE into array of long long ROW in size and returns
*/
void readKeys(FILE *dataFile, long long *keys){
    // Read keys, separeted by whitespace
    for(int k = 0; k < ROWS; ++k){
        fscanf(dataFile, "%lld", &keys[k]);
    }
}

/*
    countRowsCols

    Counts the rows and columns in the CSV file provided.
    A single line has n floating point numbers which equates to the number of
    columns.
    The number of lines in the file equates to the number of rows required by
    the data.
*/
void countRowsCols(FILE *dataFile){
    // Count columns
    // Count columns in first line
    COLS = 0;
    char buffer[10000];
    char *line = buffer;
    if(fgets(buffer, sizeof(buffer), dataFile) != NULL){
        double temp;
        int offset = 0;
        while(sscanf(line, "%lf%*c%n", &temp, &offset) == 1){
            COLS++;
            line += offset;
        }
    }
    // Rewind file stream
    rewind(dataFile);
    // Check number of rows in case incorrect or formatting errors in file
    int dataRows = 0;
    while(fgets(buffer, sizeof(buffer), dataFile) != NULL){
        dataRows++;
    }
    // Error check to see if proper key file is being used with data file
    if(dataRows != ROWS){
        fprintf(stderr, "%s Error: Number of rows counted in data file (%d) \
        does not equal the number of keys counted (%d). Exiting.", programName,
        dataRows, ROWS);
        exit(EXIT_FAILURE);
    }
    // Rewind file stream
    rewind(dataFile);
}

/*
	readData

	input char pointer to file containing comma separated floating point numbers
	Reads data from data.txt into 2d array of doubles ROWS x COLS in size
*/
void readData(FILE *dataFile, double **dataMatrix){
	for(int r = 0; r < ROWS; ++r){
		for(int c = 0; c < COLS; ++c){
			//  Use fscanf to skip the comma characters
			fscanf(dataFile, "%lf%*c", &dataMatrix[r][c]);
            //TEST OUTPUT
			//printf("%lf found at index %d : %d\n", dataMatrix[r][c], r, c);
		}
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
