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
    transposeMatrix
    
    input data matrix
    Tranposes data matrix
*/
double **transposeMatrix(double **dataMatrix) {
    //Allocate memory for transposed matrix
    double **transposed = (double **) malloc(sizeof(double *) * COLS);
    for(int i=0; i<COLS; i++) {
        transposed[i] = (double *) malloc(sizeof(double) * ROWS);
    }
    
    //Store data in transposed matrix
    for(int col=0; col<COLS; col++) {
        for(int row=0; row<ROWS; row++) {
            transposed[col][row] = dataMatrix[row][col];
        }
    }
    
    //Return pointer to transposed matrix
    return transposed;
}

/*
    freeData
 
    input dynamically allocated matrix
    Frees dynamically allocated memory
*/
void freeData(double **mat) {
    //Free all rows of matrix
    for(int i=0; i<ROWS; i++) {
        free(mat[i]);
    }
    //Free matrix
    free(mat);
}

/*
 freeTransposedData
 
 input dynamically allocated transposed
 Frees dynamically allocated memory
 */
void freeTransposedData(double **mat) {
    //Free all rows of matrix
    for(int i=0; i<COLS; i++) {
        free(mat[i]);
    }
    //Free matrix
    free(mat);
}

/*
 freeBD
 
 input blockDatabase: dynamically allocated database of Blocks
 Frees memory allocated to block database
 
 */
void freeBD(Block *bd, int numBlocks) {
    //Loop through all blocks
    for(int i=0; i<numBlocks; i++) {
        //Free each block's row database
        free(bd[i].rows);
    }
    //Free block database
    free(bd);
}

/*
 freeCollisionDB
 
 input collision database and number of collisions
 Frees all memory allocated for collision database, and individual collisions
 */
void freeCollisionDB(Collision *cdb, int numCollisions) {
    //Free each collision
    for(int i=0; i<numCollisions; i++) {
        //First, free each collision's column database
        free(cdb[i].columns);
        //Then, free each collision's Block database
        free(cdb[i].blocks);
    }
    //Now, free database memory
    free(cdb);
}

/*
 freeMergedDB
 
 Input MergedCollision database and number of merged collisions
 Frees all memory allocated for merged collision database
 */
void freeMergedDB(MergedCollision *mdb, int numMerged) {
    //Loop over all merged collisions
    for(int i=0; i<numMerged; i++) {
        //For each merged collision, free the column and row databases
        free(mdb[i].columns);
        free(mdb[i].rows);
    }
    //Free database
    free(mdb);
}

/*
 printBlock

 input Block struct
 Prints out all elements of block, along with its column number and its signature
 */
void printBlock(Block b) {
    //For clean spacing
    printf("\n");
    //First print block's signature and column numbers
    printf("Signature:          %16lld\n", b.signature);
    printf("Column:             %16d\n", b.column);
    //Next, print block's elements
    printf("Row 1:              %16d\n", b.rows[0]);
    printf("Row 2:              %16d\n", b.rows[1]);
    printf("Row 3:              %16d\n", b.rows[2]);
    printf("Row 4:              %16d\n", b.rows[3]);
    //For clean spacing
    printf("\n");
}
