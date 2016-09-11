/*
	CITS3402 Project 1 2016
	Name:			Ammar Abu Shamleh, Pradyumn Vij 
	Student number: 21521274, 21469477
    Date:           September 2016
*/
#include "blocks.h"
//	Program Name
const char* programName;

    
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
			printf("%lf found at index %d : %d\n", dataMatrix[r][c], r, c);
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
		//printf("key found = %lld at index %d\n", keyDatabase[k], k);
	}

	//	Close file pointer if open
	if(dataFile != NULL){
		fclose(dataFile);
	}
}

/*
    freeData
    
    input dynamically allocated matrix and keyDatabase
    Frees dynamically allocated data
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


int main(int argc, char** argv){
	// Program name without /, cast to constant for file
	programName = (const char*) strrchr(argv[0], '/') + 1;

	// Allocate memory for keys
	long long *keyDatabase = malloc(ROWS*sizeof(long long));
	// Read in keys
	readKeys(KEY_FILE, keyDatabase);

	// Allocate memory for database
	double **dataMatrix = (double**)malloc(ROWS*sizeof(double*));
	for(int r = 0; r < ROWS; ++r){
		dataMatrix[r] = (double*)malloc(COLS*sizeof(double));
	}

    //Read in matrix, and record number of rows and columns in ROWS and COLS
	//readData(DATA_FILE, dataMatrix);
    dataMatrix = readMatrix(DATA_FILE, dataMatrix);
    
    //Free all dynamically allocated data
    freeData(dataMatrix, keyDatabase);

	return EXIT_SUCCESS;
}
