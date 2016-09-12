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
void freeCollisionDB(Collision **cdb, int numCols) {
    //Free each collision
    for(int i=0; i<numCols; i++) {
        //First, free each collision's block database
        freeBD(cdb[i]->collidingBlocks, cdb[i]->numBlocksInCollision);
        //Now, free collision
        free(cdb[i]);
    }
    //Now, free database memory
    free(cdb);
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
//NOT WORKING PROPERLY (SEGFAULT, WHICH I THINK OCCURS HERE)
Collision **findCollisions(Block **blockDB, int numBlocks, int *numberCollisionsFound) {
    //Set up collision database
    int numCollisions = 0;
    Collision **collisions = (Collision **) malloc((numCollisions+1) * sizeof(Collision *));
    //Allocate memory for first entry of collision database
    collisions[0] = (Collision *) malloc(sizeof(Collision));
    
    //Loop over all blocks
    for(int i=0; i<numBlocks; i++) {
        Block *curBlock = blockDB[i];
        long long curSig = curBlock->signature;
        int curCollisions = 0;
        printf("cur block index = %d\n", i);
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
                    collisions[numCollisions] = (Collision *) malloc(sizeof(Collision));
                    collisions[numCollisions]->collidingBlocks = (Block **) malloc(sizeof(Block *));
                    collisions[numCollisions]->numBlocksInCollision = 0;
                }
                collisions[numCollisions]->collidingBlocks[curCollisions] = (Block *) malloc(sizeof(Block));
                collisions[numCollisions]->collidingBlocks[curCollisions] = compBlock;
                collisions[numCollisions]->numBlocksInCollision += 1;
                curCollisions++;
                numCollisions++;
                //Allocate more memory for collision database
                collisions = (Collision **) realloc(collisions, (numCollisions+1)*sizeof(Collision *));
                
                //TEST OUTPUT
                printf("Found collision at block %d and %d. Sigs are: %lld and %lld\n", i, j, curSig, compBlock->signature);
                
            }
        }
    }
    *numberCollisionsFound = numCollisions;
    return collisions;
}

int main(int argc, char** argv){
	// Program name without /, cast to constant for file
	programName = (const char*) strrchr(argv[0], '/') + 1;

	// Allocate memory for keys
	long long *keyDatabase = malloc(ROWS*sizeof(long long));
	// Read in keys
	readKeys(KEY_FILE, keyDatabase);
    
    //Allocate memory for matrix database
	double **dataMatrix = (double**)malloc(1*sizeof(double*));
    //Read in matrix, and record number of rows and columns in ROWS and COLS
    dataMatrix = readMatrix(DATA_FILE, dataMatrix);
    
    //Create pool of blocks
    Block **blockDatabase = malloc(1*sizeof(Block *));;
    
    //Find all blocks in matrix
    int numBlocks = 0;
    blockDatabase = findBlocks(blockDatabase, dataMatrix, keyDatabase, &numBlocks);
    
    //Find all collisions among blocks
    int numCollisions = 0;
    Collision **collisions = findCollisions(blockDatabase, numBlocks, &numCollisions);
    
    //Free all dynamically allocated memory for key and matrix databases
    freeData(dataMatrix, keyDatabase);
    //Free dynamically allocated memory for block database
    freeBD(blockDatabase, numBlocks);
    //Free dynamically allocated memory for collision database
    freeCollisionDB(collisions, numCollisions);
    //Exit program
	return EXIT_SUCCESS;
}
