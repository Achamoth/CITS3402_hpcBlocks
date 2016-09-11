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
void readData(char *filename, double dataMatrix[ROWS][COLS]){
	FILE *dataFile = fopen(filename, "r");

	//	Null file pointer check
	if(dataFile == NULL){
		fprintf(stderr, "%s - The file \"%s\" could not be opened\n", 
			programName, filename);

		exit(EXIT_FAILURE);
	}

	for(int r = 0; r < ROWS; ++r){
		for(int c = 0; c < COLS; ++c){
			fscanf(dataFile, "%lf[^,\n]", &dataMatrix[r][c]);
			printf("%lf\n", dataMatrix[r][c]);
		}
	}



	//	Close file pointer if open
	if(dataFile != NULL){
		fclose(dataFile);
	}
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


int main(int argc, char** argv){
	// Program name without /, cast to constant for file
	programName = (const char*) strrchr(argv[0], '/') + 1;

	double dataMatrix[ROWS][COLS] = {{0}};
	readData(DATA_FILE, dataMatrix);

	long long keyDatabase[ROWS];
	readKeys(KEY_FILE, keyDatabase);

	return EXIT_SUCCESS;
}