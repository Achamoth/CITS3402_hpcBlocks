#include <stdio.h>
#include <stdlib.h>

float **readMatrix(int *, int *, float **);
void printMatrix(float **, int,  int);
void freeMatrix(float **, int);

int main(void) {
    int rows;
    int cols;
    
    //Set up matrix
    float **matrix = (float **) malloc(sizeof(float *) * 1);
    
    //Load data into matrix
    matrix = readMatrix(&rows, &cols, matrix);
    
    
    //Do the project here
    //TODO: PROJCET
    
    
    //Print matrix (test output to confirm reading is working properly)
    //printMatrix(matrix, rows, cols);
    
    //Free matrix
    freeMatrix(matrix, rows);
    
    //Exit program
    exit(EXIT_SUCCESS);
}

//Function receives matrix as argument, and prints it to stdout
void printMatrix(float **mat, int rows, int cols) {
    //Loop through rows
    for(int i=0; i<rows; i++) {
        //Loop through columns
        for(int j=0; j<cols; j++) {
            printf("%3f ", mat[i][j]);
        }
        printf("\n");
    }
}

//Function receives matrix as argument, and frees memory allocated for matrix
void freeMatrix(float **mat, int rows) {
    //Loop through rows
    for(int i=0; i<rows; i++) {
        //Free current row
        free(mat[i]);
    }
    //Free matrix
    free(mat);
}

//Function opens input file and reads matrix into variable
float **readMatrix(int *rows, int *cols, float **mat) {
    //Open data file
    FILE *fp = fopen("test.txt", "r");
    
    //Ensure file could be opened successfully
    if(fp == NULL) {
        fprintf(stderr, "Failed to open data file\n");
        exit(EXIT_FAILURE);
    }
    
    //Count number of rows and columns in matrix
    int curRow = 0;
    int curCol = 0;
    
    //Read file row by row
    char line[BUFSIZ];
    while(fgets(line, BUFSIZ, fp)) {
        //Counts number of columns read on current line
        curCol = 0;
        
        //Reallocate memory (provide additional memory for another row)
        mat = (float **) realloc(mat, (curRow+1) * sizeof(float *));
        
        //Allocate memory for current row, starting with 1 column
        mat[curRow] = (float *) malloc(sizeof(float) * (curCol +1));
        
        //Store string in pointer, rather than array
        char *data = line;
        
        int offset;
        
        //Read all entries of current row into array
        while(sscanf(data, "%f , %n ", &mat[curRow][curCol], &offset) == 1) {
            //Reallocate larger memory block for current row (add another column space)
            printf("Read %f at index %d:%d\n", mat[curRow][curCol], curRow, curCol);
            curCol++;
            mat[curRow] = realloc(mat[curRow], sizeof(float) * (curCol+1));
            data += offset;
        }
        //Increment row counter
        curRow++;
    }
    
    //Calculate number of rows and columns
    *rows = curRow;
    *cols = curCol;
    
    //Return matrix
    return mat;
}
