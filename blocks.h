/*
	 CITS3402 Project 1 2016
	 Name:			Pradyumn Vij, Ammar Abu Shamleh
	 Student number: 21521274, 21469477
	 Date:           September 2016
*/
#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <string.h>
#include <math.h>
#include <omp.h>

//------------------------------------------------------------------
// Struct declaration for blocks
//------------------------------------------------------------------
typedef struct Block {
    long long signature;
    int column;
} Block;

//------------------------------------------------------------------
// Struct declaration for collisions
//------------------------------------------------------------------
typedef struct Collision {
    int *columns;
    int numBlocksInCollision;
} Collision;

//------------------------------------------------------------------
// Struct declaration for pair
//------------------------------------------------------------------
typedef struct pair{
	double value;
	long long key;
} pair;

//------------------------------------------------------------------
// Package accessible functions
//------------------------------------------------------------------
extern void readData(char *, double **);
extern void readKeys(char *, long long *);
extern double **readMatrix(char *, double **);
extern double **transposeMatrix(double **);
extern void freeData(double **);
extern void freeCollisionDB(Collision *, int);
extern void freeTransposedData(double **);
extern Block *findBlocks(Block *, double **, long long *, int *);
extern Block *findBlocksParallel(Block *, double **, long long *, int *);
extern Block *findBlocksOptimised(Block *, double **, long long *, int *);
extern Collision *findCollisions(Block *, int, int *);
extern Collision *findCollisionsParallel(Block *, int, int *);
extern Collision *findCollisionsOptimised(Block *, int, int *);
extern Collision *findCollisionsOptimisedParallel(Block *, int, int *);

//------------------------------------------------------------------
// Package accessible variables and definitions
//------------------------------------------------------------------
#define DATA_FILE "data.txt"
#define KEY_FILE "keys.txt"
#define DIA 0.0000025
#define NUM_THREADS 4
extern const char* programName;
extern int ROWS;
extern int COLS;
