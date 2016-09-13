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
    double sumOfElements;
    int column;
} Block;

//------------------------------------------------------------------
// Struct declaration for collisions
//------------------------------------------------------------------
typedef struct Collision {
    Block **collidingBlocks;
    int numBlocksInCollision;
} Collision;

//------------------------------------------------------------------
// Package accessible functions
//------------------------------------------------------------------
extern void readData(char *, double **);
extern void readKeys(char *, long long *);
extern double **readMatrix(char *, double **);
extern void freeData(double **, long long *);
extern void freeBD(Block **, int);
extern void freeCollisionDB(Collision **, int);
extern Block **findBlocks(Block **, double **, long long *, int *);
extern Collision **findCollisions(Block **, int, int *);

//------------------------------------------------------------------
// Package accessible variables and definitions
//------------------------------------------------------------------
#define DATA_FILE "data.txt"
#define KEY_FILE "keys.txt"
#define DIA 0.000001
extern const char* programName;
extern int ROWS;
extern int COLS;
