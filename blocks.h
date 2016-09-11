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

//------------------------------------------------------------------
// Package accessible functions
//------------------------------------------------------------------
//extern void readData(char*, double**);
extern void readKeys(char*, long long*);
//------------------------------------------------------------------
// Package accessible variables and definitions
//------------------------------------------------------------------
#define DATA_FILE "data.txt"
#define KEY_FILE "keys.txt"
extern const char* programName;
int ROWS = 4400;
int COLS = 500;
double DIA = 0.000001;
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
