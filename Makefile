# CITS3402 Project 1 2016
# Name:                 Pradyumn Vij, Ammar Abu Shamleh
# Student number(s):    21469477, 21521274
# Date:                 September 2016

#   Make file for CITS3402 project 1 program
EXEC = blocks
HEADERS = $(wildcard *.h)
SOURCES = $(wildcard *.c)
OBJECTS = $(wildcard *.o)

CC = gcc
CFLAGS = -std=c99 -Wall -pedantic -Werror -fopenmp

$(EXEC): $(OBJECTS)
	$(CC) $(CFLAGS) $(SOURCES) -o $(EXEC)

$(OBJECTS): $(SOURCES) $(HEADERS)
	$(CC) $(CFLAGS) -c $(SOURCES)

clean:
	rm -f $(OBJECTS) $(EXEC)