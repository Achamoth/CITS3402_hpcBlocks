/*
 CITS3402 Project 1 2016
 Name:           Ammar Abu Shamleh, Pradyumn Vij
 Student number: 21521274, 21469477
 Date:           September 2016
 */
#include "blocks.h"

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
 mergeBlockDatabases

 input partial block database and complete block database
 Copies all data (actually copies pointers) from partial block database to complete database
 */
Block *mergeBlockDatabases(Block *completeDB, Block *partialDB, int numBlockInPartial, int *numCopiedToComplete) {
    //Reallocate more memory for complete database
    completeDB = (Block *) realloc(completeDB, sizeof(Block) * (*numCopiedToComplete+numBlockInPartial));
    for(int i=0; i<numBlockInPartial; i++) {
        completeDB[i + *numCopiedToComplete] = partialDB[i];
    }
    *numCopiedToComplete += numBlockInPartial;
    return completeDB;
}

/*
 getBlockValues

 input column and four row numbers, and complete data matrix
 Finds four values corresponding to column number and four row numbers, and returns a double array containing them
 */
double *getBlockValues(int c, int r1, int r2, int r3, int r4, double **mat) {
    //First, allocate memory for element database
    double *blockElements = (double *) malloc(sizeof(double) * 4);
    //Now, find four values and store them in database
    blockElements[0] = mat[c][r1];
    blockElements[1] = mat[c][r2];
    blockElements[2] = mat[c][r3];
    blockElements[3] = mat[c][r4];
    //Return database
    return blockElements;
}

/*
 getRows

 input four row values
 Returns dynamically allocated int array consisting of four row values
 */
int *getRows(int r1, int r2, int r3, int r4) {
    //First, allocate memory for row database
    int *rows = (int *) malloc(sizeof(int) * 4);
    //Now, store four row values in database
    rows[0] = r1;
    rows[1] = r2;
    rows[2] = r3;
    rows[3] = r4;
    //Return database
    return rows;
}

/*
 findBlocks

 input blockDatabase and matrixDatabase
 Finds all blocks in matrixDatabase and stores them in blockDatabase
 */
Block *findBlocks(Block *blockDB, double **mat, long long *kd, int *numBlocks) {
    int nextBlock = 0;
    //Loop through matrix columns
    for(int col=0; col<COLS-1; col++) {
        //Loop through matrix rows
        for(int row1=0; row1<ROWS; row1++) {
            for(int row2=row1+1; row2<ROWS; row2++) {
                //Ensure row1 and row2 are unique
                if(row2 == row1) continue;
                //Check if they're in the same neighbourhood
                if(fabs(mat[col][row1] - mat[col][row2])>DIA) continue;
                for(int row3=row2+1; row3<ROWS; row3++) {
                    //Ensure rows 1, 2 and 3 are unique
                    if(row3 == row2 || row3 == row1) continue;
                    //Check they're in the same neighbourhood
                    if(fabs(mat[col][row3]-mat[col][row2])>DIA || fabs(mat[col][row3]-mat[col][row1])>DIA) continue;
                    for(int row4=row3+1; row4<ROWS; row4++) {
                        //Ensure all rows are unique
                        if(row4==row1 || row4==row2 || row4==row3) continue;
                        //Check they're in the same neighbourhood
                        if(fabs(mat[col][row4]-mat[col][row1])>DIA || fabs(mat[col][row4]-mat[col][row2])>DIA || fabs(mat[col][row4]-mat[col][row3])>DIA) continue;
                        //We have found a block, and must store it in the block database
                        blockDB = (Block *) realloc(blockDB, (nextBlock+1)*sizeof(Block));
                        blockDB[nextBlock].signature = findSig(row1, row2, row3, row4, kd);
                        blockDB[nextBlock].column = col;
                        //blockDB[nextBlock].values = getBlockValues(col, row1, row2, row3, row4, mat);
//                        blockDB[nextBlock].rows = getRows(row1, row2, row3, row4);
                        blockDB[nextBlock].rows[0] = row1;
                        blockDB[nextBlock].rows[1] = row2;
                        blockDB[nextBlock].rows[2] = row3;
                        blockDB[nextBlock].rows[3] = row4;
                        nextBlock++;
                        //TEST OUTPUT
                        //printf("Found block at column %d on rows %d, %d, %d, %d\n", col, row1, row2, row3, row4);
                    }
                }
            }
        }
    }
    *numBlocks = nextBlock;
    return blockDB;
}

/*
 findBlocksParallel

 input blockDatabase and matrixDatabase
 Finds all blocks in matrixDatabase and stores them in blockDatabase. Parallelized
 */
Block *findBlocksParallel(Block *blockDB, double **mat, long long *kd, int *numBlocks, int start, int end) {
    int nextBlock = 0;
    //Loop through matrix columns (excluding last column)
    for(int col=start; col<end; col++) {
        //Loop through matrix rows in parallel
        #pragma omp parallel
        {
            Block *partialBlockDB = (Block *) malloc(1 * sizeof(Block));
            int localNextBlock = 0; //Thread private counter for next index into block database
            int ID = omp_get_thread_num();
            int numThreads = omp_get_num_threads();
            for(int row1=ID; row1<ROWS; row1+=numThreads) {
                for(int row2=row1+1; row2<ROWS; row2++) {
                    //Check if they're in the same neighbourhood
                    if(fabs(mat[col][row1] - mat[col][row2])>DIA) continue;
                    for(int row3=row2+1; row3<ROWS; row3++) {
                        //Check they're in the same neighbourhood
                        if(fabs(mat[col][row3]-mat[col][row2])>DIA || fabs(mat[col][row3]-mat[col][row1])>DIA) continue;
                        for(int row4=row3+1; row4<ROWS; row4++) {
                            //Check they're in the same neighbourhood
                            if(fabs(mat[col][row4]-mat[col][row1])>DIA || fabs(mat[col][row4]-mat[col][row2])>DIA || fabs(mat[col][row4]-mat[col][row3])>DIA) continue;
                            //We have found a block, and must store it in the block database
                            partialBlockDB = (Block *) realloc(partialBlockDB, (localNextBlock+1)*sizeof(Block));
                            partialBlockDB[localNextBlock].signature = findSig(row1, row2, row3, row4, kd);
                            partialBlockDB[localNextBlock].column = col;
                            //partialBlockDB[localNextBlock].values = getBlockValues(col, row1, row2, row3, row4, mat);
//                            partialBlockDB[localNextBlock].rows = getRows(row1, row2, row3, row4);
                            partialBlockDB[localNextBlock].rows[0] = row1;
                            partialBlockDB[localNextBlock].rows[1] = row2;
                            partialBlockDB[localNextBlock].rows[2] = row3;
                            partialBlockDB[localNextBlock].rows[3] = row4;
                            localNextBlock++;
                            //TEST OUTPUT
                            //printf("Thread %d: Found block at column %d on rows %d, %d, %d, %d\n", ID, col, row1, row2, row3, row4);
                        }
                    }
                }
            }
            //Merge partial block database with complete block database
            #pragma omp critical
            {
                blockDB = mergeBlockDatabases(blockDB, partialBlockDB, localNextBlock, &nextBlock);
            }
            //Free partial block database
            free(partialBlockDB);
        }
    }

    *numBlocks = nextBlock;
    return blockDB;
}

/*
 findBlocksParallelMPI

 input blockDatabase and matrixDatabase
 Finds all blocks in matrixDatabase and stores them in blockDatabase. Parallelized
 */
Block *findBlocksParallelMPI(Block *blockDB, double **mat, long long *kd, int *numBlocks) {
    //Create partial block database for each node
    Block *block = malloc(sizeof(Block));
    
    //Now, set up MPI struct for block
    MPI_Datatype BlockType;
    MPI_Datatype type[3] = {MPI_LONG_LONG_INT, MPI_INT, MPI_INT};
    int blocklen[3] = {1,1,4};
    MPI_Aint disp[3];
    disp[0] = block->signature - partialBlockDatabase;
    disp[1] = block->column - partialBlockDatabase;
    disp[2] = block->&rows[0] = partialBlockDatabase;
    MPI_Type_create_struct(3, blocklen, disp, type, &BlockType);
    MPI_Type_commit(&BlockType);
    
    //Set up MPI struct for data matrix
    DataMatrixColumn *columns = malloc(sizeof(DataMatrixColumn));
    MPI_Datatype MatColumnType;
    MPI_Datatype colType[3] = {MPI_DOUBLE};
    int colLen[1] = {4400};
    MPI_Aint colDisp[1] = {column->&rows[0] - column};
    MPI_Type_create_struct(1, colLen, colDisp, colType, &MatColumnType);
    MPI_Type_commit(&MatColumnType);
    
    //Initialize MPI
    MPI_Init(NULL, NULL);
    
    //Create partial block database for current node
    Block *partialBlockDatabase = malloc(sizeof(Block));
    
    //Find world size and current rank
    int commSize;
    int commRank;
    MPI_COMM_RANK(MPI_COMM_WORLD, &commRank);
    MPI_COMM_SIZE(MPI_COMM_WORLD, &commSize);
    
    if(commSize == 1) {
        //Just compute using the omp parallel algorithm
        free(partialBlockDatabase);
        blockDB = findBlocksParallel(blockDB, mat, kd, numBlocks, 0, COLS);
        return blockDB;
    }
    
    if(commRank == MASTER && commSize > 1) {
        /* Master process. This process governs all the others */
        
        //First step is to split the data matrix into chunks of columns, using comm size
        for(int i=1; i<commSize; i++) {
            //For each process, compute its start and end values (the columns it has to process)
            int blockSize = COLS / (commSize-1);
            int start = (i-1) * blockSize;
            int end = ((i-1) * blockSize) + blockSize;
            
            //Send each process its start and end values
            MPI_SEND(&start, 1, MPI_INT, i, 0, MPI_COMM_WORLD);
            MPI_SEND(&end, 1, MPI_INT, i, 0, MPI_COMM_WORLD);
        }
        
        
        //Each process will now process its columns, and then send back an array of blocks
        
    }
    
    else {
        int start;
        int end;
        int status;
        //Master will send start and end values for columns
        MPI_RECV(&start, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);
        MPI_RECV(&end, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);
        
        //Next, this process has to compute all blocks over columns within start and end values
        
        
        //Now, the process will send back the partial blocks that it found
        
    }
}

/*
    compareDoubles
    Comparator function for qSort, to sort doubles in ascending order
*/
int compareDoubles(const void* a, const void* b){
    pair resA = *(pair*)a;
    pair resB = *(pair*)b;
    if(resA.value < resB.value) return -1;
    if(resA.value > resB.value) return 1;
    return 0;
}


/*
    findBlocksOptimised
    input blockDatabase and matrixDatabase
    Finds all blocks in matrixDatabase and stores them in blockDatabase
    Sequential optimised code
*/
Block *findBlocksOptimised(Block *blockDB, double **mat, long long *kd, int *numBlocks) {
    int nextBlock = 0;
    //Loop through matrix columns
    // double* tempContainer = malloc(ROWS * sizeof(double));
    pair* pairContainer = malloc(ROWS * sizeof(pair));
    for(int col=0; col< COLS - 1; col++) {
        // Create an array of doubles
        for(int row = 0; row < ROWS; ++row){
            //tempContainer[row] = mat[row][col];
            pairContainer[row].value = mat[col][row];
            pairContainer[row].key = kd[row];
            pairContainer[row].row = row;
        }

        // Use sliding technique to fill BlockDB
        // lower bound in array
        int lower = 0;
        // sort the entire array in ascending order O(nlgn)
        qsort(pairContainer, ROWS, sizeof(pair), compareDoubles);
        // Upper bound incrementing is growing the size of the window, default start at 1
        // Go out of bounds with inclusion of boundary ROWS
        for(int upper = 1; upper <= ROWS; ++upper){
            // Check if window contains elements in the same neighbourhood.
            // If not move the lower bound up till it does or until upper == lower
            while(pairContainer[upper-1].value - pairContainer[lower].value > DIA){
                ++lower;
                // Take note that the lower bound i.e. new window instance
            }
            // Get all combinations within the window of size 4
            for(int i = lower; i < upper-1; ++i){
                for(int j = i+1; j < upper-1; ++j){
                    for(int k = j+1; k < upper-1; ++k){
                        //  Increase memory for the new block pointer to the database
                        blockDB = (Block *) realloc(blockDB, (nextBlock+1)*sizeof(Block));
                        blockDB[nextBlock].signature = pairContainer[i].key + pairContainer[j].key + pairContainer[k].key + pairContainer[upper-1].key;
                        blockDB[nextBlock].column = col;
                        //blockDB[nextBlock].values = getBlockValues(col, i, j, k, upper-1, mat);
//                        blockDB[nextBlock].rows = getRows(pairContainer[i].row, pairContainer[j].row, pairContainer[k].row, pairContainer[upper-1].row);
                        blockDB[nextBlock].rows[0] = pairContainer[i].row;
                        blockDB[nextBlock].rows[1] = pairContainer[j].row;
                        blockDB[nextBlock].rows[2] = pairContainer[k].row;
                        blockDB[nextBlock].rows[3] = pairContainer[upper-1].row;
                        //  Increment number of blocks
                        nextBlock++;
                        // Uncomment to print all rows / indexes being found
                        //printf("Found block at column %d on rows %lld, %lld, %lld, %lld\n", col, pairContainer[i].key, pairContainer[j].key, pairContainer[k].key, pairContainer[upper-1].key);
                    }
                }
            }
        }
    }
    // Free the utility container
    *numBlocks = nextBlock;
    free(pairContainer);
    //printf("%d\n", *numBlocks);
    return blockDB;
}

/*
 findBlocksOptimisedParallel
 input blockDatabase and matrixDatabase
 Finds all blocks in matrixDatabase and stores them in blockDatabase
 Parallel optimised code
 */
Block *findBlocksOptimisedParallel(Block *collectiveBlockDB, double **mat, long long *kd, int *numBlocks) {
    int numBlocksTotal = 0;
    //Loop through matrix columns
    #pragma omp parallel for
    for(int col=0; col< COLS - 1; col++) {
        //Thread-private partial block database and counter
        Block *blockDB = malloc(sizeof(Block));
        int nextBlock = 0;
        pair* pairContainer = malloc(ROWS * sizeof(pair));
        // Create an array of doubles
        for(int row = 0; row < ROWS; ++row){
            //tempContainer[row] = mat[row][col];
            pairContainer[row].value = mat[col][row];
            pairContainer[row].key = kd[row];
            pairContainer[row].row = row;
        }

        // Use sliding technique to fill BlockDB
        // lower bound in array
        int lower = 0;
        // sort the entire array in ascending order O(nlgn)
        qsort(pairContainer, ROWS, sizeof(pair), compareDoubles);
        // Upper bound incrementing is growing the size of the window, default start at 1
        // Go out of bounds with inclusion of boundary ROWS
        for(int upper = 1; upper <= ROWS; ++upper){
            // Check if window contains elements in the same neighbourhood.
            // If not move the lower bound up till it does or until upper == lower
            while(pairContainer[upper-1].value - pairContainer[lower].value > DIA){
                ++lower;
                // Take note that the lower bound i.e. new window instance
            }
            // Get all combinations within the window of size 4
            for(int i = lower; i < upper-1; ++i){
                for(int j = i+1; j < upper-1; ++j){
                    for(int k = j+1; k < upper-1; ++k){
                        //  Increase memory for the new block pointer to the database
                        blockDB = (Block *) realloc(blockDB, (nextBlock+1)*sizeof(Block));
                        blockDB[nextBlock].signature = pairContainer[i].key + pairContainer[j].key + pairContainer[k].key + pairContainer[upper-1].key;
                        blockDB[nextBlock].column = col;
                        //blockDB[nextBlock].values = getBlockValues(col, i, j, k, upper-1, mat);
//                        blockDB[nextBlock].rows = getRows(pairContainer[i].row, pairContainer[j].row, pairContainer[k].row, pairContainer[upper-1].row);
                        blockDB[nextBlock].rows[0] = pairContainer[i].row;
                        blockDB[nextBlock].rows[1] = pairContainer[j].row;
                        blockDB[nextBlock].rows[2] = pairContainer[k].row;
                        blockDB[nextBlock].rows[3] = pairContainer[upper-1].row;
                        //  Increment number of blocks
                        nextBlock++;
                        // Uncomment to print all rows / indexes being found
                        //printf("Found block at column %d on rows %lld, %lld, %lld, %lld\n", col, pairContainer[i].key, pairContainer[j].key, pairContainer[k].key, pairContainer[upper-1].key);
                    }
                }
            }
        }
        free(pairContainer);
        //Merge block database with collective block database
        #pragma omp critical
        {
            collectiveBlockDB = mergeBlockDatabases(collectiveBlockDB, blockDB, nextBlock, &numBlocksTotal);
        }
    }
    // Free the utility container
    *numBlocks = numBlocksTotal;
    //printf("%d\n", *numBlocks);
    return collectiveBlockDB;
}

/*
 findCollisionsParallel

 input blockDatabase
 Finds all collisions between generated blocks and return collision database; also store number of collisions found. Sequential brute-force code
 */
Collision *findCollisions(Block *blockDB, int numBlocks, int *numberCollisionsFound) {
    //Set up collision database
    int numCollisions = 0;
    Collision *collisions = (Collision *) malloc((numCollisions+1) * sizeof(Collision));
    //Record whether or not block has already been detected in a collision
    bool *collided = malloc(sizeof(bool) * numBlocks);
    for(int i=0; i<numBlocks; i++) collided[i] = false;

    //Loop over all blocks
    for(int i=0; i<numBlocks; i++) {
        Block curBlock = blockDB[i];
        long long curSig = curBlock.signature;
        int curCollisions = 0;
        //If current block has already been detected in a collision, skip it
        if(collided[i]) continue;
        //Start inner loop to compare to all other blocks
        for(int j=i+1; j<numBlocks; j++) {
            Block compBlock = blockDB[j];
            //Ensure that blocks are not in same column
            if(curBlock.column == compBlock.column) {
                continue;
            }
            //Check for collision
            if(compBlock.signature == curSig) {
                //Add block to current collision
                if(curCollisions == 0) {
                    //Allocate memory for collision's column database, and store first column
                    collisions[numCollisions].columns = (int *) malloc(5 * sizeof(int));
                    collisions[numCollisions].numBlocksInCollision = 1;
                    collisions[numCollisions].columns[curCollisions] = curBlock.column;
                    //Allocate memory for collision's block database, and store first block
                    collisions[numCollisions].blocks = (Block *) malloc(2 * sizeof(Block));
                    collisions[numCollisions].blocks[curCollisions] = curBlock;
                    curCollisions++;
                    //Allocate more memory for collision database
                    numCollisions++;
                    collisions = (Collision *) realloc(collisions, (numCollisions+1)*sizeof(Collision));
                }
                //Store next block in collision
                collided[j] = true;
                collisions[numCollisions-1].numBlocksInCollision += 1;
                collisions[numCollisions-1].columns[curCollisions] = compBlock.column;
                collisions[numCollisions-1].blocks[curCollisions] = compBlock;
                curCollisions++;
                //Allocate more memory for current collision's block and column databases
                collisions[numCollisions-1].columns = (int *) realloc(collisions[numCollisions-1].columns, ((curCollisions+1)*sizeof(int)));
                collisions[numCollisions-1].blocks = (Block *) realloc(collisions[numCollisions-1].blocks, (curCollisions+1) * sizeof(Block));

                //TEST OUTPUT
                //printf("%d %d: ", numCollisions-1, curCollisions);
                //printf("Found collision at blocks %d and %d. Cols are: %d and %d. Sigs are: %lld and %lld \n", i, j, curBlock.column, collisions[numCollisions-1].columns[curCollisions-1], curSig, compBlock.signature);

            }
        }
    }
    free(collided);
    *numberCollisionsFound = numCollisions;
    return collisions;
}

/*
 findCollisionsParallel

 input blockDatabase
 Finds all collisions between generated blocks and return collision database; also store number of collisions found. Parallelized brute-force code
 */
Collision *findCollisionsParallel(Block *blockDB, int numBlocks, int *numberCollisionsFound) {
    //Set up collision database
    int totNumCollisions = 0;
    Collision *collisionsCollective = (Collision *) malloc((totNumCollisions+1) * sizeof(Collision));
    //  Record whether or not block has already been detected in a collision
    //  bool *collided = malloc(sizeof(bool) * numBlocks);
    //  for(int i=0; i<numBlocks; i++) collided[i] = false;
    #pragma omp parallel
    {
        //Thread variables
        int ID = omp_get_thread_num();
        int nthreads = omp_get_num_threads();
        //Private partial collision database
        int numCollisions = 0;
        Collision *collisions = (Collision *) malloc((numCollisions+1) * sizeof(Collision));

        //Loop over all blocks
        for(int i=ID; i<numBlocks; i+=nthreads) {
            Block curBlock = blockDB[i];
            long long curSig = curBlock.signature;
            int curCollisions = 0;
            //If current block has already been detected in a collision, skip it
//            if(collided[i]) continue;
            //Start inner loop to compare to all other blocks
            for(int j=i+1; j<numBlocks; j++) {
                Block compBlock = blockDB[j];
                //Ensure that blocks are not in same column
                if(curBlock.column == compBlock.column) {
                    continue;
                }
                //Check for collision
                if(compBlock.signature == curSig) {
                    //Add block to current collision
                    if(curCollisions == 0) {
                        //Allocate memory for collision's column database, and store first column
                        collisions[numCollisions].columns = (int *) malloc(5 * sizeof(int));
                        collisions[numCollisions].numBlocksInCollision = 1;
                        collisions[numCollisions].columns[curCollisions] = curBlock.column;
                        //Allocate memory for collision's block database, and store first block
                        collisions[numCollisions].blocks = (Block *) malloc(2 * sizeof(Block));
                        collisions[numCollisions].blocks[curCollisions] = curBlock;
                        curCollisions++;
                        //Allocate more memory for collision database
                        numCollisions++;
                        collisions = (Collision *) realloc(collisions, (numCollisions+1)*sizeof(Collision));
                    }
                    //Store next block in collision
//                    collided[j] = true;
                    collisions[numCollisions-1].numBlocksInCollision += 1;
                    collisions[numCollisions-1].columns[curCollisions] = compBlock.column;
                    collisions[numCollisions-1].blocks[curCollisions] = compBlock;
                    curCollisions++;
                    //Allocate more memory for current collision's block and column databases
                    collisions[numCollisions-1].columns = (int *) realloc(collisions[numCollisions-1].columns, ((curCollisions+1)*sizeof(int)));
                    collisions[numCollisions-1].blocks = (Block *) realloc(collisions[numCollisions-1].blocks, (curCollisions+1) * sizeof(Block));

                    //TEST OUTPUT
                    //printf("Thread %d\t:", ID);
                    //printf("%d %d: ", numCollisions-1, curCollisions);
                    //printf("Found collision at blocks %d and %d. Cols are: %d and %d. Sigs are: %lld and %lld \n", i, j, curBlock.column, collisions[numCollisions-1].columns[curCollisions-1], curSig, compBlock.signature);

                }
            }
        }
    }
//    free(collided);
    *numberCollisionsFound = totNumCollisions;
    return collisionsCollective;
}

/*
 mergeCollisionDatabases

 input collective collision databases and partial collision databases, along with sizes of each
 Copies all collisions in partial database into collective database
*/
Collision *mergeCollisionDatabases(Collision *cdb, Collision *partial, int *totCollisions, int numInPartial) {
    //Reallocate more memory for collective database
    cdb = (Collision *) realloc(cdb, (numInPartial+ (*totCollisions)) * sizeof(Collision));

    //Loop over all entries in partial database
    for(int i=0; i<numInPartial; i++) {
        /*//Copy all data over to collective database
        cdb[i + *totCollisions].numBlocksInCollision = partial[i].numBlocksInCollision;
        //Allocate memory for column database
        cdb[i + *totCollisions].columns = (int *) malloc(sizeof(int) * partial[i].numBlocksInCollision);
        //Copy all column entries over
        for(int j=0; j<partial[i].numBlocksInCollision; j++) {
            cdb[i + *totCollisions].columns[j] = partial[i].columns[j];
        }*/

        //Copy collision data over
        cdb[i + *totCollisions] = partial[i];
    }

    //Update collective count of collisions
    *totCollisions += numInPartial;

    //Return new pointer to collective collision database
    return cdb;
}

/*
 cmpfunc

 Input; two variables of same datatype
 Performs comparison for sorting function. Returns positive if first value is larger; negative if second value is larger
 */
int cmpfunc(const void *a, const void *b) {
    //Find value of block a's signature
    long sigA = 0;
    Block ba = *(Block *) a;
    sigA = ba.signature;

    //Find value of block b's signature
    long sigB = 0;
    Block bb = *(Block *) b;
    sigB = bb.signature;

    //Return comparison
    return sigA - sigB;
}

/*
 findCollisionsOptimised

 input blockDatabase
 Finds all collisions between generated blocks and return collision database; also store number of collisions found. Using sorting method instead of brute force
 */
Collision *findCollisionsOptimised(Block *blockDB, int numBlocks, int *numberCollisionsFound) {
    //Sort block database
    qsort(blockDB, numBlocks, sizeof(Block), cmpfunc);

    //Set up collision database
    Collision *collisions = (Collision *) malloc(1 * sizeof(Collision));

    //Linearly loop through block database, storing collisions as they're found
    int numCollisions = 0;
    long previousSig = blockDB[0].signature;
    int curBlocksInCollision = 0;
    for(int i=1; i<numBlocks; i++) {
        if(blockDB[i].signature == previousSig) {
            //Collision detected
            if(curBlocksInCollision == 0) {
                //New collision. Reallocate more memory for collision database
                collisions = (Collision *) realloc(collisions, sizeof(Collision) * (numCollisions + 1));
                //Increment counter
                numCollisions++;
                //Store first two blocks of current collision
                curBlocksInCollision = 2;
                collisions[numCollisions-1].numBlocksInCollision = 2;
                collisions[numCollisions-1].columns = (int *) malloc(sizeof(int) * 5);
                collisions[numCollisions-1].blocks = (Block *) malloc(sizeof(Block) * curBlocksInCollision);
                collisions[numCollisions-1].columns[0] = blockDB[i-1].column;
                collisions[numCollisions-1].columns[1] = blockDB[i].column;
                collisions[numCollisions-1].blocks[0] = blockDB[i-1];
                collisions[numCollisions-1].blocks[1] = blockDB[i];
                //printf("%d: Found collision on signature %ld with %d blocks in it\n", numCollisions-1, previousSig, curBlocksInCollision);
            }
            else {
                //Increment counter
                curBlocksInCollision++;
                collisions[numCollisions-1].numBlocksInCollision = curBlocksInCollision;
                //Store current block's column number, and rellocate collisions column array
                collisions[numCollisions-1].columns = (int *) realloc(collisions[numCollisions-1].columns, sizeof(int) * curBlocksInCollision);
                collisions[numCollisions-1].columns[curBlocksInCollision-1] = blockDB[i].column;
                //Store current block, and reallocate collision's block array
                collisions[numCollisions-1].blocks = (Block *) realloc(collisions[numCollisions-1].blocks, sizeof(Block) * curBlocksInCollision);
                collisions[numCollisions-1].blocks[curBlocksInCollision-1] = blockDB[i];
                //printf("%d: Found collision on signature %ld with %d blocks in it\n", numCollisions-1, previousSig, curBlocksInCollision);
            }
        }
        else {
            //No collision between previous block and current block. Record new signature
            previousSig = blockDB[i].signature;
            curBlocksInCollision = 0;
        }
    }

    //Record number of collisions found
    *numberCollisionsFound = numCollisions;

    //Return collision database
    return collisions;
}

/*
 findCollisionsOptimisedParallel

 input blockDatabase
 Finds all collisions between generated blocks and return collision database; also store number of collisions found. Using sorting method instead of brute force. Parallelized using tasks
 */
Collision *findCollisionsOptimisedParallel(Block *blockDB, int numBlocks, int *numberCollisionsFound) {
    //Sort block database
    qsort(blockDB, numBlocks, sizeof(Block), cmpfunc);

    //Set up collision database
    Collision *collectiveCollisionDB = (Collision *) malloc(1 * sizeof(Collision));
    int numTotalCollisions = 0;

    //Test time spent in critical region
    double timeInCritical = 0;

    //Begin parallel region
    #pragma omp parallel
    {
        #pragma omp single
        {
            //Set up index values for thread work
            int nthreads = omp_get_num_threads();
            int start = 0;
            int increment = numBlocks/nthreads;
            int finish = start + increment;
            //Move finishing index value to next collision (i.e. next signature)
            while(blockDB[finish].signature == blockDB[finish-1].signature) {
                finish++;
            }

            //Generate tasks to find all collisions
            while(start < numBlocks) {
                //Now create task to find all collisions between start and finish
                #pragma omp task
                {
                    //Get thread number
                    //int ID = omp_get_thread_num();

                    //Set up partial collision database
                    Collision *collisions = (Collision *) malloc(1 * sizeof(Collision));

                    //Linearly loop through block database, storing collisions as they're found
                    int numCollisions = 0;
                    long previousSig = blockDB[start].signature;
                    int curBlocksInCollision = 0;
                    for(int i=start+1; i<finish; i++) {
                        if(blockDB[i].signature == previousSig) {
                            //Collision detected
                            if(curBlocksInCollision == 0) {
                                //New collision. Reallocate more memory for collision database
                                collisions = (Collision *) realloc(collisions, sizeof(Collision) * (numCollisions + 1));
                                //Increment counter
                                numCollisions++;
                                //Store first two blocks of current collision
                                curBlocksInCollision = 2;
                                collisions[numCollisions-1].numBlocksInCollision = 2;
                                collisions[numCollisions-1].columns = (int *) malloc(sizeof(int) * 5);
                                collisions[numCollisions-1].columns[0] = blockDB[i-1].column;
                                collisions[numCollisions-1].columns[1] = blockDB[i].column;
                                collisions[numCollisions-1].blocks = (Block *) malloc(sizeof(Block) * curBlocksInCollision);
                                collisions[numCollisions-1].blocks[0] = blockDB[i-1];
                                collisions[numCollisions-1].blocks[1]= blockDB[i];
                                //printf("Thread %d: Found collision %d on signature %ld with %d blocks in it\n", ID, numCollisions-1, previousSig, curBlocksInCollision);
                            }
                            else {
                                //Increment counter
                                curBlocksInCollision++;
                                collisions[numCollisions-1].numBlocksInCollision = curBlocksInCollision;
                                //Store current block's column number, and rellocate collisions column array
                                collisions[numCollisions-1].columns = (int *) realloc(collisions[numCollisions-1].columns, sizeof(int) * curBlocksInCollision);
                                collisions[numCollisions-1].columns[curBlocksInCollision-1] = blockDB[i].column;
                                //Store current block, and reallocte collision's block array
                                collisions[numCollisions-1].blocks = (Block *) realloc(collisions[numCollisions-1].blocks, sizeof(Block) * curBlocksInCollision);
                                collisions[numCollisions-1].blocks[curBlocksInCollision-1] = blockDB[i];
                                //printf("Thread %d: Found collision %d on signature %ld with %d blocks in it\n", ID, numCollisions-1, previousSig, curBlocksInCollision);
                            }
                        }
                        else {
                            //No collision between previous block and current block. Record new signature
                            previousSig = blockDB[i].signature;
                            curBlocksInCollision = 0;
                        }
                    }
                    //Merge partial collision database with complete collision database
                    //double start = omp_get_wtime();
                    #pragma omp critical
                    {
                        collectiveCollisionDB = mergeCollisionDatabases(collectiveCollisionDB, collisions, &numTotalCollisions, numCollisions);
                    }
                    double duration = omp_get_wtime() - start;
                    #pragma omp atomic
                    timeInCritical += duration;
                    //Free partial collision database (column databases will be freed in freeCollisionDB)
                    free(collisions);

                }
                //Update indexes to create next task
                start = finish;
                finish += increment;
                if(finish > numBlocks) {
                    finish = numBlocks;
                }
                else {
                    //Move finishing index value to next collision (i.e. next signature)
                    while(blockDB[finish].signature == blockDB[finish-1].signature) {
                        finish++;
                    }
                }
            }
        }
    }

    //Print amount of time spent in critical region
//    printf("%5lf milli-seconds spent in critical region\n", (double) timeInCritical * 1000);

    //Record number of collisions found
    *numberCollisionsFound = numTotalCollisions;

    //Return collision database
    return collectiveCollisionDB;
}

/*
 getRowSet

 Input two Collision structures
 Output int representing the number of unique row numbers found between the two collisions
 */
int *getRowSet(Collision c1, Collision c2, int *numFound) {

    //Array of rows to return
    int *rowSet = (int *) malloc(sizeof(int));

    //Count number of unique values found
    int numUnique = 0;

    //Record whether or not a row has already been seen
    bool found[ROWS];
    for(int i=0; i<ROWS; i++) {
        found[i] = false;
    }

    //Loop through first collision's rows
    for(int i=0; i<4; i++) {
        //Get current row number
        int curRow = c1.blocks[0].rows[i];
        //Check if it's been seen
        if(found[curRow]) {
            continue;
        }
        //If it hasn't, add it to the set
        rowSet[numUnique] = curRow;
        //Reallocate rowSet
        numUnique++;
        rowSet = (int *) realloc(rowSet, (numUnique+1) * sizeof(int));
    }

    //Loop through second collision's rows
    for(int i=0; i<4; i++) {
        //Get current row number
        int curRow = c2.blocks[0].rows[i];
        //Check if it's been seen
        if(found[curRow]) {
            continue;
        }
        //If it hasn't, add it to the set
        rowSet[numUnique] = curRow;
        //Reallocate rowset
        numUnique++;
        rowSet = (int *) realloc(rowSet, (numUnique+1) * sizeof(int));
    }

    //Record number found
    *numFound = numUnique;

    //Return rowset
    return rowSet;
}


/*
 mergeCollisions

 Input collision database, number of collisions, and pointer for number of merged collisions found
 Merge all collisions over same column set, and return array of MergedCollision structures
 //NEEDS WORK
 */
MergedCollision *mergeCollisions(Collision *collisions, int numCollisions, int *numMerged) {
    //Create database of merged collisions
    int numFound = 0;
    MergedCollision *merged = (MergedCollision *) malloc(sizeof(MergedCollision));

    //Loop through all collisions in collision database
    for(int i=0; i<numCollisions; i++) {
        for(int j=i+1; j<numCollisions; j++) {
            //First check if they have the same number of blocks inside them (if not, the column sets can't be equal)
            if(collisions[i].numBlocksInCollision != collisions[j].numBlocksInCollision) {
                continue;
            }
            bool columnSetSame = true;
            //If they are the same, check that all columns in the first are in the second
            for(int k=0; k<collisions[i].numBlocksInCollision; k++) {
                //Record current column number of first collision, and whether or not it has been found in second collision
                int curColumn = collisions[i].columns[k];
                bool foundInSecond = false;
                //Now check against all column numbers in second collision
                for(int k1=0; k1<collisions[j].numBlocksInCollision; k1++) {
                    if(curColumn == collisions[j].columns[k1]) {
                        foundInSecond = true;
                        break;
                    }
                }
                //Now check if the column number was found in the second collision
                if(!foundInSecond) {
                    //These collisions don't have the same column set. Continue to next iteration of j's loop
                    columnSetSame = false;
                    break;
                }
            }
            if(!columnSetSame) {
                //Go to next iteration of j's loop
                continue;
            }
            //If we've made it to this point, the collisions have the same column set
            merged = (MergedCollision *) realloc(merged, sizeof(MergedCollision) * (numFound+1));
            //Store columns in merged collision database
            merged[numFound].numColumns = collisions[j].numBlocksInCollision;
            merged[numFound].columns = (int *) malloc(sizeof(int) * merged[numFound].numColumns);
            for(int k=0; k<collisions[j].numBlocksInCollision; k++) {
                merged[numFound].columns[k] = collisions[j].columns[k];
            }
            //Store rows in merged collision
            int numUniqueRows;
            merged[numFound].rows = getRowSet(collisions[i], collisions[j], &numUniqueRows);
            merged[numFound].numRows = numUniqueRows;
            //Increment numFound
            numFound++;
        }
    }
    //Record number of merged collisions found
    *numMerged = numFound;
    //Return merged collision database
    return merged;
}
