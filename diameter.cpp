/*
COMP90025 Project 1D: MPI and Diameter
San Kho Lin (829463) sanl1@student.unimelb.edu.au

Tested on VLSCI Snowy:
    mpicxx diameter.cpp -o diameter.exe

Tested with Intel Compiler on Windows:
    mpicxx.bat mandelbrot.cpp

Limitations:
    Due to Row-wise decomposition approach use, this program requires such that:
    number of vertices must be divisible by the number of processors.
    This can be further improved using some heuristic approaches or
    using of Block-wise decomposition or two-dimensional decomposition of the various matrices.

    REF:
    - Ian Foster. Designing and Building Parallel Programs. 1995. http://www.mcs.anl.gov/~itf/dbpp/text/node35.html
    - Quinn Michael J. Parallel Programming in C with MPI and OpenMP. McGraw Hill, 2003.
*/

#include <cstdio>
#include <cstdlib>
#include <mpi.h>

//using namespace std;

#define MAX 10000
#define NOT_CONNECTED -1

int distance[MAX][MAX];

//number of nodes
int nodesCount;
int nprocs;
int rank;

double start_time;

//initialize all distances to 
void Initialize() {
    for (int i = 0; i<MAX; ++i) {
        for (int j = 0; j<MAX; ++j) {
            distance[i][j] = NOT_CONNECTED;
            //printf("%d %d\n", i, j);
        }
        distance[i][i] = 0;
    }
}

void printArray(int *row, int nElements) {
    int i;
    for (i = 0; i<nElements; i++) {
        printf("%d ", row[i]);
    }
    printf("\n");
}

void printMatrix(int mat[], int n) {
    int i, j;
    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            if (mat[i * n + j] == NOT_CONNECTED)
                printf("%d ", NOT_CONNECTED);
            else
                printf("%d ", mat[i * n + j]);
        }
        printf("\n");
    }
}

int main(int argc, char *argv[]) {

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    if (rank == 0) {
        
        //start_time = MPI_Wtime(); // time measure        

        scanf("%d", &nodesCount);
    }
    MPI_Bcast(&nodesCount, 1, MPI_INT, 0, MPI_COMM_WORLD);
    //printf("rank: %d, nodesCount: %d\n", rank, nodesCount);

    //-- Row-wise Decomposition Limitations

    if (nprocs > nodesCount) {
        if (rank == 0) {
            printf("Number of processors is higher than number of graph vertices/nodes. Re-run with: %d processors\n", nodesCount);
        }
        MPI_Finalize();
        exit(1);
    }

    if (nodesCount % nprocs != 0) {
        if (rank == 0)
            printf("Number of nodes must be divisible by number of processors.\n");
        MPI_Finalize();
        exit(1);
    }

    int nxn_nodesCount = nodesCount * nodesCount; // N x N, 1D nodes count

    //--

    int *matrix;

    if (rank == 0) {

        if ((matrix = (int *)malloc(nxn_nodesCount * sizeof(int))) == NULL) {
            printf("Malloc error");
            exit(1);
        }

        Initialize();

        //edges count
        int m;
        scanf("%d", &m);

        while (m--) {
            //nodes - let the indexation begin from 1
            int a, b;

            //edge weight
            int c;

            scanf("%d-%d-%d", &a, &c, &b);
            distance[a][b] = c;
        }

        // Flatten 2D stack array into contiguous 1D heap array
        for (int n = 1; n <= nodesCount; n++) {
            for (int m = 1; m <= nodesCount; m++) {
                //printf("%d\t", distance[n][m]);
                matrix[(n-1)*nodesCount + (m-1)] = distance[n][m];
            }
            //printf("\n");
        }

        //printf("1D Adjacency Matrix: \t");
        //printArray(matrix, nxn_nodesCount);
        //printf("\n");
        //printf("\n");
    }

    //--

    int mychunk = nodesCount * (nodesCount / nprocs);

    //printf("rank: %d, mychunk: %d\n", rank, mychunk);
    //exit(0);

    int* local_mat = (int *)malloc(mychunk * sizeof(int));
    if (local_mat == NULL) {
        perror("Error in malloc 3");
        exit(1);
    }

    if (MPI_Scatter(matrix, mychunk, MPI_INT, local_mat, mychunk, MPI_INT, 0, MPI_COMM_WORLD) != MPI_SUCCESS) {
        perror("Scatter error");
        exit(1);
    }

    //printf("rank [%d], received [%d]elements: ", rank, mychunk);
    //printArray(local_mat, mychunk);
    //MPI_Barrier(MPI_COMM_WORLD);
    //exit(0);

    //--

    int* bcast_row;
    int offset;
    int root;

    bcast_row = (int *) malloc(nodesCount * sizeof(int));

    int k;
    for (k = 0; k < nodesCount; k++) {

        root = k / (nodesCount / nprocs); // block owner

        if (rank == root) {

            offset = k % (nodesCount / nprocs);
            //printf("rank: %d, root: %d, offset: %d, bcast_row: ", rank, root, offset);

            int h;
            for (h = 0; h < nodesCount; h++) {
                bcast_row[h] = local_mat[offset * nodesCount + h];
                //printf(" %d", bcast_row[h]);
            }
            //printf("\n");
        }

        MPI_Bcast(bcast_row, nodesCount, MPI_INT, root, MPI_COMM_WORLD);
        //MPI_Barrier(MPI_COMM_WORLD);

        int i;
        for (i = 0; i < nodesCount / nprocs; i++) {
            int idx_ik = i * nodesCount + k;
            if (local_mat[idx_ik] != NOT_CONNECTED) {
                int j;
                for (j = 0; j < nodesCount; j++) {
                    int idx_ij = i * nodesCount + j;

                    if (bcast_row[j] != NOT_CONNECTED && (local_mat[idx_ij] == NOT_CONNECTED || 
                        local_mat[idx_ik] + bcast_row[j] < local_mat[idx_ij])) {

                        local_mat[idx_ij] = local_mat[idx_ik] + bcast_row[j];
                    }
                }
                //printf("\n");
            }
        }
    }
    free(bcast_row);

    //-- Print local matrix
    //printf("rank: %d\n", rank);
    //printMatrix(local_mat, nodesCount);
    //printArray(local_mat, nxn_nodesCount);
    //printf("\n");
    //printf("\n");

    //--

    int* final_mat = (int *) malloc(nxn_nodesCount * sizeof(int));

    MPI_Gather(local_mat, mychunk, MPI_INT, final_mat, mychunk, MPI_INT, 0, MPI_COMM_WORLD);

    //-- Finale

    if (rank == 0) {
        //printf("Final matrix is:\n");
        //printf("2D\n");
        //printMatrix(final_mat, nodesCount);
        //printf("1D\n");
        //printArray(final_mat, nxn_nodesCount);
        
        int biggest = -1;
        for (int i = 0; i < nxn_nodesCount; i++) {
            if (biggest < final_mat[i])
                biggest = final_mat[i];
        }
        printf("%d\n", biggest);

        //--
        // Please uncomment the follow block if you want
        // to see the all most distant pairs output like in sequential code

        /*
        int longest_path = -1;
        //reconstruct the distance stack 2d array
        for (int i = 0; i < nodesCount; ++i) {
            for (int j = 0; j < nodesCount; ++j) {
                
                int vertex = final_mat[i * nodesCount + j];

                if (vertex == NOT_CONNECTED) {
                    distance[i + 1][j + 1] = NOT_CONNECTED;
                }
                else {
                    if (longest_path < vertex) {
                        longest_path = vertex;
                    }
                    distance[i + 1][j + 1] = vertex;
                }
            }
        }
        //printf("%d\n", longest_path);

        int diameter = -1;
        //look for the most distant pair
        for (int i = 1; i <= nodesCount; ++i) {
            for (int j = 1; j <= nodesCount; ++j) {
                if (diameter<distance[i][j]) {
                    diameter = distance[i][j];
                    printf("%d-%d-%d\n", i, diameter, j);
                }
            }
        }
        printf("%d\n", diameter);
        */

        //--
        /*
         Project 1 B output
         Posted on: Tuesday, 20 September 2016 11:02:27 AM EST
         Posted by: Aaron Harwood
         For the program output, only the diameter of the graph need to be output.
         The example program outputs all most distant pairs, but this line should have been commented out.
        */

        //printf("\nElapsed time(second): %f\n", (float)(MPI_Wtime() - start_time)); // time measure

    } // rank 0 end
    
    
    MPI_Finalize();

    return 0;
}
