/*
Project 2 - Tromino’s matrix multiplication algorithm.​
By Aravind Venkit
*/
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <math.h>

int cnt = 0;
/*
Displays the given matrix.
*/
void DisplayMatrix(int N, int matrix[][N]){
    int i,j;
    printf("\n");
    for(i = 0; i < N; ++i)
    {
       for(j = 0; j < N ; ++j)
       { 
         if(matrix[i][j] == -1){
            printf("%4c", 'X'); 
         }
         else{
            printf("%4d", matrix[i][j]); 
         }          
       }
       printf("\n");
    }
}
/*
TrominoTiling algorithm in recursive manner.
*/
void trominoTile(int N, int matrix[][N],int size, int start_R, int start_C, int hole_R, int hole_C){ 
    cnt++;  
    int end_R = (start_R + size)-1, end_C = (start_C + size)-1;

    if(size == 2){
        if(matrix[start_R][start_C] == 0){
            matrix[start_R][start_C] = cnt;
        }
        if(matrix[start_R][end_C] == 0){
        matrix[start_R][end_C] = cnt;
        }
        if(matrix[end_R][start_C] == 0){
            matrix[end_R][start_C] = cnt;
        }
        if(matrix[end_R][end_C] == 0){
            matrix[end_R][end_C] = cnt;
        }
    }
    else{
        int Mid_R = (size/2) + start_R;
        int Mid_C = (size/2) + start_C;

        int holeIn = 0;
        int hole_R11 = -1, hole_C11 = -1, hole_R12 = -1, hole_C12 = -1, hole_R21 = -1, 
            hole_C21 = -1, hole_R22 = -1, hole_C22 = -1;

        /*
        Finds wether the hole is in 11 block of the given matrix
        */
        if(start_R <= hole_R && start_C <= hole_C && hole_R <= Mid_R-1 && hole_C <= Mid_C-1){
            holeIn = 11;
        }

        /*
        Finds wether the hole is in 12 block of the given matrix
        */         
        if(start_R <= hole_R && Mid_C <= hole_C && hole_R <= Mid_R-1 && hole_C <= end_C){
            holeIn = 12;
        }

        /*
        Finds wether the hole is in 21 block of the given matrix
        */
        if(Mid_R <= hole_R && start_C <= hole_C && hole_R <= end_R && hole_C <= Mid_C-1){
            holeIn = 21;
        }

        /*
        Finds wether the hole is in 22 block of the given matrix
        */
        if(Mid_R <= hole_R && Mid_C <= hole_C && hole_R <= end_C && hole_R <= end_C){
            holeIn = 22;
        }

        /*
        Mark single tile of each 11,12,21 and 22 blocks other than the tail which has a hole i.e not 0
        */
        if(holeIn != 11){
            matrix[Mid_R-1][Mid_C-1] = cnt;
            hole_R11 = Mid_R-1;
            hole_C11 = Mid_C-1;
        }
        else{
            hole_R11 = hole_R;
            hole_C11 = hole_C;
        }

        if(holeIn != 12){
            matrix[Mid_R-1][Mid_C] = cnt;
            hole_R12 = Mid_R-1;
            hole_C12 = Mid_C;
        }
        else{
            hole_R12 = hole_R;
            hole_C12 = hole_C;
        }

        if(holeIn != 21){
            matrix[Mid_R][Mid_C-1] = cnt;
            hole_R21 = Mid_R;
            hole_C21 = Mid_C-1;
        }

        else{
            hole_R21 = hole_R;
            hole_C21 = hole_C;
        }

        if(holeIn != 22){
            matrix[Mid_R][Mid_C] = cnt;
            hole_R22 = Mid_R;
            hole_C22 = Mid_C;
        }
        else{
            hole_R22 = hole_R;

            hole_C22 = hole_C;
        }

        /*
        TrominoTiling function recursive in nature for each reapective 11, 12, 21, and 22 blocks.
        */
        trominoTile(N, matrix, size/2, start_R, start_C, hole_R11, hole_C11);
        trominoTile(N, matrix, size/2, start_R, Mid_C, hole_R12, hole_C12);
        trominoTile(N, matrix, size/2, Mid_R, start_C, hole_R21, hole_C21);
        trominoTile(N, matrix, size/2, Mid_R, Mid_C, hole_R22, hole_C22); 
    }
}
/*
main function to get the required inputs and calls the trimino tiling function and displays the result
*/
int main(int argc, char *argv[])
{   
    /*
    To get required inputs.
    */
    int K= 0, N = 0, hole_R = 0, hole_C = 0;
    
    if(argc < 4){
        printf("Arguments Missing for K, Row(hole), Collumn(hole) values\n");
        return 0;
    }

    K = atoi(argv[1]);
    hole_R = atoi(argv[2]);
    hole_C = atoi(argv[3]);
 
    N = pow(2, K) ;// calculate and set N to 2^K
    /*
    To generate matrix with all zero except for hole -1(X will be printed)
    */
	int matrix[N][N];
    memset(matrix, 0, sizeof(matrix[0][0]) * N * N);
    matrix[hole_R][hole_C] = -1;  

    /*
    trominoTile function is called
    */
	trominoTile(N, matrix,N, 0, 0, hole_R, hole_C);	

    /*
    DisplayMatrix function to display the result.
    */
    DisplayMatrix(N, matrix);

	printf("\n");
	return 0;
}

