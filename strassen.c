/*
Project 2 - Strassen’s matrix multiplication algorithm.​
By Aravind Venkit
*/
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <time.h>
#include <math.h>

/*
Displays the given matrix.
*/
void DisplayMatrix(int N, float matrix[][N], int M){
	int i,j;
	printf("\n");
	for(i = 0; i < M; ++i)
	{
	   for(j = 0; j < M ; ++j)
	   { 
	      printf("%7.2f", matrix[i][j]); 
	   }
	   printf("\n");
	}
}

/*
Generate Matrix with random number -5.0 t0 5.0.
*/
void GenerateMatrix(int N, float matrix[][N], int M){
  memset(matrix, 0, sizeof(matrix[0][0]) * N * N);
  int i, j;
  int range = 5;
  srand((int)time(NULL));
  for(i = 0; i < M; ++i)
	{
	   for(j = 0; j < M ; ++j)
	   {        
	   		float newElement = (-1+2*(float)rand()/(float)(RAND_MAX)) * range;
	      	matrix[i][j] = newElement;            
	   }
	}
	DisplayMatrix(N, matrix, M);
}

/*
StandardMultiplication to multiply two matrices
*/
void StandardMultiplication(int N, float matrixA[][N], float matrixB[][N], float result[][N]){
	int i,j,k;
	for (i = 0; i < N; i++) {
        for(j = 0; j < N; j++) {
        	result[i][j] = 0;
            for(k = 0; k < N; k++) {
               result[i][j] += matrixA[i][k] * matrixB[k][j];
            }
        } 
    }
   
}

/*
AddMatrix to add to two matrices.
*/
void AddMatrix(int N, float matrixA[][N], float matrixB[][N], float result[][N]){
	int i,j,k;
	for (i = 0; i < N; i++) {
        for(j = 0; j < N; j++) {
        	result[i][j] = 0;
            result[i][j] = matrixA[i][j] + matrixB[i][j];
        } 
    }
}

/*
SubMatrix to subtract to two matrices.
*/
void SubMatrix(int N, float matrixA[][N], float matrixB[][N], float result[][N]){
	int i,j,k;
	for (i = 0; i < N; i++) {
        for(j = 0; j < N; j++) {
        	result[i][j] = 0;
            result[i][j] = matrixA[i][j] - matrixB[i][j];
        } 
    }
}

/*
StrassensBase - Implements Strassen’s matrix multiplication algorithm for 2*2 Matrix.
*/
void StrassensBase(float basematrixA[][2], float basematrixB[][2], float baseResult[][2]){

    float m1 = (basematrixA[0][0] + basematrixA[1][1]) * (basematrixB[0][0] + basematrixB[1][1]);
    float m2 = (basematrixA[1][0] + basematrixA[1][1]) * basematrixB[0][0];
    float m3 = basematrixA[0][0] * (basematrixB[0][1] - basematrixB[1][1]);
    float m4 = basematrixA[1][1] * (basematrixB[1][0] - basematrixB[0][0]);
    float m5 = (basematrixA[0][0] + basematrixA[0][1]) * basematrixB[1][1];
    float m6 = (basematrixA[1][0] - basematrixA[0][0]) * (basematrixB[0][0] + basematrixB[0][1]);
    float m7 = (basematrixA[0][1] - basematrixA[1][1]) * (basematrixB[1][0] + basematrixB[1][1]);
    
    baseResult[0][0] = m1 + m4 - m5 + m7;
    baseResult[0][1] = m3 + m5;
    baseResult[1][0] = m2 + m4;
    baseResult[1][1] = m1 + m3 - m2 + m6;
}

/*
StrassensMultiplication - Implements Strassen’s matrix multiplication algorithm in recursive manner.
*/
void strassensMultiplication (int N, float matrixA[][N], float matrixB[][N], float result[][N]){
     memset(result, 0, sizeof(result[0][0]) * N * N);
    if (N == 2){
        StrassensBase(matrixA, matrixB, result);
    }
    else if (N == 1) {
        result[0][0] = (float)matrixA[0][0] * (float)matrixB[0][0];
    }
    else {
        int M = N/2;

        float matrixA11[M][M], matrixA12[M][M], matrixA21[M][M], matrixA22[M][M];
        float matrixB11[M][M], matrixB12[M][M], matrixB21[M][M], matrixB22[M][M];
        float result11[M][M], result12[M][M], result21[M][M], result22[M][M];
        float m1[M][M], m2[M][M], m3[M][M], m4[M][M], m5[M][M], m6[M][M], m7[M][M];
        float aResult[M][M], bResult[M][M];

        //dividing the matrices in 4 sub-matrices
        int i, j;
        for (i = 0; i < M; i++) {
            for (j = 0; j < M; j++) {
                matrixA11[i][j] = matrixA[i][j];
                matrixA12[i][j] = matrixA[i][j + M];
                matrixA21[i][j] = matrixA[i + M][j];
                matrixA22[i][j] = matrixA[i + M][j + M];
                matrixB11[i][j] = matrixB[i][j];
                matrixB12[i][j] = matrixB[i][j + M];
                matrixB21[i][j] = matrixB[i + M][j];
                matrixB22[i][j] = matrixB[i + M][j + M];
            }
        }

        // Calculating m1 to m7:
        AddMatrix(M, matrixA11, matrixA22, aResult); // matrixA11 + matrixA22
        AddMatrix(M, matrixB11, matrixB22, bResult); // matrixB11 + matrixB22
        strassensMultiplication(M, aResult, bResult, m1); // m1 = (matrixA11+matrixA22) * (matrixB11+matrixB22)

        AddMatrix(M, matrixA21, matrixA22, aResult); // matrixA21 + matrixA22
        strassensMultiplication(M, aResult, matrixB11, m2); // p2 = (matrixA21+matrixA22) * (matrixB11)

        SubMatrix(M, matrixB12, matrixB22, bResult); // matrixB12 - matrixB22
        strassensMultiplication(M, matrixA11, bResult, m3); // m3 = (matrixA11) * (matrixB12 - matrixB22)

        SubMatrix(M, matrixB21, matrixB11, bResult); // matrixB21 - matrixB11
        strassensMultiplication(M, matrixA22, bResult, m4); // m4 = (matrixA22) * (matrixB21 - matrixB11)

        AddMatrix(M, matrixA11, matrixA12, aResult); // matrixA11 + matrixA12
        strassensMultiplication(M, aResult, matrixB22, m5); // m5 = (matrixA11+matrixA12) * (matrixB22) 

        SubMatrix(M, matrixA21, matrixA11, aResult); // matrixA21 - matrixA11
        AddMatrix(M, matrixB11, matrixB12, bResult); // matrixB11 + matrixB12
        strassensMultiplication(M, aResult, bResult, m6); // m6 = (matrixA21-matrixA11) * (matrixB11+matrixB12)

        SubMatrix(M, matrixA12, matrixA22, aResult); // matrixA12 - matrixA22
        AddMatrix(M, matrixB21, matrixB22, bResult); // matrixB21 + matrixB22
        strassensMultiplication(M, aResult, bResult, m7); // m7 = (matrixA12-a22) * (matrixB21+matrixB22)

        // calculating result21, result21, result11, result22:

        AddMatrix(M, m3, m5, result12); // result12 = m3 + m5
        AddMatrix(M, m2, m4, result21); // result21 = m2 + m4

        AddMatrix(M, m1, m4, aResult); // m1 + m4
        AddMatrix(M, aResult, m7, bResult); // m1 + m4 + m7
        SubMatrix(M, bResult, m5, result11); // result11 = m1 + m4 - m5 + m7

        AddMatrix(M, m1, m3, aResult); // m1 + m3
        AddMatrix(M, aResult, m6, bResult); // m1 + m3 + m6
        SubMatrix(M, bResult, m2, result22); // result22 = m1 + m3 - m2 + m6

        // Grouping the results obtained in a single matrix
        for (i = 0; i < M ; i++) {
            for (j = 0 ; j < M ; j++) {
                result[i][j] = result11[i][j];
                result[i][j + M] = result12[i][j];
                result[i + M][j] = result21[i][j];
                result[i + M][j + M] = result22[i][j];
            }
        }

    }
}

/*
main - Implements main methode to call other methodes as given in projec 2 document.
*/
int main(int argc, char *argv[])
{   
	int M = 0, N = 0;
    if(argv[1] == NULL){
	   printf("\nEnter n*n Vertex Value for both A and B matrix: \n");
	   scanf("%d", &M);
    }
    else{
       M = atoi(argv[1]);
    }

    N = pow(2, ceil(log(M)/log(2)));

	printf("\nMatrix A: ");
	float matrixA[N][N];
	GenerateMatrix(N, matrixA, M);

	printf("\nMatrix B: ");
	float matrixB[N][N];
	GenerateMatrix(N, matrixB, M);

	float strassenResult[N][N];//Stressen Multiplication
    memset(strassenResult, 0, sizeof(strassenResult[0][0]) * N * N);
	printf("\nStrassen’s Multiplication Output: ");
	strassensMultiplication(N, matrixA, matrixB, strassenResult);
	DisplayMatrix(N, strassenResult, M);

    float standardResult[N][N];//Standard Multiplication
    memset(standardResult, 0, sizeof(standardResult[0][0]) * N * N);
	printf("\nStandard Multiplication Output: ");
	StandardMultiplication(N, matrixA, matrixB, standardResult);
    DisplayMatrix(N, standardResult, M);

	printf("\n");
	return 0;
}
