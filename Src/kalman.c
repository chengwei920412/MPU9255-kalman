#include "kalman.h"
#include "math.h"
#include "system.h"
#include "mpu9255.h"

extern float AngleGyro;


/*float CalculateAngleAccel(float AccelX,float AccelZ)
{
	float AngleAccel;
	AngleAccel = atan2(AccelX ,AccelZ);
	return AngleAccel;
}*/
float CalculateAngleAccel(float AccelX,float AccelY,float AccelZ)
{
	float temp;
	float res=0;
	temp = AccelX/sqrt((AccelY*AccelY+AccelZ*AccelZ)) ;
	res=atan(temp);
	return res*180/PI;
}
float CalculateAngleRateGyro(int16_t GyroX)
{
	float AngleRateGyro;
	AngleRateGyro = - GyroX ;//* PI/ 180.00;
	AngleGyro = AngleGyro + AngleRateGyro * dt;
	return AngleGyro;
}

void MatrixAddition(float* A, float* B, int m, int n, float* C)
{
//the matrix addition function
// A = input matrix (m x n)
// B = input matrix (m x n)
// m = number of rows in A = number of rows in B
// n = number of columns in A = number of columns in B
// C = output matrix = A+B (m x n)
	int i, j;
	for (i=0;i<m;i++)
		for(j=0;j<n;j++)
			C[n*i+j]=A[n*i+j]+B[n*i+j];
}

void MatrixSubtraction(float* A, float* B, int m, int n, float* C)
{
//the matrix subtraction function
// A = input matrix (m x n)
// B = input matrix (m x n)
// m = number of rows in A = number of rows in B
// n = number of columns in A = number of columns in B
// C = output matrix = A-B (m x n)
	int i, j;
	for (i=0;i<m;i++)
		for(j=0;j<n;j++)
			C[n*i+j]=A[n*i+j]-B[n*i+j];
}

void MatrixMultiply(float* A, float* B, int m, int p, int n, float* C)
{
//the matrix multipliction function
// A = input matrix (m x p)
// B = input matrix (p x n)
// m = number of rows in A
// p = number of columns in A = number of rows in B
// n = number of columns in B
// C = output matrix = A*B (m x n)
	int i, j, k;
	for (i=0;i<m;i++)
	for(j=0;j<n;j++)
	{
		C[n*i+j]=0;
		for (k=0;k<p;k++)
			C[n*i+j]= C[n*i+j]+A[p*i+k]*B[n*k+j];
	}
}

int MatrixInversion(float* A, int n, float* AInverse)
{
//the matrix inversion function
// A = input matrix (n x n)
// n = dimension of A
// AInverse = inverted matrix (n x n)
// This function inverts a matrix based on the Gauss Jordan method.
// The function returns 1 on success, 0 on failure.
	int i, j, iPass, imx, icol, irow;
	float det, temp, pivot, factor;
	float* ac = (float*)calloc(n*n, sizeof(float));
	det = 1;
	for (i = 0; i < n; i++)
	{
		for (j = 0; j < n; j++)
		{
			AInverse[n*i+j] = 0;
			ac[n*i+j] = A[n*i+j];
		}
		AInverse[n*i+i] = 1;
	}
	// The current pivot row is iPass.
	// For each pass, first find the maximum element in the pivot column.
	for (iPass = 0; iPass < n; iPass++)
	{
		imx = iPass;
		for (irow = iPass; irow < n; irow++)
		{
			if (fabs(A[n*irow+iPass]) > fabs(A[n*imx+iPass])) imx = irow;
		}
		// Interchange the elements of row iPass and row imx in both A and AInverse.
		if (imx != iPass)
		{
			for (icol = 0; icol < n; icol++)
			{
				temp = AInverse[n*iPass+icol];
				AInverse[n*iPass+icol] = AInverse[n*imx+icol];
				AInverse[n*imx+icol] = temp;
				if (icol >= iPass)
				{
					temp = A[n*iPass+icol];
					A[n*iPass+icol] = A[n*imx+icol];
					A[n*imx+icol] = temp;
				}
			}
		}
		// The current pivot is now A[iPass][iPass].
		// The determinant is the product of the pivot elements.
		pivot = A[n*iPass+iPass];
		det = det * pivot;
		if (det == 0)
		{
			free(ac);
			return 0;
		}
		for (icol = 0; icol < n; icol++)
		{
			// Normalize the pivot row by dividing by the pivot element.
			AInverse[n*iPass+icol] = AInverse[n*iPass+icol] / pivot;
			if (icol >= iPass) A[n*iPass+icol] = A[n*iPass+icol] / pivot;
		}
		for (irow = 0; irow < n; irow++)
		// Add a multiple of the pivot row to each row. The multiple factor
		// is chosen so that the element of A on the pivot column is 0.
		{
			if (irow != iPass) factor = A[n*irow+iPass];
			for (icol = 0; icol < n; icol++)
			{
				if (irow != iPass)
				{
					AInverse[n*irow+icol] -= factor * AInverse[n*iPass+icol];
					A[n*irow+icol] -= factor * A[n*iPass+icol];
				}
			}
		}
	}
	free(ac);
	return 1;
}
