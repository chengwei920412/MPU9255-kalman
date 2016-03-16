#ifndef __kalman_h
#define __kalman_h

#include "stdint.h"

#define PI 3.14159
#define dt 0.00035

#define GyroStanDev 2
#define BiasStanDev 2
#define GyroVar 4
#define AccelVar 4
#define p1 0.1
#define q1 0.1
#define m1 0.1
#define n1 0.1
#define e 0.1
#define f 0.1
#define g 0.1
#define h 0.1


//float CalculateAngleAccel(float AccelX,float AccelZ);
float CalculateAngleAccel(float AccelX,float AccelY,float AccelZ);
float CalculateAngleRateGyro(int16_t GyroX);

void MatrixAddition(float* A, float* B, int m, int n, float* C);
void MatrixSubtraction(float* A, float* B, int m, int n, float* C);
void MatrixMultiply(float* A, float* B, int m, int p, int n, float* C);
int MatrixInversion(float* A, int n, float* AInverse);



#endif
