#include <stdlib.h> 
#include <stdio.h>
#include "mpu9255.h"
#include "spi.h"
#include "system.h"
#include "kalman.h"


#define t 2000

float a,b,c;
float AngleGyro;

//initial matrix
float A[3][3] = {{1,0,-dt},{0,0,-1},{0,0,1}}; //A is a 3*3 matrix
float At[3][3] = {{1,0,0},{0,0,0},{-dt,-1,1}}; //At is the transposed matrix of A
float B[3][2] = {{dt,0},{1,0},{0,0}}; //B is a 3*2 matrix
float H[2][3] = {{0,1,1},{1,0,0}}; //H is a 2*3 matrix
float Ht[3][2] = {{1,0},{0,1},{0,1}}; //Ht is the transposed matrix of H
float X[3][1] = {0}; //X is a 3*1 vector, represent the estimation of the states of the system
float Xp[3][1] = {0}; //Xp is a 3*1 vector, represent an a priori estimation of the states of the system
//float Z[2][1] = {{AngleRateGyro},{AngleAccel}}; //Z is a 2*1 vector, represent the measured outputs of the system
float Z[2][1] = {0,0};
//float U[2][1] = {{AngleRateGyro},{AngleAccel}}; //U is a 2*1 vector, represent the inputs of the system
float U[2][1] = {0,0};
float Q[3][3] = {0}; //Q is the process noise covariance matrix, 3*3
//float Q1[3][1] = {{GyroStanDev*dt*p+q}, {GyroStanDev*p+q}, {BiasStanDev*m+n}}; 
float Q1[3][1] = {0,0,0};
//Q = Q1 * Q2, Q1 is a 3*1 matrix
//float Q2[1][3] = {GyroStanDev*dt*p+q, GyroStanDev*p+q, BiasStanDev*m+n }; 
float Q2[1][3] = {0,0,0};
//Q2 is the transposed matrix of Q1
//float R[2][2] = {{GyroVar*e+f,0},{0,AccelVar*g+h}}; 
float R[2][2] = {{0,0},{0,0}};
//R is the measurement noise covariance matrix, 2*2
float P[3][3] = {0}; //P is the estimation of the error covariance matrix, 3*3
float Pp[3][3] = {0}; //Pp is an a priori estimation of the error covariance matrix, 3*3
float K[3][2] = {0}; //K is the Kalman gain matrix, 3*2
float I[3][3] = {{1,0,0},{0,1,0},{0,0,1}}; //I is a 3*3 unit matrix

//derived matrix
float AX[3][1] = {0}; //AX = A * X
float BU[3][1] = {0}; //BU = B * U
float AP[3][3] = {0}; //AP = A * P
float APAt[3][3] = {0}; //APAt = AP * At
float PpHt[3][2] = {0}; //PpHt = Pp * Ht
float HPpHt[2][2] = {0}; //HPpHt = H * PpHt
float HPpHtR[2][2] = {0}; //HPpHtR = HPpHt + R
float InvHPpHtR[2][2] = {0}; //InvHPpHtR = Inv(HPpHtR)
float HXp[2][1] = {0}; //HXp = H * Xp
float ZHXp[2][1] = {0}; //ZHXp = Z - HXp
float KZHXp[3][1] = {0}; //KZHXp = K * ZHXp
float KH[3][3] = {0}; //KH = K * H
float IKH[3][3] = {0}; //IKH = I – KH


uint8_t MPU9255_DataBuffer[14];
S_INT16_XYZ MPU9255_ACC_LAST;
S_INT16_XYZ MPU9255_GYRO_LAST;
S_INT16_XYZ MPU9255_MAG_LAST;
S_INT32_XYZ MPU9255_ACC_OFFSET;
S_INT32_XYZ MPU9255_GYRO_OFFSET;
S_INT32_XYZ MPU9255_MAG_OFFSET;
int16_t MPU9255_TEMP_LAST;

uint8_t MPU9255_Init(void){
	if(MPU9255_Read_Reg(WHO_AM_I)==0x73)
	{
		MPU9255_Write_Reg(USER_CTRL,0X10); //Ê¹ÄÜMPU9255SPI
		MPU9255_Write_Reg(PWR_MGMT_1,0X80);  //µçÔ´¹ÜÀí£¬¸´Î»MPU9255
		MPU9255_Write_Reg(SMPLRT_DIV,0x07);//²ÉÑùÂÊ1000/(1+7)=125HZ
		MPU9255_Write_Reg(CONFIG,GYRO_BW_5);				//ÍÓÂÝÒÇÓëÎÂ¶ÈµÍÍ¨ÂË²¨Æ÷ 0x06 5hz
		MPU9255_Write_Reg(GYRO_CONFIG,GYRO_1000DPS);  //ÍÓÂÝÒÇ²âÁ¿·¶Î§ 0X18 +1000 dps
		MPU9255_Write_Reg(ACCEL_CONFIG,ACC_2G); //¼ÓËÙ¶È¼Æ²âÁ¿·¶Î§ 0X18 +-2g
		MPU9255_Write_Reg(ACCEL_CONFIG2,ACC_BW_5);	//¼ÓËÙ¶È¼ÆµÍÍ¨ÂË²¨Æ÷£¬5Hz
	/*	MPU9255_Write_Reg(XG_OFFSET_H,0x00);
		MPU9255_Write_Reg(XG_OFFSET_L,0x00);
		MPU9255_Write_Reg(YG_OFFSET_H,0x00);
		MPU9255_Write_Reg(YG_OFFSET_L,0x00);
		MPU9255_Write_Reg(ZG_OFFSET_H,0x00);
		MPU9255_Write_Reg(ZG_OFFSET_L,0x00);*/
		MPU9255_ACC_OFFSET.X=0;
		MPU9255_ACC_OFFSET.Y=0;
		MPU9255_ACC_OFFSET.Z=0;
		
		MPU9255_ACC_LAST.X=0;
		MPU9255_ACC_LAST.Y=0;
		MPU9255_ACC_LAST.Z=0;
		
		for(uint16_t i=0;i<t;i++)
		{
			MPU9255_ReadValue(0);	
			if(i==6)
			{
				MPU9255_ACC_OFFSET.X=0;
				MPU9255_ACC_OFFSET.Y=0;
				MPU9255_ACC_OFFSET.Z=0;
			}
		}
		MPU9255_ACC_OFFSET.X/=t;
		MPU9255_ACC_OFFSET.Y/=t;
		MPU9255_ACC_OFFSET.Z=MPU9255_ACC_OFFSET.Z/t-16384;	
		MPU9255_GYRO_OFFSET.X/=t;
		MPU9255_GYRO_OFFSET.Y/=t;
		MPU9255_GYRO_OFFSET.Z/=t;

		return 1;
	}
	return 0;
}

//SPIÐ´¼Ä´æÆ÷
//reg:Ö¸¶¨µÄ¼Ä´æÆ÷µØÖ·
//value:Ð´ÈëµÄÖµ
uint8_t MPU9255_Write_Reg(uint8_t reg,uint8_t value)
{
	uint8_t status;
	SPI_MPU9255_CS_L;											  //Ê¹ÄÜSPI´«Êä
	status = HAL_SPI_Transmit(&hspi1, &reg, 1, 0xFFFF);
	status = HAL_SPI_Transmit(&hspi1, &value, 1, 0xFFFF);
	SPI_MPU9255_CS_H;										  	//½ûÖ¹MPU9255
	Delay(0xFFF);
	return(status);													//·µ»Ø×´Ì¬Öµ
}

//SPI¶ÁÈ¡Ö¸¶¨¼Ä´æÆ÷
//reg:Ö¸¶¨¼Ä´æÆ÷µÄµØÖ·
uint8_t MPU9255_Read_Reg(uint8_t reg)
{
	uint8_t reg_val;
	SPI_MPU9255_CS_L;	
	reg = reg|0x80;
	HAL_SPI_Transmit(&hspi1, &reg, 1, 0xFFFF);	 	//·¢ËÍ¶ÁÃüÁî+¼Ä´æÆ÷ºÅ
 	HAL_SPI_Receive(&hspi1, &reg_val, 1, 0xFFFF);				//¶ÁÈ¡¼Ä´æÆ÷Öµ
	SPI_MPU9255_CS_H;																//½ûÖ¹SPI´«Êä
	Delay(0xFFF);
	return(reg_val);
}

//SPI¶ÁMPU9255Êý¾Ý
uint8_t MPU9255_ReadValue(uint8_t status)
{
	uint8_t data=ACCEL_XOUT_H|0x80;
	SPI_MPU9255_CS_L;  //Ê¹ÄÜSPI´«Êä
	HAL_SPI_Transmit(&hspi1, &data, 1, 0xFFFF);
	HAL_SPI_Receive(&hspi1, MPU9255_DataBuffer, 14, 0xFFFF); //¹²¶ÁÈ¡14×Ö½ÚÊý¾Ý
	//init status
	if(status == 0)
	{
		MPU9255_ACC_OFFSET.X += ((int16_t)(MPU9255_DataBuffer[0]<<8)) | (MPU9255_DataBuffer)[1];
		MPU9255_ACC_OFFSET.Y += ((int16_t)(MPU9255_DataBuffer[2]<<8)) | (MPU9255_DataBuffer)[3];
		MPU9255_ACC_OFFSET.Z += ((int16_t)(MPU9255_DataBuffer[4]<<8)) | (MPU9255_DataBuffer)[5];
		MPU9255_GYRO_OFFSET.X += ((int16_t)(MPU9255_DataBuffer[8]<<8)) | (MPU9255_DataBuffer)[9];
		MPU9255_GYRO_OFFSET.Y += ((int16_t)(MPU9255_DataBuffer[10]<<8)) | (MPU9255_DataBuffer)[11];
		MPU9255_GYRO_OFFSET.Z += ((int16_t)(MPU9255_DataBuffer[12]<<8)) | (MPU9255_DataBuffer)[13];
	}
	//measure status
	else if(status == 1)
	{
//		//¼ÓËÙ¶È¼Æ
		
		MPU9255_ACC_LAST.X = (((int16_t)(MPU9255_DataBuffer[0]<<8)) | (MPU9255_DataBuffer[1])) - (int16_t)MPU9255_ACC_OFFSET.X;
		MPU9255_ACC_LAST.Y = (((int16_t)(MPU9255_DataBuffer[2]<<8)) | (MPU9255_DataBuffer[3])) - (int16_t)MPU9255_ACC_OFFSET.Y;
		MPU9255_ACC_LAST.Z = (((int16_t)(MPU9255_DataBuffer[4]<<8)) | (MPU9255_DataBuffer[5])) - (int16_t)MPU9255_ACC_OFFSET.Z;
	
		//ÎÂ¶È
		MPU9255_TEMP_LAST =  ((int16_t)(MPU9255_DataBuffer[6]<<8)) | (MPU9255_DataBuffer)[7];
		//ÍÓÂÝÒÇ
		MPU9255_GYRO_LAST.X = (((int16_t)(MPU9255_DataBuffer[8]<<8)) | (MPU9255_DataBuffer[9])) - (int16_t)MPU9255_GYRO_OFFSET.X;
		MPU9255_GYRO_LAST.Y = (((int16_t)(MPU9255_DataBuffer[10]<<8)) | (MPU9255_DataBuffer[11])) - (int16_t)MPU9255_GYRO_OFFSET.Y;
		MPU9255_GYRO_LAST.Z = (((int16_t)(MPU9255_DataBuffer[12]<<8)) | (MPU9255_DataBuffer[13])) - (int16_t)MPU9255_GYRO_OFFSET.Z;
	  
/*			
		MPU9255_ACC_LAST.X /= (float)16384.0;	
		MPU9255_ACC_LAST.Y /= (float)16384.0;
		MPU9255_ACC_LAST.Z /= (float)16384.0;		//½«ÔËËã²¿·Ö·Åµ½matlbaÀïÃæ
	
	MPU9255_GYRO_LAST.X /=32.8;
		MPU9255_GYRO_LAST.Y /=32.8;
		MPU9255_GYRO_LAST.Z /=32.8;
	*/	
	
	//	MPU9255_ACC_LAST.Z+=16384;
		
		a=CalculateAngleAccel((float)MPU9255_ACC_LAST.X,(float)MPU9255_ACC_LAST.Y,(float)MPU9255_ACC_LAST.Z);
		//a=CalculateAngleAccel(MPU9255_ACC_LAST.X,MPU9255_ACC_LAST.Z);//xÖáÓë³õÊ¼×´Ì¬µÄ¼Ð½Ç
		//a=a/PI*180;
		b=CalculateAngleRateGyro(MPU9255_GYRO_LAST.Y);//ÈÆYÖá×ª¹ýµÄ¼Ð½Ç

	/*
		MPU9255_DataBuffer[0] = (int16_t)MPU9255_ACC_LAST.X >> 8;
		MPU9255_DataBuffer[1] = (int16_t)MPU9255_ACC_LAST.X;
		MPU9255_DataBuffer[2] = (int16_t)MPU9255_ACC_LAST.Y >> 8;
		MPU9255_DataBuffer[3] = (int16_t)MPU9255_ACC_LAST.Y;
		MPU9255_DataBuffer[4] = (int16_t)MPU9255_ACC_LAST.Z >> 8;
		MPU9255_DataBuffer[5] = (int16_t)MPU9255_ACC_LAST.Z;
	*/
		MPU9255_DataBuffer[0] = (int16_t)a >> 8;
		MPU9255_DataBuffer[1] = (int16_t)a;
		MPU9255_DataBuffer[2] = (int16_t)b >> 8;
		MPU9255_DataBuffer[3] = (int16_t)b;
		
		Z[0][0] = -MPU9255_GYRO_LAST.Y*PI/180;
		Z[1][0] = a;
		U[0][0] = -MPU9255_GYRO_LAST.Y*PI/180;
		U[1][0] = a;
		Q1[0][0] = GyroStanDev*dt*p1+q1;
		Q1[1][0] = GyroStanDev*p1+q1;
		Q1[2][0] = BiasStanDev*m1+n1;
		Q2[0][0] = GyroStanDev*dt*p1+q1;
		Q2[0][1] = GyroStanDev*p1+q1;
		Q2[0][2] = BiasStanDev*m1+n1;
		R[0][0] = GyroVar*e+f;
		R[0][1] = 0;
		R[1][0] = 0;
		R[1][1] = AccelVar*g+h;
		//time update or prediction
		//Xp = A * X + B * U
		MatrixMultiply((float*)A,(float*)X,3,3,1,(float*)AX); //AX = A * X
		MatrixMultiply((float*)B,(float*)U,3,2,1,(float*)BU); //BU = B * U
		MatrixAddition((float*)AX,(float*)BU,3,1,(float*)Xp); //Xp = AX + BU
		//Pp = A * P * At + Q
		MatrixMultiply((float*)A,(float*)P,3,3,3,(float*)AP); //AP = A * P
		MatrixMultiply((float*)AP,(float*)At,3,3,3,(float*)APAt); //APAt = AP * At
		MatrixMultiply((float*)Q1,(float*)Q2,3,1,3,(float*)Q); //Q = Q1 * Q2
		MatrixAddition((float*)APAt,(float*)Q,3,3,(float*)Pp); //Pp = APAt + Q
		//measurement update or correction
		//K = Pp * Ht * Inv(H*Pp*Ht + R)
		MatrixMultiply((float*)Pp,(float*)Ht,3,3,2,(float*)PpHt); //PpHt = Pp * Ht
		MatrixMultiply((float*)H,(float*)PpHt,2,3,2,(float*)HPpHt); //HPpHt = H * PpHt
		MatrixAddition((float*)HPpHt,(float*)R,2,2,(float*)HPpHtR); //HPpHtR = HPpHt + R
		MatrixInversion((float*)HPpHtR,2,(float*)InvHPpHtR); //InvHPpHtR = Inv(HPpHtR)
		MatrixMultiply((float*)PpHt,(float*)InvHPpHtR,3,2,2,(float*)K); //K = PpHt * InvHPpHtR
		//X = Xp + K * (Z - H*Xp)
		MatrixMultiply((float*)H,(float*)Xp,2,3,1,(float*)HXp); //HXp = H * Xp
		MatrixSubtraction((float*)Z,(float*)HXp,2,1,(float*)ZHXp); //ZHXp = Z - HXp
		MatrixMultiply((float*)K,(float*)ZHXp,3,2,1,(float*)KZHXp); //KZHXp = K * ZHXp
		MatrixAddition((float*)Xp,(float*)KZHXp,3,1,(float*)X); //X= Xp + KZHXp
		//P = (I - K*H) * Pp
		MatrixMultiply((float*)K,(float*)H,3,2,3,(float*)KH); //KH = K * H
		MatrixSubtraction((float*)I,(float*)KH,3,3,(float*)IKH); //IKH = I - KH
		MatrixMultiply((float*)IKH,(float*)Pp,3,3,3,(float*)P); //P = IKH * Pp
		//the kalman filter outputs
		c = X[0][0];
	//	AngleRate = X[1][0];
		
		
	//	MPU9255_DataBuffer[4]=0;
	//	MPU9255_DataBuffer[5]=0;
		
		
		MPU9255_DataBuffer[4] = (int16_t)c >> 8;
		MPU9255_DataBuffer[5] = (int16_t)c;
		
	//	printf("%f %f",a,b);
			
		MPU9255_DataBuffer[8] = (int16_t)MPU9255_GYRO_LAST.X >> 8;
		MPU9255_DataBuffer[9] = (int16_t)MPU9255_GYRO_LAST.X;
		MPU9255_DataBuffer[10] = (int16_t)MPU9255_GYRO_LAST.Y >> 8;
		MPU9255_DataBuffer[11] = (int16_t)MPU9255_GYRO_LAST.Y;
		MPU9255_DataBuffer[12] = (int16_t)MPU9255_GYRO_LAST.Z >> 8;
		MPU9255_DataBuffer[13] = (int16_t)MPU9255_GYRO_LAST.Z;
		
		
	}
	SPI_MPU9255_CS_H;  //½ûÖ¹MPU9255
	return 14;
}
