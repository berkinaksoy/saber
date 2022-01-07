#include <stdio.h>
#include "api.h"
#include "poly.h"
#include "poly_mul.h"
#include "pack_unpack.h"
#include "cbd.h"
#include "fips202.h"

void MatrixVectorMul(const uint16_t A[SABER_L][SABER_L][SABER_N], const uint16_t sw[SABER_L][7][64], uint16_t res[SABER_L][SABER_N], int16_t transpose)
{
	uint16_t acc[889] = {0};
	uint16_t bcc[889];
   	int32_t i, j, k;

	for (i = 0; i < SABER_L; i++)
	{
            	for (j = 0; j < 889; j++) {
                    bcc[j] = 0;
            			    	   }
		for (j = 0; j < SABER_L; j++){
		
		if (transpose == 1)
		{	
				toom_cook_4way_evaluate(A[j][i], sw[j], acc);
				for (k = 0; k < 889; k++)
			 {
        			bcc[k] += acc[k];
  			 }
		}
		else 
		{
				toom_cook_4way_evaluate(A[i][j], sw[j], acc);
				for (k = 0; k < 889; k++)
			 {
        			bcc[k] += acc[k];
  			 }
					
		}				
						}
	toom_cook_4way_interpol(bcc, res[i]);

	}
}

void InnerProd(const uint16_t b[SABER_L][SABER_N], const uint16_t sw[SABER_L][7][64], uint16_t res[SABER_N])
{
	uint16_t acc[889] = {0};
	uint16_t bcc[889] = {0};
   	int32_t j, k;

	for (j = 0; j < SABER_L; j++)
	{
		toom_cook_4way_evaluate(b[j], sw[j], acc);
		
		for (k = 0; k < 889; k++)
			 {
        			bcc[k] += acc[k];
  			 }
	}

	toom_cook_4way_interpol(bcc, res);					
				
}

void GenMatrix(uint16_t A[SABER_L][SABER_L][SABER_N], const uint8_t seed[SABER_SEEDBYTES])
{
	uint8_t buf[SABER_L * SABER_POLYVECBYTES];
	int i;

	shake128(buf, sizeof(buf), seed, SABER_SEEDBYTES);

	for (i = 0; i < SABER_L; i++)
	{
		BS2POLVECq(buf + i * SABER_POLYVECBYTES, A[i]);
	}
}

void GenSecret(uint16_t s[SABER_L][SABER_N], const uint8_t seed[SABER_NOISE_SEEDBYTES])
{
	uint8_t buf[SABER_L * SABER_POLYCOINBYTES];
	size_t i;

	shake128(buf, sizeof(buf), seed, SABER_NOISE_SEEDBYTES);

	for (i = 0; i < SABER_L; i++)
	{
		cbd(s[i], buf + i * SABER_POLYCOINBYTES);
	}
}
