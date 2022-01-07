#include "poly_mul.h"
#include <stdint.h>
#include <string.h>

#define SCHB_N 16

#define N_RES (SABER_N << 1)
#define N_SB (SABER_N >> 2)
#define N_SB_RES (2*N_SB-1)


#define KARATSUBA_N 64

void evaluation_single(const uint16_t *b, uint16_t sw[7][64]){ // for precomputing B

	uint16_t r0, r1, r2, r3, r4, r5, r6, r7;

	//uint16_t bw1[N_SB], bw2[N_SB], bw3[N_SB], bw4[N_SB], bw5[N_SB], bw6[N_SB], bw7[N_SB]; //can be reduced

	uint16_t *B0, *B1, *B2, *B3;

	int j;

	B0 = (uint16_t*)b;
	B1 = (uint16_t*)&b[N_SB];
	B2 = (uint16_t*)&b[2*N_SB];
	B3 = (uint16_t*)&b[3*N_SB];

	for (j = 0; j < N_SB; ++j) {
		r0 = B0[j];
		r1 = B1[j];
		r2 = B2[j];
		r3 = B3[j];
		r4 = r0 + r2;
		r5 = r1 + r3;
		r6 = r4 + r5; r7 = r4 - r5;
		sw[2][j] = r6;
		sw[3][j] = r7;
		r4 = ((r0 << 2)+r2) << 1;
		r5 = (r1 << 2) + r3;
		r6 = r4 + r5; r7 = r4 - r5;
		sw[4][j] = r6;
		sw[5][j] = r7;
		r4 = (r3 << 3) + (r2 << 2) + (r1 << 1) + r0;
		sw[1][j] = r4; 
		sw[6][j] = r0;
		sw[0][j] = r3;
	}
}
static void karatsuba_simple(const uint16_t *a_1, const uint16_t *b_1, uint16_t *result_final) {
    uint16_t d01[KARATSUBA_N / 2 - 1];
    uint16_t d0123[KARATSUBA_N / 2 - 1];
    uint16_t d23[KARATSUBA_N / 2 - 1];
    uint16_t result_d01[KARATSUBA_N - 1];

    int32_t i, j;

    memset(result_d01, 0, (KARATSUBA_N - 1)*sizeof(uint16_t));
    memset(d01, 0, (KARATSUBA_N / 2 - 1)*sizeof(uint16_t));
    memset(d0123, 0, (KARATSUBA_N / 2 - 1)*sizeof(uint16_t));
    memset(d23, 0, (KARATSUBA_N / 2 - 1)*sizeof(uint16_t));
    memset(result_final, 0, (2 * KARATSUBA_N - 1)*sizeof(uint16_t));

    uint16_t acc1, acc2, acc3, acc4, acc5, acc6, acc7, acc8, acc9, acc10;


    for (i = 0; i < KARATSUBA_N / 4; i++) {
        acc1 = a_1[i]; //a0
        acc2 = a_1[i + KARATSUBA_N / 4]; //a1
        acc3 = a_1[i + 2 * KARATSUBA_N / 4]; //a2
        acc4 = a_1[i + 3 * KARATSUBA_N / 4]; //a3
        for (j = 0; j < KARATSUBA_N / 4; j++) {

            acc5 = b_1[j]; //b0
            acc6 = b_1[j + KARATSUBA_N / 4]; //b1

            result_final[i + j + 0 * KARATSUBA_N / 4] = result_final[i + j + 0 * KARATSUBA_N / 4] + acc1 * acc5;
            result_final[i + j + 2 * KARATSUBA_N / 4] = result_final[i + j + 2 * KARATSUBA_N / 4] + acc2 * acc6;

            acc7 = acc5 + acc6; //b01
            acc8 = acc1 + acc2; //a01
            d01[i + j] = d01[i + j] + acc7 * acc8;
            //--------------------------------------------------------

            acc7 = b_1[j + 2 * KARATSUBA_N / 4]; //b2
            acc8 = b_1[j + 3 * KARATSUBA_N / 4]; //b3
            result_final[i + j + 4 * KARATSUBA_N / 4] = result_final[i + j + 4 * KARATSUBA_N / 4] + acc7 * acc3;

            result_final[i + j + 6 * KARATSUBA_N / 4] = result_final[i + j + 6 * KARATSUBA_N / 4] + acc8 * acc4;

            acc9 = acc3 + acc4;
            acc10 = acc7 + acc8;
            d23[i + j] = d23[i + j] + acc9 * acc10;
            //--------------------------------------------------------

            acc5 = acc5 + acc7; //b02
            acc7 = acc1 + acc3; //a02
            result_d01[i + j + 0 * KARATSUBA_N / 4] = result_d01[i + j + 0 * KARATSUBA_N / 4] + acc5 * acc7;

            acc6 = acc6 + acc8; //b13
            acc8 = acc2 + acc4;
            result_d01[i + j + 2 * KARATSUBA_N / 4] = result_d01[i + j + 2 * KARATSUBA_N / 4] + acc6 * acc8;

            acc5 = acc5 + acc6;
            acc7 = acc7 + acc8;
            d0123[i + j] = d0123[i + j] + acc5 * acc7;
        }
    }

    // 2nd last stage

    for (i = 0; i < KARATSUBA_N / 2 - 1; i++) {
        d0123[i] = d0123[i] - result_d01[i + 0 * KARATSUBA_N / 4] - result_d01[i + 2 * KARATSUBA_N / 4];
        d01[i] = d01[i] - result_final[i + 0 * KARATSUBA_N / 4] - result_final[i + 2 * KARATSUBA_N / 4];
        d23[i] = d23[i] - result_final[i + 4 * KARATSUBA_N / 4] - result_final[i + 6 * KARATSUBA_N / 4];
    }

    for (i = 0; i < KARATSUBA_N / 2 - 1; i++) {
        result_d01[i + 1 * KARATSUBA_N / 4] = result_d01[i + 1 * KARATSUBA_N / 4] + d0123[i];
        result_final[i + 1 * KARATSUBA_N / 4] = result_final[i + 1 * KARATSUBA_N / 4] + d01[i];
        result_final[i + 5 * KARATSUBA_N / 4] = result_final[i + 5 * KARATSUBA_N / 4] + d23[i];
    }

    // Last stage
    for (i = 0; i < KARATSUBA_N - 1; i++) {
        result_d01[i] = result_d01[i] - result_final[i] - result_final[i + KARATSUBA_N];
    }

    for (i = 0; i < KARATSUBA_N - 1; i++) {
        result_final[i + 1 * KARATSUBA_N / 2] = result_final[i + 1 * KARATSUBA_N / 2] + result_d01[i];
    }

}







//-------------------------------------//lazy-----------------------//lazy---------------------



void toom_cook_4way_evaluate(const uint16_t *a1, const uint16_t sw[7][64], uint16_t *res) {
    //uint16_t inv3 = 43691, inv9 = 36409, inv15 = 61167;

    uint16_t aw1[N_SB], aw2[N_SB], aw3[N_SB], aw4[N_SB], aw5[N_SB], aw6[N_SB], aw7[N_SB];
    //uint16_t bw1[N_SB], bw2[N_SB], bw3[N_SB], bw4[N_SB], bw5[N_SB], bw6[N_SB], bw7[N_SB];
    uint16_t w1[N_SB_RES] = {0}, w2[N_SB_RES] = {0}, w3[N_SB_RES] = {0}, w4[N_SB_RES] = {0},
                            w5[N_SB_RES] = {0}, w6[N_SB_RES] = {0}, w7[N_SB_RES] = {0};
    uint16_t r0, r1, r2, r3, r4, r5, r6, r7;
    uint16_t *A0, *A1, *A2, *A3;
    // *B0, *B1, *B2, *B3;
    A0 = (uint16_t *)a1;
    A1 = (uint16_t *)&a1[N_SB];
    A2 = (uint16_t *)&a1[2 * N_SB];
    A3 = (uint16_t *)&a1[3 * N_SB];
    //B0 = (uint16_t *)b1;
    //B1 = (uint16_t *)&b1[N_SB];
    //B2 = (uint16_t *)&b1[2 * N_SB];
    //B3 = (uint16_t *)&b1[3 * N_SB];

    int i, j;
    
    uint16_t C[889];
        for (i = 0; i < 889; i++) {
        C[i] = 0;
    }
    //C = result;



    // EVALUATION
    for (j = 0; j < N_SB; ++j) {
        r0 = A0[j];
        r1 = A1[j];
        r2 = A2[j];
        r3 = A3[j];
        r4 = r0 + r2;
        r5 = r1 + r3;
        r6 = r4 + r5;
        r7 = r4 - r5;
        aw3[j] = r6;
        aw4[j] = r7;
        r4 = ((r0 << 2) + r2) << 1;
        r5 = (r1 << 2) + r3;
        r6 = r4 + r5;
        r7 = r4 - r5;
        aw5[j] = r6;
        aw6[j] = r7;
        r4 = (r3 << 3) + (r2 << 2) + (r1 << 1) + r0;
        aw2[j] = r4;
        aw7[j] = r0;
        aw1[j] = r3;
    }


    // MULTIPLICATION

    karatsuba_simple(aw1, sw[0], w1);
    karatsuba_simple(aw2, sw[1], w2);
    karatsuba_simple(aw3, sw[2], w3);
    karatsuba_simple(aw4, sw[3], w4);
    karatsuba_simple(aw5, sw[4], w5);
    karatsuba_simple(aw6, sw[5], w6);
    karatsuba_simple(aw7, sw[6], w7);

    // INTERPOLATION
    for (i = 0; i < N_SB_RES; ++i) {
        r0 = w1[i];
        r1 = w2[i];
        r2 = w3[i];
        r3 = w4[i];
        r4 = w5[i];
        r5 = w6[i];
        r6 = w7[i];

        C[i]     = r0;
        C[i + 127]  = r1;
        C[i + 254] = r2;
        C[i + 381] = r3;
        C[i + 508] = r4;
        C[i + 635] = r5;
        C[i + 762] = r6;


    }
    
        for (i = 0; i < 889; i++) {
        res[i] = (C[i]);
    }
}

void toom_cook_4way_interpol(const uint16_t *a_i, uint16_t *res_i){

   	int32_t m, n;
   	
   	uint16_t r0, r1, r2, r3, r4, r5, r6;
	uint16_t inv3 = 43691, inv9 = 36409, inv15 = 61167;
	uint16_t c[512] = {0};

					
					for (m = 0; m < 127; ++m) {

					r0 = a_i[m];
					r1 = a_i[m+127];
					r2 = a_i[m+254];
					r3 = a_i[m+381];
					r4 = a_i[m+508];
					r5 = a_i[m+635];
					r6 = a_i[m+762];

					r1 = r1 + r4;
					r5 = r5 - r4;
					r3 = ((r3 - r2) >> 1);
					r4 = r4 - r0;
					r4 = r4 - (r6 << 6);
					r4 = (r4 << 1) + r5;
					r2 = r2 + r3;
					r1 = r1 - (r2 << 6) - r2;
					r2 = r2 - r6;
					r2 = r2 - r0;
					r1 = r1 + 45 * r2;
					r4 = (((r4 - (r2 << 3)) * inv3) >> 3);
					r5 = r5 + r1;
					r1 = (((r1 + (r3 << 4)) * inv9) >> 1);
					r3 = -(r3 + r1);
					r5 = (((30 * r1 - r5) * inv15) >> 2);
					r2 = r2 - r4;
					r1 = r1 - r5;

					c[m]     += r6;
					c[m + 64]  += r5;
					c[m + 128] += r4;
					c[m + 192] += r3;
					c[m + 256] += r2;
					c[m + 320] += r1;
					c[m + 384] += r0;
				}
			 		for (n = 256; n < 512; n++) {
				res_i[n - 256] = (c[n - 256] - c[n]) & 8191;
							}
							
}

