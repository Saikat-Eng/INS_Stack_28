/**
 * @file matrix.c
 * @brief Implementation of matrix operations for the INS system.
 *
 * This file contains functions for matrix initialization, copying, and multiplication.
 *
 * @author saikat
 * @date 2023-12-01
 * @version 1.0
 */

#include "matrix.h"
#include "math.h"

// size in byte
void mat_reset(float *d, int size)
{
  do
  {
    size--;
    d[size] = 0.0f;
  }while(size);
}

void mat_fill(float *B, int size)
{
    while (size--)
    {
        *B++ = 0.0f;
    }
}

void mat_copy28(float *B, const float *A)
{
  for (int i = 0; i < 784; i += 28)
  {
      B[i]  = A[i];
      B[i+1]  = A[i+1];
      B[i+2]  = A[i+2];
      B[i+3]  = A[i+3];
      B[i+4]  = A[i+4];
      B[i+5]  = A[i+5];
      B[i+6]  = A[i+6];
      B[i+7]  = A[i+7];
      B[i+8]  = A[i+8];
      B[i+9]  = A[i+9];
      B[i+10] = A[i+10];
      B[i+11] = A[i+11];
      B[i+12] = A[i+12];
      B[i+13] = A[i+13];
      B[i+14] = A[i+14];
      B[i+15] = A[i+15];
      B[i+16] = A[i+16];
      B[i+17] = A[i+17];
      B[i+18] = A[i+18];
      B[i+19] = A[i+19];
      B[i+20] = A[i+20];
      B[i+21] = A[i+21];
      B[i+22] = A[i+22];
      B[i+23] = A[i+23];
      B[i+24] = A[i+24];
      B[i+25] = A[i+25];
      B[i+26] = A[i+26];
      B[i+27] = A[i+27];
  }
}

void mat_cpy9(float *A, float *B)
{
    B[0] = A[0];
    B[1] = A[1];
    B[2] = A[2];
    B[3] = A[3];
    B[4] = A[4];
    B[5] = A[5];
    B[6] = A[6];
    B[7] = A[7];
    B[8] = A[8];
}

inline void mat_mul9_3(const float *A, const float *B, float *C)
{
    C[0] = A[0] * B[0] + A[1] * B[1] + A[2] * B[2];
    C[1] = A[3] * B[0] + A[4] * B[1] + A[5] * B[2];
    C[2] = A[6] * B[0] + A[7] * B[1] + A[8] * B[2];
}

inline void mat_mul3_9(const float *A, const float *B, float *C)
{
    C[0] = A[0] * B[0] + A[1] * B[3] + A[2] * B[6];
    C[1] = A[0] * B[1] + A[1] * B[4] + A[2] * B[7];
    C[2] = A[0] * B[2] + A[1] * B[5] + A[2] * B[8];
}


void mat_mul(const float *a, const float *b, float *c)
{
    c[0] = a[0] * b[0] + a[1] * b[3] + a[2] * b[6];
    c[1] = a[0] * b[1] + a[1] * b[4] + a[2] * b[7];
    c[2] = a[0] * b[2] + a[1] * b[5] + a[2] * b[8];

    c[3] = a[3] * b[0] + a[4] * b[3] + a[5] * b[6];
    c[4] = a[3] * b[1] + a[4] * b[4] + a[5] * b[7];
    c[5] = a[3] * b[2] + a[4] * b[5] + a[5] * b[8];

    c[6] = a[6] * b[0] + a[7] * b[3] + a[8] * b[6];
    c[7] = a[6] * b[1] + a[7] * b[4] + a[8] * b[7];
    c[8] = a[6] * b[2] + a[7] * b[5] + a[8] * b[8];
}

void mat3_transpose(float *R)
{
    float tmp;
    tmp = R[1]; R[1] = R[3]; R[3] = tmp;
    tmp = R[2]; R[2] = R[6]; R[6] = tmp;
    tmp = R[5]; R[5] = R[7]; R[7] = tmp;
}

char invert_3x3(const float *A, float *Ainv)
{
    const float a00 = A[0];
    const float a01 = A[1];
    const float a02 = A[2];

    const float a10 = A[3];
    const float a11 = A[4];
    const float a12 = A[5];

    const float a20 = A[6];
    const float a21 = A[7];
    const float a22 = A[8];

    const float c00 =  (a11 * a22 - a12 * a21);
    const float c01 = -(a10 * a22 - a12 * a20);
    const float c02 =  (a10 * a21 - a11 * a20);

    const float c10 = -(a01 * a22 - a02 * a21);
    const float c11 =  (a00 * a22 - a02 * a20);
    const float c12 = -(a00 * a21 - a01 * a20);

    const float c20 =  (a01 * a12 - a02 * a11);
    const float c21 = -(a00 * a12 - a02 * a10);
    const float c22 =  (a00 * a11 - a01 * a10);

    const float det = a00 * c00 + a01 * c01 + a02 * c02;

    if (!isfinite(det) || fabsf(det) < 1.0e-12f) {
        return 0;
    }

    const float inv_det = 1.0f / det;

    /*
     * Ainv = adj(A) / det
     * Note transpose of cofactor matrix.
     */
    Ainv[0] = c00 * inv_det;
    Ainv[1] = c10 * inv_det;
    Ainv[2] = c20 * inv_det;

    Ainv[3] = c01 * inv_det;
    Ainv[4] = c11 * inv_det;
    Ainv[5] = c21 * inv_det;

    Ainv[6] = c02 * inv_det;
    Ainv[7] = c12 * inv_det;
    Ainv[8] = c22 * inv_det;

//    for (int i = 0; i < 9; i++) {
//        if (!isfinite(Ainv[i])) {
//            return 0;
//        }
//    }

    return 1;
}

void skewm(const float *v, float *S)
{
    /*  S = [ 0 -z   y; */
    /*        z  0  -x; */
    /*       -y  x   0; ]; */
    S[0] = 0.0;
    S[1] = -v[2];
    S[2] = v[1];
    S[3] = v[2];
    S[4] = 0.0;
    S[5] = -v[0];
    S[6] = -v[1];
    S[7] = v[0];
    S[8] = 0.0;
}
void skewm_inv(const float *S, float *v)
{
    v[0] = S[7];
    v[1] = S[2];
    v[2] = S[3];
}
