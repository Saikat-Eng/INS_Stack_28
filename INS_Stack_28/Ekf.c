/**
 * @file Ekf.c  
 * @brief   Implementation of the Extended Kalman Filter (EKF) for INS/IMU data processing.
 *
 * This file contains the implementation of the Extended Kalman Filter (EKF) used for state estimation
 * in inertial navigation systems, including prediction and correction steps.
 *
 * @author saikat
 * @date 2023-12-01
 * @version 1.0
 */
#include <math.h>
#include "util.h"
#include "Ekf.h"
#include "matrix.h"
#include <stdio.h>


/* Function Definitions 

function [innov, innovCov] = equationInnovation(P, h, H, z, R)
  innovCov = H*P*(H.') + R; //S
  innov = (z(:) - h(:));

        
function [x, P, innov, innovCov] = equationCorrect(x, P, h, H, z, R)
  [innov, innovCov] = positioning.internal.EKF.equationInnovation(P, h, H, z, R);
  
  W = P*(H.') / innovCov; //Kalman gain
  x = x + W*(innov);
  P = P - W*H*P;
 
// with Joseph version
    W = P*(H.') / innovCov;
    x = x + W*innov;
    I = eye(size(P));
    A = I - W*H;
    P = A*P*A.' + W*R*W.';
end

*/

static const signed char iv[784] = {
      1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,
      0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
      1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 1};

//covariance symmetrization
void ekf_fix_covariance(float *P, int n)
{
  for (int r = 0; r < n; r++) {
          for (int c = r + 1; c < n; c++) {
              float s = 0.5f * (P[r + n*c] + P[c + n*r]);
              P[r + n*c] = s;
              P[c + n*r] = s;
          }
          float *prr = &P[r + n*r];
          if (!isfinite(*prr) || *prr < 1e-9f) *prr = 1e-9f;
      }
}


/*
 * Arguments    : float x[28]
 *                float P[784]
 *                const float h[3]
 *                const float H[84]
 *                const float z[3]
 *                const float R[9]
 * Return Type  : void
 */

char IMUBasicEKF_correctEqn(float *x, float *P, const float *H, float *innov, float *b_innovCov, const float *R)
{
  static float A[784] __attribute__((section(".DTCM_RAM")));
  static float b_P[784]__attribute__((section(".DTCM_RAM")));
  static float W[84] __attribute__((section(".DTCM_RAM")));
  static float y[84] __attribute__((section(".DTCM_RAM")));

  float a21;
  float f;
  float maxval;
  int b_innovCov_tmp;
  int c_innovCov_tmp;
  int d_innovCov_tmp;
  int e_innovCov_tmp;
  int innovCov_tmp;
  int k;
  int r1;
  int r2;
  int r3;
  int rtemp;

   for (innovCov_tmp = 0; innovCov_tmp < 3; innovCov_tmp++) {
    for (b_innovCov_tmp = 0; b_innovCov_tmp < 28; b_innovCov_tmp++) {
      maxval = 0.0F;
      for (rtemp = 0; rtemp < 28; rtemp++) {
        maxval += H[rtemp + 28 * innovCov_tmp] * P[rtemp + 28 * b_innovCov_tmp];
      }
      y[innovCov_tmp + 3 * b_innovCov_tmp] = maxval;
    }
  }
  r1 = 0;
  r2 = 1;
  r3 = 2;
  maxval = fabsf(b_innovCov[0]);
  a21 = fabsf(b_innovCov[3]);
  if (a21 > maxval) {
    maxval = a21;
    r1 = 1;
    r2 = 0;
  }
  if (fabsf(b_innovCov[6]) > maxval) {
    r1 = 2;
    r2 = 1;
    r3 = 0;
  }
  b_innovCov[3 * r2] /= b_innovCov[3 * r1];
  b_innovCov[3 * r3] /= b_innovCov[3 * r1];
  innovCov_tmp = 3 * r2 + 1;
  c_innovCov_tmp = 3 * r1 + 1;
  b_innovCov[innovCov_tmp] -= b_innovCov[3 * r2] * b_innovCov[c_innovCov_tmp];
  b_innovCov_tmp = 3 * r3 + 1;
  b_innovCov[b_innovCov_tmp] -= b_innovCov[3 * r3] * b_innovCov[c_innovCov_tmp];
  rtemp = 3 * r2 + 2;
  d_innovCov_tmp = 3 * r1 + 2;
  b_innovCov[rtemp] -= b_innovCov[3 * r2] * b_innovCov[d_innovCov_tmp];
  rtemp = 3 * r3 + 2;
  b_innovCov[rtemp] -= b_innovCov[3 * r3] * b_innovCov[d_innovCov_tmp];
  if (fabsf(b_innovCov[b_innovCov_tmp]) > fabsf(b_innovCov[innovCov_tmp])) {
    rtemp = r2;
    r2 = r3;
    r3 = rtemp;
  }
  innovCov_tmp = 3 * r3 + 1;
  b_innovCov_tmp = 3 * r2 + 1;
  b_innovCov[innovCov_tmp] /= b_innovCov[b_innovCov_tmp];
  rtemp = 3 * r3 + 2;
  e_innovCov_tmp = 3 * r2 + 2;
  b_innovCov[rtemp] -= b_innovCov[innovCov_tmp] * b_innovCov[e_innovCov_tmp];
  maxval = innov[0];
  a21 = innov[1];
  f = innov[2];
  for (k = 0; k < 28; k++) {
    int W_tmp;
    int b_W_tmp;
    int c_W_tmp;
    int d_W_tmp;
    int e_W_tmp;
    W_tmp = r1 + 3 * k;
    W[W_tmp] = y[3 * k] / b_innovCov[3 * r1];
    b_W_tmp = r2 + 3 * k;
    c_W_tmp = 3 * k + 1;
    W[b_W_tmp] = y[c_W_tmp] - W[W_tmp] * b_innovCov[c_innovCov_tmp];
    d_W_tmp = r3 + 3 * k;
    e_W_tmp = 3 * k + 2;
    W[d_W_tmp] = y[e_W_tmp] - W[W_tmp] * b_innovCov[d_innovCov_tmp];
    W[b_W_tmp] /= b_innovCov[b_innovCov_tmp];
    W[d_W_tmp] -= W[b_W_tmp] * b_innovCov[e_innovCov_tmp];
    W[d_W_tmp] /= b_innovCov[rtemp];
    W[b_W_tmp] -= W[d_W_tmp] * b_innovCov[innovCov_tmp];
    W[W_tmp] -= W[d_W_tmp] * b_innovCov[3 * r3];
    W[W_tmp] -= W[b_W_tmp] * b_innovCov[3 * r2];
    x[k] += (maxval * W[3 * k] + a21 * W[c_W_tmp]) + f * W[e_W_tmp];
  }
  for (innovCov_tmp = 0; innovCov_tmp < 28; innovCov_tmp++) {
    maxval = H[innovCov_tmp];
    a21 = H[innovCov_tmp + 28];
    f = H[innovCov_tmp + 56];
    for (b_innovCov_tmp = 0; b_innovCov_tmp < 28; b_innovCov_tmp++) {
      rtemp = innovCov_tmp + 28 * b_innovCov_tmp;
      A[rtemp] =
          (float)iv[rtemp] -
          ((maxval * W[3 * b_innovCov_tmp] + a21 * W[3 * b_innovCov_tmp + 1]) +
           f * W[3 * b_innovCov_tmp + 2]);
    }
  }
  for (innovCov_tmp = 0; innovCov_tmp < 28; innovCov_tmp++) {
    for (b_innovCov_tmp = 0; b_innovCov_tmp < 28; b_innovCov_tmp++) {
      maxval = 0.0F;
      for (rtemp = 0; rtemp < 28; rtemp++) {
        maxval += P[innovCov_tmp + 28 * rtemp] * A[rtemp + 28 * b_innovCov_tmp];
      }
      b_P[innovCov_tmp + 28 * b_innovCov_tmp] = maxval;
    }
  }
  for (innovCov_tmp = 0; innovCov_tmp < 3; innovCov_tmp++) {
    maxval = R[innovCov_tmp];
    a21 = R[innovCov_tmp + 3];
    f = R[innovCov_tmp + 6];
    for (b_innovCov_tmp = 0; b_innovCov_tmp < 28; b_innovCov_tmp++) {
      y[innovCov_tmp + 3 * b_innovCov_tmp] =
          (maxval * W[3 * b_innovCov_tmp] + a21 * W[3 * b_innovCov_tmp + 1]) +
          f * W[3 * b_innovCov_tmp + 2];
    }
  }
  for (innovCov_tmp = 0; innovCov_tmp < 28; innovCov_tmp++) {
    for (b_innovCov_tmp = 0; b_innovCov_tmp < 28; b_innovCov_tmp++) {
      maxval = 0.0F;
      for (rtemp = 0; rtemp < 28; rtemp++) {
        maxval +=
            A[rtemp + 28 * innovCov_tmp] * b_P[rtemp + 28 * b_innovCov_tmp];
      }
      P[innovCov_tmp + 28 * b_innovCov_tmp] = maxval;
    }
  }
  for (innovCov_tmp = 0; innovCov_tmp < 28; innovCov_tmp++) {
    maxval = W[3 * innovCov_tmp];
    a21 = W[3 * innovCov_tmp + 1];
    f = W[3 * innovCov_tmp + 2];
    for (b_innovCov_tmp = 0; b_innovCov_tmp < 28; b_innovCov_tmp++) {
      A[innovCov_tmp + 28 * b_innovCov_tmp] =
          (maxval * y[3 * b_innovCov_tmp] + a21 * y[3 * b_innovCov_tmp + 1]) +
          f * y[3 * b_innovCov_tmp + 2];
    }
  }
//  for (innovCov_tmp = 0; innovCov_tmp < 784; innovCov_tmp++) {
//    P[innovCov_tmp] += A[innovCov_tmp];
//  }  
     for (k = 0; k < 28; k++)
    {
        r3 = 0+k;
        P[r3] += A[r3];
        r3 = 28+k;
        P[r3] += A[r3];
        r3 = 56+k;
        P[r3] += A[r3];
        r3 = 84+k;
        P[r3] += A[r3];
        r3 = 112+k;
        P[r3] += A[r3];
        r3 = 140+k;
        P[r3] += A[r3];
        r3 = 168+k;
        P[r3] += A[r3];
        r3 = 196+k;
        P[r3] += A[r3];
        r3 = 224+k;
        P[r3] += A[r3];
        r3 = 252+k;
        P[r3] += A[r3];
        r3 = 280+k;
        P[r3] += A[r3];
        r3 = 308+k;
        P[r3] += A[r3];
        r3 = 336+k;
        P[r3] += A[r3];
        r3 = 364+k;
        P[r3] += A[r3];
        r3 = 392+k;
        P[r3] += A[r3];
        r3 = 420+k;
        P[r3] += A[r3];
        r3 = 448+k;
        P[r3] += A[r3];
        r3 = 476+k;
        P[r3] += A[r3];
        r3 = 504+k;
        P[r3] += A[r3];
        r3 = 532+k;
        P[r3] += A[r3];
        r3 = 560+k;
        P[r3] += A[r3];
        r3 = 588+k;
        P[r3] += A[r3];
        r3 = 616+k;
        P[r3] += A[r3];
        r3 = 644+k;
        P[r3] += A[r3];
        r3 = 672+k;
        P[r3] += A[r3];
        r3 = 700+k;
        P[r3] += A[r3];
        r3 = 728+k;
        P[r3] += A[r3];
        r3 = 756+k;
        P[r3] += A[r3];
    }    
    
  return 1;
}


/*
 * Arguments    : float x[28]
 *                float P[784]
 *                const float h[6]
 *                const float H[168]
 *                const float z[6]
 *                const float R[36]
 * Return Type  : void
 */
 
char GPSBasicEKF_correctEqn(float *x, float *P, const float *H, float *innov, float *innovCov, const float *R)
{
  static float b_A[784] __attribute__((section(".DTCM_RAM")));
  static float b_P[784] __attribute__((section(".DTCM_RAM")));
  static float W[168] __attribute__((section(".DTCM_RAM")));
  static float X[168] __attribute__((section(".DTCM_RAM")));
  static float A[36] __attribute__((section(".DTCM_RAM")));
 
  float f;
  float f1;
  float f2;
  float f3;
  float s;
  float smax;
  int b_i;
  int i;
  int ijA;
  int j;
  int jA;
  int jBcol;
  int jp1j;
  int k;
  signed char ipiv[6];
  
  A[0] = innovCov[0];
  A[1] = innovCov[6];
  A[2] = innovCov[12];
  A[3] = innovCov[18];
  A[4] = innovCov[24];
  A[5] = innovCov[30];
  A[6] = innovCov[1];
  A[7] = innovCov[7];
  A[8] = innovCov[13];
  A[9] = innovCov[19];
  A[10] = innovCov[25];
  A[11] = innovCov[31];
  A[12] = innovCov[2];
  A[13] = innovCov[8];
  A[14] = innovCov[14];
  A[15] = innovCov[20];
  A[16] = innovCov[26];
  A[17] = innovCov[32];
  A[18] = innovCov[3];
  A[19] = innovCov[9];
  A[20] = innovCov[15];
  A[21] = innovCov[21];
  A[22] = innovCov[27];
  A[23] = innovCov[33];
  A[24] = innovCov[4];
  A[25] = innovCov[10];
  A[26] = innovCov[16];
  A[27] = innovCov[22];
  A[28] = innovCov[28];
  A[29] = innovCov[34];
  A[30] = innovCov[5];
  A[31] = innovCov[11];
  A[32] = innovCov[17];
  A[33] = innovCov[23];
  A[34] = innovCov[29];
  A[35] = innovCov[35];
  
//  for (i = 0; i < 28; i++) {
//    for (jp1j = 0; jp1j < 6; jp1j++) {
//      smax = 0.0F;
//      for (jA = 0; jA < 28; jA++) {
//        smax += H[jA + 28 * jp1j] * P[jA + 28 * i];
//      }
//      X[i + 28 * jp1j] = smax;
//    }
//  }
//as H is a identity matrix
  X[0] = 1.0f*P[7];
  X[28] = 1.0f*P[8];
  X[56] = 1.0f*P[9];
  X[84] = 1.0f*P[10];
  X[112] = 1.0f*P[11];
  X[140] = 1.0f*P[12];
  X[1] = 1.0f*P[35];
  X[29] = 1.0f*P[36];
  X[57] = 1.0f*P[37];
  X[85] = 1.0f*P[38];
  X[113] = 1.0f*P[39];
  X[141] = 1.0f*P[40];
  X[2] = 1.0f*P[63];
  X[30] = 1.0f*P[64];
  X[58] = 1.0f*P[65];
  X[86] = 1.0f*P[66];
  X[114] = 1.0f*P[67];
  X[142] = 1.0f*P[68];
  X[3] = 1.0f*P[91];
  X[31] = 1.0f*P[92];
  X[59] = 1.0f*P[93];
  X[87] = 1.0f*P[94];
  X[115] = 1.0f*P[95];
  X[143] = 1.0f*P[96];
  X[4] = 1.0f*P[119];
  X[32] = 1.0f*P[120];
  X[60] = 1.0f*P[121];
  X[88] = 1.0f*P[122];
  X[116] = 1.0f*P[123];
  X[144] = 1.0f*P[124];
  X[5] = 1.0f*P[147];
  X[33] = 1.0f*P[148];
  X[61] = 1.0f*P[149];
  X[89] = 1.0f*P[150];
  X[117] = 1.0f*P[151];
  X[145] = 1.0f*P[152];
  X[6] = 1.0f*P[175];
  X[34] = 1.0f*P[176];
  X[62] = 1.0f*P[177];
  X[90] = 1.0f*P[178];
  X[118] = 1.0f*P[179];
  X[146] = 1.0f*P[180];
  X[7] = 1.0f*P[203];
  X[35] = 1.0f*P[204];
  X[63] = 1.0f*P[205];
  X[91] = 1.0f*P[206];
  X[119] = 1.0f*P[207];
  X[147] = 1.0f*P[208];
  X[8] = 1.0f*P[231];
  X[36] = 1.0f*P[232];
  X[64] = 1.0f*P[233];
  X[92] = 1.0f*P[234];
  X[120] = 1.0f*P[235];
  X[148] = 1.0f*P[236];
  X[9] = 1.0f*P[259];
  X[37] = 1.0f*P[260];
  X[65] = 1.0f*P[261];
  X[93] = 1.0f*P[262];
  X[121] = 1.0f*P[263];
  X[149] = 1.0f*P[264];
  X[10] = 1.0f*P[287];
  X[38] = 1.0f*P[288];
  X[66] = 1.0f*P[289];
  X[94] = 1.0f*P[290];
  X[122] = 1.0f*P[291];
  X[150] = 1.0f*P[292];
  X[11] = 1.0f*P[315];
  X[39] = 1.0f*P[316];
  X[67] = 1.0f*P[317];
  X[95] = 1.0f*P[318];
  X[123] = 1.0f*P[319];
  X[151] = 1.0f*P[320];
  X[12] = 1.0f*P[343];
  X[40] = 1.0f*P[344];
  X[68] = 1.0f*P[345];
  X[96] = 1.0f*P[346];
  X[124] = 1.0f*P[347];
  X[152] = 1.0f*P[348];
  X[13] = 1.0f*P[371];
  X[41] = 1.0f*P[372];
  X[69] = 1.0f*P[373];
  X[97] = 1.0f*P[374];
  X[125] = 1.0f*P[375];
  X[153] = 1.0f*P[376];
  X[14] = 1.0f*P[399];
  X[42] = 1.0f*P[400];
  X[70] = 1.0f*P[401];
  X[98] = 1.0f*P[402];
  X[126] = 1.0f*P[403];
  X[154] = 1.0f*P[404];
  X[15] = 1.0f*P[427];
  X[43] = 1.0f*P[428];
  X[71] = 1.0f*P[429];
  X[99] = 1.0f*P[430];
  X[127] = 1.0f*P[431];
  X[155] = 1.0f*P[432];
  X[16] = 1.0f*P[455];
  X[44] = 1.0f*P[456];
  X[72] = 1.0f*P[457];
  X[100] = 1.0f*P[458];
  X[128] = 1.0f*P[459];
  X[156] = 1.0f*P[460];
  X[17] = 1.0f*P[483];
  X[45] = 1.0f*P[484];
  X[73] = 1.0f*P[485];
  X[101] = 1.0f*P[486];
  X[129] = 1.0f*P[487];
  X[157] = 1.0f*P[488];
  X[18] = 1.0f*P[511];
  X[46] = 1.0f*P[512];
  X[74] = 1.0f*P[513];
  X[102] = 1.0f*P[514];
  X[130] = 1.0f*P[515];
  X[158] = 1.0f*P[516];
  X[19] = 1.0f*P[539];
  X[47] = 1.0f*P[540];
  X[75] = 1.0f*P[541];
  X[103] = 1.0f*P[542];
  X[131] = 1.0f*P[543];
  X[159] = 1.0f*P[544];
  X[20] = 1.0f*P[567];
  X[48] = 1.0f*P[568];
  X[76] = 1.0f*P[569];
  X[104] = 1.0f*P[570];
  X[132] = 1.0f*P[571];
  X[160] = 1.0f*P[572];
  X[21] = 1.0f*P[595];
  X[49] = 1.0f*P[596];
  X[77] = 1.0f*P[597];
  X[105] = 1.0f*P[598];
  X[133] = 1.0f*P[599];
  X[161] = 1.0f*P[600];
  X[22] = 1.0f*P[623];
  X[50] = 1.0f*P[624];
  X[78] = 1.0f*P[625];
  X[106] = 1.0f*P[626];
  X[134] = 1.0f*P[627];
  X[162] = 1.0f*P[628];
  X[23] = 1.0f*P[651];
  X[51] = 1.0f*P[652];
  X[79] = 1.0f*P[653];
  X[107] = 1.0f*P[654];
  X[135] = 1.0f*P[655];
  X[163] = 1.0f*P[656];
  X[24] = 1.0f*P[679];
  X[52] = 1.0f*P[680];
  X[80] = 1.0f*P[681];
  X[108] = 1.0f*P[682];
  X[136] = 1.0f*P[683];
  X[164] = 1.0f*P[684];
  X[25] = 1.0f*P[707];
  X[53] = 1.0f*P[708];
  X[81] = 1.0f*P[709];
  X[109] = 1.0f*P[710];
  X[137] = 1.0f*P[711];
  X[165] = 1.0f*P[712];
  X[26] = 1.0f*P[735];
  X[54] = 1.0f*P[736];
  X[82] = 1.0f*P[737];
  X[110] = 1.0f*P[738];
  X[138] = 1.0f*P[739];
  X[166] = 1.0f*P[740];
  X[27] = 1.0f*P[763];
  X[55] = 1.0f*P[764];
  X[83] = 1.0f*P[765];
  X[111] = 1.0f*P[766];
  X[139] = 1.0f*P[767];
  X[167] = 1.0f*P[768];
  
  ipiv[0] = 1;
  ipiv[1] = 2;
  ipiv[2] = 3;
  ipiv[3] = 4;
  ipiv[4] = 5;
  ipiv[5] = 6;
  for (j = 0; j < 5; j++) {
    int b_tmp;
    int mmj_tmp;
    mmj_tmp = 4 - j;
    b_tmp = j * 7;
    jp1j = b_tmp + 2;
    jA = 6 - j;
    jBcol = 0;
    smax = fabsf(A[b_tmp]);
    for (k = 2; k <= jA; k++) {
      s = fabsf(A[(b_tmp + k) - 1]);
      if (s > smax) {
        jBcol = k - 1;
        smax = s;
      }
    }
    if (A[b_tmp + jBcol] != 0.0F) {
      if (jBcol != 0) {
        jA = j + jBcol;
        ipiv[j] = (signed char)(jA + 1);
        smax = A[j];
        A[j] = A[jA];
        A[jA] = smax;
        smax = A[j + 6];
        A[j + 6] = A[jA + 6];
        A[jA + 6] = smax;
        smax = A[j + 12];
        A[j + 12] = A[jA + 12];
        A[jA + 12] = smax;
        smax = A[j + 18];
        A[j + 18] = A[jA + 18];
        A[jA + 18] = smax;
        smax = A[j + 24];
        A[j + 24] = A[jA + 24];
        A[jA + 24] = smax;
        smax = A[j + 30];
        A[j + 30] = A[jA + 30];
        A[jA + 30] = smax;
      }
      i = (b_tmp - j) + 6;
      for (b_i = jp1j; b_i <= i; b_i++) {
        A[b_i - 1] /= A[b_tmp];
      }
    }
    jA = b_tmp;
    for (jBcol = 0; jBcol <= mmj_tmp; jBcol++) {
      smax = A[(b_tmp + jBcol * 6) + 6];
      if (smax != 0.0F) {
        i = jA + 8;
        jp1j = (jA - j) + 12;
        for (ijA = i; ijA <= jp1j; ijA++) {
          A[ijA - 1] += A[((b_tmp + ijA) - jA) - 7] * -smax;
        }
      }
      jA += 6;
    }
  }
  for (j = 0; j < 6; j++) {
    jBcol = 28 * j - 1;
    jA = 6 * j;
    for (k = 0; k < j; k++) {
      jp1j = 28 * k;
      smax = A[k + jA];
      if (smax != 0.0F) {
        for (b_i = 0; b_i < 28; b_i++) {
          ijA = (b_i + jBcol) + 1;
          X[ijA] -= smax * X[b_i + jp1j];
        }
      }
    }
    smax = 1.0F / A[j + jA];
    for (b_i = 0; b_i < 28; b_i++) {
      ijA = (b_i + jBcol) + 1;
      X[ijA] *= smax;
    }
  }
  for (j = 5; j >= 0; j--) {
    jBcol = 28 * j - 1;
    jA = 6 * j - 1;
    i = j + 2;
    for (k = i; k < 7; k++) {
      jp1j = 28 * (k - 1);
      smax = A[k + jA];
      if (smax != 0.0F) {
        for (b_i = 0; b_i < 28; b_i++) {
          ijA = (b_i + jBcol) + 1;
          X[ijA] -= smax * X[b_i + jp1j];
        }
      }
    }
  }
  for (j = 4; j >= 0; j--) {
    signed char i1;
    i1 = ipiv[j];
    if (i1 != j + 1) {
      for (b_i = 0; b_i < 28; b_i++) {
        jA = b_i + 28 * j;
        smax = X[jA];
        ijA = b_i + 28 * (i1 - 1);
        X[jA] = X[ijA];
        X[ijA] = smax;
      }
    }
  }
  smax = innov[0];
  s = innov[1];
  f = innov[2];
  f1 = innov[3];
  f2 = innov[4];
  f3 = innov[5];
  for (i = 0; i < 28; i++) {
    float f4;
    float f5;
    f4 = X[i];
    W[6 * i] = f4;
    f5 = smax * f4;
    f4 = X[i + 28];
    W[6 * i + 1] = f4;
    f5 += s * f4;
    f4 = X[i + 56];
    W[6 * i + 2] = f4;
    f5 += f * f4;
    f4 = X[i + 84];
    W[6 * i + 3] = f4;
    f5 += f1 * f4;
    f4 = X[i + 112];
    W[6 * i + 4] = f4;
    f5 += f2 * f4;
    f4 = X[i + 140];
    W[6 * i + 5] = f4;
    f5 += f3 * f4;
    x[i] += f5;
  }
  for (i = 0; i < 28; i++) {
    smax = H[i];
    s = H[i + 28];
    f = H[i + 56];
    f1 = H[i + 84];
    f2 = H[i + 112];
    f3 = H[i + 140];
    for (jp1j = 0; jp1j < 28; jp1j++) {
      jA = i + 28 * jp1j;
      b_A[jA] = (float)iv[jA] - (((((smax * W[6 * jp1j] + s * W[6 * jp1j + 1]) +
                                    f * W[6 * jp1j + 2]) +
                                   f1 * W[6 * jp1j + 3]) +
                                  f2 * W[6 * jp1j + 4]) +
                                 f3 * W[6 * jp1j + 5]);
    }
  }
  for (i = 0; i < 28; i++) {
    for (jp1j = 0; jp1j < 28; jp1j++) {
      smax = 0.0F;
      for (jA = 0; jA < 28; jA++) {
        smax += P[i + 28 * jA] * b_A[jA + 28 * jp1j];
      }
      b_P[i + 28 * jp1j] = smax;
    }
  }
  for (i = 0; i < 6; i++) {
    smax = R[i];
    s = R[i + 6];
    f = R[i + 12];
    f1 = R[i + 18];
    f2 = R[i + 24];
    f3 = R[i + 30];
    for (jp1j = 0; jp1j < 28; jp1j++) {
      X[i + 6 * jp1j] =
          ((((smax * W[6 * jp1j] + s * W[6 * jp1j + 1]) + f * W[6 * jp1j + 2]) +
            f1 * W[6 * jp1j + 3]) +
           f2 * W[6 * jp1j + 4]) +
          f3 * W[6 * jp1j + 5];
    }
  }
  for (i = 0; i < 28; i++) {
    for (jp1j = 0; jp1j < 28; jp1j++) {
      smax = 0.0F;
      for (jA = 0; jA < 28; jA++) {
        smax += b_A[jA + 28 * i] * b_P[jA + 28 * jp1j];
      }
      P[i + 28 * jp1j] = smax;
    }
  }
  for (i = 0; i < 28; i++) {
    smax = W[6 * i];
    s = W[6 * i + 1];
    f = W[6 * i + 2];
    f1 = W[6 * i + 3];
    f2 = W[6 * i + 4];
    f3 = W[6 * i + 5];
    for (jp1j = 0; jp1j < 28; jp1j++) {
      b_A[i + 28 * jp1j] =
          ((((smax * X[6 * jp1j] + s * X[6 * jp1j + 1]) + f * X[6 * jp1j + 2]) +
            f1 * X[6 * jp1j + 3]) +
           f2 * X[6 * jp1j + 4]) +
          f3 * X[6 * jp1j + 5];
    }
  }
//  for (i = 0; i < 784; i++) {
//    P[i] += b_A[i];
//  }
      for (k = 0; k < 28; k++)
    {
        j = 0+k;
        P[j] += b_A[j];
        j = 28+k;
        P[j] += b_A[j];
        j = 56+k;
        P[j] += b_A[j];
        j = 84+k;
        P[j] += b_A[j];
        j = 112+k;
        P[j] += b_A[j];
        j = 140+k;
        P[j] += b_A[j];
        j = 168+k;
        P[j] += b_A[j];
        j = 196+k;
        P[j] += b_A[j];
        j = 224+k;
        P[j] += b_A[j];
        j = 252+k;
        P[j] += b_A[j];
        j = 280+k;
        P[j] += b_A[j];
        j = 308+k;
        P[j] += b_A[j];
        j = 336+k;
        P[j] += b_A[j];
        j = 364+k;
        P[j] += b_A[j];
        j = 392+k;
        P[j] += b_A[j];
        j = 420+k;
        P[j] += b_A[j];
        j = 448+k;
        P[j] += b_A[j];
        j = 476+k;
        P[j] += b_A[j];
        j = 504+k;
        P[j] += b_A[j];
        j = 532+k;
        P[j] += b_A[j];
        j = 560+k;
        P[j] += b_A[j];
        j = 588+k;
        P[j] += b_A[j];
        j = 616+k;
        P[j] += b_A[j];
        j = 644+k;
        P[j] += b_A[j];
        j = 672+k;
        P[j] += b_A[j];
        j = 700+k;
        P[j] += b_A[j];
        j = 728+k;
        P[j] += b_A[j];
        j = 756+k;
        P[j] += b_A[j];
    }
   return 1;
}



/*
 * File trailer for IMUBasicEKF.c
 *
 * [EOF]
 */
