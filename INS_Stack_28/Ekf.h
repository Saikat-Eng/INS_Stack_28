/*
 * File: IMUBasicEKF.h
 *
 * MATLAB Coder version            : 4.2
 * C/C++ source code generated on  : 14-Jun-2023 18:12:47
 */

#ifndef IMUBASICEKF_H
#define IMUBASICEKF_H


/* Function Declarations */
char IMUBasicEKF_correctEqn(float *x, float *P, const float *H, float *innov, float *b_innovCov, const float *R);
char GPSBasicEKF_correctEqn(float *x, float *P, const float *H, float *innov, float *innovCov, const float *R);
void ekf_fix_covariance(float *P, int n);
#endif

/*
 * File trailer for IMUBasicEKF.h
 *
 * [EOF]
 */
