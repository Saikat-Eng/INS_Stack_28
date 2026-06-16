#include "zupt.h"
#include <math.h>
#include "util.h"
#include "Ekf.h"

#define RM(row, col, ncols) ((row) * (ncols) + (col))

static float nis_3d(const float *innov, const float *Sinv);
                                       
float nis_3d(const float *innov, const float *Sinv)
{
    float tmp[3];

    for (int r = 0; r < 3; r++) {
        tmp[r] = 0.0f;

        for (int c = 0; c < 3; c++) {
            tmp[r] += Sinv[RM(r, c, 3)] * innov[c];
        }
    }

    return innov[0] * tmp[0] +
           innov[1] * tmp[1] +
           innov[2] * tmp[2]; 
}

void INS_Filt_fusezupt(insfilterAsync *obj)
{
    float H[84]  = {0.0f};      // row-major: 3 rows, 28 columns
    float innov[3];
    float S[9] = {0.0f};       // row-major: 3 x 3

    H[RM(0, 10, 28)] = 1.0f;
    H[RM(1, 11, 28)] = 1.0f;
    H[RM(2, 12, 28)] = 1.0f;

    innov[0] = -obj->pState[10];
    innov[1] = -obj->pState[11];
    innov[2] = -obj->pState[12];

    /*
     * S = H P H' + R
     *
     * Since H selects only velocity states,
     * S(i,j) = P[vel_i, vel_j] + Rzupt(i,j)
     */
    const float *R = ekf_cfg_instance()->Rzupt;
    
    const int v[3] = {10, 11, 12};
    for (int r = 0; r < 3; r++) {
        for (int c = 0; c < 3; c++) {
            S[RM(r, c, 3)] = obj->pStateCovariance[RM(v[r], v[c], 28)] + R[RM(r, c, 3)];
        }
    }

    float nis =
      innov[0] * innov[0] / S[0] +
      innov[1] * innov[1] / S[4] +
      innov[2] * innov[2] / S[8];

    if(nis < 10.81f){
      if (!IMUBasicEKF_correctEqn(obj->pState,obj->pStateCovariance,H,innov,S,R))
      {
          reset_ekf();
          return;
      }
    }
    get_norm_quat(obj->pState);
}

char zupt_detector(const float *accel, const float *gyro, const float *gps_vel, float dt, char position_hold_mode, char landed)
{
    static float static_time = 0.0f;

    float a_norm = sqrtf(
        accel[0] * accel[0] +
        accel[1] * accel[1] +
        accel[2] * accel[2]);

    float w_norm = sqrtf(
        gyro[0] * gyro[0] +
        gyro[1] * gyro[1] +
        gyro[2] * gyro[2]);

    float v_norm = sqrtf(
        gps_vel[0] * gps_vel[0] +
        gps_vel[1] * gps_vel[1] +
        gps_vel[2] * gps_vel[2]);

    char imu_static =
        (fabsf(a_norm - Actual_gravity) < 0.25f) && (fabsf(w_norm )< 0.05f);

    char gps_static =
        fabsf(v_norm) < 0.25f;

    char allowed =
        landed || position_hold_mode;

    if (allowed && imu_static && gps_static) {
        static_time += dt;
    } else {
        static_time = 0.0f;
    }

    return static_time>0.5f;
}
