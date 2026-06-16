/*
 * File: PoseEstimatior.h
 *
 * MATLAB Coder version            : 4.2
 * C/C++ source code generated on  : 14-Jun-2023 18:12:47
 */

#ifndef RUN_POSEESTIMATIONFROMASYNCHRONOUSSENSORSEXAMPLE_H
#define RUN_POSEESTIMATIONFROMASYNCHRONOUSSENSORSEXAMPLE_H

typedef struct {
    double ReferenceLocation[3];
    float pState[28];
    float pStateCovariance[784];
    float addProcNoise[784];
    float accl_innov[3];
    float gyro_innov[3];
    float mag_innov[3];
    float gps_innov[6];
} insfilterAsync;


typedef struct {

    float Racc[9];
    float Rgyro[9];
    float Rmag[9];
    float Rgps[36];
    float Rzupt[9];

    float QuaternionNoise[4];
    float AngularVelocityNoise[3];
    float PositionNoise[3];
    float VelocityNoise[3];
    float AccelerationNoise[3];
    float GyroscopeBiasNoise[3];
    float AccelerometerBiasNoise[3];
    float GeomagneticVectorNoise[3];
    float MagnetometerBiasNoise[3];
}insfilterAsync_Cfg;

extern insfilterAsync* get_ekf_instance(void);
extern insfilterAsync_Cfg* ekf_cfg_instance(void);
extern void insfilt_init(const double *ref_loc, const float *mag_field, float *q_init);
void update_process_noise(void);
extern void resister_ekf_reset_callback(void (*fptr)(void));
extern float *get_state(void);
extern void reset_ekf(void);
extern void set_ref_pos_ekf(const double *ref_loc);
extern void set_ref_alt_ekf(const float alt);
extern void set_zupt_pos_hold(char set);
#endif

/*
 * File trailer for run_PoseEstimationFromAsynchronousSensorsExample.h
 *
 * [EOF]
 */
