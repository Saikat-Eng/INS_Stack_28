/**
 * @file INSEstimator.c 
 * @brief Implementation of the INS estimator and EKF fusion logic.
 *
 * This file contains the implementation of the INS estimator, including initialization,
 * state updates, and configuration management for the Extended Kalman Filter (EKF).
 *
 * @author saikat
 * @date 2023-12-01
 * @version 1.0
 */
#include <string.h>
#include <INSEstimator.h>
#include <INSFuserBase.h>
#include <Ekf.h>
#include <util.h>
#include <drv_utils.h>

insfilterAsync fusionfilt __attribute__((section(".DTCM_RAM")));
insfilterAsync_Cfg fusionfilt_cfg __attribute__((section(".DTCM_RAM")));

void (*ekf_reset_callback)(void) = NULL;

/**
 * @brief Get pointer to singleton EKF instance.
 *
 * Returns a pointer to the global insfilterAsync instance used throughout the
 * application. Use this to access EKF state and internal structures.
 *
 * @return insfilterAsync* pointer to the EKF instance.
 */
insfilterAsync* get_ekf_instance(void)
{
    return &fusionfilt;
}

/**
 * @brief Get pointer to EKF configuration structure.
 *
 * Returns a pointer to the global insfilterAsync_Cfg structure which holds
 * measurement and process noise configuration used by the EKF.
 *
 * @return insfilterAsync_Cfg* pointer to the EKF configuration.
 */
insfilterAsync_Cfg* ekf_cfg_instance(void)
{
    return &fusionfilt_cfg;
}

/**
 * @brief Register a callback to be invoked when EKF is reset.
 *
 * Stores the provided function pointer and will call it from reset_ekf().
 * Only the first registered callback is stored; subsequent calls will be ignored.
 *
 * @param fptr Function pointer to call on EKF reset (void(void)).
 */
void resister_ekf_reset_callback(void (*fptr)(void))
{
    if(ekf_reset_callback == NULL)
        ekf_reset_callback = fptr;
}

/**
 * @brief Initialize EKF state and configuration.
 *
 * Initializes the EKF state vector, covariance, and process/measurement noise
 * parameters. Sets the reference location, initial quaternion, estimated biases,
 * geomagnetic field and related tuning parameters.
 *
 * @param ref_loc Pointer to reference location [lat, lon, alt] (double array).
 * @param mag_field Pointer to geomagnetic field vector (float[3], NED, uT).
 * @param q_init Initial orientation quaternion (float[4], order [w,x,y,z]).
 */


void insfilt_init(const double *ref_loc, const float *mag_field, float *q_init)
{
    memset(&fusionfilt,0,sizeof(fusionfilt));
    memset(&fusionfilt_cfg,0,sizeof(fusionfilt_cfg));

    fusionfilt.ReferenceLocation[0] = ref_loc[0];
    fusionfilt.ReferenceLocation[1] = ref_loc[1];
    fusionfilt.ReferenceLocation[2] = ref_loc[2];

    get_norm_quat(q_init);

        /*
        %    States                            Units    Index
        %    Orientation (quaternion parts)             1:4  
        %    Angular Velocity (XYZ)            rad/s    5:7  
        %    Position (NED)                    m        8:10 
        %    Velocity (NED)                    m/s      11:13
        %    Acceleration (NED)                m/s^2    14:16
        %    Accelerometer Bias (XYZ)          m/s^2    17:19
        %    Gyroscope Bias (XYZ)              rad/s    20:22
        %    Geomagnetic Field Vector (NED)    uT       23:25
        %    Magnetometer Bias (XYZ)           uT       26:28
        */

    fusionfilt.pState[0] = q_init[0];
    fusionfilt.pState[1] = q_init[1];
    fusionfilt.pState[2] = q_init[2];
    fusionfilt.pState[3] = q_init[3];
    
    fusionfilt.pState[22] = mag_field[0];
    fusionfilt.pState[23] = mag_field[1];
    fusionfilt.pState[24] = mag_field[2];   

 /*
>Increase in Sensor Bias Noise:

Reduced Accuracy: The state estimates become less accurate because the filter has to account for higher uncertainty in the measurements.
Higher Innovation: The difference between the predicted and actual measurements (innovation) will generally increase,
leading to larger corrections in the state estimates.
Increased Covariance: The covariance of the state estimate increases, indicating higher uncertainty in the estimated states.

>Decrease in Sensor Bias Noise:

Improved Accuracy: The state estimates become more accurate as the filter can rely more on the measurements.
Lower Innovation: The innovation decreases, leading to smaller corrections in the state estimates.
Decreased Covariance: The covariance of the state estimate decreases, indicating lower uncertainty in the estimated states.

Process Noise: This represents the uncertainty in the system model.
If the process noise covariance is too low, the filter may not account for all the dynamics of the system,
leading to drift over time. Increasing the process noise covariance can help the filter adapt better to changes and reduce drift.

Measurement Noise: This represents the uncertainty in the sensor measurements.
If the measurement noise covariance is not accurately set, the filter might give too much weight to noisy measurements,
causing drift. Ensuring that the measurement noise covariance accurately reflects the sensor's noise characteristics is crucial.

if velocity process noise is less it's trusting model less so more more in vartical measurment and less
accel based velocity measurment
*/
//**************************************************************************************************************************************
    
    fusionfilt_cfg.Rmag[0] = 65.0f;
    fusionfilt_cfg.Rmag[4] = 65.0f;
    fusionfilt_cfg.Rmag[8] = 65.0f;                                                                      
    
    fusionfilt_cfg.Racc[0] = 0.061f;
    fusionfilt_cfg.Racc[4] = 0.061f;
    fusionfilt_cfg.Racc[8] = 0.95f;

    fusionfilt_cfg.Rgyro[0] = 5.0e-4f;
    fusionfilt_cfg.Rgyro[4] = 5.0e-4f;
    fusionfilt_cfg.Rgyro[8] = 5.0e-4f;

    fusionfilt_cfg.Rgps[0] = 2.5f;
    fusionfilt_cfg.Rgps[7] = 2.5f;
    fusionfilt_cfg.Rgps[14]= 2.5f;   // posD (baro height): trust baro more to rein in altitude drift

    fusionfilt_cfg.Rgps[21] = 0.05f;
    fusionfilt_cfg.Rgps[28] = 0.05f;
    fusionfilt_cfg.Rgps[35] = 0.25f; // velD (baro climb rate): trust more to pull back vertical velocity drift

  //Strong ZUPT, only when definitely stationary/landed:
    fusionfilt_cfg.Rzupt[0] = 0.025f;  //strong when decreasing , weak when increasing
    fusionfilt_cfg.Rzupt[4] = 0.025f;  //
    fusionfilt_cfg.Rzupt[8] = 0.45f;   //

    //Process Noise
    fusionfilt_cfg.PositionNoise[0] = 5.0e-3f;
    fusionfilt_cfg.PositionNoise[1] = 5.0e-3f;
    fusionfilt_cfg.PositionNoise[2] = 0.25f;	//

    fusionfilt_cfg.VelocityNoise[0] = 1.0e-3f;
    fusionfilt_cfg.VelocityNoise[1] = 1.0e-3f;
    fusionfilt_cfg.VelocityNoise[2] = 0.5f;

    fusionfilt_cfg.QuaternionNoise[0] = 2.5e-7f;
    fusionfilt_cfg.QuaternionNoise[1] = 2.5e-7f;
    fusionfilt_cfg.QuaternionNoise[2] = 2.5e-7f;
    fusionfilt_cfg.QuaternionNoise[3] = 2.5e-7f;
    
    fusionfilt_cfg.AccelerationNoise[0] = 25.0f;
    fusionfilt_cfg.AccelerationNoise[1] = 25.0f;
    fusionfilt_cfg.AccelerationNoise[2] = 25.0f;

    fusionfilt_cfg.AngularVelocityNoise[0] = 0.5f;
    fusionfilt_cfg.AngularVelocityNoise[1] = 0.5f;
    fusionfilt_cfg.AngularVelocityNoise[2] = 0.5f;

    //Increasing it makes earth magnetic field bias estimation faster and noisier.
    fusionfilt_cfg.GeomagneticVectorNoise[0] = 1.0e-8f;
    fusionfilt_cfg.GeomagneticVectorNoise[1] = 1.0e-8f;
    fusionfilt_cfg.GeomagneticVectorNoise[2] = 1.0e-8f;

    //Increasing it makes compass offset estimation faster and noisier..
    fusionfilt_cfg.MagnetometerBiasNoise[0] = 1.0e-3f;
    fusionfilt_cfg.MagnetometerBiasNoise[1] = 1.0e-3f;
    fusionfilt_cfg.MagnetometerBiasNoise[2] = 1.0e-3f;

    //Increasing it makes accelerometer bias estimation faster and noisier.
    fusionfilt_cfg.AccelerometerBiasNoise[0] = 1.0e-6f;
    fusionfilt_cfg.AccelerometerBiasNoise[1] = 1.0e-6f;
    fusionfilt_cfg.AccelerometerBiasNoise[2] = 1.0e-2f;

    //Increasing it makes rate gyro bias estimation faster and noisier.
    fusionfilt_cfg.GyroscopeBiasNoise[0] = 5.0e-12f;
    fusionfilt_cfg.GyroscopeBiasNoise[1] = 5.0e-12f;
    fusionfilt_cfg.GyroscopeBiasNoise[2] = 5.0e-12f;

    //**************************************************************************************************************************************
    int i=0;

    for(i=0;i<28;i++)
        fusionfilt.pStateCovariance[i*28+i] = 1.0e-5f;
 
    fusionfilt.addProcNoise[0]   = fusionfilt_cfg.QuaternionNoise[0];
    fusionfilt.addProcNoise[29]  = fusionfilt_cfg.QuaternionNoise[1];
    fusionfilt.addProcNoise[58]  = fusionfilt_cfg.QuaternionNoise[2];
    fusionfilt.addProcNoise[87]  = fusionfilt_cfg.QuaternionNoise[3];
    fusionfilt.addProcNoise[116] = fusionfilt_cfg.AngularVelocityNoise[0];
    fusionfilt.addProcNoise[145] = fusionfilt_cfg.AngularVelocityNoise[1];
    fusionfilt.addProcNoise[174] = fusionfilt_cfg.AngularVelocityNoise[2];
    fusionfilt.addProcNoise[203] = fusionfilt_cfg.PositionNoise[0];
    fusionfilt.addProcNoise[232] = fusionfilt_cfg.PositionNoise[1];
    fusionfilt.addProcNoise[261] = fusionfilt_cfg.PositionNoise[2];
    fusionfilt.addProcNoise[290] = fusionfilt_cfg.VelocityNoise[0];
    fusionfilt.addProcNoise[319] = fusionfilt_cfg.VelocityNoise[1];
    fusionfilt.addProcNoise[348] = fusionfilt_cfg.VelocityNoise[2];
    fusionfilt.addProcNoise[377] = fusionfilt_cfg.AccelerationNoise[0];
    fusionfilt.addProcNoise[406] = fusionfilt_cfg.AccelerationNoise[1];
    fusionfilt.addProcNoise[435] = fusionfilt_cfg.AccelerationNoise[2];
    fusionfilt.addProcNoise[464] = fusionfilt_cfg.AccelerometerBiasNoise[0];
    fusionfilt.addProcNoise[493] = fusionfilt_cfg.AccelerometerBiasNoise[1];
    fusionfilt.addProcNoise[522] = fusionfilt_cfg.AccelerometerBiasNoise[2];
    fusionfilt.addProcNoise[551] = fusionfilt_cfg.GyroscopeBiasNoise[0];
    fusionfilt.addProcNoise[580] = fusionfilt_cfg.GyroscopeBiasNoise[1];
    fusionfilt.addProcNoise[609] = fusionfilt_cfg.GyroscopeBiasNoise[2];
    fusionfilt.addProcNoise[638] = fusionfilt_cfg.GeomagneticVectorNoise[0];
    fusionfilt.addProcNoise[667] = fusionfilt_cfg.GeomagneticVectorNoise[1];
    fusionfilt.addProcNoise[696] = fusionfilt_cfg.GeomagneticVectorNoise[2];
    fusionfilt.addProcNoise[725] = fusionfilt_cfg.MagnetometerBiasNoise[0];
    fusionfilt.addProcNoise[754] = fusionfilt_cfg.MagnetometerBiasNoise[1];
    fusionfilt.addProcNoise[783] = fusionfilt_cfg.MagnetometerBiasNoise[2];
}

/**
 * @brief Update EKF process noise vector from configuration.
 *
 * Copies the current noise values from fusionfilt_cfg into the EKF's addProcNoise
 * array. Call this after modifying fusionfilt_cfg to apply new tuning values.
 */
void update_process_noise(void)
{
    fusionfilt.addProcNoise[0]   = fusionfilt_cfg.QuaternionNoise[0];
    fusionfilt.addProcNoise[29]  = fusionfilt_cfg.QuaternionNoise[1];
    fusionfilt.addProcNoise[58]  = fusionfilt_cfg.QuaternionNoise[2];
    fusionfilt.addProcNoise[87]  = fusionfilt_cfg.QuaternionNoise[3];
    fusionfilt.addProcNoise[116] = fusionfilt_cfg.AngularVelocityNoise[0];
    fusionfilt.addProcNoise[145] = fusionfilt_cfg.AngularVelocityNoise[1];
    fusionfilt.addProcNoise[174] = fusionfilt_cfg.AngularVelocityNoise[2];
    fusionfilt.addProcNoise[203] = fusionfilt_cfg.PositionNoise[0];
    fusionfilt.addProcNoise[232] = fusionfilt_cfg.PositionNoise[1];
    fusionfilt.addProcNoise[261] = fusionfilt_cfg.PositionNoise[2];
    fusionfilt.addProcNoise[290] = fusionfilt_cfg.VelocityNoise[0];
    fusionfilt.addProcNoise[319] = fusionfilt_cfg.VelocityNoise[1];
    fusionfilt.addProcNoise[348] = fusionfilt_cfg.VelocityNoise[2];
    fusionfilt.addProcNoise[377] = fusionfilt_cfg.AccelerationNoise[0];
    fusionfilt.addProcNoise[406] = fusionfilt_cfg.AccelerationNoise[1];
    fusionfilt.addProcNoise[435] = fusionfilt_cfg.AccelerationNoise[2];
    fusionfilt.addProcNoise[464] = fusionfilt_cfg.AccelerometerBiasNoise[0];
    fusionfilt.addProcNoise[493] = fusionfilt_cfg.AccelerometerBiasNoise[1];
    fusionfilt.addProcNoise[522] = fusionfilt_cfg.AccelerometerBiasNoise[2];
    fusionfilt.addProcNoise[551] = fusionfilt_cfg.GyroscopeBiasNoise[0];
    fusionfilt.addProcNoise[580] = fusionfilt_cfg.GyroscopeBiasNoise[1];
    fusionfilt.addProcNoise[609] = fusionfilt_cfg.GyroscopeBiasNoise[2];
    fusionfilt.addProcNoise[638] = fusionfilt_cfg.GeomagneticVectorNoise[0];
    fusionfilt.addProcNoise[667] = fusionfilt_cfg.GeomagneticVectorNoise[1];
    fusionfilt.addProcNoise[696] = fusionfilt_cfg.GeomagneticVectorNoise[2];
    fusionfilt.addProcNoise[725] = fusionfilt_cfg.MagnetometerBiasNoise[0];
    fusionfilt.addProcNoise[754] = fusionfilt_cfg.MagnetometerBiasNoise[1];
    fusionfilt.addProcNoise[783] = fusionfilt_cfg.MagnetometerBiasNoise[2];
}

/**
 * @brief Return pointer to EKF state vector.
 *
 * Exposes the EKF internal state vector (pState). Use only for read access or
 * controlled diagnostics. The returned pointer references fusionfilt.pState.
 *
 * @return float* pointer to 28-element state array.
 */
float *get_state(void)
{
    return fusionfilt.pState;
}

/**
 * @brief Trigger EKF reset and invoke registered callback.
 *
 * If a reset callback was previously registered with resister_ekf_reset_callback,
 * this function will call it. Use reset_ekf() to reinitialize EKF-related state
 * from higher level tasks.
 */
void reset_ekf(void)
{
    if(ekf_reset_callback != NULL)
    {
        print_Debug("Reseting EKF\n");
        ekf_reset_callback();
    }
}

/**
 * @brief Update EKF reference position.
 *
 * Sets the EKF internal reference location (latitude, longitude, altitude).
 * This is the origin used for converting between geodetic and local frames.
 *
 * @param ref_loc Pointer to double[3] with {lat, lon, alt}.
 */
void set_ref_pos_ekf(const double *ref_loc)
{
    fusionfilt.ReferenceLocation[0] = ref_loc[0];
    fusionfilt.ReferenceLocation[1] = ref_loc[1];
    fusionfilt.ReferenceLocation[2] = ref_loc[2];
    
    print_Debug("Reference Position Update Done\n");
}

/**
 * @brief Update EKF reference altitude.
 *
 * Convenience function to update only the altitude component of the EKF
 * reference location.
 *
 * @param alt Reference altitude in same units used by the EKF (meters).
 */
void set_ref_alt_ekf(const float alt)
{
   fusionfilt.ReferenceLocation[2] = alt;
    
   print_Debug("Reference Altitude Update Done\n");
}

void set_zupt_pos_hold(char set)
{
  if(set)
  {
    fusionfilt_cfg.Rzupt[0] = 0.0005f;  //strong when decreasing , weak when increasing
    fusionfilt_cfg.Rzupt[4] = 0.0005f; 
    fusionfilt_cfg.Rzupt[8] = 0.45f; 
  }
  else
  {
    fusionfilt_cfg.Rzupt[0] = 0.025f;  //strong when decreasing , weak when increasing
    fusionfilt_cfg.Rzupt[4] = 0.025f;   
    fusionfilt_cfg.Rzupt[8] = 0.45f; 
  }
}

/*
 * File trailer for run_PoseEstimationFromAsynchronousSensorsExample.c
 *
 * [EOF]
 */
