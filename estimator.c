/**
 * @file estimator.c
 * @brief Implementation of the estimator module for the NAV_SYS_M7 project.
 *
 * This file contains the implementation of the estimator module, which processes sensor data
 * (IMU, GPS, magnetometer) to estimate the state of the system (position, velocity, attitude).
 *
 * @author Saikat
 * @date 2023-11-25
 * @version 1.0
 */
#include <stm32f7xx_hal.h>
#include "cmsis_os2.h"
#include <drv_utils.h>
#include <task.h>
#include <imu.h>
#include <nav.h>
#include <nav_sys_config.h>
#include <flight_control.h>
#include <matrix.h>
#include <system.h>
#include <estimator.h>
#include <string.h>
#include <math.h>
#include <INSEstimator.h>
#include <INSFuserBase.h>
#include <GeoBase.h>
#include <GeoMag.h>
#include <util.h>
#include <matrix.h>
#include <imu_calibration.h>
#include <quternion.h>
#include <zupt.h>


#define mag_fusing_us 20000U  //~50 hz
#define gps_fusing_us 50000U //~20  hz

static volatile e_data estimated_data __attribute__((section(".DTCM_RAM")));
volatile float dt_ahrs;
volatile char rst_flag = 0;
static osMutexId_t ekf_mutex;
static float mMag_sf = 0.0f;
static alpha_3d lpf_accel, lpf_gyro, lpf_mag;

void estimator_init(void)
{
	int i = 0;
	double pos_ref[3] = {0, 0, 0};
	double pos_ecef[3] = {0, 0, 0};
	float q[4] = {0, 0, 0, 0};
	float q_sum[4] = {0, 0, 0, 0};
	float accel[3] = {0, 0, 0};
  float gyro[3] = {0, 0, 0};
	float mag[3] = {0, 0, 0};
	float mag_field[3] = {0, 0, 0};

	ekf_mutex = osMutexNew(NULL);

	/*resetting EKF*/
	system_state.ekf_health = 0;

	pos_ref[0] = struct_home_data.lla_origin[N];
	pos_ref[1] = struct_home_data.lla_origin[E];
	pos_ref[2] = struct_home_data.lla_origin[D];

	memset(struct_mag_cal.Mag_vectors, 0, sizeof(struct_mag_cal.Mag_vectors));

	lla2ecef(pos_ref, pos_ecef);
	get_geomag_vector(2025.5, pos_ecef, mag_field);
	magField2Elements(mag_field, (float)pos_ref[0], (float)pos_ref[1], struct_mag_cal.Mag_vectors);

	// converting to uT
	// 41.6799011
	// 3.01447105
	// 24.4759007

	struct_mag_cal.Mag_vectors[0] *= 0.001f; // north
	struct_mag_cal.Mag_vectors[1] *= 0.001f; // east
	struct_mag_cal.Mag_vectors[2] *= 0.001f; // down
  
	for(i = 0; i<50; i++)
	{
		get_accel_mss(accel);
    get_gyro_rad(gyro);
		get_mag_ut(mag);
    run_alpha_lpf(&lpf_mag, mag, mag);
		get_initial_orientation(accel, mag, q);
      
		q_sum[0] += q[0];
		q_sum[1] += q[1];
		q_sum[2] += q[2];
    q_sum[3] += q[3];

		osDelay(10);
	}
	q_sum[0] *= 0.02f;
	q_sum[1] *= 0.02f;
	q_sum[2] *= 0.02f;
  q_sum[3] *= 0.02f;
 
  //mahony_init(q_sum);
  
  const float mNorm = norm_xyz(mag[0],mag[1],mag[2]);
  mMag_sf = struct_mag_cal.Mag_vectors[3]/fmaxf(mNorm,10.0f);
  
  insfilt_init(pos_ref, struct_mag_cal.Mag_vectors, q_sum);
  
	system_state.ekf_health = 1;
	set_target_course(struct_home_data.Heading_home);
	print_Debug("Origin LLA:\n");
	printf("%f\n%f\n%f\n", pos_ref[0], pos_ref[1], pos_ref[2]);
	print_Debug("North Magnetic field:\n");
	printf("%f\n%f\n%f\n", struct_mag_cal.Mag_vectors[0], struct_mag_cal.Mag_vectors[1], struct_mag_cal.Mag_vectors[2]);
  printf("Initial Heading = %f\n",struct_home_data.Heading_home);
	print_Debug("Estimatior Init Done\n");
}

__NO_RETURN void task_estimator(void *argument)
{
	(void)argument;
	static uint64_t ahrs_us, rate_mag_fusing, rate_gps_fusing;
  static uint32_t diag_ms;
  static float accel[3] __attribute__((section(".DTCM_RAM")));
  static float gyro[3] __attribute__((section(".DTCM_RAM")));
  static float mag[3] __attribute__((section(".DTCM_RAM")));
  static double gps_lla[3] __attribute__((section(".DTCM_RAM")));
	static float gps_vel[3] __attribute__((section(".DTCM_RAM")));
  float accel_mss[3];
  float gyro_rad[3];
	float dt = 1.0f / 200.0f;
	float R_pos;

	rst_flag = 0;
	dt_ahrs = dt;

	static insfilterAsync *m_ekf;
	static insfilterAsync_Cfg *m_cfg;
  
  m_ekf = get_ekf_instance();
  m_cfg = ekf_cfg_instance();

	while (!system_state.state_init)
	{
		ahrs_us = get_us64();
		rate_gps_fusing = get_us64();
		rate_mag_fusing = get_us64();
		// system_state.ekf_helth = 0;
		get_accel_mss(accel_mss);
		get_gyro_rad(gyro);
		get_mag_ut(mag);

    accel[0] = -accel_mss[0];
    accel[1] = accel_mss[1];
    accel[2] = accel_mss[2];
    
    run_alpha_lpf(&lpf_mag, mag, mag);
    run_alpha_lpf(&lpf_accel,accel,accel);
    run_alpha_lpf(&lpf_gyro,gyro,gyro);      
		osDelay(10);
	}

	while (1)
	{
		if (rst_flag)
		{
			rst_flag = 0;
			estimator_init();
		}

		if (system_state.ekf_health)
		{
			if (system_state.accl_health)
			{
				//in body frame

				get_accel_mss(accel_mss);
				get_gyro_rad(gyro_rad);
				
				float accel_rate[3];
				get_angular_accel(gyro_rad, accel_rate, dt_ahrs);
				// apply lever arm to cg in body frame
				imu_lever_arm_correction(accel_rate, gyro_rad, accel_mss);  

        accel[0] = accel_mss[0];
        accel[1] = accel_mss[1];
        accel[2] = accel_mss[2];
        
        gyro[0] = gyro_rad[0];
        gyro[1] = gyro_rad[1];
        gyro[2] = gyro_rad[2];
        
        run_alpha_lpf(&lpf_accel,accel,accel);
        run_alpha_lpf(&lpf_gyro,gyro,gyro);  
			}
			if (system_state.baro_health)
			{
				gps_lla[D] = get_baro_alt();
				gps_vel[D] = -get_alt_rate();
			}
			if (system_state.mag_health)
			{
				float mag_temp[3];
				get_mag_ut(mag_temp);
        
        mag[0] = mMag_sf*mag_temp[0];
        mag[1] = mMag_sf*mag_temp[1];
        mag[2] = mMag_sf*mag_temp[2];

        run_alpha_lpf(&lpf_mag, mag, mag);
			}
			if (system_state.gps_health)
			{
				R_pos = get_gps_Hdev();
				m_cfg->Rgps[0] = sqr(R_pos);
				m_cfg->Rgps[7] = sqr(R_pos);
        get_gps_vel(gps_vel);
				get_gps_pos(gps_lla);
			}
      
			//-------------------------------------------------------------------------------------
			const uint64_t now = get_us64();
			dt = (float)(now - ahrs_us) * 1e-6f;
			ahrs_us = now;
			dt_ahrs += 0.999f * (dt - dt_ahrs);
          
      update_process_noise();
      
      INS_Filt_predict(m_ekf,dt_ahrs);
              
      osDelay(1);
      
      INS_Filt_fuseaccel(m_ekf,accel);
      INS_Filt_fusegyro(m_ekf,gyro);
      
      quat2_rotm(m_ekf->pState, estimated_data.e_rotmat);
      quat2_eular(estimated_data.e_eular,estimated_data.e_rotmat);	
      
      if ((now - rate_mag_fusing) > mag_fusing_us) // 62hz
      {   
        rate_mag_fusing = now;
        INS_Filt_fusemag(m_ekf, mag);
      }
      if ((now - rate_gps_fusing) > gps_fusing_us) // 20hz
      {
        rate_gps_fusing = now;
        INS_Filt_fusegps(m_ekf,gps_lla,gps_vel);
      }
      if (zupt_detector(accel, gyro, gps_vel, dt, 1, 1)) 
      {
        INS_Filt_fusezupt(m_ekf);
      }
     
      osDelay(2);
			osMutexAcquire(ekf_mutex, 5);		   
      /*
             State22                            Units        State Vector Index
             Orientation as a quaternion                      1:4
             Position (NED)                      m            5:7
             Velocity (NED)                      m/s          8:10 
             Delta Angle Bias (XYZ)              rad          11:13
             Delta Velocity Bias (XYZ)           m/s          14:16
             Geomagnetic Field Vector (NED)      uT           17:19
             Magnetometer Bias (XYZ)             uT           20:22
      
        %    States28                          Units    Index
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

      estimated_data.quat[0] = m_ekf->pState[0];
      estimated_data.quat[1] = m_ekf->pState[1];
      estimated_data.quat[2] = m_ekf->pState[2];
      estimated_data.quat[3] = m_ekf->pState[3];
      estimated_data.e_angXYZ[0] = m_ekf->pState[4];
      estimated_data.e_angXYZ[1] = m_ekf->pState[5];
      estimated_data.e_angXYZ[2] = m_ekf->pState[6];
      estimated_data.e_posNED[0] = m_ekf->pState[7];
      estimated_data.e_posNED[1] = m_ekf->pState[8];
      estimated_data.e_posNED[2] = m_ekf->pState[9];
			estimated_data.e_velNED[0] = m_ekf->pState[10];
			estimated_data.e_velNED[1] = m_ekf->pState[11];
			estimated_data.e_velNED[2] = m_ekf->pState[12];
      estimated_data.e_accNED[0] = m_ekf->pState[13];
      estimated_data.e_accNED[1] = m_ekf->pState[14];
      estimated_data.e_accNED[2] = m_ekf->pState[15];
      
      //get_linear_ned_accel(accel, &m_ekf->pState[13], estimated_data.quat, estimated_data.e_accNED, dt_ahrs);

			double ned_pos[3] = {estimated_data.e_posNED[0], estimated_data.e_posNED[1], estimated_data.e_posNED[2]};
			ned2lla(ned_pos, m_ekf->ReferenceLocation, estimated_data.e_pos_lla);
			osMutexRelease(ekf_mutex);
		}
		osDelay(2);
	}
}

void NED_2Body(const float *ned, float *body)
{
	static float m[9];
	// Correct 3D Transformation (NED to Body)
	if (osMutexAcquire(ekf_mutex, 0) == osOK)
	{
		for (int i = 0; i < 9; i++)
			m[i] = estimated_data.e_rotmat[i];
		osMutexRelease(ekf_mutex);
	}
	body[0] = m[0] * ned[0] + m[3] * ned[1] + m[6] * ned[2];
	body[1] = m[1] * ned[0] + m[4] * ned[1] + m[7] * ned[2];
	body[2] = m[2] * ned[0] + m[5] * ned[1] + m[8] * ned[2];
}

void get_initial_orientation(const float *m_accl, const float *m_mag, float *quat)
{
	/*correct orientation*/
	/*ZYX frame*/
	//	roll right +ve
	//	roll left -ve
	//	yaw right +ve
	//	yaw left -ve
	//	pitch up +ve
	//	pitch dn -ve
	float eul[3]={0};
	
	const float pitch = atan2f(-m_accl[0], sqrtf(m_accl[1] * m_accl[1] + m_accl[2] * m_accl[2]));
	const float roll = atan2f(m_accl[1], m_accl[2]);

  // tilt compansested yaw;
  float cosRoll  = cosf(roll);
  float sinRoll  = sinf(roll);
  float cosPitch = cosf(pitch);
  float sinPitch = sinf(pitch);

  float mx_h = m_mag[0] * cosPitch + m_mag[2] * sinPitch;
  float my_h = m_mag[0] * sinRoll * sinPitch
               + m_mag[1] * cosRoll
               - m_mag[2] * sinRoll * cosPitch;

  const float yaw = atan2f(-my_h, mx_h);
  
  eul[0] = yaw;
  eul[1] = -pitch;
  eul[2] = roll;
 
	// ZYX
	eul2_quat(eul[2], eul[1], eul[0], quat);
	get_norm_quat(quat);

	struct_home_data.Heading_home = yaw;
	/*90 deg diff from mag ned to body frame*/

	if (struct_home_data.Heading_home < 0.0f)
		struct_home_data.Heading_home += mPi_2;

	struct_home_data.Heading_home = degrees_f(struct_home_data.Heading_home);
  //printf("Heading_home = %f\n",struct_home_data.Heading_home);
}

void gps_lever_arm_correction()
{
  
}

void get_linear_ned_accel(const float *accel_body, const float *dv, const float *quat, float *accel_ned, const float dt) 
{
    //remove bias
  const float dt_inv = 1.0f/dt;
  const float fbx = accel_body[0];  // + (dv[0]*dt_inv);
  const float fby = accel_body[1];  // + (dv[1]*dt_inv);
  const float fbz = accel_body[2];  // + (dv[2]*dt_inv);
    
  float u[3];
  float o_a;
  float o_b;
  float o_c;
  float q0[4];
  float qv_a;
  float qv_b;
  float qv_c;
  float qv_d;

  q0[0] = quat[0];
  q0[1] = quat[1];
  q0[2] = quat[2];
  q0[3] = quat[3];
  
  get_norm_quat(q0);

  qv_a = 0.0f;
  qv_b = fbx;
  qv_c = fby;
  qv_d = fbz;

  o_a = ((q0[0] * qv_a - q0[1] * qv_b) - q0[2] * qv_c) - q0[3] * qv_d;
  o_b = ((q0[0] * qv_b + q0[1] * qv_a) + q0[2] * qv_d) - q0[3] * qv_c;
  o_c = ((q0[0] * qv_c - q0[1] * qv_d) + q0[2] * qv_a) + q0[3] * qv_b;
  qv_a = ((q0[0] * qv_d + q0[1] * qv_c) - q0[2] * qv_b) + q0[3] * qv_a;
  q0[1] = -q0[1];
  q0[2] = -q0[2];
  q0[3] = -q0[3];
  u[0] = ((o_a * q0[1] + o_b * q0[0]) + o_c * q0[3]) - qv_a * q0[2];
  u[1] = ((o_a * q0[2] - o_b * q0[3]) + o_c * q0[0]) + qv_a * q0[1];
  u[2] = ((o_a * q0[3] + o_b * q0[2]) - o_c * q0[1]) + qv_a * q0[0];
 
  accel_ned[0] = u[0];
  accel_ned[1] = u[1];
  accel_ned[2] = u[2] - Actual_gravity;
}

void get_est_quat(float *q)
{
  q[0] = estimated_data.quat[0];
  q[1] = estimated_data.quat[1];
  q[2] = estimated_data.quat[2];
  q[3] = estimated_data.quat[3];
}
float get_est_roll(void)
{
	static float roll;
	if (osMutexAcquire(ekf_mutex, 0) == osOK)
	{
		roll = degrees_f(estimated_data.e_eular[0]);
		osMutexRelease(ekf_mutex);
	}
	return roll;
}
float get_est_pitch(void)
{
	static float pitch;
	if (osMutexAcquire(ekf_mutex, 0) == osOK)
	{
		pitch = degrees_f(estimated_data.e_eular[1]);
		osMutexRelease(ekf_mutex);
	}
	return pitch;
}
float get_est_yaw(void)
{
	static float yaw;
	if (osMutexAcquire(ekf_mutex, 0) == osOK)
	{
    yaw = estimated_data.e_eular[2];
    if(yaw < 0.0f)
      yaw += mPi_2;
    
		yaw = degrees_f(yaw);
		osMutexRelease(ekf_mutex);
	}
	return yaw;
}
double get_est_latD(void)
{
	static double lat;
	if (osMutexAcquire(ekf_mutex, 0) == osOK)
	{
		lat = estimated_data.e_pos_lla[0];
		osMutexRelease(ekf_mutex);
	}
	return lat;
}
double get_est_lonD(void)
{
	static double lon;
	if (osMutexAcquire(ekf_mutex, 0) == osOK)
	{
		lon = estimated_data.e_pos_lla[1];
		osMutexRelease(ekf_mutex);
	}
	return lon;
}
float get_est_alt_m(void)
{
	static float altm;
	if (osMutexAcquire(ekf_mutex, 0) == osOK)
	{
		altm += 0.75f * (estimated_data.e_pos_lla[2] - altm);
		osMutexRelease(ekf_mutex);
	}
	return altm;
}
float get_est_posN(void)
{
	static float posN;
	if (osMutexAcquire(ekf_mutex, 0) == osOK)
	{
		posN += 0.75f * (estimated_data.e_posNED[0] - posN);
		osMutexRelease(ekf_mutex);
	}
	return posN;
}
float get_est_posE(void)
{
	static float posE;
	if (osMutexAcquire(ekf_mutex, 0) == osOK)
	{
		posE += 0.75f * (estimated_data.e_posNED[1] - posE);
		osMutexRelease(ekf_mutex);
	}
	return posE;
}

float get_est_posD(void)
{
	static float posD;
	if (osMutexAcquire(ekf_mutex, 0) == osOK)
	{
    posD += 0.75f * (estimated_data.e_posNED[2] - posD);
		osMutexRelease(ekf_mutex);
	}
	return posD;
}

float get_est_velN(void)
{
	static float velN;
	if (osMutexAcquire(ekf_mutex, 0) == osOK)
	{
		velN += 0.65f * (estimated_data.e_velNED[N] - velN);
		osMutexRelease(ekf_mutex);
	}
	return velN;
}
float get_est_velE(void)
{
	static float velE;
	if (osMutexAcquire(ekf_mutex, 0) == osOK)
	{
		velE += 0.65f * (estimated_data.e_velNED[E] - velE);
		osMutexRelease(ekf_mutex);
	}
	return velE;
}
float get_est_velU(void)
{
	static float velU;
	if (osMutexAcquire(ekf_mutex, 0) == osOK)
	{
    const float velu = -estimated_data.e_velNED[D];
		velU += 0.65f * (velu - velU);
		osMutexRelease(ekf_mutex);
	}
	return velU;
}

float get_est_velR(void)
{
	static float vel_r;
	if (osMutexAcquire(ekf_mutex, 0) == osOK)
	{
		float brg = 0;
		brg = degrees_f(estimated_data.e_eular[2]);
		float v_north = estimated_data.e_velNED[0] * cosf(brg);
		float v_east = estimated_data.e_velNED[1] * sinf(brg);
		vel_r = sqrtf(fabsf(v_north * v_north) + fabsf(v_east * v_east));
		osMutexRelease(ekf_mutex);
	}
	return vel_r;
}

float get_est_acclN(void)
{
	static float accel_outN;
	if (osMutexAcquire(ekf_mutex, 0) == osOK)
	{
		accel_outN += 0.72f * (estimated_data.e_accNED[N] - accel_outN);
		osMutexRelease(ekf_mutex);
	}
	return accel_outN;
}
float get_est_acclE(void)
{
	static float accel_outE;
	if (osMutexAcquire(ekf_mutex, 0) == osOK)
	{
		accel_outE += 0.72f * (estimated_data.e_accNED[E] - accel_outE);
		osMutexRelease(ekf_mutex);
	}
	return accel_outE;
}
float get_est_acclD(void)
{
	static float accel_outD;
	if (osMutexAcquire(ekf_mutex, 0) == osOK)
	{
		accel_outD += 0.72f * (estimated_data.e_accNED[D] - accel_outD);
		osMutexRelease(ekf_mutex);
	}
	return accel_outD;
}
float get_est_turnR()
{
	static float turnRate;
	if (osMutexAcquire(ekf_mutex, 0) == osOK)
	{
		turnRate = estimated_data.e_turnRate;
		osMutexRelease(ekf_mutex);
	}
	return turnRate;
}
float get_est_windN(void)
{
	return 0;
}
float get_est_windE(void)
{

	return 0;
}
