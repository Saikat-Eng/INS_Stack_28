/**
 * @file pos_control.c
 * @brief Position controller implementation.
 *
 * This file contains the position controller implementation for the navigation system.
 *
 * @author saikat
 * @date 2023-12-01
 * @version 1.0
 */
#include <pos_control.h>
#include <math.h>
#include <util.h>
#include <flight_control.h>
#include <nav_sys_config.h>
#include <system.h>
#include <nav.h>
#include <stdio.h>
#include <estimator.h>
#include <vel_control.h>
#include <INSEstimator.h>

static unsigned char alt_lock = 0, pos_lock = 0, course_lock = 0, circle_lock = 0;
//static float integral_anti_windup(float *i_term_state,float error,float Ki,float non_i_term,float out_limit,float dt);

void pos_controller_init()
{
	/*Stable controller*/
  PID_config.PI_Salt_Kp = 0.35f;	
  PID_config.PI_Spos_Kp = 0.0001f;
    
  PID_config.PI_Spos_Ki = 0.085f;
  PID_config.PI_Salt_Ki = 0.05f;
    
	set_target_passive_alt(Min_ground_alt);
}

void pos_controller_reset()
{
  control_data.integral_posN = 0.0f;
  control_data.integral_posE = 0.0f;
  control_data.integral_posD = 0.0f;
  
  control_data.pid_out_posN = 0.0f;
  control_data.pid_out_posE = 0.0f;
  control_data.pid_out_posD = 0.0f;  
}

inline void set_alt_lock(unsigned char lock)
{
	control_data.pid_out_velD  = 0.0f;
	control_data.integral_velD = 0.0f;
	control_data.integral_posD = 0.0f;
	alt_lock = lock;
}
inline void set_pos_lock(unsigned char lock)
{
	pos_lock = lock;
}
inline void set_course_lock(unsigned char lock)
{
	course_lock = lock;
}
inline void set_circle_lock(unsigned char lock)
{
	circle_lock = lock;
}

inline unsigned char get_alt_lock(void)
{
	return alt_lock;
}
inline unsigned char get_pos_lock(void)
{
	return pos_lock;
}
inline unsigned char get_course_lock(void)
{
	return course_lock;
}
inline unsigned char get_circle_lock(void)
{
	return circle_lock;
}
/*
set target altitude from baro altitude
*/
char set_target_active_alt(const float alt)
{
	static float alt_avg = 0;
	static int idx = 0;
	alt_avg += 0.72f * (alt - alt_avg);
	idx++;									
	
	if((idx>20))
	{	
		idx = 0;
		//if(control_data.set_point_pwm[THR]>RC_min || (fabsf(control_data.target_alt - alt_avg)<Min_ground_alt))
		{
			//if diffrence is < Min_ground_alt, set target alt from baro alt
			control_data.target_alt = alt_avg;
		}
		control_data.ground_alt = alt_avg;
    alt_avg = 0.0f;
		set_alt_lock(1);
		return 1;
	}
	return 0;
}

/*
set target altitude from other source
*/

char set_target_passive_alt(const float alt)
{
	if(alt>=Min_ground_alt)
	{
		control_data.target_alt = alt;
		set_alt_lock(1);
		return 1;
	}
	return 0;
}

char set_target_active_pos(const float pos_n, const float pos_e)
{
	static float pos_avg[2];
	static int idx = 0;
  
	pos_avg[0] += 0.72f * (pos_n - pos_avg[0]);
  pos_avg[1] += 0.72f * (pos_e - pos_avg[1]);
  
	idx++;									
	
	if(idx>20)
	{	
    set_zupt_pos_hold(1);
		idx = 0;
		control_data.target_posN = 	pos_avg[0];
    control_data.target_posE = 	pos_avg[1];
    pos_avg[0] = 0.0f;
    pos_avg[1] = 0.0f;
		set_pos_lock(1);
		return 1;
	}
	return 0;
}

float integral_anti_windup(float *i_term_state,float error,float Ki,float non_i_term,float out_limit,float dt)
{
    const float Kaw = 2.0f;
    float i_term = *i_term_state;

    float unsat = non_i_term + i_term;

    float sat = unsat;

    if (out_limit > 0.0f) {
        if (sat > out_limit) {
            sat = out_limit;
        } else if (sat < -out_limit) {
            sat = -out_limit;
        }
    }
    
    bool saturated_high = false;
    bool saturated_low  = false;

    if (out_limit > 0.0f) {
        saturated_high = unsat > out_limit;
        saturated_low  = unsat < -out_limit;
    }

    float error_for_integration = error;

    if (saturated_high && error > 0.0f) {
        error_for_integration = 0.0f;
    } else if (saturated_low && error < 0.0f) {
        error_for_integration = 0.0f;
    }

    i_term += (Ki * error_for_integration + Kaw * (sat - unsat)) * dt;

    if (out_limit > 0.0f) {
        if (i_term > out_limit) {
            i_term = out_limit;
        } else if (i_term < -out_limit) {
            i_term = -out_limit;
        }
    }
    *i_term_state = i_term;

    return i_term;
}


void run_pos_controller(const float lat_diff, const float lon_diff,const float dt)
{
	const float DEADZONE_M = 0.005f;   //0.5cm       /* position deadband (meters) */
	const float MAX_POS_VEL = ctrl_config.rate_Hspd;

	const float Kp = PID_config.PI_Spos_Kp;
  const float Ki = PID_config.PI_Spos_Ki;
  
  float err_n = 0.0f;
  if (fabsf(lat_diff) > DEADZONE_M) {
        err_n = (lat_diff > 0) ? (lat_diff - DEADZONE_M) : (lat_diff + DEADZONE_M);
  }
    
  float err_e = 0.0f;
  if (fabsf(lon_diff) > DEADZONE_M) {
      err_e = (lon_diff > 0) ? (lon_diff - DEADZONE_M) : (lon_diff + DEADZONE_M);
  }
    /*P*/
	float target_vel_n = Kp * err_n;
	float target_vel_e = Kp * err_e;

 // float integral_anti_windup(float *i_term_state,float error,float Ki,float non_i_term,float out_limit,float dt)
	target_vel_n += integral_anti_windup(&control_data.integral_posN, err_n, Ki, target_vel_n, MAX_POS_VEL, dt);
	target_vel_e += integral_anti_windup(&control_data.integral_posE, err_e, Ki, target_vel_e, MAX_POS_VEL, dt);
	/* P + I */

	control_data.pid_out_velN = constrain_f(target_vel_n, -MAX_POS_VEL, MAX_POS_VEL);
	control_data.pid_out_velE = constrain_f(target_vel_e, -MAX_POS_VEL, MAX_POS_VEL);
}

void run_alt_controller(const float alt_diff, const float dt)
{
    const float Kp = PID_config.PI_Salt_Kp;
    const float Ki = PID_config.PI_Salt_Ki;
    
    const float max_velU = ctrl_config.rate_Vspd;

    float velU_sp = Kp * alt_diff;
  
    velU_sp += integral_anti_windup(&control_data.integral_posD,
                                    alt_diff,
                                    Ki,
                                    velU_sp,
                                    max_velU,
                                    dt);

    control_data.pid_out_velD = constrain_f(velU_sp, -max_velU, max_velU);
}
