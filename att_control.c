/**
 * @file att_control.c
 * @brief Attitude control implementation.
 *
 * This file contains the implementation of the attitude controller, including angle and rate control loops.
 *
 * @author Saikat
 * @date 2023-11-10
 * @version 1.0
 */
#include <att_control.h>
#include <math.h>
#include <util.h>
#include <flight_control.h>
#include <nav_sys_config.h>
#include <system.h>
#include <nav.h>
#include <stdio.h>
#include <imu.h>
#include <estimator.h>

static float ExpKp(float k,float error);
static float ExpKd(float rate);
static float integral_anti_windup(float *i_term_state,float error,float Ki,float non_i_term,float out_limit,float dt);
static float filt_derivitive[3];
static float I_lim[3];

#define pid_kb 0.5f

float ExpKp(float k,float error)
{
 return k *(1 + P_Constant*error*error);
}

float ExpKd(float rate)
{
  return rate *(1 + D_Constant*rate*rate);
}
#include <stdbool.h>
#include <math.h>

float integral_anti_windup(float *i_term_state,float error,float Ki,float non_i_term,float out_limit,float dt)
{
    const float Kaw = 0.5f;
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

void att_controller_init()
{	
	/*pid init */
  
  I_lim[0] = 50.0f;
  I_lim[1] = 50.0f;
  I_lim[2] = 50.0f;
  
	PID_config.PI_Sroll_Kp = 4.65f;      
	PID_config.PI_Sroll_Ki = 0.85f;    
																					
	PID_config.PI_Spitch_Kp = 4.65f;     
	PID_config.PI_Spitch_Ki = 0.85f;   
																					
	PID_config.PI_Syaw_Kp = 4.25f;       
	PID_config.PI_Syaw_Ki = 0.85f;      
																					
	PID_config.PID_Rroll_Kp = 1.2f;     
	PID_config.PID_Rroll_Ki = 0.05f;   
	PID_config.PID_Rroll_Kd = 0.045f;   
																					
	PID_config.PID_Rpitch_Kp = 1.2f;    
	PID_config.PID_Rpitch_Ki = 0.05f;  
	PID_config.PID_Rpitch_Kd = 0.045f;
	
	PID_config.PID_Ryaw_Kp = 1.2f;
	PID_config.PID_Ryaw_Ki = 0.05f;
	PID_config.PID_Ryaw_Kd = 0.05f;
}

void att_controller_reset()
{
  memset(filt_derivitive,0,sizeof(filt_derivitive)); 
	memset((void*)control_data.set_point_att,0,sizeof(control_data.set_point_att));
	memset((void*)control_data.s_integral_att,0,sizeof(control_data.s_integral_att));
	memset((void*)control_data.r_integral_att,0,sizeof(control_data.r_integral_att));
	memset((void*)control_data.s_pid_out_att,0,sizeof(control_data.s_pid_out_att));
  memset((void*)control_data.r_pid_out_att,0,sizeof(control_data.r_pid_out_att));
}

void run_angle_controller(float dt)
{
  float error[3] = {0,0,0};
  float pi_out[3] = {0,0,0};
  float cfg_turn_rate = (float)ctrl_config.rate_turn;  //deg / sec

  float est_yaw = get_est_yaw();
  float est_pitch = get_est_pitch();
  float est_roll = get_est_roll();
  
	error[Yaw] = control_data.set_point_att[Yaw] - est_yaw;
	if(error[Yaw]<=-180.0f) error[Yaw] += 360.0f;
	if(error[Yaw]>= 180.0f) error[Yaw] -= 360.0f;

  if(error[Yaw] > cfg_turn_rate)
    error[Yaw] = cfg_turn_rate;
  if(error[Yaw] < -cfg_turn_rate)
    error[Yaw] = -cfg_turn_rate;
      
	error[Roll] = control_data.set_point_att[Roll] - est_roll; // * 1.5f; //sub deg presition
	error[Pitch] = control_data.set_point_att[Pitch] - est_pitch;
	
	pi_out[Roll]  = error[Roll] * PID_config.PI_Sroll_Kp; 
	pi_out[Pitch] = error[Pitch] * PID_config.PI_Spitch_Kp;
	pi_out[Yaw]   = error[Yaw] * PID_config.PI_Syaw_Kp; 
			//float integral_anti_windup(float *i_term_state,float error,float Ki,float non_i_term,float out_limit,float dt)
	pi_out[Roll] += integral_anti_windup((float*)&control_data.s_integral_att[Roll],error[Roll],PID_config.PI_Sroll_Ki,pi_out[Roll],I_lim[0], dt);
  pi_out[Pitch]+= integral_anti_windup((float*)&control_data.s_integral_att[Pitch],error[Pitch],PID_config.PI_Spitch_Ki,pi_out[Pitch],I_lim[1], dt);
  pi_out[Yaw]  += integral_anti_windup((float*)&control_data.s_integral_att[Yaw],error[Yaw],PID_config.PI_Syaw_Ki,pi_out[Yaw],I_lim[2], dt);
  	
  pi_out[Roll]  += (control_data.set_point_att[Roll]*K_ff);   //feed forward
  pi_out[Pitch] += (control_data.set_point_att[Pitch]*K_ff);
  pi_out[Yaw]   += pi_out[Yaw];

  control_data.s_pid_out_att[Roll]  = constrain_f(pi_out[Roll],-Att_PIDout_Lim,Att_PIDout_Lim);
  control_data.s_pid_out_att[Pitch] = constrain_f(pi_out[Pitch],-Att_PIDout_Lim,Att_PIDout_Lim);
  control_data.s_pid_out_att[Yaw]   = constrain_f(pi_out[Yaw],-Att_PIDout_Lim,Att_PIDout_Lim);

}

void run_rate_controller(float dt)
{
	static float r_last_gyro[3];

	float derivative[3]={0,0,0};
  float gyro[3];
	float r_error_att[3];
	float pid_out[3]={0,0,0};
	
	get_gyro_deg(gyro);

	r_error_att[Roll]  = control_data.s_pid_out_att[Roll] - gyro[X];
	r_error_att[Pitch] = control_data.s_pid_out_att[Pitch] - gyro[Y];
	r_error_att[Yaw]   = control_data.s_pid_out_att[Yaw] - gyro[Z];
	
//P
	pid_out[Roll]  += r_error_att[Roll]  * PID_config.PID_Rroll_Kp;
	pid_out[Pitch] += r_error_att[Pitch] * PID_config.PID_Rpitch_Kp;
	pid_out[Yaw]   += r_error_att[Yaw]   * PID_config.PID_Ryaw_Kp;
		
//I
  if(r_error_att[Roll] < 1.25f){    	
    control_data.r_integral_att[Roll]	+= (r_error_att[Roll]*dt*dt);
    pid_out[Roll]  += control_data.r_integral_att[Roll]*PID_config.PID_Rroll_Ki;
  }
  if(r_error_att[Pitch] < 1.25){
      control_data.r_integral_att[Pitch]+= (r_error_att[Pitch]*dt*dt);
      pid_out[Pitch] += control_data.r_integral_att[Pitch]*PID_config.PID_Rpitch_Ki;
  }
  if(r_error_att[Yaw] < 1.25f){
    control_data.r_integral_att[Yaw]	+= (r_error_att[Yaw]*dt*dt);
    pid_out[Yaw] += control_data.r_integral_att[Yaw]*PID_config.PID_Ryaw_Ki;
  }
  
//D
  float inv_dt = 1.0f/fmaxf(dt, 1e-4f);
	derivative[Roll]  = (r_error_att[Roll]  - r_last_gyro[Roll])*inv_dt; 
	derivative[Pitch] = (r_error_att[Pitch] - r_last_gyro[Pitch])*inv_dt; 
	derivative[Yaw]   = (r_error_att[Yaw]   - r_last_gyro[Yaw])*inv_dt; 
	
	r_last_gyro[Roll]	= r_error_att[Roll];
	r_last_gyro[Pitch]= r_error_att[Pitch];
	r_last_gyro[Yaw]	= r_error_att[Yaw];
  	   
  filt_derivitive[Roll]  = 0.86f*filt_derivitive[Roll] + derivative[Roll] * 0.14f;
	filt_derivitive[Pitch] = 0.86f*filt_derivitive[Pitch]+ derivative[Pitch]* 0.14f;
	filt_derivitive[Yaw]   = 0.86f*filt_derivitive[Yaw]  + derivative[Yaw]  * 0.14f;
    
	pid_out[Roll]  += (filt_derivitive[Roll] * PID_config.PID_Rroll_Kd);
	pid_out[Pitch] += (filt_derivitive[Pitch]* PID_config.PID_Rpitch_Kd);
	pid_out[Yaw]   += (filt_derivitive[Yaw]  * PID_config.PID_Ryaw_Kd);
	
	control_data.r_pid_out_att[Roll]  = constrain_f(pid_out[Roll], -Att_PIDout_Lim, Att_PIDout_Lim);
	control_data.r_pid_out_att[Pitch] = constrain_f(pid_out[Pitch],-Att_PIDout_Lim, Att_PIDout_Lim);
	control_data.r_pid_out_att[Yaw]   = constrain_f(pid_out[Yaw],  -Att_PIDout_Lim, Att_PIDout_Lim);
}

