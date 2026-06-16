/**
 * @file vel_control.c
 * @brief Velocity control implementation for the navigation system.
 *
 * This file contains the implementation of velocity control algorithms, including filtering and PID control.
 * 
 * @author saikat
 * @date 2023-12-01
 * @version 1.0
 */
#include <vel_control.h>
#include <math.h>
#include <util.h>
#include <flight_control.h>
#include <nav_sys_config.h>
#include <system.h>
#include <nav.h>
#include <stdio.h>
#include <estimator.h>
#include <calibration.h>
#include <imu.h>
#include <drv_utils.h>
#include <pos_control.h>

typedef struct f
{
	float m_vel;
	float alpha;
	
}filt_alpha;

#define velNE_ff_gain 0.05f
#define velD_ff_gain 0.15f

//filt_alpha filt_velD;

extern float dbg_var1,dbg_var2,dbg_var3,dbg_var4;

static float integral_anti_windup(float *i_term_state,float error,float Ki,float non_i_term,float out_limit,float dt);

static inline float vnorm3_(const float v[3]) 
{
  return sqrtf(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
}
static inline void vscale3_(float v[3], float s) 
{ 
  v[0]*=s; v[1]*=s; v[2]*=s; 
}

void vel_controller_init()
{
	// Vertical velocity loop (1-3 kg multirotor starting point):
	// lower Kp/Kd to reduce overshoot, higher Ki to hold altitude (no drift).

  PID_config.PID_Ralt_Kp = 2.5f;
  PID_config.PID_Ralt_Ki = 0.15f;
  PID_config.PID_Ralt_Kd = 1.2f;


	// Horizontal velocity loop: reduce Kp (less bounce) and Kd (D is on noisy
	// measured accel), raise Ki to kill steady-state velocity/position drift.
	PID_config.PID_Rpos_Kp = 1.80f;
	PID_config.PID_Rpos_Ki = 0.005f;
	PID_config.PID_Rpos_Kd = 0.425f;
	
//	filt_velD.alpha = 0.005f;
}

void vel_controller_reset()
{ 
  control_data.integral_velN = 0.0f;
  control_data.integral_velE = 0.0f;
  control_data.integral_velD = 0.0f;
  
  control_data.pid_out_velN = 0.0f;
  control_data.pid_out_velE = 0.0f;
  control_data.pid_out_velD = 0.0f; 
}


void velocity_2lean_angle(float vel)
{
		if(vel>2.0f && get_est_turnR()>5.0f)
		{
			control_data.lean_angle = f_atan2f(vel*vel,Actual_gravity-max_turn_redius);
			control_data.lean_angle = degrees_f(control_data.lean_angle);
		}
		else
			control_data.lean_angle = 0.0f;
}

void update_velocity_setpoints(float forward_speed_cmd, float lateral_speed_cmd)
{
    // 1. Get current heading from EKF
    float yaw = radians_f(get_est_yaw()); 

 // from controller 
    control_data.target_velN_cmd = (forward_speed_cmd * cosf(yaw) - lateral_speed_cmd * sinf(yaw))*Gravity_inv;
    control_data.target_velE_cmd = (forward_speed_cmd * sinf(yaw) + lateral_speed_cmd * cosf(yaw))*Gravity_inv;
}

float integral_anti_windup(float *i_term_state,float error,float Ki, float non_i_term,float out_limit,float dt)
{
    const float Kaw = 2.5f;
    float i_term = *i_term_state;
  
  //no integrator if error is more than 2
//    if(fabsf(error)>2.5f)
//    {
//      *i_term_state = 0.0f;
//      return 0.0f;
//    }
//  
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

/*
Vertical velocity controller
*/

void run_velD_controller(const float target_velD, float dt)
{
	float velU_pid_out = 0.0f;
	static float last_target_vz;
	static float filt_accel_ff;
	const float Kp = PID_config.PID_Ralt_Kp;
	const float Ki = PID_config.PID_Ralt_Ki;
	const float Kd = PID_config.PID_Ralt_Kd;
	
  //in NEU frame,so Down acceleration and Down velocity need to be negetive
	const float accD_est = -get_est_acclD();
	float velD_est = get_est_velU();
  
  dbg_var1 = target_velD;
  //alt rate PID 
	float error = target_velD - velD_est;
  
	velU_pid_out = Kp*error; 
 	velU_pid_out += integral_anti_windup(&control_data.integral_velD,error,Ki,velU_pid_out,A_Integral_Lim,dt);
	velU_pid_out -= (Kd*accD_est);
    
  control_data.pid_out_posD = constrain_f(velU_pid_out,-A_POSPIDout_Lim,A_POSPIDout_Lim);	  
}

/*
direction    |
controller   |
____\/_______|

position	
   \/
velocity 		
   \/
accelaration				
   \/
NED to body translation						

roll command , pitch command

//  accel east +ve
//  accel west -ve
//  accel north +ve
//  accel south -ve
//  accel up +ve
//  accel dn -ve

//	vel east +ve
//	vel west -ve
//	vel north +ve
//	vel south -ve
//	vel up +ve
//	vel dn -ve
  
//	pos east +ve
//	pos west -ve
//	pos north +ve
//	pos south -ve
//	pos up +ve
//	pos dn -ve  

//	roll right +ve  ||   left vel +ve roll
//	roll left -ve
//	pitch up -ve    || fw vel -pitch
//	pitch dn +ve
//	yaw clockwise +ve
//	yaw anticlockwise -ve
********************************/

void accel_2lean_angle(const float *a_des_ned, float yaw_des_deg, float *euler_deg)
{
    const float e3[3] = {0.0f, 0.0f, 1.0f};
    static float body_rate[3];

    // 1) Desired body Z (down) in NED: z_b = normalize(g*e3 - a_des)
    float zb[3] = {
        Actual_gravity*e3[0] - a_des_ned[0],
        Actual_gravity*e3[1] - a_des_ned[1],
        Actual_gravity*e3[2] - a_des_ned[2]
    };
    float zbn = vnorm3_(zb);
    const float zbn_min = 0.2f * Actual_gravity;          // floor to avoid singularity near free fall
    if (zbn < zbn_min) zbn = zbn_min;
    vscale3_(zb, 1.0f / zbn);

    // 2) Desired heading basis in NED (horizontal)
    const float yaw = radians_f(yaw_des_deg);
    float xc[3] = { cosf(yaw), sinf(yaw), 0.0f };

    // 3) Project xc onto plane orthogonal to z_b to get x_b (in NED)
    float dot = xc[0]*zb[0] + xc[1]*zb[1] + xc[2]*zb[2];
    float xb[3] = { xc[0] - dot*zb[0], xc[1] - dot*zb[1], -dot*zb[2] };
    float xb_n = vnorm3_(xb);
    if (xb_n < 1e-6f) {  // degenerate: z_b aligned with xc
        xb[0] = 1.0f; xb[1] = 0.0f; xb[2] = 0.0f;
    } else {
        vscale3_(xb, 1.0f / xb_n);
    }

    // 4) y_b = z_b × x_b
    float yb[3] = {
        zb[1]*xb[2] - zb[2]*xb[1],
        zb[2]*xb[0] - zb[0]*xb[2],
        zb[0]*xb[1] - zb[1]*xb[0]
    };

    // 5) Rotation with columns as body axes in NED (R_wb)
    float R[9] = {
        xb[0], yb[0], zb[0],
        xb[1], yb[1], zb[1],
        xb[2], yb[2], zb[2]
    };

    // 6) Extract ZYX Euler from R_wb (roll, pitch, yaw)
    // Clamp for numerical safety
    float r21 = constrain_f(R[7], -1.0f, 1.0f);
    float r22 = constrain_f(R[8], -1.0f, 1.0f);
    float r20 = constrain_f(R[6], -1.0f, 1.0f);
    float r10 = constrain_f(R[3], -1.0f, 1.0f);
    float r00 = constrain_f(R[0], -1.0f, 1.0f);

    const float phi   = atan2f(r21, r22);   // roll
    const float theta = -asinf(r20);        // pitch
    const float psi   = atan2f(r10, r00);   // yaw

    // Optional tilt limit (e.g., 35 deg)
    const float max_tilt_deg = ctrl_config.s_div;
    float phi_d  = degrees_f(phi);
    float theta_d= degrees_f(theta);
    float tilt_d = fmaxf(fabsf(phi_d), fabsf(theta_d));
    if (tilt_d > max_tilt_deg) {
        float s = max_tilt_deg / tilt_d;
        phi_d   *= s;
        theta_d *= s;
    }
    
    body_rate[0] += 0.72f*(phi_d - body_rate[0]);
    body_rate[1] += 0.72f*(theta_d - body_rate[1]);
    body_rate[2] += 0.72f*(psi - body_rate[2]);
         
//    filt_vel_e = body_rate[0];
//    filt_vel_n = body_rate[1];
    
    euler_deg[Roll]  = body_rate[0];
    euler_deg[Pitch] = -body_rate[1];
    euler_deg[Yaw]   = degrees_f(body_rate[2]);
    
//    if(console_buf[1]=='n')
//      printf("%f,%f,%f\n",euler_deg[Roll],euler_deg[Pitch],euler_deg[Yaw]);
}

void run_velNE_controller(float target_velN, float target_velE, float dt)
{
	float error_n,error_e;
	float vel_n,vel_e;
  float accel_n,accel_e;
	float pid_velocty[2] = {0,0};
	float ff_accel_n = 0.0f;
	float ff_accel_e = 0.0f;
	static float filt_accel_n,filt_accel_e;
	static float last_target_velN,last_target_velE;

	const float Kp = PID_config.PID_Rpos_Kp;
	const float Ki = PID_config.PID_Rpos_Ki;
	const float Kd = PID_config.PID_Rpos_Kd;

	vel_n = get_est_velN();
	vel_e = get_est_velE();
  	
	accel_n = get_est_acclN();
	accel_e = get_est_acclE();
  
	error_n = target_velN - vel_n;
	error_e = target_velE - vel_e;
    
	/*Propotional + Integral*/
	pid_velocty[N] = error_n * Kp;
  pid_velocty[E] = error_e * Kp;
    
    //float integral_anti_windup(float *i_term_state,float error,float Ki,float non_i_term,float out_limit,float dt)

  if(fabsf(error_n)>0.25f)
    control_data.integral_velN = 0.0f;
  
  if(fabsf(error_e)>0.25f)
    control_data.integral_velE = 0.0f;
  
  pid_velocty[N] += integral_anti_windup(&control_data.integral_velN,error_n,Ki,pid_velocty[N],P_Integral_Lim,dt);
  pid_velocty[E] += integral_anti_windup(&control_data.integral_velE,error_e,Ki,pid_velocty[E],P_Integral_Lim,dt);
  
	pid_velocty[N] -= (accel_n * Kd);
	pid_velocty[E] -= (accel_e * Kd);
  //e = +ve
  //n = +ve
   
  float momentum_ff_n = (target_velN - last_target_velN) / dt;
  last_target_velN = target_velN;
  
  float momentum_ff_e = (target_velE - last_target_velE) / dt;
  last_target_velE = target_velE;
  
  filt_accel_n += 0.82f * (momentum_ff_n - filt_accel_n);
  filt_accel_e += 0.82f * (momentum_ff_e - filt_accel_e);
       
  ff_accel_n = (ctrl_config.v_mass * filt_accel_n) * Gravity_inv;
  ff_accel_e = (ctrl_config.v_mass * filt_accel_e) * Gravity_inv;
        
  //out pid from velocity ne controller 
  control_data.pid_out_posN = pid_velocty[N]+(ff_accel_n*velNE_ff_gain);
  control_data.pid_out_posE = pid_velocty[E]+(ff_accel_e*velNE_ff_gain);
  
//  filt_vel_n = control_data.pid_out_pos[N];
//  filt_vel_e = control_data.pid_out_pos[E];
}

