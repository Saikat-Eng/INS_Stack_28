/**
 * @file flight_control.c
 * @brief Implementation of the flight control module for the NAV_SYS_M7 project.
 *
 * This file contains the implementation of the flight control module, which manages
 * the control inputs and outputs for the aircraft.
 *
 * @author saikat
 * @date 2023-11-30
 * @version 1.0
 */
#include "stm32f7xx_hal.h" // Keil::Device:STM32Cube HAL:Common
#include "cmsis_os2.h"
#include <drv_utils.h>
#include <gpio.h>
#include <pwm.h>
#include <timer.h>
#include <math.h>
#include <util.h>
#include <flight_control.h>
#include <nav_sys_config.h>
#include <system.h>
#include <nav.h>
#include <imu.h>
#include <task.h>
#include <telemetry.h>
#include <e2prom.h>
#include <combiner.h>
#include <calibration.h>
#include <estimator.h>
#include <att_control.h>
#include <pos_control.h>
#include <vel_control.h>
#include <thr_control.h>
#include <wp_control.h>
#include <GeoBase.h>
#include <INSEstimator.h>

struct_input ext_control __attribute__((section(".DTCM_RAM")));
struct_PID PID_config __attribute__((section(".DTCM_RAM")));
control_params ctrl_config __attribute__((section(".DTCM_RAM")));
volatile struct_control control_data __attribute__((section(".DTCM_RAM")));
static float idle_thr;
static float delta_f;
static hpf_3D lpf_att;
static hpf_1D lpf_alt;
static hpf_1D lpf_pos[2];
static char arming_init;
volatile float velocityN_sp, velocityE_sp, velocityD_sp;
volatile float dt_control;
volatile float dt_process;
volatile float dt_navctrl;
volatile float dt_stable_pid;



void Update_Controls(const float thr)
{
	float att_out[3] = {0.0f, 0.0f, 0.0f};
  doResonantLPF(&lpf_att, control_data.r_pid_out_att, att_out);
  
#ifdef Heli

	f[0] = _thr;
	f[5] = PID_OUT[Yaw];

	f[1] = PID_OUT[Pitch] + PID_OUT[Roll];
	f[2] = PID_OUT[Pitch] - PID_OUT[Roll];
	f[3] = PID_OUT[Pitch] - PID_OUT[Roll];
	f[4] = PID_OUT[Pitch] - PID_OUT[Roll];

	Mt1 = constrain_f(pulse_out[0], RC_min, RC_max);
	Mt2 = constrain_f(pulse_out[1], RC_min, RC_max);
	Mt3 = constrain_f(pulse_out[2], RC_min, RC_max);
	Mt4 = constrain_f(pulse_out[3], RC_min, RC_max);
	Mt5 = constrain_f(pulse_out[4], RC_min, RC_max);
	Mt6 = constrain_f(pulse_out[5], RC_min, RC_max);
#endif

#ifdef Quad
/*
		roll right +ve
		roll left -ve
		yaw right +ve
		yaw left -ve
		pitch up -ve
		pitch dn +ve

	  accel up -ve
	  accel dn +ve
	  accel left -ve
	  accel right +ve
	  accel forward -ve
	  accel backward +ve
*/

	static float pulse_out0,pulse_out1,pulse_out2,pulse_out3;

	pulse_out0 = thr + att_out[Pitch] - att_out[Roll] - att_out[Yaw];
	pulse_out1 = thr - att_out[Pitch] - att_out[Roll] + att_out[Yaw];
	pulse_out2 = thr + att_out[Pitch] + att_out[Roll] + att_out[Yaw];
	pulse_out3 = thr - att_out[Pitch] + att_out[Roll] - att_out[Yaw];

	pulse_out0 = constrain_f(pulse_out0, RC_min, RC_max);
	pulse_out1 = constrain_f(pulse_out1, RC_min, RC_max);
	pulse_out2 = constrain_f(pulse_out2, RC_min, RC_max);
	pulse_out3 = constrain_f(pulse_out3, RC_min, RC_max);

	set_pulse_ch0(pulse_out0);
	set_pulse_ch1(pulse_out1);
	set_pulse_ch2(pulse_out2);
	set_pulse_ch3(pulse_out3);
#endif

#ifdef Hexa
  static float pulse_out0,pulse_out1,pulse_out2,pulse_out3,pulse_out4,pulse_out5;
  
  raw_motors[0] = throttle + att_out[Roll] + att_out[Pitch] - att_out[Yaw];
  raw_motors[1] = throttle + att_out[Roll] + 0.0f + yaw;
  raw_motors[2] = throttle + att_out[Roll] - att_out[Pitch] - att_out[Yaw];
  raw_motors[3] = throttle - att_out[Roll] + att_out[Pitch] + att_out[Yaw];
  raw_motors[4] = throttle - att_out[Roll] + 0.0f - yaw;
  raw_motors[5] = throttle - att_out[Roll] - att_out[Pitch] + att_out[Yaw];
  
  pulse_out0 = constrain_f(pulse_out0, RC_min, RC_max);
  pulse_out1 = constrain_f(pulse_out1, RC_min, RC_max);
  pulse_out2 = constrain_f(pulse_out2, RC_min, RC_max);
  pulse_out3 = constrain_f(pulse_out3, RC_min, RC_max);
  pulse_out4 = constrain_f(pulse_out4, RC_min, RC_max);
  pulse_out5 = constrain_f(pulse_out5, RC_min, RC_max);

	set_pulse_ch0(pulse_out0);
	set_pulse_ch1(pulse_out1);
	set_pulse_ch2(pulse_out2);
	set_pulse_ch3(pulse_out3);
  set_pulse_ch4(pulse_out4);
	set_pulse_ch5(pulse_out5);

#endif

}

__NO_RETURN void task_control(void *argument)
{
	(void)argument;
	static uint64_t us;
	static uint64_t us_stable_pid;
	static float thr_alt;

	while (!system_state.state_init)
	{
		update_imu_raw();
		iwdg_refresh();
		Update_Controls(RC_min);
		osDelay(10);
	}

	while (1)
	{
		const uint64_t now = get_us64();
		float dt_ = (float)(now - us) * 1e-6f;
		us = now;
		dt_control += 0.999f * (dt_ - dt_control);

		update_imu_raw();

		switch (system_state.Armed_mode)
		{
		case ARM_Auto:
		case ARM_Normal:
		case ARM_GCS:
		{
			// run stable pid @100hz
			uint64_t us_ms = now - us_stable_pid;
			if (us_ms > 9000 && us_ms < 10000)
			{
				us_stable_pid = now;

				const float dt_spid = (float)us_ms * 1.0e-6f;
				dt_stable_pid = 0.998f * dt_stable_pid + dt_spid * 0.002f;

				run_angle_controller(dt_stable_pid);
				run_thr_controller(control_data.thr_base, dt_stable_pid);
        thr_alt = control_data.final_thr_out;
				
        if (thr_alt < idle_thr)
				{
					thr_alt = idle_thr;
				}
				strobe_led('a');
			}
			run_rate_controller(dt_control);
			/*Final Thr Out*/
			Update_Controls(thr_alt);
		}
		break;
		default:
			if (!m_calState.cal_esc)
				Update_Controls(RC_min);
			us_stable_pid = now;
			break;
		}
		//	iwdg_refresh();
		osDelay(1);
	}
}

__NO_RETURN void task_navprocess(void *argument)
{
	(void)argument;
	static uint64_t us;
	float alt_hold_thr = RC_min;
	float thr_level = 0.0f;
	float diff_latN, diff_lonE;
	static float est_pos_n, est_pos_e, vel_n, vel_e, est_pos_u, est_yaw;
	static float current_throttle;
	static float circle_bearing;
	float dt;
	float alt_err = 0.0f;
	float max_geofence_alt;

	// sould be for running position controller, outer loop
	// to genarate velocity set points
	while (!system_state.state_init)
	{
		osDelay(2);
	}

	while (1)
	{
		const uint64_t now = get_us64();
		dt = (float)(now - us) * 1e-6f;
		us = now;
		dt_process += 0.999f * (dt - dt_process);

		// low gain lpf to reduce jitter noise
		// using altitde in body frame to support configurations
		// as there no diffrence apart from sign between NED and Body
		est_pos_n = get_est_posN();	 // in m
		est_pos_e = get_est_posE();	 // in m
		est_pos_u = get_est_alt_m(); // in m
		vel_n = get_est_velN();
		vel_e = get_est_velE();
		est_yaw = get_est_yaw();
		current_throttle = control_data.set_point_pwm[THR];
		max_geofence_alt = (float)ctrl_config.max_geofence_alt;
		const float rate_vertical_spd = ctrl_config.rate_Vspd;
		const float rate_horizontal_spd = ctrl_config.rate_Hspd;

		switch (system_state.Ctl_Filt_mode)
		{
		case Stabilize:
		{
			set_pos_lock(0);
			set_alt_lock(0);
			set_course_lock(0);
			set_circle_lock(0);
      set_zupt_pos_hold(0);
			alt_err = 0.0f;
			thr_level = 0.0f;
			control_data.integral_posN = 0.0f;
			control_data.integral_posE = 0.0f;
			control_data.integral_posD = 0.0f;
			control_data.pid_out_velN = 0.0f;
			control_data.pid_out_velE = 0.0f;
			control_data.pid_out_velD = 0.0f;
			velocityN_sp = 0.0f;
			velocityE_sp = 0.0f;
			velocityD_sp = 0.0f;
			control_data.thr_base = current_throttle; // if none of flight mode set
		}
		break;
		case Heading_hold:
		{
			set_course_lock(1);
		}
		break;
		case Auto_takeoff:
		{
			char alt_l = get_alt_lock();
			// only if system not airboarn
			if (!alt_l && !system_state.Airborne)
			{
				// Take avg before setting altitude
				// if target altitude set befor triggering althold.
				// set altitude should be more than +-Min_ground_alt mtr
				if (set_target_active_alt(est_pos_u) > 0)
				{
          if(current_throttle<=RC_min)
            alt_hold_thr = current_throttle;
				}
				if (fabsf(control_data.target_alt) < Min_takeoff_alt)
				{
					control_data.target_alt = Min_takeoff_alt;
				}
			}
			else if (alt_l)
			{
				// make the system airboarn untill system have minimum ground alt
				// increase the throttle
				detect_airborne(est_pos_u);
				if (!system_state.Airborne)
				{
					alt_hold_thr += (25.0f * dt_process); // per second
					alt_hold_thr = (alt_hold_thr > RC_max) ? RC_max : alt_hold_thr;
				}
				else
				{
					// put to altitude hold
					system_state.Ctl_Filt_mode = Altitude_hold;
				}
			}
		}
		break;
		case Altitude_hold:
		{
			// for first time altitude hold trigger
			if (!get_alt_lock())
			{
				// Take avg before setting altitude
				// if target altitude set befor triggering alt hold.
				// set altitude should be more than +-Min_ground_alt mtr
				if (set_target_active_alt(est_pos_u) > 0)
				{
          velocityD_sp = 0;
					alt_hold_thr = current_throttle;
				}
			}
			if (get_alt_lock())
			{
				// 1 / (RC_max - RC_min)
				float thr_sp = (current_throttle - alt_hold_thr) * 9.090909e-4f;

				if (fabsf(thr_sp) > THR_ALT_DB)
				{ // NEU frame
					// climbing (Up positive)
					// descending (Down negetive)
					thr_level = thr_sp * rate_vertical_spd;
					control_data.target_alt = est_pos_u;
				}
				else
				{
					thr_sp = 0.0f;
					thr_level = 0.0f;
				}

				alt_err = control_data.target_alt - est_pos_u;

				run_alt_controller(alt_err, dt_process);
				velocityD_sp = control_data.pid_out_velD + thr_level;

				control_data.thr_base = alt_hold_thr;
			}
		}
		break;
		case Position_hold:
		{
			if (fabsf(control_data.set_point_pwm[Roll]) > 2.0f || fabsf(control_data.set_point_pwm[Pitch]) > 2.0f)
			{
				set_pos_lock(0);
			}
			else if (!get_pos_lock())
			{
				if (sqrtf(vel_n * vel_n + vel_e * vel_e) < 0.2f)
				{
          velocityN_sp = 0;
          velocityE_sp = 0;
					set_target_active_pos(est_pos_n, est_pos_e);
				}
			}
			if (get_pos_lock())
			{
				diff_latN = control_data.target_posN - est_pos_n;
				diff_lonE = control_data.target_posE - est_pos_e;

				/* distance to hold point (m) */
				float d = sqrtf(diff_latN * diff_latN + diff_lonE * diff_lonE);

				if (d > pos_guid_range && d > 1e-6f)
				{
					control_data.target_posN = est_pos_n + (diff_latN / d) * pos_guid_range;
					control_data.target_posE = est_pos_e + (diff_lonE / d) * pos_guid_range;
					diff_latN = control_data.target_posN - est_pos_n; // Recalculate
					diff_lonE = control_data.target_posE - est_pos_e;
				}
				run_pos_controller(diff_latN, diff_lonE, dt_process);

				velocityN_sp = control_data.pid_out_velN;
				velocityE_sp = control_data.pid_out_velE;
				float magv = sqrtf(velocityN_sp * velocityN_sp + velocityE_sp * velocityE_sp);
				if (magv > rate_horizontal_spd)
				{
					float s = rate_horizontal_spd / magv;
					velocityN_sp *= s;
					velocityE_sp *= s;
				}
				if (!get_alt_lock())
				{
					control_data.thr_base = current_throttle;
				}
			}
		}
		break;
		case Loiter_mode:
		{
			// alt+position
			if (!get_alt_lock())
			{
				if (set_target_active_alt(est_pos_u) > 0)
				{
					alt_hold_thr = current_throttle;
          velocityD_sp = 0;
				}
			}
			if (fabsf(control_data.set_point_pwm[Roll]) > 2.0f || fabsf(control_data.set_point_pwm[Pitch]) > 2.0f)
			{
				set_pos_lock(0);
			}
			else if (!get_pos_lock())
			{
				if (sqrtf(vel_n * vel_n + vel_e * vel_e) < 0.2f)
				{
          velocityN_sp = 0;
          velocityE_sp = 0;
					set_target_active_pos(est_pos_n, est_pos_e);
				}
			}
			if (get_alt_lock())
			{
				float thr_sp = (current_throttle - alt_hold_thr) * 9.090909e-4f;

				if (fabsf(thr_sp) > THR_ALT_DB)
				{ // NEU frame
					thr_level = thr_sp * rate_vertical_spd;
					control_data.target_alt = est_pos_u;
				}
				else
				{
					thr_sp = 0.0f;
					thr_level = 0.0f;
				}

				alt_err = control_data.target_alt - est_pos_u;
				run_alt_controller(alt_err, dt_process);
				velocityD_sp = control_data.pid_out_velD + thr_level;
				// geting control_data.pid_out_pos[2] from vel D controller
				control_data.thr_base = alt_hold_thr;
			}
			if (get_pos_lock())
			{
				diff_latN = control_data.target_posN - est_pos_n;
				diff_lonE = control_data.target_posE - est_pos_e;

				/* distance to hold point (m) */
				float d = sqrtf(diff_latN * diff_latN + diff_lonE * diff_lonE);

				if (d > pos_guid_range && d > 1e-6f)
				{
					control_data.target_posN = est_pos_n + (diff_latN / d) * pos_guid_range;
					control_data.target_posE = est_pos_e + (diff_lonE / d) * pos_guid_range;
					diff_latN = control_data.target_posN - est_pos_n; // Recalculate
					diff_lonE = control_data.target_posE - est_pos_e;
				}
				run_pos_controller(diff_latN, diff_lonE, dt_process);

				velocityN_sp = control_data.pid_out_velN;
				velocityE_sp = control_data.pid_out_velE;
				float magv = sqrtf(velocityN_sp * velocityN_sp + velocityE_sp * velocityE_sp);
				if (magv > rate_horizontal_spd)
				{
					float s = rate_horizontal_spd / magv;
					velocityN_sp *= s;
					velocityE_sp *= s;
				}
				if (!get_alt_lock())
				{
					control_data.thr_base = current_throttle;
				}
			}
		}
		break;
		case RTL_mode:
		{
			float dist_2origin;
			double origin_ned[3];
			double lla_current[3] = {0, 0, 0};
			// convert origin pos to ned
			lla_current[0] = get_est_latD();
			lla_current[1] = get_est_lonD();

			lla2ned(struct_home_data.lla_origin, lla_current, origin_ned);

			diff_latN = origin_ned[0] - est_pos_n;
			diff_lonE = origin_ned[1] - est_pos_e;
			// distance to origin
			dist_2origin = dx_dy2_length(origin_ned[0], est_pos_n, origin_ned[1], est_pos_e);
			// heading to origin
			control_data.target_course = dx_dy2_bearing(origin_ned[0], est_pos_n, origin_ned[1], est_pos_e);

			if (dist_2origin > 5.0f)
			{
				control_data.target_alt = RTL_Altitude;
				control_data.target_Hvelocity = RTL_Velocity;
			}
			else
			{
				if (control_data.target_Hvelocity > 0.0f)
					control_data.target_Hvelocity -= 0.2f;

				if (fabsf(dist_2origin) < 1.0f || system_state.Airborne)
				{
					if (control_data.target_alt > 0.0f)
						control_data.target_alt -= 1.0f;
				}
			}
			float alt_dff = control_data.target_alt - est_pos_u;
			run_alt_controller(alt_dff, dt_process);
			run_pos_controller(diff_latN, diff_lonE, dt_process);
			detect_landing();
		}
		break;
		case Circle_mode:
		{
			float edg_ned[2];
			if (!get_circle_lock())
			{ // init the circle mode if lock is not set
				circle_bearing = est_yaw;
				// create circle edge coordinate with same heading and cicle dia
				if ((struct_circle.flags & 0x1) == 0) // origin not set
				{									  // ned
					struct_circle.circle_origin[0] = est_pos_n;
					struct_circle.circle_origin[1] = est_pos_e;
					struct_circle.flags = 1;
				}

				set_circle_lock(1);
			}
			if (get_circle_lock())
			{
				// get the edge pos in meters
				bearing_length2_xy(struct_circle.circle_origin[0], struct_circle.circle_origin[1], &edg_ned[0], &edg_ned[1], circle_bearing, struct_circle.circle_dia);

				float dist_2edge = dx_dy2_length(struct_circle.circle_origin[0], struct_circle.circle_origin[1], edg_ned[0], edg_ned[1]);
				if (dist_2edge < 1.5f)
				{
					circle_bearing++;
					if (circle_bearing > 360.0f)
						circle_bearing = 0.0f;

					float cir_course = circle_bearing + 90.0f;

					if (cir_course > 360.0f)
						control_data.target_course = cir_course - 360.0f;
					else
						control_data.target_course = cir_course;
				}

				diff_latN = edg_ned[0] - est_pos_n;
				diff_lonE = edg_ned[1] - est_pos_e;

				if (control_data.target_alt != struct_circle.circle_alt)
					control_data.target_alt = struct_circle.circle_alt;

				control_data.target_Hvelocity = struct_circle.circle_speed;
				float alt_diff = control_data.target_alt - est_pos_u;
				run_alt_controller(alt_diff, dt_process);
				run_pos_controller(diff_latN, diff_lonE, dt_process);
			}
		}
		break;
		case Mission_mode:
			run_wp_controller(est_pos_n, est_pos_e, est_yaw, est_pos_u, dt_process);
			break;
		}

		if (est_pos_u > max_geofence_alt)
			control_data.thr_base -= 0.5f;

		osDelay(18);
	}
}

__NO_RETURN void task_navcontrol(void *argument)
{
	(void)argument;
	static uint64_t us;
	static float set_point_yaw;
	static float set_point_body[3];
	static float velocity_sp_n, velocity_sp_e, velocity_sp_d;

	// sould be for running velocity controller, inner loop
	while (!system_state.state_init)
	{
		osDelay(10);
	}

	while (1)
	{
		uint64_t now = get_us64();
		float dt = (float)(now - us) * 1e-6f;
		us = now;

		dt_navctrl += 0.999f * (dt - dt_navctrl);

		calculate_alt_data();
		float est_yaw = get_est_yaw();
		float yaw_sp = control_data.set_point_pwm[Yaw];

		// control_data.pid_out_vel[N] , coming from position controller
		// control_data.target_vel_ned[N] , coming from ext control
		// first create velocity set point

		velocity_sp_n = 0;
		velocity_sp_e = 0;
		velocity_sp_d = 0;

		// horizontal position part
		velocity_sp_n = velocityN_sp + control_data.target_vel_n + control_data.target_velN_cmd;
		velocity_sp_e = velocityE_sp + control_data.target_vel_e + control_data.target_velE_cmd;
		if (get_pos_lock())
		{
			run_velNE_controller(velocity_sp_n, velocity_sp_e, dt_navctrl);  
		}
		else
		{
			run_velNE_controller(velocity_sp_n, velocity_sp_e, dt_navctrl);
//      control_data.pid_out_pos[N] = 0.0f;
//      control_data.pid_out_pos[E] = 0.0f;
		}
		// vertical position part
		velocity_sp_d = velocityD_sp + control_data.target_vel_d;
		if (get_alt_lock())
		{
			run_velD_controller(velocity_sp_d, dt_navctrl);
		}
		else
			control_data.pid_out_posD = 0.0f;

		// control_data.pid_out_pos , coming from velocity controller
		float accel_sp_ned[3] = {0, 0, 0};
		// translating NED frame to bofy frame
		accel_sp_ned[0] = control_data.pid_out_posN; // from velocity controller
		accel_sp_ned[1] = control_data.pid_out_posE; // from velocity controller
		
		accel_2lean_angle(accel_sp_ned, est_yaw, set_point_body);

//		if (fabsf(control_data.set_point_pwm[Roll]) > 0.01f && system_state.RC_Link)
//			control_data.set_point_att[Roll] = control_data.set_point_pwm[Roll];
//		else
			control_data.set_point_att[Roll] = set_point_body[Roll];

//		if (fabsf(control_data.set_point_pwm[Pitch]) > 0.01f && system_state.RC_Link)
//			control_data.set_point_att[Pitch] = control_data.set_point_pwm[Pitch];
//		else
			control_data.set_point_att[Pitch] = set_point_body[Pitch];

		if (!get_course_lock())
		{
			//*Yaw control, genarating yaw setpoint
			/*In this modes course will control by self*/
			set_point_yaw = 0.0f;
			if (yaw_sp > 5.0f)
			{
				control_data.target_course = set_point_body[Yaw];
				set_point_yaw = yaw_sp - 2.5f;
			}
			if (yaw_sp < -5.0f)
			{
				control_data.target_course = set_point_body[Yaw];
				set_point_yaw = yaw_sp + 2.5f;
			}
		}
		else
			set_point_yaw = 0.0f;

		control_data.set_point_att[Yaw] = control_data.target_course + set_point_yaw;
	
		osDelay(5);
	}
}

void update_control_setpoint(void)
{
	float delta[6] = {0.0f, 0.0f, 0.0f, RC_min, 0.0f, 0.0f};
	static float wait_t;

	if (Init_Cfg.control_bits & 0x1)
	{
		delta[Roll] = ext_control.primary[Roll] * 0.01f;
		delta[Pitch] = ext_control.primary[Pitch] * 0.01f;
		delta[Yaw] = ext_control.primary[Yaw] * 0.01f;
		delta[THR] = map(ext_control.primary[THR], struct_data_e2p.cal_pwm_min[THR], struct_data_e2p.cal_pwm_max[THR], RC_min, RC_max);
		delta[AUX1] = ext_control.primary[AUX1] * 0.01f;
		delta[AUX2] = ext_control.primary[AUX2] * 0.01f;
	}
	else
	{ // if control inputs are coming from rc transmitter
		if (m_calState.cal_rc == 0)
		{
			short pwm_data[6];
			get_pwm_data(pwm_data);
			delta[Roll] = pwm_2deg(AIL, pwm_data);
			delta[Pitch] = pwm_2deg(ELE, pwm_data);
			delta[Yaw] = pwm_2deg(RUD, pwm_data);
			delta[THR] = map(pwm_data[THR], struct_data_e2p.cal_pwm_min[THR], struct_data_e2p.cal_pwm_max[THR], RC_min, RC_max);
			delta[AUX1] = pwm_data[AUX1];
			delta[AUX2] = pwm_data[AUX2];
		}
	}
	if (system_state.Armed_mode > DIS_ARM)
	{
		// yes sytem is armed
		// failsafe
		detect_airborne(get_est_alt_m());
		if (system_state.Airborne)
		{
			if (!system_state.RC_Link)
			{
				if (wait_t > airborne_wait_ms)
				{
					wait_t = 0;
					print_Debug("No Signal from Transmitter\n");
					// set_ctl_flight_mode(ctrl_config.breach_flag);
					return;
				}
				wait_t++;
			}
		}
	}

	delta[Roll] = fabsf(delta[Roll]) > pwm_deadband ? delta[Roll] : 0.0f;
	delta[Pitch] = fabsf(delta[Pitch]) > pwm_deadband ? delta[Pitch] : 0.0f;
	delta[Yaw] = fabsf(delta[Yaw]) > pwm_deadband ? delta[Yaw] : 0.0f;
  delta[THR] = (delta[THR]<RC_min) ? RC_min : delta[THR];
	delta[AUX1] = fabsf(delta[AUX1]) > pwm_deadband ? delta[AUX1] : 0.0f;
	delta[AUX2] = fabsf(delta[AUX2]) > pwm_deadband ? delta[AUX2] : 0.0f;

	control_data.set_point_pwm[Roll] += 0.05f * (delta[Roll] - control_data.set_point_pwm[Roll]);
	control_data.set_point_pwm[Pitch] += 0.05f * (delta[Pitch] - control_data.set_point_pwm[Pitch]);
	control_data.set_point_pwm[Yaw] += 0.05f * (delta[Yaw] - control_data.set_point_pwm[Yaw]);
	control_data.set_point_pwm[THR] += 0.05f * (delta[THR] - control_data.set_point_pwm[THR]);
	control_data.set_point_pwm[AUX1] = delta[AUX1];
	control_data.set_point_pwm[AUX2] = delta[AUX2];

	// convert stick input to NED frame
	// if control bypass from angle controller to velocity control
	update_velocity_setpoints(control_data.set_point_pwm[Pitch],control_data.set_point_pwm[Roll]);
}

float pwm_2deg(int control, const short *pwm_in)
{
	float delta = 0.0f;

	if (pwm_in[control] < (RC_min - 100))
	{
		system_state.RC_Link = 0;
		return 0;
	}
	else
		system_state.RC_Link = 1;

	if (pwm_in[control] > cal_RCusec_mid[control])
		delta = (pwm_in[control] - cal_RCusec_mid[control]) * delta_f;
	else if (pwm_in[control] < cal_RCusec_mid[control])
		delta = (pwm_in[control] - cal_RCusec_mid[control]) * delta_f;
	else
	{
		delta = 0.0f;
	}
	return delta * Scale_factor_att;
}

void set_sysARM_mode(void *mode)
{
	unsigned char arm_mode = DIS_ARM;
	if (mode != NULL)
	{
		arm_mode = (*(unsigned char *)mode);
	}

	switch (arm_mode)
	{
	case ARM_Normal:
	case ARM_Auto:
	case ARM_GCS:
		if (system_state.ekf_health && system_state.mag_health && system_state.accl_health && system_state.gyro_health)
		{
			system_state.Armed_mode = arm_mode & 0b11;
			control_data.target_course = get_est_yaw();
			control_data.set_point_att[Yaw] = get_est_yaw();
			control_data.ground_alt = get_est_alt_m();
      thr_controller_init();
			print_Debug("System ARMED\n");
      arming_init = 0;
		}
		else
		{
			print_Debug("System Bad Sensor Health, ARMing Skiped\n");
		}
		break;
	default:
		system_state.Armed_mode = DIS_ARM;
		break;
	}
}

void set_preARM_data(const unsigned char *arm_mode)
{
	if (*arm_mode == DIS_ARM)
	{
		system_state.Armed_mode = DIS_ARM;
		gpioG_rstPin(status_led_green);
		// print_Debug("System DisARMED\n");
	}
	else
	{
		if (m_calState.cal_rc == 0)
		{
			if (!(system_state.Armed_mode > DIS_ARM))
			{
        if(!arming_init)
        {
          arming_init = 1;
          initiate_caibration(Cal_Gps, &struct_home_data, NULL);
          initiate_caibration(Cal_Baro, &struct_home_data, (void *)arm_mode);
        }
			}
		}
	}

	system_state.Airborne = 0;
	att_controller_reset();
	memset(lpf_att.momentum, 0, 12);
	memset(lpf_att.lastInput, 0, 12);
	control_data.target_posN = 0.0f;
	control_data.target_posE = 0.0f;
	control_data.target_course = get_est_yaw();
	control_data.set_point_att[Yaw] = get_est_yaw();
}

char detect_landing(void)
{
	static uint32_t wait_t;
	static float accel_norm_change;

	if (wait_t == 0)
		wait_t = HAL_GetTick();

	float baro_alt_err = struct_home_data.lla_origin[2] - get_est_alt_m();
	float accelN = get_est_acclN();
	float accelE = get_est_acclE();
	float accelD = get_est_acclD();

	float accel_norm = sqrtf(sqr(accelN) + sqr(accelE) + sqr(accelD));
	float change = (accel_norm - accel_norm_change) * dt_process;
	accel_norm_change = accel_norm;

	if (baro_alt_err <= 1.0f && control_data.thr_base <= (RC_min + 100) && fabsf(change) < 0.5f)
	{
		if ((HAL_GetTick() - wait_t) > Landing_wait_sec)
		{
			wait_t = HAL_GetTick();
			system_state.Airborne = 0;
			control_data.thr_base = RC_min;
			return 1;
		}
	}
	return 0;
}

void detect_airborne(float est_alt_m)
{
	static unsigned int wait_t;

	if (wait_t == 0)
		wait_t = HAL_GetTick();

	if (!system_state.Airborne)
	{
		// validate the ground diffrence
		if (fabsf(est_alt_m) - fabsf(control_data.ground_alt) >= Min_ground_alt)
		{
			// wait 5 sec
			if ((HAL_GetTick() - wait_t) > airborne_wait_ms)
			{
				control_data.target_alt = est_alt_m;
				wait_t = 0;
				system_state.Airborne = 1;
			}
		}
	}
}

void set_ctl_flight_mode(char ctl_mode)
{
	if (ctl_mode != system_state.Ctl_Filt_mode)
	{
		set_pos_lock(0);
		set_alt_lock(0);
		set_course_lock(0);
		set_circle_lock(0);
		control_data.integral_posN = 0.0f;
		control_data.integral_posE = 0.0f;
		control_data.integral_posD = 0.0f;
		control_data.pid_out_velN = 0.0f;
		control_data.pid_out_velE = 0.0f;
		control_data.pid_out_velD = 0.0f;
		velocityN_sp = 0.0f;
		velocityE_sp = 0.0f;
		velocityD_sp = 0.0f;
		system_state.Ctl_Filt_mode = ctl_mode & 0b1111;
	}
}

void set_nav_flight_mode(const unsigned char ctl_mode, const unsigned char nav_mode)
{
	//	if(ext_control.pwm_in[AUX1]<struct_data_e2p.flight_mode_us[index]+50 && ext_control.pwm_in[AUX1]>struct_data_e2p.flight_mode_us[index]-50)
	//	{
	//		mode = struct_data_e2p.flight_mode_index[index];
	//	}

	if (nav_mode != system_state.Nav_Filt_mode)
	{
		switch (nav_mode)
		{
		case Mission_mode:
		case Circle_mode:
			system_state.Nav_Filt_mode = nav_mode & 0b11;
			break;
		}
	}
	set_ctl_flight_mode(ctl_mode);
}

void set_target_course(float course)
{
	control_data.target_course = course;
}

void strobe_led(char pattern)
{
	// should be called at task frequency rate
	static uint32_t on_ms = 0;	// = HAL_GetTick();
	static uint32_t off_ms = 0; // = HAL_GetTick();
	static int8_t flag_tg = 0;

	switch (pattern)
	{
	case 'a':
		if ((HAL_GetTick() - on_ms > 75) && flag_tg) // 200ms
		{
			GPIOG->ODR ^= status_led_green;

			on_ms = HAL_GetTick();
		}
		if (HAL_GetTick() - off_ms > 500)
		{
			flag_tg = !flag_tg;
			GPIOG->BSRR |= status_led_green << 16U;
			off_ms = HAL_GetTick();
		}
		break;
	case 'b':

		break;
	case 'c':
		if (HAL_GetTick() - on_ms > 200) // 200ms
		{
			gpioG_togPin(status_led_green | status_led_blue | status_led_green);
			on_ms = HAL_GetTick();
		}
		break;
	}
}

void fcs_init(void)
{
	dt_control = 0.001f;
	dt_process = 1.0f / 50.0f;
	dt_navctrl = 1.0f / 200.0f;
	dt_stable_pid = 0.004f;

	lpf_att.c1 = 0.185f;
	lpf_att.c2 = 1.25e-10f;

	lpf_alt.c1 = 0.350f;
	lpf_alt.c2 = 1.25e-7f;

	lpf_pos[N].c1 = 0.35f;
	lpf_pos[N].c2 = 1.0e-7f;

	lpf_pos[E].c1 = 0.35f;
	lpf_pos[E].c2 = 1.0e-7f;

	delta_f = 15.0f/(float)ctrl_config.s_div;

	idle_thr = RC_min + ((float)RC_min * ((float)ctrl_config.idle_thr * 0.01f));

	cal_RCusec_mid[AIL] = (struct_data_e2p.cal_pwm_max[AIL] + struct_data_e2p.cal_pwm_min[AIL]) * 0.5f;
	cal_RCusec_mid[ELE] = (struct_data_e2p.cal_pwm_max[ELE] + struct_data_e2p.cal_pwm_min[ELE]) * 0.5f;
	cal_RCusec_mid[RUD] = (struct_data_e2p.cal_pwm_max[RUD] + struct_data_e2p.cal_pwm_min[RUD]) * 0.5f;
	cal_RCusec_mid[THR] = (struct_data_e2p.cal_pwm_max[THR] + struct_data_e2p.cal_pwm_min[THR]) * 0.5f;
	cal_RCusec_mid[AUX1] = (struct_data_e2p.cal_pwm_min[AUX1] + struct_data_e2p.cal_pwm_max[AUX1]) * 0.5f;
	cal_RCusec_mid[AUX2] = (struct_data_e2p.cal_pwm_min[AUX2] + struct_data_e2p.cal_pwm_max[AUX2]) * 0.5f;

	control_data.set_point_pwm[Roll] = 0.0f;
	control_data.set_point_pwm[Pitch] = 0.0f;
	control_data.set_point_pwm[Yaw] = 0.0f;
	control_data.set_point_pwm[THR] = RC_min;
	control_data.set_point_pwm[AUX1] = 0.0f;
	control_data.set_point_pwm[AUX2] = 0.0f;

	att_controller_init();
	vel_controller_init();
	pos_controller_init();
  thr_controller_init();
	print_Debug("FCS Init Done\n");
}
