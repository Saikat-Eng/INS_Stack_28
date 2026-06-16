#include <drv_utils.h>
#include <math.h>
#include <thr_control.h>
#include <flight_control.h>
#include <estimator.h>
#include <pos_control.h>
#include <util.h>
#include <nav_sys_config.h>
#include <calibration.h>

// Relaxed bounds for 100Hz EKF-fused states on high-vibration frames
#define HOVER_EST_ACCEL_THRESH  1.2f  
#define HOVER_EST_VEL_THRESH    0.4f 

// Time constant in seconds (Independent of frequency, calibrated for 100Hz execution)
#define HOVER_LEARN_TIME_CONST  3.0f 

extern float dbg_var1,dbg_var2,dbg_var3,dbg_var4;
static float thr_centre,thr_Kp;

static float estimated_hover_ratio; 
static float accel_2thrust_cmd(float accel_des, float cos_tilt);
static void hover_thrust_estimator(float velU_est, float accU_est, float current_pwm, float dt);
static float get_clamped_cos_tilt(void);

void thr_controller_init()
{
  thr_Kp = 0.8f;
  estimated_hover_ratio = 0.5f;  
  thr_centre = (RC_max+RC_min)*0.5f;
}

void run_thr_controller(float thr_base, float dt)
{
	float thr_out_pwm = RC_min;
	float accel_z = -get_est_acclD();
	float vel_z = get_est_velU();
  
	float cos_tilt = get_clamped_cos_tilt();

	if (get_alt_lock())
	{  
     float target_accel = control_data.pid_out_posD;
     dbg_var4 = target_accel;
     thr_out_pwm = accel_2thrust_cmd(target_accel, cos_tilt);
     hover_thrust_estimator(vel_z, accel_z, thr_out_pwm, dt);
	}
	else
	{
     // Manual Throttle Mode
     float active_thr_request = thr_base - RC_min;
     if (active_thr_request < 0.0f) active_thr_request = 0.0f;
	 
     float compensated_thr = active_thr_request / cos_tilt;

     // Active velocity-dampening feedback loop for manual flight stability
     float momentum_damping = accel_z * thr_Kp;
     thr_out_pwm = RC_min + compensated_thr - momentum_damping;
	}
  dbg_var1 = thr_out_pwm;
  control_data.final_thr_out = thr_out_pwm;
}

void hover_thrust_estimator(float velU_est, float accU_est, float current_pwm, float dt) 
{
    if (fabsf(velU_est) < HOVER_EST_VEL_THRESH && fabsf(accU_est) < HOVER_EST_ACCEL_THRESH) {
        
        float pwm_span = thr_centre - RC_min;
        if (pwm_span <= 0.0f) return;

        // Normal PWM ESC Inverse Mapping (quadratic relationship verified)
        float normalized_pwm = (current_pwm - RC_min) / pwm_span;
        float actual_ratio = normalized_pwm * normalized_pwm; 

        actual_ratio = constrain_f(actual_ratio, 0.25f, 0.85f);

        // Mathematical Time-Constant conversion at 100Hz execution 
        float alpha = dt / (HOVER_LEARN_TIME_CONST + dt);
        
        estimated_hover_ratio = ((1.0f - alpha) * estimated_hover_ratio) + (alpha * actual_ratio);
        dbg_var2 = estimated_hover_ratio;
    }
}

float accel_2thrust_cmd(float accel_des, float cos_tilt)
{        
    float total_accel = accel_des + Actual_gravity;
    float thrust_ratio = (estimated_hover_ratio * total_accel * Gravity_inv) / cos_tilt;
    
    thrust_ratio = constrain_f(thrust_ratio, 0.15f, 0.95f);
    float pwm_span = thr_centre - RC_min;
    
    // Output curve mapping for normal PWM ESCs
    float desired_pwm = RC_min + (pwm_span * sqrtf(thrust_ratio));
    dbg_var3 = desired_pwm;
    return desired_pwm;
}

float get_clamped_cos_tilt(void) 
{
    const float roll_rad = radians_f(get_est_roll());
    const float pitch_rad = radians_f(get_est_pitch());
    
    float cos_tilt = cosf(pitch_rad) * cosf(roll_rad);
    // Clamp to 45 degrees max lean compensation (0.707) to preserve roll/pitch authority
    return constrain_f(cos_tilt, 0.707f, 1.0f);
}