/**
 * @file imu.c
 * @brief Implementation of the IMU data processing and calibration routines.
 *
 * This file contains the implementation of the IMU data processing and calibration routines.
 *
 * @author saikat
 * @date 2023-12-01
 * @version 1.0
 */
#include <math.h>
#include <spi.h>
#include <gpio.h>
#include <e2prom.h>
#include <imu.h>
#include <imu_defines.h>
#include <nav_sys_config.h>
#include <nav.h>
#include <system.h>
#include <combiner.h>
#include <estimator.h>
#include <util.h>
#include <drv_utils.h>
#include <dma.h>
#include "cmsis_os2.h"
#include <task.h>
#include <matrix.h>
#include <stdio.h>
#include <imu_calibration.h>
//#include <INSEstimator.h>
#include <gyro_fft_peak.h>

static volatile imu_data struct_imu_raw __attribute__((section(".DTCM_RAM")));
static volatile air_data struct_air_raw __attribute__((section(".DTCM_RAM")));
static volatile mag_data struct_mag_raw __attribute__((section(".DTCM_RAM")));
static volatile imu_ring_buffer m_accel_ring __attribute__((section(".DTCM_RAM")));
static volatile imu_ring_buffer m_gyro_ring __attribute__((section(".DTCM_RAM")));
static volatile imu_ring_buffer m_mag_ring __attribute__((section(".DTCM_RAM")));
static volatile int16_t A_16x,A_16y,A_16z,G_16x,G_16y,G_16z;
static volatile uint8_t mpu_flag = 0;
static volatile uint8_t buf_mpu9250[15];
static void accumulate_imu_data(void);
static void imu_exti_callback(void);
static m_avg lpf_press;
static hpf_1D lpf_alt, lpf_alt_rate;
static hpf_3D lpf_gyro;
static float dt_gfft;

static osMutexId_t imu_mutex;
static osMutexId_t mag_mutex;
static osMutexId_t baro_mutex;

float sim_alt_rate = 0.0f;

void imu_exti_callback(void)
{
	static const uint8_t push[15] = {ACCEL_OUT_Read};

	GPIOE->ODR &= ~spi4_cs;
	//	accumulate_imu_data();
	spi4_read_dma(buf_mpu9250, 15);
	spi4_write_dma(push, 15);
}

void accumulate_imu_data(void)
{
	mpu_flag = 0;
	int16_t imu_data_0,imu_data_1,imu_data_2,imu_data_3,imu_data_4,imu_data_5;
	// updating from interrupt
	imu_data_0 = ((int16_t)buf_mpu9250[1] << 8 | (int16_t)buf_mpu9250[2]);
	imu_data_1 = ((int16_t)buf_mpu9250[3] << 8 | (int16_t)buf_mpu9250[4]);
	imu_data_2 = ((int16_t)buf_mpu9250[5] << 8 | (int16_t)buf_mpu9250[6]);
	imu_data_3 = ((int16_t)buf_mpu9250[9] << 8 | (int16_t)buf_mpu9250[10]);
	imu_data_4 = ((int16_t)buf_mpu9250[11] << 8 | (int16_t)buf_mpu9250[12]);
	imu_data_5 = -((int16_t)buf_mpu9250[13] << 8 | (int16_t)buf_mpu9250[14]);

  GPIOE->ODR |= spi4_cs;
	// added 1st order lpf
	A_16x = (A_16x * 7 + imu_data_0) >> 3;
	A_16y = (A_16y * 7 + imu_data_1) >> 3;
	A_16z = (A_16z * 7 + imu_data_2) >> 3;
                               
	G_16x = (G_16x * 7 + imu_data_3) >> 3;
	G_16y = (G_16y * 7 + imu_data_4) >> 3;
	G_16z = (G_16z * 7 + imu_data_5) >> 3;


	mpu_flag = 1;
}

char get_gyroXYZ_raw(short int *gyro)
{
	if (mpu_flag)
	{
		mpu_flag = 0;
		gyro[X] = G_16x;
		gyro[Y] = G_16y;
		gyro[Z] = G_16z;
		return 1;
	}
	return 0;
}

char get_accelXYZ_raw(short int *accel)
{
	if (mpu_flag)
	{
		mpu_flag = 0;
		accel[X] = A_16x;
		accel[Y] = A_16y;
		accel[Z] = A_16z;
		return 1;
	}
	return 0;
}

void update_imu_raw(void)
{
	static int16_t accel_x,accel_y,accel_z;
	static int16_t gyro_x,gyro_y,gyro_z;

	if (mpu_flag)
	{ // first copy the imu data then process
		accel_x = A_16x;
		accel_y = A_16y;
		accel_z = A_16z;
		gyro_x = G_16x;
		gyro_y = G_16y;
		gyro_z = G_16z;

		mpu_flag = 0;
		/*Considering accel and gyro have same time stamp*/
		osMutexAcquire(imu_mutex, 1);
		process_gyro_data(gyro_x,gyro_y,gyro_z);
		process_accel_data(accel_x,accel_y,accel_z);

		osMutexRelease(imu_mutex);
		struct_imu_raw.time_us = get_us64();
		gpioG_togPin(status_led_blue);
	}
}
/*
update accel raw to g and mss
*/
void process_accel_data(const short int accel_x, const short int accel_y, const short int accel_z)
{
	float accl_correct[XYZ] = {0, 0, 0};
	static int idx;
	// raw data in g
	accl_correct[X] = accel_x * Res_Accl;
	accl_correct[Y] = accel_y * Res_Accl;
	accl_correct[Z] = accel_z * Res_Accl;

	accl_correct[X] -= struct_imu_cal.A_bias[X];
	accl_correct[Y] -= struct_imu_cal.A_bias[Y];
	accl_correct[Z] -= struct_imu_cal.A_bias[Z];

	mat_mul9_3(struct_imu_cal.acclTm, accl_correct, struct_imu_raw.Accl);
	// get_norm_xyz(struct_imu_raw.Accl,struct_imu_raw.Accl);

	m_accel_ring.x_data[idx] = struct_imu_raw.Accl[0] * -Actual_gravity;
	m_accel_ring.y_data[idx] = struct_imu_raw.Accl[1] * Actual_gravity;
	m_accel_ring.z_data[idx] = struct_imu_raw.Accl[2] * Actual_gravity;

	idx = (idx + 1) % imu_buf_len;

	system_state.accl_health = 1;
}
/*
update gyro raw to degrees and radians
*/
void process_gyro_data(const short int gyro_x, const short int gyro_y, const short int gyro_z)
{
	float gyro[XYZ] = {0, 0, 0};
	float filt_gyro[XYZ] = {0, 0, 0};
	static int idx;

	gyro[X] = gyro_x * Res_Gyro;
	gyro[Y] = gyro_y * Res_Gyro;
	gyro[Z] = gyro_z * Res_Gyro;

	gyro[X] -= struct_imu_cal.G_bias[X];
	gyro[Y] -= struct_imu_cal.G_bias[Y];
	gyro[Z] -= struct_imu_cal.G_bias[Z];

	run_dynamic_notch(gyro, filt_gyro);

	doResonantLPF(&lpf_gyro, filt_gyro, struct_imu_raw.Gyro);

	m_gyro_ring.x_data[idx] = radians_f(filt_gyro[X]);
	m_gyro_ring.y_data[idx] = radians_f(-filt_gyro[Y]);
	m_gyro_ring.z_data[idx] = radians_f(filt_gyro[Z]);

	idx = (idx + 1) % imu_buf_len;

	system_state.gyro_health = 1;
}
/*
update mag raw to micro tesla
*/

void process_mag_data(const short int mag_x, const short int mag_y, const short int mag_z)
{
	float mag_correct[3] = {0, 0, 0};
	float mag[3] = {0, 0, 0};
	static int idx;

	mag[0] = (float)mag_x * 0.01f;
	mag[1] = (float)mag_y * 0.01f;
	mag[2] = (float)mag_z * 0.01f;

	mag_correct[X] = mag[X] - struct_mag_cal.M_bias[X];
	mag_correct[Y] = mag[Y] - struct_mag_cal.M_bias[Y];
	mag_correct[Z] = mag[Z] - struct_mag_cal.M_bias[Z];

	// A*(m-b)
	osMutexAcquire(mag_mutex, 5);
	mat_mul9_3(struct_mag_cal.magTm, mag_correct, struct_mag_raw.mag_ut);

	m_mag_ring.x_data[idx] = struct_mag_raw.mag_ut[X];
	m_mag_ring.y_data[idx] = struct_mag_raw.mag_ut[Y];
	m_mag_ring.z_data[idx] = struct_mag_raw.mag_ut[Z];

	idx = (idx + 1) % 2;

	osMutexRelease(mag_mutex);
	struct_mag_raw.time_us = get_us64();
	system_state.mag_health = 1;
}

void process_air_data(const float P_mbr, const float T_deg)
{
  static float baro_dt;
	struct_air_raw.t_deg += 0.999f * (T_deg - struct_air_raw.t_deg);
	/* Add lpf to raw pressure measurements */
	struct_air_raw.p_mbr = doMovingAvgLPF(&lpf_press, P_mbr);
  
  uint64_t now = get_us64();
 // float dt = (float)(now - struct_air_raw.time_us) * 1e-6f;
	struct_air_raw.time_us = now;
//	
//  baro_dt += 0.999f * (dt - baro_dt);
//  if(baro_dt<0.012f)
//    struct_air_raw.p_mbr = sim_baro_pressure(sim_alt_rate, baro_dt);
}

float sim_baro_pressure(float vel_z_mps, float dt)
{
   const float P0 = 1013.25f;      // sea-level pressure, mbar
   const float T0 = 288.15f;       // sea-level temperature, K
   const float L  = 0.0065f;       // temperature lapse rate, K/m
   static float alt_m;
  
   float altitude_m = alt_m + vel_z_mps * dt;

   if (altitude_m < 0.0f)
   {
       altitude_m = 0.0f;
       alt_m = 0.0f;
   }
   if (altitude_m > 11000.0f)
   {
       altitude_m = 11000.0f;
       alt_m = 11000.0f;
   }
   float pressure_mbar = P0 * powf(1.0f - (L * altitude_m) / T0, 5.255787740552169866026f);
   
   alt_m += vel_z_mps * dt;
   
   return pressure_mbar;
}

void calculate_alt_data()
{
	static float alt_next;
	static float baro_sample_rate = 0.01f;
	static uint64_t baro_us;
	/*should be more than 80 hz*/
	static float alt_rate, altitude;

	uint64_t time_diff = struct_air_raw.time_us - baro_us;
	baro_us = struct_air_raw.time_us;

	if (time_diff > 7000U && time_diff < 12000U)
	{
		float _dt = time_diff * 1.0e-6f;

		// lpf to baro_sample_rate
		baro_sample_rate += 0.9995f * (_dt - baro_sample_rate);

		const float a = -0.0065f; // temperature laps rate
		// current pressure at MSL in kPa
		const float p1 = struct_air_raw.ref_press * 0.1f;
		// measured pressure in kPa
		const float p = struct_air_raw.p_mbr * 0.1f;
		// measured temperature in k
		const float T1 = struct_air_raw.t_deg + 273.15f;

		altitude = (((powf((p / p1), (-(a * AIR_GAS_CONST) * Gravity_inv))) * T1) - T1) / a;
		alt_rate = (altitude - alt_next) / fmaxf(baro_sample_rate,1.0e-4f);
		alt_next = altitude;
	}

	osMutexAcquire(baro_mutex, 1);
	struct_air_raw.alt_mtr = altitude;
	struct_air_raw.alt_rate = alt_rate;
	osMutexRelease(baro_mutex);
}

void calculate_as_data()
{
}

void get_angular_accel(const float *gyro_rad, float *accel_rate,const float dt)
{
	static float last_gyro[3] __attribute__((section(".DTCM_RAM")));
	static float filt_ang_accel[3] __attribute__((section(".DTCM_RAM")));
	float inv_dt = 1.0f / dt;

	float ang_accel_x = (gyro_rad[0] - last_gyro[0]) * inv_dt;
	float ang_accel_y = (gyro_rad[1] - last_gyro[1]) * inv_dt;
	float ang_accel_z = (gyro_rad[2] - last_gyro[2]) * inv_dt;

	filt_ang_accel[0] += 0.82f * (ang_accel_x - filt_ang_accel[0]);
	filt_ang_accel[1] += 0.82f * (ang_accel_y - filt_ang_accel[1]);
	filt_ang_accel[2] += 0.82f * (ang_accel_z - filt_ang_accel[2]);

	last_gyro[0] = gyro_rad[0];
	last_gyro[1] = gyro_rad[1];
	last_gyro[2] = gyro_rad[2];

	accel_rate[0] = filt_ang_accel[0];
	accel_rate[1] = filt_ang_accel[1];
	accel_rate[2] = filt_ang_accel[2];
}

void imu_lever_arm_correction(const float *accel_rate,const float *gyro, float *accel_cg)
{
	float r_off[3];
	float w_x_r[3];

	r_off[0] = (float)(Init_Cfg.imu_x_offset_mm+imu_board_X) * 0.001f;
	r_off[1] = (float)(Init_Cfg.imu_y_offset_mm+imu_board_Y) * 0.001f;
	r_off[2] = (float)(Init_Cfg.imu_z_offset_mm+imu_board_Z) * 0.001f;

	w_x_r[0] = gyro[1] * r_off[2] - gyro[2] * r_off[1];
	w_x_r[1] = gyro[2] * r_off[0] - gyro[0] * r_off[2];
	w_x_r[2] = gyro[0] * r_off[1] - gyro[1] * r_off[0];

	// Cross product: w x (w x r)
	float centripetal_x = gyro[1] * w_x_r[2] - gyro[2] * w_x_r[1];
	float centripetal_y = gyro[2] * w_x_r[0] - gyro[0] * w_x_r[2];
	float centripetal_z = gyro[0] * w_x_r[1] - gyro[1] * w_x_r[0];

	// Tangential Term: (dw/dt) x r
	// gyro_accel is the angular acceleration (rad/s^2)
	float tangential_x = accel_rate[1] * r_off[2] - accel_rate[2] * r_off[1];
	float tangential_y = accel_rate[2] * r_off[0] - accel_rate[0] * r_off[2];
	float tangential_z = accel_rate[0] * r_off[1] - accel_rate[1] * r_off[0];

	// Apply Correction: a_cg = a_imu - centripetal - tangential
	accel_cg[0] -= (centripetal_x + tangential_x);
	accel_cg[1] -= (centripetal_y + tangential_y);
	accel_cg[2] -= (centripetal_z + tangential_z);
}

char mpu9250_init(void)
{
	char temp[10];

	uint8_t wmi = 0;

	write_imu_reg(PWR_MGMNT_1, CLOCK_SEL_PLL);

	write_imu_reg(USER_CTRL, I2C_MST_EN);

	write_imu_reg(PWR_MGMNT_1, PWR_RESET);

	osDelay(40);
	write_imu_reg(PWR_MGMNT_1, CLOCK_SEL_PLL);

	wmi = read_imu_reg(MPU_WHOAMI);
	osDelay(10);

	if (wmi == 0x73 || wmi == 0x71)
	{
		write_imu_reg(PWR_MGMNT_2, SEN_ENABLE);

		write_imu_reg(GYRO_CONFIG, GYRO_FS_SEL_250DPS);
		write_imu_reg(CONFIG, GYRO_DLPF_99);

		write_imu_reg(ACCEL_CONFIG, ACCEL_FS_SEL_16G);
		write_imu_reg(ACCEL_CONFIG2, ACCEL_DLPF_99);

		write_imu_reg(SMPDIV, SMPL_RATE);

		write_imu_reg(INT_PIN_CFG, INT_PULSE_50US | 0x10);
		write_imu_reg(INT_ENABLE, INT_RAW_RDY_EN);
		osDelay(10);

		system_state.gyro_health = 1;
		system_state.accl_health = 1;

		print_Debug("MPU9250 OK\n");
		return 1;
	}
	else
	{
		dto_str(temp, wmi, '\n');
		print_Debug("mpu9250 who am i = ");
		print_Debug(temp);
		print_Debug("Can't communicate with mpu9250\n");
		system_state.gyro_health = 0;
		system_state.accl_health = 0;
		return 0;
	}
}

void write_imu_reg(unsigned char reg, unsigned char data)
{
	GPIOE->ODR &= ~spi4_cs;

	spi4_transfer(reg, 0x100);
	spi4_transfer(data, 0x100);
	GPIOE->ODR |= spi4_cs;

	osDelay(5);
}

unsigned char read_imu_reg(unsigned char reg)
{
	uint8_t x;
	GPIOE->ODR &= ~spi4_cs;
	spi4_transfer(reg | 0x80, 0x100);
	x = spi4_transfer(0x00, 0x100);
	GPIOE->ODR |= spi4_cs;
	return x;
}

char ms4525_init(void)
{
	return 1;
}

void imu_init(void)
{
	uint8_t c = try_counter;
	lpf_gyro.c1 = 0.165f;
	lpf_gyro.c2 = 1.0e-10f;
	lpf_alt.c1 = 0.185f;
	lpf_alt.c2 = 1.0e-10f;
	lpf_alt_rate.c1 = 0.00985f;
	lpf_alt_rate.c2 = 1.0e-10f;

	/*
	at cg = > 0 0 0
	fw from cg => +ve
	right from cg => +ve
	down from cg => +ve
	*/
//	Init_Cfg.imu_x_offset_mm = +10;
//	Init_Cfg.imu_y_offset_mm = 0;
//	Init_Cfg.imu_z_offset_mm = -30;

	imu_mutex = osMutexNew(NULL);
	mag_mutex = osMutexNew(NULL);
	baro_mutex = osMutexNew(NULL);

	resister_exti_callback(imu_exti_callback, 10);
	DMA2_Stream_Callback(accumulate_imu_data, 3);
	HAL_NVIC_EnableIRQ(EXTI15_10_IRQn);

	if (mpu9250_init() == 1)
	{
		spi4_init(4, dma_enable);
		osDelay(10);
		gpio_exti_init(GPIOE, gpio_pin_10, reset);
		system_state.sens_fault = 0;
		dynamic_notch_init();
		print_Debug("IMU Init Done\n");
	}
	else
	{
		while (c && !(mpu9250_init() == 1))
		{
			gpioG_setPin(status_led_error);
			print_Debug("Can't Init the IMU\n");
			system_state.sens_fault = 1;
			c--;
			osDelay(100);
		}
		if (c == 0)
			HAL_NVIC_SystemReset();
	}

	c = try_counter;

	if (ms4525_init() == 1)
	{
		osDelay(2);
		print_Debug("AS Init Done\n");
		system_state.sens_fault = 0;
	}
	else
	{
		while (!(ms4525_init() == 1) && c--)
		{
			gpioG_setPin(status_led_error);
			print_Debug("Can't Init the AS\n");
			system_state.sens_fault = 1;
			osDelay(100);
		}
	}
	if (c == 0)
		return;
}

void get_mag_ut(float *data)
{
	if (osMutexAcquire(mag_mutex, 1) == osOK)
	{
		float sum_x = 0;
		float sum_y = 0;
		float sum_z = 0;
		for (int i = 0; i < 2; i++)
		{
			sum_x += m_mag_ring.x_data[i];
			sum_y += m_mag_ring.y_data[i];
			sum_z += m_mag_ring.z_data[i];
		}
		data[0] = sum_y * 0.5f;
		data[1] = -sum_x * 0.5f;
		data[2] = sum_z * 0.5f;

		osMutexRelease(mag_mutex);
	}
}
void get_accel_mss(float *data)
{
	if (osMutexAcquire(imu_mutex, 1) == osOK)
	{ // axis changed for NED allignment

		float sum_x = 0;
		float sum_y = 0;
		float sum_z = 0;
		for (int i = 0; i < imu_buf_len; i++)
		{
			sum_x += m_accel_ring.x_data[i];
			sum_y += m_accel_ring.y_data[i];
			sum_z += m_accel_ring.z_data[i];
		}
		data[0] = sum_x * 0.2f;
		data[1] = sum_y * 0.2f;
		data[2] = sum_z * 0.2f;
		osMutexRelease(imu_mutex);
	}
}
void get_gyro_rad(float *data)
{
	if (osMutexAcquire(imu_mutex, 1) == osOK)
	{ // axis changed for NED allignment
		float sum_x = 0;
		float sum_y = 0;
		float sum_z = 0;
		for (int i = 0; i < imu_buf_len; i++)
		{
			sum_x += m_gyro_ring.x_data[i];
			sum_y += m_gyro_ring.y_data[i];
			sum_z += m_gyro_ring.z_data[i];
		}
		data[0] = sum_x * 0.2f;
		data[1] = sum_y * 0.2f;
		data[2] = sum_z * 0.2f;
		osMutexRelease(imu_mutex);
	}
}
void get_gyro_deg(float *data)
{
	if (osMutexAcquire(imu_mutex, 0) == osOK)
	{
		data[0] = struct_imu_raw.Gyro[0];
		data[1] = struct_imu_raw.Gyro[1];
		data[2] = struct_imu_raw.Gyro[2];
		osMutexRelease(imu_mutex);
	}
}

void set_baro_press_ref(const float p)
{
	struct_air_raw.ref_press = p;
}
float get_baro_press_ref(void)
{
	return struct_air_raw.ref_press;
}

float get_temperature_C(void)
{
	static float temp;
	if (osMutexAcquire(baro_mutex, 5) == osOK)
	{
		temp = struct_air_raw.t_deg;
		osMutexRelease(baro_mutex);
	}
	return temp;
}
float get_baro_pressure(void)
{
	static float press;
	if (osMutexAcquire(baro_mutex, 5) == osOK)
	{
		press = struct_air_raw.p_mbr;
		osMutexRelease(baro_mutex);
	}
	return press;
}
float get_baro_alt(void)
{
	static float altitude;
	if (osMutexAcquire(baro_mutex, 1) == osOK)
	{
		altitude = doResonantLPF_1D(&lpf_alt, struct_air_raw.alt_mtr);
		osMutexRelease(baro_mutex);
	}
	return altitude;
}
float get_alt_rate(void)
{
	static float alt_rate;
	if (osMutexAcquire(baro_mutex, 1) == osOK)
	{
		alt_rate = doResonantLPF_1D(&lpf_alt_rate, struct_air_raw.alt_rate);
		osMutexRelease(baro_mutex);
	}
	return alt_rate;
}
float get_air_speed(void)
{
	return struct_air_raw.air_speed;
}

__NO_RETURN void task_gyro_fft(void *argument)
{
	(void)argument;
	/*use to monitor the gyro frequency components
	  notch filter centre frequency */
	static uint64_t us;
	while (1)
	{
		uint64_t now = get_us64();
		dt_gfft = (now - us) * 1.0e-6f;
		us = now;

		check_notch_center_freq();
		osDelay(120);
	}
}
