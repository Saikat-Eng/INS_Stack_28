#ifndef _ZUPT_H
#define _ZUPT_H

#include "INSEstimator.h"

extern char zupt_detector(const float *accel, const float *gyro, const float *gps_vel, float dt, char position_hold_mode, char landed);
extern void INS_Filt_fusezupt(insfilterAsync *obj);


#endif
