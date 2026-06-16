/**
 * @file INSFuserBase.c 
 * @brief Implementation of the INS EKF fusion base functions.
 *
 * This file contains the core EKF fusion functions for integrating inertial and GPS measurements.  
 * 
 * @author saikat
 * @date 2023-12-01
 * @version 1.0
 */

#include "INSEstimator.h"
#include "INSFuserBase.h"
#include "Ekf.h"
#include "GeoBase.h"
#include "matrix.h"
#include "util.h"

/* Function Declarations */
static void accelMeans(const float *x, float *z);
static void accelJacobian(const float *x, float *dhdx);
static void magMeans(const float *x, float *z);
static void stateTransJacobianFcn(const float *x, float *dfdx);
static void magJakobian(const float *x, float *dhdx);
static void stateTransFcn(float *x);
static void ContinuousEKFPredictor_x(const float *x, const float *xdot, float *e, float dt);
static void ContinuousEKFPredictor_p(const float *P, const float *dfdx, float *pStateCovariance, float dt);

static void AccelEKF_equationInnovation(const float *P, const float *h,
                                        const float *H, const float *z,
                                        const float *R, float *innov,
                                        float *innovCov);
static void GyroEKF_equationInnovation(const float *P, const float *h,
                                       const float *H, const float *z,
                                       const float *R, float *innov,
                                       float *innovCov);
static void MagEKF_equationInnovation(const float *P, const float *h,
                                      const float *H, const float *z,
                                      const float *R, float *innov,
                                      float *innovCov);

static void GpsEKF_equationInnovation(const float *P, const float *h,
                                      const float *z, const float *R,
                                      float *innov, float *innovCov);


                                   
/* Function Definitions */
void AccelEKF_equationInnovation(const float *P, const float *h,
                                 const float *H, const float *z,
                                 const float *R, float *innov,
                                 float *innovCov)
{
    const float b_P_0=P[0]*H[0]+P[28]*H[1]+P[56]*H[2]+P[84]*H[3]+P[364]*H[13]+P[392]*H[14]+P[420]*H[15]+P[448]*1.0f;
    const float b_P_28=P[0]*H[28]+P[28]*H[29]+P[56]*H[30]+P[84]*H[31]+P[364]*H[41]+P[392]*H[42]+P[420]*H[43]+P[476]*1.0f;
    const float b_P_56=P[0]*H[56]+P[28]*H[57]+P[56]*H[58]+P[84]*H[59]+P[364]*H[69]+P[392]*H[70]+P[420]*H[71]+P[504]*1.0f;
    const float b_P_1=P[1]*H[0]+P[29]*H[1]+P[57]*H[2]+P[85]*H[3]+P[365]*H[13]+P[393]*H[14]+P[421]*H[15]+P[449]*H[16]*1.0f;
    const float b_P_29=P[1]*H[28]+P[29]*H[29]+P[57]*H[30]+P[85]*H[31]+P[365]*H[41]+P[393]*H[42]+P[421]*H[43]+P[477]*1.0f;
    const float b_P_57=P[1]*H[56]+P[29]*H[57]+P[57]*H[58]+P[85]*H[59]+P[365]*H[69]+P[393]*H[70]+P[421]*H[71]+P[505]*1.0f;
    const float b_P_2=P[2]*H[0]+P[30]*H[1]+P[58]*H[2]+P[86]*H[3]+P[366]*H[13]+P[394]*H[14]+P[422]*H[15]+P[450]*H[16]*1.0f;
    const float b_P_30=P[2]*H[28]+P[30]*H[29]+P[58]*H[30]+P[86]*H[31]+P[366]*H[41]+P[394]*H[42]+P[422]*H[43]+P[478]*1.0f;
    const float b_P_58=P[2]*H[56]+P[30]*H[57]+P[58]*H[58]+P[86]*H[59]+P[366]*H[69]+P[394]*H[70]+P[422]*H[71]+P[506]*1.0f;
    const float b_P_3=P[3]*H[0]+P[31]*H[1]+P[59]*H[2]+P[87]*H[3]+P[367]*H[13]+P[395]*H[14]+P[423]*H[15]+P[451]*H[16]*1.0f;
    const float b_P_31=P[3]*H[28]+P[31]*H[29]+P[59]*H[30]+P[87]*H[31]+P[367]*H[41]+P[395]*H[42]+P[423]*H[43]+P[479]*1.0f;
    const float b_P_59=P[3]*H[56]+P[31]*H[57]+P[59]*H[58]+P[87]*H[59]+P[367]*H[69]+P[395]*H[70]+P[423]*H[71]+P[507]*1.0f;
    const float b_P_13=P[13]*H[0]+P[41]*H[1]+P[69]*H[2]+P[97]*H[3]+P[377]*H[13]+P[405]*H[14]+P[433]*H[15]+P[461]*1.0f;
    const float b_P_41=P[13]*H[28]+P[41]*H[29]+P[69]*H[30]+P[97]*H[31]+P[377]*H[41]+P[405]*H[42]+P[433]*H[43]+P[489]*1.0f;
    const float b_P_69=P[13]*H[56]+P[41]*H[57]+P[69]*H[58]+P[97]*H[59]+P[377]*H[69]+P[405]*H[70]+P[433]*H[71]+P[517]*1.0f;
    const float b_P_14=P[14]*H[0]+P[42]*H[1]+P[70]*H[2]+P[98]*H[3]+P[378]*H[13]+P[406]*H[14]+P[434]*H[15]+P[462]*1.0f;
    const float b_P_42=P[14]*H[28]+P[42]*H[29]+P[70]*H[30]+P[98]*H[31]+P[378]*H[41]+P[406]*H[42]+P[434]*H[43]+P[490]*1.0f;
    const float b_P_70=P[14]*H[56]+P[42]*H[57]+P[70]*H[58]+P[98]*H[59]+P[378]*H[69]+P[406]*H[70]+P[434]*H[71]+P[518]*1.0f;
    const float b_P_15=P[15]*H[0]+P[43]*H[1]+P[71]*H[2]+P[99]*H[3]+P[379]*H[13]+P[407]*H[14]+P[435]*H[15]+P[463]*1.0f;
    const float b_P_43=P[15]*H[28]+P[43]*H[29]+P[71]*H[30]+P[99]*H[31]+P[379]*H[41]+P[407]*H[42]+P[435]*H[43]+P[491]*1.0f;
    const float b_P_71=P[15]*H[56]+P[43]*H[57]+P[71]*H[58]+P[99]*H[59]+P[379]*H[69]+P[407]*H[70]+P[435]*H[71]+P[519]*1.0f;
    const float b_P_16=P[16]*H[0]+P[44]*H[1]+P[72]*H[2]+P[100]*H[3]+P[380]*H[13]+P[408]*H[14]+P[436]*H[15]+P[464]*1.0f;
    const float b_P_44=P[16]*H[28]+P[44]*H[29]+P[72]*H[30]+P[100]*H[31]+P[380]*H[41]+P[408]*H[42]+P[436]*H[43]+P[492]*1.0f;
    const float b_P_72=P[16]*H[56]+P[44]*H[57]+P[72]*H[58]+P[100]*H[59]+P[380]*H[69]+P[408]*H[70]+P[436]*H[71]+P[520]*1.0f;
    const float b_P_17=P[17]*H[0]+P[45]*H[1]+P[73]*H[2]+P[101]*H[3]+P[381]*H[13]+P[409]*H[14]+P[437]*H[15]+P[465]*1.0f;
    const float b_P_45=P[17]*H[28]+P[45]*H[29]+P[73]*H[30]+P[101]*H[31]+P[381]*H[41]+P[409]*H[42]+P[437]*H[43]+P[493]*1.0f;
    const float b_P_73=P[17]*H[56]+P[45]*H[57]+P[73]*H[58]+P[101]*H[59]+P[381]*H[69]+P[409]*H[70]+P[437]*H[71]+P[521]*1.0f;
    const float b_P_18=P[18]*H[0]+P[46]*H[1]+P[74]*H[2]+P[102]*H[3]+P[382]*H[13]+P[410]*H[14]+P[438]*H[15]+P[466]*1.0f;
    const float b_P_46=P[18]*H[28]+P[46]*H[29]+P[74]*H[30]+P[102]*H[31]+P[382]*H[41]+P[410]*H[42]+P[438]*H[43]+P[494]*1.0f;
    const float b_P_74=P[18]*H[56]+P[46]*H[57]+P[74]*H[58]+P[102]*H[59]+P[382]*H[69]+P[410]*H[70]+P[438]*H[71]+P[522]*1.0f;


    //    for (b_i = 0; b_i < 3; b_i++) {
    //        for (i = 0; i < 3; i++) {
    //            d = 0.0f;
    //            for (i1 = 0; i1 < 28; i1++) {
    //                d += H[i1 + 28 * b_i] * b_P[i1 + 28 * i];
    //            }
    //            innovCov_tmp = b_i + 3 * i;
    //            innovCov[innovCov_tmp] = d + R[innovCov_tmp];
    //        }
    //        innov[b_i] = z[b_i] - h[b_i];
    //    }

    innovCov[0] = H[0] * b_P_0 +H[1] * b_P_1+H[2] * b_P_2+H[3] * b_P_3+H[13] * b_P_13+H[14] * b_P_14+H[15] * b_P_15+H[16] * b_P_16+ R[0];
    innovCov[3] = H[0] * b_P_28+H[1] * b_P_29+H[2] * b_P_30+H[3] * b_P_31+H[13] * b_P_41+H[14] * b_P_42+H[15] * b_P_43+H[16] * b_P_44+ R[3];
    innovCov[6] = H[0] * b_P_56+H[1] * b_P_57+H[2] * b_P_58+H[3] * b_P_59+H[13] * b_P_69+H[14] * b_P_70+H[15] * b_P_71+H[16] * b_P_72+ R[6];
    innov[0] = constrain_f(z[0] - h[0], -10.0f, 10.0f);
    innovCov[1] = H[28] * b_P_0 +H[29] * b_P_1+H[30] * b_P_2+H[31] * b_P_3+H[41] * b_P_13+H[42] * b_P_14+H[43] * b_P_15+H[45] * b_P_17+ R[1];
    innovCov[4] = H[28] * b_P_28+H[29] * b_P_29+H[30] * b_P_30+H[31] * b_P_31+H[41] * b_P_41+H[42] * b_P_42+H[43] * b_P_43+H[45] * b_P_45+ R[4];
    innovCov[7] = H[28] * b_P_56+H[29] * b_P_57+H[30] * b_P_58+H[31] * b_P_59+H[41] * b_P_69+H[42] * b_P_70+H[43] * b_P_71+H[45] * b_P_73+ R[7];
    innov[1] = constrain_f(z[1] - h[1], -10.0f, 10.0f);
    innovCov[2] = H[56] * b_P_0 +H[57] * b_P_1+H[58] * b_P_2+H[59] * b_P_3+H[69] * b_P_13+H[70] * b_P_14+H[71] * b_P_15+H[74] * b_P_18+ R[2];
    innovCov[5] = H[56] * b_P_28+H[57] * b_P_29+H[58] * b_P_30+H[59] * b_P_31+H[69] * b_P_41+H[70] * b_P_42+H[71] * b_P_43+H[74] * b_P_46+ R[5];
    innovCov[8] = H[56] * b_P_56+H[57] * b_P_57+H[58] * b_P_58+H[59] * b_P_59+H[69] * b_P_69+H[70] * b_P_70+H[71] * b_P_71+H[74] * b_P_74+ R[8];
    innov[2] = constrain_f(z[2] - h[2], -10.0f, 10.0f);

//    float PHt[28 * 3];
//    for (int i = 0; i < 28; i++) {
//        for (int j = 0; j < 3; j++) {
//            float sum = 0.0f;
//            for (int k = 0; k < 28; k++) {
//                sum += P[i + 28*k] * H[k + 28*j];
//            }
//            PHt[i + 28*j] = sum;
//        }
//    }

//    for (int i = 0; i < 3; i++) {
//        innov[i] = z[i] - h[i];
//        for (int j = 0; j < 3; j++) {
//            float sum = R[i + 3*j];
//            for (int k = 0; k < 28; k++) {
//                sum += H[k + 28*i] * PHt[k + 28*j];
//            }
//            innovCov[i + 3*j] = sum;
//        }
//    }
}

void GyroEKF_equationInnovation(const float *P, const float *h,
                                const float *H, const float *z,
                                const float *R, float *innov,
                                float *innovCov)
{
  //  float b_P[84];
    //    float d;
    //    //    register int b_i;
    //    int i;
    //    int i1;
    //    int innovCov_tmp;

    //    for (i = 0; i < 28; i++) {
    //        for (i1 = 0; i1 < 3; i1++) {
    //            d = 0.0f;
    //            for (innovCov_tmp = 0; innovCov_tmp < 28; innovCov_tmp++) {
    //                d += P[i + 28 * innovCov_tmp] * H[innovCov_tmp + 28 * i1];
    //            }
    //            b_P[i + 28 * i1] = d;
    //        }
    //    }

    const float b_P_4=P[116]*H[4]+P[536]*H[19];
    const float b_P_32=P[144]*H[33]+P[564]*H[48];
    const float b_P_60=P[172]*H[62]+P[592]*H[77];
    const float b_P_5=P[117]*H[4]+P[537]*H[19];
    const float b_P_33=P[145]*H[33]+P[565]*H[48];
    const float b_P_61=P[173]*H[62]+P[593]*H[77];
    const float b_P_6=P[118]*H[4]+P[538]*H[19];
    const float b_P_34=P[146]*H[33]+P[566]*H[48];
    const float b_P_62=P[174]*H[62]+P[594]*H[77];
    const float b_P_19=P[131]*H[4]+P[551]*H[19];
    const float b_P_47=P[159]*H[33]+P[579]*H[48];
    const float b_P_75=P[187]*H[62]+P[607]*H[77];
    const float b_P_20=P[132]*H[4]+P[552]*H[19];
    const float b_P_48=P[160]*H[33]+P[580]*H[48];
    const float b_P_76=P[188]*H[62]+P[608]*H[77];
    const float b_P_21=P[133]*H[4]+P[553]*H[19];
    const float b_P_49=P[161]*H[33]+P[581]*H[48];
    const float b_P_77=P[189]*H[62]+P[609]*H[77];



    //    for (b_i = 0; b_i < 3; b_i++) {
    //        for (i = 0; i < 3; i++) {
    //            d = 0.0f;
    //            for (i1 = 0; i1 < 28; i1++) {
    //                d += H[i1 + 28 * b_i] * b_P[i1 + 28 * i];
    //            }
    //            innovCov_tmp = b_i + 3 * i;
    //            innovCov[innovCov_tmp] = d + R[innovCov_tmp];
    //        }
    //        innov[b_i] = z[b_i] - h[b_i];
    //    }

    innovCov[0] = H[4] * b_P_4   + H[19] * b_P_19 + R[0];
    innovCov[3] = H[4] * b_P_32  + H[19] * b_P_47 + R[3];
    innovCov[6] = H[4] * b_P_60  + H[19] * b_P_75 + R[6];
    innov[0] = z[0] - h[0];
    innovCov[1] = H[33] * b_P_5  + H[48] * b_P_20 + R[1];
    innovCov[4] = H[33] * b_P_33 + H[48] * b_P_48 + R[4];
    innovCov[7] = H[33] * b_P_61 + H[48] * b_P_76 + R[7];
    innov[1] = z[1] - h[1];
    innovCov[2] = H[62] * b_P_6  + H[77] * b_P_21 + R[2];
    innovCov[5] = H[62] * b_P_34 + H[77] * b_P_49 + R[5];
    innovCov[8] = H[62] * b_P_62 + H[77] * b_P_77 + R[8];
    innov[2] = z[2] - h[2];
}

void MagEKF_equationInnovation(const float *P, const float *h,
                               const float *H, const float *z,
                               const float *R, float *innov,
                               float *innovCov)
{
    //float b_P[84];
    //    float d;
    //    //    register int b_i;
    //    int i;
    //    int i1;
    //    int innovCov_tmp;

    //    for (i = 0; i < 28; i++) {
    //        for (i1 = 0; i1 < 3; i1++) {
    //            d = 0.0f;
    //            for (innovCov_tmp = 0; innovCov_tmp < 28; innovCov_tmp++) {
    //                d += P[i + 28 * innovCov_tmp] * H[innovCov_tmp + 28 * i1];
    //            }
    //            b_P[i + 28 * i1] = d;
    //        }
    //    }
    const float b_P_0=P[0]*H[0]+P[28]*H[1]+P[56]*H[2]+P[84]*H[3]+P[616]*H[22]+P[644]*H[23]+P[672]*H[24]+P[700]*1.0f;
    const float b_P_28=P[0]*H[28]+P[28]*H[29]+P[56]*H[30]+P[84]*H[31]+P[616]*H[50]+P[644]*H[51]+P[672]*H[52]+P[728]*1.0f;
    const float b_P_56=P[0]*H[56]+P[28]*H[57]+P[56]*H[58]+P[84]*H[59]+P[616]*H[78]+P[644]*H[79]+P[672]*H[80]+P[756]*1.0f;
    const float b_P_1=P[1]*H[0]+P[29]*H[1]+P[57]*H[2]+P[85]*H[3]+P[617]*H[22]+P[645]*H[23]+P[673]*H[24]+P[701]*H[25]*1.0f;
    const float b_P_29=P[1]*H[28]+P[29]*H[29]+P[57]*H[30]+P[85]*H[31]+P[617]*H[50]+P[645]*H[51]+P[673]*H[52]+P[729]*1.0f;
    const float b_P_57=P[1]*H[56]+P[29]*H[57]+P[57]*H[58]+P[85]*H[59]+P[617]*H[78]+P[645]*H[79]+P[673]*H[80]+P[757]*1.0f;
    const float b_P_2=P[2]*H[0]+P[30]*H[1]+P[58]*H[2]+P[86]*H[3]+P[618]*H[22]+P[646]*H[23]+P[674]*H[24]+P[702]*H[25]*1.0f;
    const float b_P_30=P[2]*H[28]+P[30]*H[29]+P[58]*H[30]+P[86]*H[31]+P[618]*H[50]+P[646]*H[51]+P[674]*H[52]+P[730]*1.0f;
    const float b_P_58=P[2]*H[56]+P[30]*H[57]+P[58]*H[58]+P[86]*H[59]+P[618]*H[78]+P[646]*H[79]+P[674]*H[80]+P[758]*1.0f;
    const float b_P_3=P[3]*H[0]+P[31]*H[1]+P[59]*H[2]+P[87]*H[3]+P[619]*H[22]+P[647]*H[23]+P[675]*H[24]+P[703]*H[25]*1.0f;
    const float b_P_31=P[3]*H[28]+P[31]*H[29]+P[59]*H[30]+P[87]*H[31]+P[619]*H[50]+P[647]*H[51]+P[675]*H[52]+P[731]*1.0f;
    const float b_P_59=P[3]*H[56]+P[31]*H[57]+P[59]*H[58]+P[87]*H[59]+P[619]*H[78]+P[647]*H[79]+P[675]*H[80]+P[759]*1.0f;
    const float b_P_22=P[22]*H[0]+P[50]*H[1]+P[78]*H[2]+P[106]*H[3]+P[638]*H[22]+P[666]*H[23]+P[694]*H[24]+P[722]*1.0f;
    const float b_P_50=P[22]*H[28]+P[50]*H[29]+P[78]*H[30]+P[106]*H[31]+P[638]*H[50]+P[666]*H[51]+P[694]*H[52]+P[750]*1.0f;
    const float b_P_78=P[22]*H[56]+P[50]*H[57]+P[78]*H[58]+P[106]*H[59]+P[638]*H[78]+P[666]*H[79]+P[694]*H[80]+P[778]*1.0f;
    const float b_P_23=P[23]*H[0]+P[51]*H[1]+P[79]*H[2]+P[107]*H[3]+P[639]*H[22]+P[667]*H[23]+P[695]*H[24]+P[723]*1.0f;
    const float b_P_51=P[23]*H[28]+P[51]*H[29]+P[79]*H[30]+P[107]*H[31]+P[639]*H[50]+P[667]*H[51]+P[695]*H[52]+P[751]*1.0f;
    const float b_P_79=P[23]*H[56]+P[51]*H[57]+P[79]*H[58]+P[107]*H[59]+P[639]*H[78]+P[667]*H[79]+P[695]*H[80]+P[779]*1.0f;
    const float b_P_24=P[24]*H[0]+P[52]*H[1]+P[80]*H[2]+P[108]*H[3]+P[640]*H[22]+P[668]*H[23]+P[696]*H[24]+P[724]*1.0f;
    const float b_P_52=P[24]*H[28]+P[52]*H[29]+P[80]*H[30]+P[108]*H[31]+P[640]*H[50]+P[668]*H[51]+P[696]*H[52]+P[752]*1.0f;
    const float b_P_80=P[24]*H[56]+P[52]*H[57]+P[80]*H[58]+P[108]*H[59]+P[640]*H[78]+P[668]*H[79]+P[696]*H[80]+P[780]*1.0f;
    const float b_P_25=P[25]*H[0]+P[53]*H[1]+P[81]*H[2]+P[109]*H[3]+P[641]*H[22]+P[669]*H[23]+P[697]*H[24]+P[725]*1.0f;
    const float b_P_53=P[25]*H[28]+P[53]*H[29]+P[81]*H[30]+P[109]*H[31]+P[641]*H[50]+P[669]*H[51]+P[697]*H[52]+P[753]*1.0f;
    const float b_P_81=P[25]*H[56]+P[53]*H[57]+P[81]*H[58]+P[109]*H[59]+P[641]*H[78]+P[669]*H[79]+P[697]*H[80]+P[781]*1.0f;
    const float b_P_26=P[26]*H[0]+P[54]*H[1]+P[82]*H[2]+P[110]*H[3]+P[642]*H[22]+P[670]*H[23]+P[698]*H[24]+P[726]*1.0f;
    const float b_P_54=P[26]*H[28]+P[54]*H[29]+P[82]*H[30]+P[110]*H[31]+P[642]*H[50]+P[670]*H[51]+P[698]*H[52]+P[754]*1.0f;
    const float b_P_82=P[26]*H[56]+P[54]*H[57]+P[82]*H[58]+P[110]*H[59]+P[642]*H[78]+P[670]*H[79]+P[698]*H[80]+P[782]*1.0f;
    const float b_P_27=P[27]*H[0]+P[55]*H[1]+P[83]*H[2]+P[111]*H[3]+P[643]*H[22]+P[671]*H[23]+P[699]*H[24]+P[727]*1.0f;
    const float b_P_55=P[27]*H[28]+P[55]*H[29]+P[83]*H[30]+P[111]*H[31]+P[643]*H[50]+P[671]*H[51]+P[699]*H[52]+P[755]*1.0f;
    const float b_P_83=P[27]*H[56]+P[55]*H[57]+P[83]*H[58]+P[111]*H[59]+P[643]*H[78]+P[671]*H[79]+P[699]*H[80]+P[783]*1.0f;
     
    //    for (b_i = 0; b_i < 3; b_i++) {
    //        for (i = 0; i < 3; i++) {
    //            d = 0.0f;
    //            for (i1 = 0; i1 < 28; i1++) {
    //                d += H[i1 + 28 * b_i] * b_P[i1 + 28 * i];
    //            }
    //            innovCov_tmp = b_i + 3 * i;
    //            innovCov[innovCov_tmp] = d + R[innovCov_tmp];
    //        }
    //        innov[b_i] = z[b_i] - h[b_i];
    //    }

    innovCov[0] = H[0] * b_P_0+H[1] * b_P_1+H[2] * b_P_2+H[3] * b_P_3+H[22] * b_P_22+H[23] * b_P_23+H[24] * b_P_24+H[25] * b_P_25+ R[0];
    innovCov[3] = H[0] * b_P_28+H[1] * b_P_29+H[2] * b_P_30+H[3] * b_P_31+H[22] * b_P_50+H[23] * b_P_51+H[24] * b_P_52+H[25] * b_P_53+ R[3];
    innovCov[6] = H[0] * b_P_56+H[1] * b_P_57+H[2] * b_P_58+H[3] * b_P_59+H[22] * b_P_78+H[23] * b_P_79+H[24] * b_P_80+H[25] * b_P_81+ R[6];
    innov[0] = z[0] - h[0];
    innovCov[1] = H[28] * b_P_0+H[29] * b_P_1+H[30] * b_P_2+H[31] * b_P_3+H[50] * b_P_22+H[51] * b_P_23+H[52] * b_P_24+H[54] * b_P_26+ R[1];
    innovCov[4] = H[28] * b_P_28+H[29] * b_P_29+H[30] * b_P_30+H[31] * b_P_31+H[50] * b_P_50+H[51] * b_P_51+H[52] * b_P_52+H[54] * b_P_54+ R[4];
    innovCov[7] = H[28] * b_P_56+H[29] * b_P_57+H[30] * b_P_58+H[31] * b_P_59+H[50] * b_P_78+H[51] * b_P_79+H[52] * b_P_80+H[54] * b_P_82+ R[7];
    innov[1] = z[1] - h[1];
    innovCov[2] = H[56] * b_P_0+H[57] * b_P_1+H[58] * b_P_2+H[59] * b_P_3+H[78] * b_P_22+H[79] * b_P_23+H[80] * b_P_24+H[83] * b_P_27+ R[2];
    innovCov[5] = H[56] * b_P_28+H[57] * b_P_29+H[58] * b_P_30+H[59] * b_P_31+H[78] * b_P_50+H[79] * b_P_51+H[80] * b_P_52+H[83] * b_P_55+ R[5];
    innovCov[8] = H[56] * b_P_56+H[57] * b_P_57+H[58] * b_P_58+H[59] * b_P_59+H[78] * b_P_78+H[79] * b_P_79+H[80] * b_P_80+H[83] * b_P_83+ R[8];
    innov[2] = z[2] - h[2];
}

static void GpsEKF_equationInnovation(const float *P, const float *h,
                                      const float *z, const float *R,
                                      float *innov, float *innovCov)
{
    //float b_P[168];
    //    float d;
    //    register int b_i;
    //    register int i;
    //    register int i1;
    //    register int innovCov_tmp;

    //    for (i = 0; i < 28; i++) {
    //        for (i1 = 0; i1 < 6; i1++) {
    //            d = 0.0f;
    //            for (innovCov_tmp = 0; innovCov_tmp < 28; innovCov_tmp++) {
    //                d += P[i + 28 * innovCov_tmp] * (float)iv[innovCov_tmp + 28 * i1];
    //            }
    //            b_P[i + 28 * i1] = d;
    //        }
    //    }

    //   for (b_i = 0; b_i < 6; b_i++) {
    //       for (i = 0; i < 6; i++) {
    //           d = 0.0f;
    //           for (i1 = 0; i1 < 28; i1++) {
    //               d += (float)iv1[b_i + 6 * i1] * b_P[i1 + 28 * i];
    //           }
    //           innovCov_tmp = b_i + 6 * i;
    //           innovCov[innovCov_tmp] = d + R[innovCov_tmp];
    //       }
    //       innov[b_i] = z[b_i] - h[b_i];
    //   }

    innovCov[0] = 1.0f*P[203]+ R[0];
    innovCov[6] = 1.0f*P[231]+ R[6];
    innovCov[12] = 1.0f*P[259]+ R[12];
    innovCov[18] = 1.0f*P[287]+ R[18];
    innovCov[24] = 1.0f*P[351]+ R[24];
    innovCov[30] = 1.0f*P[343]+ R[30];
    innov[0] = z[0] - h[0];
    innovCov[1] = 1.0f*P[204]+ R[1];
    innovCov[7] = 1.0f*P[232]+ R[7];
    innovCov[13] = 1.0f*P[260]+ R[13];
    innovCov[19] = 1.0f*P[288]+ R[19];
    innovCov[25] = 1.0f*P[316]+ R[25];
    innovCov[31] = 1.0f*P[344]+ R[31];
    innov[1] = z[1] - h[1];
    innovCov[2] = 1.0f*P[205]+ R[2];
    innovCov[8] = 1.0f*P[233]+ R[8];
    innovCov[14] = 1.0f*P[261]+ R[14];
    innovCov[20] = 1.0f*P[289]+ R[20];
    innovCov[26] = 1.0f*P[317]+ R[26];
    innovCov[32] = 1.0f*P[345]+ R[32];
    innov[2] = z[2] - h[2];
    innovCov[3] = 1.0f*P[206]+ R[3];
    innovCov[9] = 1.0f*P[234]+ R[9];
    innovCov[15] = 1.0f*P[262]+ R[15];
    innovCov[21] = 1.0f*P[290]+ R[21];
    innovCov[27] = 1.0f*P[318]+ R[27];
    innovCov[33] = 1.0f*P[346]+ R[33];
    innov[3] = z[3] - h[3];
    innovCov[4] = 1.0f*P[207]+ R[4];
    innovCov[10] = 1.0f*P[235]+ R[10];
    innovCov[16] = 1.0f*P[263]+ R[16];
    innovCov[22] = 1.0f*P[291]+ R[22];
    innovCov[28] = 1.0f*P[319]+ R[28];
    innovCov[34] = 1.0f*P[347]+ R[34];
    innov[4] = z[4] - h[4];
    innovCov[5] = 1.0f*P[208]+ R[5];
    innovCov[11] = 1.0f*P[236]+ R[11];
    innovCov[17] = 1.0f*P[264]+ R[17];
    innovCov[23] = 1.0f*P[292]+ R[23];
    innovCov[29] = 1.0f*P[320]+ R[29];
    innovCov[35] = 1.0f*P[348]+ R[35];
    innov[5] = z[5] - h[5];
}

char accel_ok_for_fusion(const float *accel, const float *gyro)
{
    const float a_norm = sqrtf(accel[0]*accel[0] + accel[1]*accel[1] + accel[2]*accel[2]);
    const float w_norm = sqrtf(gyro[0]*gyro[0] + gyro[1]*gyro[1] + gyro[2]*gyro[2]);
    return (fabsf(a_norm - Actual_gravity) < 1.0f) && (w_norm < 1.5f);
}

void INS_Filt_predict(insfilterAsync *obj, const float dt)
{
    static float *P;
    static float xdot[28] __attribute__((section(".DTCM_RAM")));
    static float dfdx[784] __attribute__((section(".DTCM_RAM")));
    static float pDot[784] __attribute__((section(".DTCM_RAM")));
    static float b_dfdx[784] __attribute__((section(".DTCM_RAM")));
    int i0;
    int i1;
    int i2;
    int dfdx_tmp;

    xdot[0]  = obj->pState[0];	xdot[1]  = obj->pState[1];	xdot[2]  = obj->pState[2];	xdot[3]  = obj->pState[3];
    xdot[4]  = obj->pState[4];	xdot[5]  = obj->pState[5];	xdot[6]  = obj->pState[6];	xdot[7]  = obj->pState[7];
    xdot[8]  = obj->pState[8];	xdot[9]  = obj->pState[9];	xdot[10] = obj->pState[10];	xdot[11] = obj->pState[11];
    xdot[12] = obj->pState[12];	xdot[13] = obj->pState[13];	xdot[14] = obj->pState[14];	xdot[15] = obj->pState[15];
    xdot[16] = obj->pState[16];	xdot[17] = obj->pState[17];	xdot[18] = obj->pState[18];	xdot[19] = obj->pState[19];
    xdot[20] = obj->pState[20];	xdot[21] = obj->pState[21];	xdot[22] = obj->pState[22];	xdot[23] = obj->pState[23];
    xdot[24] = obj->pState[24];	xdot[25] = obj->pState[25];	xdot[26] = obj->pState[26];	xdot[27] = obj->pState[27];

    P = obj->pStateCovariance;

    //state transition
    stateTransFcn(xdot);
    //state transition jacobian
    stateTransJacobianFcn(obj->pState, dfdx);

    for (i0 = 0; i0 < 28; i0++) {
        for (i1 = 0; i1 < 28; i1++) {
            float d0 = 0.0f;
            float d1 = 0.0f;
            for (dfdx_tmp = 0; dfdx_tmp < 28; dfdx_tmp++) {
                i2 = dfdx_tmp + 28 * i1;
                d0 += P[i0 + 28 * dfdx_tmp] * dfdx[i2];
                d1 += dfdx[dfdx_tmp + 28 * i0] * P[i2];
            }
            dfdx_tmp = i0 + 28 * i1;
            b_dfdx[dfdx_tmp] = d1;
            pDot[dfdx_tmp] = d0;
        }
    }

    //  for (i0 = 0; i0 < 784; i0++)
    //  {
    //    pDot[i0] = (pDot[i0] + b_dfdx[i0]) + obj->addProcNoise[i0];
    //  }

    for (i0 = 0; i0 < 28; i0++)
    {
        i1 = 0+i0;
        pDot[i1] = (pDot[i1] + b_dfdx[i1]) + obj->addProcNoise[i1];
        i1 = 28+i0;
        pDot[i1] = (pDot[i1] + b_dfdx[i1]) + obj->addProcNoise[i1];
        i1 = 56+i0;
        pDot[i1] = (pDot[i1] + b_dfdx[i1]) + obj->addProcNoise[i1];
        i1 = 84+i0;
        pDot[i1] = (pDot[i1] + b_dfdx[i1]) + obj->addProcNoise[i1];
        i1 = 112+i0;
        pDot[i1] = (pDot[i1] + b_dfdx[i1]) + obj->addProcNoise[i1];
        i1 = 140+i0;
        pDot[i1] = (pDot[i1] + b_dfdx[i1]) + obj->addProcNoise[i1];
        i1 = 168+i0;
        pDot[i1] = (pDot[i1] + b_dfdx[i1]) + obj->addProcNoise[i1];
        i1 = 196+i0;
        pDot[i1] = (pDot[i1] + b_dfdx[i1]) + obj->addProcNoise[i1];
        i1 = 224+i0;
        pDot[i1] = (pDot[i1] + b_dfdx[i1]) + obj->addProcNoise[i1];
        i1 = 252+i0;
        pDot[i1] = (pDot[i1] + b_dfdx[i1]) + obj->addProcNoise[i1];
        i1 = 280+i0;
        pDot[i1] = (pDot[i1] + b_dfdx[i1]) + obj->addProcNoise[i1];
        i1 = 308+i0;
        pDot[i1] = (pDot[i1] + b_dfdx[i1]) + obj->addProcNoise[i1];
        i1 = 336+i0;
        pDot[i1] = (pDot[i1] + b_dfdx[i1]) + obj->addProcNoise[i1];
        i1 = 364+i0;
        pDot[i1] = (pDot[i1] + b_dfdx[i1]) + obj->addProcNoise[i1];
        i1 = 392+i0;
        pDot[i1] = (pDot[i1] + b_dfdx[i1]) + obj->addProcNoise[i1];
        i1 = 420+i0;
        pDot[i1] = (pDot[i1] + b_dfdx[i1]) + obj->addProcNoise[i1];
        i1 = 448+i0;
        pDot[i1] = (pDot[i1] + b_dfdx[i1]) + obj->addProcNoise[i1];
        i1 = 476+i0;
        pDot[i1] = (pDot[i1] + b_dfdx[i1]) + obj->addProcNoise[i1];
        i1 = 504+i0;
        pDot[i1] = (pDot[i1] + b_dfdx[i1]) + obj->addProcNoise[i1];
        i1 = 532+i0;
        pDot[i1] = (pDot[i1] + b_dfdx[i1]) + obj->addProcNoise[i1];
        i1 = 560+i0;
        pDot[i1] = (pDot[i1] + b_dfdx[i1]) + obj->addProcNoise[i1];
        i1 = 588+i0;
        pDot[i1] = (pDot[i1] + b_dfdx[i1]) + obj->addProcNoise[i1];
        i1 = 616+i0;
        pDot[i1] = (pDot[i1] + b_dfdx[i1]) + obj->addProcNoise[i1];
        i1 = 644+i0;
        pDot[i1] = (pDot[i1] + b_dfdx[i1]) + obj->addProcNoise[i1];
        i1 = 672+i0;
        pDot[i1] = (pDot[i1] + b_dfdx[i1]) + obj->addProcNoise[i1];
        i1 = 700+i0;
        pDot[i1] = (pDot[i1] + b_dfdx[i1]) + obj->addProcNoise[i1];
        i1 = 728+i0;
        pDot[i1] = (pDot[i1] + b_dfdx[i1]) + obj->addProcNoise[i1];
        i1 = 756+i0;
        pDot[i1] = (pDot[i1] + b_dfdx[i1]) + obj->addProcNoise[i1];
    }

    //  for (i0 = 0; i0 < 28; i0++) {
    //    for (i1 = 0; i1 < 28; i1++) {
    //      dfdx_tmp = i1 + 28 * i0;
    //      b_dfdx[dfdx_tmp] = 0.5f * (pDot[dfdx_tmp] + pDot[i0 + 28 * i1]);
    //    }
    //  }
    // b_dfdx:pDot = 0.5 * (pDot + pDot.'); % ensure symmetry
    for (i0 = 0; i0 < 28; i0++)
    {
        dfdx_tmp = 0 + 28 * i0;
        b_dfdx[dfdx_tmp] = 0.5f * (pDot[dfdx_tmp] + pDot[i0 + 0]);
        dfdx_tmp = 1 + 28 * i0;
        b_dfdx[dfdx_tmp] = 0.5f * (pDot[dfdx_tmp] + pDot[i0 + 28]);
        dfdx_tmp = 2 + 28 * i0;
        b_dfdx[dfdx_tmp] = 0.5f * (pDot[dfdx_tmp] + pDot[i0 + 56]);
        dfdx_tmp = 3 + 28 * i0;
        b_dfdx[dfdx_tmp] = 0.5f * (pDot[dfdx_tmp] + pDot[i0 + 84]);
        dfdx_tmp = 4 + 28 * i0;
        b_dfdx[dfdx_tmp] = 0.5f * (pDot[dfdx_tmp] + pDot[i0 + 112]);
        dfdx_tmp = 5 + 28 * i0;
        b_dfdx[dfdx_tmp] = 0.5f * (pDot[dfdx_tmp] + pDot[i0 + 140]);
        dfdx_tmp = 6 + 28 * i0;
        b_dfdx[dfdx_tmp] = 0.5f * (pDot[dfdx_tmp] + pDot[i0 + 168]);
        dfdx_tmp = 7 + 28 * i0;
        b_dfdx[dfdx_tmp] = 0.5f * (pDot[dfdx_tmp] + pDot[i0 + 196]);
        dfdx_tmp = 8 + 28 * i0;
        b_dfdx[dfdx_tmp] = 0.5f * (pDot[dfdx_tmp] + pDot[i0 + 224]);
        dfdx_tmp = 9 + 28 * i0;
        b_dfdx[dfdx_tmp] = 0.5f * (pDot[dfdx_tmp] + pDot[i0 + 252]);
        dfdx_tmp = 10 + 28 * i0;
        b_dfdx[dfdx_tmp] = 0.5f * (pDot[dfdx_tmp] + pDot[i0 + 280]);
        dfdx_tmp = 11 + 28 * i0;
        b_dfdx[dfdx_tmp] = 0.5f * (pDot[dfdx_tmp] + pDot[i0 + 308]);
        dfdx_tmp = 12 + 28 * i0;
        b_dfdx[dfdx_tmp] = 0.5f * (pDot[dfdx_tmp] + pDot[i0 + 336]);
        dfdx_tmp = 13 + 28 * i0;
        b_dfdx[dfdx_tmp] = 0.5f * (pDot[dfdx_tmp] + pDot[i0 + 364]);
        dfdx_tmp = 14 + 28 * i0;
        b_dfdx[dfdx_tmp] = 0.5f * (pDot[dfdx_tmp] + pDot[i0 + 392]);
        dfdx_tmp = 15 + 28 * i0;
        b_dfdx[dfdx_tmp] = 0.5f * (pDot[dfdx_tmp] + pDot[i0 + 420]);
        dfdx_tmp = 16 + 28 * i0;
        b_dfdx[dfdx_tmp] = 0.5f * (pDot[dfdx_tmp] + pDot[i0 + 448]);
        dfdx_tmp = 17 + 28 * i0;
        b_dfdx[dfdx_tmp] = 0.5f * (pDot[dfdx_tmp] + pDot[i0 + 476]);
        dfdx_tmp = 18 + 28 * i0;
        b_dfdx[dfdx_tmp] = 0.5f * (pDot[dfdx_tmp] + pDot[i0 + 504]);
        dfdx_tmp = 19 + 28 * i0;
        b_dfdx[dfdx_tmp] = 0.5f * (pDot[dfdx_tmp] + pDot[i0 + 532]);
        dfdx_tmp = 20 + 28 * i0;
        b_dfdx[dfdx_tmp] = 0.5f * (pDot[dfdx_tmp] + pDot[i0 + 560]);
        dfdx_tmp = 21 + 28 * i0;
        b_dfdx[dfdx_tmp] = 0.5f * (pDot[dfdx_tmp] + pDot[i0 + 588]);
        dfdx_tmp = 22 + 28 * i0;
        b_dfdx[dfdx_tmp] = 0.5f * (pDot[dfdx_tmp] + pDot[i0 + 616]);
        dfdx_tmp = 23 + 28 * i0;
        b_dfdx[dfdx_tmp] = 0.5f * (pDot[dfdx_tmp] + pDot[i0 + 644]);
        dfdx_tmp = 24 + 28 * i0;
        b_dfdx[dfdx_tmp] = 0.5f * (pDot[dfdx_tmp] + pDot[i0 + 672]);
        dfdx_tmp = 25 + 28 * i0;
        b_dfdx[dfdx_tmp] = 0.5f * (pDot[dfdx_tmp] + pDot[i0 + 700]);
        dfdx_tmp = 26 + 28 * i0;
        b_dfdx[dfdx_tmp] = 0.5f * (pDot[dfdx_tmp] + pDot[i0 + 728]);
        dfdx_tmp = 27 + 28 * i0;
        b_dfdx[dfdx_tmp] = 0.5f * (pDot[dfdx_tmp] + pDot[i0 + 756]);

    }

    ContinuousEKFPredictor_x(obj->pState, xdot, obj->pState, dt);
    get_norm_quat(obj->pState);
    //update covarience
    ContinuousEKFPredictor_p(P, b_dfdx, obj->pStateCovariance, dt);
    
//    ekf_fix_covariance(P,28);
//    for (int r=0; r<28; ++r) 
//    {
//      for (int c=r+1; c<28; ++c) 
//      {   
//        //pStateCovariance[i1] = P[i1] + dt * dfdx[i1];
//        float s = 0.5f*(P[r + 28*c] + P[c + 28*r]);
//        P[r + 28*c] = P[c + 28*r] = s;
//      }
//      if (P[r + 28*r] < 1e-12f) P[r + 28*r] = 1e-12f;
//   }
}


/*
 * Arguments    : const float x[28]
 *                float z[3]
 * Return Type  : void
 */
void accelMeans(const float *x, float *z)
{  
    const float q0 = x[0];
    const float q1 = x[1];
    const float q2 = x[2];
    const float q3 = x[3];
    
    // Extract accelerations in body frame
    const float accx_b = x[16];
    const float accy_b = x[17];
    const float accz_b = x[18];
    
    // Precompute reused terms for optimization
    const float q0q0 = q0 * q0;
    const float q1q1 = q1 * q1;
    const float q2q2 = q2 * q2;
    const float q3q3 = q3 * q3;
    
    const float two_q0q1 = 2.0f * q0 * q1;
    const float two_q0q2 = 2.0f * q0 * q2;
    const float two_q0q3 = 2.0f * q0 * q3;
    const float two_q1q2 = 2.0f * q1 * q2;
    const float two_q1q3 = 2.0f * q1 * q3;
    const float two_q2q3 = 2.0f * q2 * q3;
    
    const float an_gx = x[13];
    const float ae_gy = x[14];
    const float ad_gz = x[15] - Actual_gravity;

    // Compute z vector components
    z[0] = accx_b - an_gx * (q0q0 + q1q1 - q2q2 - q3q3) + ad_gz * (two_q0q2 - two_q1q3) - ae_gy * (two_q0q3 + two_q1q2);
    z[1] = accy_b - ae_gy * (q0q0 - q1q1 + q2q2 - q3q3) - ad_gz * (two_q0q1 + two_q2q3) + an_gx * (two_q0q3 - two_q1q2);
    z[2] = accz_b - ad_gz * (q0q0 - q1q1 - q2q2 + q3q3) + ae_gy * (two_q0q1 - two_q2q3) - an_gx * (two_q0q2 + two_q1q3);
}


/*
 * Arguments    : const float x[28]
 *                float z[3]
 * Return Type  : void
 */
void magMeans(const float *x, float *z)
{
    const float q0 = x[0];
    const float q1 = x[1];
    const float q2 = x[2];
    const float q3 = x[3];
    const float q0q0  = q0*q0;
    const float q1q1  = q1*q1;
    const float q2q2  = q2*q2;
    const float q3q3  = q3*q3;
    const float q1q2 = q1*q2;
    const float q1q3 = q1*q3;
    const float q2q3 = q2*q3;
    const float q0q1 = q0*q1;
    const float q0q2 = q0*q2;
    const float q0q3 = q0*q3;
    const float magNavX = x[22];
    const float magNavY = x[23];
    const float magNavZ = x[24];

    z[0] = x[25] + magNavX*(q0q0 + q1q1 - q2q2 - q3q3) - magNavZ*(2*(q0q2 - q1q3)) + magNavY*(2*(q0q3 + q1q2));
    z[1] = x[26] + magNavY*(q0q0 - q1q1 + q2q2 - q3q3) + magNavZ*(2*(q0q1 + q2q3)) - magNavX*(2*(q0q3 - q1q2));
    z[2] = x[27] + magNavZ*(q0q0 - q1q1 - q2q2 + q3q3) - magNavY*(2*(q0q1 - q2q3)) + magNavX*(2*(q0q2 + q1q3));

}


/*
 * Arguments    : const float x[28]
 *                float dfdx[784]
 * Return Type  : void
 */
void stateTransJacobianFcn(const float *x, float *dfdx)
{
    const float q0   = x[0]*0.5f;
    const float q1   = x[1]*0.5f;
    const float q2   = x[2]*0.5f;
    const float q3   = x[3]*0.5f;
    const float wx   = x[4]*0.5f;
    const float wy   = x[5]*0.5f;
    const float wz   = x[6]*0.5f;

    int i;

    const char b_iv[28] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                           1, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                           0, 0, 0, 0, 0, 0, 0, 0};
    const char b_iv1[28] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                            0, 1, 0, 0, 0, 0, 0, 0, 0, 0,
                            0, 0, 0, 0, 0, 0, 0, 0};
    const char iv2[28] =   {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0,
                            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    const char iv3[28] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,
                          0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    const char iv4[28] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                          1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    const char iv5[28] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                          0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};

    dfdx[0]=0;
    dfdx[1]=-wx;
    dfdx[2]=-wy;
    dfdx[3]=-wz;
    dfdx[4]=-q1;
    dfdx[5]=-q2;
    dfdx[6]=-q3;
    dfdx[7]=0;
    dfdx[8]=0;
    dfdx[9]=0;
    dfdx[10]=0;
    dfdx[11]=0;
    dfdx[12]=0;
    dfdx[13]=0;
    dfdx[14]=0;
    dfdx[15]=0;
    dfdx[16]=0;
    dfdx[17]=0;
    dfdx[18]=0;
    dfdx[19]=0;
    dfdx[20]=0;
    dfdx[21]=0;
    dfdx[22]=0;
    dfdx[23]=0;
    dfdx[24]=0;
    dfdx[25]=0;
    dfdx[26]=0;
    dfdx[27]=0;
    dfdx[28]=wx;
    dfdx[29]=0;
    dfdx[30]= wz;
    dfdx[31]=-wy;
    dfdx[32]= q0;
    dfdx[33]=-q3;
    dfdx[34]= q2;
    dfdx[35]=0;
    dfdx[36]=0;
    dfdx[37]=0;
    dfdx[38]=0;
    dfdx[39]=0;
    dfdx[40]=0;
    dfdx[41]=0;
    dfdx[42]=0;
    dfdx[43]=0;
    dfdx[44]=0;
    dfdx[45]=0;
    dfdx[46]=0;
    dfdx[47]=0;
    dfdx[48]=0;
    dfdx[49]=0;
    dfdx[50]=0;
    dfdx[51]=0;
    dfdx[52]=0;
    dfdx[53]=0;
    dfdx[54]=0;
    dfdx[55]=0;
    dfdx[56]=wy;
    dfdx[57]=-wz;
    dfdx[58]=0;
    dfdx[59]= wx;
    dfdx[60]= q3;
    dfdx[61]= q0;
    dfdx[62]=-q1;
    dfdx[63]=0;
    dfdx[64]=0;
    dfdx[65]=0;
    dfdx[66]=0;
    dfdx[67]=0;
    dfdx[68]=0;
    dfdx[69]=0;
    dfdx[70]=0;
    dfdx[71]=0;
    dfdx[72]=0;
    dfdx[73]=0;
    dfdx[74]=0;
    dfdx[75]=0;
    dfdx[76]=0;
    dfdx[77]=0;
    dfdx[78]=0;
    dfdx[79]=0;
    dfdx[80]=0;
    dfdx[81]=0;
    dfdx[82]=0;
    dfdx[83]=0;
    dfdx[84]=wz;
    dfdx[85]= wy;
    dfdx[86]=-wx;
    dfdx[87]=0;
    dfdx[88]=-q2;
    dfdx[89]= q1;
    dfdx[90]= q0;
    dfdx[91]=0;
    dfdx[92]=0;
    dfdx[93]=0;
    dfdx[94]=0;
    dfdx[95]=0;
    dfdx[96]=0;
    dfdx[97]=0;
    dfdx[98]=0;
    dfdx[99]=0;
    dfdx[100]=0;
    dfdx[101]=0;
    dfdx[102]=0;
    dfdx[103]=0;
    dfdx[104]=0;
    dfdx[105]=0;
    dfdx[106]=0;
    dfdx[107]=0;
    dfdx[108]=0;
    dfdx[109]=0;
    dfdx[110]=0;
    dfdx[111]=0;
    
    for (i = 0; i < 28; i++)
    {
        dfdx[i + 112] = 0.0f;
        dfdx[i + 140] = 0.0f;
        dfdx[i + 168] = 0.0f;
        dfdx[i + 196] = b_iv[i];
        dfdx[i + 224] = b_iv1[i];
        dfdx[i + 252] = iv2[i];
        dfdx[i + 280] = iv3[i];
        dfdx[i + 308] = iv4[i];
        dfdx[i + 336] = iv5[i];
        dfdx[i + 364] = 0.0f;
        dfdx[i + 392] = 0.0f;
        dfdx[i + 420] = 0.0f;
        dfdx[i + 448] = 0.0f;
        dfdx[i + 476] = 0.0f;
        dfdx[i + 504] = 0.0f;
        dfdx[i + 532] = 0.0f;
        dfdx[i + 560] = 0.0f;
        dfdx[i + 588] = 0.0f;
        dfdx[i + 616] = 0.0f;
        dfdx[i + 644] = 0.0f;
        dfdx[i + 672] = 0.0f;
        dfdx[i + 700] = 0.0f;
        dfdx[i + 728] = 0.0f;
        dfdx[i + 756] = 0.0f;
    }
}

/*
 * Arguments    : const float x[28]
 *                float dhdx[84]
 * Return Type  : void
 */
void accelJacobian(const float *x, float *dhdx)
{
    const float q0 = x[0];
    const float q1 = x[1];
    const float q2 = x[2];
    const float q3 = x[3];
     
    // Precompute reused terms
    const float an_g = x[13];
    const float ae_g = x[14];
    const float ad_g = x[15] - Actual_gravity;
    
    const float q0q0 = q0 * q0;
    const float q1q1 = q1 * q1;
    const float q2q2 = q2 * q2;
    const float q3q3 = q3 * q3;
    
    const float q0q1 = q0 * q1;
    const float q0q2 = q0 * q2;
    const float q0q3 = q0 * q3;
    const float q1q2 = q1 * q2;
    const float q1q3 = q1 * q3;
    const float q2q3 = q2 * q3;
//	 2*q2*(ad - gnavz) - 2*q3*(ae - gnavy) - 2*q0*(an - gnavx), - 2*q3*(ad - gnavz) - 2*q2*(ae - gnavy) - 2*q1*(an - gnavx),   2*q0*(ad - gnavz) - 2*q1*(ae - gnavy) + 2*q2*(an - gnavx),   2*q3*(an - gnavx) - 2*q0*(ae - gnavy) - 2*q1*(ad - gnavz), 0, 0, 0, 0, 0, 0, 0, 0, 0, - q0^2 - q1^2 + q2^2 + q3^2,         - 2*q0*q3 - 2*q1*q2,           2*q0*q2 - 2*q1*q3, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
//   2*q3*(an - gnavx) - 2*q0*(ae - gnavy) - 2*q1*(ad - gnavz),   2*q1*(ae - gnavy) - 2*q0*(ad - gnavz) - 2*q2*(an - gnavx), - 2*q3*(ad - gnavz) - 2*q2*(ae - gnavy) - 2*q1*(an - gnavx),   2*q3*(ae - gnavy) - 2*q2*(ad - gnavz) + 2*q0*(an - gnavx), 0, 0, 0, 0, 0, 0, 0, 0, 0,           2*q0*q3 - 2*q1*q2, - q0^2 + q1^2 - q2^2 + q3^2,         - 2*q0*q1 - 2*q2*q3, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
//   2*q1*(ae - gnavy) - 2*q0*(ad - gnavz) - 2*q2*(an - gnavx),   2*q1*(ad - gnavz) + 2*q0*(ae - gnavy) - 2*q3*(an - gnavx),   2*q2*(ad - gnavz) - 2*q3*(ae - gnavy) - 2*q0*(an - gnavx), - 2*q3*(ad - gnavz) - 2*q2*(ae - gnavy) - 2*q1*(an - gnavx), 0, 0, 0, 0, 0, 0, 0, 0, 0,         - 2*q0*q2 - 2*q1*q3,           2*q0*q1 - 2*q2*q3, - q0^2 + q1^2 + q2^2 - q3^2, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0];
// 

    dhdx[0] =  2.0 * q2 * ad_g - 2.0 * q3 * ae_g - 2.0 * q0 * an_g;
    dhdx[1] = -2.0 * q3 * ad_g - 2.0 * q2 * ae_g - 2.0 * q1 * an_g;
    dhdx[2] =  2.0 * q0 * ad_g - 2.0 * q1 * ae_g + 2.0 * q2 * an_g;
    dhdx[3] =  2.0 * q3 * an_g - 2.0 * q0 * ae_g - 2.0 * q1 * ad_g;
    dhdx[4] =  0;
    dhdx[5] =  0;
    dhdx[6] =  0;
    dhdx[7] =  0;
    dhdx[8] =  0;
    dhdx[9] =  0;
    dhdx[10] =  0;
    dhdx[11] =  0;
    dhdx[12] =  0;
    dhdx[13] = -q0q0 - q1q1 + q2q2 + q3q3;
    dhdx[14] = -2.0 * q0q3 - 2.0 * q1q2;
    dhdx[15] =  2.0 * q0q2 - 2.0 * q1q3;
    dhdx[16] = 1.0;
    dhdx[17] =  0;
    dhdx[18] =  0;
    dhdx[19] =  0;
    dhdx[20] =  0;
    dhdx[21] =  0;
    dhdx[22] =  0;
    dhdx[23] =  0;
    dhdx[24] =  0;
    dhdx[25] =  0;
    dhdx[26] =  0;
    dhdx[27] =  0;
    dhdx[28] =  2.0 * q3 * an_g - 2.0 * q0 * ae_g - 2.0 * q1 * ad_g;
    dhdx[29] =  2.0 * q1 * ae_g - 2.0 * q0 * ad_g - 2.0 * q2 * an_g;
    dhdx[30] =  -2.0 * q3 * ad_g - 2.0 * q2 * ae_g - 2.0 * q1 * an_g;
    dhdx[31] =  2.0 * q3 * ae_g - 2.0 * q2 * ad_g + 2.0 * q0 * an_g;
    dhdx[32] =  0;
    dhdx[33] =  0;
    dhdx[34] =  0;
    dhdx[35] =  0;
    dhdx[36] =  0;
    dhdx[37] =  0;
    dhdx[38] =  0;
    dhdx[39] =  0;
    dhdx[40] =  0;
    dhdx[41] =  2.0 * q0q3 - 2.0 * q1q2;
    dhdx[42] =  -q0q0 + q1q1 - q2q2 + q3q3;
    dhdx[43] =  -2.0 * q0q1 - 2.0 * q2q3;
    dhdx[44] =  0;
    dhdx[45] =  1;
    dhdx[46] =  0;
    dhdx[47] =  0;
    dhdx[48] =  0;
    dhdx[49] =  0;
    dhdx[50] =  0;
    dhdx[51] =  0;
    dhdx[52] =  0;
    dhdx[53] =  0;
    dhdx[54] =  0;
    dhdx[55] =  0;
    dhdx[56] =  2.0 * q1 * ae_g - 2.0 * q0 * ad_g - 2.0 * q2 * an_g;
    dhdx[57] =  2.0 * q1 * ad_g + 2.0 * q0 * ae_g - 2.0 * q3 * an_g;
    dhdx[58] =  2.0 * q2 * ad_g - 2.0 * q3 * ae_g - 2.0 * q0 * an_g;
    dhdx[59] =  -2.0 * q3 * ad_g - 2.0 * q2 * ae_g - 2.0 * q1 * an_g;
    dhdx[60] =  0;
    dhdx[61] =  0;
    dhdx[62] =  0;
    dhdx[63] =  0;
    dhdx[64] =  0;
    dhdx[65] =  0;
    dhdx[66] =  0;
    dhdx[67] =  0;
    dhdx[68] =  0;
    dhdx[69] = -2.0 * q0q2 - 2.0 * q1q3;
    dhdx[70] = 2.0 * q0q1 - 2.0 * q2q3;
    dhdx[71] = -q0q0 + q1q1 + q2q2 - q3q3;
    dhdx[72] =  0;
    dhdx[73] =  0;
    dhdx[74] =  1;
    dhdx[75] =  0;
    dhdx[76] =  0;
    dhdx[77] =  0;
    dhdx[78] =  0;
    dhdx[79] =  0;
    dhdx[80] =  0;
    dhdx[81] =  0;
    dhdx[82] =  0;
    dhdx[83] =  0;
}

/*
 * Arguments    : const float x[28]
 *                float dhdx[84]
 * Return Type  : void
 */
void magJakobian(const float *x, float *dhdx)
{  
    const float q0 = x[0];
    const float q1 = x[1];
    const float q2 = x[2];
    const float q3 = x[3];
    const float q0q0   = q0*q0;
    const float q1q1   = q1*q1;
    const float q2q2   = q2*q2;
    const float q3q3   = q3*q3;
    const float q0q3   = 2.0f*q0*q3;
    const float q0q2   = 2.0f*q0*q2;
    const float q0q1   = 2.0f*q0*q1;
    const float q1q3   = 2.0f*q1*q3;
    const float q1q2   = 2.0f*q1*q2;
    const float q2q3   = 2.0f*q2*q3;
    const float magNavYq3 = 2*x[23]*q3;
    const float magNavZq2 = 2*x[24]*q2;
    const float magNavXq0 = 2*x[22]*q0;
    const float magNavZq3 = 2*x[24]*q3;
    const float magNavYq2 = 2*x[23]*q2;
    const float magNavXq1 = 2*x[22]*q1;
    const float magNavYq1 = 2*x[23]*q1;
    const float magNavZq0 = 2*x[24]*q0;
    const float magNavXq2 = 2*x[22]*q2;
    const float magNavZq1 = 2*x[24]*q1;
    const float magNavYq0 = 2*x[23]*q0;
    const float magNavXq3 = 2*x[22]*q3;

    dhdx[0] = magNavYq3 - magNavZq2 + magNavXq0;
    dhdx[1] = magNavZq3 + magNavYq2 + magNavXq1;
    dhdx[2] = magNavYq1 - magNavZq0 - magNavXq2;
    dhdx[3] = magNavZq1 + magNavYq0 - magNavXq3;
    dhdx[4] = 0;
    dhdx[5] = 0;
    dhdx[6] = 0;
    dhdx[7] = 0;
    dhdx[8] = 0;
    dhdx[9] = 0;
    dhdx[10] = 0;
    dhdx[11] = 0;
    dhdx[12] = 0;
    dhdx[13] = 0;
    dhdx[14] = 0;
    dhdx[15] = 0;
    dhdx[16] = 0;
    dhdx[17] = 0;
    dhdx[18] = 0;
    dhdx[19] = 0;
    dhdx[20] = 0;
    dhdx[21] = 0;
    dhdx[22] = q0q0 + q1q1 - q2q2 - q3q3;
    dhdx[23] = q0q3 + q1q2;
    dhdx[24] = q1q3 - q0q2;
    dhdx[25] = 1;
    dhdx[26] = 0;
    dhdx[27] = 0;
    dhdx[28] = magNavZq1 + magNavYq0 - magNavXq3;
    dhdx[29] = magNavZq0 - magNavYq1 + magNavXq2;
    dhdx[30] = magNavZq3 + magNavYq2 + magNavXq1;
    dhdx[31] = magNavZq2 - magNavYq3 - magNavXq0;
    dhdx[32] = 0;
    dhdx[33] = 0;
    dhdx[34] = 0;
    dhdx[35] = 0;
    dhdx[36] = 0;
    dhdx[37] = 0;
    dhdx[38] = 0;
    dhdx[39] = 0;
    dhdx[40] = 0;
    dhdx[41] = 0;
    dhdx[42] = 0;
    dhdx[43] = 0;
    dhdx[44] = 0;
    dhdx[45] = 0;
    dhdx[46] = 0;
    dhdx[47] = 0;
    dhdx[48] = 0;
    dhdx[49] = 0;
    dhdx[50] = q1q2 - q0q3;
    dhdx[51] = q0q0 - q1q1 + q2q2 - q3q3;
    dhdx[52] = q0q1 + q2q3;
    dhdx[53] = 0;
    dhdx[54] = 1;
    dhdx[55] = 0;
    dhdx[56] = magNavZq0 - magNavYq1 + magNavXq2;
    dhdx[57] = magNavXq3 - magNavYq0 - magNavZq1;
    dhdx[58] = magNavYq3 - magNavZq2 + magNavXq0;
    dhdx[59] = magNavZq3 + magNavYq2 + magNavXq1;
    dhdx[60] = 0;
    dhdx[61] = 0;
    dhdx[62] = 0;
    dhdx[63] = 0;
    dhdx[64] = 0;
    dhdx[65] = 0;
    dhdx[66] = 0;
    dhdx[67] = 0;
    dhdx[68] = 0;
    dhdx[69] = 0;
    dhdx[70] = 0;
    dhdx[71] = 0;
    dhdx[72] = 0;
    dhdx[73] = 0;
    dhdx[74] = 0;
    dhdx[75] = 0;
    dhdx[76] = 0;
    dhdx[77] = 0;
    dhdx[78] = q0q2 + q1q3;
    dhdx[79] = q2q3 - q0q1;
    dhdx[80] = q0q0 - q1q1 - q2q2 + q3q3;
    dhdx[81] = 0;
    dhdx[82] = 0;
    dhdx[83] = 1;
}


/*
 * Arguments    : float x[28]
 * Return Type  : void
 */
void stateTransFcn(float *x)
{
    const float q0   = x[0];
    const float q1   = x[1];
    const float q2   = x[2];
    const float q3   = x[3];
    const float wx   = x[4];
    const float wy   = x[5];
    const float wz   = x[6];

    x[0] = - (q1*wx)*0.5f - (q2*wy)*0.5f - (q3*wz)*0.5f;
    x[1] = (q0*wx)*0.5f - (q3*wy)*0.5f + (q2*wz)*0.5f;
    x[2] = (q3*wx)*0.5f + (q0*wy)*0.5f - (q1*wz)*0.5f;
    x[3] = (q1*wy)*0.5f - (q2*wx)*0.5f + (q0*wz)*0.5f;
    x[4] = 0.0f;
    x[5] = 0.0f;
    x[6] = 0.0f;
    x[7] = x[10];
    x[8] = x[11];
    x[9] = x[12];
    x[10] = x[13];
    x[11] = x[14];
    x[12] = x[15];
    x[13] = 0;
    x[14] = 0;
    x[15] = 0;
    x[16] = 0;
    x[17] = 0;
    x[18] = 0;
    x[19] = 0;
    x[20] = 0;
    x[21] = 0;
    x[22] = 0;
    x[23] = 0;
    x[24] = 0;
    x[25] = 0;
    x[26] = 0;
    x[27] = 0;
}

/*
 * Arguments    : c_fusion_internal_coder_AsyncMA *obj
 *                const float accel[3]
 *                float Raccel
 * Return Type  : void
 */
void INS_Filt_fuseaccel(insfilterAsync *obj, const float *accel)
{
    float dv0[3];
    float innov_covariance[9] = {0,0,0,
                                 0,0,0,
                                 0,0,0};
    static float dv1[84] __attribute__((section(".DTCM_RAM")));

                                 
    //means
    accelMeans(obj->pState, dv0);
    //jacobian
    accelJacobian(obj->pState, dv1);
                                 
    const float *R = ekf_cfg_instance()->Racc;                              

    AccelEKF_equationInnovation(obj->pStateCovariance, dv0, dv1, accel, R, obj->accl_innov,innov_covariance);

    if(!IMUBasicEKF_correctEqn(obj->pState, obj->pStateCovariance, dv1, obj->accl_innov, innov_covariance, R))
      reset_ekf();
    get_norm_quat(obj->pState);
}

/*
 * Arguments    : c_fusion_internal_coder_AsyncMA *obj
 *                const float lla[3]
 *                float Rpos
 *                const float vel[3]
 *                float Rvel
 * Return Type  : void
 */
void INS_Filt_fusegps(insfilterAsync *obj, const double *lla, const float *vel)
{
    float dv1[6];
    double dv0[3];
    float h[6];
    float innov_covariance[36] = {0,0,0,0,0,0,
                                  0,0,0,0,0,0,
                                  0,0,0,0,0,0,
                                  0,0,0,0,0,0,
                                  0,0,0,0,0,0,
                                  0,0,0,0,0,0};
    const float dhdx[168] = { 
                0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
                                  
    lla2ned(lla,obj->ReferenceLocation,dv0);

    dv1[0] = (float)dv0[0];
    dv1[1] = (float)dv0[1];
    dv1[2] = (float)dv0[2];
    dv1[3] = vel[0];
    dv1[4] = vel[1];
    dv1[5] = vel[2];

    h[0] = obj->pState[7];
    h[1] = obj->pState[8];
    h[2] = obj->pState[9];
    h[3] = obj->pState[10];
    h[4] = obj->pState[11];
    h[5] = obj->pState[12];

    const float *R = ekf_cfg_instance()->Rgps;                              
                                  
    GpsEKF_equationInnovation(obj->pStateCovariance, h, dv1, R, obj->gps_innov, innov_covariance);

    GPSBasicEKF_correctEqn(obj->pState, obj->pStateCovariance, dhdx, obj->gps_innov, innov_covariance, R);

    get_norm_quat(obj->pState);
}

/*
 * Arguments    : c_fusion_internal_coder_AsyncMA *obj
 *                const float gyro[3]
 *                float Rgyro
 * Return Type  : void
 */
void INS_Filt_fusegyro(insfilterAsync *obj, const float *gyro)
{
    float dv0[3];
    float innov_covariance[9] = {0,0,0,
                                 0,0,0,
                                 0,0,0};

    const float dv1[84] = {
        0.0f, 0.0f, 0.0f, 0.0f, 1.0, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f,
        0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 1.0, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f,
        0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 1.0, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f,
        0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 1.0, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f,
        0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 1.0, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f,
        0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 1.0, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f};

    dv0[0] = obj->pState[19] + obj->pState[4];
    dv0[1] = obj->pState[20] + obj->pState[5];
    dv0[2] = obj->pState[21] + obj->pState[6];

    const float *R = ekf_cfg_instance()->Rgyro;
        
    GyroEKF_equationInnovation(obj->pStateCovariance, dv0, dv1, gyro, R, obj->gyro_innov,innov_covariance);

    if(!IMUBasicEKF_correctEqn(obj->pState, obj->pStateCovariance, dv1, obj->gyro_innov, innov_covariance, R))
      reset_ekf();
    
    get_norm_quat(obj->pState);

}

/*
 * Arguments    : c_fusion_internal_coder_AsyncMA *obj
 *                const float mag[3]
 *                float Rmag
 * Return Type  : void
 */
void INS_Filt_fusemag(insfilterAsync *obj, float *mag)
{
    float dv0[3];
    float innov_covariance[9] = {0,0,0,
                                 0,0,0,
                                 0,0,0};
    static float dv1[84] __attribute__((section(".DTCM_RAM")));

    //Means
    magMeans(obj->pState, dv0);
    //jakobian
    magJakobian(obj->pState, dv1);
  
    mag[2] = dv0[2];                                     
               
    const float *R = ekf_cfg_instance()->Rmag;                              
    MagEKF_equationInnovation(obj->pStateCovariance, dv0, dv1, mag, R, obj->mag_innov, innov_covariance);

    if(!IMUBasicEKF_correctEqn(obj->pState, obj->pStateCovariance, dv1, obj->mag_innov, innov_covariance, R))
      reset_ekf();

    get_norm_quat(obj->pState);
}

void INS_Filt_fusebaro_height(insfilterAsync *obj, float posD_baro, float R_baro)
{
    const int n = 28;
    const int idx = 9; // posD

    float innov = posD_baro - obj->pState[idx];
    float S = obj->pStateCovariance[idx + n*idx] + R_baro;

    if (!isfinite(S) || S < 1e-6f) {
        return;
    }

    float K[28];

    for (int i = 0; i < n; i++) {
        K[i] = obj->pStateCovariance[i + n*idx] / S;
    }

    for (int i = 0; i < n; i++) {
        obj->pState[i] += K[i] * innov;
    }

    // Joseph scalar update:
    // P_new = P - K*H*P - P*H'*K' + K*S*K'
    // H selects state idx only.

    float Pold[28 * 28];
    //memcpy(Pold, obj->pStateCovariance, sizeof(Pold));
    mat_copy28(Pold, obj->pStateCovariance);

    for (int r = 0; r < n; r++) {
        for (int c = 0; c < n; c++) {
            obj->pStateCovariance[r + n*c] =
                Pold[r + n*c]
              - K[r] * Pold[idx + n*c]
              - Pold[r + n*idx] * K[c]
              + K[r] * S * K[c];
        }
    }

    ekf_fix_covariance(obj->pStateCovariance, 28);
}

void ContinuousEKFPredictor_x(const float *x, const float *xdot, float *e, float dt)
{
    e[0]  = x[0]  + dt * xdot[0] ;
    e[1]  = x[1]  + dt * xdot[1] ;
    e[2]  = x[2]  + dt * xdot[2] ;
    e[3]  = x[3]  + dt * xdot[3] ;
    e[4]  = x[4]  + dt * xdot[4] ;
    e[5]  = x[5]  + dt * xdot[5] ;
    e[6]  = x[6]  + dt * xdot[6] ;
    e[7]  = x[7]  + dt * xdot[7] ;
    e[8]  = x[8]  + dt * xdot[8] ;
    e[9]  = x[9]  + dt * xdot[9] ;
    e[10] = x[10] + dt * xdot[10];
    e[11] = x[11] + dt * xdot[11];
    e[12] = x[12] + dt * xdot[12];
    e[13] = x[13] + dt * xdot[13];
    e[14] = x[14] + dt * xdot[14];
    e[15] = x[15] + dt * xdot[15];
    e[16] = x[16] + dt * xdot[16];
    e[17] = x[17] + dt * xdot[17];
    e[18] = x[18] + dt * xdot[18];
    e[19] = x[19] + dt * xdot[19];
    e[20] = x[20] + dt * xdot[20];
    e[21] = x[21] + dt * xdot[21];
    e[22] = x[22] + dt * xdot[22];
    e[23] = x[23] + dt * xdot[23];
    e[24] = x[24] + dt * xdot[24];
    e[25] = x[25] + dt * xdot[25];
    e[26] = x[26] + dt * xdot[26];
    e[27] = x[27] + dt * xdot[27];
}


void ContinuousEKFPredictor_p(const float *P,const float *dfdx, float *pStateCovariance, float dt)
{
    int i0,i1;
    for (i0 = 0; i0 < 28; i0++)
    {
        i1 = 0+i0;
        pStateCovariance[i1] = P[i1] + dt * dfdx[i1];
        i1 = 28+i0;
        pStateCovariance[i1] = P[i1] + dt * dfdx[i1];
        i1 = 56+i0;
        pStateCovariance[i1] = P[i1] + dt * dfdx[i1];
        i1 = 84+i0;
        pStateCovariance[i1] = P[i1] + dt * dfdx[i1];
        i1 = 112+i0;
        pStateCovariance[i1] = P[i1] + dt * dfdx[i1];
        i1 = 140+i0;
        pStateCovariance[i1] = P[i1] + dt * dfdx[i1];
        i1 = 168+i0;
        pStateCovariance[i1] = P[i1] + dt * dfdx[i1];
        i1 = 196+i0;
        pStateCovariance[i1] = P[i1] + dt * dfdx[i1];
        i1 = 224+i0;
        pStateCovariance[i1] = P[i1] + dt * dfdx[i1];
        i1 = 252+i0;
        pStateCovariance[i1] = P[i1] + dt * dfdx[i1];
        i1 = 280+i0;
        pStateCovariance[i1] = P[i1] + dt * dfdx[i1];
        i1 = 308+i0;
        pStateCovariance[i1] = P[i1] + dt * dfdx[i1];
        i1 = 336+i0;
        pStateCovariance[i1] = P[i1] + dt * dfdx[i1];
        i1 = 364+i0;
        pStateCovariance[i1] = P[i1] + dt * dfdx[i1];
        i1 = 392+i0;
        pStateCovariance[i1] = P[i1] + dt * dfdx[i1];
        i1 = 420+i0;
        pStateCovariance[i1] = P[i1] + dt * dfdx[i1];
        i1 = 448+i0;
        pStateCovariance[i1] = P[i1] + dt * dfdx[i1];
        i1 = 476+i0;
        pStateCovariance[i1] = P[i1] + dt * dfdx[i1];
        i1 = 504+i0;
        pStateCovariance[i1] = P[i1] + dt * dfdx[i1];
        i1 = 532+i0;
        pStateCovariance[i1] = P[i1] + dt * dfdx[i1];
        i1 = 560+i0;
        pStateCovariance[i1] = P[i1] + dt * dfdx[i1];
        i1 = 588+i0;
        pStateCovariance[i1] = P[i1] + dt * dfdx[i1];
        i1 = 616+i0;
        pStateCovariance[i1] = P[i1] + dt * dfdx[i1];
        i1 = 644+i0;
        pStateCovariance[i1] = P[i1] + dt * dfdx[i1];
        i1 = 672+i0;
        pStateCovariance[i1] = P[i1] + dt * dfdx[i1];
        i1 = 700+i0;
        pStateCovariance[i1] = P[i1] + dt * dfdx[i1];
        i1 = 728+i0;
        pStateCovariance[i1] = P[i1] + dt * dfdx[i1];
        i1 = 756+i0;
        pStateCovariance[i1] = P[i1] + dt * dfdx[i1];
    }
}

