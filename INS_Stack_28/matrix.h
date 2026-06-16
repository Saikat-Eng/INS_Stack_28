#ifndef _MATRIX_h
#define _MATRIX_h


extern void mat_reset(float *d, int size);
extern void mat_fill(float *m, int size);
extern void mat_copy28(float *dst, const float *src);
extern void mat_cpy9(float *A, float *B);
extern void mat_mul3_9(const float *A,const float *B, float *C);
extern void mat_mul9_3(const float *A,const float *B, float *C);
extern void mat_mul(const float *a, const float *b, float *c);
extern void mat3_transpose(float *R);
extern char invert_3x3(const float *A, float *Ainv);
extern void skewm(const float *v, float *S);
extern void skewm_inv(const float *S, float *v);

#endif
