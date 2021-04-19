#ifndef BSR_UTIL_H
#define BSR_UTIL_H

int waitForWorkerThreads(bsr_state_t *bsr_state, int min_status);
int waitForMainThread(bsr_state_t *bsr_state, int min_status);
int checkExceptions(bsr_state_t *bsr_state);
int limitIntensity(double *pixel_r, double *pixel_g, double *pixel_b);
int limitIntensityPreserveColor(double *pixel_r, double *pixel_g, double *pixel_b);

#endif // BSR_UTIL_H
