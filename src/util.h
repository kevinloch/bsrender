#ifndef BSR_UTIL_H
#define BSR_UTIL_H

int waitForWorkerThreads(int *status_array, int num_worker_threads);
int limitIntensity(double *pixel_r, double *pixel_g, double *pixel_b);
int limitIntensityPreserveColor(double *pixel_r, double *pixel_g, double *pixel_b);

#endif // BSR_UTIL_H
