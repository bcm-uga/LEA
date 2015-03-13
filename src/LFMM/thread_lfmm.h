/**
 * @addtogroup thread_lfmm
 * @ingroup LFMM
 * @{
 * @file thread_lfmm.h
 *
 * @brief general function and structure to manage multithreading
 */

#ifndef THREAD_LFMM_H
#define THREAD_LFMM_H

#ifndef WIN32

typedef struct _multithreading_lfmm *Multithreading_lfmm;

/**
 * @brief structure containing generic parameters for multithreading
 * 	in LFMM functions.
 */
typedef struct _multithreading_lfmm {
        float *R;
        double *A;
        double *B;
        double *C;
        double *m;
        double *inv_cov;
        double *L;
        int J;
        int N;
        int M;
        int K;
        int mode;
        double *alpha;
        double alpha_R;
        int slice;
        int c;
        int num_thrd;
} multithreading_lfmm;

/**
 * general multithreading function manager. Some parameters can be NULL
 *
 */
void thread_fct_lfmm(float *R, double *A, double *B, double *C, double *m,
                     double *inv_cov, double *L, int J, int K, int N, int M,
                     double *alpha, double alpha_R, int num_thrd, int mode,
                     void (*fct) ());

#endif

#endif                          // THREAD_H

/** @} */
