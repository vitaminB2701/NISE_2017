#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
#ifdef _WIN32
    #include <Windows.h>
#else
    #include <unistd.h>
#endif
#include "omp.h"
#include "types.h"
#include "NISE_subs.h"
#include "polar.h"
#include "calc_CG_2DES.h"
#include "propagate.h"
#include <stdarg.h>
#include "project.h"

/* This is the main routine which will control the flow of the
 * coarse grained 2DES calculation. It takes the general NISE
 * input as input.
 */

void calc_CG_2DES(t_non *non){
    printf("Hello world!\n");
    float *P_DA;
    float K[] = {-0.1,0.1,0.01,-0.011};
    float P0[] = {1,0,0,1};
    int pro_dim=2;
    FILE *outone;

    P_DA=(float *)calloc(non->tmax2*sizeof(K)/sizeof(K[0]),sizeof(float));
    CG_2DES_P_DA(non,P_DA,K);
    // Write to file
    outone=fopen("KPop.dat","w");
    for (int t1=0;t1<non->tmax2;t1+=non->dt1){
        fprintf(outone,"%f ",t1*non->deltat);
        for (int a=0;a<pro_dim;a++){
            for (int b=0;b<pro_dim;b++){
                fprintf(outone,"%f ",P_DA[t1+(non->singles*b+a)*non->tmax2]);
            }
        }
        fprintf(outone,"\n"); 
    }
    fclose(outone);
    free(P_DA);
    return;
};

void CG_2DES_doorway(t_non *non,float *re_doorway,float *im_doorway){
    printf("Calculate doorway part!\n");
    return;
};

void CG_2DES_P_DA(t_non *non,float *P_DA,float K[]){
    printf("Calculate population transfer!\n");
    // int index, N;
    float *eigK_re, *eigK_im; // eigenvalues of K
    float *evecL, *evecR; // eigenvectors of K
    float *ivecR, *ivecL; //inverse eigenvectors of K
    float *EKt; // expm(K*dt)
    float *cnr;
    float *crr;
    // float re, im;
    int a, b, c;
    int N = 2;
    
    eigK_re = (float *)calloc(N*N,sizeof(float));
    eigK_im = (float *)calloc(N*N,sizeof(float));
    evecL = (float *)calloc(N*N,sizeof(float));
    evecR = (float *)calloc(N*N,sizeof(float));
    ivecL = (float *)calloc(N*N,sizeof(float));
    ivecR = (float *)calloc(N*N,sizeof(float));
    EKt = (float *)calloc(N,sizeof(float));
    cnr = (float *)calloc(N * N, sizeof(float));
    crr = (float *)calloc(N * N, sizeof(float));
    diagonalize_real_nonsym(K, eigK_re, eigK_im, evecL, evecR, ivecL, ivecR, N);
    // printf("%f %f %f %f\n",evecR[0],evecR[1],evecR[2],evecR[3]);
    // printf("%f %f\n",eigK_re[0],eigK_re[1]);
    // printf("%f %f %f %f\n",ivecR[0],ivecR[1],ivecR[2],ivecR[3]);
    for (int a = 0; a<N; a++) {
        if (eigK_im[a]!=0) {
            printf("Transfer lifetime is not real!\n");
            exit(0);
        }
    }

    // P(t2) = expm(-K*t2) = V*exp(EK*t2)/V

        for (a = 0; a < N; a++) {
            EKt[a] = exp(eigK_re[a] * non->deltat);
        }

        /* Multiply with eigenvector */
        for (a = 0; a < N; a++) {
            for (b = 0; b < N; b++) {
                cnr[b + a * N] += EKt[b] * Kt[b];
            }
        }
        for (a = 0; a < N; a++) {
            for (b = 0; b < N; b++) {
                for (c = 0; c < N; c++) {
                    P_DA[nt2+(a + c * N)*non->tmax2] += evecR[b + a * N] * cnr[b + c * N];
                }
            }
        }
        /* The one exciton propagator has been calculated */

        // for (a = 0; a < N; a++) {
        //     cnr[a] = 0;
        //     for (b = 0; b < N; b++) {
        //         cnr[a] += crr[a + b * N] * cr[b];
        //     }
        // }

        // for (a = 0; a < N; a++) {
        //     cr[a] = cnr[a];
        // }
        // Loop over t2
    for (int nt2 = 0; nt2<non->tmax2; nt2++) {}

    free(cnr);
    free(crr);
    free(eigK_im);
    free(eigK_re);
    free(Kt);
    free(evecL);
    free(evecR);
    free(ivecL);
    free(ivecR);

    return;
};
void CG_2DES_window_GB(t_non *non,float *re_window_GB,float *im_window_GB);
void CG_2DES_window_SE(t_non *non,float *re_window_SE,float *im_window_SE);
void CG_2DES_window_EA(t_non *non,float *re_window_EA,float *im_window_EA);
void CG_full_2DES_segments(t_non *non,float *re_2DES,float *im_2DES);
void combine_CG_2DES(t_non *non,float *re_doorway,float *im_doorway,
    float *P_DA,float *re_window_GB,float *im_window_GB,
    float *re_window_SE,float *im_window_SE,float *re_window_EA,float *im_window_EA,
    float *re_2DES,float *im_2DES);

void diagonalize_real_nonsym(float* K, float* eig_re, float* eig_im, float* evecL, float* evecR, float* ivecL, float* ivecR, int N) {
    int INFO, lwork;
    float *work, *Kcopy;
    int i, j;
    float *pivot;
    /* Find lwork for diagonalization */
    lwork = -1;
    work = (float *)calloc(1, sizeof(float));
    sgeev_("V", "V", &N, Kcopy, &N, eig_re, eig_im, evecL, &N, evecR, &N, work, &lwork, &INFO);
    lwork = work[0];
    free(work);
    work = (float *)calloc(lwork, sizeof(float));
    Kcopy = (float *)calloc(N * N, sizeof(float));
    /* Copy matrix */
    for (i = 0; i < N; i++) {
        for (j = 0; j < N; j++) {
            Kcopy[i * N + j] = K[i * N + j];
        }
    }
    /* Do diagonalization*/
    sgeev_("V", "V", &N, Kcopy, &N, eig_re, eig_im, evecL, &N, evecR, &N, work, &lwork, &INFO);
    if (INFO != 0) {
        printf("Something went wrong trying to diagonalize a matrix...\nExit code %d\n",INFO);
        exit(0);
    }
    free(work);
    printf("eigenvalue completed %p\n",eig_re);

    /* Copy matrix */
    for (i = 0; i < N; i++) {
        for (j = 0; j < N; j++) {
            ivecL[i * N + j] = evecL[i * N + j];
            ivecR[i * N + j] = evecR[i * N + j];
        }
    }

    /* Inverse right eigenvectors*/
    pivot = (float *)calloc(N,sizeof(float));
    sgetrf_(&N, &N, ivecR, &N, pivot, &INFO);
    if (INFO != 0) {
        printf("Something went wrong trying to factorize right eigenvector matrix...\nExit code %d\n",INFO);
        exit(0);
    }    
    lwork = -1; /* Find lwork for diagonalization */
    work = (float *)calloc(1, sizeof(float));
    sgetri_(&N, ivecR, &N, pivot, work, &lwork, &INFO);
    lwork = work[0];
    free(work);
    work = (float *)calloc(lwork, sizeof(float));
    sgetri_(&N, ivecR, &N, pivot, work, &lwork, &INFO);
    if (INFO != 0) {
        printf("Something went wrong trying to inverse right eigenvector matrix...\nExit code %d\n",INFO);
        exit(0);
    }
    free(work), free(pivot);

        /* Inverse left eigenvectors*/
    pivot = (float *)calloc(N,sizeof(float));
    sgetrf_(&N, &N, ivecL, &N, pivot, &INFO);
    if (INFO != 0) {
        printf("Something went wrong trying to factorize left eigenvector matrix...\nExit code %d\n",INFO);
        exit(0);
    }    
    lwork = -1; /* Find lwork for diagonalization */
    work = (float *)calloc(1, sizeof(float));
    sgetri_(&N, ivecL, &N, pivot, work, &lwork, &INFO);
    lwork = work[0];
    free(work);
    work = (float *)calloc(lwork, sizeof(float));
    sgetri_(&N, ivecL, &N, pivot, work, &lwork, &INFO);
    if (INFO != 0) {
        printf("Something went wrong trying to inverse left eigenvector matrix...\nExit code %d\n",INFO);
        exit(0);
    }

    /* Free space */
    free(Kcopy), free(work), free(pivot);
    return;
}
