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
    float K[] = {0.1,0,-0.1,0.01};
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
    float *eigK, *Kt;
    float *cnr;
    float *crr;
    // float re, im;
    int a, b, c;
    int N = 2;
    
    cnr = (float *)calloc(N * N, sizeof(float));
    crr = (float *)calloc(N * N, sizeof(float));
    
    diagonalizeLPD(&K[0], eigK, N);
    printf("%d\n",N);
    // P(t2) = exp(-K*t2)
    // Loop over t2
    for (int nt2 = 0; nt2<non->tmax2; nt2++) {
    for (a = 0; a < N; a++) {
        Kt[a] = exp(-eigK[a] * nt2 * non->deltat);
    }

    /* Transform to site basis */
    for (a = 0; a < N; a++) {
        for (b = 0; b < N; b++) {
            cnr[b + a * N] += K[b + a * N] * Kt[b];
        }
    }
    for (a = 0; a < N; a++) {
        for (b = 0; b < N; b++) {
            for (c = 0; c < N; c++) {
                P_DA[nt2+(a + c * N)*non->tmax2] += K[b + a * N] * cnr[b + c * N];
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
    }

    free(cnr);
    free(crr);

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
