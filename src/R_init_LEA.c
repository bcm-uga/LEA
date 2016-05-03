#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

#include "R_ancestrymap2geno.h"
#include "R_ancestrymap2lfmm.h"
#include "R_lfmm2geno.h"
#include "R_geno2lfmm.h"
#include "R_ped2geno.h"
#include "R_ped2lfmm.h"
#include "R_vcf2geno.h"
#include "R_createDataSet.h"
#include "R_crossEntropy.h"
#include "R_LFMM.h"
#include "R_pca.h"
#include "R_sNMF.h"
#include "R_tracyWidom.h"

// register the routines
static R_CMethodDef cMethods[] = {
        {"R_ancestrymap2geno", (DL_FUNC) &R_ancestrymap2geno, 4, 
                (R_NativePrimitiveArgType[4]) {STRSXP, STRSXP,
                                                                   INTSXP,
                                                                   INTSXP}},
        {"R_ancestrymap2lfmm", (DL_FUNC) &R_ancestrymap2lfmm, 4, 
                (R_NativePrimitiveArgType[4]) {STRSXP,
                                                                   STRSXP,
                                                                   INTSXP,
                                                                   INTSXP}},
        {"R_lfmm2geno", (DL_FUNC) &R_lfmm2geno, 4, 
                (R_NativePrimitiveArgType[4]) {STRSXP,
                                                     STRSXP, INTSXP, INTSXP}},
        {"R_geno2lfmm", (DL_FUNC) &R_geno2lfmm, 4, 
                (R_NativePrimitiveArgType[4]) {STRSXP,
                                                     STRSXP, INTSXP, INTSXP}},
        {"R_vcf2geno", (DL_FUNC) &R_vcf2geno, 4, 
                (R_NativePrimitiveArgType[4]) {STRSXP,
                                                   STRSXP, INTSXP, INTSXP}},
        {"R_ped2geno", (DL_FUNC) &R_ped2geno, 4, 
                (R_NativePrimitiveArgType[4]) {STRSXP,
                                                   STRSXP, INTSXP, INTSXP}},
        {"R_ped2lfmm", (DL_FUNC) &R_ped2lfmm, 4, 
                (R_NativePrimitiveArgType[4]) {STRSXP,
                                                   STRSXP, INTSXP, INTSXP}},
        {"R_createDataSet", (DL_FUNC) &R_createDataSet, 4, 
                (R_NativePrimitiveArgType[4]) {STRSXP,
                                                             INTSXP, REALSXP,
                                                             STRSXP}}, 
        {"R_crossEntropy", (DL_FUNC) &R_crossEntropy, 8, 
                (R_NativePrimitiveArgType[8]) {STRSXP, STRSXP,
                                                        STRSXP, STRSXP,
                                                           INTSXP, INTSXP,
                                                           REALSXP, REALSXP}},
        {"R_LFMM", (DL_FUNC) &R_LFMM, 19, 
                (R_NativePrimitiveArgType[19]) {STRSXP, STRSXP, STRSXP, INTSXP,
                                            INTSXP, INTSXP, INTSXP, INTSXP,
                                            INTSXP, INTSXP, INTSXP, INTSXP,
                                            INTSXP,
                                            INTSXP, REALSXP, REALSXP, REALSXP,
                                            REALSXP, LGLSXP}},
        {"R_pca", (DL_FUNC) &R_pca, 10, (R_NativePrimitiveArgType[10]) 
                                            {STRSXP, STRSXP, STRSXP, STRSXP,
                                            STRSXP, INTSXP, INTSXP, INTSXP,
                                            INTSXP, INTSXP}},
        {"R_sNMF", (DL_FUNC) &R_sNMF, 17, (R_NativePrimitiveArgType[17])
                                            {STRSXP, INTSXP, REALSXP, REALSXP,
                                            REALSXP, INTSXP, INTSXP, INTSXP,
                                            INTSXP, STRSXP, STRSXP, STRSXP,
                                            INTSXP, REALSXP, REALSXP, INTSXP,
                                            INTSXP}}, 
        {"R_tracyWidom", (DL_FUNC) &R_tracyWidom, 2, 
                (R_NativePrimitiveArgType[4]){STRSXP, STRSXP}},
        {NULL, NULL, 0}
};

void R_init_LEA(DllInfo * info)
{
           R_registerRoutines(info, cMethods, NULL, NULL, NULL);
}
