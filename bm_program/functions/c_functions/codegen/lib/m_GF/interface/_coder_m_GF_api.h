/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * _coder_m_GF_api.h
 *
 * Code generation for function '_coder_m_GF_api'
 *
 */

#ifndef _CODER_M_GF_API_H
#define _CODER_M_GF_API_H

/* Include files */
#include "tmwtypes.h"
#include "mex.h"
#include "emlrt.h"
#include <stddef.h>
#include <stdlib.h>
#include "_coder_m_GF_api.h"

/* Type Definitions */
#ifndef typedef_emxArray_creal_T
#define typedef_emxArray_creal_T

typedef struct {
  creal_T *data;
  int32_T *size;
  int32_T allocatedSize;
  int32_T numDimensions;
  boolean_T canFreeData;
} emxArray_creal_T;

#endif                                 /*typedef_emxArray_creal_T*/

#ifndef struct_emxArray_real_T
#define struct_emxArray_real_T

struct emxArray_real_T
{
  real_T *data;
  int32_T *size;
  int32_T allocatedSize;
  int32_T numDimensions;
  boolean_T canFreeData;
};

#endif                                 /*struct_emxArray_real_T*/

#ifndef typedef_emxArray_real_T
#define typedef_emxArray_real_T

typedef struct emxArray_real_T emxArray_real_T;

#endif                                 /*typedef_emxArray_real_T*/

/* Variable Declarations */
extern emlrtCTX emlrtRootTLSGlobal;
extern emlrtContext emlrtContextGlobal;

/* Function Declarations */
extern void m_GF(emxArray_real_T *points, real_T i0p, emxArray_real_T *E_p,
                 emxArray_real_T *E_n, emxArray_real_T *E_q, emxArray_creal_T
                 *C1, emxArray_creal_T *C2, emxArray_creal_T *C3,
                 emxArray_creal_T *C4, emxArray_creal_T *rho, emxArray_creal_T
                 *GR, emxArray_creal_T *GK);
extern void m_GF_api(const mxArray * const prhs[10], int32_T nlhs, const mxArray
                     *plhs[2]);
extern void m_GF_atexit(void);
extern void m_GF_initialize(void);
extern void m_GF_terminate(void);
extern void m_GF_xil_terminate(void);

#endif

/* End of code generation (_coder_m_GF_api.h) */
