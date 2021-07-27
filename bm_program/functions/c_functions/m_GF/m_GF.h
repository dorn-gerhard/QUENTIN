/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * m_GF.h
 *
 * Code generation for function 'm_GF'
 *
 */

#ifndef M_GF_H
#define M_GF_H

/* Include files */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "tmwtypes.h"
#include "mex.h"
#include "emlrt.h"
#include "rtwtypes.h"
#include "m_GF_types.h"

/* Function Declarations */
extern void m_GF(const emlrtStack *sp, const emxArray_real_T *points, real_T i0p,
                 const emxArray_real_T *E_p, const emxArray_real_T *E_n, const
                 emxArray_real_T *E_q, const emxArray_creal_T *C1, const
                 emxArray_creal_T *C2, const emxArray_creal_T *C3, const
                 emxArray_creal_T *C4, const emxArray_creal_T *rho,
                 emxArray_creal_T *GR, emxArray_creal_T *GK);

#endif

/* End of code generation (m_GF.h) */
