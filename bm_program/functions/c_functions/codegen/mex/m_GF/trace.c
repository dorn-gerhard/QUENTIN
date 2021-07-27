/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * trace.c
 *
 * Code generation for function 'trace'
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "m_GF.h"
#include "trace.h"

/* Variable Definitions */
static emlrtRTEInfo ab_emlrtRTEI = { 11,/* lineNo */
  15,                                  /* colNo */
  "trace",                             /* fName */
  "/afs/itp.tugraz.at/opt/matlab/R2018b/toolbox/eml/lib/matlab/matfun/trace.m"/* pName */
};

/* Function Definitions */
creal_T trace(const emlrtStack *sp, const emxArray_creal_T *a)
{
  creal_T t;
  boolean_T b0;
  int32_T i5;
  int32_T k;
  b0 = (a->size[0] == a->size[1]);
  if (!b0) {
    emlrtErrorWithMessageIdR2018a(sp, &ab_emlrtRTEI, "Coder:MATLAB:square",
      "Coder:MATLAB:square", 0);
  }

  t.re = 0.0;
  t.im = 0.0;
  i5 = a->size[0];
  for (k = 0; k < i5; k++) {
    t.re += a->data[k + a->size[0] * k].re;
    t.im += a->data[k + a->size[0] * k].im;
  }

  return t;
}

/* End of code generation (trace.c) */
