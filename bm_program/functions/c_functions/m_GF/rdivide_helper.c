/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * rdivide_helper.c
 *
 * Code generation for function 'rdivide_helper'
 *
 */

/* Include files */
#include "mwmathutil.h"
#include "rt_nonfinite.h"
#include "m_GF.h"
#include "rdivide_helper.h"
#include "m_GF_emxutil.h"

/* Variable Definitions */
static emlrtRTEInfo t_emlrtRTEI = { 28,/* lineNo */
  1,                                   /* colNo */
  "rdivide_helper",                    /* fName */
  "/afs/itp.tugraz.at/opt/matlab/R2018b/toolbox/eml/eml/+coder/+internal/rdivide_helper.m"/* pName */
};

/* Function Definitions */
void b_rdivide_helper(const emlrtStack *sp, const emxArray_real_T *y,
                      emxArray_real_T *z)
{
  int32_T i6;
  int32_T loop_ub;
  i6 = z->size[0] * z->size[1];
  z->size[0] = y->size[0];
  z->size[1] = y->size[1];
  emxEnsureCapacity_real_T(sp, z, i6, &t_emlrtRTEI);
  loop_ub = y->size[0] * y->size[1];
  for (i6 = 0; i6 < loop_ub; i6++) {
    z->data[i6] = 1.0 / y->data[i6];
  }
}

void rdivide_helper(const emlrtStack *sp, const emxArray_creal_T *y,
                    emxArray_creal_T *z)
{
  int32_T i4;
  int32_T loop_ub;
  real_T y_re;
  real_T y_im;
  real_T brm;
  real_T bim;
  real_T s;
  i4 = z->size[0] * z->size[1];
  z->size[0] = y->size[0];
  z->size[1] = y->size[1];
  emxEnsureCapacity_creal_T(sp, z, i4, &t_emlrtRTEI);
  loop_ub = y->size[0] * y->size[1];
  for (i4 = 0; i4 < loop_ub; i4++) {
    y_re = y->data[i4].re;
    y_im = y->data[i4].im;
    if (y_im == 0.0) {
      z->data[i4].re = 1.0 / y_re;
      z->data[i4].im = 0.0;
    } else if (y_re == 0.0) {
      z->data[i4].re = 0.0;
      z->data[i4].im = -(1.0 / y_im);
    } else {
      brm = muDoubleScalarAbs(y_re);
      bim = muDoubleScalarAbs(y_im);
      if (brm > bim) {
        s = y_im / y_re;
        bim = y_re + s * y_im;
        z->data[i4].re = (1.0 + s * 0.0) / bim;
        z->data[i4].im = (0.0 - s) / bim;
      } else if (bim == brm) {
        if (y_re > 0.0) {
          s = 0.5;
        } else {
          s = -0.5;
        }

        if (y_im > 0.0) {
          bim = 0.5;
        } else {
          bim = -0.5;
        }

        z->data[i4].re = (s + 0.0 * bim) / brm;
        z->data[i4].im = (0.0 * s - bim) / brm;
      } else {
        s = y_re / y_im;
        bim = y_im + s * y_re;
        z->data[i4].re = s / bim;
        z->data[i4].im = (s * 0.0 - 1.0) / bim;
      }
    }
  }
}

/* End of code generation (rdivide_helper.c) */
