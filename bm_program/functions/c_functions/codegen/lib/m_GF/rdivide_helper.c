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
#include <math.h>
#include "m_GF.h"
#include "rdivide_helper.h"
#include "m_GF_emxutil.h"

/* Function Definitions */
void b_rdivide_helper(const emxArray_real_T *y, emxArray_real_T *z)
{
  int i3;
  int loop_ub;
  i3 = z->size[0] * z->size[1];
  z->size[0] = y->size[0];
  z->size[1] = y->size[1];
  emxEnsureCapacity_real_T(z, i3);
  loop_ub = y->size[0] * y->size[1];
  for (i3 = 0; i3 < loop_ub; i3++) {
    z->data[i3] = 1.0 / y->data[i3];
  }
}

void rdivide_helper(const emxArray_creal_T *y, emxArray_creal_T *z)
{
  int i2;
  int loop_ub;
  double y_re;
  double y_im;
  double brm;
  double bim;
  double s;
  i2 = z->size[0] * z->size[1];
  z->size[0] = y->size[0];
  z->size[1] = y->size[1];
  emxEnsureCapacity_creal_T(z, i2);
  loop_ub = y->size[0] * y->size[1];
  for (i2 = 0; i2 < loop_ub; i2++) {
    y_re = y->data[i2].re;
    y_im = y->data[i2].im;
    if (y_im == 0.0) {
      z->data[i2].re = 1.0 / y_re;
      z->data[i2].im = 0.0;
    } else if (y_re == 0.0) {
      z->data[i2].re = 0.0;
      z->data[i2].im = -(1.0 / y_im);
    } else {
      brm = fabs(y_re);
      bim = fabs(y_im);
      if (brm > bim) {
        s = y_im / y_re;
        bim = y_re + s * y_im;
        z->data[i2].re = (1.0 + s * 0.0) / bim;
        z->data[i2].im = (0.0 - s) / bim;
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

        z->data[i2].re = (s + 0.0 * bim) / brm;
        z->data[i2].im = (0.0 * s - bim) / brm;
      } else {
        s = y_re / y_im;
        bim = y_im + s * y_re;
        z->data[i2].re = s / bim;
        z->data[i2].im = (s * 0.0 - 1.0) / bim;
      }
    }
  }
}

/* End of code generation (rdivide_helper.c) */
