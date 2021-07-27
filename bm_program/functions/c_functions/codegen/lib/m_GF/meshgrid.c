/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * meshgrid.c
 *
 * Code generation for function 'meshgrid'
 *
 */

/* Include files */
#include "m_GF.h"
#include "meshgrid.h"
#include "m_GF_emxutil.h"

/* Function Definitions */
void meshgrid(const emxArray_real_T *x, const emxArray_real_T *y,
              emxArray_real_T *xx, emxArray_real_T *yy)
{
  int nx;
  int ny;
  int j;
  int i;
  nx = x->size[0] * x->size[1];
  ny = y->size[0] * y->size[1];
  j = xx->size[0] * xx->size[1];
  xx->size[0] = ny;
  xx->size[1] = nx;
  emxEnsureCapacity_real_T(xx, j);
  j = yy->size[0] * yy->size[1];
  yy->size[0] = ny;
  yy->size[1] = nx;
  emxEnsureCapacity_real_T(yy, j);
  if ((nx == 0) || (ny == 0)) {
  } else {
    for (j = 0; j < nx; j++) {
      for (i = 0; i < ny; i++) {
        xx->data[i + xx->size[0] * j] = x->data[j];
        yy->data[i + yy->size[0] * j] = y->data[i];
      }
    }
  }
}

/* End of code generation (meshgrid.c) */
