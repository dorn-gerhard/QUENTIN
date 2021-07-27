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
#include "rt_nonfinite.h"
#include "m_GF.h"
#include "meshgrid.h"
#include "m_GF_emxutil.h"
#include "eml_int_forloop_overflow_check.h"
#include "m_GF_data.h"

/* Variable Definitions */
static emlrtRSInfo e_emlrtRSI = { 31,  /* lineNo */
  "meshgrid",                          /* fcnName */
  "/afs/itp.tugraz.at/opt/matlab/R2018b/toolbox/eml/lib/matlab/elmat/meshgrid.m"/* pathName */
};

static emlrtRSInfo f_emlrtRSI = { 32,  /* lineNo */
  "meshgrid",                          /* fcnName */
  "/afs/itp.tugraz.at/opt/matlab/R2018b/toolbox/eml/lib/matlab/elmat/meshgrid.m"/* pathName */
};

static emlrtRTEInfo s_emlrtRTEI = { 1, /* lineNo */
  23,                                  /* colNo */
  "meshgrid",                          /* fName */
  "/afs/itp.tugraz.at/opt/matlab/R2018b/toolbox/eml/lib/matlab/elmat/meshgrid.m"/* pName */
};

/* Function Definitions */
void meshgrid(const emlrtStack *sp, const emxArray_real_T *x, const
              emxArray_real_T *y, emxArray_real_T *xx, emxArray_real_T *yy)
{
  int32_T nx;
  int32_T ny;
  int32_T j;
  int32_T i;
  emlrtStack st;
  emlrtStack b_st;
  st.prev = sp;
  st.tls = sp->tls;
  b_st.prev = &st;
  b_st.tls = st.tls;
  nx = x->size[0] * x->size[1];
  ny = y->size[0] * y->size[1];
  j = xx->size[0] * xx->size[1];
  xx->size[0] = ny;
  xx->size[1] = nx;
  emxEnsureCapacity_real_T(sp, xx, j, &s_emlrtRTEI);
  j = yy->size[0] * yy->size[1];
  yy->size[0] = ny;
  yy->size[1] = nx;
  emxEnsureCapacity_real_T(sp, yy, j, &s_emlrtRTEI);
  if ((nx == 0) || (ny == 0)) {
  } else {
    st.site = &e_emlrtRSI;
    if ((1 <= nx) && (nx > 2147483646)) {
      b_st.site = &g_emlrtRSI;
      check_forloop_overflow_error(&b_st);
    }

    for (j = 0; j < nx; j++) {
      st.site = &f_emlrtRSI;
      if ((1 <= ny) && (ny > 2147483646)) {
        b_st.site = &g_emlrtRSI;
        check_forloop_overflow_error(&b_st);
      }

      for (i = 0; i < ny; i++) {
        xx->data[i + xx->size[0] * j] = x->data[j];
        yy->data[i + yy->size[0] * j] = y->data[i];
      }
    }
  }
}

/* End of code generation (meshgrid.c) */
