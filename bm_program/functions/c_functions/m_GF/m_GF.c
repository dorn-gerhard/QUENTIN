/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * m_GF.c
 *
 * Code generation for function 'm_GF'
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "m_GF.h"
#include "rdivide_helper.h"
#include "m_GF_emxutil.h"
#include "trace.h"
#include "power.h"
#include "meshgrid.h"
#include "m_GF_data.h"
#include "blas.h"

/* Variable Definitions */
static emlrtRSInfo emlrtRSI = { 4,     /* lineNo */
  "m_GF",                              /* fcnName */
  "/itp/MooseFS/dorn/run/bm_program/functions/c_functions/m_GF.m"/* pathName */
};

static emlrtRSInfo b_emlrtRSI = { 5,   /* lineNo */
  "m_GF",                              /* fcnName */
  "/itp/MooseFS/dorn/run/bm_program/functions/c_functions/m_GF.m"/* pathName */
};

static emlrtRSInfo c_emlrtRSI = { 9,   /* lineNo */
  "m_GF",                              /* fcnName */
  "/itp/MooseFS/dorn/run/bm_program/functions/c_functions/m_GF.m"/* pathName */
};

static emlrtRSInfo d_emlrtRSI = { 10,  /* lineNo */
  "m_GF",                              /* fcnName */
  "/itp/MooseFS/dorn/run/bm_program/functions/c_functions/m_GF.m"/* pathName */
};

static emlrtRSInfo i_emlrtRSI = { 52,  /* lineNo */
  "eml_mtimes_helper",                 /* fcnName */
  "/afs/itp.tugraz.at/opt/matlab/R2018b/toolbox/eml/lib/matlab/ops/eml_mtimes_helper.m"/* pathName */
};

static emlrtRSInfo j_emlrtRSI = { 21,  /* lineNo */
  "eml_mtimes_helper",                 /* fcnName */
  "/afs/itp.tugraz.at/opt/matlab/R2018b/toolbox/eml/lib/matlab/ops/eml_mtimes_helper.m"/* pathName */
};

static emlrtRSInfo k_emlrtRSI = { 114, /* lineNo */
  "mtimes",                            /* fcnName */
  "/afs/itp.tugraz.at/opt/matlab/R2018b/toolbox/eml/eml/+coder/+internal/+blas/mtimes.m"/* pathName */
};

static emlrtRSInfo l_emlrtRSI = { 118, /* lineNo */
  "mtimes",                            /* fcnName */
  "/afs/itp.tugraz.at/opt/matlab/R2018b/toolbox/eml/eml/+coder/+internal/+blas/mtimes.m"/* pathName */
};

static emlrtRSInfo p_emlrtRSI = { 40,  /* lineNo */
  "mpower",                            /* fcnName */
  "/afs/itp.tugraz.at/opt/matlab/R2018b/toolbox/eml/lib/matlab/ops/mpower.m"/* pathName */
};

static emlrtRTEInfo emlrtRTEI = { 6,   /* lineNo */
  1,                                   /* colNo */
  "m_GF",                              /* fName */
  "/itp/MooseFS/dorn/run/bm_program/functions/c_functions/m_GF.m"/* pName */
};

static emlrtRTEInfo b_emlrtRTEI = { 7, /* lineNo */
  1,                                   /* colNo */
  "m_GF",                              /* fName */
  "/itp/MooseFS/dorn/run/bm_program/functions/c_functions/m_GF.m"/* pName */
};

static emlrtRTEInfo c_emlrtRTEI = { 9, /* lineNo */
  31,                                  /* colNo */
  "m_GF",                              /* fName */
  "/itp/MooseFS/dorn/run/bm_program/functions/c_functions/m_GF.m"/* pName */
};

static emlrtRTEInfo d_emlrtRTEI = { 9, /* lineNo */
  93,                                  /* colNo */
  "m_GF",                              /* fName */
  "/itp/MooseFS/dorn/run/bm_program/functions/c_functions/m_GF.m"/* pName */
};

static emlrtRTEInfo e_emlrtRTEI = { 9, /* lineNo */
  27,                                  /* colNo */
  "m_GF",                              /* fName */
  "/itp/MooseFS/dorn/run/bm_program/functions/c_functions/m_GF.m"/* pName */
};

static emlrtRTEInfo f_emlrtRTEI = { 9, /* lineNo */
  20,                                  /* colNo */
  "m_GF",                              /* fName */
  "/itp/MooseFS/dorn/run/bm_program/functions/c_functions/m_GF.m"/* pName */
};

static emlrtRTEInfo g_emlrtRTEI = { 118,/* lineNo */
  13,                                  /* colNo */
  "mtimes",                            /* fName */
  "/afs/itp.tugraz.at/opt/matlab/R2018b/toolbox/eml/eml/+coder/+internal/+blas/mtimes.m"/* pName */
};

static emlrtRTEInfo h_emlrtRTEI = { 9, /* lineNo */
  77,                                  /* colNo */
  "m_GF",                              /* fName */
  "/itp/MooseFS/dorn/run/bm_program/functions/c_functions/m_GF.m"/* pName */
};

static emlrtRTEInfo i_emlrtRTEI = { 9, /* lineNo */
  89,                                  /* colNo */
  "m_GF",                              /* fName */
  "/itp/MooseFS/dorn/run/bm_program/functions/c_functions/m_GF.m"/* pName */
};

static emlrtRTEInfo j_emlrtRTEI = { 10,/* lineNo */
  43,                                  /* colNo */
  "m_GF",                              /* fName */
  "/itp/MooseFS/dorn/run/bm_program/functions/c_functions/m_GF.m"/* pName */
};

static emlrtRTEInfo k_emlrtRTEI = { 10,/* lineNo */
  42,                                  /* colNo */
  "m_GF",                              /* fName */
  "/itp/MooseFS/dorn/run/bm_program/functions/c_functions/m_GF.m"/* pName */
};

static emlrtRTEInfo l_emlrtRTEI = { 10,/* lineNo */
  124,                                 /* colNo */
  "m_GF",                              /* fName */
  "/itp/MooseFS/dorn/run/bm_program/functions/c_functions/m_GF.m"/* pName */
};

static emlrtRTEInfo m_emlrtRTEI = { 10,/* lineNo */
  123,                                 /* colNo */
  "m_GF",                              /* fName */
  "/itp/MooseFS/dorn/run/bm_program/functions/c_functions/m_GF.m"/* pName */
};

static emlrtRTEInfo n_emlrtRTEI = { 10,/* lineNo */
  20,                                  /* colNo */
  "m_GF",                              /* fName */
  "/itp/MooseFS/dorn/run/bm_program/functions/c_functions/m_GF.m"/* pName */
};

static emlrtRTEInfo o_emlrtRTEI = { 10,/* lineNo */
  37,                                  /* colNo */
  "m_GF",                              /* fName */
  "/itp/MooseFS/dorn/run/bm_program/functions/c_functions/m_GF.m"/* pName */
};

static emlrtRTEInfo p_emlrtRTEI = { 10,/* lineNo */
  97,                                  /* colNo */
  "m_GF",                              /* fName */
  "/itp/MooseFS/dorn/run/bm_program/functions/c_functions/m_GF.m"/* pName */
};

static emlrtRTEInfo q_emlrtRTEI = { 10,/* lineNo */
  118,                                 /* colNo */
  "m_GF",                              /* fName */
  "/itp/MooseFS/dorn/run/bm_program/functions/c_functions/m_GF.m"/* pName */
};

static emlrtRTEInfo r_emlrtRTEI = { 1, /* lineNo */
  21,                                  /* colNo */
  "m_GF",                              /* fName */
  "/itp/MooseFS/dorn/run/bm_program/functions/c_functions/m_GF.m"/* pName */
};

static emlrtECInfo emlrtECI = { 2,     /* nDims */
  9,                                   /* lineNo */
  44,                                  /* colNo */
  "m_GF",                              /* fName */
  "/itp/MooseFS/dorn/run/bm_program/functions/c_functions/m_GF.m"/* pName */
};

static emlrtECInfo b_emlrtECI = { 2,   /* nDims */
  9,                                   /* lineNo */
  27,                                  /* colNo */
  "m_GF",                              /* fName */
  "/itp/MooseFS/dorn/run/bm_program/functions/c_functions/m_GF.m"/* pName */
};

static emlrtECInfo c_emlrtECI = { 2,   /* nDims */
  9,                                   /* lineNo */
  106,                                 /* colNo */
  "m_GF",                              /* fName */
  "/itp/MooseFS/dorn/run/bm_program/functions/c_functions/m_GF.m"/* pName */
};

static emlrtECInfo d_emlrtECI = { 2,   /* nDims */
  9,                                   /* lineNo */
  89,                                  /* colNo */
  "m_GF",                              /* fName */
  "/itp/MooseFS/dorn/run/bm_program/functions/c_functions/m_GF.m"/* pName */
};

static emlrtECInfo e_emlrtECI = { 2,   /* nDims */
  9,                                   /* lineNo */
  20,                                  /* colNo */
  "m_GF",                              /* fName */
  "/itp/MooseFS/dorn/run/bm_program/functions/c_functions/m_GF.m"/* pName */
};

static emlrtECInfo f_emlrtECI = { 2,   /* nDims */
  10,                                  /* lineNo */
  56,                                  /* colNo */
  "m_GF",                              /* fName */
  "/itp/MooseFS/dorn/run/bm_program/functions/c_functions/m_GF.m"/* pName */
};

static emlrtECInfo g_emlrtECI = { 2,   /* nDims */
  10,                                  /* lineNo */
  37,                                  /* colNo */
  "m_GF",                              /* fName */
  "/itp/MooseFS/dorn/run/bm_program/functions/c_functions/m_GF.m"/* pName */
};

static emlrtECInfo h_emlrtECI = { 2,   /* nDims */
  10,                                  /* lineNo */
  137,                                 /* colNo */
  "m_GF",                              /* fName */
  "/itp/MooseFS/dorn/run/bm_program/functions/c_functions/m_GF.m"/* pName */
};

static emlrtECInfo i_emlrtECI = { 2,   /* nDims */
  10,                                  /* lineNo */
  118,                                 /* colNo */
  "m_GF",                              /* fName */
  "/itp/MooseFS/dorn/run/bm_program/functions/c_functions/m_GF.m"/* pName */
};

static emlrtECInfo j_emlrtECI = { 2,   /* nDims */
  10,                                  /* lineNo */
  20,                                  /* colNo */
  "m_GF",                              /* fName */
  "/itp/MooseFS/dorn/run/bm_program/functions/c_functions/m_GF.m"/* pName */
};

static emlrtRTEInfo w_emlrtRTEI = { 88,/* lineNo */
  23,                                  /* colNo */
  "eml_mtimes_helper",                 /* fName */
  "/afs/itp.tugraz.at/opt/matlab/R2018b/toolbox/eml/lib/matlab/ops/eml_mtimes_helper.m"/* pName */
};

static emlrtRTEInfo x_emlrtRTEI = { 83,/* lineNo */
  23,                                  /* colNo */
  "eml_mtimes_helper",                 /* fName */
  "/afs/itp.tugraz.at/opt/matlab/R2018b/toolbox/eml/lib/matlab/ops/eml_mtimes_helper.m"/* pName */
};

static emlrtBCInfo emlrtBCI = { -1,    /* iFirst */
  -1,                                  /* iLast */
  9,                                   /* lineNo */
  31,                                  /* colNo */
  "points",                            /* aName */
  "m_GF",                              /* fName */
  "/itp/MooseFS/dorn/run/bm_program/functions/c_functions/m_GF.m",/* pName */
  0                                    /* checkKind */
};

static emlrtBCInfo b_emlrtBCI = { -1,  /* iFirst */
  -1,                                  /* iLast */
  9,                                   /* lineNo */
  93,                                  /* colNo */
  "points",                            /* aName */
  "m_GF",                              /* fName */
  "/itp/MooseFS/dorn/run/bm_program/functions/c_functions/m_GF.m",/* pName */
  0                                    /* checkKind */
};

static emlrtBCInfo c_emlrtBCI = { -1,  /* iFirst */
  -1,                                  /* iLast */
  9,                                   /* lineNo */
  5,                                   /* colNo */
  "GR",                                /* aName */
  "m_GF",                              /* fName */
  "/itp/MooseFS/dorn/run/bm_program/functions/c_functions/m_GF.m",/* pName */
  0                                    /* checkKind */
};

static emlrtBCInfo d_emlrtBCI = { -1,  /* iFirst */
  -1,                                  /* iLast */
  10,                                  /* lineNo */
  43,                                  /* colNo */
  "points",                            /* aName */
  "m_GF",                              /* fName */
  "/itp/MooseFS/dorn/run/bm_program/functions/c_functions/m_GF.m",/* pName */
  0                                    /* checkKind */
};

static emlrtBCInfo e_emlrtBCI = { -1,  /* iFirst */
  -1,                                  /* iLast */
  10,                                  /* lineNo */
  124,                                 /* colNo */
  "points",                            /* aName */
  "m_GF",                              /* fName */
  "/itp/MooseFS/dorn/run/bm_program/functions/c_functions/m_GF.m",/* pName */
  0                                    /* checkKind */
};

static emlrtBCInfo f_emlrtBCI = { -1,  /* iFirst */
  -1,                                  /* iLast */
  10,                                  /* lineNo */
  5,                                   /* colNo */
  "GK",                                /* aName */
  "m_GF",                              /* fName */
  "/itp/MooseFS/dorn/run/bm_program/functions/c_functions/m_GF.m",/* pName */
  0                                    /* checkKind */
};

/* Function Definitions */
void m_GF(const emlrtStack *sp, const emxArray_real_T *points, real_T i0p, const
          emxArray_real_T *E_p, const emxArray_real_T *E_n, const
          emxArray_real_T *E_q, const emxArray_creal_T *C1, const
          emxArray_creal_T *C2, const emxArray_creal_T *C3, const
          emxArray_creal_T *C4, const emxArray_creal_T *rho, emxArray_creal_T
          *GR, emxArray_creal_T *GK)
{
  emxArray_real_T *E_pp;
  emxArray_real_T *E_nn;
  emxArray_real_T *E_mm;
  emxArray_real_T *E_qq;
  int32_T i0;
  int32_T loop_ub;
  emxArray_creal_T *r0;
  emxArray_creal_T *r1;
  emxArray_creal_T *b;
  emxArray_creal_T *a;
  emxArray_creal_T *b_b;
  emxArray_real_T *r2;
  emxArray_real_T *r3;
  emxArray_real_T *b_points;
  emxArray_creal_T *c_points;
  int32_T k;
  int32_T b_E_pp[2];
  int32_T b_E_nn[2];
  real_T y_re;
  int32_T i1;
  int32_T i2;
  real_T d_points;
  real_T b_re;
  real_T y_im;
  ptrdiff_t m_t;
  ptrdiff_t n_t;
  int32_T b_loop_ub;
  ptrdiff_t k_t;
  ptrdiff_t lda_t;
  static const creal_T beta1 = { 0.0,  /* re */
    0.0                                /* im */
  };

  ptrdiff_t ldb_t;
  ptrdiff_t ldc_t;
  int32_T c_loop_ub;
  int32_T i3;
  char_T TRANSA;
  char_T TRANSB;
  static const creal_T alpha1 = { 1.0, /* re */
    0.0                                /* im */
  };

  int32_T iv0[2];
  emlrtStack *r4;
  emlrtStack st;
  emlrtStack b_st;
  emlrtStack c_st;
  emlrtStack d_st;
  st.prev = sp;
  st.tls = sp->tls;
  b_st.prev = &st;
  b_st.tls = st.tls;
  c_st.prev = &b_st;
  c_st.tls = b_st.tls;
  d_st.prev = &c_st;
  d_st.tls = c_st.tls;
  emlrtHeapReferenceStackEnterFcnR2012b(sp);
  emxInit_real_T(sp, &E_pp, 2, &r_emlrtRTEI, true);
  emxInit_real_T(sp, &E_nn, 2, &r_emlrtRTEI, true);
  emxInit_real_T(sp, &E_mm, 2, &r_emlrtRTEI, true);
  emxInit_real_T(sp, &E_qq, 2, &r_emlrtRTEI, true);
  st.site = &emlrtRSI;
  meshgrid(&st, E_p, E_n, E_pp, E_nn);
  st.site = &b_emlrtRSI;
  meshgrid(&st, E_n, E_q, E_mm, E_qq);
  i0 = GR->size[0];
  GR->size[0] = points->size[0];
  emxEnsureCapacity_creal_T(sp, GR, i0, &emlrtRTEI);
  loop_ub = points->size[0];
  for (i0 = 0; i0 < loop_ub; i0++) {
    GR->data[i0].re = 0.0;
    GR->data[i0].im = 0.0;
  }

  i0 = GK->size[0];
  GK->size[0] = points->size[0];
  emxEnsureCapacity_creal_T(sp, GK, i0, &b_emlrtRTEI);
  loop_ub = points->size[0];
  for (i0 = 0; i0 < loop_ub; i0++) {
    GK->data[i0].re = 0.0;
    GK->data[i0].im = 0.0;
  }

  i0 = points->size[0];
  emxInit_creal_T(sp, &r0, 2, &r_emlrtRTEI, true);
  emxInit_creal_T(sp, &r1, 2, &r_emlrtRTEI, true);
  emxInit_creal_T(sp, &b, 2, &e_emlrtRTEI, true);
  emxInit_creal_T(sp, &a, 2, &f_emlrtRTEI, true);
  emxInit_creal_T(sp, &b_b, 2, &i_emlrtRTEI, true);
  emxInit_real_T(sp, &r2, 2, &k_emlrtRTEI, true);
  emxInit_real_T(sp, &r3, 2, &k_emlrtRTEI, true);
  emxInit_real_T(sp, &b_points, 2, &j_emlrtRTEI, true);
  emxInit_creal_T(sp, &c_points, 2, &c_emlrtRTEI, true);
  for (k = 0; k < i0; k++) {
    b_E_pp[0] = E_pp->size[0];
    b_E_pp[1] = E_pp->size[1];
    b_E_nn[0] = E_nn->size[0];
    b_E_nn[1] = E_nn->size[1];
    if ((b_E_pp[0] != b_E_nn[0]) || (b_E_pp[1] != b_E_nn[1])) {
      emlrtSizeEqCheckNDR2012b(&b_E_pp[0], &b_E_nn[0], &emlrtECI, sp);
    }

    y_re = i0p * 0.0;
    i1 = c_points->size[0] * c_points->size[1];
    c_points->size[0] = E_pp->size[0];
    c_points->size[1] = E_pp->size[1];
    emxEnsureCapacity_creal_T(sp, c_points, i1, &c_emlrtRTEI);
    i1 = points->size[0];
    i2 = 1 + k;
    if ((i2 < 1) || (i2 > i1)) {
      emlrtDynamicBoundsCheckR2012b(i2, 1, i1, &emlrtBCI, sp);
    }

    d_points = points->data[i2 - 1];
    loop_ub = E_pp->size[0] * E_pp->size[1];
    for (i1 = 0; i1 < loop_ub; i1++) {
      c_points->data[i1].re = (d_points - (E_pp->data[i1] - E_nn->data[i1])) +
        y_re;
      c_points->data[i1].im = i0p;
    }

    rdivide_helper(sp, c_points, b);
    b_E_pp[0] = b->size[0];
    b_E_pp[1] = b->size[1];
    b_E_nn[0] = C1->size[0];
    b_E_nn[1] = C1->size[1];
    if ((b_E_pp[0] != b_E_nn[0]) || (b_E_pp[1] != b_E_nn[1])) {
      emlrtSizeEqCheckNDR2012b(&b_E_pp[0], &b_E_nn[0], &b_emlrtECI, sp);
    }

    b_E_pp[0] = E_mm->size[0];
    b_E_pp[1] = E_mm->size[1];
    b_E_nn[0] = E_qq->size[0];
    b_E_nn[1] = E_qq->size[1];
    if ((b_E_pp[0] != b_E_nn[0]) || (b_E_pp[1] != b_E_nn[1])) {
      emlrtSizeEqCheckNDR2012b(&b_E_pp[0], &b_E_nn[0], &c_emlrtECI, sp);
    }

    y_re = i0p * 0.0;
    i1 = c_points->size[0] * c_points->size[1];
    c_points->size[0] = E_mm->size[0];
    c_points->size[1] = E_mm->size[1];
    emxEnsureCapacity_creal_T(sp, c_points, i1, &d_emlrtRTEI);
    i1 = points->size[0];
    i2 = 1 + k;
    if ((i2 < 1) || (i2 > i1)) {
      emlrtDynamicBoundsCheckR2012b(i2, 1, i1, &b_emlrtBCI, sp);
    }

    d_points = points->data[i2 - 1];
    loop_ub = E_mm->size[0] * E_mm->size[1];
    for (i1 = 0; i1 < loop_ub; i1++) {
      c_points->data[i1].re = (d_points - (E_mm->data[i1] - E_qq->data[i1])) +
        y_re;
      c_points->data[i1].im = i0p;
    }

    rdivide_helper(sp, c_points, b_b);
    b_E_pp[0] = b_b->size[0];
    b_E_pp[1] = b_b->size[1];
    b_E_nn[0] = C4->size[0];
    b_E_nn[1] = C4->size[1];
    if ((b_E_pp[0] != b_E_nn[0]) || (b_E_pp[1] != b_E_nn[1])) {
      emlrtSizeEqCheckNDR2012b(&b_E_pp[0], &b_E_nn[0], &d_emlrtECI, sp);
    }

    st.site = &c_emlrtRSI;
    i1 = b->size[0] * b->size[1];
    i2 = b->size[0] * b->size[1];
    emxEnsureCapacity_creal_T(&st, b, i2, &e_emlrtRTEI);
    loop_ub = i1 - 1;
    for (i1 = 0; i1 <= loop_ub; i1++) {
      b_re = b->data[i1].re;
      y_re = b->data[i1].im;
      d_points = C1->data[i1].re;
      y_im = C1->data[i1].im;
      b->data[i1].re = b_re * d_points - y_re * y_im;
      b->data[i1].im = b_re * y_im + y_re * d_points;
    }

    b_st.site = &j_emlrtRSI;
    if (rho->size[1] != b->size[0]) {
      if (((rho->size[0] == 1) && (rho->size[1] == 1)) || ((b->size[0] == 1) &&
           (b->size[1] == 1))) {
        emlrtErrorWithMessageIdR2018a(&b_st, &x_emlrtRTEI,
          "Coder:toolbox:mtimes_noDynamicScalarExpansion",
          "Coder:toolbox:mtimes_noDynamicScalarExpansion", 0);
      } else {
        emlrtErrorWithMessageIdR2018a(&b_st, &w_emlrtRTEI, "MATLAB:innerdim",
          "MATLAB:innerdim", 0);
      }
    }

    if ((rho->size[1] == 1) || (b->size[0] == 1)) {
      i1 = a->size[0] * a->size[1];
      a->size[0] = rho->size[0];
      a->size[1] = b->size[1];
      emxEnsureCapacity_creal_T(&st, a, i1, &f_emlrtRTEI);
      loop_ub = rho->size[0];
      for (i1 = 0; i1 < loop_ub; i1++) {
        b_loop_ub = b->size[1];
        for (i2 = 0; i2 < b_loop_ub; i2++) {
          a->data[i1 + a->size[0] * i2].re = 0.0;
          a->data[i1 + a->size[0] * i2].im = 0.0;
          c_loop_ub = rho->size[1];
          for (i3 = 0; i3 < c_loop_ub; i3++) {
            d_points = rho->data[i1 + rho->size[0] * i3].re * b->data[i3 +
              b->size[0] * i2].re - rho->data[i1 + rho->size[0] * i3].im *
              b->data[i3 + b->size[0] * i2].im;
            b_re = rho->data[i1 + rho->size[0] * i3].re * b->data[i3 + b->size[0]
              * i2].im + rho->data[i1 + rho->size[0] * i3].im * b->data[i3 +
              b->size[0] * i2].re;
            a->data[i1 + a->size[0] * i2].re += d_points;
            a->data[i1 + a->size[0] * i2].im += b_re;
          }
        }
      }
    } else {
      b_st.site = &i_emlrtRSI;
      if ((rho->size[0] == 0) || (rho->size[1] == 0) || (b->size[0] == 0) ||
          (b->size[1] == 0)) {
        i1 = a->size[0] * a->size[1];
        a->size[0] = rho->size[0];
        a->size[1] = b->size[1];
        emxEnsureCapacity_creal_T(&b_st, a, i1, &f_emlrtRTEI);
        loop_ub = rho->size[0] * b->size[1];
        for (i1 = 0; i1 < loop_ub; i1++) {
          a->data[i1] = beta1;
        }
      } else {
        c_st.site = &k_emlrtRSI;
        c_st.site = &l_emlrtRSI;
        m_t = (ptrdiff_t)rho->size[0];
        n_t = (ptrdiff_t)b->size[1];
        k_t = (ptrdiff_t)rho->size[1];
        lda_t = (ptrdiff_t)rho->size[0];
        ldb_t = (ptrdiff_t)rho->size[1];
        ldc_t = (ptrdiff_t)rho->size[0];
        i1 = a->size[0] * a->size[1];
        a->size[0] = rho->size[0];
        a->size[1] = b->size[1];
        emxEnsureCapacity_creal_T(&c_st, a, i1, &g_emlrtRTEI);
        TRANSA = 'N';
        TRANSB = 'N';
        zgemm(&TRANSA, &TRANSB, &m_t, &n_t, &k_t, (real_T *)&alpha1, (real_T *)
              &rho->data[0], &lda_t, (real_T *)&b->data[0], &ldb_t, (real_T *)
              &beta1, (real_T *)&a->data[0], &ldc_t);
      }
    }

    st.site = &c_emlrtRSI;
    b_st.site = &j_emlrtRSI;
    if (a->size[1] != C2->size[0]) {
      if (((a->size[0] == 1) && (a->size[1] == 1)) || ((C2->size[0] == 1) &&
           (C2->size[1] == 1))) {
        emlrtErrorWithMessageIdR2018a(&b_st, &x_emlrtRTEI,
          "Coder:toolbox:mtimes_noDynamicScalarExpansion",
          "Coder:toolbox:mtimes_noDynamicScalarExpansion", 0);
      } else {
        emlrtErrorWithMessageIdR2018a(&b_st, &w_emlrtRTEI, "MATLAB:innerdim",
          "MATLAB:innerdim", 0);
      }
    }

    if ((a->size[1] == 1) || (C2->size[0] == 1)) {
      i1 = r0->size[0] * r0->size[1];
      r0->size[0] = a->size[0];
      r0->size[1] = C2->size[1];
      emxEnsureCapacity_creal_T(&st, r0, i1, &f_emlrtRTEI);
      loop_ub = a->size[0];
      for (i1 = 0; i1 < loop_ub; i1++) {
        b_loop_ub = C2->size[1];
        for (i2 = 0; i2 < b_loop_ub; i2++) {
          r0->data[i1 + r0->size[0] * i2].re = 0.0;
          r0->data[i1 + r0->size[0] * i2].im = 0.0;
          c_loop_ub = a->size[1];
          for (i3 = 0; i3 < c_loop_ub; i3++) {
            y_im = a->data[i1 + a->size[0] * i3].re * C2->data[i3 + C2->size[0] *
              i2].re - a->data[i1 + a->size[0] * i3].im * C2->data[i3 + C2->
              size[0] * i2].im;
            d_points = a->data[i1 + a->size[0] * i3].re * C2->data[i3 + C2->
              size[0] * i2].im + a->data[i1 + a->size[0] * i3].im * C2->data[i3
              + C2->size[0] * i2].re;
            r0->data[i1 + r0->size[0] * i2].re += y_im;
            r0->data[i1 + r0->size[0] * i2].im += d_points;
          }
        }
      }
    } else {
      b_st.site = &i_emlrtRSI;
      if ((a->size[0] == 0) || (a->size[1] == 0) || (C2->size[0] == 0) ||
          (C2->size[1] == 0)) {
        i1 = r0->size[0] * r0->size[1];
        r0->size[0] = a->size[0];
        r0->size[1] = C2->size[1];
        emxEnsureCapacity_creal_T(&b_st, r0, i1, &f_emlrtRTEI);
        loop_ub = a->size[0] * C2->size[1];
        for (i1 = 0; i1 < loop_ub; i1++) {
          r0->data[i1] = beta1;
        }
      } else {
        c_st.site = &k_emlrtRSI;
        c_st.site = &l_emlrtRSI;
        m_t = (ptrdiff_t)a->size[0];
        n_t = (ptrdiff_t)C2->size[1];
        k_t = (ptrdiff_t)a->size[1];
        lda_t = (ptrdiff_t)a->size[0];
        ldb_t = (ptrdiff_t)a->size[1];
        ldc_t = (ptrdiff_t)a->size[0];
        i1 = r0->size[0] * r0->size[1];
        r0->size[0] = a->size[0];
        r0->size[1] = C2->size[1];
        emxEnsureCapacity_creal_T(&c_st, r0, i1, &g_emlrtRTEI);
        TRANSA = 'N';
        TRANSB = 'N';
        zgemm(&TRANSA, &TRANSB, &m_t, &n_t, &k_t, (real_T *)&alpha1, (real_T *)
              &a->data[0], &lda_t, (real_T *)&C2->data[0], &ldb_t, (real_T *)
              &beta1, (real_T *)&r0->data[0], &ldc_t);
      }
    }

    st.site = &c_emlrtRSI;
    b_st.site = &j_emlrtRSI;
    if (rho->size[1] != C3->size[0]) {
      if (((rho->size[0] == 1) && (rho->size[1] == 1)) || ((C3->size[0] == 1) &&
           (C3->size[1] == 1))) {
        emlrtErrorWithMessageIdR2018a(&b_st, &x_emlrtRTEI,
          "Coder:toolbox:mtimes_noDynamicScalarExpansion",
          "Coder:toolbox:mtimes_noDynamicScalarExpansion", 0);
      } else {
        emlrtErrorWithMessageIdR2018a(&b_st, &w_emlrtRTEI, "MATLAB:innerdim",
          "MATLAB:innerdim", 0);
      }
    }

    if ((rho->size[1] == 1) || (C3->size[0] == 1)) {
      i1 = a->size[0] * a->size[1];
      a->size[0] = rho->size[0];
      a->size[1] = C3->size[1];
      emxEnsureCapacity_creal_T(&st, a, i1, &h_emlrtRTEI);
      loop_ub = rho->size[0];
      for (i1 = 0; i1 < loop_ub; i1++) {
        b_loop_ub = C3->size[1];
        for (i2 = 0; i2 < b_loop_ub; i2++) {
          a->data[i1 + a->size[0] * i2].re = 0.0;
          a->data[i1 + a->size[0] * i2].im = 0.0;
          c_loop_ub = rho->size[1];
          for (i3 = 0; i3 < c_loop_ub; i3++) {
            d_points = rho->data[i1 + rho->size[0] * i3].re * C3->data[i3 +
              C3->size[0] * i2].re - rho->data[i1 + rho->size[0] * i3].im *
              C3->data[i3 + C3->size[0] * i2].im;
            b_re = rho->data[i1 + rho->size[0] * i3].re * C3->data[i3 + C3->
              size[0] * i2].im + rho->data[i1 + rho->size[0] * i3].im * C3->
              data[i3 + C3->size[0] * i2].re;
            a->data[i1 + a->size[0] * i2].re += d_points;
            a->data[i1 + a->size[0] * i2].im += b_re;
          }
        }
      }
    } else {
      b_st.site = &i_emlrtRSI;
      if ((rho->size[0] == 0) || (rho->size[1] == 0) || (C3->size[0] == 0) ||
          (C3->size[1] == 0)) {
        i1 = a->size[0] * a->size[1];
        a->size[0] = rho->size[0];
        a->size[1] = C3->size[1];
        emxEnsureCapacity_creal_T(&b_st, a, i1, &h_emlrtRTEI);
        loop_ub = rho->size[0] * C3->size[1];
        for (i1 = 0; i1 < loop_ub; i1++) {
          a->data[i1] = beta1;
        }
      } else {
        c_st.site = &k_emlrtRSI;
        c_st.site = &l_emlrtRSI;
        m_t = (ptrdiff_t)rho->size[0];
        n_t = (ptrdiff_t)C3->size[1];
        k_t = (ptrdiff_t)rho->size[1];
        lda_t = (ptrdiff_t)rho->size[0];
        ldb_t = (ptrdiff_t)rho->size[1];
        ldc_t = (ptrdiff_t)rho->size[0];
        i1 = a->size[0] * a->size[1];
        a->size[0] = rho->size[0];
        a->size[1] = C3->size[1];
        emxEnsureCapacity_creal_T(&c_st, a, i1, &g_emlrtRTEI);
        TRANSA = 'N';
        TRANSB = 'N';
        zgemm(&TRANSA, &TRANSB, &m_t, &n_t, &k_t, (real_T *)&alpha1, (real_T *)
              &rho->data[0], &lda_t, (real_T *)&C3->data[0], &ldb_t, (real_T *)
              &beta1, (real_T *)&a->data[0], &ldc_t);
      }
    }

    st.site = &c_emlrtRSI;
    i1 = b_b->size[0] * b_b->size[1];
    i2 = b_b->size[0] * b_b->size[1];
    emxEnsureCapacity_creal_T(&st, b_b, i2, &i_emlrtRTEI);
    loop_ub = i1 - 1;
    for (i1 = 0; i1 <= loop_ub; i1++) {
      b_re = b_b->data[i1].re;
      y_re = b_b->data[i1].im;
      d_points = C4->data[i1].re;
      y_im = C4->data[i1].im;
      b_b->data[i1].re = b_re * d_points - y_re * y_im;
      b_b->data[i1].im = b_re * y_im + y_re * d_points;
    }

    b_st.site = &j_emlrtRSI;
    if (a->size[1] != b_b->size[0]) {
      if (((a->size[0] == 1) && (a->size[1] == 1)) || ((b_b->size[0] == 1) &&
           (b_b->size[1] == 1))) {
        emlrtErrorWithMessageIdR2018a(&b_st, &x_emlrtRTEI,
          "Coder:toolbox:mtimes_noDynamicScalarExpansion",
          "Coder:toolbox:mtimes_noDynamicScalarExpansion", 0);
      } else {
        emlrtErrorWithMessageIdR2018a(&b_st, &w_emlrtRTEI, "MATLAB:innerdim",
          "MATLAB:innerdim", 0);
      }
    }

    if ((a->size[1] == 1) || (b_b->size[0] == 1)) {
      i1 = r1->size[0] * r1->size[1];
      r1->size[0] = a->size[0];
      r1->size[1] = b_b->size[1];
      emxEnsureCapacity_creal_T(&st, r1, i1, &h_emlrtRTEI);
      loop_ub = a->size[0];
      for (i1 = 0; i1 < loop_ub; i1++) {
        b_loop_ub = b_b->size[1];
        for (i2 = 0; i2 < b_loop_ub; i2++) {
          r1->data[i1 + r1->size[0] * i2].re = 0.0;
          r1->data[i1 + r1->size[0] * i2].im = 0.0;
          c_loop_ub = a->size[1];
          for (i3 = 0; i3 < c_loop_ub; i3++) {
            y_im = a->data[i1 + a->size[0] * i3].re * b_b->data[i3 + b_b->size[0]
              * i2].re - a->data[i1 + a->size[0] * i3].im * b_b->data[i3 +
              b_b->size[0] * i2].im;
            d_points = a->data[i1 + a->size[0] * i3].re * b_b->data[i3 +
              b_b->size[0] * i2].im + a->data[i1 + a->size[0] * i3].im *
              b_b->data[i3 + b_b->size[0] * i2].re;
            r1->data[i1 + r1->size[0] * i2].re += y_im;
            r1->data[i1 + r1->size[0] * i2].im += d_points;
          }
        }
      }
    } else {
      b_st.site = &i_emlrtRSI;
      if ((a->size[0] == 0) || (a->size[1] == 0) || (b_b->size[0] == 0) ||
          (b_b->size[1] == 0)) {
        i1 = r1->size[0] * r1->size[1];
        r1->size[0] = a->size[0];
        r1->size[1] = b_b->size[1];
        emxEnsureCapacity_creal_T(&b_st, r1, i1, &h_emlrtRTEI);
        loop_ub = a->size[0] * b_b->size[1];
        for (i1 = 0; i1 < loop_ub; i1++) {
          r1->data[i1] = beta1;
        }
      } else {
        c_st.site = &k_emlrtRSI;
        c_st.site = &l_emlrtRSI;
        m_t = (ptrdiff_t)a->size[0];
        n_t = (ptrdiff_t)b_b->size[1];
        k_t = (ptrdiff_t)a->size[1];
        lda_t = (ptrdiff_t)a->size[0];
        ldb_t = (ptrdiff_t)a->size[1];
        ldc_t = (ptrdiff_t)a->size[0];
        i1 = r1->size[0] * r1->size[1];
        r1->size[0] = a->size[0];
        r1->size[1] = b_b->size[1];
        emxEnsureCapacity_creal_T(&c_st, r1, i1, &g_emlrtRTEI);
        TRANSA = 'N';
        TRANSB = 'N';
        zgemm(&TRANSA, &TRANSB, &m_t, &n_t, &k_t, (real_T *)&alpha1, (real_T *)
              &a->data[0], &lda_t, (real_T *)&b_b->data[0], &ldb_t, (real_T *)
              &beta1, (real_T *)&r1->data[0], &ldc_t);
      }
    }

    iv0[0] = r0->size[0];
    iv0[1] = r0->size[1];
    b_E_pp[0] = r1->size[0];
    b_E_pp[1] = r1->size[1];
    if ((iv0[0] != b_E_pp[0]) || (iv0[1] != b_E_pp[1])) {
      emlrtSizeEqCheckNDR2012b(&iv0[0], &b_E_pp[0], &e_emlrtECI, sp);
    }

    i1 = c_points->size[0] * c_points->size[1];
    c_points->size[0] = r0->size[0];
    c_points->size[1] = r0->size[1];
    emxEnsureCapacity_creal_T(sp, c_points, i1, &f_emlrtRTEI);
    loop_ub = r0->size[0] * r0->size[1];
    for (i1 = 0; i1 < loop_ub; i1++) {
      c_points->data[i1].re = r0->data[i1].re + r1->data[i1].re;
      c_points->data[i1].im = r0->data[i1].im + r1->data[i1].im;
    }

    st.site = &c_emlrtRSI;
    r4 = &st;
    i1 = GR->size[0];
    i2 = 1 + k;
    if ((i2 < 1) || (i2 > i1)) {
      emlrtDynamicBoundsCheckR2012b(i2, 1, i1, &c_emlrtBCI, &st);
    }

    GR->data[i2 - 1] = trace(r4, c_points);
    b_E_pp[0] = E_pp->size[0];
    b_E_pp[1] = E_pp->size[1];
    b_E_nn[0] = E_nn->size[0];
    b_E_nn[1] = E_nn->size[1];
    if ((b_E_pp[0] != b_E_nn[0]) || (b_E_pp[1] != b_E_nn[1])) {
      emlrtSizeEqCheckNDR2012b(&b_E_pp[0], &b_E_nn[0], &f_emlrtECI, &st);
    }

    b_st.site = &d_emlrtRSI;
    c_st.site = &p_emlrtRSI;
    y_im = i0p * i0p;
    i1 = b_points->size[0] * b_points->size[1];
    b_points->size[0] = E_pp->size[0];
    b_points->size[1] = E_pp->size[1];
    emxEnsureCapacity_real_T(&st, b_points, i1, &j_emlrtRTEI);
    i1 = points->size[0];
    i2 = 1 + k;
    if ((i2 < 1) || (i2 > i1)) {
      emlrtDynamicBoundsCheckR2012b(i2, 1, i1, &d_emlrtBCI, &st);
    }

    d_points = points->data[i2 - 1];
    loop_ub = E_pp->size[0] * E_pp->size[1];
    for (i1 = 0; i1 < loop_ub; i1++) {
      b_points->data[i1] = d_points - (E_pp->data[i1] - E_nn->data[i1]);
    }

    b_st.site = &d_emlrtRSI;
    power(&b_st, b_points, r2);
    i1 = r3->size[0] * r3->size[1];
    r3->size[0] = r2->size[0];
    r3->size[1] = r2->size[1];
    emxEnsureCapacity_real_T(&st, r3, i1, &k_emlrtRTEI);
    loop_ub = r2->size[0] * r2->size[1];
    for (i1 = 0; i1 < loop_ub; i1++) {
      r3->data[i1] = r2->data[i1] + y_im;
    }

    b_rdivide_helper(&st, r3, r2);
    iv0[0] = r2->size[0];
    iv0[1] = r2->size[1];
    b_E_nn[0] = C1->size[0];
    b_E_nn[1] = C1->size[1];
    if ((iv0[0] != b_E_nn[0]) || (iv0[1] != b_E_nn[1])) {
      emlrtSizeEqCheckNDR2012b(&iv0[0], &b_E_nn[0], &g_emlrtECI, &st);
    }

    b_E_pp[0] = E_mm->size[0];
    b_E_pp[1] = E_mm->size[1];
    b_E_nn[0] = E_qq->size[0];
    b_E_nn[1] = E_qq->size[1];
    if ((b_E_pp[0] != b_E_nn[0]) || (b_E_pp[1] != b_E_nn[1])) {
      emlrtSizeEqCheckNDR2012b(&b_E_pp[0], &b_E_nn[0], &h_emlrtECI, &st);
    }

    b_st.site = &d_emlrtRSI;
    c_st.site = &p_emlrtRSI;
    i1 = b_points->size[0] * b_points->size[1];
    b_points->size[0] = E_mm->size[0];
    b_points->size[1] = E_mm->size[1];
    emxEnsureCapacity_real_T(&st, b_points, i1, &l_emlrtRTEI);
    i1 = points->size[0];
    i2 = 1 + k;
    if ((i2 < 1) || (i2 > i1)) {
      emlrtDynamicBoundsCheckR2012b(i2, 1, i1, &e_emlrtBCI, &st);
    }

    d_points = points->data[i2 - 1];
    loop_ub = E_mm->size[0] * E_mm->size[1];
    for (i1 = 0; i1 < loop_ub; i1++) {
      b_points->data[i1] = d_points - (E_mm->data[i1] - E_qq->data[i1]);
    }

    b_st.site = &d_emlrtRSI;
    power(&b_st, b_points, r3);
    i1 = b_points->size[0] * b_points->size[1];
    b_points->size[0] = r3->size[0];
    b_points->size[1] = r3->size[1];
    emxEnsureCapacity_real_T(&st, b_points, i1, &m_emlrtRTEI);
    loop_ub = r3->size[0] * r3->size[1];
    for (i1 = 0; i1 < loop_ub; i1++) {
      b_points->data[i1] = r3->data[i1] + y_im;
    }

    b_rdivide_helper(&st, b_points, r3);
    iv0[0] = r3->size[0];
    iv0[1] = r3->size[1];
    b_E_nn[0] = C4->size[0];
    b_E_nn[1] = C4->size[1];
    if ((iv0[0] != b_E_nn[0]) || (iv0[1] != b_E_nn[1])) {
      emlrtSizeEqCheckNDR2012b(&iv0[0], &b_E_nn[0], &i_emlrtECI, &st);
    }

    y_re = i0p * 0.0;
    y_im = i0p * -2.0;
    i1 = a->size[0] * a->size[1];
    a->size[0] = rho->size[0];
    a->size[1] = rho->size[1];
    emxEnsureCapacity_creal_T(&st, a, i1, &n_emlrtRTEI);
    loop_ub = rho->size[0] * rho->size[1];
    for (i1 = 0; i1 < loop_ub; i1++) {
      d_points = rho->data[i1].re;
      b_re = rho->data[i1].im;
      a->data[i1].re = y_re * d_points - y_im * b_re;
      a->data[i1].im = y_re * b_re + y_im * d_points;
    }

    b_st.site = &d_emlrtRSI;
    i1 = b->size[0] * b->size[1];
    b->size[0] = r2->size[0];
    b->size[1] = r2->size[1];
    emxEnsureCapacity_creal_T(&b_st, b, i1, &o_emlrtRTEI);
    loop_ub = r2->size[0] * r2->size[1];
    for (i1 = 0; i1 < loop_ub; i1++) {
      b->data[i1].re = r2->data[i1] * C1->data[i1].re;
      b->data[i1].im = r2->data[i1] * C1->data[i1].im;
    }

    c_st.site = &j_emlrtRSI;
    if (a->size[1] != b->size[0]) {
      if (((a->size[0] == 1) && (a->size[1] == 1)) || ((b->size[0] == 1) &&
           (b->size[1] == 1))) {
        emlrtErrorWithMessageIdR2018a(&c_st, &x_emlrtRTEI,
          "Coder:toolbox:mtimes_noDynamicScalarExpansion",
          "Coder:toolbox:mtimes_noDynamicScalarExpansion", 0);
      } else {
        emlrtErrorWithMessageIdR2018a(&c_st, &w_emlrtRTEI, "MATLAB:innerdim",
          "MATLAB:innerdim", 0);
      }
    }

    if ((a->size[1] == 1) || (b->size[0] == 1)) {
      i1 = b_b->size[0] * b_b->size[1];
      b_b->size[0] = a->size[0];
      b_b->size[1] = b->size[1];
      emxEnsureCapacity_creal_T(&b_st, b_b, i1, &n_emlrtRTEI);
      loop_ub = a->size[0];
      for (i1 = 0; i1 < loop_ub; i1++) {
        b_loop_ub = b->size[1];
        for (i2 = 0; i2 < b_loop_ub; i2++) {
          b_b->data[i1 + b_b->size[0] * i2].re = 0.0;
          b_b->data[i1 + b_b->size[0] * i2].im = 0.0;
          c_loop_ub = a->size[1];
          for (i3 = 0; i3 < c_loop_ub; i3++) {
            y_im = a->data[i1 + a->size[0] * i3].re * b->data[i3 + b->size[0] *
              i2].re - a->data[i1 + a->size[0] * i3].im * b->data[i3 + b->size[0]
              * i2].im;
            d_points = a->data[i1 + a->size[0] * i3].re * b->data[i3 + b->size[0]
              * i2].im + a->data[i1 + a->size[0] * i3].im * b->data[i3 + b->
              size[0] * i2].re;
            b_b->data[i1 + b_b->size[0] * i2].re += y_im;
            b_b->data[i1 + b_b->size[0] * i2].im += d_points;
          }
        }
      }
    } else {
      c_st.site = &i_emlrtRSI;
      if ((a->size[0] == 0) || (a->size[1] == 0) || (b->size[0] == 0) ||
          (b->size[1] == 0)) {
        i1 = b_b->size[0] * b_b->size[1];
        b_b->size[0] = a->size[0];
        b_b->size[1] = b->size[1];
        emxEnsureCapacity_creal_T(&c_st, b_b, i1, &n_emlrtRTEI);
        loop_ub = a->size[0] * b->size[1];
        for (i1 = 0; i1 < loop_ub; i1++) {
          b_b->data[i1] = beta1;
        }
      } else {
        d_st.site = &l_emlrtRSI;
        m_t = (ptrdiff_t)a->size[0];
        n_t = (ptrdiff_t)b->size[1];
        k_t = (ptrdiff_t)a->size[1];
        lda_t = (ptrdiff_t)a->size[0];
        ldb_t = (ptrdiff_t)a->size[1];
        ldc_t = (ptrdiff_t)a->size[0];
        i1 = b_b->size[0] * b_b->size[1];
        b_b->size[0] = a->size[0];
        b_b->size[1] = b->size[1];
        emxEnsureCapacity_creal_T(&d_st, b_b, i1, &g_emlrtRTEI);
        TRANSA = 'N';
        TRANSB = 'N';
        zgemm(&TRANSA, &TRANSB, &m_t, &n_t, &k_t, (real_T *)&alpha1, (real_T *)
              &a->data[0], &lda_t, (real_T *)&b->data[0], &ldb_t, (real_T *)
              &beta1, (real_T *)&b_b->data[0], &ldc_t);
      }
    }

    b_st.site = &d_emlrtRSI;
    c_st.site = &j_emlrtRSI;
    if (b_b->size[1] != C2->size[0]) {
      if (((b_b->size[0] == 1) && (b_b->size[1] == 1)) || ((C2->size[0] == 1) &&
           (C2->size[1] == 1))) {
        emlrtErrorWithMessageIdR2018a(&c_st, &x_emlrtRTEI,
          "Coder:toolbox:mtimes_noDynamicScalarExpansion",
          "Coder:toolbox:mtimes_noDynamicScalarExpansion", 0);
      } else {
        emlrtErrorWithMessageIdR2018a(&c_st, &w_emlrtRTEI, "MATLAB:innerdim",
          "MATLAB:innerdim", 0);
      }
    }

    if ((b_b->size[1] == 1) || (C2->size[0] == 1)) {
      i1 = r0->size[0] * r0->size[1];
      r0->size[0] = b_b->size[0];
      r0->size[1] = C2->size[1];
      emxEnsureCapacity_creal_T(&b_st, r0, i1, &n_emlrtRTEI);
      loop_ub = b_b->size[0];
      for (i1 = 0; i1 < loop_ub; i1++) {
        b_loop_ub = C2->size[1];
        for (i2 = 0; i2 < b_loop_ub; i2++) {
          r0->data[i1 + r0->size[0] * i2].re = 0.0;
          r0->data[i1 + r0->size[0] * i2].im = 0.0;
          c_loop_ub = b_b->size[1];
          for (i3 = 0; i3 < c_loop_ub; i3++) {
            b_re = b_b->data[i1 + b_b->size[0] * i3].re * C2->data[i3 + C2->
              size[0] * i2].re - b_b->data[i1 + b_b->size[0] * i3].im * C2->
              data[i3 + C2->size[0] * i2].im;
            y_re = b_b->data[i1 + b_b->size[0] * i3].re * C2->data[i3 + C2->
              size[0] * i2].im + b_b->data[i1 + b_b->size[0] * i3].im * C2->
              data[i3 + C2->size[0] * i2].re;
            r0->data[i1 + r0->size[0] * i2].re += b_re;
            r0->data[i1 + r0->size[0] * i2].im += y_re;
          }
        }
      }
    } else {
      c_st.site = &i_emlrtRSI;
      if ((b_b->size[0] == 0) || (b_b->size[1] == 0) || (C2->size[0] == 0) ||
          (C2->size[1] == 0)) {
        i1 = r0->size[0] * r0->size[1];
        r0->size[0] = b_b->size[0];
        r0->size[1] = C2->size[1];
        emxEnsureCapacity_creal_T(&c_st, r0, i1, &n_emlrtRTEI);
        loop_ub = b_b->size[0] * C2->size[1];
        for (i1 = 0; i1 < loop_ub; i1++) {
          r0->data[i1] = beta1;
        }
      } else {
        d_st.site = &l_emlrtRSI;
        m_t = (ptrdiff_t)b_b->size[0];
        n_t = (ptrdiff_t)C2->size[1];
        k_t = (ptrdiff_t)b_b->size[1];
        lda_t = (ptrdiff_t)b_b->size[0];
        ldb_t = (ptrdiff_t)b_b->size[1];
        ldc_t = (ptrdiff_t)b_b->size[0];
        i1 = r0->size[0] * r0->size[1];
        r0->size[0] = b_b->size[0];
        r0->size[1] = C2->size[1];
        emxEnsureCapacity_creal_T(&d_st, r0, i1, &g_emlrtRTEI);
        TRANSA = 'N';
        TRANSB = 'N';
        zgemm(&TRANSA, &TRANSB, &m_t, &n_t, &k_t, (real_T *)&alpha1, (real_T *)
              &b_b->data[0], &lda_t, (real_T *)&C2->data[0], &ldb_t, (real_T *)
              &beta1, (real_T *)&r0->data[0], &ldc_t);
      }
    }

    y_re = i0p * 0.0;
    y_im = i0p * 2.0;
    i1 = a->size[0] * a->size[1];
    a->size[0] = rho->size[0];
    a->size[1] = rho->size[1];
    emxEnsureCapacity_creal_T(&st, a, i1, &p_emlrtRTEI);
    loop_ub = rho->size[0] * rho->size[1];
    for (i1 = 0; i1 < loop_ub; i1++) {
      d_points = rho->data[i1].re;
      b_re = rho->data[i1].im;
      a->data[i1].re = y_re * d_points - y_im * b_re;
      a->data[i1].im = y_re * b_re + y_im * d_points;
    }

    b_st.site = &d_emlrtRSI;
    c_st.site = &j_emlrtRSI;
    if (a->size[1] != C3->size[0]) {
      if (((a->size[0] == 1) && (a->size[1] == 1)) || ((C3->size[0] == 1) &&
           (C3->size[1] == 1))) {
        emlrtErrorWithMessageIdR2018a(&c_st, &x_emlrtRTEI,
          "Coder:toolbox:mtimes_noDynamicScalarExpansion",
          "Coder:toolbox:mtimes_noDynamicScalarExpansion", 0);
      } else {
        emlrtErrorWithMessageIdR2018a(&c_st, &w_emlrtRTEI, "MATLAB:innerdim",
          "MATLAB:innerdim", 0);
      }
    }

    if ((a->size[1] == 1) || (C3->size[0] == 1)) {
      i1 = b_b->size[0] * b_b->size[1];
      b_b->size[0] = a->size[0];
      b_b->size[1] = C3->size[1];
      emxEnsureCapacity_creal_T(&b_st, b_b, i1, &p_emlrtRTEI);
      loop_ub = a->size[0];
      for (i1 = 0; i1 < loop_ub; i1++) {
        b_loop_ub = C3->size[1];
        for (i2 = 0; i2 < b_loop_ub; i2++) {
          b_b->data[i1 + b_b->size[0] * i2].re = 0.0;
          b_b->data[i1 + b_b->size[0] * i2].im = 0.0;
          c_loop_ub = a->size[1];
          for (i3 = 0; i3 < c_loop_ub; i3++) {
            y_im = a->data[i1 + a->size[0] * i3].re * C3->data[i3 + C3->size[0] *
              i2].re - a->data[i1 + a->size[0] * i3].im * C3->data[i3 + C3->
              size[0] * i2].im;
            d_points = a->data[i1 + a->size[0] * i3].re * C3->data[i3 + C3->
              size[0] * i2].im + a->data[i1 + a->size[0] * i3].im * C3->data[i3
              + C3->size[0] * i2].re;
            b_b->data[i1 + b_b->size[0] * i2].re += y_im;
            b_b->data[i1 + b_b->size[0] * i2].im += d_points;
          }
        }
      }
    } else {
      c_st.site = &i_emlrtRSI;
      if ((a->size[0] == 0) || (a->size[1] == 0) || (C3->size[0] == 0) ||
          (C3->size[1] == 0)) {
        i1 = b_b->size[0] * b_b->size[1];
        b_b->size[0] = a->size[0];
        b_b->size[1] = C3->size[1];
        emxEnsureCapacity_creal_T(&c_st, b_b, i1, &p_emlrtRTEI);
        loop_ub = a->size[0] * C3->size[1];
        for (i1 = 0; i1 < loop_ub; i1++) {
          b_b->data[i1] = beta1;
        }
      } else {
        d_st.site = &l_emlrtRSI;
        m_t = (ptrdiff_t)a->size[0];
        n_t = (ptrdiff_t)C3->size[1];
        k_t = (ptrdiff_t)a->size[1];
        lda_t = (ptrdiff_t)a->size[0];
        ldb_t = (ptrdiff_t)a->size[1];
        ldc_t = (ptrdiff_t)a->size[0];
        i1 = b_b->size[0] * b_b->size[1];
        b_b->size[0] = a->size[0];
        b_b->size[1] = C3->size[1];
        emxEnsureCapacity_creal_T(&d_st, b_b, i1, &g_emlrtRTEI);
        TRANSA = 'N';
        TRANSB = 'N';
        zgemm(&TRANSA, &TRANSB, &m_t, &n_t, &k_t, (real_T *)&alpha1, (real_T *)
              &a->data[0], &lda_t, (real_T *)&C3->data[0], &ldb_t, (real_T *)
              &beta1, (real_T *)&b_b->data[0], &ldc_t);
      }
    }

    b_st.site = &d_emlrtRSI;
    i1 = b->size[0] * b->size[1];
    b->size[0] = r3->size[0];
    b->size[1] = r3->size[1];
    emxEnsureCapacity_creal_T(&b_st, b, i1, &q_emlrtRTEI);
    loop_ub = r3->size[0] * r3->size[1];
    for (i1 = 0; i1 < loop_ub; i1++) {
      b->data[i1].re = r3->data[i1] * C4->data[i1].re;
      b->data[i1].im = r3->data[i1] * C4->data[i1].im;
    }

    c_st.site = &j_emlrtRSI;
    if (b_b->size[1] != b->size[0]) {
      if (((b_b->size[0] == 1) && (b_b->size[1] == 1)) || ((b->size[0] == 1) &&
           (b->size[1] == 1))) {
        emlrtErrorWithMessageIdR2018a(&c_st, &x_emlrtRTEI,
          "Coder:toolbox:mtimes_noDynamicScalarExpansion",
          "Coder:toolbox:mtimes_noDynamicScalarExpansion", 0);
      } else {
        emlrtErrorWithMessageIdR2018a(&c_st, &w_emlrtRTEI, "MATLAB:innerdim",
          "MATLAB:innerdim", 0);
      }
    }

    if ((b_b->size[1] == 1) || (b->size[0] == 1)) {
      i1 = r1->size[0] * r1->size[1];
      r1->size[0] = b_b->size[0];
      r1->size[1] = b->size[1];
      emxEnsureCapacity_creal_T(&b_st, r1, i1, &p_emlrtRTEI);
      loop_ub = b_b->size[0];
      for (i1 = 0; i1 < loop_ub; i1++) {
        b_loop_ub = b->size[1];
        for (i2 = 0; i2 < b_loop_ub; i2++) {
          r1->data[i1 + r1->size[0] * i2].re = 0.0;
          r1->data[i1 + r1->size[0] * i2].im = 0.0;
          c_loop_ub = b_b->size[1];
          for (i3 = 0; i3 < c_loop_ub; i3++) {
            b_re = b_b->data[i1 + b_b->size[0] * i3].re * b->data[i3 + b->size[0]
              * i2].re - b_b->data[i1 + b_b->size[0] * i3].im * b->data[i3 +
              b->size[0] * i2].im;
            y_re = b_b->data[i1 + b_b->size[0] * i3].re * b->data[i3 + b->size[0]
              * i2].im + b_b->data[i1 + b_b->size[0] * i3].im * b->data[i3 +
              b->size[0] * i2].re;
            r1->data[i1 + r1->size[0] * i2].re += b_re;
            r1->data[i1 + r1->size[0] * i2].im += y_re;
          }
        }
      }
    } else {
      c_st.site = &i_emlrtRSI;
      if ((b_b->size[0] == 0) || (b_b->size[1] == 0) || (b->size[0] == 0) ||
          (b->size[1] == 0)) {
        i1 = r1->size[0] * r1->size[1];
        r1->size[0] = b_b->size[0];
        r1->size[1] = b->size[1];
        emxEnsureCapacity_creal_T(&c_st, r1, i1, &p_emlrtRTEI);
        loop_ub = b_b->size[0] * b->size[1];
        for (i1 = 0; i1 < loop_ub; i1++) {
          r1->data[i1] = beta1;
        }
      } else {
        d_st.site = &l_emlrtRSI;
        m_t = (ptrdiff_t)b_b->size[0];
        n_t = (ptrdiff_t)b->size[1];
        k_t = (ptrdiff_t)b_b->size[1];
        lda_t = (ptrdiff_t)b_b->size[0];
        ldb_t = (ptrdiff_t)b_b->size[1];
        ldc_t = (ptrdiff_t)b_b->size[0];
        i1 = r1->size[0] * r1->size[1];
        r1->size[0] = b_b->size[0];
        r1->size[1] = b->size[1];
        emxEnsureCapacity_creal_T(&d_st, r1, i1, &g_emlrtRTEI);
        TRANSA = 'N';
        TRANSB = 'N';
        zgemm(&TRANSA, &TRANSB, &m_t, &n_t, &k_t, (real_T *)&alpha1, (real_T *)
              &b_b->data[0], &lda_t, (real_T *)&b->data[0], &ldb_t, (real_T *)
              &beta1, (real_T *)&r1->data[0], &ldc_t);
      }
    }

    iv0[0] = r0->size[0];
    iv0[1] = r0->size[1];
    b_E_pp[0] = r1->size[0];
    b_E_pp[1] = r1->size[1];
    if ((iv0[0] != b_E_pp[0]) || (iv0[1] != b_E_pp[1])) {
      emlrtSizeEqCheckNDR2012b(&iv0[0], &b_E_pp[0], &j_emlrtECI, &st);
    }

    i1 = c_points->size[0] * c_points->size[1];
    c_points->size[0] = r0->size[0];
    c_points->size[1] = r0->size[1];
    emxEnsureCapacity_creal_T(&st, c_points, i1, &n_emlrtRTEI);
    loop_ub = r0->size[0] * r0->size[1];
    for (i1 = 0; i1 < loop_ub; i1++) {
      c_points->data[i1].re = r0->data[i1].re + r1->data[i1].re;
      c_points->data[i1].im = r0->data[i1].im + r1->data[i1].im;
    }

    b_st.site = &d_emlrtRSI;
    r4 = &b_st;
    i1 = GK->size[0];
    i2 = 1 + k;
    if ((i2 < 1) || (i2 > i1)) {
      emlrtDynamicBoundsCheckR2012b(i2, 1, i1, &f_emlrtBCI, &b_st);
    }

    GK->data[i2 - 1] = trace(r4, c_points);
    if (*emlrtBreakCheckR2012bFlagVar != 0) {
      emlrtBreakCheckR2012b(&b_st);
    }
  }

  emxFree_creal_T(&c_points);
  emxFree_real_T(&b_points);
  emxFree_real_T(&r3);
  emxFree_real_T(&r2);
  emxFree_creal_T(&b_b);
  emxFree_creal_T(&a);
  emxFree_creal_T(&b);
  emxFree_creal_T(&r1);
  emxFree_creal_T(&r0);
  emxFree_real_T(&E_qq);
  emxFree_real_T(&E_mm);
  emxFree_real_T(&E_nn);
  emxFree_real_T(&E_pp);

  /*   */
  /*  GR = mmat(mmat(rho, bsxfun(@times, (1./(permute(points,[2,3,1]) - (E_p' - E_n) + i0p*1i)),C1)), C2) + ... */
  /*      mmat(mmat(rho, C3), bsxfun(@times, (1./(permute(points, [2,3,1]) - (E_n' - E_q) + i0p*1i)), C4)); */
  /*  n = size(GR,3); */
  /*  [GR_real, GR_imag, GK_real, GK_imag] = deal(zeros(n,1)); */
  /*  for k = 1:n */
  /*      GR_real(k,1) = real(trace(GR(:,:,k))); */
  /*      GR_imag(k,1) = imag(trace(GR(:,:,k))); */
  /*  end */
  /*   */
  /*  GK = -2i*i0p*   mmat(rho, mmat(bsxfun(@times, (1./((permute(points,[2,3,1]) - (E_p' - E_n)).^2 + i0p^2)), C1), C2)) + ... */
  /*      2i*i0p*   mmat(rho, mmat(C3, bsxfun(@times, (1./((permute(points,[2,3,1]) - (E_n' - E_q)).^2 + i0p^2)), C4))); */
  /*   */
  /*   */
  /*  for k = 1:n */
  /*      GK_real(k,1) = real(trace(GK(:,:,k))); */
  /*      GK_imag(k,1) = imag(trace(GK(:,:,k))); */
  /*  end */
  /*   */
  /*   */
  /*   */
  /*  end */
  /*   */
  /*  function C = mmat(A,B) */
  /*   */
  /*  na = ndims(A); */
  /*  nb = ndims(B); */
  /*  C = zeros(size(A,1), size(B,2), max(size(A,3), size(B,3))); */
  /*  if na == nb */
  /*      %equal size case */
  /*      if na == 2 */
  /*          %matrix multiplication */
  /*          C = A*B; */
  /*      elseif na == 3 */
  /*          n = size(A,3); */
  /*          for k = 1:n */
  /*              C(:,:,k) = A(:,:,k) * B(:,:,k); */
  /*          end */
  /*      end */
  /*  else */
  /*      if na == 3 && nb == 2 */
  /*          n = size(A,3); */
  /*          for k = 1:n */
  /*              C(:,:,k) = A(:,:,k) * B(:,:); */
  /*          end */
  /*      elseif na == 2 && nb == 3 */
  /*          n = size(B,3); */
  /*           */
  /*          for k = 1:n */
  /*              C(:,:,k) = A(:,:) * B(:,:,k); */
  /*          end */
  /*      end */
  /*  end */
  /*   */
  /*  end */
  /*   */
  emlrtHeapReferenceStackLeaveFcnR2012b(sp);
}

/* End of code generation (m_GF.c) */
