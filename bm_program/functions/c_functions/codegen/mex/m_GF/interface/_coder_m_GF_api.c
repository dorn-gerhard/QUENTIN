/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * _coder_m_GF_api.c
 *
 * Code generation for function '_coder_m_GF_api'
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "m_GF.h"
#include "_coder_m_GF_api.h"
#include "m_GF_emxutil.h"
#include "m_GF_data.h"

/* Variable Definitions */
static emlrtRTEInfo v_emlrtRTEI = { 1, /* lineNo */
  1,                                   /* colNo */
  "_coder_m_GF_api",                   /* fName */
  ""                                   /* pName */
};

/* Function Declarations */
static void b_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u, const
  emlrtMsgIdentifier *parentId, emxArray_real_T *y);
static real_T c_emlrt_marshallIn(const emlrtStack *sp, const mxArray *i0p, const
  char_T *identifier);
static real_T d_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u, const
  emlrtMsgIdentifier *parentId);
static void e_emlrt_marshallIn(const emlrtStack *sp, const mxArray *E_p, const
  char_T *identifier, emxArray_real_T *y);
static void emlrt_marshallIn(const emlrtStack *sp, const mxArray *points, const
  char_T *identifier, emxArray_real_T *y);
static const mxArray *emlrt_marshallOut(const emlrtStack *sp, const
  emxArray_creal_T *u);
static void f_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u, const
  emlrtMsgIdentifier *parentId, emxArray_real_T *y);
static void g_emlrt_marshallIn(const emlrtStack *sp, const mxArray *C1, const
  char_T *identifier, emxArray_creal_T *y);
static void h_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u, const
  emlrtMsgIdentifier *parentId, emxArray_creal_T *y);
static void i_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src, const
  emlrtMsgIdentifier *msgId, emxArray_real_T *ret);
static real_T j_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src, const
  emlrtMsgIdentifier *msgId);
static void k_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src, const
  emlrtMsgIdentifier *msgId, emxArray_real_T *ret);
static void l_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src, const
  emlrtMsgIdentifier *msgId, emxArray_creal_T *ret);

/* Function Definitions */
static void b_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u, const
  emlrtMsgIdentifier *parentId, emxArray_real_T *y)
{
  i_emlrt_marshallIn(sp, emlrtAlias(u), parentId, y);
  emlrtDestroyArray(&u);
}

static real_T c_emlrt_marshallIn(const emlrtStack *sp, const mxArray *i0p, const
  char_T *identifier)
{
  real_T y;
  emlrtMsgIdentifier thisId;
  thisId.fIdentifier = (const char *)identifier;
  thisId.fParent = NULL;
  thisId.bParentIsCell = false;
  y = d_emlrt_marshallIn(sp, emlrtAlias(i0p), &thisId);
  emlrtDestroyArray(&i0p);
  return y;
}

static real_T d_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u, const
  emlrtMsgIdentifier *parentId)
{
  real_T y;
  y = j_emlrt_marshallIn(sp, emlrtAlias(u), parentId);
  emlrtDestroyArray(&u);
  return y;
}

static void e_emlrt_marshallIn(const emlrtStack *sp, const mxArray *E_p, const
  char_T *identifier, emxArray_real_T *y)
{
  emlrtMsgIdentifier thisId;
  thisId.fIdentifier = (const char *)identifier;
  thisId.fParent = NULL;
  thisId.bParentIsCell = false;
  f_emlrt_marshallIn(sp, emlrtAlias(E_p), &thisId, y);
  emlrtDestroyArray(&E_p);
}

static void emlrt_marshallIn(const emlrtStack *sp, const mxArray *points, const
  char_T *identifier, emxArray_real_T *y)
{
  emlrtMsgIdentifier thisId;
  thisId.fIdentifier = (const char *)identifier;
  thisId.fParent = NULL;
  thisId.bParentIsCell = false;
  b_emlrt_marshallIn(sp, emlrtAlias(points), &thisId, y);
  emlrtDestroyArray(&points);
}

static const mxArray *emlrt_marshallOut(const emlrtStack *sp, const
  emxArray_creal_T *u)
{
  const mxArray *y;
  const mxArray *m0;
  y = NULL;
  m0 = emlrtCreateNumericArray(1, *(int32_T (*)[1])u->size, mxDOUBLE_CLASS,
    mxCOMPLEX);
  emlrtExportNumericArrayR2013b(sp, m0, &u->data[0], 8);
  emlrtAssign(&y, m0);
  return y;
}

static void f_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u, const
  emlrtMsgIdentifier *parentId, emxArray_real_T *y)
{
  k_emlrt_marshallIn(sp, emlrtAlias(u), parentId, y);
  emlrtDestroyArray(&u);
}

static void g_emlrt_marshallIn(const emlrtStack *sp, const mxArray *C1, const
  char_T *identifier, emxArray_creal_T *y)
{
  emlrtMsgIdentifier thisId;
  thisId.fIdentifier = (const char *)identifier;
  thisId.fParent = NULL;
  thisId.bParentIsCell = false;
  h_emlrt_marshallIn(sp, emlrtAlias(C1), &thisId, y);
  emlrtDestroyArray(&C1);
}

static void h_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u, const
  emlrtMsgIdentifier *parentId, emxArray_creal_T *y)
{
  l_emlrt_marshallIn(sp, emlrtAlias(u), parentId, y);
  emlrtDestroyArray(&u);
}

static void i_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src, const
  emlrtMsgIdentifier *msgId, emxArray_real_T *ret)
{
  static const int32_T dims[1] = { -1 };

  const boolean_T bv0[1] = { true };

  int32_T iv1[1];
  int32_T i7;
  emlrtCheckVsBuiltInR2012b(sp, msgId, src, "double", false, 1U, dims, &bv0[0],
    iv1);
  ret->allocatedSize = iv1[0];
  i7 = ret->size[0];
  ret->size[0] = iv1[0];
  emxEnsureCapacity_real_T(sp, ret, i7, (emlrtRTEInfo *)NULL);
  ret->data = (real_T *)emlrtMxGetData(src);
  ret->canFreeData = false;
  ret->canFreeData = false;
  emlrtDestroyArray(&src);
}

static real_T j_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src, const
  emlrtMsgIdentifier *msgId)
{
  real_T ret;
  static const int32_T dims = 0;
  emlrtCheckBuiltInR2012b(sp, msgId, src, "double", false, 0U, &dims);
  ret = *(real_T *)emlrtMxGetData(src);
  emlrtDestroyArray(&src);
  return ret;
}

static void k_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src, const
  emlrtMsgIdentifier *msgId, emxArray_real_T *ret)
{
  static const int32_T dims[2] = { -1, -1 };

  const boolean_T bv1[2] = { true, true };

  int32_T iv2[2];
  int32_T i8;
  emlrtCheckVsBuiltInR2012b(sp, msgId, src, "double", false, 2U, dims, &bv1[0],
    iv2);
  ret->allocatedSize = iv2[0] * iv2[1];
  i8 = ret->size[0] * ret->size[1];
  ret->size[0] = iv2[0];
  ret->size[1] = iv2[1];
  emxEnsureCapacity_real_T(sp, ret, i8, (emlrtRTEInfo *)NULL);
  ret->data = (real_T *)emlrtMxGetData(src);
  ret->canFreeData = false;
  ret->canFreeData = false;
  emlrtDestroyArray(&src);
}

static void l_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src, const
  emlrtMsgIdentifier *msgId, emxArray_creal_T *ret)
{
  static const int32_T dims[2] = { -1, -1 };

  const boolean_T bv2[2] = { true, true };

  int32_T iv3[2];
  int32_T i9;
  emlrtCheckVsBuiltInR2012b(sp, msgId, src, "double", true, 2U, dims, &bv2[0],
    iv3);
  i9 = ret->size[0] * ret->size[1];
  ret->size[0] = iv3[0];
  ret->size[1] = iv3[1];
  emxEnsureCapacity_creal_T(sp, ret, i9, (emlrtRTEInfo *)NULL);
  emlrtImportArrayR2015b(sp, src, ret->data, 8, true);
  emlrtDestroyArray(&src);
}

void m_GF_api(const mxArray * const prhs[10], int32_T nlhs, const mxArray *plhs
              [2])
{
  emxArray_real_T *points;
  emxArray_real_T *E_p;
  emxArray_real_T *E_n;
  emxArray_real_T *E_q;
  emxArray_creal_T *C1;
  emxArray_creal_T *C2;
  emxArray_creal_T *C3;
  emxArray_creal_T *C4;
  emxArray_creal_T *rho;
  emxArray_creal_T *GR;
  emxArray_creal_T *GK;
  real_T i0p;
  emlrtStack st = { NULL,              /* site */
    NULL,                              /* tls */
    NULL                               /* prev */
  };

  st.tls = emlrtRootTLSGlobal;
  emlrtHeapReferenceStackEnterFcnR2012b(&st);
  emxInit_real_T(&st, &points, 1, &v_emlrtRTEI, true);
  emxInit_real_T(&st, &E_p, 2, &v_emlrtRTEI, true);
  emxInit_real_T(&st, &E_n, 2, &v_emlrtRTEI, true);
  emxInit_real_T(&st, &E_q, 2, &v_emlrtRTEI, true);
  emxInit_creal_T(&st, &C1, 2, &v_emlrtRTEI, true);
  emxInit_creal_T(&st, &C2, 2, &v_emlrtRTEI, true);
  emxInit_creal_T(&st, &C3, 2, &v_emlrtRTEI, true);
  emxInit_creal_T(&st, &C4, 2, &v_emlrtRTEI, true);
  emxInit_creal_T(&st, &rho, 2, &v_emlrtRTEI, true);
  emxInit_creal_T(&st, &GR, 1, &v_emlrtRTEI, true);
  emxInit_creal_T(&st, &GK, 1, &v_emlrtRTEI, true);

  /* Marshall function inputs */
  points->canFreeData = false;
  emlrt_marshallIn(&st, emlrtAlias(prhs[0]), "points", points);
  i0p = c_emlrt_marshallIn(&st, emlrtAliasP(prhs[1]), "i0p");
  E_p->canFreeData = false;
  e_emlrt_marshallIn(&st, emlrtAlias(prhs[2]), "E_p", E_p);
  E_n->canFreeData = false;
  e_emlrt_marshallIn(&st, emlrtAlias(prhs[3]), "E_n", E_n);
  E_q->canFreeData = false;
  e_emlrt_marshallIn(&st, emlrtAlias(prhs[4]), "E_q", E_q);
  g_emlrt_marshallIn(&st, emlrtAliasP(prhs[5]), "C1", C1);
  g_emlrt_marshallIn(&st, emlrtAliasP(prhs[6]), "C2", C2);
  g_emlrt_marshallIn(&st, emlrtAliasP(prhs[7]), "C3", C3);
  g_emlrt_marshallIn(&st, emlrtAliasP(prhs[8]), "C4", C4);
  g_emlrt_marshallIn(&st, emlrtAliasP(prhs[9]), "rho", rho);

  /* Invoke the target function */
  m_GF(&st, points, i0p, E_p, E_n, E_q, C1, C2, C3, C4, rho, GR, GK);

  /* Marshall function outputs */
  plhs[0] = emlrt_marshallOut(&st, GR);
  emxFree_creal_T(&GR);
  emxFree_creal_T(&rho);
  emxFree_creal_T(&C4);
  emxFree_creal_T(&C3);
  emxFree_creal_T(&C2);
  emxFree_creal_T(&C1);
  emxFree_real_T(&E_q);
  emxFree_real_T(&E_n);
  emxFree_real_T(&E_p);
  emxFree_real_T(&points);
  if (nlhs > 1) {
    plhs[1] = emlrt_marshallOut(&st, GK);
  }

  emxFree_creal_T(&GK);
  emlrtHeapReferenceStackLeaveFcnR2012b(&st);
}

/* End of code generation (_coder_m_GF_api.c) */
