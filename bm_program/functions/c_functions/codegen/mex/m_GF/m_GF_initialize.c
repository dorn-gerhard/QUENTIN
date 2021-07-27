/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * m_GF_initialize.c
 *
 * Code generation for function 'm_GF_initialize'
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "m_GF.h"
#include "m_GF_initialize.h"
#include "_coder_m_GF_mex.h"
#include "m_GF_data.h"

/* Function Definitions */
void m_GF_initialize(void)
{
  emlrtStack st = { NULL,              /* site */
    NULL,                              /* tls */
    NULL                               /* prev */
  };

  mexFunctionCreateRootTLS();
  emlrtBreakCheckR2012bFlagVar = emlrtGetBreakCheckFlagAddressR2012b();
  st.tls = emlrtRootTLSGlobal;
  emlrtClearAllocCountR2012b(&st, false, 0U, 0);
  emlrtEnterRtStackR2012b(&st);
  emlrtFirstTimeR2012b(emlrtRootTLSGlobal);
}

/* End of code generation (m_GF_initialize.c) */
