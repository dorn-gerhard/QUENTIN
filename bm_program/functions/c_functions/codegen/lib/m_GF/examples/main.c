/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * main.c
 *
 * Code generation for function 'main'
 *
 */

/*************************************************************************/
/* This automatically generated example C main file shows how to call    */
/* entry-point functions that MATLAB Coder generated. You must customize */
/* this file for your application. Do not modify this file directly.     */
/* Instead, make a copy of this file, modify it, and integrate it into   */
/* your development environment.                                         */
/*                                                                       */
/* This file initializes entry-point function arguments to a default     */
/* size and value before calling the entry-point functions. It does      */
/* not store or use any values returned from the entry-point functions.  */
/* If necessary, it does pre-allocate memory for returned values.        */
/* You can use this file as a starting point for a main function that    */
/* you can deploy in your application.                                   */
/*                                                                       */
/* After you copy the file, and before you deploy it, you must make the  */
/* following changes:                                                    */
/* * For variable-size function arguments, change the example sizes to   */
/* the sizes that your application requires.                             */
/* * Change the example values of function arguments to the values that  */
/* your application requires.                                            */
/* * If the entry-point functions return values, store these values or   */
/* otherwise use them as required by your application.                   */
/*                                                                       */
/*************************************************************************/
/* Include files */
#include "m_GF.h"
#include "main.h"
#include "m_GF_terminate.h"
#include "m_GF_emxAPI.h"
#include "m_GF_initialize.h"

/* Function Declarations */
static emxArray_real_T *argInit_Unboundedx1_real_T(void);
static creal_T argInit_creal_T(void);
static double argInit_real_T(void);
static emxArray_creal_T *c_argInit_UnboundedxUnbounded_c(void);
static emxArray_real_T *c_argInit_UnboundedxUnbounded_r(void);
static void main_m_GF(void);

/* Function Definitions */
static emxArray_real_T *argInit_Unboundedx1_real_T(void)
{
  emxArray_real_T *result;
  static int iv0[1] = { 2 };

  int idx0;

  /* Set the size of the array.
     Change this size to the value that the application requires. */
  result = emxCreateND_real_T(1, iv0);

  /* Loop over the array to initialize each element. */
  for (idx0 = 0; idx0 < result->size[0U]; idx0++) {
    /* Set the value of the array element.
       Change this value to the value that the application requires. */
    result->data[idx0] = argInit_real_T();
  }

  return result;
}

static creal_T argInit_creal_T(void)
{
  creal_T result;
  double result_tmp;

  /* Set the value of the complex variable.
     Change this value to the value that the application requires. */
  result_tmp = argInit_real_T();
  result.re = result_tmp;
  result.im = result_tmp;
  return result;
}

static double argInit_real_T(void)
{
  return 0.0;
}

static emxArray_creal_T *c_argInit_UnboundedxUnbounded_c(void)
{
  emxArray_creal_T *result;
  int idx0;
  int idx1;

  /* Set the size of the array.
     Change this size to the value that the application requires. */
  result = emxCreate_creal_T(2, 2);

  /* Loop over the array to initialize each element. */
  for (idx0 = 0; idx0 < result->size[0U]; idx0++) {
    for (idx1 = 0; idx1 < result->size[1U]; idx1++) {
      /* Set the value of the array element.
         Change this value to the value that the application requires. */
      result->data[idx0 + result->size[0] * idx1] = argInit_creal_T();
    }
  }

  return result;
}

static emxArray_real_T *c_argInit_UnboundedxUnbounded_r(void)
{
  emxArray_real_T *result;
  int idx0;
  int idx1;

  /* Set the size of the array.
     Change this size to the value that the application requires. */
  result = emxCreate_real_T(2, 2);

  /* Loop over the array to initialize each element. */
  for (idx0 = 0; idx0 < result->size[0U]; idx0++) {
    for (idx1 = 0; idx1 < result->size[1U]; idx1++) {
      /* Set the value of the array element.
         Change this value to the value that the application requires. */
      result->data[idx0 + result->size[0] * idx1] = argInit_real_T();
    }
  }

  return result;
}

static void main_m_GF(void)
{
  emxArray_creal_T *GR;
  emxArray_creal_T *GK;
  emxArray_real_T *points;
  double i0p;
  emxArray_real_T *E_p;
  emxArray_real_T *E_n;
  emxArray_real_T *E_q;
  emxArray_creal_T *C1;
  emxArray_creal_T *C2;
  emxArray_creal_T *C3;
  emxArray_creal_T *C4;
  emxArray_creal_T *rho;
  emxInitArray_creal_T(&GR, 1);
  emxInitArray_creal_T(&GK, 1);

  /* Initialize function 'm_GF' input arguments. */
  /* Initialize function input argument 'points'. */
  points = argInit_Unboundedx1_real_T();
  i0p = argInit_real_T();

  /* Initialize function input argument 'E_p'. */
  E_p = c_argInit_UnboundedxUnbounded_r();

  /* Initialize function input argument 'E_n'. */
  E_n = c_argInit_UnboundedxUnbounded_r();

  /* Initialize function input argument 'E_q'. */
  E_q = c_argInit_UnboundedxUnbounded_r();

  /* Initialize function input argument 'C1'. */
  C1 = c_argInit_UnboundedxUnbounded_c();

  /* Initialize function input argument 'C2'. */
  C2 = c_argInit_UnboundedxUnbounded_c();

  /* Initialize function input argument 'C3'. */
  C3 = c_argInit_UnboundedxUnbounded_c();

  /* Initialize function input argument 'C4'. */
  C4 = c_argInit_UnboundedxUnbounded_c();

  /* Initialize function input argument 'rho'. */
  rho = c_argInit_UnboundedxUnbounded_c();

  /* Call the entry-point 'm_GF'. */
  m_GF(points, i0p, E_p, E_n, E_q, C1, C2, C3, C4, rho, GR, GK);
  emxDestroyArray_creal_T(GK);
  emxDestroyArray_creal_T(GR);
  emxDestroyArray_creal_T(rho);
  emxDestroyArray_creal_T(C4);
  emxDestroyArray_creal_T(C3);
  emxDestroyArray_creal_T(C2);
  emxDestroyArray_creal_T(C1);
  emxDestroyArray_real_T(E_q);
  emxDestroyArray_real_T(E_n);
  emxDestroyArray_real_T(E_p);
  emxDestroyArray_real_T(points);
}

int main(int argc, const char * const argv[])
{
  (void)argc;
  (void)argv;

  /* Initialize the application.
     You do not need to do this more than one time. */
  m_GF_initialize();

  /* Invoke the entry-point functions.
     You can call entry-point functions multiple times. */
  main_m_GF();

  /* Terminate the application.
     You do not need to do this more than one time. */
  m_GF_terminate();
  return 0;
}

/* End of code generation (main.c) */
