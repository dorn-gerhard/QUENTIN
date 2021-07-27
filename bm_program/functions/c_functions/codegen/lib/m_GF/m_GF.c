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
#include "m_GF.h"
#include "m_GF_emxutil.h"
#include "rdivide_helper.h"
#include "power.h"
#include "meshgrid.h"

/* Function Definitions */
void m_GF(const emxArray_real_T *points, double i0p, const emxArray_real_T *E_p,
          const emxArray_real_T *E_n, const emxArray_real_T *E_q, const
          emxArray_creal_T *C1, const emxArray_creal_T *C2, const
          emxArray_creal_T *C3, const emxArray_creal_T *C4, const
          emxArray_creal_T *rho, emxArray_creal_T *GR, emxArray_creal_T *GK)
{
  emxArray_real_T *E_pp;
  emxArray_real_T *E_nn;
  emxArray_real_T *E_mm;
  emxArray_real_T *E_qq;
  int i0;
  int loop_ub;
  emxArray_creal_T *r0;
  emxArray_creal_T *y;
  emxArray_creal_T *b;
  emxArray_creal_T *a;
  emxArray_real_T *r1;
  emxArray_real_T *b_points;
  emxArray_creal_T *c_points;
  int k;
  double temp_re;
  double y_re;
  int i1;
  double d_points;
  int m;
  double temp_im;
  int inner;
  int n;
  double C1_im;
  int j;
  int coffset;
  int boffset;
  int i;
  int b_k;
  int aoffset;
  double c_tmp;
  emxInit_real_T(&E_pp, 2);
  emxInit_real_T(&E_nn, 2);
  emxInit_real_T(&E_mm, 2);
  emxInit_real_T(&E_qq, 2);
  meshgrid(E_p, E_n, E_pp, E_nn);
  meshgrid(E_n, E_q, E_mm, E_qq);
  i0 = GR->size[0];
  GR->size[0] = points->size[0];
  emxEnsureCapacity_creal_T(GR, i0);
  loop_ub = points->size[0];
  for (i0 = 0; i0 < loop_ub; i0++) {
    GR->data[i0].re = 0.0;
    GR->data[i0].im = 0.0;
  }

  i0 = GK->size[0];
  GK->size[0] = points->size[0];
  emxEnsureCapacity_creal_T(GK, i0);
  loop_ub = points->size[0];
  for (i0 = 0; i0 < loop_ub; i0++) {
    GK->data[i0].re = 0.0;
    GK->data[i0].im = 0.0;
  }

  i0 = points->size[0];
  emxInit_creal_T(&r0, 2);
  emxInit_creal_T(&y, 2);
  emxInit_creal_T(&b, 2);
  emxInit_creal_T(&a, 2);
  emxInit_real_T(&r1, 2);
  emxInit_real_T(&b_points, 2);
  emxInit_creal_T(&c_points, 2);
  for (k = 0; k < i0; k++) {
    temp_re = i0p * 0.0;
    y_re = i0p * 0.0;
    i1 = c_points->size[0] * c_points->size[1];
    c_points->size[0] = E_pp->size[0];
    c_points->size[1] = E_pp->size[1];
    emxEnsureCapacity_creal_T(c_points, i1);
    d_points = points->data[k];
    loop_ub = E_pp->size[0] * E_pp->size[1];
    for (i1 = 0; i1 < loop_ub; i1++) {
      c_points->data[i1].re = (d_points - (E_pp->data[i1] - E_nn->data[i1])) +
        temp_re;
      c_points->data[i1].im = i0p;
    }

    rdivide_helper(c_points, b);
    i1 = b->size[0] * b->size[1];
    m = b->size[0] * b->size[1];
    emxEnsureCapacity_creal_T(b, m);
    loop_ub = i1 - 1;
    for (i1 = 0; i1 <= loop_ub; i1++) {
      temp_im = b->data[i1].re;
      temp_re = b->data[i1].im;
      d_points = C1->data[i1].re;
      C1_im = C1->data[i1].im;
      b->data[i1].re = temp_im * d_points - temp_re * C1_im;
      b->data[i1].im = temp_im * C1_im + temp_re * d_points;
    }

    if ((rho->size[1] == 1) || (b->size[0] == 1)) {
      i1 = y->size[0] * y->size[1];
      y->size[0] = rho->size[0];
      y->size[1] = b->size[1];
      emxEnsureCapacity_creal_T(y, i1);
      loop_ub = rho->size[0];
      for (i1 = 0; i1 < loop_ub; i1++) {
        inner = b->size[1];
        for (m = 0; m < inner; m++) {
          y->data[i1 + y->size[0] * m].re = 0.0;
          y->data[i1 + y->size[0] * m].im = 0.0;
          n = rho->size[1];
          for (j = 0; j < n; j++) {
            d_points = rho->data[i1 + rho->size[0] * j].re * b->data[j + b->
              size[0] * m].re - rho->data[i1 + rho->size[0] * j].im * b->data[j
              + b->size[0] * m].im;
            C1_im = rho->data[i1 + rho->size[0] * j].re * b->data[j + b->size[0]
              * m].im + rho->data[i1 + rho->size[0] * j].im * b->data[j +
              b->size[0] * m].re;
            y->data[i1 + y->size[0] * m].re += d_points;
            y->data[i1 + y->size[0] * m].im += C1_im;
          }
        }
      }
    } else {
      m = rho->size[0];
      inner = rho->size[1];
      n = b->size[1];
      i1 = y->size[0] * y->size[1];
      y->size[0] = rho->size[0];
      y->size[1] = b->size[1];
      emxEnsureCapacity_creal_T(y, i1);
      for (j = 0; j < n; j++) {
        coffset = j * m;
        boffset = j * inner;
        for (i = 0; i < m; i++) {
          i1 = coffset + i;
          y->data[i1].re = 0.0;
          y->data[i1].im = 0.0;
        }

        for (b_k = 0; b_k < inner; b_k++) {
          aoffset = b_k * m;
          loop_ub = boffset + b_k;
          temp_re = b->data[loop_ub].re;
          temp_im = b->data[loop_ub].im;
          for (i = 0; i < m; i++) {
            loop_ub = aoffset + i;
            d_points = temp_re * rho->data[loop_ub].re - temp_im * rho->
              data[loop_ub].im;
            C1_im = temp_re * rho->data[loop_ub].im + temp_im * rho->
              data[loop_ub].re;
            i1 = coffset + i;
            y->data[i1].re += d_points;
            y->data[i1].im += C1_im;
          }
        }
      }
    }

    if ((y->size[1] == 1) || (C2->size[0] == 1)) {
      i1 = a->size[0] * a->size[1];
      a->size[0] = y->size[0];
      a->size[1] = C2->size[1];
      emxEnsureCapacity_creal_T(a, i1);
      loop_ub = y->size[0];
      for (i1 = 0; i1 < loop_ub; i1++) {
        inner = C2->size[1];
        for (m = 0; m < inner; m++) {
          a->data[i1 + a->size[0] * m].re = 0.0;
          a->data[i1 + a->size[0] * m].im = 0.0;
          n = y->size[1];
          for (j = 0; j < n; j++) {
            d_points = y->data[i1 + y->size[0] * j].re * C2->data[j + C2->size[0]
              * m].re - y->data[i1 + y->size[0] * j].im * C2->data[j + C2->size
              [0] * m].im;
            C1_im = y->data[i1 + y->size[0] * j].re * C2->data[j + C2->size[0] *
              m].im + y->data[i1 + y->size[0] * j].im * C2->data[j + C2->size[0]
              * m].re;
            a->data[i1 + a->size[0] * m].re += d_points;
            a->data[i1 + a->size[0] * m].im += C1_im;
          }
        }
      }
    } else {
      m = y->size[0];
      inner = y->size[1];
      n = C2->size[1];
      i1 = a->size[0] * a->size[1];
      a->size[0] = y->size[0];
      a->size[1] = C2->size[1];
      emxEnsureCapacity_creal_T(a, i1);
      for (j = 0; j < n; j++) {
        coffset = j * m;
        boffset = j * inner;
        for (i = 0; i < m; i++) {
          i1 = coffset + i;
          a->data[i1].re = 0.0;
          a->data[i1].im = 0.0;
        }

        for (b_k = 0; b_k < inner; b_k++) {
          aoffset = b_k * m;
          loop_ub = boffset + b_k;
          temp_re = C2->data[loop_ub].re;
          temp_im = C2->data[loop_ub].im;
          for (i = 0; i < m; i++) {
            loop_ub = aoffset + i;
            d_points = temp_re * y->data[loop_ub].re - temp_im * y->data[loop_ub]
              .im;
            C1_im = temp_re * y->data[loop_ub].im + temp_im * y->data[loop_ub].
              re;
            i1 = coffset + i;
            a->data[i1].re += d_points;
            a->data[i1].im += C1_im;
          }
        }
      }
    }

    if ((rho->size[1] == 1) || (C3->size[0] == 1)) {
      i1 = y->size[0] * y->size[1];
      y->size[0] = rho->size[0];
      y->size[1] = C3->size[1];
      emxEnsureCapacity_creal_T(y, i1);
      loop_ub = rho->size[0];
      for (i1 = 0; i1 < loop_ub; i1++) {
        inner = C3->size[1];
        for (m = 0; m < inner; m++) {
          y->data[i1 + y->size[0] * m].re = 0.0;
          y->data[i1 + y->size[0] * m].im = 0.0;
          n = rho->size[1];
          for (j = 0; j < n; j++) {
            d_points = rho->data[i1 + rho->size[0] * j].re * C3->data[j +
              C3->size[0] * m].re - rho->data[i1 + rho->size[0] * j].im *
              C3->data[j + C3->size[0] * m].im;
            C1_im = rho->data[i1 + rho->size[0] * j].re * C3->data[j + C3->size
              [0] * m].im + rho->data[i1 + rho->size[0] * j].im * C3->data[j +
              C3->size[0] * m].re;
            y->data[i1 + y->size[0] * m].re += d_points;
            y->data[i1 + y->size[0] * m].im += C1_im;
          }
        }
      }
    } else {
      m = rho->size[0];
      inner = rho->size[1];
      n = C3->size[1];
      i1 = y->size[0] * y->size[1];
      y->size[0] = rho->size[0];
      y->size[1] = C3->size[1];
      emxEnsureCapacity_creal_T(y, i1);
      for (j = 0; j < n; j++) {
        coffset = j * m;
        boffset = j * inner;
        for (i = 0; i < m; i++) {
          i1 = coffset + i;
          y->data[i1].re = 0.0;
          y->data[i1].im = 0.0;
        }

        for (b_k = 0; b_k < inner; b_k++) {
          aoffset = b_k * m;
          loop_ub = boffset + b_k;
          temp_re = C3->data[loop_ub].re;
          temp_im = C3->data[loop_ub].im;
          for (i = 0; i < m; i++) {
            loop_ub = aoffset + i;
            d_points = temp_re * rho->data[loop_ub].re - temp_im * rho->
              data[loop_ub].im;
            C1_im = temp_re * rho->data[loop_ub].im + temp_im * rho->
              data[loop_ub].re;
            i1 = coffset + i;
            y->data[i1].re += d_points;
            y->data[i1].im += C1_im;
          }
        }
      }
    }

    i1 = c_points->size[0] * c_points->size[1];
    c_points->size[0] = E_mm->size[0];
    c_points->size[1] = E_mm->size[1];
    emxEnsureCapacity_creal_T(c_points, i1);
    d_points = points->data[k];
    loop_ub = E_mm->size[0] * E_mm->size[1];
    for (i1 = 0; i1 < loop_ub; i1++) {
      c_points->data[i1].re = (d_points - (E_mm->data[i1] - E_qq->data[i1])) +
        y_re;
      c_points->data[i1].im = i0p;
    }

    rdivide_helper(c_points, b);
    i1 = b->size[0] * b->size[1];
    m = b->size[0] * b->size[1];
    emxEnsureCapacity_creal_T(b, m);
    loop_ub = i1 - 1;
    for (i1 = 0; i1 <= loop_ub; i1++) {
      temp_im = b->data[i1].re;
      temp_re = b->data[i1].im;
      d_points = C4->data[i1].re;
      C1_im = C4->data[i1].im;
      b->data[i1].re = temp_im * d_points - temp_re * C1_im;
      b->data[i1].im = temp_im * C1_im + temp_re * d_points;
    }

    if ((y->size[1] == 1) || (b->size[0] == 1)) {
      i1 = r0->size[0] * r0->size[1];
      r0->size[0] = y->size[0];
      r0->size[1] = b->size[1];
      emxEnsureCapacity_creal_T(r0, i1);
      loop_ub = y->size[0];
      for (i1 = 0; i1 < loop_ub; i1++) {
        inner = b->size[1];
        for (m = 0; m < inner; m++) {
          r0->data[i1 + r0->size[0] * m].re = 0.0;
          r0->data[i1 + r0->size[0] * m].im = 0.0;
          n = y->size[1];
          for (j = 0; j < n; j++) {
            y_re = y->data[i1 + y->size[0] * j].re * b->data[j + b->size[0] * m]
              .re - y->data[i1 + y->size[0] * j].im * b->data[j + b->size[0] * m]
              .im;
            C1_im = y->data[i1 + y->size[0] * j].re * b->data[j + b->size[0] * m]
              .im + y->data[i1 + y->size[0] * j].im * b->data[j + b->size[0] * m]
              .re;
            r0->data[i1 + r0->size[0] * m].re += y_re;
            r0->data[i1 + r0->size[0] * m].im += C1_im;
          }
        }
      }
    } else {
      m = y->size[0];
      inner = y->size[1];
      n = b->size[1];
      i1 = r0->size[0] * r0->size[1];
      r0->size[0] = y->size[0];
      r0->size[1] = b->size[1];
      emxEnsureCapacity_creal_T(r0, i1);
      for (j = 0; j < n; j++) {
        coffset = j * m;
        boffset = j * inner;
        for (i = 0; i < m; i++) {
          i1 = coffset + i;
          r0->data[i1].re = 0.0;
          r0->data[i1].im = 0.0;
        }

        for (b_k = 0; b_k < inner; b_k++) {
          aoffset = b_k * m;
          loop_ub = boffset + b_k;
          temp_re = b->data[loop_ub].re;
          temp_im = b->data[loop_ub].im;
          for (i = 0; i < m; i++) {
            loop_ub = aoffset + i;
            d_points = temp_re * y->data[loop_ub].re - temp_im * y->data[loop_ub]
              .im;
            C1_im = temp_re * y->data[loop_ub].im + temp_im * y->data[loop_ub].
              re;
            i1 = coffset + i;
            r0->data[i1].re += d_points;
            r0->data[i1].im += C1_im;
          }
        }
      }
    }

    i1 = a->size[0] * a->size[1];
    m = a->size[0] * a->size[1];
    emxEnsureCapacity_creal_T(a, m);
    loop_ub = i1 - 1;
    for (i1 = 0; i1 <= loop_ub; i1++) {
      a->data[i1].re += r0->data[i1].re;
      a->data[i1].im += r0->data[i1].im;
    }

    temp_re = 0.0;
    temp_im = 0.0;
    i1 = a->size[0];
    for (b_k = 0; b_k < i1; b_k++) {
      temp_re += a->data[b_k + a->size[0] * b_k].re;
      temp_im += a->data[b_k + a->size[0] * b_k].im;
    }

    GR->data[k].re = temp_re;
    GR->data[k].im = temp_im;
    c_tmp = i0p * i0p;
    temp_re = i0p * 0.0;
    temp_im = i0p * -2.0;
    i1 = y->size[0] * y->size[1];
    y->size[0] = rho->size[0];
    y->size[1] = rho->size[1];
    emxEnsureCapacity_creal_T(y, i1);
    loop_ub = rho->size[0] * rho->size[1];
    for (i1 = 0; i1 < loop_ub; i1++) {
      d_points = rho->data[i1].re;
      C1_im = rho->data[i1].im;
      y->data[i1].re = temp_re * d_points - temp_im * C1_im;
      y->data[i1].im = temp_re * C1_im + temp_im * d_points;
    }

    i1 = b_points->size[0] * b_points->size[1];
    b_points->size[0] = E_pp->size[0];
    b_points->size[1] = E_pp->size[1];
    emxEnsureCapacity_real_T(b_points, i1);
    d_points = points->data[k];
    loop_ub = E_pp->size[0] * E_pp->size[1];
    for (i1 = 0; i1 < loop_ub; i1++) {
      b_points->data[i1] = d_points - (E_pp->data[i1] - E_nn->data[i1]);
    }

    power(b_points, r1);
    i1 = b_points->size[0] * b_points->size[1];
    b_points->size[0] = r1->size[0];
    b_points->size[1] = r1->size[1];
    emxEnsureCapacity_real_T(b_points, i1);
    loop_ub = r1->size[0] * r1->size[1];
    for (i1 = 0; i1 < loop_ub; i1++) {
      b_points->data[i1] = r1->data[i1] + c_tmp;
    }

    b_rdivide_helper(b_points, r1);
    i1 = b->size[0] * b->size[1];
    b->size[0] = r1->size[0];
    b->size[1] = r1->size[1];
    emxEnsureCapacity_creal_T(b, i1);
    loop_ub = r1->size[0] * r1->size[1];
    for (i1 = 0; i1 < loop_ub; i1++) {
      b->data[i1].re = r1->data[i1] * C1->data[i1].re;
      b->data[i1].im = r1->data[i1] * C1->data[i1].im;
    }

    if ((y->size[1] == 1) || (b->size[0] == 1)) {
      i1 = c_points->size[0] * c_points->size[1];
      c_points->size[0] = y->size[0];
      c_points->size[1] = b->size[1];
      emxEnsureCapacity_creal_T(c_points, i1);
      loop_ub = y->size[0];
      for (i1 = 0; i1 < loop_ub; i1++) {
        inner = b->size[1];
        for (m = 0; m < inner; m++) {
          c_points->data[i1 + c_points->size[0] * m].re = 0.0;
          c_points->data[i1 + c_points->size[0] * m].im = 0.0;
          n = y->size[1];
          for (j = 0; j < n; j++) {
            y_re = y->data[i1 + y->size[0] * j].re * b->data[j + b->size[0] * m]
              .re - y->data[i1 + y->size[0] * j].im * b->data[j + b->size[0] * m]
              .im;
            C1_im = y->data[i1 + y->size[0] * j].re * b->data[j + b->size[0] * m]
              .im + y->data[i1 + y->size[0] * j].im * b->data[j + b->size[0] * m]
              .re;
            c_points->data[i1 + c_points->size[0] * m].re += y_re;
            c_points->data[i1 + c_points->size[0] * m].im += C1_im;
          }
        }
      }
    } else {
      m = y->size[0];
      inner = y->size[1];
      n = b->size[1];
      i1 = c_points->size[0] * c_points->size[1];
      c_points->size[0] = y->size[0];
      c_points->size[1] = b->size[1];
      emxEnsureCapacity_creal_T(c_points, i1);
      for (j = 0; j < n; j++) {
        coffset = j * m;
        boffset = j * inner;
        for (i = 0; i < m; i++) {
          i1 = coffset + i;
          c_points->data[i1].re = 0.0;
          c_points->data[i1].im = 0.0;
        }

        for (b_k = 0; b_k < inner; b_k++) {
          aoffset = b_k * m;
          loop_ub = boffset + b_k;
          temp_re = b->data[loop_ub].re;
          temp_im = b->data[loop_ub].im;
          for (i = 0; i < m; i++) {
            loop_ub = aoffset + i;
            d_points = temp_re * y->data[loop_ub].re - temp_im * y->data[loop_ub]
              .im;
            C1_im = temp_re * y->data[loop_ub].im + temp_im * y->data[loop_ub].
              re;
            i1 = coffset + i;
            c_points->data[i1].re += d_points;
            c_points->data[i1].im += C1_im;
          }
        }
      }
    }

    if ((c_points->size[1] == 1) || (C2->size[0] == 1)) {
      i1 = a->size[0] * a->size[1];
      a->size[0] = c_points->size[0];
      a->size[1] = C2->size[1];
      emxEnsureCapacity_creal_T(a, i1);
      loop_ub = c_points->size[0];
      for (i1 = 0; i1 < loop_ub; i1++) {
        inner = C2->size[1];
        for (m = 0; m < inner; m++) {
          a->data[i1 + a->size[0] * m].re = 0.0;
          a->data[i1 + a->size[0] * m].im = 0.0;
          n = c_points->size[1];
          for (j = 0; j < n; j++) {
            d_points = c_points->data[i1 + c_points->size[0] * j].re * C2->
              data[j + C2->size[0] * m].re - c_points->data[i1 + c_points->size
              [0] * j].im * C2->data[j + C2->size[0] * m].im;
            C1_im = c_points->data[i1 + c_points->size[0] * j].re * C2->data[j +
              C2->size[0] * m].im + c_points->data[i1 + c_points->size[0] * j].
              im * C2->data[j + C2->size[0] * m].re;
            a->data[i1 + a->size[0] * m].re += d_points;
            a->data[i1 + a->size[0] * m].im += C1_im;
          }
        }
      }
    } else {
      m = c_points->size[0];
      inner = c_points->size[1];
      n = C2->size[1];
      i1 = a->size[0] * a->size[1];
      a->size[0] = c_points->size[0];
      a->size[1] = C2->size[1];
      emxEnsureCapacity_creal_T(a, i1);
      for (j = 0; j < n; j++) {
        coffset = j * m;
        boffset = j * inner;
        for (i = 0; i < m; i++) {
          i1 = coffset + i;
          a->data[i1].re = 0.0;
          a->data[i1].im = 0.0;
        }

        for (b_k = 0; b_k < inner; b_k++) {
          aoffset = b_k * m;
          loop_ub = boffset + b_k;
          temp_re = C2->data[loop_ub].re;
          temp_im = C2->data[loop_ub].im;
          for (i = 0; i < m; i++) {
            loop_ub = aoffset + i;
            d_points = temp_re * c_points->data[loop_ub].re - temp_im *
              c_points->data[loop_ub].im;
            C1_im = temp_re * c_points->data[loop_ub].im + temp_im *
              c_points->data[loop_ub].re;
            i1 = coffset + i;
            a->data[i1].re += d_points;
            a->data[i1].im += C1_im;
          }
        }
      }
    }

    temp_re = i0p * 0.0;
    temp_im = i0p * 2.0;
    i1 = y->size[0] * y->size[1];
    y->size[0] = rho->size[0];
    y->size[1] = rho->size[1];
    emxEnsureCapacity_creal_T(y, i1);
    loop_ub = rho->size[0] * rho->size[1];
    for (i1 = 0; i1 < loop_ub; i1++) {
      d_points = rho->data[i1].re;
      C1_im = rho->data[i1].im;
      y->data[i1].re = temp_re * d_points - temp_im * C1_im;
      y->data[i1].im = temp_re * C1_im + temp_im * d_points;
    }

    if ((y->size[1] == 1) || (C3->size[0] == 1)) {
      i1 = c_points->size[0] * c_points->size[1];
      c_points->size[0] = y->size[0];
      c_points->size[1] = C3->size[1];
      emxEnsureCapacity_creal_T(c_points, i1);
      loop_ub = y->size[0];
      for (i1 = 0; i1 < loop_ub; i1++) {
        inner = C3->size[1];
        for (m = 0; m < inner; m++) {
          c_points->data[i1 + c_points->size[0] * m].re = 0.0;
          c_points->data[i1 + c_points->size[0] * m].im = 0.0;
          n = y->size[1];
          for (j = 0; j < n; j++) {
            y_re = y->data[i1 + y->size[0] * j].re * C3->data[j + C3->size[0] *
              m].re - y->data[i1 + y->size[0] * j].im * C3->data[j + C3->size[0]
              * m].im;
            C1_im = y->data[i1 + y->size[0] * j].re * C3->data[j + C3->size[0] *
              m].im + y->data[i1 + y->size[0] * j].im * C3->data[j + C3->size[0]
              * m].re;
            c_points->data[i1 + c_points->size[0] * m].re += y_re;
            c_points->data[i1 + c_points->size[0] * m].im += C1_im;
          }
        }
      }
    } else {
      m = y->size[0];
      inner = y->size[1];
      n = C3->size[1];
      i1 = c_points->size[0] * c_points->size[1];
      c_points->size[0] = y->size[0];
      c_points->size[1] = C3->size[1];
      emxEnsureCapacity_creal_T(c_points, i1);
      for (j = 0; j < n; j++) {
        coffset = j * m;
        boffset = j * inner;
        for (i = 0; i < m; i++) {
          i1 = coffset + i;
          c_points->data[i1].re = 0.0;
          c_points->data[i1].im = 0.0;
        }

        for (b_k = 0; b_k < inner; b_k++) {
          aoffset = b_k * m;
          loop_ub = boffset + b_k;
          temp_re = C3->data[loop_ub].re;
          temp_im = C3->data[loop_ub].im;
          for (i = 0; i < m; i++) {
            loop_ub = aoffset + i;
            d_points = temp_re * y->data[loop_ub].re - temp_im * y->data[loop_ub]
              .im;
            C1_im = temp_re * y->data[loop_ub].im + temp_im * y->data[loop_ub].
              re;
            i1 = coffset + i;
            c_points->data[i1].re += d_points;
            c_points->data[i1].im += C1_im;
          }
        }
      }
    }

    i1 = b_points->size[0] * b_points->size[1];
    b_points->size[0] = E_mm->size[0];
    b_points->size[1] = E_mm->size[1];
    emxEnsureCapacity_real_T(b_points, i1);
    d_points = points->data[k];
    loop_ub = E_mm->size[0] * E_mm->size[1];
    for (i1 = 0; i1 < loop_ub; i1++) {
      b_points->data[i1] = d_points - (E_mm->data[i1] - E_qq->data[i1]);
    }

    power(b_points, r1);
    i1 = b_points->size[0] * b_points->size[1];
    b_points->size[0] = r1->size[0];
    b_points->size[1] = r1->size[1];
    emxEnsureCapacity_real_T(b_points, i1);
    loop_ub = r1->size[0] * r1->size[1];
    for (i1 = 0; i1 < loop_ub; i1++) {
      b_points->data[i1] = r1->data[i1] + c_tmp;
    }

    b_rdivide_helper(b_points, r1);
    i1 = b->size[0] * b->size[1];
    b->size[0] = r1->size[0];
    b->size[1] = r1->size[1];
    emxEnsureCapacity_creal_T(b, i1);
    loop_ub = r1->size[0] * r1->size[1];
    for (i1 = 0; i1 < loop_ub; i1++) {
      b->data[i1].re = r1->data[i1] * C4->data[i1].re;
      b->data[i1].im = r1->data[i1] * C4->data[i1].im;
    }

    if ((c_points->size[1] == 1) || (b->size[0] == 1)) {
      i1 = r0->size[0] * r0->size[1];
      r0->size[0] = c_points->size[0];
      r0->size[1] = b->size[1];
      emxEnsureCapacity_creal_T(r0, i1);
      loop_ub = c_points->size[0];
      for (i1 = 0; i1 < loop_ub; i1++) {
        inner = b->size[1];
        for (m = 0; m < inner; m++) {
          r0->data[i1 + r0->size[0] * m].re = 0.0;
          r0->data[i1 + r0->size[0] * m].im = 0.0;
          n = c_points->size[1];
          for (j = 0; j < n; j++) {
            d_points = c_points->data[i1 + c_points->size[0] * j].re * b->data[j
              + b->size[0] * m].re - c_points->data[i1 + c_points->size[0] * j].
              im * b->data[j + b->size[0] * m].im;
            C1_im = c_points->data[i1 + c_points->size[0] * j].re * b->data[j +
              b->size[0] * m].im + c_points->data[i1 + c_points->size[0] * j].im
              * b->data[j + b->size[0] * m].re;
            r0->data[i1 + r0->size[0] * m].re += d_points;
            r0->data[i1 + r0->size[0] * m].im += C1_im;
          }
        }
      }
    } else {
      m = c_points->size[0];
      inner = c_points->size[1];
      n = b->size[1];
      i1 = r0->size[0] * r0->size[1];
      r0->size[0] = c_points->size[0];
      r0->size[1] = b->size[1];
      emxEnsureCapacity_creal_T(r0, i1);
      for (j = 0; j < n; j++) {
        coffset = j * m;
        boffset = j * inner;
        for (i = 0; i < m; i++) {
          i1 = coffset + i;
          r0->data[i1].re = 0.0;
          r0->data[i1].im = 0.0;
        }

        for (b_k = 0; b_k < inner; b_k++) {
          aoffset = b_k * m;
          loop_ub = boffset + b_k;
          temp_re = b->data[loop_ub].re;
          temp_im = b->data[loop_ub].im;
          for (i = 0; i < m; i++) {
            loop_ub = aoffset + i;
            d_points = temp_re * c_points->data[loop_ub].re - temp_im *
              c_points->data[loop_ub].im;
            C1_im = temp_re * c_points->data[loop_ub].im + temp_im *
              c_points->data[loop_ub].re;
            i1 = coffset + i;
            r0->data[i1].re += d_points;
            r0->data[i1].im += C1_im;
          }
        }
      }
    }

    i1 = a->size[0] * a->size[1];
    m = a->size[0] * a->size[1];
    emxEnsureCapacity_creal_T(a, m);
    loop_ub = i1 - 1;
    for (i1 = 0; i1 <= loop_ub; i1++) {
      a->data[i1].re += r0->data[i1].re;
      a->data[i1].im += r0->data[i1].im;
    }

    temp_re = 0.0;
    temp_im = 0.0;
    i1 = a->size[0];
    for (b_k = 0; b_k < i1; b_k++) {
      temp_re += a->data[b_k + a->size[0] * b_k].re;
      temp_im += a->data[b_k + a->size[0] * b_k].im;
    }

    GK->data[k].re = temp_re;
    GK->data[k].im = temp_im;
  }

  emxFree_creal_T(&c_points);
  emxFree_real_T(&b_points);
  emxFree_real_T(&r1);
  emxFree_creal_T(&a);
  emxFree_creal_T(&b);
  emxFree_creal_T(&y);
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
}

/* End of code generation (m_GF.c) */
