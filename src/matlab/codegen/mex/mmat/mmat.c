/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * mmat.c
 *
 * Code generation for function 'mmat'
 *
 */

/* Include files */
#include "mmat.h"
#include "mmat_types.h"
#include "rt_nonfinite.h"
#include "mmat_automatic.h"
#include "mwmathutil.h"

/* Variable Definitions */
static emlrtRSInfo emlrtRSI = {
    20,     /* lineNo */
    "mmat", /* fcnName */
    "/Users/idesjard/Documents/Planetary Surface and Spacecraft Lab/Plasma "
    "Debris Wake Research/KVIST/src/matlab/mmat.m" /* pathName */
};

static emlrtMCInfo emlrtMCI = {
    27,      /* lineNo */
    5,       /* colNo */
    "error", /* fName */
    "/Applications/MATLAB_R2021b.app/toolbox/eml/lib/matlab/lang/error.m" /* pName
                                                                           */
};

static emlrtRSInfo b_emlrtRSI = {
    27,      /* lineNo */
    "error", /* fcnName */
    "/Applications/MATLAB_R2021b.app/toolbox/eml/lib/matlab/lang/error.m" /* pathName
                                                                           */
};

/* Function Declarations */
static void b_error(const emlrtStack *sp, const mxArray *m, const mxArray *m1,
                    emlrtMCInfo *location);

/* Function Definitions */
static void b_error(const emlrtStack *sp, const mxArray *m, const mxArray *m1,
                    emlrtMCInfo *location)
{
  const mxArray *pArrays[2];
  pArrays[0] = m;
  pArrays[1] = m1;
  emlrtCallMATLABR2012b((emlrtCTX)sp, 0, NULL, 2, &pArrays[0],
                        (const char_T *)"error", true, location);
}

void mmat(const emlrtStack *sp, real_T E, const emxArray_real_T *xz,
          const emxArray_real_T *u, real_T newM[4], creal_T *k, real_T *delx)
{
  static const int32_T iv[2] = {1, 33};
  static const char_T varargin_1[33] = {
      'S', 'o', 'm', 'e', 't', 'h', 'i', 'n', 'g', ' ',  'w',
      'e', 'n', 't', ' ', 'w', 'r', 'o', 'n', 'g', ',',  ' ',
      'f', 'l', 'a', 'g', ' ', '=', ' ', '%', 'g', '\\', 'n'};
  emlrtStack b_st;
  emlrtStack st;
  const mxArray *b_y;
  const mxArray *m;
  const mxArray *y;
  real_T M[4];
  const real_T *u_data;
  const real_T *xz_data;
  real_T yi;
  real_T yr;
  int32_T flag;
  st.prev = sp;
  st.tls = sp->tls;
  b_st.prev = &st;
  b_st.tls = st.tls;
  u_data = u->data;
  xz_data = xz->data;
  /*  for code generation, preinitialize the output variable */
  /*  data type, size, and complexity  */
  /*  generate an include in the C code */
  if (E < 0.0) {
    yr = 0.0;
    yi = muDoubleScalarSqrt(-E);
  } else {
    yr = muDoubleScalarSqrt(E);
    yi = 0.0;
  }
  k->re = yr;
  k->im = yi;
  /*  evaluate the C function */
  flag = mmat_automatic(E, &xz_data[0], &u_data[0], xz->size[1], &M[0], delx);
  if (flag < 0) {
    st.site = &emlrtRSI;
    y = NULL;
    m = emlrtCreateCharArray(2, &iv[0]);
    emlrtInitCharArrayR2013a(&st, 33, m, &varargin_1[0]);
    emlrtAssign(&y, m);
    b_y = NULL;
    m = emlrtCreateNumericMatrix(1, 1, mxINT32_CLASS, mxREAL);
    *(int32_T *)emlrtMxGetData(m) = flag;
    emlrtAssign(&b_y, m);
    b_st.site = &b_emlrtRSI;
    b_error(&b_st, y, b_y, &emlrtMCI);
  }
  newM[0] = M[0];
  newM[2] = M[1];
  newM[1] = M[2];
  newM[3] = M[3];
}

/* End of code generation (mmat.c) */
