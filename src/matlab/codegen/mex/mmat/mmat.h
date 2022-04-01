/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * mmat.h
 *
 * Code generation for function 'mmat'
 *
 */

#pragma once

/* Include files */
#include "mmat_types.h"
#include "rtwtypes.h"
#include "emlrt.h"
#include "mex.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* Function Declarations */
void mmat(const emlrtStack *sp, real_T E, const emxArray_real_T *xz,
          const emxArray_real_T *u, real_T newM[4], creal_T *k, real_T *delx);

/* End of code generation (mmat.h) */
