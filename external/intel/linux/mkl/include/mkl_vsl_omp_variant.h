/*******************************************************************************
* Copyright 2020-2021 Intel Corporation.
*
* This software and the related documents are Intel copyrighted  materials,  and
* your use of  them is  governed by the  express license  under which  they were
* provided to you (License).  Unless the License provides otherwise, you may not
* use, modify, copy, publish, distribute,  disclose or transmit this software or
* the related documents without Intel's prior written permission.
*
* This software and the related documents  are provided as  is,  with no express
* or implied  warranties,  other  than those  that are  expressly stated  in the
* License.
*******************************************************************************/

/*
!  Content:
!      Intel(R) Math Kernel Library (Intel(R) MKL) for OpenMP compiler offload
!      interface
!******************************************************************************/

#ifndef _MKL_VSL_OMP_VARIANT_H_
#define _MKL_VSL_OMP_VARIANT_H_

#include "mkl_vsl_types.h"
#include "mkl_omp_variant.h"

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

int MKL_VARIANT_NAME(vsl, vsRngUniform)(const MKL_INT, VSLStreamStatePtr, const MKL_INT, float [], const float, const float) NOTHROW;

int MKL_VARIANT_NAME(vsl, vdRngUniform)(const MKL_INT, VSLStreamStatePtr, const MKL_INT, double [], const double, const double) NOTHROW;

int MKL_VARIANT_NAME(vsl, viRngUniform)(const MKL_INT, VSLStreamStatePtr, const MKL_INT, int [], const int, const int) NOTHROW;

int MKL_VARIANT_NAME(vsl, vsRngGaussian)(const MKL_INT, VSLStreamStatePtr, const MKL_INT, float [], const float, const float) NOTHROW;

int MKL_VARIANT_NAME(vsl, vdRngGaussian)(const MKL_INT, VSLStreamStatePtr, const MKL_INT, double [], const double, const double) NOTHROW;

int MKL_VARIANT_NAME(vsl, vsRngGaussianMV)(const MKL_INT, VSLStreamStatePtr, const MKL_INT, float [], const MKL_INT, const MKL_INT, const float*, const float*) NOTHROW;

int MKL_VARIANT_NAME(vsl, vdRngGaussianMV)(const MKL_INT, VSLStreamStatePtr, const MKL_INT, double [], const MKL_INT, const MKL_INT, const double*, const double*) NOTHROW;

int MKL_VARIANT_NAME(vsl, vsRngLognormal)(const MKL_INT, VSLStreamStatePtr, const MKL_INT, float [], const float, const float, const float, const float) NOTHROW;

int MKL_VARIANT_NAME(vsl, vdRngLognormal)(const MKL_INT, VSLStreamStatePtr, const MKL_INT, double [], const double, const double, const double, const double) NOTHROW;

int MKL_VARIANT_NAME(vsl, vsRngCauchy)(const MKL_INT, VSLStreamStatePtr, const MKL_INT, float [], const float, const float) NOTHROW;

int MKL_VARIANT_NAME(vsl, vdRngCauchy)(const MKL_INT, VSLStreamStatePtr, const MKL_INT, double [], const double, const double) NOTHROW;

int MKL_VARIANT_NAME(vsl, vsRngExponential)(const MKL_INT, VSLStreamStatePtr, const MKL_INT, float [], const float, const float) NOTHROW;

int MKL_VARIANT_NAME(vsl, vdRngExponential)(const MKL_INT, VSLStreamStatePtr, const MKL_INT, double [], const double, const double) NOTHROW;

int MKL_VARIANT_NAME(vsl, vsRngGumbel)(const MKL_INT, VSLStreamStatePtr, const MKL_INT, float [], const float, const float) NOTHROW;

int MKL_VARIANT_NAME(vsl, vdRngGumbel)(const MKL_INT, VSLStreamStatePtr, const MKL_INT, double [], const double, const double) NOTHROW;

int MKL_VARIANT_NAME(vsl, vsRngLaplace)(const MKL_INT, VSLStreamStatePtr, const MKL_INT, float [], const float, const float) NOTHROW;

int MKL_VARIANT_NAME(vsl, vdRngLaplace)(const MKL_INT, VSLStreamStatePtr, const MKL_INT, double [], const double, const double) NOTHROW;

int MKL_VARIANT_NAME(vsl, vsRngRayleigh)(const MKL_INT, VSLStreamStatePtr, const MKL_INT, float [], const float, const float) NOTHROW;

int MKL_VARIANT_NAME(vsl, vdRngRayleigh)(const MKL_INT, VSLStreamStatePtr, const MKL_INT, double [], const double, const double) NOTHROW;

int MKL_VARIANT_NAME(vsl, vsRngWeibull)(const MKL_INT, VSLStreamStatePtr, const MKL_INT, float [], const float, const float, const float) NOTHROW;

int MKL_VARIANT_NAME(vsl, vdRngWeibull)(const MKL_INT, VSLStreamStatePtr, const MKL_INT, double [], const double, const double, const double) NOTHROW;

int MKL_VARIANT_NAME(vsl, viRngBernoulli)(const MKL_INT, VSLStreamStatePtr, const MKL_INT, int [], const double) NOTHROW;

int MKL_VARIANT_NAME(vsl, viRngGeometric)(const MKL_INT, VSLStreamStatePtr, const MKL_INT, int [], const double) NOTHROW;

int MKL_VARIANT_NAME(vsl, viRngPoisson)(const MKL_INT, VSLStreamStatePtr, const MKL_INT, int [], const double) NOTHROW;

int MKL_VARIANT_NAME(vsl, viRngUniformBits)(const MKL_INT, VSLStreamStatePtr, const MKL_INT, unsigned int []) NOTHROW;

int MKL_VARIANT_NAME(vsl, viRngUniformBits32)(const MKL_INT, VSLStreamStatePtr, const MKL_INT, unsigned int []) NOTHROW;

int MKL_VARIANT_NAME(vsl, viRngUniformBits64)(const MKL_INT, VSLStreamStatePtr, const MKL_INT, unsigned MKL_INT64 []) NOTHROW;

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif // _MKL_VSL_OMP_VARIANT_H_
