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

#ifndef _MKL_VSL_OMP_OFFLOAD_H_
#define _MKL_VSL_OMP_OFFLOAD_H_

#include "mkl_vsl_types.h"
#include "mkl_vsl_omp_variant.h"

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

#pragma omp declare variant (MKL_VARIANT_NAME(vsl, vsRngUniform)) match(construct={target variant dispatch}, device={arch(gen)})
int vsRngUniform(const MKL_INT, VSLStreamStatePtr, const MKL_INT, float [], const float, const float) NOTHROW;

#pragma omp declare variant (MKL_VARIANT_NAME(vsl, vdRngUniform)) match(construct={target variant dispatch}, device={arch(gen)})
int vdRngUniform(const MKL_INT, VSLStreamStatePtr, const MKL_INT, double [], const double, const double) NOTHROW;

#pragma omp declare variant (MKL_VARIANT_NAME(vsl, viRngUniform)) match(construct={target variant dispatch}, device={arch(gen)})
int viRngUniform(const MKL_INT, VSLStreamStatePtr, const MKL_INT, int [], const int, const int) NOTHROW;

#pragma omp declare variant (MKL_VARIANT_NAME(vsl, vsRngGaussian)) match(construct={target variant dispatch}, device={arch(gen)})
int vsRngGaussian(const MKL_INT, VSLStreamStatePtr, const MKL_INT, float [], const float, const float) NOTHROW;

#pragma omp declare variant (MKL_VARIANT_NAME(vsl, vdRngGaussian)) match(construct={target variant dispatch}, device={arch(gen)})
int vdRngGaussian(const MKL_INT, VSLStreamStatePtr, const MKL_INT, double [], const double, const double) NOTHROW;

#pragma omp declare variant (MKL_VARIANT_NAME(vsl, vsRngGaussianMV)) match(construct={target variant dispatch}, device={arch(gen)})
int vsRngGaussianMV(const MKL_INT, VSLStreamStatePtr, const MKL_INT, float [], const MKL_INT, const MKL_INT, const float*, const float*) NOTHROW;

#pragma omp declare variant (MKL_VARIANT_NAME(vsl, vdRngGaussianMV)) match(construct={target variant dispatch}, device={arch(gen)})
int vdRngGaussianMV(const MKL_INT, VSLStreamStatePtr, const MKL_INT, double [], const MKL_INT, const MKL_INT, const double*, const double*) NOTHROW;

#pragma omp declare variant (MKL_VARIANT_NAME(vsl, vsRngLognormal)) match(construct={target variant dispatch}, device={arch(gen)})
int vsRngLognormal(const MKL_INT, VSLStreamStatePtr, const MKL_INT, float [], const float, const float, const float, const float) NOTHROW;

#pragma omp declare variant (MKL_VARIANT_NAME(vsl, vdRngLognormal)) match(construct={target variant dispatch}, device={arch(gen)})
int vdRngLognormal(const MKL_INT, VSLStreamStatePtr, const MKL_INT, double [], const double, const double, const double, const double) NOTHROW;

#pragma omp declare variant (MKL_VARIANT_NAME(vsl, vsRngCauchy)) match(construct={target variant dispatch}, device={arch(gen)})
int vsRngCauchy(const MKL_INT, VSLStreamStatePtr, const MKL_INT, float [], const float, const float) NOTHROW;

#pragma omp declare variant (MKL_VARIANT_NAME(vsl, vdRngCauchy)) match(construct={target variant dispatch}, device={arch(gen)})
int vdRngCauchy(const MKL_INT, VSLStreamStatePtr, const MKL_INT, double [], const double, const double) NOTHROW;

#pragma omp declare variant (MKL_VARIANT_NAME(vsl, vsRngExponential)) match(construct={target variant dispatch}, device={arch(gen)})
int vsRngExponential(const MKL_INT, VSLStreamStatePtr, const MKL_INT, float [], const float, const float) NOTHROW;

#pragma omp declare variant (MKL_VARIANT_NAME(vsl, vdRngExponential)) match(construct={target variant dispatch}, device={arch(gen)})
int vdRngExponential(const MKL_INT, VSLStreamStatePtr, const MKL_INT, double [], const double, const double) NOTHROW;

#pragma omp declare variant (MKL_VARIANT_NAME(vsl, vsRngGumbel)) match(construct={target variant dispatch}, device={arch(gen)})
int vsRngGumbel(const MKL_INT, VSLStreamStatePtr, const MKL_INT, float [], const float, const float) NOTHROW;

#pragma omp declare variant (MKL_VARIANT_NAME(vsl, vdRngGumbel)) match(construct={target variant dispatch}, device={arch(gen)})
int vdRngGumbel(const MKL_INT, VSLStreamStatePtr, const MKL_INT, double [], const double, const double) NOTHROW;

#pragma omp declare variant (MKL_VARIANT_NAME(vsl, vsRngLaplace)) match(construct={target variant dispatch}, device={arch(gen)})
int vsRngLaplace(const MKL_INT, VSLStreamStatePtr, const MKL_INT, float [], const float, const float) NOTHROW;

#pragma omp declare variant (MKL_VARIANT_NAME(vsl, vdRngLaplace)) match(construct={target variant dispatch}, device={arch(gen)})
int vdRngLaplace(const MKL_INT, VSLStreamStatePtr, const MKL_INT, double [], const double, const double) NOTHROW;

#pragma omp declare variant (MKL_VARIANT_NAME(vsl, vsRngRayleigh)) match(construct={target variant dispatch}, device={arch(gen)})
int vsRngRayleigh(const MKL_INT, VSLStreamStatePtr, const MKL_INT, float [], const float, const float) NOTHROW;

#pragma omp declare variant (MKL_VARIANT_NAME(vsl, vdRngRayleigh)) match(construct={target variant dispatch}, device={arch(gen)})
int vdRngRayleigh(const MKL_INT, VSLStreamStatePtr, const MKL_INT, double [], const double, const double) NOTHROW;

#pragma omp declare variant (MKL_VARIANT_NAME(vsl, vsRngWeibull)) match(construct={target variant dispatch}, device={arch(gen)})
int vsRngWeibull(const MKL_INT, VSLStreamStatePtr, const MKL_INT, float [], const float, const float, const float) NOTHROW;

#pragma omp declare variant (MKL_VARIANT_NAME(vsl, vdRngWeibull)) match(construct={target variant dispatch}, device={arch(gen)})
int vdRngWeibull(const MKL_INT, VSLStreamStatePtr, const MKL_INT, double [], const double, const double, const double) NOTHROW;

#pragma omp declare variant (MKL_VARIANT_NAME(vsl, viRngBernoulli)) match(construct={target variant dispatch}, device={arch(gen)})
int viRngBernoulli(const MKL_INT, VSLStreamStatePtr, const MKL_INT, int [], const double) NOTHROW;

#pragma omp declare variant (MKL_VARIANT_NAME(vsl, viRngGeometric)) match(construct={target variant dispatch}, device={arch(gen)})
int viRngGeometric(const MKL_INT, VSLStreamStatePtr, const MKL_INT, int [], const double) NOTHROW;

#pragma omp declare variant (MKL_VARIANT_NAME(vsl, viRngPoisson)) match(construct={target variant dispatch}, device={arch(gen)})
int viRngPoisson(const MKL_INT, VSLStreamStatePtr, const MKL_INT, int [], const double) NOTHROW;

#pragma omp declare variant (MKL_VARIANT_NAME(vsl, viRngUniformBits)) match(construct={target variant dispatch}, device={arch(gen)})
int viRngUniformBits(const MKL_INT, VSLStreamStatePtr, const MKL_INT, unsigned int []) NOTHROW;

#pragma omp declare variant (MKL_VARIANT_NAME(vsl, viRngUniformBits32)) match(construct={target variant dispatch}, device={arch(gen)})
int viRngUniformBits32(const MKL_INT, VSLStreamStatePtr, const MKL_INT, unsigned int []) NOTHROW;

#pragma omp declare variant (MKL_VARIANT_NAME(vsl, viRngUniformBits64)) match(construct={target variant dispatch}, device={arch(gen)})
int viRngUniformBits64(const MKL_INT, VSLStreamStatePtr, const MKL_INT, unsigned MKL_INT64 []) NOTHROW;

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif // _MKL_VSL_OMP_OFFLOAD_H_
