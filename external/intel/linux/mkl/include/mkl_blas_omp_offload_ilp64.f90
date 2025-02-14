!===============================================================================
! Copyright 2020-2021 Intel Corporation.
!
! This software and the related documents are Intel copyrighted  materials,  and
! your use of  them is  governed by the  express license  under which  they were
! provided to you (License).  Unless the License provides otherwise, you may not
! use, modify, copy, publish, distribute,  disclose or transmit this software or
! the related documents without Intel's prior written permission.
!
! This software and the related documents  are provided as  is,  with no express
! or implied  warranties,  other  than those  that are  expressly stated  in the
! License.
!===============================================================================

!  Content:
!      Intel(R) oneAPI Math Kernel Library (oneMKL) FORTRAN interface for BLAS
!      OpenMP offload
!*******************************************************************************

module onemkl_blas_omp_offload_ilp64

  include "mkl_blas_omp_variant_ilp64.f90"

  interface

     ! BLAS Level3

     subroutine dgemm ( transa, transb, m, n, k, alpha, a, lda,        &
          &b, ldb, beta, c, ldc ) BIND(C)
       character*1,intent(in)             :: transa, transb
       integer,intent(in)                 :: m, n, k, lda, ldb, ldc
       double precision,intent(in)        :: alpha, beta
       double precision,intent(in)        :: a( lda, * ), b( ldb, * )
       double precision,intent(inout)     :: c( ldc, * )
       !$omp  declare variant( dgemm:mkl_blas_dgemm_omp_offload_ilp64 ) match( construct={target variant dispatch}, device={arch(gen)} )
     end subroutine dgemm

     subroutine sgemm ( transa, transb, m, n, k, alpha, a, lda,        &
          &b, ldb, beta, c, ldc ) BIND(C)
       character*1,intent(in)              :: transa, transb
       integer,intent(in)                  :: m, n, k, lda, ldb, ldc
       real,intent(in)                     :: alpha, beta
       real,intent(in)                     :: a( lda, * ), b( ldb, * )
       real,intent(in)                     :: c( ldc, * )
       !$omp  declare variant( sgemm:mkl_blas_sgemm_omp_offload_ilp64 ) match( construct={target variant dispatch}, device={arch(gen)} )
     end subroutine sgemm

     subroutine zgemm ( transa, transb, m, n, k, alpha, a, lda,        &
          &b, ldb, beta, c, ldc ) BIND(C)
       character*1,intent(in)             :: transa, transb
       integer,intent(in)                 :: m, n, k, lda, ldb, ldc
       complex*16,intent(in)              :: alpha, beta
       complex*16,intent(in)              :: a( lda, * ), b( ldb, * )
       complex*16,intent(inout)           :: c( ldc, * )
       !$omp  declare variant( zgemm:mkl_blas_zgemm_omp_offload_ilp64 ) match( construct={target variant dispatch}, device={arch(gen)} )
     end subroutine zgemm

     subroutine cgemm ( transa, transb, m, n, k, alpha, a, lda,        &
          &b, ldb, beta, c, ldc ) BIND(C)
       character*1,intent(in)              :: transa, transb
       integer,intent(in)                  :: m, n, k, lda, ldb, ldc
       complex,intent(in)                  :: alpha, beta
       complex,intent(in)                  :: a( lda, * ), b( ldb, * )
       complex,intent(in)                  :: c( ldc, * )
       !$omp  declare variant( cgemm:mkl_blas_cgemm_omp_offload_ilp64 ) match( construct={target variant dispatch}, device={arch(gen)} )
     end subroutine cgemm

     subroutine dsymm ( side, uplo, m, n, alpha, a, lda,        &
          &b, ldb, beta, c, ldc ) BIND(C)
       character*1,intent(in)             :: side, uplo
       integer,intent(in)                 :: m, n, lda, ldb, ldc
       double precision,intent(in)        :: alpha, beta
       double precision,intent(in)        :: a( lda, * ), b( ldb, * )
       double precision,intent(inout)     :: c( ldc, * )
       !$omp  declare variant( dsymm:mkl_blas_dsymm_omp_offload_ilp64 ) match( construct={target variant dispatch}, device={arch(gen)} )
     end subroutine dsymm

     subroutine ssymm ( side, uplo, m, n, alpha, a, lda,        &
          &b, ldb, beta, c, ldc ) BIND(C)
       character*1,intent(in)              :: side, uplo
       integer,intent(in)                  :: m, n, lda, ldb, ldc
       real,intent(in)                     :: alpha, beta
       real,intent(in)                     :: a( lda, * ), b( ldb, * )
       real,intent(in)                     :: c( ldc, * )
       !$omp  declare variant( ssymm:mkl_blas_ssymm_omp_offload_ilp64 ) match( construct={target variant dispatch}, device={arch(gen)} )
     end subroutine ssymm

     subroutine zsymm ( side, uplo, m, n, alpha, a, lda,        &
          &b, ldb, beta, c, ldc ) BIND(C)
       character*1,intent(in)             :: side, uplo
       integer,intent(in)                 :: m, n, lda, ldb, ldc
       complex*16,intent(in)              :: alpha, beta
       complex*16,intent(in)              :: a( lda, * ), b( ldb, * )
       complex*16,intent(inout)           :: c( ldc, * )
       !$omp  declare variant( zsymm:mkl_blas_zsymm_omp_offload_ilp64 ) match( construct={target variant dispatch}, device={arch(gen)} )
     end subroutine zsymm

     subroutine csymm ( side, uplo, m, n, alpha, a, lda,        &
          &b, ldb, beta, c, ldc ) BIND(C)
       character*1,intent(in)              :: side, uplo
       integer,intent(in)                  :: m, n, lda, ldb, ldc
       complex,intent(in)                  :: alpha, beta
       complex,intent(in)                  :: a( lda, * ), b( ldb, * )
       complex,intent(in)                  :: c( ldc, * )
       !$omp  declare variant( csymm:mkl_blas_csymm_omp_offload_ilp64 ) match( construct={target variant dispatch}, device={arch(gen)} )
     end subroutine csymm

     subroutine zhemm ( side, uplo, m, n, alpha, a, lda,        &
          &b, ldb, beta, c, ldc ) BIND(C)
       character*1,intent(in)             :: side, uplo
       integer,intent(in)                 :: m, n, lda, ldb, ldc
       complex*16,intent(in)              :: alpha, beta
       complex*16,intent(in)              :: a( lda, * ), b( ldb, * )
       complex*16,intent(inout)           :: c( ldc, * )
       !$omp  declare variant( zhemm:mkl_blas_zhemm_omp_offload_ilp64 ) match( construct={target variant dispatch}, device={arch(gen)} )
     end subroutine zhemm

     subroutine chemm ( side, uplo, m, n, alpha, a, lda,        &
          &b, ldb, beta, c, ldc ) BIND(C)
       character*1,intent(in)              :: side, uplo
       integer,intent(in)                  :: m, n, lda, ldb, ldc
       complex,intent(in)                  :: alpha, beta
       complex,intent(in)                  :: a( lda, * ), b( ldb, * )
       complex,intent(in)                  :: c( ldc, * )
       !$omp  declare variant( chemm:mkl_blas_chemm_omp_offload_ilp64 ) match( construct={target variant dispatch}, device={arch(gen)} )
     end subroutine chemm

     subroutine dsyrk ( uplo, trans, n, k, alpha, a, lda,        &
          &beta, c, ldc ) BIND(C)
       character*1,intent(in)             :: uplo, trans
       integer,intent(in)                 :: n, k, lda, ldc
       double precision,intent(in)        :: alpha, beta
       double precision,intent(in)        :: a( lda, * )
       double precision,intent(inout)     :: c( ldc, * )
       !$omp  declare variant( dsyrk:mkl_blas_dsyrk_omp_offload_ilp64 ) match( construct={target variant dispatch}, device={arch(gen)} )
     end subroutine dsyrk

     subroutine ssyrk ( uplo, trans, n, k, alpha, a, lda,        &
          &beta, c, ldc ) BIND(C)
       character*1,intent(in)              :: uplo, trans
       integer,intent(in)                  :: n, k, lda, ldc
       real,intent(in)                     :: alpha, beta
       real,intent(in)                     :: a( lda, * )
       real,intent(in)                     :: c( ldc, * )
       !$omp  declare variant( ssyrk:mkl_blas_ssyrk_omp_offload_ilp64 ) match( construct={target variant dispatch}, device={arch(gen)} )
     end subroutine ssyrk

     subroutine zsyrk ( uplo, trans, n, k, alpha, a, lda,        &
          &beta, c, ldc ) BIND(C)
       character*1,intent(in)             :: uplo, trans
       integer,intent(in)                 :: n, k, lda, ldc
       complex*16,intent(in)              :: alpha, beta
       complex*16,intent(in)              :: a( lda, * )
       complex*16,intent(inout)           :: c( ldc, * )
       !$omp  declare variant( zsyrk:mkl_blas_zsyrk_omp_offload_ilp64 ) match( construct={target variant dispatch}, device={arch(gen)} )
     end subroutine zsyrk

     subroutine csyrk ( uplo, trans, n, k, alpha, a, lda,        &
          &beta, c, ldc ) BIND(C)
       character*1,intent(in)              :: uplo, trans
       integer,intent(in)                  :: n, k, lda, ldc
       complex,intent(in)                  :: alpha, beta
       complex,intent(in)                  :: a( lda, * )
       complex,intent(in)                  :: c( ldc, * )
       !$omp  declare variant( csyrk:mkl_blas_csyrk_omp_offload_ilp64 ) match( construct={target variant dispatch}, device={arch(gen)} )
     end subroutine csyrk

     subroutine zherk ( uplo, trans, n, k, alpha, a, lda,        &
          &beta, c, ldc ) BIND(C)
       character*1,intent(in)             :: uplo, trans
       integer,intent(in)                 :: n, k, lda, ldc
       double precision,intent(in)        :: alpha, beta
       complex*16,intent(in)              :: a( lda, * )
       complex*16,intent(inout)           :: c( ldc, * )
       !$omp  declare variant( zherk:mkl_blas_zherk_omp_offload_ilp64 ) match( construct={target variant dispatch}, device={arch(gen)} )
     end subroutine zherk

     subroutine cherk ( uplo, trans, n, k, alpha, a, lda,        &
          &beta, c, ldc ) BIND(C)
       character*1,intent(in)              :: uplo, trans
       integer,intent(in)                  :: n, k, lda, ldc
       real,intent(in)                     :: alpha, beta
       complex,intent(in)                  :: a( lda, * )
       complex,intent(in)                  :: c( ldc, * )
       !$omp  declare variant( cherk:mkl_blas_cherk_omp_offload_ilp64 ) match( construct={target variant dispatch}, device={arch(gen)} )
     end subroutine cherk

     subroutine dsyr2k ( uplo, trans, n, k, alpha, a, lda,        &
          &b, ldb, beta, c, ldc ) BIND(C)
       character*1,intent(in)             :: uplo, trans
       integer,intent(in)                 :: n, k, lda, ldb, ldc
       double precision,intent(in)        :: alpha, beta
       double precision,intent(in)        :: a( lda, * ), b( ldb, * )
       double precision,intent(inout)     :: c( ldc, * )
       !$omp  declare variant( dsyr2k:mkl_blas_dsyr2k_omp_offload_ilp64 ) match( construct={target variant dispatch}, device={arch(gen)} )
     end subroutine dsyr2k

     subroutine ssyr2k ( uplo, trans, n, k, alpha, a, lda,        &
          &b, ldb, beta, c, ldc ) BIND(C)
       character*1,intent(in)              :: uplo, trans
       integer,intent(in)                  :: n, k, lda, ldb, ldc
       real,intent(in)                     :: alpha, beta
       real,intent(in)                     :: a( lda, * ), b( ldb, * )
       real,intent(in)                     :: c( ldc, * )
       !$omp  declare variant( ssyr2k:mkl_blas_ssyr2k_omp_offload_ilp64 ) match( construct={target variant dispatch}, device={arch(gen)} )
     end subroutine ssyr2k

     subroutine zsyr2k ( uplo, trans, n, k, alpha, a, lda,        &
          &b, ldb, beta, c, ldc ) BIND(C)
       character*1,intent(in)             :: uplo, trans
       integer,intent(in)                 :: n, k, lda, ldb, ldc
       complex*16,intent(in)              :: alpha, beta
       complex*16,intent(in)              :: a( lda, * ), b( ldb, * )
       complex*16,intent(inout)           :: c( ldc, * )
       !$omp  declare variant( zsyr2k:mkl_blas_zsyr2k_omp_offload_ilp64 ) match( construct={target variant dispatch}, device={arch(gen)} )
     end subroutine zsyr2k

     subroutine csyr2k ( uplo, trans, n, k, alpha, a, lda,        &
          &b, ldb, beta, c, ldc ) BIND(C)
       character*1,intent(in)              :: uplo, trans
       integer,intent(in)                  :: n, k, lda, ldb, ldc
       complex,intent(in)                  :: alpha, beta
       complex,intent(in)                  :: a( lda, * ), b( ldb, * )
       complex,intent(in)                  :: c( ldc, * )
       !$omp  declare variant( csyr2k:mkl_blas_csyr2k_omp_offload_ilp64 ) match( construct={target variant dispatch}, device={arch(gen)} )
     end subroutine csyr2k

     subroutine zher2k ( uplo, trans, n, k, alpha, a, lda,        &
          &b, ldb, beta, c, ldc ) BIND(C)
       character*1,intent(in)             :: uplo, trans
       integer,intent(in)                 :: n, k, lda, ldb, ldc
       double precision,intent(in)        :: beta
       complex*16,intent(in)              :: alpha
       complex*16,intent(in)              :: a( lda, * ), b( ldb, * )
       complex*16,intent(inout)           :: c( ldc, * )
       !$omp  declare variant( zher2k:mkl_blas_zher2k_omp_offload_ilp64 ) match( construct={target variant dispatch}, device={arch(gen)} )
     end subroutine zher2k

     subroutine cher2k ( uplo, trans, n, k, alpha, a, lda,        &
          &b, ldb, beta, c, ldc ) BIND(C)
       character*1,intent(in)              :: uplo, trans
       integer,intent(in)                  :: n, k, lda, ldb, ldc
       real,intent(in)                     :: beta
       complex,intent(in)                  :: alpha
       complex,intent(in)                  :: a( lda, * ), b( ldb, * )
       complex,intent(in)                  :: c( ldc, * )
       !$omp  declare variant( cher2k:mkl_blas_cher2k_omp_offload_ilp64 ) match( construct={target variant dispatch}, device={arch(gen)} )
     end subroutine cher2k

     subroutine dtrmm ( side, uplo, trans, diag, m, n, alpha, a, lda,        &
          &b, ldb ) BIND(C)
       character*1,intent(in)             :: side, uplo, trans, diag
       integer,intent(in)                 :: m, n, lda, ldb
       double precision,intent(in)        :: alpha
       double precision,intent(in)        :: a( lda, * )
       double precision,intent(inout)     :: b( ldb, * )
       !$omp  declare variant( dtrmm:mkl_blas_dtrmm_omp_offload_ilp64 ) match( construct={target variant dispatch}, device={arch(gen)} )
     end subroutine dtrmm

     subroutine strmm ( side, uplo, trans, diag, m, n, alpha, a, lda,        &
          &b, ldb ) BIND(C)
       character*1,intent(in)              :: side, uplo, trans, diag
       integer,intent(in)                  :: m, n, lda, ldb
       real,intent(in)                     :: alpha
       real,intent(in)                     :: a( lda, * )
       real,intent(in)                     :: b( ldb, * )
       !$omp  declare variant( strmm:mkl_blas_strmm_omp_offload_ilp64 ) match( construct={target variant dispatch}, device={arch(gen)} )
     end subroutine strmm

     subroutine ztrmm ( side, uplo, trans, diag, m, n, alpha, a, lda,        &
          &b, ldb ) BIND(C)
       character*1,intent(in)             :: side, uplo, trans, diag
       integer,intent(in)                 :: m, n, lda, ldb
       complex*16,intent(in)              :: alpha
       complex*16,intent(in)              :: a( lda, * )
       complex*16,intent(inout)           :: b( ldb, * )
       !$omp  declare variant( ztrmm:mkl_blas_ztrmm_omp_offload_ilp64 ) match( construct={target variant dispatch}, device={arch(gen)} )
     end subroutine ztrmm

     subroutine ctrmm ( side, uplo, trans, diag, m, n, alpha, a, lda,        &
          &b, ldb ) BIND(C)
       character*1,intent(in)              :: side, uplo, trans, diag
       integer,intent(in)                  :: m, n, lda, ldb
       complex,intent(in)                  :: alpha
       complex,intent(in)                  :: a( lda, * )
       complex,intent(in)                  :: b( ldb, * )
       !$omp  declare variant( ctrmm:mkl_blas_ctrmm_omp_offload_ilp64 ) match( construct={target variant dispatch}, device={arch(gen)} )
     end subroutine ctrmm

     subroutine dtrsm ( side, uplo, trans, diag, m, n, alpha, a, lda,        &
          &b, ldb ) BIND(C)
       character*1,intent(in)             :: side, uplo, trans, diag
       integer,intent(in)                 :: m, n, lda, ldb
       double precision,intent(in)        :: alpha
       double precision,intent(in)        :: a( lda, * )
       double precision,intent(inout)     :: b( ldb, * )
       !$omp  declare variant( dtrsm:mkl_blas_dtrsm_omp_offload_ilp64 ) match( construct={target variant dispatch}, device={arch(gen)} )
     end subroutine dtrsm

     subroutine strsm ( side, uplo, trans, diag, m, n, alpha, a, lda,        &
          &b, ldb ) BIND(C)
       character*1,intent(in)              :: side, uplo, trans, diag
       integer,intent(in)                  :: m, n, lda, ldb
       real,intent(in)                     :: alpha
       real,intent(in)                     :: a( lda, * )
       real,intent(in)                     :: b( ldb, * )
       !$omp  declare variant( strsm:mkl_blas_strsm_omp_offload_ilp64 ) match( construct={target variant dispatch}, device={arch(gen)} )
     end subroutine strsm

     subroutine ztrsm ( side, uplo, trans, diag, m, n, alpha, a, lda,        &
          &b, ldb ) BIND(C)
       character*1,intent(in)             :: side, uplo, trans, diag
       integer,intent(in)                 :: m, n, lda, ldb
       complex*16,intent(in)              :: alpha
       complex*16,intent(in)              :: a( lda, * )
       complex*16,intent(inout)           :: b( ldb, * )
       !$omp  declare variant( ztrsm:mkl_blas_ztrsm_omp_offload_ilp64 ) match( construct={target variant dispatch}, device={arch(gen)} )
     end subroutine ztrsm

     subroutine ctrsm ( side, uplo, trans, diag, m, n, alpha, a, lda,        &
          &b, ldb ) BIND(C)
       character*1,intent(in)              :: side, uplo, trans, diag
       integer,intent(in)                  :: m, n, lda, ldb
       complex,intent(in)                  :: alpha
       complex,intent(in)                  :: a( lda, * )
       complex,intent(in)                  :: b( ldb, * )
       !$omp  declare variant( ctrsm:mkl_blas_ctrsm_omp_offload_ilp64 ) match( construct={target variant dispatch}, device={arch(gen)} )
     end subroutine ctrsm

     subroutine dgemmt ( uplo, transa, transb, m, n, alpha, a, lda,        &
          &b, ldb, beta, c, ldc ) BIND(C)
       character*1,intent(in)             :: uplo, transa, transb
       integer,intent(in)                 :: m, n, lda, ldb, ldc
       double precision,intent(in)        :: alpha, beta
       double precision,intent(in)        :: a( lda, * ), b( ldb, * )
       double precision,intent(inout)     :: c( ldc, * )
       !$omp  declare variant( dgemmt:mkl_blas_dgemmt_omp_offload_ilp64 ) match( construct={target variant dispatch}, device={arch(gen)} )
     end subroutine dgemmt

     subroutine sgemmt ( uplo, transa, transb, m, n, alpha, a, lda,        &
          &b, ldb, beta, c, ldc ) BIND(C)
       character*1,intent(in)              :: uplo, transa, transb
       integer,intent(in)                  :: m, n, lda, ldb, ldc
       real,intent(in)                     :: alpha, beta
       real,intent(in)                     :: a( lda, * ), b( ldb, * )
       real,intent(in)                     :: c( ldc, * )
       !$omp  declare variant( sgemmt:mkl_blas_sgemmt_omp_offload_ilp64 ) match( construct={target variant dispatch}, device={arch(gen)} )
     end subroutine sgemmt

     subroutine zgemmt ( uplo, transa, transb, m, n, alpha, a, lda,        &
          &b, ldb, beta, c, ldc ) BIND(C)
       character*1,intent(in)             :: uplo, transa, transb
       integer,intent(in)                 :: m, n, lda, ldb, ldc
       complex*16,intent(in)              :: alpha, beta
       complex*16,intent(in)              :: a( lda, * ), b( ldb, * )
       complex*16,intent(inout)           :: c( ldc, * )
       !$omp  declare variant( zgemmt:mkl_blas_zgemmt_omp_offload_ilp64 ) match( construct={target variant dispatch}, device={arch(gen)} )
     end subroutine zgemmt

     subroutine cgemmt ( uplo, transa, transb, m, n, alpha, a, lda,        &
          &b, ldb, beta, c, ldc ) BIND(C)
       character*1,intent(in)              :: uplo, transa, transb
       integer,intent(in)                  :: m, n, lda, ldb, ldc
       complex,intent(in)                  :: alpha, beta
       complex,intent(in)                  :: a( lda, * ), b( ldb, * )
       complex,intent(in)                  :: c( ldc, * )
       !$omp  declare variant( cgemmt:mkl_blas_cgemmt_omp_offload_ilp64 ) match( construct={target variant dispatch}, device={arch(gen)} )
     end subroutine cgemmt

     ! BLAS Level2

     subroutine dgbmv ( trans, m, n, kl, ku, alpha, a, lda, x, incx,   &
          &beta, y, incy ) BIND(C)
       double precision,intent(in)         :: alpha, beta
       integer,intent(in)                  :: incx, incy, kl, ku, lda, m, n
       character*1,intent(in)              :: trans
       double precision,intent(in)         :: a( lda, * ), x( * )
       double precision,intent(inout)      :: y( * )
       !$omp  declare variant( dgbmv:mkl_blas_dgbmv_omp_offload_ilp64 ) match( construct={target variant dispatch}, device={arch(gen)} )
     end subroutine dgbmv

     subroutine sgbmv ( trans, m, n, kl, ku, alpha, a, lda, x, incx,   &
          &beta, y, incy ) BIND(C)
       real,intent(in)                     :: alpha, beta
       integer,intent(in)                  :: incx, incy, kl, ku, lda, m, n
       character*1,intent(in)              :: trans
       real,intent(in)                     :: a( lda, * ), x( * )
       real,intent(inout)                  :: y( * )
       !$omp  declare variant( sgbmv:mkl_blas_sgbmv_omp_offload_ilp64 ) match( construct={target variant dispatch}, device={arch(gen)} )
     end subroutine sgbmv

     subroutine zgbmv ( trans, m, n, kl, ku, alpha, a, lda, x, incx,   &
          &beta, y, incy ) BIND(C)
       complex*16,intent(in)               :: alpha, beta
       integer,intent(in)                  :: incx, incy, kl, ku, lda, m, n
       character*1,intent(in)              :: trans
       complex*16,intent(in)               :: a( lda, * ), x( * )
       complex*16,intent(inout)            :: y( * )
       !$omp  declare variant( zgbmv:mkl_blas_zgbmv_omp_offload_ilp64 ) match( construct={target variant dispatch}, device={arch(gen)} )
     end subroutine zgbmv

     subroutine cgbmv ( trans, m, n, kl, ku, alpha, a, lda, x, incx,   &
          &beta, y, incy ) BIND(C)
       complex,intent(in)                  :: alpha, beta
       integer,intent(in)                  :: incx, incy, kl, ku, lda, m, n
       character*1,intent(in)              :: trans
       complex,intent(in)                  :: a( lda, * ), x( * )
       complex,intent(inout)               :: y( * )
       !$omp  declare variant( cgbmv:mkl_blas_cgbmv_omp_offload_ilp64 ) match( construct={target variant dispatch}, device={arch(gen)} )
     end subroutine cgbmv

     subroutine dgemv ( trans, m, n, alpha, a, lda, x, incx,   &
          &beta, y, incy ) BIND(C)
       double precision,intent(in)         :: alpha, beta
       integer,intent(in)                  :: incx, incy, lda, m, n
       character*1,intent(in)              :: trans
       double precision,intent(in)         :: a( lda, * ), x( * )
       double precision,intent(inout)      :: y( * )
       !$omp  declare variant( dgemv:mkl_blas_dgemv_omp_offload_ilp64 ) match( construct={target variant dispatch}, device={arch(gen)} )
     end subroutine dgemv

     subroutine sgemv ( trans, m, n, alpha, a, lda, x, incx,   &
          &beta, y, incy ) BIND(C)
       real,intent(in)                     :: alpha, beta
       integer,intent(in)                  :: incx, incy, lda, m, n
       character*1,intent(in)              :: trans
       real,intent(in)                     :: a( lda, * ), x( * )
       real,intent(inout)                  :: y( * )
       !$omp  declare variant( sgemv:mkl_blas_sgemv_omp_offload_ilp64 ) match( construct={target variant dispatch}, device={arch(gen)} )
     end subroutine sgemv

     subroutine zgemv ( trans, m, n, alpha, a, lda, x, incx,   &
          &beta, y, incy ) BIND(C)
       complex*16,intent(in)               :: alpha, beta
       integer,intent(in)                  :: incx, incy, lda, m, n
       character*1,intent(in)              :: trans
       complex*16,intent(in)               :: a( lda, * ), x( * )
       complex*16,intent(inout)            :: y( * )
       !$omp  declare variant( zgemv:mkl_blas_zgemv_omp_offload_ilp64 ) match( construct={target variant dispatch}, device={arch(gen)} )
     end subroutine zgemv

     subroutine cgemv ( trans, m, n, alpha, a, lda, x, incx,   &
          &beta, y, incy ) BIND(C)
       complex,intent(in)                  :: alpha, beta
       integer,intent(in)                  :: incx, incy, lda, m, n
       character*1,intent(in)              :: trans
       complex,intent(in)                  :: a( lda, * ), x( * )
       complex,intent(inout)               :: y( * )
       !$omp  declare variant( cgemv:mkl_blas_cgemv_omp_offload_ilp64 ) match( construct={target variant dispatch}, device={arch(gen)} )
     end subroutine cgemv

     subroutine dsbmv ( uplo, n, k, alpha, a, lda, x, incx,   &
          &beta, y, incy ) BIND(C)
       double precision,intent(in)         :: alpha, beta
       integer,intent(in)                  :: incx, incy, k, lda, n
       character*1,intent(in)              :: uplo
       double precision,intent(in)         :: a( lda, * ), x( * )
       double precision,intent(inout)      :: y( * )
       !$omp  declare variant( dsbmv:mkl_blas_dsbmv_omp_offload_ilp64 ) match( construct={target variant dispatch}, device={arch(gen)} )
     end subroutine dsbmv

     subroutine ssbmv ( uplo, n, k, alpha, a, lda, x, incx,   &
          &beta, y, incy ) BIND(C)
       real,intent(in)                     :: alpha, beta
       integer,intent(in)                  :: incx, incy, k, lda, n
       character*1,intent(in)              :: uplo
       real,intent(in)                     :: a( lda, * ), x( * )
       real,intent(inout)                  :: y( * )
       !$omp  declare variant( ssbmv:mkl_blas_ssbmv_omp_offload_ilp64 ) match( construct={target variant dispatch}, device={arch(gen)} )
     end subroutine ssbmv

     subroutine zhbmv ( uplo, n, k, alpha, a, lda, x, incx,   &
          &beta, y, incy ) BIND(C)
       complex*16,intent(in)               :: alpha, beta
       integer,intent(in)                  :: incx, incy, k, lda, n
       character*1,intent(in)              :: uplo
       complex*16,intent(in)               :: a( lda, * ), x( * )
       complex*16,intent(inout)            :: y( * )
       !$omp  declare variant( zhbmv:mkl_blas_zhbmv_omp_offload_ilp64 ) match( construct={target variant dispatch}, device={arch(gen)} )
     end subroutine zhbmv

     subroutine chbmv ( uplo, n, k, alpha, a, lda, x, incx,   &
          &beta, y, incy ) BIND(C)
       complex,intent(in)                  :: alpha, beta
       integer,intent(in)                  :: incx, incy, k, lda, n
       character*1,intent(in)              :: uplo
       complex,intent(in)                  :: a( lda, * ), x( * )
       complex,intent(inout)               :: y( * )
       !$omp  declare variant( chbmv:mkl_blas_chbmv_omp_offload_ilp64 ) match( construct={target variant dispatch}, device={arch(gen)} )
     end subroutine chbmv

     subroutine dsymv ( uplo, n, alpha, a, lda, x, incx,   &
          &beta, y, incy ) BIND(C)
       double precision,intent(in)         :: alpha, beta
       integer,intent(in)                  :: incx, incy, lda, n
       character*1,intent(in)              :: uplo
       double precision,intent(in)         :: a( lda, * ), x( * )
       double precision,intent(inout)      :: y( * )
       !$omp  declare variant( dsymv:mkl_blas_dsymv_omp_offload_ilp64 ) match( construct={target variant dispatch}, device={arch(gen)} )
     end subroutine dsymv

     subroutine ssymv ( uplo, n, alpha, a, lda, x, incx,   &
          &beta, y, incy ) BIND(C)
       real,intent(in)                     :: alpha, beta
       integer,intent(in)                  :: incx, incy, lda, n
       character*1,intent(in)              :: uplo
       real,intent(in)                     :: a( lda, * ), x( * )
       real,intent(inout)                  :: y( * )
       !$omp  declare variant( ssymv:mkl_blas_ssymv_omp_offload_ilp64 ) match( construct={target variant dispatch}, device={arch(gen)} )
     end subroutine ssymv

     subroutine zhemv ( uplo, n, alpha, a, lda, x, incx,   &
          &beta, y, incy ) BIND(C)
       complex*16,intent(in)               :: alpha, beta
       integer,intent(in)                  :: incx, incy, lda, n
       character*1,intent(in)              :: uplo
       complex*16,intent(in)               :: a( lda, * ), x( * )
       complex*16,intent(inout)            :: y( * )
       !$omp  declare variant( zhemv:mkl_blas_zhemv_omp_offload_ilp64 ) match( construct={target variant dispatch}, device={arch(gen)} )
     end subroutine zhemv

     subroutine chemv ( uplo, n, alpha, a, lda, x, incx,   &
          &beta, y, incy ) BIND(C)
       complex,intent(in)                  :: alpha, beta
       integer,intent(in)                  :: incx, incy, lda, n
       character*1,intent(in)              :: uplo
       complex,intent(in)                  :: a( lda, * ), x( * )
       complex,intent(inout)               :: y( * )
       !$omp  declare variant( chemv:mkl_blas_chemv_omp_offload_ilp64 ) match( construct={target variant dispatch}, device={arch(gen)} )
     end subroutine chemv

     subroutine dspmv ( uplo, n, alpha, a, x, incx,   &
          &beta, y, incy ) BIND(C)
       double precision,intent(in)         :: alpha, beta
       integer,intent(in)                  :: incx, incy, n
       character*1,intent(in)              :: uplo
       double precision,intent(in)         :: a( * ), x( * )
       double precision,intent(inout)      :: y( * )
       !$omp  declare variant( dspmv:mkl_blas_dspmv_omp_offload_ilp64 ) match( construct={target variant dispatch}, device={arch(gen)} )
     end subroutine dspmv

     subroutine sspmv ( uplo, n, alpha, a, x, incx,   &
          &beta, y, incy ) BIND(C)
       real,intent(in)                     :: alpha, beta
       integer,intent(in)                  :: incx, incy, n
       character*1,intent(in)              :: uplo
       real,intent(in)                     :: a( * ), x( * )
       real,intent(inout)                  :: y( * )
       !$omp  declare variant( sspmv:mkl_blas_sspmv_omp_offload_ilp64 ) match( construct={target variant dispatch}, device={arch(gen)} )
     end subroutine sspmv

     subroutine zhpmv ( uplo, n, alpha, a, x, incx,   &
          &beta, y, incy ) BIND(C)
       complex*16,intent(in)               :: alpha, beta
       integer,intent(in)                  :: incx, incy, n
       character*1,intent(in)              :: uplo
       complex*16,intent(in)               :: a( * ), x( * )
       complex*16,intent(inout)            :: y( * )
       !$omp  declare variant( zhpmv:mkl_blas_zhpmv_omp_offload_ilp64 ) match( construct={target variant dispatch}, device={arch(gen)} )
     end subroutine zhpmv

     subroutine chpmv ( uplo, n, alpha, a, x, incx,   &
          &beta, y, incy ) BIND(C)
       complex,intent(in)                  :: alpha, beta
       integer,intent(in)                  :: incx, incy, n
       character*1,intent(in)              :: uplo
       complex,intent(in)                  :: a( * ), x( * )
       complex,intent(inout)               :: y( * )
       !$omp  declare variant( chpmv:mkl_blas_chpmv_omp_offload_ilp64 ) match( construct={target variant dispatch}, device={arch(gen)} )
     end subroutine chpmv

     subroutine dger ( m, n, alpha, x, incx,   &
          &y, incy, a, lda ) BIND(C)
       double precision,intent(in)         :: alpha
       integer,intent(in)                  :: incx, incy, lda, m, n
       double precision,intent(in)         :: y( * ), x( * )
       double precision,intent(inout)      :: a( lda, * )
       !$omp  declare variant( dger:mkl_blas_dger_omp_offload_ilp64 ) match( construct={target variant dispatch}, device={arch(gen)} )
     end subroutine dger

     subroutine sger ( m, n, alpha, x, incx,   &
          &y, incy, a, lda ) BIND(C)
       real,intent(in)                     :: alpha
       integer,intent(in)                  :: incx, incy, lda, m, n
       real,intent(in)                     :: y( * ), x( * )
       real,intent(inout)                  :: a( lda, * )
       !$omp  declare variant( sger:mkl_blas_sger_omp_offload_ilp64 ) match( construct={target variant dispatch}, device={arch(gen)} )
     end subroutine sger

     subroutine zgerc ( m, n, alpha, x, incx,   &
          &y, incy, a, lda ) BIND(C)
       complex*16,intent(in)               :: alpha
       integer,intent(in)                  :: incx, incy, lda, m, n
       complex*16,intent(in)               :: y( * ), x( * )
       complex*16,intent(inout)            :: a( lda, * )
       !$omp  declare variant( zgerc:mkl_blas_zgerc_omp_offload_ilp64 ) match( construct={target variant dispatch}, device={arch(gen)} )
     end subroutine zgerc

     subroutine cgerc ( m, n, alpha, x, incx,   &
          &y, incy, a, lda ) BIND(C)
       complex,intent(in)                  :: alpha
       integer,intent(in)                  :: incx, incy, lda, m, n
       complex,intent(in)                  :: y( * ), x( * )
       complex,intent(inout)               :: a( lda, * )
       !$omp  declare variant( cgerc:mkl_blas_cgerc_omp_offload_ilp64 ) match( construct={target variant dispatch}, device={arch(gen)} )
     end subroutine cgerc

     subroutine zgeru ( m, n, alpha, x, incx,   &
          &y, incy, a, lda ) BIND(C)
       complex*16,intent(in)               :: alpha
       integer,intent(in)                  :: incx, incy, lda, m, n
       complex*16,intent(in)               :: y( * ), x( * )
       complex*16,intent(inout)            :: a( lda, * )
       !$omp  declare variant( zgeru:mkl_blas_zgeru_omp_offload_ilp64 ) match( construct={target variant dispatch}, device={arch(gen)} )
     end subroutine zgeru

     subroutine cgeru ( m, n, alpha, x, incx,   &
          &y, incy, a, lda ) BIND(C)
       complex,intent(in)                  :: alpha
       integer,intent(in)                  :: incx, incy, lda, m, n
       complex,intent(in)                  :: y( * ), x( * )
       complex,intent(inout)               :: a( lda, * )
       !$omp  declare variant( cgeru:mkl_blas_cgeru_omp_offload_ilp64 ) match( construct={target variant dispatch}, device={arch(gen)} )
     end subroutine cgeru

     subroutine dsyr ( uplo, n, alpha, x, incx,   &
          &a, lda ) BIND(C)
       character*1,intent(in)              :: uplo
       double precision,intent(in)         :: alpha
       integer,intent(in)                  :: incx, lda, n
       double precision,intent(in)         :: x( * )
       double precision,intent(inout)      :: a( lda, * )
       !$omp  declare variant( dsyr:mkl_blas_dsyr_omp_offload_ilp64 ) match( construct={target variant dispatch}, device={arch(gen)} )
     end subroutine dsyr

     subroutine ssyr ( uplo, n, alpha, x, incx,   &
          &a, lda ) BIND(C)
       character*1,intent(in)              :: uplo
       real,intent(in)                     :: alpha
       integer,intent(in)                  :: incx, lda, n
       real,intent(in)                     :: x( * )
       real,intent(inout)                  :: a( lda, * )
       !$omp  declare variant( ssyr:mkl_blas_ssyr_omp_offload_ilp64 ) match( construct={target variant dispatch}, device={arch(gen)} )
     end subroutine ssyr

     subroutine zher ( uplo, n, alpha, x, incx,   &
          &a, lda ) BIND(C)
       character*1,intent(in)              :: uplo
       double precision,intent(in)         :: alpha
       integer,intent(in)                  :: incx, lda, n
       complex*16,intent(in)               :: x( * )
       complex*16,intent(inout)            :: a( lda, * )
       !$omp  declare variant( zher:mkl_blas_zher_omp_offload_ilp64 ) match( construct={target variant dispatch}, device={arch(gen)} )
     end subroutine zher

     subroutine cher ( uplo, n, alpha, x, incx,   &
          &a, lda ) BIND(C)
       character*1,intent(in)              :: uplo
       real,intent(in)                     :: alpha
       integer,intent(in)                  :: incx, lda, n
       complex,intent(in)                  :: x( * )
       complex,intent(inout)               :: a( lda, * )
       !$omp  declare variant( cher:mkl_blas_cher_omp_offload_ilp64 ) match( construct={target variant dispatch}, device={arch(gen)} )
     end subroutine cher

     subroutine dspr ( uplo, n, alpha, x, incx,   &
          &a ) BIND(C)
       character*1,intent(in)              :: uplo
       double precision,intent(in)         :: alpha
       integer,intent(in)                  :: incx, n
       double precision,intent(in)         :: x( * )
       double precision,intent(inout)      :: a( * )
       !$omp  declare variant( dspr:mkl_blas_dspr_omp_offload_ilp64 ) match( construct={target variant dispatch}, device={arch(gen)} )
     end subroutine dspr

     subroutine sspr ( uplo, n, alpha, x, incx,   &
          &a ) BIND(C)
       character*1,intent(in)              :: uplo
       real,intent(in)                     :: alpha
       integer,intent(in)                  :: incx, n
       real,intent(in)                     :: x( * )
       real,intent(inout)                  :: a( * )
       !$omp  declare variant( sspr:mkl_blas_sspr_omp_offload_ilp64 ) match( construct={target variant dispatch}, device={arch(gen)} )
     end subroutine sspr

     subroutine zhpr ( uplo, n, alpha, x, incx,   &
          &a ) BIND(C)
       character*1,intent(in)              :: uplo
       double precision,intent(in)         :: alpha
       integer,intent(in)                  :: incx, n
       complex*16,intent(in)               :: x( * )
       complex*16,intent(inout)            :: a( * )
       !$omp  declare variant( zhpr:mkl_blas_zhpr_omp_offload_ilp64 ) match( construct={target variant dispatch}, device={arch(gen)} )
     end subroutine zhpr

     subroutine chpr ( uplo, n, alpha, x, incx,   &
          &a ) BIND(C)
       character*1,intent(in)              :: uplo
       real,intent(in)                     :: alpha
       integer,intent(in)                  :: incx, n
       complex,intent(in)                  :: x( * )
       complex,intent(inout)               :: a( * )
       !$omp  declare variant( chpr:mkl_blas_chpr_omp_offload_ilp64 ) match( construct={target variant dispatch}, device={arch(gen)} )
     end subroutine chpr

     subroutine dsyr2 ( uplo, n, alpha, x, incx, y, incy,   &
          &a, lda ) BIND(C)
       character*1,intent(in)              :: uplo
       double precision,intent(in)         :: alpha
       integer,intent(in)                  :: incx, incy, lda, n
       double precision,intent(in)         :: x( * ), y( * )
       double precision,intent(inout)      :: a( lda, * )
       !$omp  declare variant( dsyr2:mkl_blas_dsyr2_omp_offload_ilp64 ) match( construct={target variant dispatch}, device={arch(gen)} )
     end subroutine dsyr2

     subroutine ssyr2 ( uplo, n, alpha, x, incx, y, incy,   &
          &a, lda ) BIND(C)
       character*1,intent(in)              :: uplo
       real,intent(in)                     :: alpha
       integer,intent(in)                  :: incx, incy, lda, n
       real,intent(in)                     :: x( * ), y( * )
       real,intent(inout)                  :: a( lda, * )
       !$omp  declare variant( ssyr2:mkl_blas_ssyr2_omp_offload_ilp64 ) match( construct={target variant dispatch}, device={arch(gen)} )
     end subroutine ssyr2

     subroutine zher2 ( uplo, n, alpha, x, incx, y, incy,   &
          &a, lda ) BIND(C)
       character*1,intent(in)              :: uplo
       complex*16,intent(in)               :: alpha
       integer,intent(in)                  :: incx, incy, lda, n
       complex*16,intent(in)               :: x( * ), y( * )
       complex*16,intent(inout)            :: a( lda, * )
       !$omp  declare variant( zher2:mkl_blas_zher2_omp_offload_ilp64 ) match( construct={target variant dispatch}, device={arch(gen)} )
     end subroutine zher2

     subroutine cher2 ( uplo, n, alpha, x, incx, y, incy,   &
          &a, lda ) BIND(C)
       character*1,intent(in)              :: uplo
       complex,intent(in)                  :: alpha
       integer,intent(in)                  :: incx, incy, lda, n
       complex,intent(in)                  :: x( * ), y( * )
       complex,intent(inout)               :: a( lda, * )
       !$omp  declare variant( cher2:mkl_blas_cher2_omp_offload_ilp64 ) match( construct={target variant dispatch}, device={arch(gen)} )
     end subroutine cher2

     subroutine dspr2 ( uplo, n, alpha, x, incx, y, incy,   &
          &a ) BIND(C)
       character*1,intent(in)              :: uplo
       double precision,intent(in)         :: alpha
       integer,intent(in)                  :: incx, incy, n
       double precision,intent(in)         :: x( * ), y( * )
       double precision,intent(inout)      :: a( * )
       !$omp  declare variant( dspr2:mkl_blas_dspr2_omp_offload_ilp64 ) match( construct={target variant dispatch}, device={arch(gen)} )
     end subroutine dspr2

     subroutine sspr2 ( uplo, n, alpha, x, incx, y, incy,   &
          &a ) BIND(C)
       character*1,intent(in)              :: uplo
       real,intent(in)                     :: alpha
       integer,intent(in)                  :: incx, incy, n
       real,intent(in)                     :: x( * ), y( * )
       real,intent(inout)                  :: a( * )
       !$omp  declare variant( sspr2:mkl_blas_sspr2_omp_offload_ilp64 ) match( construct={target variant dispatch}, device={arch(gen)} )
     end subroutine sspr2

     subroutine zhpr2 ( uplo, n, alpha, x, incx, y, incy,   &
          &a ) BIND(C)
       character*1,intent(in)              :: uplo
       complex*16,intent(in)               :: alpha
       integer,intent(in)                  :: incx, incy, n
       complex*16,intent(in)               :: x( * ), y( * )
       complex*16,intent(inout)            :: a( * )
       !$omp  declare variant( zhpr2:mkl_blas_zhpr2_omp_offload_ilp64 ) match( construct={target variant dispatch}, device={arch(gen)} )
     end subroutine zhpr2

     subroutine chpr2 ( uplo, n, alpha, x, incx, y, incy,   &
          &a ) BIND(C)
       character*1,intent(in)              :: uplo
       complex,intent(in)                  :: alpha
       integer,intent(in)                  :: incx, incy, n
       complex,intent(in)                  :: x( * ), y( * )
       complex,intent(inout)               :: a( * )
       !$omp  declare variant( chpr2:mkl_blas_chpr2_omp_offload_ilp64 ) match( construct={target variant dispatch}, device={arch(gen)} )
     end subroutine chpr2

     subroutine dtbmv ( uplo, trans, diag, n, k, a, lda, x, incx ) BIND(C)
       integer,intent(in)                  :: incx, k, lda, n
       character*1,intent(in)              :: uplo, trans, diag
       double precision,intent(in)         :: a( lda, * )
       double precision,intent(inout)      :: x( * )
       !$omp  declare variant( dtbmv:mkl_blas_dtbmv_omp_offload_ilp64 ) match( construct={target variant dispatch}, device={arch(gen)} )
     end subroutine dtbmv

     subroutine stbmv ( uplo, trans, diag, n, k, a, lda, x, incx ) BIND(C)
       integer,intent(in)                  :: incx, k, lda, n
       character*1,intent(in)              :: uplo, trans, diag
       real,intent(in)                     :: a( lda, * )
       real,intent(inout)                  :: x( * )
       !$omp  declare variant( stbmv:mkl_blas_stbmv_omp_offload_ilp64 ) match( construct={target variant dispatch}, device={arch(gen)} )
     end subroutine stbmv

     subroutine ztbmv ( uplo, trans, diag, n, k, a, lda, x, incx ) BIND(C)
       integer,intent(in)                  :: incx, k, lda, n
       character*1,intent(in)              :: uplo, trans, diag
       complex*16,intent(in)               :: a( lda, * )
       complex*16,intent(inout)            :: x( * )
       !$omp  declare variant( ztbmv:mkl_blas_ztbmv_omp_offload_ilp64 ) match( construct={target variant dispatch}, device={arch(gen)} )
     end subroutine ztbmv

     subroutine ctbmv ( uplo, trans, diag, n, k, a, lda, x, incx ) BIND(C)
       integer,intent(in)                  :: incx, k, lda, n
       character*1,intent(in)              :: uplo, trans, diag
       complex,intent(in)                  :: a( lda, * )
       complex,intent(inout)               :: x( * )
       !$omp  declare variant( ctbmv:mkl_blas_ctbmv_omp_offload_ilp64 ) match( construct={target variant dispatch}, device={arch(gen)} )
     end subroutine ctbmv

     subroutine dtrmv ( uplo, trans, diag, n, a, lda, x, incx ) BIND(C)
       integer,intent(in)                  :: incx, lda, n
       character*1,intent(in)              :: uplo, trans, diag
       double precision,intent(in)         :: a( lda, * )
       double precision,intent(inout)      :: x( * )
       !$omp  declare variant( dtrmv:mkl_blas_dtrmv_omp_offload_ilp64 ) match( construct={target variant dispatch}, device={arch(gen)} )
     end subroutine dtrmv

     subroutine strmv ( uplo, trans, diag, n, a, lda, x, incx ) BIND(C)
       integer,intent(in)                  :: incx, lda, n
       character*1,intent(in)              :: uplo, trans, diag
       real,intent(in)                     :: a( lda, * )
       real,intent(inout)                  :: x( * )
       !$omp  declare variant( strmv:mkl_blas_strmv_omp_offload_ilp64 ) match( construct={target variant dispatch}, device={arch(gen)} )
     end subroutine strmv

     subroutine ztrmv ( uplo, trans, diag, n, a, lda, x, incx ) BIND(C)
       integer,intent(in)                  :: incx, lda, n
       character*1,intent(in)              :: uplo, trans, diag
       complex*16,intent(in)               :: a( lda, * )
       complex*16,intent(inout)            :: x( * )
       !$omp  declare variant( ztrmv:mkl_blas_ztrmv_omp_offload_ilp64 ) match( construct={target variant dispatch}, device={arch(gen)} )
     end subroutine ztrmv

     subroutine ctrmv ( uplo, trans, diag, n, a, lda, x, incx ) BIND(C)
       integer,intent(in)                  :: incx, lda, n
       character*1,intent(in)              :: uplo, trans, diag
       complex,intent(in)                  :: a( lda, * )
       complex,intent(inout)               :: x( * )
       !$omp  declare variant( ctrmv:mkl_blas_ctrmv_omp_offload_ilp64 ) match( construct={target variant dispatch}, device={arch(gen)} )
     end subroutine ctrmv

     subroutine dtpmv ( uplo, trans, diag, n, a, x, incx ) BIND(C)
       integer,intent(in)                  :: incx, n
       character*1,intent(in)              :: uplo, trans, diag
       double precision,intent(in)         :: a( * )
       double precision,intent(inout)      :: x( * )
       !$omp  declare variant( dtpmv:mkl_blas_dtpmv_omp_offload_ilp64 ) match( construct={target variant dispatch}, device={arch(gen)} )
     end subroutine dtpmv

     subroutine stpmv ( uplo, trans, diag, n, a, x, incx ) BIND(C)
       integer,intent(in)                  :: incx, n
       character*1,intent(in)              :: uplo, trans, diag
       real,intent(in)                     :: a( * )
       real,intent(inout)                  :: x( * )
       !$omp  declare variant( stpmv:mkl_blas_stpmv_omp_offload_ilp64 ) match( construct={target variant dispatch}, device={arch(gen)} )
     end subroutine stpmv

     subroutine ztpmv ( uplo, trans, diag, n, a, x, incx ) BIND(C)
       integer,intent(in)                  :: incx, n
       character*1,intent(in)              :: uplo, trans, diag
       complex*16,intent(in)               :: a( * )
       complex*16,intent(inout)            :: x( * )
       !$omp  declare variant( ztpmv:mkl_blas_ztpmv_omp_offload_ilp64 ) match( construct={target variant dispatch}, device={arch(gen)} )
     end subroutine ztpmv

     subroutine ctpmv ( uplo, trans, diag, n, a, x, incx ) BIND(C)
       integer,intent(in)                  :: incx, n
       character*1,intent(in)              :: uplo, trans, diag
       complex,intent(in)                  :: a( * )
       complex,intent(inout)               :: x( * )
       !$omp  declare variant( ctpmv:mkl_blas_ctpmv_omp_offload_ilp64 ) match( construct={target variant dispatch}, device={arch(gen)} )
     end subroutine ctpmv
     
     subroutine dtbsv ( uplo, trans, diag, n, k, a, lda, x, incx ) BIND(C)
       integer,intent(in)                  :: incx, k, lda, n
       character*1,intent(in)              :: uplo, trans, diag
       double precision,intent(in)         :: a( lda, * )
       double precision,intent(inout)      :: x( * )
       !$omp  declare variant( dtbsv:mkl_blas_dtbsv_omp_offload_ilp64 ) match( construct={target variant dispatch}, device={arch(gen)} )
     end subroutine dtbsv

     subroutine stbsv ( uplo, trans, diag, n, k, a, lda, x, incx ) BIND(C)
       integer,intent(in)                  :: incx, k, lda, n
       character*1,intent(in)              :: uplo, trans, diag
       real,intent(in)                     :: a( lda, * )
       real,intent(inout)                  :: x( * )
       !$omp  declare variant( stbsv:mkl_blas_stbsv_omp_offload_ilp64 ) match( construct={target variant dispatch}, device={arch(gen)} )
     end subroutine stbsv

     subroutine ztbsv ( uplo, trans, diag, n, k, a, lda, x, incx ) BIND(C)
       integer,intent(in)                  :: incx, k, lda, n
       character*1,intent(in)              :: uplo, trans, diag
       complex*16,intent(in)               :: a( lda, * )
       complex*16,intent(inout)            :: x( * )
       !$omp  declare variant( ztbsv:mkl_blas_ztbsv_omp_offload_ilp64 ) match( construct={target variant dispatch}, device={arch(gen)} )
     end subroutine ztbsv

     subroutine ctbsv ( uplo, trans, diag, n, k, a, lda, x, incx ) BIND(C)
       integer,intent(in)                  :: incx, k, lda, n
       character*1,intent(in)              :: uplo, trans, diag
       complex,intent(in)                  :: a( lda, * )
       complex,intent(inout)               :: x( * )
       !$omp  declare variant( ctbsv:mkl_blas_ctbsv_omp_offload_ilp64 ) match( construct={target variant dispatch}, device={arch(gen)} )
     end subroutine ctbsv

     subroutine dtrsv ( uplo, trans, diag, n, a, lda, x, incx ) BIND(C)
       integer,intent(in)                  :: incx, lda, n
       character*1,intent(in)              :: uplo, trans, diag
       double precision,intent(in)         :: a( lda, * )
       double precision,intent(inout)      :: x( * )
       !$omp  declare variant( dtrsv:mkl_blas_dtrsv_omp_offload_ilp64 ) match( construct={target variant dispatch}, device={arch(gen)} )
     end subroutine dtrsv

     subroutine strsv ( uplo, trans, diag, n, a, lda, x, incx ) BIND(C)
       integer,intent(in)                  :: incx, lda, n
       character*1,intent(in)              :: uplo, trans, diag
       real,intent(in)                     :: a( lda, * )
       real,intent(inout)                  :: x( * )
       !$omp  declare variant( strsv:mkl_blas_strsv_omp_offload_ilp64 ) match( construct={target variant dispatch}, device={arch(gen)} )
     end subroutine strsv

     subroutine ztrsv ( uplo, trans, diag, n, a, lda, x, incx ) BIND(C)
       integer,intent(in)                  :: incx, lda, n
       character*1,intent(in)              :: uplo, trans, diag
       complex*16,intent(in)               :: a( lda, * )
       complex*16,intent(inout)            :: x( * )
       !$omp  declare variant( ztrsv:mkl_blas_ztrsv_omp_offload_ilp64 ) match( construct={target variant dispatch}, device={arch(gen)} )
     end subroutine ztrsv

     subroutine ctrsv ( uplo, trans, diag, n, a, lda, x, incx ) BIND(C)
       integer,intent(in)                  :: incx, lda, n
       character*1,intent(in)              :: uplo, trans, diag
       complex,intent(in)                  :: a( lda, * )
       complex,intent(inout)               :: x( * )
       !$omp  declare variant( ctrsv:mkl_blas_ctrsv_omp_offload_ilp64 ) match( construct={target variant dispatch}, device={arch(gen)} )
     end subroutine ctrsv

     subroutine dtpsv ( uplo, trans, diag, n, a, x, incx ) BIND(C)
       integer,intent(in)                  :: incx, n
       character*1,intent(in)              :: uplo, trans, diag
       double precision,intent(in)         :: a( * )
       double precision,intent(inout)      :: x( * )
       !$omp  declare variant( dtpsv:mkl_blas_dtpsv_omp_offload_ilp64 ) match( construct={target variant dispatch}, device={arch(gen)} )
     end subroutine dtpsv

     subroutine stpsv ( uplo, trans, diag, n, a, x, incx ) BIND(C)
       integer,intent(in)                  :: incx, n
       character*1,intent(in)              :: uplo, trans, diag
       real,intent(in)                     :: a( * )
       real,intent(inout)                  :: x( * )
       !$omp  declare variant( stpsv:mkl_blas_stpsv_omp_offload_ilp64 ) match( construct={target variant dispatch}, device={arch(gen)} )
     end subroutine stpsv

     subroutine ztpsv ( uplo, trans, diag, n, a, x, incx ) BIND(C)
       integer,intent(in)                  :: incx, n
       character*1,intent(in)              :: uplo, trans, diag
       complex*16,intent(in)               :: a( * )
       complex*16,intent(inout)            :: x( * )
       !$omp  declare variant( ztpsv:mkl_blas_ztpsv_omp_offload_ilp64 ) match( construct={target variant dispatch}, device={arch(gen)} )
     end subroutine ztpsv

     subroutine ctpsv ( uplo, trans, diag, n, a, x, incx ) BIND(C)
       integer,intent(in)                  :: incx, n
       character*1,intent(in)              :: uplo, trans, diag
       complex,intent(in)                  :: a( * )
       complex,intent(inout)               :: x( * )
       !$omp  declare variant( ctpsv:mkl_blas_ctpsv_omp_offload_ilp64 ) match( construct={target variant dispatch}, device={arch(gen)} )
     end subroutine ctpsv

     ! BLAS Level1

     function dasum ( n, x, incx ) BIND(C)
       double precision                    :: dasum
       integer,intent(in)                  :: incx, n
       double precision,intent(in)         :: x( * )
       !$omp  declare variant( dasum:mkl_blas_dasum_omp_offload_ilp64 ) match( construct={target variant dispatch}, device={arch(gen)} )
     end function dasum

     function sasum ( n, x, incx ) BIND(C)
       real                                :: sasum
       integer,intent(in)                  :: incx, n
       real,intent(in)                     :: x( * )
       !$omp  declare variant( sasum:mkl_blas_sasum_omp_offload_ilp64 ) match( construct={target variant dispatch}, device={arch(gen)} )
     end function sasum

     function dzasum ( n, x, incx ) BIND(C)
       double precision                    :: dzasum
       integer,intent(in)                  :: incx, n
       complex*16,intent(in)               :: x( * )
       !$omp  declare variant( dzasum:mkl_blas_dzasum_omp_offload_ilp64 ) match( construct={target variant dispatch}, device={arch(gen)} )
     end function dzasum

     function scasum ( n, x, incx ) BIND(C)
       real                                :: scasum
       integer,intent(in)                  :: incx, n
       complex,intent(in)                  :: x( * )
       !$omp  declare variant( scasum:mkl_blas_scasum_omp_offload_ilp64 ) match( construct={target variant dispatch}, device={arch(gen)} )
     end function scasum

     function dnrm2 ( n, x, incx ) BIND(C)
       double precision                    :: dnrm2
       integer,intent(in)                  :: incx, n
       double precision,intent(in)         :: x( * )
       !$omp  declare variant( dnrm2:mkl_blas_dnrm2_omp_offload_ilp64 ) match( construct={target variant dispatch}, device={arch(gen)} )
     end function dnrm2

     function snrm2 ( n, x, incx ) BIND(C)
       real                                :: snrm2
       integer,intent(in)                  :: incx, n
       real,intent(in)                     :: x( * )
       !$omp  declare variant( snrm2:mkl_blas_snrm2_omp_offload_ilp64 ) match( construct={target variant dispatch}, device={arch(gen)} )
     end function snrm2

     function dznrm2 ( n, x, incx ) BIND(C)
       double precision                    :: dznrm2
       integer,intent(in)                  :: incx, n
       complex*16,intent(in)               :: x( * )
       !$omp  declare variant( dznrm2:mkl_blas_dznrm2_omp_offload_ilp64 ) match( construct={target variant dispatch}, device={arch(gen)} )
     end function dznrm2

     function scnrm2 ( n, x, incx ) BIND(C)
       real                                :: scnrm2
       integer,intent(in)                  :: incx, n
       complex,intent(in)                  :: x( * )
       !$omp  declare variant( scnrm2:mkl_blas_scnrm2_omp_offload_ilp64 ) match( construct={target variant dispatch}, device={arch(gen)} )
     end function scnrm2

     subroutine daxpy ( n, alpha, x, incx,   &
          &y, incy ) BIND(C)
       double precision,intent(in)         :: alpha
       integer,intent(in)                  :: incx, incy, n
       double precision,intent(in)         :: x( * )
       double precision,intent(inout)      :: y( * )
       !$omp  declare variant( daxpy:mkl_blas_daxpy_omp_offload_ilp64 ) match( construct={target variant dispatch}, device={arch(gen)} )
     end subroutine daxpy

     subroutine saxpy ( n, alpha, x, incx,   &
          &y, incy ) BIND(C)
       real,intent(in)                     :: alpha
       integer,intent(in)                  :: incx, incy, n
       real,intent(in)                     :: x( * )
       real,intent(inout)                  :: y( * )
       !$omp  declare variant( saxpy:mkl_blas_saxpy_omp_offload_ilp64 ) match( construct={target variant dispatch}, device={arch(gen)} )
     end subroutine saxpy

     subroutine zaxpy ( n, alpha, x, incx,   &
          &y, incy ) BIND(C)
       complex*16,intent(in)               :: alpha
       integer,intent(in)                  :: incx, incy, n
       complex*16,intent(in)               :: x( * )
       complex*16,intent(inout)            :: y( * )
       !$omp  declare variant( zaxpy:mkl_blas_zaxpy_omp_offload_ilp64 ) match( construct={target variant dispatch}, device={arch(gen)} )
     end subroutine zaxpy

     subroutine caxpy ( n, alpha, x, incx,   &
          &y, incy ) BIND(C)
       complex,intent(in)                  :: alpha
       integer,intent(in)                  :: incx, incy, n
       complex,intent(in)                  :: x( * )
       complex,intent(inout)               :: y( * )
       !$omp  declare variant( caxpy:mkl_blas_caxpy_omp_offload_ilp64 ) match( construct={target variant dispatch}, device={arch(gen)} )
     end subroutine caxpy

     subroutine daxpby ( n, alpha, x, incx,   &
          &beta, y, incy ) BIND(C)
       double precision,intent(in)         :: alpha, beta
       integer,intent(in)                  :: incx, incy, n
       double precision,intent(in)         :: x( * )
       double precision,intent(inout)      :: y( * )
       !$omp  declare variant( daxpby:mkl_blas_daxpby_omp_offload_ilp64 ) match( construct={target variant dispatch}, device={arch(gen)} )
     end subroutine daxpby

     subroutine saxpby ( n, alpha, x, incx,   &
          &beta, y, incy ) BIND(C)
       real,intent(in)                     :: alpha, beta
       integer,intent(in)                  :: incx, incy, n
       real,intent(in)                     :: x( * )
       real,intent(inout)                  :: y( * )
       !$omp  declare variant( saxpby:mkl_blas_saxpby_omp_offload_ilp64 ) match( construct={target variant dispatch}, device={arch(gen)} )
     end subroutine saxpby

     subroutine zaxpby ( n, alpha, x, incx,   &
          &beta, y, incy ) BIND(C)
       complex*16,intent(in)               :: alpha, beta
       integer,intent(in)                  :: incx, incy, n
       complex*16,intent(in)               :: x( * )
       complex*16,intent(inout)            :: y( * )
       !$omp  declare variant( zaxpby:mkl_blas_zaxpby_omp_offload_ilp64 ) match( construct={target variant dispatch}, device={arch(gen)} )
     end subroutine zaxpby

     subroutine caxpby ( n, alpha, x, incx,   &
          &beta, y, incy ) BIND(C)
       complex,intent(in)                  :: alpha, beta
       integer,intent(in)                  :: incx, incy, n
       complex,intent(in)                  :: x( * )
       complex,intent(inout)               :: y( * )
       !$omp  declare variant( caxpby:mkl_blas_caxpby_omp_offload_ilp64 ) match( construct={target variant dispatch}, device={arch(gen)} )
     end subroutine caxpby

     subroutine dcopy ( n, x, incx,   &
          &y, incy ) BIND(C)
       integer,intent(in)                  :: incx, incy, n
       double precision,intent(in)         :: x( * )
       double precision,intent(inout)      :: y( * )
       !$omp  declare variant( dcopy:mkl_blas_dcopy_omp_offload_ilp64 ) match( construct={target variant dispatch}, device={arch(gen)} )
     end subroutine dcopy

     subroutine scopy ( n, x, incx,   &
          &y, incy ) BIND(C)
       integer,intent(in)                  :: incx, incy, n
       real,intent(in)                     :: x( * )
       real,intent(inout)                  :: y( * )
       !$omp  declare variant( scopy:mkl_blas_scopy_omp_offload_ilp64 ) match( construct={target variant dispatch}, device={arch(gen)} )
     end subroutine scopy

     subroutine zcopy ( n, x, incx,   &
          &y, incy ) BIND(C)
       integer,intent(in)                  :: incx, incy, n
       complex*16,intent(in)               :: x( * )
       complex*16,intent(inout)            :: y( * )
       !$omp  declare variant( zcopy:mkl_blas_zcopy_omp_offload_ilp64 ) match( construct={target variant dispatch}, device={arch(gen)} )
     end subroutine zcopy

     subroutine ccopy ( n, x, incx,   &
          &y, incy ) BIND(C)
       integer,intent(in)                  :: incx, incy, n
       complex,intent(in)                  :: x( * )
       complex,intent(inout)               :: y( * )
       !$omp  declare variant( ccopy:mkl_blas_ccopy_omp_offload_ilp64 ) match( construct={target variant dispatch}, device={arch(gen)} )
     end subroutine ccopy

     function ddot ( n, x, incx, y, incy ) BIND(C)
       double precision                    :: ddot
       integer,intent(in)                  :: incx, incy, n
       double precision,intent(in)         :: x( * ), y( * )
       !$omp  declare variant( ddot:mkl_blas_ddot_omp_offload_ilp64 ) match( construct={target variant dispatch}, device={arch(gen)} )
     end function ddot

     function sdot ( n, x, incx, y, incy ) BIND(C)
       real                                :: sdot
       integer,intent(in)                  :: incx, incy, n
       real,intent(in)                     :: x( * ), y( * )
       !$omp  declare variant( sdot:mkl_blas_sdot_omp_offload_ilp64 ) match( construct={target variant dispatch}, device={arch(gen)} )
     end function sdot

     function dsdot ( n, x, incx, y, incy ) BIND(C)
       double precision                    :: dsdot
       integer,intent(in)                  :: incx, incy, n
       real,intent(in)                     :: x( * ), y( * )
       !$omp  declare variant( dsdot:mkl_blas_dsdot_omp_offload_ilp64 ) match( construct={target variant dispatch}, device={arch(gen)} )
     end function dsdot

     function sdsdot ( n, sb, x, incx, y, incy ) BIND(C)
       real                                :: sdsdot
       integer,intent(in)                  :: incx, incy, n
       real,intent(in)                     :: sb, x( * ), y( * )
       !$omp  declare variant( sdsdot:mkl_blas_sdsdot_omp_offload_ilp64 ) match( construct={target variant dispatch}, device={arch(gen)} )
     end function sdsdot

     function zdotc ( n, x, incx, y, incy )
       complex*16                          :: zdotc
       integer,intent(in)                  :: incx, incy, n
       complex*16,intent(in)               :: x( * ), y( * )
       !$omp  declare variant( zdotc:mkl_blas_zdotc_omp_offload_ilp64 ) match( construct={target variant dispatch}, device={arch(gen)} )
     end function zdotc

     function cdotc ( n, x, incx, y, incy )
       complex                             :: cdotc
       integer,intent(in)                  :: incx, incy, n
       complex,intent(in)                  :: x( * ), y( * )
       !$omp  declare variant( cdotc:mkl_blas_cdotc_omp_offload_ilp64 ) match( construct={target variant dispatch}, device={arch(gen)} )
     end function cdotc

     function zdotu ( n, x, incx, y, incy )
       complex*16                          :: zdotu
       integer,intent(in)                  :: incx, incy, n
       complex*16,intent(in)               :: x( * ), y( * )
       !$omp  declare variant( zdotu:mkl_blas_zdotu_omp_offload_ilp64 ) match( construct={target variant dispatch}, device={arch(gen)} )
     end function zdotu

     function cdotu ( n, x, incx, y, incy )
       complex                             :: cdotu
       integer,intent(in)                  :: incx, incy, n
       complex,intent(in)                  :: x( * ), y( * )
       !$omp  declare variant( cdotu:mkl_blas_cdotu_omp_offload_ilp64 ) match( construct={target variant dispatch}, device={arch(gen)} )
     end function cdotu

     subroutine drot ( n, x, incx,   &
          &y, incy, c, s ) BIND(C)
       double precision,intent(in)         :: c, s         
       integer,intent(in)                  :: incx, incy, n
       double precision,intent(inout)      :: x( * ), y( * )
       !$omp  declare variant( drot:mkl_blas_drot_omp_offload_ilp64 ) match( construct={target variant dispatch}, device={arch(gen)} )
     end subroutine drot

     subroutine srot ( n, x, incx,   &
          &y, incy, c, s ) BIND(C)
       real,intent(in)                     :: c, s         
       integer,intent(in)                  :: incx, incy, n
       real,intent(inout)                  :: x( * ), y( * )
       !$omp  declare variant( srot:mkl_blas_srot_omp_offload_ilp64 ) match( construct={target variant dispatch}, device={arch(gen)} )
     end subroutine srot

     subroutine zdrot ( n, x, incx,   &
          &y, incy, c, s ) BIND(C)
       double precision,intent(in)         :: c, s         
       integer,intent(in)                  :: incx, incy, n
       complex*16,intent(inout)            :: x( * ), y( * )
       !$omp  declare variant( zdrot:mkl_blas_zdrot_omp_offload_ilp64 ) match( construct={target variant dispatch}, device={arch(gen)} )
     end subroutine zdrot

     subroutine csrot ( n, x, incx,   &
          &y, incy, c, s ) BIND(C)
       real,intent(in)                     :: c, s         
       integer,intent(in)                  :: incx, incy, n
       complex,intent(inout)               :: x( * ), y( * )
       !$omp  declare variant( csrot:mkl_blas_csrot_omp_offload_ilp64 ) match( construct={target variant dispatch}, device={arch(gen)} )
     end subroutine csrot

     subroutine drotg( a, b, c, s ) BIND(C)
       double precision,intent(inout)      :: a( * ), b( * ), c( * ), s( * )
       !$omp  declare variant( drotg:mkl_blas_drotg_omp_offload_ilp64 ) match( construct={target variant dispatch}, device={arch(gen)} )
     end subroutine drotg

     subroutine srotg( a, b, c, s ) BIND(C)
       real,intent(inout)                  :: a( * ), b( * ), c( * ), s( * )
       !$omp  declare variant( srotg:mkl_blas_srotg_omp_offload_ilp64 ) match( construct={target variant dispatch}, device={arch(gen)} )
     end subroutine srotg

     subroutine zrotg( a, b, c, s ) BIND(C)
       complex*16,intent(inout)            :: a( * ), b( * ), s( * )
       double precision,intent(inout)      :: c( * )
       !$omp  declare variant( zrotg:mkl_blas_zrotg_omp_offload_ilp64 ) match( construct={target variant dispatch}, device={arch(gen)} )
     end subroutine zrotg

     subroutine crotg( a, b, c, s ) BIND(C)
       complex,intent(inout)               :: a( * ), b( * ), s( * )
       real,intent(inout)                  :: c( * )
       !$omp  declare variant( crotg:mkl_blas_crotg_omp_offload_ilp64 ) match( construct={target variant dispatch}, device={arch(gen)} )
     end subroutine crotg

     subroutine drotm ( n, x, incx,   &
          &y, incy, param ) BIND(C)
       integer,intent(in)                  :: incx, incy, n
       double precision,intent(in)         :: param( * )
       double precision,intent(inout)      :: x( * ), y( * )
       !$omp  declare variant( drotm:mkl_blas_drotm_omp_offload_ilp64 ) match( construct={target variant dispatch}, device={arch(gen)} )
     end subroutine drotm

     subroutine srotm ( n, x, incx,   &
          &y, incy, param ) BIND(C)
       integer,intent(in)                  :: incx, incy, n
       real,intent(in)                     :: param( * )
       real,intent(inout)                  :: x( * ), y( * )
       !$omp  declare variant( srotm:mkl_blas_srotm_omp_offload_ilp64 ) match( construct={target variant dispatch}, device={arch(gen)} )
     end subroutine srotm

     subroutine drotmg ( d1, d2, x1, x2, param ) BIND(C)
       double precision,intent(in)         :: x2
       double precision,intent(inout)      :: d1( * ), d2( * ), x1( * ), param( * )
       !$omp  declare variant( drotmg:mkl_blas_drotmg_omp_offload_ilp64 ) match( construct={target variant dispatch}, device={arch(gen)} )
     end subroutine drotmg

     subroutine srotmg ( d1, d2, x1, x2, param ) BIND(C)
       real,intent(in)                     :: x2
       real,intent(inout)                  :: d1( * ), d2( * ), x1( * ), param( * )
       !$omp  declare variant( srotmg:mkl_blas_srotmg_omp_offload_ilp64 ) match( construct={target variant dispatch}, device={arch(gen)} )
     end subroutine srotmg

     subroutine dscal ( n, alpha, x, incx ) BIND(C)
       double precision,intent(in)         :: alpha
       integer,intent(in)                  :: incx, n
       double precision,intent(inout)      :: x( * )
       !$omp  declare variant( dscal:mkl_blas_dscal_omp_offload_ilp64 ) match( construct={target variant dispatch}, device={arch(gen)} )
     end subroutine dscal

     subroutine sscal ( n, alpha, x, incx ) BIND(C)
       real,intent(in)                     :: alpha
       integer,intent(in)                  :: incx, n
       real,intent(inout)                  :: x( * )
       !$omp  declare variant( sscal:mkl_blas_sscal_omp_offload_ilp64 ) match( construct={target variant dispatch}, device={arch(gen)} )
     end subroutine sscal

     subroutine zscal ( n, alpha, x, incx ) BIND(C)
       complex*16,intent(in)               :: alpha
       integer,intent(in)                  :: incx, n
       complex*16,intent(inout)            :: x( * )
       !$omp  declare variant( zscal:mkl_blas_zscal_omp_offload_ilp64 ) match( construct={target variant dispatch}, device={arch(gen)} )
     end subroutine zscal

     subroutine cscal ( n, alpha, x, incx ) BIND(C)
       complex,intent(in)                  :: alpha
       integer,intent(in)                  :: incx, n
       complex,intent(inout)               :: x( * )
       !$omp  declare variant( cscal:mkl_blas_cscal_omp_offload_ilp64 ) match( construct={target variant dispatch}, device={arch(gen)} )
     end subroutine cscal

     subroutine zdscal ( n, alpha, x, incx ) BIND(C)
       double precision,intent(in)         :: alpha
       integer,intent(in)                  :: incx, n
       complex*16,intent(inout)            :: x( * )
       !$omp  declare variant( zdscal:mkl_blas_zdscal_omp_offload_ilp64 ) match( construct={target variant dispatch}, device={arch(gen)} )
     end subroutine zdscal

     subroutine csscal ( n, alpha, x, incx ) BIND(C)
       real,intent(in)                     :: alpha
       integer,intent(in)                  :: incx, n
       complex,intent(inout)               :: x( * )
       !$omp  declare variant( csscal:mkl_blas_csscal_omp_offload_ilp64 ) match( construct={target variant dispatch}, device={arch(gen)} )
     end subroutine csscal

     subroutine dswap ( n, x, incx,   &
          &y, incy ) BIND(C)
       integer,intent(in)                  :: incx, incy, n
       double precision,intent(inout)      :: x( * ), y( * )
       !$omp  declare variant( dswap:mkl_blas_dswap_omp_offload_ilp64 ) match( construct={target variant dispatch}, device={arch(gen)} )
     end subroutine dswap

     subroutine sswap ( n, x, incx,   &
          &y, incy ) BIND(C)
       integer,intent(in)                  :: incx, incy, n
       real,intent(inout)                  :: x( * ), y( * )
       !$omp  declare variant( sswap:mkl_blas_sswap_omp_offload_ilp64 ) match( construct={target variant dispatch}, device={arch(gen)} )
     end subroutine sswap

     subroutine zswap ( n, x, incx,   &
          &y, incy ) BIND(C)
       integer,intent(in)                  :: incx, incy, n
       complex*16,intent(inout)            :: x( * ), y( * )
       !$omp  declare variant( zswap:mkl_blas_zswap_omp_offload_ilp64 ) match( construct={target variant dispatch}, device={arch(gen)} )
     end subroutine zswap

     subroutine cswap ( n, x, incx,   &
          &y, incy ) BIND(C)
       integer,intent(in)                  :: incx, incy, n
       complex,intent(inout)               :: x( * ), y( * )
       !$omp  declare variant( cswap:mkl_blas_cswap_omp_offload_ilp64 ) match( construct={target variant dispatch}, device={arch(gen)} )
     end subroutine cswap

     function idamax ( n, x, incx ) BIND(C)
       integer                             :: idamax
       integer,intent(in)                  :: incx, n
       double precision,intent(in)         :: x( * )
       !$omp  declare variant( idamax:mkl_blas_idamax_omp_offload_ilp64 ) match( construct={target variant dispatch}, device={arch(gen)} )
     end function idamax

     function isamax ( n, x, incx ) BIND(C)
       integer                             :: isamax
       integer,intent(in)                  :: incx, n
       real,intent(in)                     :: x( * )
       !$omp  declare variant( isamax:mkl_blas_isamax_omp_offload_ilp64 ) match( construct={target variant dispatch}, device={arch(gen)} )
     end function isamax

     function izamax ( n, x, incx ) BIND(C)
       integer                             :: izamax
       integer,intent(in)                  :: incx, n
       complex*16,intent(in)               :: x( * )
       !$omp  declare variant( izamax:mkl_blas_izamax_omp_offload_ilp64 ) match( construct={target variant dispatch}, device={arch(gen)} )
     end function izamax

     function icamax ( n, x, incx ) BIND(C)
       integer                             :: icamax
       integer,intent(in)                  :: incx, n
       complex,intent(in)                  :: x( * )
       !$omp  declare variant( icamax:mkl_blas_icamax_omp_offload_ilp64 ) match( construct={target variant dispatch}, device={arch(gen)} )
     end function icamax

     function idamin ( n, x, incx ) BIND(C)
       integer                             :: idamin
       integer,intent(in)                  :: incx, n
       double precision,intent(in)         :: x( * )
       !$omp  declare variant( idamin:mkl_blas_idamin_omp_offload_ilp64 ) match( construct={target variant dispatch}, device={arch(gen)} )
     end function idamin

     function isamin ( n, x, incx ) BIND(C)
       integer                             :: isamin
       integer,intent(in)                  :: incx, n
       real,intent(in)                     :: x( * )
       !$omp  declare variant( isamin:mkl_blas_isamin_omp_offload_ilp64 ) match( construct={target variant dispatch}, device={arch(gen)} )
     end function isamin

     function izamin ( n, x, incx ) BIND(C)
       integer                             :: izamin
       integer,intent(in)                  :: incx, n
       complex*16,intent(in)               :: x( * )
       !$omp  declare variant( izamin:mkl_blas_izamin_omp_offload_ilp64 ) match( construct={target variant dispatch}, device={arch(gen)} )
     end function izamin

     function icamin ( n, x, incx ) BIND(C)
       integer                             :: icamin
       integer,intent(in)                  :: incx, n
       complex,intent(in)                  :: x( * )
       !$omp  declare variant( icamin:mkl_blas_icamin_omp_offload_ilp64 ) match( construct={target variant dispatch}, device={arch(gen)} )
     end function icamin

     ! BLAS batch

     subroutine dtrsm_batch_strided ( side, uplo, trans, diag, m, n, alpha, a, lda, stridea,        &
          &b, ldb, strideb, batch_size ) BIND(C)
       character*1,intent(in)             :: side, uplo, trans, diag
       integer,intent(in)                 :: m, n, lda, stridea, ldb, strideb, batch_size
       double precision,intent(in)        :: alpha
       double precision,intent(in)        :: a( lda, * )
       double precision,intent(inout)     :: b( ldb, * )
       !$omp  declare variant( dtrsm_batch_strided:mkl_blas_dtrsm_batch_strided_omp_offload_ilp64 ) match( construct={target variant dispatch}, device={arch(gen)} )
     end subroutine dtrsm_batch_strided

     subroutine strsm_batch_strided ( side, uplo, trans, diag, m, n, alpha, a, lda, stridea,        &
          &b, ldb, strideb, batch_size ) BIND(C)
       character*1,intent(in)              :: side, uplo, trans, diag
       integer,intent(in)                  :: m, n, lda, stridea, ldb, strideb, batch_size
       real,intent(in)                     :: alpha
       real,intent(in)                     :: a( lda, * )
       real,intent(in)                     :: b( ldb, * )
       !$omp  declare variant( strsm_batch_strided:mkl_blas_strsm_batch_strided_omp_offload_ilp64 ) match( construct={target variant dispatch}, device={arch(gen)} )
     end subroutine strsm_batch_strided

     subroutine ztrsm_batch_strided ( side, uplo, trans, diag, m, n, alpha, a, lda, stridea,        &
          &b, ldb, strideb, batch_size ) BIND(C)
       character*1,intent(in)             :: side, uplo, trans, diag
       integer,intent(in)                 :: m, n, lda, stridea, ldb, strideb, batch_size
       complex*16,intent(in)              :: alpha
       complex*16,intent(in)              :: a( lda, * )
       complex*16,intent(inout)           :: b( ldb, * )
       !$omp  declare variant( ztrsm_batch_strided:mkl_blas_ztrsm_batch_strided_omp_offload_ilp64 ) match( construct={target variant dispatch}, device={arch(gen)} )
     end subroutine ztrsm_batch_strided

     subroutine ctrsm_batch_strided ( side, uplo, trans, diag, m, n, alpha, a, lda, stridea,        &
          &b, ldb, strideb, batch_size ) BIND(C)
       character*1,intent(in)              :: side, uplo, trans, diag
       integer,intent(in)                  :: m, n, lda, stridea, ldb, strideb, batch_size
       complex,intent(in)                  :: alpha
       complex,intent(in)                  :: a( lda, * )
       complex,intent(in)                  :: b( ldb, * )
       !$omp  declare variant( ctrsm_batch_strided:mkl_blas_ctrsm_batch_strided_omp_offload_ilp64 ) match( construct={target variant dispatch}, device={arch(gen)} )
     end subroutine ctrsm_batch_strided

     subroutine dgemm_batch_strided ( transa, transb, m, n, k, alpha, a, lda, stridea,        &
          &b, ldb, strideb, beta, c, ldc, stridec, batch_size ) BIND(C)
       character*1,intent(in)             :: transa, transb
       integer,intent(in)                 :: m, n, k, lda, stridea, ldb, strideb, ldc, stridec, batch_size
       double precision,intent(in)        :: alpha, beta
       double precision,intent(in)        :: a( lda, * ), b( ldb, * )
       double precision,intent(inout)     :: c( ldc, * )
       !$omp  declare variant( dgemm_batch_strided:mkl_blas_dgemm_batch_strided_omp_offload_ilp64 ) match( construct={target variant dispatch}, device={arch(gen)} )
     end subroutine dgemm_batch_strided

     subroutine sgemm_batch_strided ( transa, transb, m, n, k, alpha, a, lda, stridea,        &
          &b, ldb, strideb, beta, c, ldc, stridec, batch_size ) BIND(C)
       character*1,intent(in)              :: transa, transb
       integer,intent(in)                  :: m, n, k, lda, stridea, ldb, strideb, ldc, stridec, batch_size
       real,intent(in)                     :: alpha, beta
       real,intent(in)                     :: a( lda, * ), b( ldb, * )
       real,intent(in)                     :: c( ldc, * )
       !$omp  declare variant( sgemm_batch_strided:mkl_blas_sgemm_batch_strided_omp_offload_ilp64 ) match( construct={target variant dispatch}, device={arch(gen)} )
     end subroutine sgemm_batch_strided

     subroutine zgemm_batch_strided ( transa, transb, m, n, k, alpha, a, lda, stridea,        &
          &b, ldb, strideb, beta, c, ldc, stridec, batch_size ) BIND(C)
       character*1,intent(in)             :: transa, transb
       integer,intent(in)                 :: m, n, k, lda, stridea, ldb, strideb, ldc, stridec, batch_size
       complex*16,intent(in)              :: alpha, beta
       complex*16,intent(in)              :: a( lda, * ), b( ldb, * )
       complex*16,intent(inout)           :: c( ldc, * )
       !$omp  declare variant( zgemm_batch_strided:mkl_blas_zgemm_batch_strided_omp_offload_ilp64 ) match( construct={target variant dispatch}, device={arch(gen)} )
     end subroutine zgemm_batch_strided

     subroutine cgemm_batch_strided ( transa, transb, m, n, k, alpha, a, lda, stridea,        &
          &b, ldb, strideb, beta, c, ldc, stridec, batch_size ) BIND(C)
       character*1,intent(in)              :: transa, transb
       integer,intent(in)                  :: m, n, k, lda, stridea, ldb, strideb, ldc, stridec, batch_size
       complex,intent(in)                  :: alpha, beta
       complex,intent(in)                  :: a( lda, * ), b( ldb, * )
       complex,intent(in)                  :: c( ldc, * )
       !$omp  declare variant( cgemm_batch_strided:mkl_blas_cgemm_batch_strided_omp_offload_ilp64 ) match( construct={target variant dispatch}, device={arch(gen)} )
     end subroutine cgemm_batch_strided

    subroutine daxpy_batch_strided ( n, alpha, x, incx, stridex,   &
          &y, incy, stridey, batch_size ) BIND(C)
       double precision,intent(in)         :: alpha
       integer,intent(in)                  :: incx, stridex, incy, stridey, batch_size, n
       double precision,intent(in)         :: x( * )
       double precision,intent(inout)      :: y( * )
       !$omp  declare variant( daxpy_batch_strided:mkl_blas_daxpy_batch_strided_omp_offload_ilp64 ) match( construct={target variant dispatch}, device={arch(gen)} )
     end subroutine daxpy_batch_strided

     subroutine saxpy_batch_strided ( n, alpha, x, incx, stridex,   &
          &y, incy, stridey, batch_size ) BIND(C)
       real,intent(in)                     :: alpha
       integer,intent(in)                  :: incx, stridex, incy, stridey, batch_size, n
       real,intent(in)                     :: x( * )
       real,intent(inout)                  :: y( * )
       !$omp  declare variant( saxpy_batch_strided:mkl_blas_saxpy_batch_strided_omp_offload_ilp64 ) match( construct={target variant dispatch}, device={arch(gen)} )
     end subroutine saxpy_batch_strided

     subroutine zaxpy_batch_strided ( n, alpha, x, incx, stridex,   &
          &y, incy, stridey, batch_size ) BIND(C)
       complex*16,intent(in)               :: alpha
       integer,intent(in)                  :: incx, stridex, incy, stridey, batch_size, n
       complex*16,intent(in)               :: x( * )
       complex*16,intent(inout)            :: y( * )
       !$omp  declare variant( zaxpy_batch_strided:mkl_blas_zaxpy_batch_strided_omp_offload_ilp64 ) match( construct={target variant dispatch}, device={arch(gen)} )
     end subroutine zaxpy_batch_strided

     subroutine caxpy_batch_strided ( n, alpha, x, incx, stridex,   &
          &y, incy, stridey, batch_size ) BIND(C)
       complex,intent(in)                  :: alpha
       integer,intent(in)                  :: incx, stridex, incy, stridey, batch_size, n
       complex,intent(in)                  :: x( * )
       complex,intent(inout)               :: y( * )
       !$omp  declare variant( caxpy_batch_strided:mkl_blas_caxpy_batch_strided_omp_offload_ilp64 ) match( construct={target variant dispatch}, device={arch(gen)} )
     end subroutine caxpy_batch_strided

     subroutine dgemv_batch_strided ( trans, m, n, alpha, a, lda, stridea, x, incx,   &
          &stridex, beta, y, incy, stridey, batch_size ) BIND(C)
       double precision,intent(in)         :: alpha, beta
       integer,intent(in)                  :: incx, incy, lda, m, n
       integer,intent(in)                  :: stridea, stridex, stridey, batch_size
       character*1,intent(in)              :: trans
       double precision,intent(in)         :: a( lda, * ), x( * )
       double precision,intent(inout)      :: y( * )
       !$omp  declare variant( dgemv_batch_strided:mkl_blas_dgemv_batch_strided_omp_offload_ilp64 ) match( construct={target variant dispatch}, device={arch(gen)} )
     end subroutine dgemv_batch_strided

     subroutine sgemv_batch_strided ( trans, m, n, alpha, a, lda, stridea, x, incx,   &
          &stridex, beta, y, incy, stridey, batch_size ) BIND(C)
       real,intent(in)                     :: alpha, beta
       integer,intent(in)                  :: incx, incy, lda, m, n
       integer,intent(in)                  :: stridea, stridex, stridey, batch_size
       character*1,intent(in)              :: trans
       real,intent(in)                     :: a( lda, * ), x( * )
       real,intent(inout)                  :: y( * )
       !$omp  declare variant( sgemv_batch_strided:mkl_blas_sgemv_batch_strided_omp_offload_ilp64 ) match( construct={target variant dispatch}, device={arch(gen)} )
     end subroutine sgemv_batch_strided

     subroutine zgemv_batch_strided ( trans, m, n, alpha, a, lda, stridea, x, incx,   &
          &stridex, beta, y, incy, stridey, batch_size ) BIND(C)
       complex*16,intent(in)               :: alpha, beta
       integer,intent(in)                  :: incx, incy, lda, m, n
       integer,intent(in)                  :: stridea, stridex, stridey, batch_size
       character*1,intent(in)              :: trans
       complex*16,intent(in)               :: a( lda, * ), x( * )
       complex*16,intent(inout)            :: y( * )
       !$omp  declare variant( zgemv_batch_strided:mkl_blas_zgemv_batch_strided_omp_offload_ilp64 ) match( construct={target variant dispatch}, device={arch(gen)} )
     end subroutine zgemv_batch_strided

     subroutine cgemv_batch_strided ( trans, m, n, alpha, a, lda, stridea, x, incx,   &
          &stridex, beta, y, incy, stridey, batch_size ) BIND(C)
       complex,intent(in)                  :: alpha, beta
       integer,intent(in)                  :: incx, incy, lda, m, n
       integer,intent(in)                  :: stridea, stridex, stridey, batch_size
       character*1,intent(in)              :: trans
       complex,intent(in)                  :: a( lda, * ), x( * )
       complex,intent(inout)               :: y( * )
       !$omp  declare variant( cgemv_batch_strided:mkl_blas_cgemv_batch_strided_omp_offload_ilp64 ) match( construct={target variant dispatch}, device={arch(gen)} )
     end subroutine cgemv_batch_strided

     subroutine ddgmm_batch_strided ( side, m, n, a, lda, stridea, x, incx,   &
          &stridex, c, ldc, stridec, batch_size ) BIND(C)
       integer,intent(in)                  :: incx, ldc, lda, m, n
       integer,intent(in)                  :: stridea, stridex, stridec, batch_size
       character*1,intent(in)              :: side
       double precision,intent(in)         :: a( lda, * ), x( * )
       double precision,intent(inout)      :: c( ldc, * )
       !$omp  declare variant( ddgmm_batch_strided:mkl_blas_ddgmm_batch_strided_omp_offload_ilp64 ) match( construct={target variant dispatch}, device={arch(gen)} )
     end subroutine ddgmm_batch_strided

     subroutine sdgmm_batch_strided ( side, m, n, a, lda, stridea, x, incx,   &
          &stridex, c, ldc, stridec, batch_size ) BIND(C)
       integer,intent(in)                  :: incx, ldc, lda, m, n
       integer,intent(in)                  :: stridea, stridex, stridec, batch_size
       character*1,intent(in)              :: side
       real,intent(in)                     :: a( lda, * ), x( * )
       real,intent(inout)                  :: c( ldc, * )
       !$omp  declare variant( sdgmm_batch_strided:mkl_blas_sdgmm_batch_strided_omp_offload_ilp64 ) match( construct={target variant dispatch}, device={arch(gen)} )
     end subroutine sdgmm_batch_strided

     subroutine zdgmm_batch_strided ( side, m, n, a, lda, stridea, x, incx,   &
          &stridex, c, ldc, stridec, batch_size ) BIND(C)
       integer,intent(in)                  :: incx, ldc, lda, m, n
       integer,intent(in)                  :: stridea, stridex, stridec, batch_size
       character*1,intent(in)              :: side
       complex*16,intent(in)               :: a( lda, * ), x( * )
       complex*16,intent(inout)            :: c( ldc, * )
       !$omp  declare variant( zdgmm_batch_strided:mkl_blas_zdgmm_batch_strided_omp_offload_ilp64 ) match( construct={target variant dispatch}, device={arch(gen)} )
     end subroutine zdgmm_batch_strided

     subroutine cdgmm_batch_strided ( side, m, n, a, lda, stridea, x, incx,   &
          &stridex, c, ldc, stridec, batch_size ) BIND(C)
       integer,intent(in)                  :: incx, ldc, lda, m, n
       integer,intent(in)                  :: stridea, stridex, stridec, batch_size
       character*1,intent(in)              :: side
       complex,intent(in)                  :: a( lda, * ), x( * )
       complex,intent(inout)               :: c( ldc, * )
       !$omp  declare variant( cdgmm_batch_strided:mkl_blas_cdgmm_batch_strided_omp_offload_ilp64 ) match( construct={target variant dispatch}, device={arch(gen)} )
     end subroutine cdgmm_batch_strided

  end interface

end module onemkl_blas_omp_offload_ilp64
