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
  
  interface

     ! BLAS Level3

     subroutine mkl_blas_dgemm_omp_offload_lp64 ( transa, transb, m, n, k, alpha, a, lda,        &
          &b, ldb, beta, c, ldc ) BIND(C)
       character*1,intent(in)             :: transa, transb
       integer,intent(in)                 :: m, n, k, lda, ldb, ldc
       double precision,intent(in)        :: alpha, beta
       double precision,intent(in)        :: a( lda, * ), b( ldb, * )
       double precision,intent(inout)     :: c( ldc, * )
     end subroutine mkl_blas_dgemm_omp_offload_lp64

     subroutine mkl_blas_sgemm_omp_offload_lp64 ( transa, transb, m, n, k, alpha, a, lda,        &
          &b, ldb, beta, c, ldc ) BIND(C)
       character*1,intent(in)              :: transa, transb
       integer,intent(in)                  :: m, n, k, lda, ldb, ldc
       real,intent(in)                     :: alpha, beta
       real,intent(in)                     :: a( lda, * ), b( ldb, * )
       real,intent(in)                     :: c( ldc, * )
     end subroutine mkl_blas_sgemm_omp_offload_lp64

     subroutine mkl_blas_zgemm_omp_offload_lp64 ( transa, transb, m, n, k, alpha, a, lda,        &
          &b, ldb, beta, c, ldc ) BIND(C)
       character*1,intent(in)             :: transa, transb
       integer,intent(in)                 :: m, n, k, lda, ldb, ldc
       complex*16,intent(in)              :: alpha, beta
       complex*16,intent(in)              :: a( lda, * ), b( ldb, * )
       complex*16,intent(inout)           :: c( ldc, * )
     end subroutine mkl_blas_zgemm_omp_offload_lp64

     subroutine mkl_blas_cgemm_omp_offload_lp64 ( transa, transb, m, n, k, alpha, a, lda,        &
          &b, ldb, beta, c, ldc ) BIND(C)
       character*1,intent(in)              :: transa, transb
       integer,intent(in)                  :: m, n, k, lda, ldb, ldc
       complex,intent(in)                  :: alpha, beta
       complex,intent(in)                  :: a( lda, * ), b( ldb, * )
       complex,intent(in)                  :: c( ldc, * )
     end subroutine mkl_blas_cgemm_omp_offload_lp64

     subroutine mkl_blas_dsymm_omp_offload_lp64 ( side, uplo, m, n, alpha, a, lda,        &
          &b, ldb, beta, c, ldc ) BIND(C)
       character*1,intent(in)             :: side, uplo
       integer,intent(in)                 :: m, n, lda, ldb, ldc
       double precision,intent(in)        :: alpha, beta
       double precision,intent(in)        :: a( lda, * ), b( ldb, * )
       double precision,intent(inout)     :: c( ldc, * )
     end subroutine mkl_blas_dsymm_omp_offload_lp64

     subroutine mkl_blas_ssymm_omp_offload_lp64 ( side, uplo, m, n, alpha, a, lda,        &
          &b, ldb, beta, c, ldc ) BIND(C)
       character*1,intent(in)              :: side, uplo
       integer,intent(in)                  :: m, n, lda, ldb, ldc
       real,intent(in)                     :: alpha, beta
       real,intent(in)                     :: a( lda, * ), b( ldb, * )
       real,intent(in)                     :: c( ldc, * )
     end subroutine mkl_blas_ssymm_omp_offload_lp64

     subroutine mkl_blas_zsymm_omp_offload_lp64 ( side, uplo, m, n, alpha, a, lda,        &
          &b, ldb, beta, c, ldc ) BIND(C)
       character*1,intent(in)             :: side, uplo
       integer,intent(in)                 :: m, n, lda, ldb, ldc
       complex*16,intent(in)              :: alpha, beta
       complex*16,intent(in)              :: a( lda, * ), b( ldb, * )
       complex*16,intent(inout)           :: c( ldc, * )
     end subroutine mkl_blas_zsymm_omp_offload_lp64

     subroutine mkl_blas_csymm_omp_offload_lp64 ( side, uplo, m, n, alpha, a, lda,        &
          &b, ldb, beta, c, ldc ) BIND(C)
       character*1,intent(in)              :: side, uplo
       integer,intent(in)                  :: m, n, lda, ldb, ldc
       complex,intent(in)                  :: alpha, beta
       complex,intent(in)                  :: a( lda, * ), b( ldb, * )
       complex,intent(in)                  :: c( ldc, * )
     end subroutine mkl_blas_csymm_omp_offload_lp64

     subroutine mkl_blas_zhemm_omp_offload_lp64 ( side, uplo, m, n, alpha, a, lda,        &
          &b, ldb, beta, c, ldc ) BIND(C)
       character*1,intent(in)             :: side, uplo
       integer,intent(in)                 :: m, n, lda, ldb, ldc
       complex*16,intent(in)              :: alpha, beta
       complex*16,intent(in)              :: a( lda, * ), b( ldb, * )
       complex*16,intent(inout)           :: c( ldc, * )
     end subroutine mkl_blas_zhemm_omp_offload_lp64

     subroutine mkl_blas_chemm_omp_offload_lp64 ( side, uplo, m, n, alpha, a, lda,        &
          &b, ldb, beta, c, ldc ) BIND(C)
       character*1,intent(in)              :: side, uplo
       integer,intent(in)                  :: m, n, lda, ldb, ldc
       complex,intent(in)                  :: alpha, beta
       complex,intent(in)                  :: a( lda, * ), b( ldb, * )
       complex,intent(in)                  :: c( ldc, * )
     end subroutine mkl_blas_chemm_omp_offload_lp64

     subroutine mkl_blas_dsyrk_omp_offload_lp64 ( uplo, trans, n, k, alpha, a, lda,        &
          &beta, c, ldc ) BIND(C)
       character*1,intent(in)             :: uplo, trans
       integer,intent(in)                 :: n, k, lda, ldc
       double precision,intent(in)        :: alpha, beta
       double precision,intent(in)        :: a( lda, * )
       double precision,intent(inout)     :: c( ldc, * )
     end subroutine mkl_blas_dsyrk_omp_offload_lp64

     subroutine mkl_blas_ssyrk_omp_offload_lp64 ( uplo, trans, n, k, alpha, a, lda,        &
          &beta, c, ldc ) BIND(C)
       character*1,intent(in)              :: uplo, trans
       integer,intent(in)                  :: n, k, lda, ldc
       real,intent(in)                     :: alpha, beta
       real,intent(in)                     :: a( lda, * )
       real,intent(in)                     :: c( ldc, * )
     end subroutine mkl_blas_ssyrk_omp_offload_lp64

     subroutine mkl_blas_zsyrk_omp_offload_lp64 ( uplo, trans, n, k, alpha, a, lda,        &
          &beta, c, ldc ) BIND(C)
       character*1,intent(in)             :: uplo, trans
       integer,intent(in)                 :: n, k, lda, ldc
       complex*16,intent(in)              :: alpha, beta
       complex*16,intent(in)              :: a( lda, * )
       complex*16,intent(inout)           :: c( ldc, * )
     end subroutine mkl_blas_zsyrk_omp_offload_lp64

     subroutine mkl_blas_csyrk_omp_offload_lp64 ( uplo, trans, n, k, alpha, a, lda,        &
          &beta, c, ldc ) BIND(C)
       character*1,intent(in)              :: uplo, trans
       integer,intent(in)                  :: n, k, lda, ldc
       complex,intent(in)                  :: alpha, beta
       complex,intent(in)                  :: a( lda, * )
       complex,intent(in)                  :: c( ldc, * )
     end subroutine mkl_blas_csyrk_omp_offload_lp64

     subroutine mkl_blas_zherk_omp_offload_lp64 ( uplo, trans, n, k, alpha, a, lda,        &
          &beta, c, ldc ) BIND(C)
       character*1,intent(in)             :: uplo, trans
       integer,intent(in)                 :: n, k, lda, ldc
       complex*16,intent(in)              :: alpha, beta
       complex*16,intent(in)              :: a( lda, * )
       complex*16,intent(inout)           :: c( ldc, * )
     end subroutine mkl_blas_zherk_omp_offload_lp64

     subroutine mkl_blas_cherk_omp_offload_lp64 ( uplo, trans, n, k, alpha, a, lda,        &
          &beta, c, ldc ) BIND(C)
       character*1,intent(in)              :: uplo, trans
       integer,intent(in)                  :: n, k, lda, ldc
       complex,intent(in)                  :: alpha, beta
       complex,intent(in)                  :: a( lda, * )
       complex,intent(in)                  :: c( ldc, * )
     end subroutine mkl_blas_cherk_omp_offload_lp64

     subroutine mkl_blas_dsyr2k_omp_offload_lp64 ( uplo, trans, n, k, alpha, a, lda,        &
          &b, ldb, beta, c, ldc ) BIND(C)
       character*1,intent(in)             :: uplo, trans
       integer,intent(in)                 :: n, k, lda, ldb, ldc
       double precision,intent(in)        :: alpha, beta
       double precision,intent(in)        :: a( lda, * ), b( ldb, * )
       double precision,intent(inout)     :: c( ldc, * )
     end subroutine mkl_blas_dsyr2k_omp_offload_lp64

     subroutine mkl_blas_ssyr2k_omp_offload_lp64 ( uplo, trans, n, k, alpha, a, lda,        &
          &b, ldb, beta, c, ldc ) BIND(C)
       character*1,intent(in)              :: uplo, trans
       integer,intent(in)                  :: n, k, lda, ldb, ldc
       real,intent(in)                     :: alpha, beta
       real,intent(in)                     :: a( lda, * ), b( ldb, * )
       real,intent(in)                     :: c( ldc, * )
     end subroutine mkl_blas_ssyr2k_omp_offload_lp64

     subroutine mkl_blas_zsyr2k_omp_offload_lp64 ( uplo, trans, n, k, alpha, a, lda,        &
          &b, ldb, beta, c, ldc ) BIND(C)
       character*1,intent(in)             :: uplo, trans
       integer,intent(in)                 :: n, k, lda, ldb, ldc
       complex*16,intent(in)              :: alpha, beta
       complex*16,intent(in)              :: a( lda, * ), b( ldb, * )
       complex*16,intent(inout)           :: c( ldc, * )
     end subroutine mkl_blas_zsyr2k_omp_offload_lp64

     subroutine mkl_blas_csyr2k_omp_offload_lp64 ( uplo, trans, n, k, alpha, a, lda,        &
          &b, ldb, beta, c, ldc ) BIND(C)
       character*1,intent(in)              :: uplo, trans
       integer,intent(in)                  :: n, k, lda, ldb, ldc
       complex,intent(in)                  :: alpha, beta
       complex,intent(in)                  :: a( lda, * ), b( ldb, * )
       complex,intent(in)                  :: c( ldc, * )
     end subroutine mkl_blas_csyr2k_omp_offload_lp64

     subroutine mkl_blas_zher2k_omp_offload_lp64 ( uplo, trans, n, k, alpha, a, lda,        &
          &b, ldb, beta, c, ldc ) BIND(C)
       character*1,intent(in)             :: uplo, trans
       integer,intent(in)                 :: n, k, lda, ldb, ldc
       complex*16,intent(in)              :: alpha, beta
       complex*16,intent(in)              :: a( lda, * ), b( ldb, * )
       complex*16,intent(inout)           :: c( ldc, * )
     end subroutine mkl_blas_zher2k_omp_offload_lp64

     subroutine mkl_blas_cher2k_omp_offload_lp64 ( uplo, trans, n, k, alpha, a, lda,        &
          &b, ldb, beta, c, ldc ) BIND(C)
       character*1,intent(in)              :: uplo, trans
       integer,intent(in)                  :: n, k, lda, ldb, ldc
       complex,intent(in)                  :: alpha, beta
       complex,intent(in)                  :: a( lda, * ), b( ldb, * )
       complex,intent(in)                  :: c( ldc, * )
     end subroutine mkl_blas_cher2k_omp_offload_lp64

     subroutine mkl_blas_dtrmm_omp_offload_lp64 ( side, uplo, trans, diag, m, n, alpha, a, lda,        &
          &b, ldb ) BIND(C)
       character*1,intent(in)             :: side, uplo, trans, diag
       integer,intent(in)                 :: m, n, lda, ldb
       double precision,intent(in)        :: alpha
       double precision,intent(in)        :: a( lda, * )
       double precision,intent(inout)     :: b( ldb, * )
     end subroutine mkl_blas_dtrmm_omp_offload_lp64

     subroutine mkl_blas_strmm_omp_offload_lp64 ( side, uplo, trans, diag, m, n, alpha, a, lda,        &
          &b, ldb ) BIND(C)
       character*1,intent(in)              :: side, uplo, trans, diag
       integer,intent(in)                  :: m, n, lda, ldb
       real,intent(in)                     :: alpha
       real,intent(in)                     :: a( lda, * )
       real,intent(in)                     :: b( ldb, * )
     end subroutine mkl_blas_strmm_omp_offload_lp64

     subroutine mkl_blas_ztrmm_omp_offload_lp64 ( side, uplo, trans, diag, m, n, alpha, a, lda,        &
          &b, ldb ) BIND(C)
       character*1,intent(in)             :: side, uplo, trans, diag
       integer,intent(in)                 :: m, n, lda, ldb
       complex*16,intent(in)              :: alpha
       complex*16,intent(in)              :: a( lda, * )
       complex*16,intent(inout)           :: b( ldb, * )
     end subroutine mkl_blas_ztrmm_omp_offload_lp64

     subroutine mkl_blas_ctrmm_omp_offload_lp64 ( side, uplo, trans, diag, m, n, alpha, a, lda,        &
          &b, ldb ) BIND(C)
       character*1,intent(in)              :: side, uplo, trans, diag
       integer,intent(in)                  :: m, n, lda, ldb
       complex,intent(in)                  :: alpha
       complex,intent(in)                  :: a( lda, * )
       complex,intent(in)                  :: b( ldb, * )
     end subroutine mkl_blas_ctrmm_omp_offload_lp64

     subroutine mkl_blas_dtrsm_omp_offload_lp64 ( side, uplo, trans, diag, m, n, alpha, a, lda,        &
          &b, ldb ) BIND(C)
       character*1,intent(in)             :: side, uplo, trans, diag
       integer,intent(in)                 :: m, n, lda, ldb
       double precision,intent(in)        :: alpha
       double precision,intent(in)        :: a( lda, * )
       double precision,intent(inout)     :: b( ldb, * )
     end subroutine mkl_blas_dtrsm_omp_offload_lp64

     subroutine mkl_blas_strsm_omp_offload_lp64 ( side, uplo, trans, diag, m, n, alpha, a, lda,        &
          &b, ldb ) BIND(C)
       character*1,intent(in)              :: side, uplo, trans, diag
       integer,intent(in)                  :: m, n, lda, ldb
       real,intent(in)                     :: alpha
       real,intent(in)                     :: a( lda, * )
       real,intent(in)                     :: b( ldb, * )
     end subroutine mkl_blas_strsm_omp_offload_lp64

     subroutine mkl_blas_ztrsm_omp_offload_lp64 ( side, uplo, trans, diag, m, n, alpha, a, lda,        &
          &b, ldb ) BIND(C)
       character*1,intent(in)             :: side, uplo, trans, diag
       integer,intent(in)                 :: m, n, lda, ldb
       complex*16,intent(in)              :: alpha
       complex*16,intent(in)              :: a( lda, * )
       complex*16,intent(inout)           :: b( ldb, * )
     end subroutine mkl_blas_ztrsm_omp_offload_lp64

     subroutine mkl_blas_ctrsm_omp_offload_lp64 ( side, uplo, trans, diag, m, n, alpha, a, lda,        &
          &b, ldb ) BIND(C)
       character*1,intent(in)              :: side, uplo, trans, diag
       integer,intent(in)                  :: m, n, lda, ldb
       complex,intent(in)                  :: alpha
       complex,intent(in)                  :: a( lda, * )
       complex,intent(in)                  :: b( ldb, * )
     end subroutine mkl_blas_ctrsm_omp_offload_lp64

     subroutine mkl_blas_dgemmt_omp_offload_lp64 ( uplo, transa, transb, m, n, alpha, a, lda,        &
          &b, ldb, beta, c, ldc ) BIND(C)
       character*1,intent(in)             :: uplo, transa, transb
       integer,intent(in)                 :: m, n, lda, ldb, ldc
       double precision,intent(in)        :: alpha, beta
       double precision,intent(in)        :: a( lda, * ), b( ldb, * )
       double precision,intent(inout)     :: c( ldc, * )
     end subroutine mkl_blas_dgemmt_omp_offload_lp64

     subroutine mkl_blas_sgemmt_omp_offload_lp64 ( uplo, transa, transb, m, n, alpha, a, lda,        &
          &b, ldb, beta, c, ldc ) BIND(C)
       character*1,intent(in)              :: uplo, transa, transb
       integer,intent(in)                  :: m, n, lda, ldb, ldc
       real,intent(in)                     :: alpha, beta
       real,intent(in)                     :: a( lda, * ), b( ldb, * )
       real,intent(in)                     :: c( ldc, * )
     end subroutine mkl_blas_sgemmt_omp_offload_lp64

     subroutine mkl_blas_zgemmt_omp_offload_lp64 ( uplo, transa, transb, m, n, alpha, a, lda,        &
          &b, ldb, beta, c, ldc ) BIND(C)
       character*1,intent(in)             :: uplo, transa, transb
       integer,intent(in)                 :: m, n, lda, ldb, ldc
       complex*16,intent(in)              :: alpha, beta
       complex*16,intent(in)              :: a( lda, * ), b( ldb, * )
       complex*16,intent(inout)           :: c( ldc, * )
     end subroutine mkl_blas_zgemmt_omp_offload_lp64

     subroutine mkl_blas_cgemmt_omp_offload_lp64 ( uplo, transa, transb, m, n, alpha, a, lda,        &
          &b, ldb, beta, c, ldc ) BIND(C)
       character*1,intent(in)              :: uplo, transa, transb
       integer,intent(in)                  :: m, n, lda, ldb, ldc
       complex,intent(in)                  :: alpha, beta
       complex,intent(in)                  :: a( lda, * ), b( ldb, * )
       complex,intent(in)                  :: c( ldc, * )
     end subroutine mkl_blas_cgemmt_omp_offload_lp64

     ! BLAS Level2

     subroutine mkl_blas_dgbmv_omp_offload_lp64 ( trans, m, n, kl, ku, alpha, a, lda, x, incx,   &
          &beta, y, incy ) BIND(C)
       double precision,intent(in)         :: alpha, beta
       integer,intent(in)                  :: incx, incy, kl, ku, lda, m, n
       character*1,intent(in)              :: trans
       double precision,intent(in)         :: a( lda, * ), x( * )
       double precision,intent(inout)      :: y( * )
     end subroutine mkl_blas_dgbmv_omp_offload_lp64

     subroutine mkl_blas_sgbmv_omp_offload_lp64 ( trans, m, n, kl, ku, alpha, a, lda, x, incx,   &
          &beta, y, incy ) BIND(C)
       real,intent(in)                     :: alpha, beta
       integer,intent(in)                  :: incx, incy, kl, ku, lda, m, n
       character*1,intent(in)              :: trans
       real,intent(in)                     :: a( lda, * ), x( * )
       real,intent(inout)                  :: y( * )
     end subroutine mkl_blas_sgbmv_omp_offload_lp64

     subroutine mkl_blas_zgbmv_omp_offload_lp64 ( trans, m, n, kl, ku, alpha, a, lda, x, incx,   &
          &beta, y, incy ) BIND(C)
       complex*16,intent(in)               :: alpha, beta
       integer,intent(in)                  :: incx, incy, kl, ku, lda, m, n
       character*1,intent(in)              :: trans
       complex*16,intent(in)               :: a( lda, * ), x( * )
       complex*16,intent(inout)            :: y( * )
     end subroutine mkl_blas_zgbmv_omp_offload_lp64

     subroutine mkl_blas_cgbmv_omp_offload_lp64 ( trans, m, n, kl, ku, alpha, a, lda, x, incx,   &
          &beta, y, incy ) BIND(C)
       complex,intent(in)                  :: alpha, beta
       integer,intent(in)                  :: incx, incy, kl, ku, lda, m, n
       character*1,intent(in)              :: trans
       complex,intent(in)                  :: a( lda, * ), x( * )
       complex,intent(inout)               :: y( * )
     end subroutine mkl_blas_cgbmv_omp_offload_lp64

     subroutine mkl_blas_dgemv_omp_offload_lp64 ( trans, m, n, alpha, a, lda, x, incx,   &
          &beta, y, incy ) BIND(C)
       double precision,intent(in)         :: alpha, beta
       integer,intent(in)                  :: incx, incy, lda, m, n
       character*1,intent(in)              :: trans
       double precision,intent(in)         :: a( lda, * ), x( * )
       double precision,intent(inout)      :: y( * )
     end subroutine mkl_blas_dgemv_omp_offload_lp64

     subroutine mkl_blas_sgemv_omp_offload_lp64 ( trans, m, n, alpha, a, lda, x, incx,   &
          &beta, y, incy ) BIND(C)
       real,intent(in)                     :: alpha, beta
       integer,intent(in)                  :: incx, incy, lda, m, n
       character*1,intent(in)              :: trans
       real,intent(in)                     :: a( lda, * ), x( * )
       real,intent(inout)                  :: y( * )
     end subroutine mkl_blas_sgemv_omp_offload_lp64

     subroutine mkl_blas_zgemv_omp_offload_lp64 ( trans, m, n, alpha, a, lda, x, incx,   &
          &beta, y, incy ) BIND(C)
       complex*16,intent(in)               :: alpha, beta
       integer,intent(in)                  :: incx, incy, lda, m, n
       character*1,intent(in)              :: trans
       complex*16,intent(in)               :: a( lda, * ), x( * )
       complex*16,intent(inout)            :: y( * )
     end subroutine mkl_blas_zgemv_omp_offload_lp64

     subroutine mkl_blas_cgemv_omp_offload_lp64 ( trans, m, n, alpha, a, lda, x, incx,   &
          &beta, y, incy ) BIND(C)
       complex,intent(in)                  :: alpha, beta
       integer,intent(in)                  :: incx, incy, lda, m, n
       character*1,intent(in)              :: trans
       complex,intent(in)                  :: a( lda, * ), x( * )
       complex,intent(inout)               :: y( * )
     end subroutine mkl_blas_cgemv_omp_offload_lp64

     subroutine mkl_blas_dsbmv_omp_offload_lp64 ( uplo, n, k, alpha, a, lda, x, incx,   &
          &beta, y, incy ) BIND(C)
       double precision,intent(in)         :: alpha, beta
       integer,intent(in)                  :: incx, incy, k, lda, n
       character*1,intent(in)              :: uplo
       double precision,intent(in)         :: a( lda, * ), x( * )
       double precision,intent(inout)      :: y( * )
     end subroutine mkl_blas_dsbmv_omp_offload_lp64

     subroutine mkl_blas_ssbmv_omp_offload_lp64 ( uplo, n, k, alpha, a, lda, x, incx,   &
          &beta, y, incy ) BIND(C)
       real,intent(in)                     :: alpha, beta
       integer,intent(in)                  :: incx, incy, k, lda, n
       character*1,intent(in)              :: uplo
       real,intent(in)                     :: a( lda, * ), x( * )
       real,intent(inout)                  :: y( * )
     end subroutine mkl_blas_ssbmv_omp_offload_lp64

     subroutine mkl_blas_zhbmv_omp_offload_lp64 ( uplo, n, k, alpha, a, lda, x, incx,   &
          &beta, y, incy ) BIND(C)
       complex*16,intent(in)               :: alpha, beta
       integer,intent(in)                  :: incx, incy, k, lda, n
       character*1,intent(in)              :: uplo
       complex*16,intent(in)               :: a( lda, * ), x( * )
       complex*16,intent(inout)            :: y( * )
     end subroutine mkl_blas_zhbmv_omp_offload_lp64

     subroutine mkl_blas_chbmv_omp_offload_lp64 ( uplo, n, k, alpha, a, lda, x, incx,   &
          &beta, y, incy ) BIND(C)
       complex,intent(in)                  :: alpha, beta
       integer,intent(in)                  :: incx, incy, k, lda, n
       character*1,intent(in)              :: uplo
       complex,intent(in)                  :: a( lda, * ), x( * )
       complex,intent(inout)               :: y( * )
     end subroutine mkl_blas_chbmv_omp_offload_lp64

     subroutine mkl_blas_dsymv_omp_offload_lp64 ( uplo, n, alpha, a, lda, x, incx,   &
          &beta, y, incy ) BIND(C)
       double precision,intent(in)         :: alpha, beta
       integer,intent(in)                  :: incx, incy, lda, n
       character*1,intent(in)              :: uplo
       double precision,intent(in)         :: a( lda, * ), x( * )
       double precision,intent(inout)      :: y( * )
     end subroutine mkl_blas_dsymv_omp_offload_lp64

     subroutine mkl_blas_ssymv_omp_offload_lp64 ( uplo, n, alpha, a, lda, x, incx,   &
          &beta, y, incy ) BIND(C)
       real,intent(in)                     :: alpha, beta
       integer,intent(in)                  :: incx, incy, lda, n
       character*1,intent(in)              :: uplo
       real,intent(in)                     :: a( lda, * ), x( * )
       real,intent(inout)                  :: y( * )
     end subroutine mkl_blas_ssymv_omp_offload_lp64

     subroutine mkl_blas_zhemv_omp_offload_lp64 ( uplo, n, alpha, a, lda, x, incx,   &
          &beta, y, incy ) BIND(C)
       complex*16,intent(in)               :: alpha, beta
       integer,intent(in)                  :: incx, incy, lda, n
       character*1,intent(in)              :: uplo
       complex*16,intent(in)               :: a( lda, * ), x( * )
       complex*16,intent(inout)            :: y( * )
     end subroutine mkl_blas_zhemv_omp_offload_lp64

     subroutine mkl_blas_chemv_omp_offload_lp64 ( uplo, n, alpha, a, lda, x, incx,   &
          &beta, y, incy ) BIND(C)
       complex,intent(in)                  :: alpha, beta
       integer,intent(in)                  :: incx, incy, lda, n
       character*1,intent(in)              :: uplo
       complex,intent(in)                  :: a( lda, * ), x( * )
       complex,intent(inout)               :: y( * )
     end subroutine mkl_blas_chemv_omp_offload_lp64

     subroutine mkl_blas_dspmv_omp_offload_lp64 ( uplo, n, alpha, a, x, incx,   &
          &beta, y, incy ) BIND(C)
       double precision,intent(in)         :: alpha, beta
       integer,intent(in)                  :: incx, incy, n
       character*1,intent(in)              :: uplo
       double precision,intent(in)         :: a( * ), x( * )
       double precision,intent(inout)      :: y( * )
     end subroutine mkl_blas_dspmv_omp_offload_lp64

     subroutine mkl_blas_sspmv_omp_offload_lp64 ( uplo, n, alpha, a, x, incx,   &
          &beta, y, incy ) BIND(C)
       real,intent(in)                     :: alpha, beta
       integer,intent(in)                  :: incx, incy, n
       character*1,intent(in)              :: uplo
       real,intent(in)                     :: a( * ), x( * )
       real,intent(inout)                  :: y( * )
     end subroutine mkl_blas_sspmv_omp_offload_lp64

     subroutine mkl_blas_zhpmv_omp_offload_lp64 ( uplo, n, alpha, a, x, incx,   &
          &beta, y, incy ) BIND(C)
       complex*16,intent(in)               :: alpha, beta
       integer,intent(in)                  :: incx, incy, n
       character*1,intent(in)              :: uplo
       complex*16,intent(in)               :: a( * ), x( * )
       complex*16,intent(inout)            :: y( * )
     end subroutine mkl_blas_zhpmv_omp_offload_lp64

     subroutine mkl_blas_chpmv_omp_offload_lp64 ( uplo, n, alpha, a, x, incx,   &
          &beta, y, incy ) BIND(C)
       complex,intent(in)                  :: alpha, beta
       integer,intent(in)                  :: incx, incy, n
       character*1,intent(in)              :: uplo
       complex,intent(in)                  :: a( * ), x( * )
       complex,intent(inout)               :: y( * )
     end subroutine mkl_blas_chpmv_omp_offload_lp64

     subroutine mkl_blas_dger_omp_offload_lp64 ( m, n, alpha, x, incx,   &
          &y, incy, a, lda ) BIND(C)
       double precision,intent(in)         :: alpha
       integer,intent(in)                  :: incx, incy, lda, m, n
       double precision,intent(in)         :: y( * ), x( * )
       double precision,intent(inout)      :: a( lda, * )
     end subroutine mkl_blas_dger_omp_offload_lp64

     subroutine mkl_blas_sger_omp_offload_lp64 ( m, n, alpha, x, incx,   &
          &y, incy, a, lda ) BIND(C)
       real,intent(in)                     :: alpha
       integer,intent(in)                  :: incx, incy, lda, m, n
       real,intent(in)                     :: y( * ), x( * )
       real,intent(inout)                  :: a( lda, * )
     end subroutine mkl_blas_sger_omp_offload_lp64

     subroutine mkl_blas_zgerc_omp_offload_lp64 ( m, n, alpha, x, incx,   &
          &y, incy, a, lda ) BIND(C)
       complex*16,intent(in)               :: alpha
       integer,intent(in)                  :: incx, incy, lda, m, n
       complex*16,intent(in)               :: y( * ), x( * )
       complex*16,intent(inout)            :: a( lda, * )
     end subroutine mkl_blas_zgerc_omp_offload_lp64

     subroutine mkl_blas_cgerc_omp_offload_lp64 ( m, n, alpha, x, incx,   &
          &y, incy, a, lda ) BIND(C)
       complex,intent(in)                  :: alpha
       integer,intent(in)                  :: incx, incy, lda, m, n
       complex,intent(in)                  :: y( * ), x( * )
       complex,intent(inout)               :: a( lda, * )
     end subroutine mkl_blas_cgerc_omp_offload_lp64

     subroutine mkl_blas_zgeru_omp_offload_lp64 ( m, n, alpha, x, incx,   &
          &y, incy, a, lda ) BIND(C)
       complex*16,intent(in)               :: alpha
       integer,intent(in)                  :: incx, incy, lda, m, n
       complex*16,intent(in)               :: y( * ), x( * )
       complex*16,intent(inout)            :: a( lda, * )
     end subroutine mkl_blas_zgeru_omp_offload_lp64

     subroutine mkl_blas_cgeru_omp_offload_lp64 ( m, n, alpha, x, incx,   &
          &y, incy, a, lda ) BIND(C)
       complex,intent(in)                  :: alpha
       integer,intent(in)                  :: incx, incy, lda, m, n
       complex,intent(in)                  :: y( * ), x( * )
       complex,intent(inout)               :: a( lda, * )
     end subroutine mkl_blas_cgeru_omp_offload_lp64

     subroutine mkl_blas_dsyr_omp_offload_lp64 ( uplo, n, alpha, x, incx,   &
          &a, lda ) BIND(C)
       character*1,intent(in)              :: uplo
       double precision,intent(in)         :: alpha
       integer,intent(in)                  :: incx, lda, n
       double precision,intent(in)         :: x( * )
       double precision,intent(inout)      :: a( lda, * )
     end subroutine mkl_blas_dsyr_omp_offload_lp64

     subroutine mkl_blas_ssyr_omp_offload_lp64 ( uplo, n, alpha, x, incx,   &
          &a, lda ) BIND(C)
       character*1,intent(in)              :: uplo
       real,intent(in)                     :: alpha
       integer,intent(in)                  :: incx, lda, n
       real,intent(in)                     :: x( * )
       real,intent(inout)                  :: a( lda, * )
     end subroutine mkl_blas_ssyr_omp_offload_lp64

     subroutine mkl_blas_zher_omp_offload_lp64 ( uplo, n, alpha, x, incx,   &
          &a, lda ) BIND(C)
       character*1,intent(in)              :: uplo
       double precision,intent(in)         :: alpha
       integer,intent(in)                  :: incx, lda, n
       complex*16,intent(in)               :: x( * )
       complex*16,intent(inout)            :: a( lda, * )
     end subroutine mkl_blas_zher_omp_offload_lp64

     subroutine mkl_blas_cher_omp_offload_lp64 ( uplo, n, alpha, x, incx,   &
          &a, lda ) BIND(C)
       character*1,intent(in)              :: uplo
       real,intent(in)                     :: alpha
       integer,intent(in)                  :: incx, lda, n
       complex,intent(in)                  :: x( * )
       complex,intent(inout)               :: a( lda, * )
     end subroutine mkl_blas_cher_omp_offload_lp64

     subroutine mkl_blas_dspr_omp_offload_lp64 ( uplo, n, alpha, x, incx,   &
          &a ) BIND(C)
       character*1,intent(in)              :: uplo
       double precision,intent(in)         :: alpha
       integer,intent(in)                  :: incx, n
       double precision,intent(in)         :: x( * )
       double precision,intent(inout)      :: a( * )
     end subroutine mkl_blas_dspr_omp_offload_lp64

     subroutine mkl_blas_sspr_omp_offload_lp64 ( uplo, n, alpha, x, incx,   &
          &a ) BIND(C)
       character*1,intent(in)              :: uplo
       real,intent(in)                     :: alpha
       integer,intent(in)                  :: incx, n
       real,intent(in)                     :: x( * )
       real,intent(inout)                  :: a( * )
     end subroutine mkl_blas_sspr_omp_offload_lp64

     subroutine mkl_blas_zhpr_omp_offload_lp64 ( uplo, n, alpha, x, incx,   &
          &a ) BIND(C)
       character*1,intent(in)              :: uplo
       double precision,intent(in)         :: alpha
       integer,intent(in)                  :: incx, n
       complex*16,intent(in)               :: x( * )
       complex*16,intent(inout)            :: a( * )
     end subroutine mkl_blas_zhpr_omp_offload_lp64

     subroutine mkl_blas_chpr_omp_offload_lp64 ( uplo, n, alpha, x, incx,   &
          &a ) BIND(C)
       character*1,intent(in)              :: uplo
       real,intent(in)                  :: alpha
       integer,intent(in)                  :: incx, n
       complex,intent(in)                  :: x( * )
       complex,intent(inout)               :: a( * )
     end subroutine mkl_blas_chpr_omp_offload_lp64

     subroutine mkl_blas_dsyr2_omp_offload_lp64 ( uplo, n, alpha, x, incx, y, incy,   &
          &a, lda ) BIND(C)
       character*1,intent(in)              :: uplo
       double precision,intent(in)         :: alpha
       integer,intent(in)                  :: incx, incy, lda, n
       double precision,intent(in)         :: x( * ), y( * )
       double precision,intent(inout)      :: a( lda, * )
     end subroutine mkl_blas_dsyr2_omp_offload_lp64

     subroutine mkl_blas_ssyr2_omp_offload_lp64 ( uplo, n, alpha, x, incx, y, incy,   &
          &a, lda ) BIND(C)
       character*1,intent(in)              :: uplo
       real,intent(in)                     :: alpha
       integer,intent(in)                  :: incx, incy, lda, n
       real,intent(in)                     :: x( * ), y( * )
       real,intent(inout)                  :: a( lda, * )
     end subroutine mkl_blas_ssyr2_omp_offload_lp64

     subroutine mkl_blas_zher2_omp_offload_lp64 ( uplo, n, alpha, x, incx, y, incy,   &
          &a, lda ) BIND(C)
       character*1,intent(in)              :: uplo
       complex*16,intent(in)               :: alpha
       integer,intent(in)                  :: incx, incy, lda, n
       complex*16,intent(in)               :: x( * ), y( * )
       complex*16,intent(inout)            :: a( lda, * )
     end subroutine mkl_blas_zher2_omp_offload_lp64

     subroutine mkl_blas_cher2_omp_offload_lp64 ( uplo, n, alpha, x, incx, y, incy,   &
          &a, lda ) BIND(C)
       character*1,intent(in)              :: uplo
       complex,intent(in)                  :: alpha
       integer,intent(in)                  :: incx, incy, lda, n
       complex,intent(in)                  :: x( * ), y( * )
       complex,intent(inout)               :: a( lda, * )
     end subroutine mkl_blas_cher2_omp_offload_lp64

     subroutine mkl_blas_dspr2_omp_offload_lp64 ( uplo, n, alpha, x, incx, y, incy,   &
          &a ) BIND(C)
       character*1,intent(in)              :: uplo
       double precision,intent(in)         :: alpha
       integer,intent(in)                  :: incx, incy, n
       double precision,intent(in)         :: x( * ), y( * )
       double precision,intent(inout)      :: a( * )
     end subroutine mkl_blas_dspr2_omp_offload_lp64

     subroutine mkl_blas_sspr2_omp_offload_lp64 ( uplo, n, alpha, x, incx, y, incy,   &
          &a ) BIND(C)
       character*1,intent(in)              :: uplo
       real,intent(in)                     :: alpha
       integer,intent(in)                  :: incx, incy, n
       real,intent(in)                     :: x( * ), y( * )
       real,intent(inout)                  :: a( * )
     end subroutine mkl_blas_sspr2_omp_offload_lp64

     subroutine mkl_blas_zhpr2_omp_offload_lp64 ( uplo, n, alpha, x, incx, y, incy,   &
          &a ) BIND(C)
       character*1,intent(in)              :: uplo
       complex*16,intent(in)               :: alpha
       integer,intent(in)                  :: incx, incy, n
       complex*16,intent(in)               :: x( * ), y( * )
       complex*16,intent(inout)            :: a( * )
     end subroutine mkl_blas_zhpr2_omp_offload_lp64

     subroutine mkl_blas_chpr2_omp_offload_lp64 ( uplo, n, alpha, x, incx, y, incy,   &
          &a ) BIND(C)
       character*1,intent(in)              :: uplo
       complex,intent(in)                  :: alpha
       integer,intent(in)                  :: incx, incy, n
       complex,intent(in)                  :: x( * ), y( * )
       complex,intent(inout)               :: a( * )
     end subroutine mkl_blas_chpr2_omp_offload_lp64

     subroutine mkl_blas_dtbmv_omp_offload_lp64 ( uplo, trans, diag, n, k, a, lda, x, incx ) BIND(C)
       integer,intent(in)                  :: incx, k, lda, n
       character*1,intent(in)              :: uplo, trans, diag
       double precision,intent(in)         :: a( lda, * )
       double precision,intent(inout)      :: x( * )
     end subroutine mkl_blas_dtbmv_omp_offload_lp64

     subroutine mkl_blas_stbmv_omp_offload_lp64 ( uplo, trans, diag, n, k, a, lda, x, incx ) BIND(C)
       integer,intent(in)                  :: incx, k, lda, n
       character*1,intent(in)              :: uplo, trans, diag
       real,intent(in)                     :: a( lda, * )
       real,intent(inout)                  :: x( * )
     end subroutine mkl_blas_stbmv_omp_offload_lp64

     subroutine mkl_blas_ztbmv_omp_offload_lp64 ( uplo, trans, diag, n, k, a, lda, x, incx ) BIND(C)
       integer,intent(in)                  :: incx, k, lda, n
       character*1,intent(in)              :: uplo, trans, diag
       complex*16,intent(in)               :: a( lda, * )
       complex*16,intent(inout)            :: x( * )
     end subroutine mkl_blas_ztbmv_omp_offload_lp64

     subroutine mkl_blas_ctbmv_omp_offload_lp64 ( uplo, trans, diag, n, k, a, lda, x, incx ) BIND(C)
       integer,intent(in)                  :: incx, k, lda, n
       character*1,intent(in)              :: uplo, trans, diag
       complex,intent(in)                  :: a( lda, * )
       complex,intent(inout)               :: x( * )
     end subroutine mkl_blas_ctbmv_omp_offload_lp64

     subroutine mkl_blas_dtrmv_omp_offload_lp64 ( uplo, trans, diag, n, a, lda, x, incx ) BIND(C)
       integer,intent(in)                  :: incx, lda, n
       character*1,intent(in)              :: uplo, trans, diag
       double precision,intent(in)         :: a( lda, * )
       double precision,intent(inout)      :: x( * )
     end subroutine mkl_blas_dtrmv_omp_offload_lp64

     subroutine mkl_blas_strmv_omp_offload_lp64 ( uplo, trans, diag, n, a, lda, x, incx ) BIND(C)
       integer,intent(in)                  :: incx, lda, n
       character*1,intent(in)              :: uplo, trans, diag
       real,intent(in)                     :: a( lda, * )
       real,intent(inout)                  :: x( * )
     end subroutine mkl_blas_strmv_omp_offload_lp64

     subroutine mkl_blas_ztrmv_omp_offload_lp64 ( uplo, trans, diag, n, a, lda, x, incx ) BIND(C)
       integer,intent(in)                  :: incx, lda, n
       character*1,intent(in)              :: uplo, trans, diag
       complex*16,intent(in)               :: a( lda, * )
       complex*16,intent(inout)            :: x( * )
     end subroutine mkl_blas_ztrmv_omp_offload_lp64

     subroutine mkl_blas_ctrmv_omp_offload_lp64 ( uplo, trans, diag, n, a, lda, x, incx ) BIND(C)
       integer,intent(in)                  :: incx, lda, n
       character*1,intent(in)              :: uplo, trans, diag
       complex,intent(in)                  :: a( lda, * )
       complex,intent(inout)               :: x( * )
     end subroutine mkl_blas_ctrmv_omp_offload_lp64

     subroutine mkl_blas_dtpmv_omp_offload_lp64 ( uplo, trans, diag, n, a, x, incx ) BIND(C)
       integer,intent(in)                  :: incx, n
       character*1,intent(in)              :: uplo, trans, diag
       double precision,intent(in)         :: a( * )
       double precision,intent(inout)      :: x( * )
     end subroutine mkl_blas_dtpmv_omp_offload_lp64

     subroutine mkl_blas_stpmv_omp_offload_lp64 ( uplo, trans, diag, n, a, x, incx ) BIND(C)
       integer,intent(in)                  :: incx, n
       character*1,intent(in)              :: uplo, trans, diag
       real,intent(in)                     :: a( * )
       real,intent(inout)                  :: x( * )
     end subroutine mkl_blas_stpmv_omp_offload_lp64

     subroutine mkl_blas_ztpmv_omp_offload_lp64 ( uplo, trans, diag, n, a, x, incx ) BIND(C)
       integer,intent(in)                  :: incx, n
       character*1,intent(in)              :: uplo, trans, diag
       complex*16,intent(in)               :: a( * )
       complex*16,intent(inout)            :: x( * )
     end subroutine mkl_blas_ztpmv_omp_offload_lp64

     subroutine mkl_blas_ctpmv_omp_offload_lp64 ( uplo, trans, diag, n, a, x, incx ) BIND(C)
       integer,intent(in)                  :: incx, n
       character*1,intent(in)              :: uplo, trans, diag
       complex,intent(in)                  :: a( * )
       complex,intent(inout)               :: x( * )
     end subroutine mkl_blas_ctpmv_omp_offload_lp64
     
     subroutine mkl_blas_dtbsv_omp_offload_lp64 ( uplo, trans, diag, n, k, a, lda, x, incx ) BIND(C)
       integer,intent(in)                  :: incx, k, lda, n
       character*1,intent(in)              :: uplo, trans, diag
       double precision,intent(in)         :: a( lda, * )
       double precision,intent(inout)      :: x( * )
     end subroutine mkl_blas_dtbsv_omp_offload_lp64

     subroutine mkl_blas_stbsv_omp_offload_lp64 ( uplo, trans, diag, n, k, a, lda, x, incx ) BIND(C)
       integer,intent(in)                  :: incx, k, lda, n
       character*1,intent(in)              :: uplo, trans, diag
       real,intent(in)                     :: a( lda, * )
       real,intent(inout)                  :: x( * )
     end subroutine mkl_blas_stbsv_omp_offload_lp64

     subroutine mkl_blas_ztbsv_omp_offload_lp64 ( uplo, trans, diag, n, k, a, lda, x, incx ) BIND(C)
       integer,intent(in)                  :: incx, k, lda, n
       character*1,intent(in)              :: uplo, trans, diag
       complex*16,intent(in)               :: a( lda, * )
       complex*16,intent(inout)            :: x( * )
     end subroutine mkl_blas_ztbsv_omp_offload_lp64

     subroutine mkl_blas_ctbsv_omp_offload_lp64 ( uplo, trans, diag, n, k, a, lda, x, incx ) BIND(C)
       integer,intent(in)                  :: incx, k, lda, n
       character*1,intent(in)              :: uplo, trans, diag
       complex,intent(in)                  :: a( lda, * )
       complex,intent(inout)               :: x( * )
     end subroutine mkl_blas_ctbsv_omp_offload_lp64

     subroutine mkl_blas_dtrsv_omp_offload_lp64 ( uplo, trans, diag, n, a, lda, x, incx ) BIND(C)
       integer,intent(in)                  :: incx, lda, n
       character*1,intent(in)              :: uplo, trans, diag
       double precision,intent(in)         :: a( lda, * )
       double precision,intent(inout)      :: x( * )
     end subroutine mkl_blas_dtrsv_omp_offload_lp64

     subroutine mkl_blas_strsv_omp_offload_lp64 ( uplo, trans, diag, n, a, lda, x, incx ) BIND(C)
       integer,intent(in)                  :: incx, lda, n
       character*1,intent(in)              :: uplo, trans, diag
       real,intent(in)                     :: a( lda, * )
       real,intent(inout)                  :: x( * )
     end subroutine mkl_blas_strsv_omp_offload_lp64

     subroutine mkl_blas_ztrsv_omp_offload_lp64 ( uplo, trans, diag, n, a, lda, x, incx ) BIND(C)
       integer,intent(in)                  :: incx, lda, n
       character*1,intent(in)              :: uplo, trans, diag
       complex*16,intent(in)               :: a( lda, * )
       complex*16,intent(inout)            :: x( * )
     end subroutine mkl_blas_ztrsv_omp_offload_lp64

     subroutine mkl_blas_ctrsv_omp_offload_lp64 ( uplo, trans, diag, n, a, lda, x, incx ) BIND(C)
       integer,intent(in)                  :: incx, lda, n
       character*1,intent(in)              :: uplo, trans, diag
       complex,intent(in)                  :: a( lda, * )
       complex,intent(inout)               :: x( * )
     end subroutine mkl_blas_ctrsv_omp_offload_lp64

     subroutine mkl_blas_dtpsv_omp_offload_lp64 ( uplo, trans, diag, n, a, x, incx ) BIND(C)
       integer,intent(in)                  :: incx, n
       character*1,intent(in)              :: uplo, trans, diag
       double precision,intent(in)         :: a( * )
       double precision,intent(inout)      :: x( * )
     end subroutine mkl_blas_dtpsv_omp_offload_lp64

     subroutine mkl_blas_stpsv_omp_offload_lp64 ( uplo, trans, diag, n, a, x, incx ) BIND(C)
       integer,intent(in)                  :: incx, n
       character*1,intent(in)              :: uplo, trans, diag
       real,intent(in)                     :: a( * )
       real,intent(inout)                  :: x( * )
     end subroutine mkl_blas_stpsv_omp_offload_lp64

     subroutine mkl_blas_ztpsv_omp_offload_lp64 ( uplo, trans, diag, n, a, x, incx ) BIND(C)
       integer,intent(in)                  :: incx, n
       character*1,intent(in)              :: uplo, trans, diag
       complex*16,intent(in)               :: a( * )
       complex*16,intent(inout)            :: x( * )
     end subroutine mkl_blas_ztpsv_omp_offload_lp64

     subroutine mkl_blas_ctpsv_omp_offload_lp64 ( uplo, trans, diag, n, a, x, incx ) BIND(C)
       integer,intent(in)                  :: incx, n
       character*1,intent(in)              :: uplo, trans, diag
       complex,intent(in)                  :: a( * )
       complex,intent(inout)               :: x( * )
     end subroutine mkl_blas_ctpsv_omp_offload_lp64

     ! BLAS Level1

     function mkl_blas_dasum_omp_offload_lp64 ( n, x, incx ) BIND(C)
       double precision                    :: mkl_blas_dasum_omp_offload_lp64
       integer,intent(in)                  :: incx, n
       double precision,intent(in)         :: x( * )
     end function mkl_blas_dasum_omp_offload_lp64

     function mkl_blas_sasum_omp_offload_lp64 ( n, x, incx ) BIND(C)
       real                                :: mkl_blas_sasum_omp_offload_lp64
       integer,intent(in)                  :: incx, n
       real,intent(in)                     :: x( * )
     end function mkl_blas_sasum_omp_offload_lp64

     function mkl_blas_dzasum_omp_offload_lp64 ( n, x, incx ) BIND(C)
       double precision                    :: mkl_blas_dzasum_omp_offload_lp64
       integer,intent(in)                  :: incx, n
       complex*16,intent(in)               :: x( * )
     end function mkl_blas_dzasum_omp_offload_lp64

     function mkl_blas_scasum_omp_offload_lp64 ( n, x, incx ) BIND(C)
       real                                :: mkl_blas_scasum_omp_offload_lp64
       integer,intent(in)                  :: incx, n
       complex,intent(in)                  :: x( * )
     end function mkl_blas_scasum_omp_offload_lp64

     function mkl_blas_dnrm2_omp_offload_lp64 ( n, x, incx ) BIND(C)
       double precision                    :: mkl_blas_dnrm2_omp_offload_lp64
       integer,intent(in)                  :: incx, n
       double precision,intent(in)         :: x( * )
     end function mkl_blas_dnrm2_omp_offload_lp64

     function mkl_blas_snrm2_omp_offload_lp64 ( n, x, incx ) BIND(C)
       real                                :: mkl_blas_snrm2_omp_offload_lp64
       integer,intent(in)                  :: incx, n
       real,intent(in)                     :: x( * )
     end function mkl_blas_snrm2_omp_offload_lp64

     function mkl_blas_dznrm2_omp_offload_lp64 ( n, x, incx ) BIND(C)
       double precision                    :: mkl_blas_dznrm2_omp_offload_lp64
       integer,intent(in)                  :: incx, n
       complex*16,intent(in)               :: x( * )
     end function mkl_blas_dznrm2_omp_offload_lp64

     function mkl_blas_scnrm2_omp_offload_lp64 ( n, x, incx ) BIND(C)
       real                                :: mkl_blas_scnrm2_omp_offload_lp64
       integer,intent(in)                  :: incx, n
       complex,intent(in)                  :: x( * )
     end function mkl_blas_scnrm2_omp_offload_lp64

     subroutine mkl_blas_daxpy_omp_offload_lp64 ( n, alpha, x, incx,   &
          &y, incy ) BIND(C)
       double precision,intent(in)         :: alpha
       integer,intent(in)                  :: incx, incy, n
       double precision,intent(in)         :: x( * )
       double precision,intent(inout)      :: y( * )
     end subroutine mkl_blas_daxpy_omp_offload_lp64

     subroutine mkl_blas_saxpy_omp_offload_lp64 ( n, alpha, x, incx,   &
          &y, incy ) BIND(C)
       real,intent(in)                     :: alpha
       integer,intent(in)                  :: incx, incy, n
       real,intent(in)                     :: x( * )
       real,intent(inout)                  :: y( * )
     end subroutine mkl_blas_saxpy_omp_offload_lp64

     subroutine mkl_blas_zaxpy_omp_offload_lp64 ( n, alpha, x, incx,   &
          &y, incy ) BIND(C)
       complex*16,intent(in)               :: alpha
       integer,intent(in)                  :: incx, incy, n
       complex*16,intent(in)               :: x( * )
       complex*16,intent(inout)            :: y( * )
     end subroutine mkl_blas_zaxpy_omp_offload_lp64

     subroutine mkl_blas_caxpy_omp_offload_lp64 ( n, alpha, x, incx,   &
          &y, incy ) BIND(C)
       complex,intent(in)                  :: alpha
       integer,intent(in)                  :: incx, incy, n
       complex,intent(in)                  :: x( * )
       complex,intent(inout)               :: y( * )
     end subroutine mkl_blas_caxpy_omp_offload_lp64

     subroutine mkl_blas_daxpby_omp_offload_lp64 ( n, alpha, x, incx,   &
          &beta, y, incy ) BIND(C)
       double precision,intent(in)         :: alpha, beta
       integer,intent(in)                  :: incx, incy, n
       double precision,intent(in)         :: x( * )
       double precision,intent(inout)      :: y( * )
     end subroutine mkl_blas_daxpby_omp_offload_lp64

     subroutine mkl_blas_saxpby_omp_offload_lp64 ( n, alpha, x, incx,   &
          &beta, y, incy ) BIND(C)
       real,intent(in)                     :: alpha, beta
       integer,intent(in)                  :: incx, incy, n
       real,intent(in)                     :: x( * )
       real,intent(inout)                  :: y( * )
     end subroutine mkl_blas_saxpby_omp_offload_lp64

     subroutine mkl_blas_zaxpby_omp_offload_lp64 ( n, alpha, x, incx,   &
          &beta, y, incy ) BIND(C)
       complex*16,intent(in)               :: alpha, beta
       integer,intent(in)                  :: incx, incy, n
       complex*16,intent(in)               :: x( * )
       complex*16,intent(inout)            :: y( * )
     end subroutine mkl_blas_zaxpby_omp_offload_lp64

     subroutine mkl_blas_caxpby_omp_offload_lp64 ( n, alpha, x, incx,   &
          &beta, y, incy ) BIND(C)
       complex,intent(in)                  :: alpha, beta
       integer,intent(in)                  :: incx, incy, n
       complex,intent(in)                  :: x( * )
       complex,intent(inout)               :: y( * )
     end subroutine mkl_blas_caxpby_omp_offload_lp64

     subroutine mkl_blas_dcopy_omp_offload_lp64 ( n, x, incx,   &
          &y, incy ) BIND(C)
       integer,intent(in)                  :: incx, incy, n
       double precision,intent(in)         :: x( * )
       double precision,intent(inout)      :: y( * )
     end subroutine mkl_blas_dcopy_omp_offload_lp64

     subroutine mkl_blas_scopy_omp_offload_lp64 ( n, x, incx,   &
          &y, incy ) BIND(C)
       integer,intent(in)                  :: incx, incy, n
       real,intent(in)                     :: x( * )
       real,intent(inout)                  :: y( * )
     end subroutine mkl_blas_scopy_omp_offload_lp64

     subroutine mkl_blas_zcopy_omp_offload_lp64 ( n, x, incx,   &
          &y, incy ) BIND(C)
       integer,intent(in)                  :: incx, incy, n
       complex*16,intent(in)               :: x( * )
       complex*16,intent(inout)            :: y( * )
     end subroutine mkl_blas_zcopy_omp_offload_lp64

     subroutine mkl_blas_ccopy_omp_offload_lp64 ( n, x, incx,   &
          &y, incy ) BIND(C)
       integer,intent(in)                  :: incx, incy, n
       complex,intent(in)                  :: x( * )
       complex,intent(inout)               :: y( * )
     end subroutine mkl_blas_ccopy_omp_offload_lp64

     function mkl_blas_ddot_omp_offload_lp64 ( n, x, incx, y, incy ) BIND(C)
       double precision                    :: mkl_blas_ddot_omp_offload_lp64
       integer,intent(in)                  :: incx, incy, n
       double precision,intent(in)         :: x( * ), y( * )
     end function mkl_blas_ddot_omp_offload_lp64

     function mkl_blas_sdot_omp_offload_lp64 ( n, x, incx, y, incy ) BIND(C)
       real                                :: mkl_blas_sdot_omp_offload_lp64
       integer,intent(in)                  :: incx, incy, n
       real,intent(in)                     :: x( * ), y( * )
     end function mkl_blas_sdot_omp_offload_lp64

     function mkl_blas_dsdot_omp_offload_lp64 ( n, x, incx, y, incy ) BIND(C)
       double precision                    :: mkl_blas_dsdot_omp_offload_lp64
       integer,intent(in)                  :: incx, incy, n
       real,intent(in)                     :: x( * ), y( * )
     end function mkl_blas_dsdot_omp_offload_lp64

     function mkl_blas_sdsdot_omp_offload_lp64 ( n, sb, x, incx, y, incy ) BIND(C)
       real                                :: mkl_blas_sdsdot_omp_offload_lp64
       integer,intent(in)                  :: incx, incy, n
       real,intent(in)                     :: sb, x( * ), y( * )
     end function mkl_blas_sdsdot_omp_offload_lp64

     function mkl_blas_zdotc_omp_offload_lp64 ( n, x, incx, y, incy )
       complex*16                          :: mkl_blas_zdotc_omp_offload_lp64
       integer,intent(in)                  :: incx, incy, n
       complex*16,intent(in)               :: x( * ), y( * )
     end function mkl_blas_zdotc_omp_offload_lp64

     function mkl_blas_cdotc_omp_offload_lp64 ( n, x, incx, y, incy )
       complex                             :: mkl_blas_cdotc_omp_offload_lp64
       integer,intent(in)                  :: incx, incy, n
       complex,intent(in)                  :: x( * ), y( * )
     end function mkl_blas_cdotc_omp_offload_lp64

     function mkl_blas_zdotu_omp_offload_lp64 ( n, x, incx, y, incy )
       complex*16                          :: mkl_blas_zdotu_omp_offload_lp64
       integer,intent(in)                  :: incx, incy, n
       complex*16,intent(in)               :: x( * ), y( * )
     end function mkl_blas_zdotu_omp_offload_lp64

     function mkl_blas_cdotu_omp_offload_lp64 ( n, x, incx, y, incy )
       complex                             :: mkl_blas_cdotu_omp_offload_lp64
       integer,intent(in)                  :: incx, incy, n
       complex,intent(in)                  :: x( * ), y( * )
     end function mkl_blas_cdotu_omp_offload_lp64

     subroutine mkl_blas_drot_omp_offload_lp64 ( n, x, incx,   &
          &y, incy, c, s ) BIND(C)
       double precision,intent(in)         :: c, s         
       integer,intent(in)                  :: incx, incy, n
       double precision,intent(inout)      :: x( * ), y( * )
     end subroutine mkl_blas_drot_omp_offload_lp64

     subroutine mkl_blas_srot_omp_offload_lp64 ( n, x, incx,   &
          &y, incy, c, s ) BIND(C)
       real,intent(in)                     :: c, s         
       integer,intent(in)                  :: incx, incy, n
       real,intent(inout)                  :: x( * ), y( * )
     end subroutine mkl_blas_srot_omp_offload_lp64

     subroutine mkl_blas_zdrot_omp_offload_lp64 ( n, x, incx,   &
          &y, incy, c, s ) BIND(C)
       double precision,intent(in)         :: c, s         
       integer,intent(in)                  :: incx, incy, n
       complex*16,intent(inout)            :: x( * ), y( * )
     end subroutine mkl_blas_zdrot_omp_offload_lp64

     subroutine mkl_blas_csrot_omp_offload_lp64 ( n, x, incx,   &
          &y, incy, c, s ) BIND(C)
       real,intent(in)                     :: c, s         
       integer,intent(in)                  :: incx, incy, n
       complex,intent(inout)               :: x( * ), y( * )
     end subroutine mkl_blas_csrot_omp_offload_lp64

     subroutine mkl_blas_drotg_omp_offload_lp64 ( a, b, c, s ) BIND(C)
       double precision,intent(inout)      :: a( * ), b( * ), c( * ), s( * )
     end subroutine mkl_blas_drotg_omp_offload_lp64

     subroutine mkl_blas_srotg_omp_offload_lp64 ( a, b, c, s ) BIND(C)
       real,intent(inout)                  :: a( * ), b( * ), c( * ), s( * )
     end subroutine mkl_blas_srotg_omp_offload_lp64

     subroutine mkl_blas_zrotg_omp_offload_lp64 ( a, b, c, s ) BIND(C)
       complex*16,intent(inout)            :: a( * ), b( * ), s( * )
       double precision,intent(inout)      :: c( * )
     end subroutine mkl_blas_zrotg_omp_offload_lp64

     subroutine mkl_blas_crotg_omp_offload_lp64 ( a, b, c, s ) BIND(C)
       complex,intent(inout)               :: a( * ), b( * ), s( * )
       real,intent(inout)                  :: c( * )
     end subroutine mkl_blas_crotg_omp_offload_lp64

     subroutine mkl_blas_drotm_omp_offload_lp64 ( n, x, incx,   &
          &y, incy, param ) BIND(C)
       integer,intent(in)                  :: incx, incy, n
       double precision,intent(in)         :: param( * )
       double precision,intent(inout)      :: x( * ), y( * )
     end subroutine mkl_blas_drotm_omp_offload_lp64

     subroutine mkl_blas_srotm_omp_offload_lp64 ( n, x, incx,   &
          &y, incy, param ) BIND(C)
       integer,intent(in)                  :: incx, incy, n
       real,intent(in)                     :: param( * )
       real,intent(inout)                  :: x( * ), y( * )
     end subroutine mkl_blas_srotm_omp_offload_lp64

     subroutine mkl_blas_drotmg_omp_offload_lp64 ( d1, d2, x1, x2, param ) BIND(C)
       double precision,intent(in)         :: x2
       double precision,intent(inout)      :: d1( * ), d2( * ), x1( * ), param( * )
     end subroutine mkl_blas_drotmg_omp_offload_lp64

     subroutine mkl_blas_srotmg_omp_offload_lp64 ( d1, d2, x1, x2, param ) BIND(C)
       real,intent(in)                     :: x2
       real,intent(inout)                  :: d1( * ), d2( * ), x1( * ), param( * )
     end subroutine mkl_blas_srotmg_omp_offload_lp64

     subroutine mkl_blas_dscal_omp_offload_lp64 ( n, alpha, x, incx ) BIND(C)
       double precision,intent(in)         :: alpha
       integer,intent(in)                  :: incx, n
       double precision,intent(inout)      :: x( * )
     end subroutine mkl_blas_dscal_omp_offload_lp64

     subroutine mkl_blas_sscal_omp_offload_lp64 ( n, alpha, x, incx ) BIND(C)
       real,intent(in)                     :: alpha
       integer,intent(in)                  :: incx, n
       real,intent(inout)                  :: x( * )
     end subroutine mkl_blas_sscal_omp_offload_lp64

     subroutine mkl_blas_zscal_omp_offload_lp64 ( n, alpha, x, incx ) BIND(C)
       complex*16,intent(in)               :: alpha
       integer,intent(in)                  :: incx, n
       complex*16,intent(inout)            :: x( * )
     end subroutine mkl_blas_zscal_omp_offload_lp64

     subroutine mkl_blas_cscal_omp_offload_lp64 ( n, alpha, x, incx ) BIND(C)
       complex,intent(in)                  :: alpha
       integer,intent(in)                  :: incx, n
       complex,intent(inout)               :: x( * )
     end subroutine mkl_blas_cscal_omp_offload_lp64

     subroutine mkl_blas_zdscal_omp_offload_lp64 ( n, alpha, x, incx ) BIND(C)
       double precision,intent(in)         :: alpha
       integer,intent(in)                  :: incx, n
       complex*16,intent(inout)            :: x( * )
     end subroutine mkl_blas_zdscal_omp_offload_lp64

     subroutine mkl_blas_csscal_omp_offload_lp64 ( n, alpha, x, incx ) BIND(C)
       real,intent(in)                     :: alpha
       integer,intent(in)                  :: incx, n
       complex,intent(inout)               :: x( * )
     end subroutine mkl_blas_csscal_omp_offload_lp64

     subroutine mkl_blas_dswap_omp_offload_lp64 ( n, x, incx,   &
          &y, incy ) BIND(C)
       integer,intent(in)                  :: incx, incy, n
       double precision,intent(inout)      :: x( * ), y( * )
     end subroutine mkl_blas_dswap_omp_offload_lp64

     subroutine mkl_blas_sswap_omp_offload_lp64 ( n, x, incx,   &
          &y, incy ) BIND(C)
       integer,intent(in)                  :: incx, incy, n
       real,intent(inout)                  :: x( * ), y( * )
     end subroutine mkl_blas_sswap_omp_offload_lp64

     subroutine mkl_blas_zswap_omp_offload_lp64 ( n, x, incx,   &
          &y, incy ) BIND(C)
       integer,intent(in)                  :: incx, incy, n
       complex*16,intent(inout)            :: x( * ), y( * )
     end subroutine mkl_blas_zswap_omp_offload_lp64

     subroutine mkl_blas_cswap_omp_offload_lp64 ( n, x, incx,   &
          &y, incy ) BIND(C)
       integer,intent(in)                  :: incx, incy, n
       complex,intent(inout)               :: x( * ), y( * )
     end subroutine mkl_blas_cswap_omp_offload_lp64

     function mkl_blas_idamax_omp_offload_lp64 ( n, x, incx ) BIND(C)
       integer                             :: mkl_blas_idamax_omp_offload_lp64
       integer,intent(in)                  :: incx, n
       double precision,intent(in)         :: x( * )
     end function mkl_blas_idamax_omp_offload_lp64

     function mkl_blas_isamax_omp_offload_lp64 ( n, x, incx ) BIND(C)
       integer                             :: mkl_blas_isamax_omp_offload_lp64
       integer,intent(in)                  :: incx, n
       real,intent(in)                     :: x( * )
     end function mkl_blas_isamax_omp_offload_lp64

     function mkl_blas_izamax_omp_offload_lp64 ( n, x, incx ) BIND(C)
       integer                             :: mkl_blas_izamax_omp_offload_lp64
       integer,intent(in)                  :: incx, n
       complex*16,intent(in)               :: x( * )
     end function mkl_blas_izamax_omp_offload_lp64

     function mkl_blas_icamax_omp_offload_lp64 ( n, x, incx ) BIND(C)
       integer                             :: mkl_blas_icamax_omp_offload_lp64
       integer,intent(in)                  :: incx, n
       complex,intent(in)                  :: x( * )
     end function mkl_blas_icamax_omp_offload_lp64

     function mkl_blas_idamin_omp_offload_lp64 ( n, x, incx ) BIND(C)
       integer                             :: mkl_blas_idamin_omp_offload_lp64
       integer,intent(in)                  :: incx, n
       double precision,intent(in)         :: x( * )
     end function mkl_blas_idamin_omp_offload_lp64

     function mkl_blas_isamin_omp_offload_lp64 ( n, x, incx ) BIND(C)
       integer                             :: mkl_blas_isamin_omp_offload_lp64
       integer,intent(in)                  :: incx, n
       real,intent(in)                     :: x( * )
     end function mkl_blas_isamin_omp_offload_lp64

     function mkl_blas_izamin_omp_offload_lp64 ( n, x, incx ) BIND(C)
       integer                             :: mkl_blas_izamin_omp_offload_lp64
       integer,intent(in)                  :: incx, n
       complex*16,intent(in)               :: x( * )
     end function mkl_blas_izamin_omp_offload_lp64

     function mkl_blas_icamin_omp_offload_lp64 ( n, x, incx ) BIND(C)
       integer                             :: mkl_blas_icamin_omp_offload_lp64
       integer,intent(in)                  :: incx, n
       complex,intent(in)                  :: x( * )
     end function mkl_blas_icamin_omp_offload_lp64

     ! BLAS batch

     subroutine mkl_blas_dtrsm_batch_strided_omp_offload_lp64 ( side, uplo, trans, diag, m, n, alpha, a, lda, stridea,        &
          &b, ldb, strideb, batch_size ) BIND(C)
       character*1,intent(in)             :: side, uplo, trans, diag
       integer,intent(in)                 :: m, n, lda, stridea, ldb, strideb, batch_size
       double precision,intent(in)        :: alpha
       double precision,intent(in)        :: a( lda, * )
       double precision,intent(inout)     :: b( ldb, * )
     end subroutine mkl_blas_dtrsm_batch_strided_omp_offload_lp64

     subroutine mkl_blas_strsm_batch_strided_omp_offload_lp64 ( side, uplo, trans, diag, m, n, alpha, a, lda, stridea,        &
          &b, ldb, strideb, batch_size ) BIND(C)
       character*1,intent(in)              :: side, uplo, trans, diag
       integer,intent(in)                  :: m, n, lda, stridea, ldb, strideb, batch_size
       real,intent(in)                     :: alpha
       real,intent(in)                     :: a( lda, * )
       real,intent(in)                     :: b( ldb, * )
     end subroutine mkl_blas_strsm_batch_strided_omp_offload_lp64

     subroutine mkl_blas_ztrsm_batch_strided_omp_offload_lp64 ( side, uplo, trans, diag, m, n, alpha, a, lda, stridea,        &
          &b, ldb, strideb, batch_size ) BIND(C)
       character*1,intent(in)             :: side, uplo, trans, diag
       integer,intent(in)                 :: m, n, lda, stridea, ldb, strideb, batch_size
       complex*16,intent(in)              :: alpha
       complex*16,intent(in)              :: a( lda, * )
       complex*16,intent(inout)           :: b( ldb, * )
     end subroutine mkl_blas_ztrsm_batch_strided_omp_offload_lp64

     subroutine mkl_blas_ctrsm_batch_strided_omp_offload_lp64 ( side, uplo, trans, diag, m, n, alpha, a, lda, stridea,        &
          &b, ldb, strideb, batch_size ) BIND(C)
       character*1,intent(in)              :: side, uplo, trans, diag
       integer,intent(in)                  :: m, n, lda, stridea, ldb, strideb, batch_size
       complex,intent(in)                  :: alpha
       complex,intent(in)                  :: a( lda, * )
       complex,intent(in)                  :: b( ldb, * )
     end subroutine mkl_blas_ctrsm_batch_strided_omp_offload_lp64

     subroutine mkl_blas_dgemm_batch_strided_omp_offload_lp64 ( transa, transb, m, n, k, alpha, a, lda, stridea,        &
          &b, ldb, strideb, beta, c, ldc, stridec, batch_size ) BIND(C)
       character*1,intent(in)             :: transa, transb
       integer,intent(in)                 :: m, n, k, lda, stridea, ldb, strideb, ldc, stridec, batch_size
       double precision,intent(in)        :: alpha, beta
       double precision,intent(in)        :: a( lda, * ), b( ldb, * )
       double precision,intent(inout)     :: c( ldc, * )
     end subroutine mkl_blas_dgemm_batch_strided_omp_offload_lp64

     subroutine mkl_blas_sgemm_batch_strided_omp_offload_lp64 ( transa, transb, m, n, k, alpha, a, lda, stridea,        &
          &b, ldb, strideb, beta, c, ldc, stridec, batch_size ) BIND(C)
       character*1,intent(in)              :: transa, transb
       integer,intent(in)                  :: m, n, k, lda, stridea, ldb, strideb, ldc, stridec, batch_size
       real,intent(in)                     :: alpha, beta
       real,intent(in)                     :: a( lda, * ), b( ldb, * )
       real,intent(in)                     :: c( ldc, * )
     end subroutine mkl_blas_sgemm_batch_strided_omp_offload_lp64

     subroutine mkl_blas_zgemm_batch_strided_omp_offload_lp64 ( transa, transb, m, n, k, alpha, a, lda, stridea,        &
          &b, ldb, strideb, beta, c, ldc, stridec, batch_size ) BIND(C)
       character*1,intent(in)             :: transa, transb
       integer,intent(in)                 :: m, n, k, lda, stridea, ldb, strideb, ldc, stridec, batch_size
       complex*16,intent(in)              :: alpha, beta
       complex*16,intent(in)              :: a( lda, * ), b( ldb, * )
       complex*16,intent(inout)           :: c( ldc, * )
     end subroutine mkl_blas_zgemm_batch_strided_omp_offload_lp64

     subroutine mkl_blas_cgemm_batch_strided_omp_offload_lp64 ( transa, transb, m, n, k, alpha, a, lda, stridea,        &
          &b, ldb, strideb, beta, c, ldc, stridec, batch_size ) BIND(C)
       character*1,intent(in)              :: transa, transb
       integer,intent(in)                  :: m, n, k, lda, stridea, ldb, strideb, ldc, stridec, batch_size
       complex,intent(in)                  :: alpha, beta
       complex,intent(in)                  :: a( lda, * ), b( ldb, * )
       complex,intent(in)                  :: c( ldc, * )
     end subroutine mkl_blas_cgemm_batch_strided_omp_offload_lp64

    subroutine mkl_blas_daxpy_batch_strided_omp_offload_lp64 ( n, alpha, x, incx, stridex,   &
          &y, incy, stridey, batch_size ) BIND(C)
       double precision,intent(in)         :: alpha
       integer,intent(in)                  :: incx, stridex, incy, stridey, batch_size, n
       double precision,intent(in)         :: x( * )
       double precision,intent(inout)      :: y( * )
     end subroutine mkl_blas_daxpy_batch_strided_omp_offload_lp64

     subroutine mkl_blas_saxpy_batch_strided_omp_offload_lp64 ( n, alpha, x, incx, stridex,   &
          &y, incy, stridey, batch_size ) BIND(C)
       real,intent(in)                     :: alpha
       integer,intent(in)                  :: incx, stridex, incy, stridey, batch_size, n
       real,intent(in)                     :: x( * )
       real,intent(inout)                  :: y( * )
     end subroutine mkl_blas_saxpy_batch_strided_omp_offload_lp64

     subroutine mkl_blas_zaxpy_batch_strided_omp_offload_lp64 ( n, alpha, x, incx, stridex,   &
          &y, incy, stridey, batch_size ) BIND(C)
       complex*16,intent(in)               :: alpha
       integer,intent(in)                  :: incx, stridex, incy, stridey, batch_size, n
       complex*16,intent(in)               :: x( * )
       complex*16,intent(inout)            :: y( * )
     end subroutine mkl_blas_zaxpy_batch_strided_omp_offload_lp64

     subroutine mkl_blas_caxpy_batch_strided_omp_offload_lp64 ( n, alpha, x, incx, stridex,   &
          &y, incy, stridey, batch_size ) BIND(C)
       complex,intent(in)                  :: alpha
       integer,intent(in)                  :: incx, stridex, incy, stridey, batch_size, n
       complex,intent(in)                  :: x( * )
       complex,intent(inout)               :: y( * )
     end subroutine mkl_blas_caxpy_batch_strided_omp_offload_lp64

     subroutine mkl_blas_dgemv_batch_strided_omp_offload_lp64 ( trans, m, n, alpha, a, lda, stridea, x, incx,   &
          &stridex, beta, y, incy, stridey, batch_size ) BIND(C)
       double precision,intent(in)         :: alpha, beta
       integer,intent(in)                  :: incx, incy, lda, m, n
       integer,intent(in)                  :: stridea, stridex, stridey, batch_size
       character*1,intent(in)              :: trans
       double precision,intent(in)         :: a( lda, * ), x( * )
       double precision,intent(inout)      :: y( * )
     end subroutine mkl_blas_dgemv_batch_strided_omp_offload_lp64

     subroutine mkl_blas_sgemv_batch_strided_omp_offload_lp64 ( trans, m, n, alpha, a, lda, stridea, x, incx,   &
          &stridex, beta, y, incy, stridey, batch_size ) BIND(C)
       real,intent(in)                     :: alpha, beta
       integer,intent(in)                  :: incx, incy, lda, m, n
       integer,intent(in)                  :: stridea, stridex, stridey, batch_size
       character*1,intent(in)              :: trans
       real,intent(in)                     :: a( lda, * ), x( * )
       real,intent(inout)                  :: y( * )
     end subroutine mkl_blas_sgemv_batch_strided_omp_offload_lp64

     subroutine mkl_blas_zgemv_batch_strided_omp_offload_lp64 ( trans, m, n, alpha, a, lda, stridea, x, incx,   &
          &stridex, beta, y, incy, stridey, batch_size ) BIND(C)
       complex*16,intent(in)               :: alpha, beta
       integer,intent(in)                  :: incx, incy, lda, m, n
       integer,intent(in)                  :: stridea, stridex, stridey, batch_size
       character*1,intent(in)              :: trans
       complex*16,intent(in)               :: a( lda, * ), x( * )
       complex*16,intent(inout)            :: y( * )
     end subroutine mkl_blas_zgemv_batch_strided_omp_offload_lp64

     subroutine mkl_blas_cgemv_batch_strided_omp_offload_lp64 ( trans, m, n, alpha, a, lda, stridea, x, incx,   &
          &stridex, beta, y, incy, stridey, batch_size ) BIND(C)
       complex,intent(in)                  :: alpha, beta
       integer,intent(in)                  :: incx, incy, lda, m, n
       integer,intent(in)                  :: stridea, stridex, stridey, batch_size
       character*1,intent(in)              :: trans
       complex,intent(in)                  :: a( lda, * ), x( * )
       complex,intent(inout)               :: y( * )
     end subroutine mkl_blas_cgemv_batch_strided_omp_offload_lp64

     subroutine mkl_blas_ddgmm_batch_strided_omp_offload_lp64 ( side, m, n, a, lda, stridea, x, incx,   &
          &stridex, c, ldc, stridec, batch_size ) BIND(C)
       integer,intent(in)                  :: incx, ldc, lda, m, n
       integer,intent(in)                  :: stridea, stridex, stridec, batch_size
       character*1,intent(in)              :: side
       double precision,intent(in)         :: a( lda, * ), x( * )
       double precision,intent(inout)      :: c( ldc, * )
     end subroutine mkl_blas_ddgmm_batch_strided_omp_offload_lp64

     subroutine mkl_blas_sdgmm_batch_strided_omp_offload_lp64 ( side, m, n, a, lda, stridea, x, incx,   &
          &stridex, c, ldc, stridec, batch_size ) BIND(C)
       integer,intent(in)                  :: incx, ldc, lda, m, n
       integer,intent(in)                  :: stridea, stridex, stridec, batch_size
       character*1,intent(in)              :: side
       real,intent(in)                     :: a( lda, * ), x( * )
       real,intent(inout)                  :: c( ldc, * )
     end subroutine mkl_blas_sdgmm_batch_strided_omp_offload_lp64

     subroutine mkl_blas_zdgmm_batch_strided_omp_offload_lp64 ( side, m, n, a, lda, stridea, x, incx,   &
          &stridex, c, ldc, stridec, batch_size ) BIND(C)
       integer,intent(in)                  :: incx, ldc, lda, m, n
       integer,intent(in)                  :: stridea, stridex, stridec, batch_size
       character*1,intent(in)              :: side
       complex*16,intent(in)               :: a( lda, * ), x( * )
       complex*16,intent(inout)            :: c( ldc, * )
     end subroutine mkl_blas_zdgmm_batch_strided_omp_offload_lp64

     subroutine mkl_blas_cdgmm_batch_strided_omp_offload_lp64 ( side, m, n, a, lda, stridea, x, incx,   &
          &stridex, c, ldc, stridec, batch_size ) BIND(C)
       integer,intent(in)                  :: incx, ldc, lda, m, n
       integer,intent(in)                  :: stridea, stridex, stridec, batch_size
       character*1,intent(in)              :: side
       complex,intent(in)                  :: a( lda, * ), x( * )
       complex,intent(inout)               :: c( ldc, * )
     end subroutine mkl_blas_cdgmm_batch_strided_omp_offload_lp64

  end interface
