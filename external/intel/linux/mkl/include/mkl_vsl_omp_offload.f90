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
!      Intel(R) oneMKL Library FORTRAN interface for VSL OpenMP offload
!*******************************************************************************

include "mkl_vsl_omp_variant.f90"

module onemkl_vsl_omp_offload
  use onemkl_vsl_omp_variant

  !++
  !  VSL CONTINUOUS DISTRIBUTION GENERATOR FUNCTION INTERFACES.
  !--

  ! Uniform distribution
  interface
     integer function vsrnguniform(method, stream, n, r, a, b)
       import :: VSL_STREAM_STATE
       integer,intent(in)       :: method
       type(VSL_STREAM_STATE)   :: stream
       integer,intent(in)       :: n
       real(kind=4),intent(out) :: r(n)
       real(kind=4),intent(in)  :: a
       real(kind=4),intent(in)  :: b
      !$omp  declare variant( vsrnguniform:mkl_vsl_vsrnguniform_omp_offload ) match( construct={target variant dispatch}, device={arch(gen)} )
     end function vsrnguniform
  end interface

  interface
     integer function vdrnguniform (method, stream, n, r, a, b)
       import :: VSL_STREAM_STATE
       integer,intent(in)       :: method
       type(VSL_STREAM_STATE)   :: stream
       integer,intent(in)       :: n
       real(kind=8),intent(out) :: r(n)
       real(kind=8),intent(in)  :: a
       real(kind=8),intent(in)  :: b
       !$omp  declare variant( vdrnguniform:mkl_vsl_vdrnguniform_omp_offload ) match( construct={target variant dispatch}, device={arch(gen)} )
     end function vdrnguniform
  end interface
  
  ! Gaussian distribution
  interface
     integer function vsrnggaussian (method, stream, n, r, a, sigma)
       import :: VSL_STREAM_STATE
       integer,intent(in)       :: method
       type(VSL_STREAM_STATE)   :: stream
       integer,intent(in)       :: n
       real(kind=4),intent(out) :: r(n)
       real(kind=4),intent(in)  :: a
       real(kind=4),intent(in)  :: sigma
       !$omp  declare variant( vsrnggaussian:mkl_vsl_vsrnggaussian_omp_offload ) match( construct={target variant dispatch}, device={arch(gen)} )
     end function vsrnggaussian
  end interface

  interface
     integer function vdrnggaussian (method, stream, n, r, a, sigma)
       import :: VSL_STREAM_STATE
       integer,intent(in)       :: method
       type(VSL_STREAM_STATE)   :: stream
       integer,intent(in)       :: n
       real(kind=8),intent(out) :: r(n)
       real(kind=8),intent(in)  :: a
       real(kind=8),intent(in)  :: sigma
       !$omp  declare variant( vdrnggaussian:mkl_vsl_vdrnggaussian_omp_offload ) match( construct={target variant dispatch}, device={arch(gen)} )
     end function vdrnggaussian
  end interface

  ! GaussianMV distribution
  interface
     integer function vsrnggaussianmv (method, stream, n, r, dimen, mstorage, a, t)
       import :: VSL_STREAM_STATE
       integer,intent(in)       :: method
       type(VSL_STREAM_STATE)   :: stream
       integer,intent(in)       :: n
       integer, intent(in)      :: dimen
       real(kind=4),intent(out) :: r(dimen,n)
       integer, intent(in)      :: mstorage
       real(kind=4),intent(in)  :: a(dimen)
       real(kind=4),intent(in)  :: t(dimen,dimen)
       !$omp  declare variant( vsrnggaussianmv:mkl_vsl_vsrnggaussianmv_omp_offload ) match( construct={target variant dispatch}, device={arch(gen)} )
     end function vsrnggaussianmv
  end interface

  interface
     integer function vdrnggaussianmv (method, stream, n, r, dimen, mstorage, a, t)
       import :: VSL_STREAM_STATE
       integer,intent(in)       :: method
       type(VSL_STREAM_STATE)   :: stream
       integer,intent(in)       :: n
       integer, intent(in)      :: dimen
       real(kind=8),intent(out) :: r(dimen,n)
       integer, intent(in)      :: mstorage
       real(kind=8),intent(in)  :: a(dimen)
       real(kind=8),intent(in)  :: t(dimen,dimen)
       !$omp  declare variant( vdrnggaussianmv:mkl_vsl_vdrnggaussianmv_omp_offload ) match( construct={target variant dispatch}, device={arch(gen)} )
     end function vdrnggaussianmv
  end interface

  ! Exponential distribution
  interface
     integer function vsrngexponential (method, stream, n, r, a, beta)
       import :: VSL_STREAM_STATE
       integer,intent(in)       :: method
       type(VSL_STREAM_STATE)   :: stream
       integer,intent(in)       :: n
       real(kind=4),intent(out) :: r(n)
       real(kind=4),intent(in)  :: a
       real(kind=4),intent(in)  :: beta
       !$omp  declare variant( vsrngexponential:mkl_vsl_vsrngexponential_omp_offload ) match( construct={target variant dispatch}, device={arch(gen)} )
     end function vsrngexponential
  end interface

  interface
     integer function vdrngexponential (method, stream, n, r, a, beta)
       import :: VSL_STREAM_STATE
       integer,intent(in)       :: method
       type(VSL_STREAM_STATE)   :: stream
       integer,intent(in)       :: n
       real(kind=8),intent(out) :: r(n)
       real(kind=8),intent(in)  :: a
       real(kind=8),intent(in)  :: beta
       !$omp  declare variant( vdrngexponential:mkl_vsl_vdrngexponential_omp_offload ) match( construct={target variant dispatch}, device={arch(gen)} )
     end function vdrngexponential
  end interface

  ! Laplace distribution
  interface
     integer function vsrnglaplace (method, stream, n, r, a, beta)
       import :: VSL_STREAM_STATE
       integer,intent(in)       :: method
       type(VSL_STREAM_STATE)   :: stream
       integer,intent(in)       :: n
       real(kind=4),intent(out) :: r(n)
       real(kind=4),intent(in)  :: a
       real(kind=4),intent(in)  :: beta
       !$omp  declare variant( vsrnglaplace:mkl_vsl_vsrnglaplace_omp_offload ) match( construct={target variant dispatch}, device={arch(gen)} )
     end function vsrnglaplace
  end interface

  interface
     integer function vdrnglaplace (method, stream, n, r, a, beta)
       import :: VSL_STREAM_STATE
       integer,intent(in)       :: method
       type(VSL_STREAM_STATE)   :: stream
       integer,intent(in)       :: n
       real(kind=8),intent(out) :: r(n)
       real(kind=8),intent(in)  :: a
       real(kind=8),intent(in)  :: beta
       !$omp  declare variant( vdrnglaplace:mkl_vsl_vdrnglaplace_omp_offload ) match( construct={target variant dispatch}, device={arch(gen)} )
     end function vdrnglaplace
  end interface

  ! Weibull distribution
  interface
     integer function vsrngweibull (method, stream, n, r, alpha, a, beta)
       import :: VSL_STREAM_STATE
       integer,intent(in)       :: method
       type(VSL_STREAM_STATE)   :: stream
       integer,intent(in)       :: n
       real(kind=4),intent(out) :: r(n)
       real(kind=4),intent(in)  :: alpha
       real(kind=4),intent(in)  :: a
       real(kind=4),intent(in)  :: beta
       !$omp  declare variant( vsrngweibull:mkl_vsl_vsrngweibull_omp_offload ) match( construct={target variant dispatch}, device={arch(gen)} )
     end function vsrngweibull
  end interface

  interface
     integer function vdrngweibull (method, stream, n, r, alpha, a, beta)
       import :: VSL_STREAM_STATE
       integer,intent(in)       :: method
       type(VSL_STREAM_STATE)   :: stream
       integer,intent(in)       :: n
       real(kind=8),intent(out) :: r(n)
       real(kind=8),intent(in)  :: alpha
       real(kind=8),intent(in)  :: a
       real(kind=8),intent(in)  :: beta
       !$omp  declare variant( vdrngweibull:mkl_vsl_vdrngweibull_omp_offload ) match( construct={target variant dispatch}, device={arch(gen)} )
     end function vdrngweibull
  end interface

  ! Cauchy distribution
  interface
     integer function vsrngcauchy (method, stream, n, r, a, beta)
       import :: VSL_STREAM_STATE
       integer,intent(in)       :: method
       type(VSL_STREAM_STATE)   :: stream
       integer,intent(in)       :: n
       real(kind=4),intent(out) :: r(n)
       real(kind=4),intent(in)  :: a
       real(kind=4),intent(in)  :: beta
       !$omp  declare variant( vsrngcauchy:mkl_vsl_vsrngcauchy_omp_offload ) match( construct={target variant dispatch}, device={arch(gen)} )
     end function vsrngcauchy
  end interface

  interface
     integer function vdrngcauchy (method, stream, n, r, a, beta)
       import :: VSL_STREAM_STATE
       integer,intent(in)       :: method
       type(VSL_STREAM_STATE)   :: stream
       integer,intent(in)       :: n
       real(kind=8),intent(out) :: r(n)
       real(kind=8),intent(in)  :: a
       real(kind=8),intent(in)  :: beta
       !$omp  declare variant( vdrngcauchy:mkl_vsl_vdrngcauchy_omp_offload ) match( construct={target variant dispatch}, device={arch(gen)} )
     end function vdrngcauchy
  end interface

  ! Rayleigh distribution
  interface
     integer function vsrngrayleigh (method, stream, n, r, a, beta)
       import :: VSL_STREAM_STATE
       integer,intent(in)       :: method
       type(VSL_STREAM_STATE)   :: stream
       integer,intent(in)       :: n
       real(kind=4),intent(out) :: r(n)
       real(kind=4),intent(in)  :: a
       real(kind=4),intent(in)  :: beta
       !$omp  declare variant( vsrngrayleigh:mkl_vsl_vsrngrayleigh_omp_offload ) match( construct={target variant dispatch}, device={arch(gen)} )
     end function vsrngrayleigh
  end interface

  interface
     integer function vdrngrayleigh (method, stream, n, r, a, beta)
       import :: VSL_STREAM_STATE
       integer,intent(in)       :: method
       type(VSL_STREAM_STATE)   :: stream
       integer,intent(in)       :: n
       real(kind=8),intent(out) :: r(n)
       real(kind=8),intent(in)  :: a
       real(kind=8),intent(in)  :: beta
       !$omp  declare variant( vdrngrayleigh:mkl_vsl_vdrngrayleigh_omp_offload ) match( construct={target variant dispatch}, device={arch(gen)} )
     end function vdrngrayleigh
  end interface

  ! Lognormal distribution
  interface
     integer function vsrnglognormal (method, stream, n, r, a, sigma, b, beta)
       import :: VSL_STREAM_STATE
       integer,intent(in)       :: method
       type(VSL_STREAM_STATE)   :: stream
       integer,intent(in)       :: n
       real(kind=4),intent(out) :: r(n)
       real(kind=4),intent(in)  :: a
       real(kind=4),intent(in)  :: sigma
       real(kind=4),intent(in)  :: b
       real(kind=4),intent(in)  :: beta
       !$omp  declare variant( vsrnglognormal:mkl_vsl_vsrnglognormal_omp_offload ) match( construct={target variant dispatch}, device={arch(gen)} )
     end function vsrnglognormal
  end interface

  interface
     integer function vdrnglognormal (method, stream, n, r, a, sigma, b, beta)
       import :: VSL_STREAM_STATE
       integer,intent(in)       :: method
       type(VSL_STREAM_STATE)   :: stream
       integer,intent(in)       :: n
       real(kind=8),intent(out) :: r(n)
       real(kind=8),intent(in)  :: a
       real(kind=8),intent(in)  :: sigma
       real(kind=8),intent(in)  :: b
       real(kind=8),intent(in)  :: beta
       !$omp  declare variant( vdrnglognormal:mkl_vsl_vdrnglognormal_omp_offload ) match( construct={target variant dispatch}, device={arch(gen)} )
     end function vdrnglognormal
  end interface

  ! Gumbel distribution
  interface
     integer function vsrnggumbel (method, stream, n, r, a, beta)
       import :: VSL_STREAM_STATE
       integer,intent(in)       :: method
       type(VSL_STREAM_STATE)   :: stream
       integer,intent(in)       :: n
       real(kind=4),intent(out) :: r(n)
       real(kind=4),intent(in)  :: a
       real(kind=4),intent(in)  :: beta
       !$omp  declare variant( vsrnggumbel:mkl_vsl_vsrnggumbel_omp_offload ) match( construct={target variant dispatch}, device={arch(gen)} )
     end function vsrnggumbel
  end interface

  interface
     integer function vdrnggumbel (method, stream, n, r, a, beta)
       import :: VSL_STREAM_STATE
       integer,intent(in)       :: method
       type(VSL_STREAM_STATE)   :: stream
       integer,intent(in)       :: n
       real(kind=8),intent(out) :: r(n)
       real(kind=8),intent(in)  :: a
       real(kind=8),intent(in)  :: beta
       !$omp  declare variant( vdrnggumbel:mkl_vsl_vdrnggumbel_omp_offload ) match( construct={target variant dispatch}, device={arch(gen)} )
     end function vdrnggumbel
  end interface

!++
!  VSL DISCRETE DISTRIBUTION GENERATOR FUNCTION INTERFACES.
!--
  ! Uniform distribution
  interface
     integer function virnguniform (method, stream, n, r, a, b)
       import :: VSL_STREAM_STATE
       integer,intent(in)          :: method
       type(VSL_STREAM_STATE)      :: stream
       integer,intent(in)          :: n
       integer(kind=4),intent(out) :: r(n)
       integer(kind=4),intent(in)  :: a
       integer(kind=4),intent(in)  :: b
       !$omp  declare variant( virnguniform:mkl_vsl_virnguniform_omp_offload ) match( construct={target variant dispatch}, device={arch(gen)} )
     end function virnguniform
  end interface

  ! UniformBits distribution
  interface
     integer function virnguniformbits (method, stream, n, r)
       import :: VSL_STREAM_STATE
       integer,intent(in)          :: method
       type(VSL_STREAM_STATE)      :: stream
       integer,intent(in)          :: n
       integer(kind=4),intent(out) :: r(n)
       !$omp  declare variant( virnguniformbits:mkl_vsl_virnguniformbits_omp_offload ) match( construct={target variant dispatch}, device={arch(gen)} )
     end function virnguniformbits
  end interface

  ! UniformBits32 distribution
  interface
     integer function virnguniformbits32 (method, stream, n, r)
       import :: VSL_STREAM_STATE
       integer,intent(in)          :: method
       type(VSL_STREAM_STATE)      :: stream
       integer,intent(in)          :: n
       integer(kind=4),intent(out) :: r(n)
       !$omp  declare variant( virnguniformbits32:mkl_vsl_virnguniformbits32_omp_offload ) match( construct={target variant dispatch}, device={arch(gen)} )
     end function virnguniformbits32
  end interface

  ! UniformBits64 distribution
  interface
     integer function virnguniformbits64 (method, stream, n, r)
       import :: VSL_STREAM_STATE
       integer,intent(in)          :: method
       type(VSL_STREAM_STATE)      :: stream
       integer,intent(in)          :: n
       integer(kind=8),intent(out) :: r(n)
       !$omp  declare variant( virnguniformbits64:mkl_vsl_virnguniformbits64_omp_offload ) match( construct={target variant dispatch}, device={arch(gen)} )
     end function virnguniformbits64
  end interface

  ! Bernoulli distribution
  interface
     integer function virngbernoulli (method, stream, n, r, p)
       import :: VSL_STREAM_STATE
       integer,intent(in)          :: method
       type(VSL_STREAM_STATE)      :: stream
       integer,intent(in)          :: n
       integer(kind=4),intent(out) :: r(n)
       real(kind=8),intent(in)     :: p
       !$omp  declare variant( virngbernoulli:mkl_vsl_virngbernoulli_omp_offload ) match( construct={target variant dispatch}, device={arch(gen)} )
     end function virngbernoulli
  end interface

  ! Geometric distribution
  interface
     integer function virnggeometric (method, stream, n, r, p)
       import :: VSL_STREAM_STATE
       integer,intent(in)          :: method
       type(VSL_STREAM_STATE)      :: stream
       integer,intent(in)          :: n
       integer(kind=4),intent(out) :: r(n)
       real(kind=8),intent(in)     :: p
       !$omp  declare variant( virnggeometric:mkl_vsl_virnggeometric_omp_offload ) match( construct={target variant dispatch}, device={arch(gen)} )
     end function virnggeometric
  end interface

  ! Poisson distribution
  interface
     integer function virngpoisson (method, stream, n, r, lambda)
       import :: VSL_STREAM_STATE
       integer,intent(in)          :: method
       type(VSL_STREAM_STATE)      :: stream
       integer,intent(in)          :: n
       integer(kind=4),intent(out) :: r(n)
       real(kind=8),intent(in)     :: lambda
       !$omp  declare variant( virngpoisson:mkl_vsl_virngpoisson_omp_offload ) match( construct={target variant dispatch}, device={arch(gen)} )
     end function virngpoisson
  end interface

  !++
  !  VSL SERVICE FUNCTION INTERFACES.
  !--

  ! NewStream - stream creation/initialization
        interface
          integer function vslnewstream( stream, brng, seed )
              import :: VSL_STREAM_STATE
            type(VSL_STREAM_STATE) :: stream
            integer,intent(in)     :: brng
            integer,intent(in)     :: seed
          end function
        end interface

  ! NewStreamEx - advanced stream creation/initialization
        interface
          integer function vslnewstreamex( stream, brng, n, params )
              import :: VSL_STREAM_STATE
            type(VSL_STREAM_STATE)     :: stream
            integer,intent(in)         :: brng
            integer,intent(in)         :: n
            integer(kind=4),intent(in) :: params(n)
          end function
        end interface

  ! DeleteStream - delete stream
        interface
          integer function vsldeletestream( stream )
              import :: VSL_STREAM_STATE
            type(VSL_STREAM_STATE) :: stream
          end function
        end interface

  ! CopyStream - copy all stream information
        interface
          integer function vslcopystream( newstream, srcstream )
              import :: VSL_STREAM_STATE
            type(VSL_STREAM_STATE) :: newstream
            type(VSL_STREAM_STATE) :: srcstream
          end function
        end interface

  ! CopyStreamState - copy stream state only
        interface
          integer function vslcopystreamstate( deststream, srcstream )
              import :: VSL_STREAM_STATE
            type(VSL_STREAM_STATE) :: deststream
            type(VSL_STREAM_STATE) :: srcstream
          end function
        end interface

  ! LeapfrogStream - leapfrog method
        interface
          integer function vslleapfrogstream( stream, k, nstreams )
              import :: VSL_STREAM_STATE
            type(VSL_STREAM_STATE) :: stream
            integer,intent(in)     :: k
            integer,intent(in)     :: nstreams
          end function
        end interface

  ! SkipAheadStream - skip-ahead method
        interface
          integer function vslskipaheadstream( stream, nskip )
              import :: VSL_STREAM_STATE
            type(VSL_STREAM_STATE)     :: stream
            integer(kind=8),intent(in) :: nskip
          end function
        end interface

  ! SkipAheadStreamEx - skip-ahead extended method
        interface
          integer function vslskipaheadstreamex( stream, n, params )
              import :: VSL_STREAM_STATE
            type(VSL_STREAM_STATE)     :: stream
            integer,intent(in)         :: n
            integer(kind=8),intent(in) :: params(n)
          end function
        end interface

  ! GetBrngProperties - get BRNG properties
        interface
          integer function vslgetbrngproperties( brng, properties )
              import :: VSL_BRNG_PROPERTIES
            integer(kind=4),intent(in)            :: brng
            type(VSL_BRNG_PROPERTIES),intent(out) :: properties
          end function
        end interface

  ! GetStreamStateBrng - get BRNG associated with given stream
        interface
          integer function vslgetstreamstatebrng( stream )
              import :: VSL_STREAM_STATE
            type(VSL_STREAM_STATE) :: stream
          end function
        end interface

  ! SaveStreamF - save stream to file
        interface
          integer function vslsavestreamf( stream, fname )
              import :: VSL_STREAM_STATE
            character(*)           :: fname
            type(VSL_STREAM_STATE) :: stream
          end function
        end interface

  ! LoadStreamF - save stream to file
        interface
          integer function vslloadstreamf( stream, fname )
              import :: VSL_STREAM_STATE
            character(*)           :: fname
            type(VSL_STREAM_STATE) :: stream
          end function
        end interface

  ! SaveStreamM - save stream to memory
        interface
          integer function vslsavestreamm( stream, memptr )
              import :: VSL_STREAM_STATE
            integer(kind=1),dimension(*),intent(out)::memptr
            type(VSL_STREAM_STATE),intent(in)       :: stream
          end function
        end interface

  ! LoadStreamM - load stream from memory
        interface
          integer function vslloadstreamm( stream, memptr )
              import :: VSL_STREAM_STATE
            integer(kind=1),dimension(*),intent(in)::memptr
            type(VSL_STREAM_STATE),intent(out)      ::stream
          end function
        end interface

  ! GetStreamSize - get size of random stream
        interface
          integer function vslgetstreamsize( stream )
              import :: VSL_STREAM_STATE
            type(VSL_STREAM_STATE),intent(in) :: stream
          end function
        end interface

end module onemkl_vsl_omp_offload


module onemkl_vsl_omp_offload_ilp64
  use onemkl_vsl_omp_offload
end module onemkl_vsl_omp_offload_ilp64


module onemkl_vsl_omp_offload_lp64
  use onemkl_vsl_omp_variant_lp64

  !++
  !  VSL CONTINUOUS DISTRIBUTION GENERATOR FUNCTION INTERFACES.
  !--

  ! Uniform distribution
  interface
     integer function vsrnguniform(method, stream, n, r, a, b)
       import :: VSL_STREAM_STATE
       integer,intent(in)       :: method
       type(VSL_STREAM_STATE)   :: stream
       integer,intent(in)       :: n
       real(kind=4),intent(out) :: r(n)
       real(kind=4),intent(in)  :: a
       real(kind=4),intent(in)  :: b
      !$omp  declare variant( vsrnguniform:mkl_vsl_vsrnguniform_omp_offload_lp64 ) match( construct={target variant dispatch}, device={arch(gen)} )
     end function vsrnguniform
  end interface

  interface
     integer function vdrnguniform (method, stream, n, r, a, b)
       import :: VSL_STREAM_STATE
       integer,intent(in)       :: method
       type(VSL_STREAM_STATE)   :: stream
       integer,intent(in)       :: n
       real(kind=8),intent(out) :: r(n)
       real(kind=8),intent(in)  :: a
       real(kind=8),intent(in)  :: b
       !$omp  declare variant( vdrnguniform:mkl_vsl_vdrnguniform_omp_offload_lp64 ) match( construct={target variant dispatch}, device={arch(gen)} )
     end function vdrnguniform
  end interface
  
  ! Gaussian distribution
  interface
     integer function vsrnggaussian (method, stream, n, r, a, sigma)
       import :: VSL_STREAM_STATE
       integer,intent(in)       :: method
       type(VSL_STREAM_STATE)   :: stream
       integer,intent(in)       :: n
       real(kind=4),intent(out) :: r(n)
       real(kind=4),intent(in)  :: a
       real(kind=4),intent(in)  :: sigma
       !$omp  declare variant( vsrnggaussian:mkl_vsl_vsrnggaussian_omp_offload_lp64 ) match( construct={target variant dispatch}, device={arch(gen)} )
     end function vsrnggaussian
  end interface

  interface
     integer function vdrnggaussian (method, stream, n, r, a, sigma)
       import :: VSL_STREAM_STATE
       integer,intent(in)       :: method
       type(VSL_STREAM_STATE)   :: stream
       integer,intent(in)       :: n
       real(kind=8),intent(out) :: r(n)
       real(kind=8),intent(in)  :: a
       real(kind=8),intent(in)  :: sigma
       !$omp  declare variant( vdrnggaussian:mkl_vsl_vdrnggaussian_omp_offload_lp64 ) match( construct={target variant dispatch}, device={arch(gen)} )
     end function vdrnggaussian
  end interface

  ! GaussianMV distribution
  interface
     integer function vsrnggaussianmv (method, stream, n, r, dimen, mstorage, a, t)
       import :: VSL_STREAM_STATE
       integer,intent(in)       :: method
       type(VSL_STREAM_STATE)   :: stream
       integer,intent(in)       :: n
       integer, intent(in)      :: dimen
       real(kind=4),intent(out) :: r(dimen,n)
       integer, intent(in)      :: mstorage
       real(kind=4),intent(in)  :: a(dimen)
       real(kind=4),intent(in)  :: t(dimen,dimen)
       !$omp  declare variant( vsrnggaussianmv:mkl_vsl_vsrnggaussianmv_omp_offload_lp64 ) match( construct={target variant dispatch}, device={arch(gen)} )
     end function vsrnggaussianmv
  end interface

  interface
     integer function vdrnggaussianmv (method, stream, n, r, dimen, mstorage, a, t)
       import :: VSL_STREAM_STATE
       integer,intent(in)       :: method
       type(VSL_STREAM_STATE)   :: stream
       integer,intent(in)       :: n
       integer, intent(in)      :: dimen
       real(kind=8),intent(out) :: r(dimen,n)
       integer, intent(in)      :: mstorage
       real(kind=8),intent(in)  :: a(dimen)
       real(kind=8),intent(in)  :: t(dimen,dimen)
       !$omp  declare variant( vdrnggaussianmv:mkl_vsl_vdrnggaussianmv_omp_offload_lp64 ) match( construct={target variant dispatch}, device={arch(gen)} )
     end function vdrnggaussianmv
  end interface

  ! Exponential distribution
  interface
     integer function vsrngexponential (method, stream, n, r, a, beta)
       import :: VSL_STREAM_STATE
       integer,intent(in)       :: method
       type(VSL_STREAM_STATE)   :: stream
       integer,intent(in)       :: n
       real(kind=4),intent(out) :: r(n)
       real(kind=4),intent(in)  :: a
       real(kind=4),intent(in)  :: beta
       !$omp  declare variant( vsrngexponential:mkl_vsl_vsrngexponential_omp_offload_lp64 ) match( construct={target variant dispatch}, device={arch(gen)} )
     end function vsrngexponential
  end interface

  interface
     integer function vdrngexponential (method, stream, n, r, a, beta)
       import :: VSL_STREAM_STATE
       integer,intent(in)       :: method
       type(VSL_STREAM_STATE)   :: stream
       integer,intent(in)       :: n
       real(kind=8),intent(out) :: r(n)
       real(kind=8),intent(in)  :: a
       real(kind=8),intent(in)  :: beta
       !$omp  declare variant( vdrngexponential:mkl_vsl_vdrngexponential_omp_offload_lp64 ) match( construct={target variant dispatch}, device={arch(gen)} )
     end function vdrngexponential
  end interface

  ! Laplace distribution
  interface
     integer function vsrnglaplace (method, stream, n, r, a, beta)
       import :: VSL_STREAM_STATE
       integer,intent(in)       :: method
       type(VSL_STREAM_STATE)   :: stream
       integer,intent(in)       :: n
       real(kind=4),intent(out) :: r(n)
       real(kind=4),intent(in)  :: a
       real(kind=4),intent(in)  :: beta
       !$omp  declare variant( vsrnglaplace:mkl_vsl_vsrnglaplace_omp_offload_lp64 ) match( construct={target variant dispatch}, device={arch(gen)} )
     end function vsrnglaplace
  end interface

  interface
     integer function vdrnglaplace (method, stream, n, r, a, beta)
       import :: VSL_STREAM_STATE
       integer,intent(in)       :: method
       type(VSL_STREAM_STATE)   :: stream
       integer,intent(in)       :: n
       real(kind=8),intent(out) :: r(n)
       real(kind=8),intent(in)  :: a
       real(kind=8),intent(in)  :: beta
       !$omp  declare variant( vdrnglaplace:mkl_vsl_vdrnglaplace_omp_offload_lp64 ) match( construct={target variant dispatch}, device={arch(gen)} )
     end function vdrnglaplace
  end interface

  ! Weibull distribution
  interface
     integer function vsrngweibull (method, stream, n, r, alpha, a, beta)
       import :: VSL_STREAM_STATE
       integer,intent(in)       :: method
       type(VSL_STREAM_STATE)   :: stream
       integer,intent(in)       :: n
       real(kind=4),intent(out) :: r(n)
       real(kind=4),intent(in)  :: alpha
       real(kind=4),intent(in)  :: a
       real(kind=4),intent(in)  :: beta
       !$omp  declare variant( vsrngweibull:mkl_vsl_vsrngweibull_omp_offload_lp64 ) match( construct={target variant dispatch}, device={arch(gen)} )
     end function vsrngweibull
  end interface

  interface
     integer function vdrngweibull (method, stream, n, r, alpha, a, beta)
       import :: VSL_STREAM_STATE
       integer,intent(in)       :: method
       type(VSL_STREAM_STATE)   :: stream
       integer,intent(in)       :: n
       real(kind=8),intent(out) :: r(n)
       real(kind=8),intent(in)  :: alpha
       real(kind=8),intent(in)  :: a
       real(kind=8),intent(in)  :: beta
       !$omp  declare variant( vdrngweibull:mkl_vsl_vdrngweibull_omp_offload_lp64 ) match( construct={target variant dispatch}, device={arch(gen)} )
     end function vdrngweibull
  end interface

  ! Cauchy distribution
  interface
     integer function vsrngcauchy (method, stream, n, r, a, beta)
       import :: VSL_STREAM_STATE
       integer,intent(in)       :: method
       type(VSL_STREAM_STATE)   :: stream
       integer,intent(in)       :: n
       real(kind=4),intent(out) :: r(n)
       real(kind=4),intent(in)  :: a
       real(kind=4),intent(in)  :: beta
       !$omp  declare variant( vsrngcauchy:mkl_vsl_vsrngcauchy_omp_offload_lp64 ) match( construct={target variant dispatch}, device={arch(gen)} )
     end function vsrngcauchy
  end interface

  interface
     integer function vdrngcauchy (method, stream, n, r, a, beta)
       import :: VSL_STREAM_STATE
       integer,intent(in)       :: method
       type(VSL_STREAM_STATE)   :: stream
       integer,intent(in)       :: n
       real(kind=8),intent(out) :: r(n)
       real(kind=8),intent(in)  :: a
       real(kind=8),intent(in)  :: beta
       !$omp  declare variant( vdrngcauchy:mkl_vsl_vdrngcauchy_omp_offload_lp64 ) match( construct={target variant dispatch}, device={arch(gen)} )
     end function vdrngcauchy
  end interface

  ! Rayleigh distribution
  interface
     integer function vsrngrayleigh (method, stream, n, r, a, beta)
       import :: VSL_STREAM_STATE
       integer,intent(in)       :: method
       type(VSL_STREAM_STATE)   :: stream
       integer,intent(in)       :: n
       real(kind=4),intent(out) :: r(n)
       real(kind=4),intent(in)  :: a
       real(kind=4),intent(in)  :: beta
       !$omp  declare variant( vsrngrayleigh:mkl_vsl_vsrngrayleigh_omp_offload_lp64 ) match( construct={target variant dispatch}, device={arch(gen)} )
     end function vsrngrayleigh
  end interface

  interface
     integer function vdrngrayleigh (method, stream, n, r, a, beta)
       import :: VSL_STREAM_STATE
       integer,intent(in)       :: method
       type(VSL_STREAM_STATE)   :: stream
       integer,intent(in)       :: n
       real(kind=8),intent(out) :: r(n)
       real(kind=8),intent(in)  :: a
       real(kind=8),intent(in)  :: beta
       !$omp  declare variant( vdrngrayleigh:mkl_vsl_vdrngrayleigh_omp_offload_lp64 ) match( construct={target variant dispatch}, device={arch(gen)} )
     end function vdrngrayleigh
  end interface

  ! Lognormal distribution
  interface
     integer function vsrnglognormal (method, stream, n, r, a, sigma, b, beta)
       import :: VSL_STREAM_STATE
       integer,intent(in)       :: method
       type(VSL_STREAM_STATE)   :: stream
       integer,intent(in)       :: n
       real(kind=4),intent(out) :: r(n)
       real(kind=4),intent(in)  :: a
       real(kind=4),intent(in)  :: sigma
       real(kind=4),intent(in)  :: b
       real(kind=4),intent(in)  :: beta
       !$omp  declare variant( vsrnglognormal:mkl_vsl_vsrnglognormal_omp_offload_lp64 ) match( construct={target variant dispatch}, device={arch(gen)} )
     end function vsrnglognormal
  end interface

  interface
     integer function vdrnglognormal (method, stream, n, r, a, sigma, b, beta)
       import :: VSL_STREAM_STATE
       integer,intent(in)       :: method
       type(VSL_STREAM_STATE)   :: stream
       integer,intent(in)       :: n
       real(kind=8),intent(out) :: r(n)
       real(kind=8),intent(in)  :: a
       real(kind=8),intent(in)  :: sigma
       real(kind=8),intent(in)  :: b
       real(kind=8),intent(in)  :: beta
       !$omp  declare variant( vdrnglognormal:mkl_vsl_vdrnglognormal_omp_offload_lp64 ) match( construct={target variant dispatch}, device={arch(gen)} )
     end function vdrnglognormal
  end interface

  ! Gumbel distribution
  interface
     integer function vsrnggumbel (method, stream, n, r, a, beta)
       import :: VSL_STREAM_STATE
       integer,intent(in)       :: method
       type(VSL_STREAM_STATE)   :: stream
       integer,intent(in)       :: n
       real(kind=4),intent(out) :: r(n)
       real(kind=4),intent(in)  :: a
       real(kind=4),intent(in)  :: beta
       !$omp  declare variant( vsrnggumbel:mkl_vsl_vsrnggumbel_omp_offload_lp64 ) match( construct={target variant dispatch}, device={arch(gen)} )
     end function vsrnggumbel
  end interface

  interface
     integer function vdrnggumbel (method, stream, n, r, a, beta)
       import :: VSL_STREAM_STATE
       integer,intent(in)       :: method
       type(VSL_STREAM_STATE)   :: stream
       integer,intent(in)       :: n
       real(kind=8),intent(out) :: r(n)
       real(kind=8),intent(in)  :: a
       real(kind=8),intent(in)  :: beta
       !$omp  declare variant( vdrnggumbel:mkl_vsl_vdrnggumbel_omp_offload_lp64 ) match( construct={target variant dispatch}, device={arch(gen)} )
     end function vdrnggumbel
  end interface

!++
!  VSL DISCRETE DISTRIBUTION GENERATOR FUNCTION INTERFACES.
!--
  ! Uniform distribution
  interface
     integer function virnguniform (method, stream, n, r, a, b)
       import :: VSL_STREAM_STATE
       integer,intent(in)          :: method
       type(VSL_STREAM_STATE)      :: stream
       integer,intent(in)          :: n
       integer(kind=4),intent(out) :: r(n)
       integer(kind=4),intent(in)  :: a
       integer(kind=4),intent(in)  :: b
       !$omp  declare variant( virnguniform:mkl_vsl_virnguniform_omp_offload_lp64 ) match( construct={target variant dispatch}, device={arch(gen)} )
     end function virnguniform
  end interface

  ! UniformBits distribution
  interface
     integer function virnguniformbits (method, stream, n, r)
       import :: VSL_STREAM_STATE
       integer,intent(in)          :: method
       type(VSL_STREAM_STATE)      :: stream
       integer,intent(in)          :: n
       integer(kind=4),intent(out) :: r(n)
       !$omp  declare variant( virnguniformbits:mkl_vsl_virnguniformbits_omp_offload_lp64 ) match( construct={target variant dispatch}, device={arch(gen)} )
     end function virnguniformbits
  end interface

  ! UniformBits32 distribution
  interface
     integer function virnguniformbits32 (method, stream, n, r)
       import :: VSL_STREAM_STATE
       integer,intent(in)          :: method
       type(VSL_STREAM_STATE)      :: stream
       integer,intent(in)          :: n
       integer(kind=4),intent(out) :: r(n)
       !$omp  declare variant( virnguniformbits32:mkl_vsl_virnguniformbits32_omp_offload_lp64 ) match( construct={target variant dispatch}, device={arch(gen)} )
     end function virnguniformbits32
  end interface

  ! UniformBits64 distribution
  interface
     integer function virnguniformbits64 (method, stream, n, r)
       import :: VSL_STREAM_STATE
       integer,intent(in)          :: method
       type(VSL_STREAM_STATE)      :: stream
       integer,intent(in)          :: n
       integer(kind=8),intent(out) :: r(n)
       !$omp  declare variant( virnguniformbits64:mkl_vsl_virnguniformbits64_omp_offload_lp64 ) match( construct={target variant dispatch}, device={arch(gen)} )
     end function virnguniformbits64
  end interface

  ! Bernoulli distribution
  interface
     integer function virngbernoulli (method, stream, n, r, p)
       import :: VSL_STREAM_STATE
       integer,intent(in)          :: method
       type(VSL_STREAM_STATE)      :: stream
       integer,intent(in)          :: n
       integer(kind=4),intent(out) :: r(n)
       real(kind=8),intent(in)     :: p
       !$omp  declare variant( virngbernoulli:mkl_vsl_virngbernoulli_omp_offload_lp64 ) match( construct={target variant dispatch}, device={arch(gen)} )
     end function virngbernoulli
  end interface

  ! Geometric distribution
  interface
     integer function virnggeometric (method, stream, n, r, p)
       import :: VSL_STREAM_STATE
       integer,intent(in)          :: method
       type(VSL_STREAM_STATE)      :: stream
       integer,intent(in)          :: n
       integer(kind=4),intent(out) :: r(n)
       real(kind=8),intent(in)     :: p
       !$omp  declare variant( virnggeometric:mkl_vsl_virnggeometric_omp_offload_lp64 ) match( construct={target variant dispatch}, device={arch(gen)} )
     end function virnggeometric
  end interface

  ! Poisson distribution
  interface
     integer function virngpoisson (method, stream, n, r, lambda)
       import :: VSL_STREAM_STATE
       integer,intent(in)          :: method
       type(VSL_STREAM_STATE)      :: stream
       integer,intent(in)          :: n
       integer(kind=4),intent(out) :: r(n)
       real(kind=8),intent(in)     :: lambda
       !$omp  declare variant( virngpoisson:mkl_vsl_virngpoisson_omp_offload_lp64 ) match( construct={target variant dispatch}, device={arch(gen)} )
     end function virngpoisson
  end interface

  !++
  !  VSL SERVICE FUNCTION INTERFACES.
  !--

  ! NewStream - stream creation/initialization
        interface
          integer function vslnewstream( stream, brng, seed )
              import :: VSL_STREAM_STATE
            type(VSL_STREAM_STATE) :: stream
            integer,intent(in)     :: brng
            integer,intent(in)     :: seed
          end function
        end interface

  ! NewStreamEx - advanced stream creation/initialization
        interface
          integer function vslnewstreamex( stream, brng, n, params )
              import :: VSL_STREAM_STATE
            type(VSL_STREAM_STATE)     :: stream
            integer,intent(in)         :: brng
            integer,intent(in)         :: n
            integer(kind=4),intent(in) :: params(n)
          end function
        end interface

  ! DeleteStream - delete stream
        interface
          integer function vsldeletestream( stream )
              import :: VSL_STREAM_STATE
            type(VSL_STREAM_STATE) :: stream
          end function
        end interface

  ! CopyStream - copy all stream information
        interface
          integer function vslcopystream( newstream, srcstream )
              import :: VSL_STREAM_STATE
            type(VSL_STREAM_STATE) :: newstream
            type(VSL_STREAM_STATE) :: srcstream
          end function
        end interface

  ! CopyStreamState - copy stream state only
        interface
          integer function vslcopystreamstate( deststream, srcstream )
              import :: VSL_STREAM_STATE
            type(VSL_STREAM_STATE) :: deststream
            type(VSL_STREAM_STATE) :: srcstream
          end function
        end interface

  ! LeapfrogStream - leapfrog method
        interface
          integer function vslleapfrogstream( stream, k, nstreams )
              import :: VSL_STREAM_STATE
            type(VSL_STREAM_STATE) :: stream
            integer,intent(in)     :: k
            integer,intent(in)     :: nstreams
          end function
        end interface

  ! SkipAheadStream - skip-ahead method
        interface
          integer function vslskipaheadstream( stream, nskip )
              import :: VSL_STREAM_STATE
            type(VSL_STREAM_STATE)     :: stream
            integer(kind=8),intent(in) :: nskip
          end function
        end interface

  ! SkipAheadStreamEx - skip-ahead extended method
        interface
          integer function vslskipaheadstreamex( stream, n, params )
              import :: VSL_STREAM_STATE
            type(VSL_STREAM_STATE)     :: stream
            integer,intent(in)         :: n
            integer(kind=8),intent(in) :: params(n)
          end function
        end interface

  ! GetBrngProperties - get BRNG properties
        interface
          integer function vslgetbrngproperties( brng, properties )
              import :: VSL_BRNG_PROPERTIES
            integer(kind=4),intent(in)            :: brng
            type(VSL_BRNG_PROPERTIES),intent(out) :: properties
          end function
        end interface

  ! GetStreamStateBrng - get BRNG associated with given stream
        interface
          integer function vslgetstreamstatebrng( stream )
              import :: VSL_STREAM_STATE
            type(VSL_STREAM_STATE) :: stream
          end function
        end interface

  ! SaveStreamF - save stream to file
        interface
          integer function vslsavestreamf( stream, fname )
              import :: VSL_STREAM_STATE
            character(*)           :: fname
            type(VSL_STREAM_STATE) :: stream
          end function
        end interface

  ! LoadStreamF - save stream to file
        interface
          integer function vslloadstreamf( stream, fname )
              import :: VSL_STREAM_STATE
            character(*)           :: fname
            type(VSL_STREAM_STATE) :: stream
          end function
        end interface

  ! SaveStreamM - save stream to memory
        interface
          integer function vslsavestreamm( stream, memptr )
              import :: VSL_STREAM_STATE
            integer(kind=1),dimension(*),intent(out)::memptr
            type(VSL_STREAM_STATE),intent(in)       :: stream
          end function
        end interface

  ! LoadStreamM - load stream from memory
        interface
          integer function vslloadstreamm( stream, memptr )
              import :: VSL_STREAM_STATE
            integer(kind=1),dimension(*),intent(in)::memptr
            type(VSL_STREAM_STATE),intent(out)      ::stream
          end function
        end interface

  ! GetStreamSize - get size of random stream
        interface
          integer function vslgetstreamsize( stream )
              import :: VSL_STREAM_STATE
            type(VSL_STREAM_STATE),intent(in) :: stream
          end function
        end interface

end module onemkl_vsl_omp_offload_lp64
