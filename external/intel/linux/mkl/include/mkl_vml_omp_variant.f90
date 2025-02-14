! -== file: mkl_vml_omp_variant.f90 ==-
!===============================================================================
! Copyright 2006-2021 Intel Corporation.
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
!      Intel(R) oneAPI Math Kernel Library (oneMKL) FORTRAN interface for
!      OpenMP offload for VM
!*******************************************************************************

    interface

        integer(kind=4) function mkl_vm_vmlsetmode_omp_offload_ilp64(n)
        integer(kind=8), intent(in) :: n
        end function mkl_vm_vmlsetmode_omp_offload_ilp64

        integer(kind=4) function mkl_vm_vmlgetmode_omp_offload_ilp64()
        end function mkl_vm_vmlgetmode_omp_offload_ilp64

        integer(kind=4) function mkl_vm_vmlseterrstatus_omp_offload_ilp64(n)
        integer(kind=8), intent(in) :: n
        end function mkl_vm_vmlseterrstatus_omp_offload_ilp64

        integer(kind=4) function mkl_vm_vmlgeterrstatus_omp_offload_ilp64()
        end function mkl_vm_vmlgeterrstatus_omp_offload_ilp64

        integer(kind=4) function mkl_vm_vmlclearerrstatus_omp_offload_ilp64()
        end function mkl_vm_vmlclearerrstatus_omp_offload_ilp64





! function: Abs, indexing: simple
        subroutine mkl_vm_vsabs_omp_offload_ilp64(n, a, y)
        integer(kind=8), intent(in) :: n
        real(kind=4), intent(in) :: a(*)
        real(kind=4), intent(inout) :: y(*)
        end subroutine mkl_vm_vsabs_omp_offload_ilp64

        subroutine mkl_vm_vmsabs_omp_offload_ilp64(n, a, y, mode)
        integer(kind=8), intent(in) :: n
        real(kind=4), intent(in) :: a(*)
        real(kind=4), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmsabs_omp_offload_ilp64

        subroutine mkl_vm_vdabs_omp_offload_ilp64(n, a, y)
        integer(kind=8), intent(in) :: n
        real(kind=8), intent(in) :: a(*)
        real(kind=8), intent(inout) :: y(*)
        end subroutine mkl_vm_vdabs_omp_offload_ilp64

        subroutine mkl_vm_vmdabs_omp_offload_ilp64(n, a, y, mode)
        integer(kind=8), intent(in) :: n
        real(kind=8), intent(in) :: a(*)
        real(kind=8), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmdabs_omp_offload_ilp64

        subroutine mkl_vm_vcabs_omp_offload_ilp64(n, a, y)
        integer(kind=8), intent(in) :: n
        complex(kind=4), intent(in) :: a(*)
        real(kind=4), intent(inout) :: y(*)
        end subroutine mkl_vm_vcabs_omp_offload_ilp64

        subroutine mkl_vm_vmcabs_omp_offload_ilp64(n, a, y, mode)
        integer(kind=8), intent(in) :: n
        complex(kind=4), intent(in) :: a(*)
        real(kind=4), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmcabs_omp_offload_ilp64

        subroutine mkl_vm_vzabs_omp_offload_ilp64(n, a, y)
        integer(kind=8), intent(in) :: n
        complex(kind=8), intent(in) :: a(*)
        real(kind=8), intent(inout) :: y(*)
        end subroutine mkl_vm_vzabs_omp_offload_ilp64

        subroutine mkl_vm_vmzabs_omp_offload_ilp64(n, a, y, mode)
        integer(kind=8), intent(in) :: n
        complex(kind=8), intent(in) :: a(*)
        real(kind=8), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmzabs_omp_offload_ilp64


! function: Abs, indexing: strided
        subroutine mkl_vm_vsabsi_omp_offload_ilp64(n, a, inca, y, incy)
        integer(kind=8), intent(in) :: n
        real(kind=4), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        real(kind=4), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        end subroutine mkl_vm_vsabsi_omp_offload_ilp64

        subroutine mkl_vm_vmsabsi_omp_offload_ilp64(n, a, inca, y, incy, mode)
        integer(kind=8), intent(in) :: n
        real(kind=4), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        real(kind=4), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmsabsi_omp_offload_ilp64

        subroutine mkl_vm_vdabsi_omp_offload_ilp64(n, a, inca, y, incy)
        integer(kind=8), intent(in) :: n
        real(kind=8), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        real(kind=8), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        end subroutine mkl_vm_vdabsi_omp_offload_ilp64

        subroutine mkl_vm_vmdabsi_omp_offload_ilp64(n, a, inca, y, incy, mode)
        integer(kind=8), intent(in) :: n
        real(kind=8), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        real(kind=8), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmdabsi_omp_offload_ilp64

        subroutine mkl_vm_vcabsi_omp_offload_ilp64(n, a, inca, y, incy)
        integer(kind=8), intent(in) :: n
        complex(kind=4), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        real(kind=4), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        end subroutine mkl_vm_vcabsi_omp_offload_ilp64

        subroutine mkl_vm_vmcabsi_omp_offload_ilp64(n, a, inca, y, incy, mode)
        integer(kind=8), intent(in) :: n
        complex(kind=4), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        real(kind=4), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmcabsi_omp_offload_ilp64

        subroutine mkl_vm_vzabsi_omp_offload_ilp64(n, a, inca, y, incy)
        integer(kind=8), intent(in) :: n
        complex(kind=8), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        real(kind=8), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        end subroutine mkl_vm_vzabsi_omp_offload_ilp64

        subroutine mkl_vm_vmzabsi_omp_offload_ilp64(n, a, inca, y, incy, mode)
        integer(kind=8), intent(in) :: n
        complex(kind=8), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        real(kind=8), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmzabsi_omp_offload_ilp64




! function: Acos, indexing: simple
        subroutine mkl_vm_vsacos_omp_offload_ilp64(n, a, y)
        integer(kind=8), intent(in) :: n
        real(kind=4), intent(in) :: a(*)
        real(kind=4), intent(inout) :: y(*)
        end subroutine mkl_vm_vsacos_omp_offload_ilp64

        subroutine mkl_vm_vmsacos_omp_offload_ilp64(n, a, y, mode)
        integer(kind=8), intent(in) :: n
        real(kind=4), intent(in) :: a(*)
        real(kind=4), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmsacos_omp_offload_ilp64

        subroutine mkl_vm_vdacos_omp_offload_ilp64(n, a, y)
        integer(kind=8), intent(in) :: n
        real(kind=8), intent(in) :: a(*)
        real(kind=8), intent(inout) :: y(*)
        end subroutine mkl_vm_vdacos_omp_offload_ilp64

        subroutine mkl_vm_vmdacos_omp_offload_ilp64(n, a, y, mode)
        integer(kind=8), intent(in) :: n
        real(kind=8), intent(in) :: a(*)
        real(kind=8), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmdacos_omp_offload_ilp64

        subroutine mkl_vm_vcacos_omp_offload_ilp64(n, a, y)
        integer(kind=8), intent(in) :: n
        complex(kind=4), intent(in) :: a(*)
        complex(kind=4), intent(inout) :: y(*)
        end subroutine mkl_vm_vcacos_omp_offload_ilp64

        subroutine mkl_vm_vmcacos_omp_offload_ilp64(n, a, y, mode)
        integer(kind=8), intent(in) :: n
        complex(kind=4), intent(in) :: a(*)
        complex(kind=4), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmcacos_omp_offload_ilp64

        subroutine mkl_vm_vzacos_omp_offload_ilp64(n, a, y)
        integer(kind=8), intent(in) :: n
        complex(kind=8), intent(in) :: a(*)
        complex(kind=8), intent(inout) :: y(*)
        end subroutine mkl_vm_vzacos_omp_offload_ilp64

        subroutine mkl_vm_vmzacos_omp_offload_ilp64(n, a, y, mode)
        integer(kind=8), intent(in) :: n
        complex(kind=8), intent(in) :: a(*)
        complex(kind=8), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmzacos_omp_offload_ilp64


! function: Acos, indexing: strided
        subroutine mkl_vm_vsacosi_omp_offload_ilp64(n, a, inca, y, incy)
        integer(kind=8), intent(in) :: n
        real(kind=4), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        real(kind=4), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        end subroutine mkl_vm_vsacosi_omp_offload_ilp64

        subroutine mkl_vm_vmsacosi_omp_offload_ilp64(n, a, inca, y, incy, mode)
        integer(kind=8), intent(in) :: n
        real(kind=4), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        real(kind=4), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmsacosi_omp_offload_ilp64

        subroutine mkl_vm_vdacosi_omp_offload_ilp64(n, a, inca, y, incy)
        integer(kind=8), intent(in) :: n
        real(kind=8), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        real(kind=8), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        end subroutine mkl_vm_vdacosi_omp_offload_ilp64

        subroutine mkl_vm_vmdacosi_omp_offload_ilp64(n, a, inca, y, incy, mode)
        integer(kind=8), intent(in) :: n
        real(kind=8), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        real(kind=8), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmdacosi_omp_offload_ilp64

        subroutine mkl_vm_vcacosi_omp_offload_ilp64(n, a, inca, y, incy)
        integer(kind=8), intent(in) :: n
        complex(kind=4), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        complex(kind=4), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        end subroutine mkl_vm_vcacosi_omp_offload_ilp64

        subroutine mkl_vm_vmcacosi_omp_offload_ilp64(n, a, inca, y, incy, mode)
        integer(kind=8), intent(in) :: n
        complex(kind=4), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        complex(kind=4), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmcacosi_omp_offload_ilp64

        subroutine mkl_vm_vzacosi_omp_offload_ilp64(n, a, inca, y, incy)
        integer(kind=8), intent(in) :: n
        complex(kind=8), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        complex(kind=8), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        end subroutine mkl_vm_vzacosi_omp_offload_ilp64

        subroutine mkl_vm_vmzacosi_omp_offload_ilp64(n, a, inca, y, incy, mode)
        integer(kind=8), intent(in) :: n
        complex(kind=8), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        complex(kind=8), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmzacosi_omp_offload_ilp64




! function: Acosh, indexing: simple
        subroutine mkl_vm_vsacosh_omp_offload_ilp64(n, a, y)
        integer(kind=8), intent(in) :: n
        real(kind=4), intent(in) :: a(*)
        real(kind=4), intent(inout) :: y(*)
        end subroutine mkl_vm_vsacosh_omp_offload_ilp64

        subroutine mkl_vm_vmsacosh_omp_offload_ilp64(n, a, y, mode)
        integer(kind=8), intent(in) :: n
        real(kind=4), intent(in) :: a(*)
        real(kind=4), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmsacosh_omp_offload_ilp64

        subroutine mkl_vm_vdacosh_omp_offload_ilp64(n, a, y)
        integer(kind=8), intent(in) :: n
        real(kind=8), intent(in) :: a(*)
        real(kind=8), intent(inout) :: y(*)
        end subroutine mkl_vm_vdacosh_omp_offload_ilp64

        subroutine mkl_vm_vmdacosh_omp_offload_ilp64(n, a, y, mode)
        integer(kind=8), intent(in) :: n
        real(kind=8), intent(in) :: a(*)
        real(kind=8), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmdacosh_omp_offload_ilp64

        subroutine mkl_vm_vcacosh_omp_offload_ilp64(n, a, y)
        integer(kind=8), intent(in) :: n
        complex(kind=4), intent(in) :: a(*)
        complex(kind=4), intent(inout) :: y(*)
        end subroutine mkl_vm_vcacosh_omp_offload_ilp64

        subroutine mkl_vm_vmcacosh_omp_offload_ilp64(n, a, y, mode)
        integer(kind=8), intent(in) :: n
        complex(kind=4), intent(in) :: a(*)
        complex(kind=4), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmcacosh_omp_offload_ilp64

        subroutine mkl_vm_vzacosh_omp_offload_ilp64(n, a, y)
        integer(kind=8), intent(in) :: n
        complex(kind=8), intent(in) :: a(*)
        complex(kind=8), intent(inout) :: y(*)
        end subroutine mkl_vm_vzacosh_omp_offload_ilp64

        subroutine mkl_vm_vmzacosh_omp_offload_ilp64(n, a, y, mode)
        integer(kind=8), intent(in) :: n
        complex(kind=8), intent(in) :: a(*)
        complex(kind=8), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmzacosh_omp_offload_ilp64


! function: Acosh, indexing: strided
        subroutine mkl_vm_vsacoshi_omp_offload_ilp64(n, a, inca, y, incy)
        integer(kind=8), intent(in) :: n
        real(kind=4), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        real(kind=4), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        end subroutine mkl_vm_vsacoshi_omp_offload_ilp64

        subroutine mkl_vm_vmsacoshi_omp_offload_ilp64(n, a, inca, y, incy, mode)
        integer(kind=8), intent(in) :: n
        real(kind=4), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        real(kind=4), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmsacoshi_omp_offload_ilp64

        subroutine mkl_vm_vdacoshi_omp_offload_ilp64(n, a, inca, y, incy)
        integer(kind=8), intent(in) :: n
        real(kind=8), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        real(kind=8), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        end subroutine mkl_vm_vdacoshi_omp_offload_ilp64

        subroutine mkl_vm_vmdacoshi_omp_offload_ilp64(n, a, inca, y, incy, mode)
        integer(kind=8), intent(in) :: n
        real(kind=8), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        real(kind=8), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmdacoshi_omp_offload_ilp64

        subroutine mkl_vm_vcacoshi_omp_offload_ilp64(n, a, inca, y, incy)
        integer(kind=8), intent(in) :: n
        complex(kind=4), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        complex(kind=4), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        end subroutine mkl_vm_vcacoshi_omp_offload_ilp64

        subroutine mkl_vm_vmcacoshi_omp_offload_ilp64(n, a, inca, y, incy, mode)
        integer(kind=8), intent(in) :: n
        complex(kind=4), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        complex(kind=4), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmcacoshi_omp_offload_ilp64

        subroutine mkl_vm_vzacoshi_omp_offload_ilp64(n, a, inca, y, incy)
        integer(kind=8), intent(in) :: n
        complex(kind=8), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        complex(kind=8), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        end subroutine mkl_vm_vzacoshi_omp_offload_ilp64

        subroutine mkl_vm_vmzacoshi_omp_offload_ilp64(n, a, inca, y, incy, mode)
        integer(kind=8), intent(in) :: n
        complex(kind=8), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        complex(kind=8), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmzacoshi_omp_offload_ilp64




! function: Acospi, indexing: simple
        subroutine mkl_vm_vsacospi_omp_offload_ilp64(n, a, y)
        integer(kind=8), intent(in) :: n
        real(kind=4), intent(in) :: a(*)
        real(kind=4), intent(inout) :: y(*)
        end subroutine mkl_vm_vsacospi_omp_offload_ilp64

        subroutine mkl_vm_vmsacospi_omp_offload_ilp64(n, a, y, mode)
        integer(kind=8), intent(in) :: n
        real(kind=4), intent(in) :: a(*)
        real(kind=4), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmsacospi_omp_offload_ilp64

        subroutine mkl_vm_vdacospi_omp_offload_ilp64(n, a, y)
        integer(kind=8), intent(in) :: n
        real(kind=8), intent(in) :: a(*)
        real(kind=8), intent(inout) :: y(*)
        end subroutine mkl_vm_vdacospi_omp_offload_ilp64

        subroutine mkl_vm_vmdacospi_omp_offload_ilp64(n, a, y, mode)
        integer(kind=8), intent(in) :: n
        real(kind=8), intent(in) :: a(*)
        real(kind=8), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmdacospi_omp_offload_ilp64


! function: Acospi, indexing: strided
        subroutine mkl_vm_vsacospii_omp_offload_ilp64(n, a, inca, y, incy)
        integer(kind=8), intent(in) :: n
        real(kind=4), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        real(kind=4), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        end subroutine mkl_vm_vsacospii_omp_offload_ilp64

        subroutine mkl_vm_vmsacospii_omp_offload_ilp64(n, a, inca, y, incy, mode)
        integer(kind=8), intent(in) :: n
        real(kind=4), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        real(kind=4), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmsacospii_omp_offload_ilp64

        subroutine mkl_vm_vdacospii_omp_offload_ilp64(n, a, inca, y, incy)
        integer(kind=8), intent(in) :: n
        real(kind=8), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        real(kind=8), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        end subroutine mkl_vm_vdacospii_omp_offload_ilp64

        subroutine mkl_vm_vmdacospii_omp_offload_ilp64(n, a, inca, y, incy, mode)
        integer(kind=8), intent(in) :: n
        real(kind=8), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        real(kind=8), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmdacospii_omp_offload_ilp64




! function: Add, indexing: simple
        subroutine mkl_vm_vsadd_omp_offload_ilp64(n, a, b, y)
        integer(kind=8), intent(in) :: n
        real(kind=4), intent(in) :: a(*)
        real(kind=4), intent(in) :: b(*)
        real(kind=4), intent(inout) :: y(*)
        end subroutine mkl_vm_vsadd_omp_offload_ilp64

        subroutine mkl_vm_vmsadd_omp_offload_ilp64(n, a, b, y, mode)
        integer(kind=8), intent(in) :: n
        real(kind=4), intent(in) :: a(*)
        real(kind=4), intent(in) :: b(*)
        real(kind=4), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmsadd_omp_offload_ilp64

        subroutine mkl_vm_vdadd_omp_offload_ilp64(n, a, b, y)
        integer(kind=8), intent(in) :: n
        real(kind=8), intent(in) :: a(*)
        real(kind=8), intent(in) :: b(*)
        real(kind=8), intent(inout) :: y(*)
        end subroutine mkl_vm_vdadd_omp_offload_ilp64

        subroutine mkl_vm_vmdadd_omp_offload_ilp64(n, a, b, y, mode)
        integer(kind=8), intent(in) :: n
        real(kind=8), intent(in) :: a(*)
        real(kind=8), intent(in) :: b(*)
        real(kind=8), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmdadd_omp_offload_ilp64

        subroutine mkl_vm_vcadd_omp_offload_ilp64(n, a, b, y)
        integer(kind=8), intent(in) :: n
        complex(kind=4), intent(in) :: a(*)
        complex(kind=4), intent(in) :: b(*)
        complex(kind=4), intent(inout) :: y(*)
        end subroutine mkl_vm_vcadd_omp_offload_ilp64

        subroutine mkl_vm_vmcadd_omp_offload_ilp64(n, a, b, y, mode)
        integer(kind=8), intent(in) :: n
        complex(kind=4), intent(in) :: a(*)
        complex(kind=4), intent(in) :: b(*)
        complex(kind=4), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmcadd_omp_offload_ilp64

        subroutine mkl_vm_vzadd_omp_offload_ilp64(n, a, b, y)
        integer(kind=8), intent(in) :: n
        complex(kind=8), intent(in) :: a(*)
        complex(kind=8), intent(in) :: b(*)
        complex(kind=8), intent(inout) :: y(*)
        end subroutine mkl_vm_vzadd_omp_offload_ilp64

        subroutine mkl_vm_vmzadd_omp_offload_ilp64(n, a, b, y, mode)
        integer(kind=8), intent(in) :: n
        complex(kind=8), intent(in) :: a(*)
        complex(kind=8), intent(in) :: b(*)
        complex(kind=8), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmzadd_omp_offload_ilp64


! function: Add, indexing: strided
        subroutine mkl_vm_vsaddi_omp_offload_ilp64(n, a, inca, b, incb, y, incy)
        integer(kind=8), intent(in) :: n
        real(kind=4), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        real(kind=4), intent(in) :: b(*)
        integer(kind=8), intent(in) :: incb
        real(kind=4), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        end subroutine mkl_vm_vsaddi_omp_offload_ilp64

        subroutine mkl_vm_vmsaddi_omp_offload_ilp64(n, a, inca, b, incb, y, incy, mode)
        integer(kind=8), intent(in) :: n
        real(kind=4), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        real(kind=4), intent(in) :: b(*)
        integer(kind=8), intent(in) :: incb
        real(kind=4), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmsaddi_omp_offload_ilp64

        subroutine mkl_vm_vdaddi_omp_offload_ilp64(n, a, inca, b, incb, y, incy)
        integer(kind=8), intent(in) :: n
        real(kind=8), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        real(kind=8), intent(in) :: b(*)
        integer(kind=8), intent(in) :: incb
        real(kind=8), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        end subroutine mkl_vm_vdaddi_omp_offload_ilp64

        subroutine mkl_vm_vmdaddi_omp_offload_ilp64(n, a, inca, b, incb, y, incy, mode)
        integer(kind=8), intent(in) :: n
        real(kind=8), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        real(kind=8), intent(in) :: b(*)
        integer(kind=8), intent(in) :: incb
        real(kind=8), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmdaddi_omp_offload_ilp64

        subroutine mkl_vm_vcaddi_omp_offload_ilp64(n, a, inca, b, incb, y, incy)
        integer(kind=8), intent(in) :: n
        complex(kind=4), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        complex(kind=4), intent(in) :: b(*)
        integer(kind=8), intent(in) :: incb
        complex(kind=4), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        end subroutine mkl_vm_vcaddi_omp_offload_ilp64

        subroutine mkl_vm_vmcaddi_omp_offload_ilp64(n, a, inca, b, incb, y, incy, mode)
        integer(kind=8), intent(in) :: n
        complex(kind=4), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        complex(kind=4), intent(in) :: b(*)
        integer(kind=8), intent(in) :: incb
        complex(kind=4), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmcaddi_omp_offload_ilp64

        subroutine mkl_vm_vzaddi_omp_offload_ilp64(n, a, inca, b, incb, y, incy)
        integer(kind=8), intent(in) :: n
        complex(kind=8), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        complex(kind=8), intent(in) :: b(*)
        integer(kind=8), intent(in) :: incb
        complex(kind=8), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        end subroutine mkl_vm_vzaddi_omp_offload_ilp64

        subroutine mkl_vm_vmzaddi_omp_offload_ilp64(n, a, inca, b, incb, y, incy, mode)
        integer(kind=8), intent(in) :: n
        complex(kind=8), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        complex(kind=8), intent(in) :: b(*)
        integer(kind=8), intent(in) :: incb
        complex(kind=8), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmzaddi_omp_offload_ilp64




! function: Arg, indexing: simple
        subroutine mkl_vm_vcarg_omp_offload_ilp64(n, a, y)
        integer(kind=8), intent(in) :: n
        complex(kind=4), intent(in) :: a(*)
        real(kind=4), intent(inout) :: y(*)
        end subroutine mkl_vm_vcarg_omp_offload_ilp64

        subroutine mkl_vm_vmcarg_omp_offload_ilp64(n, a, y, mode)
        integer(kind=8), intent(in) :: n
        complex(kind=4), intent(in) :: a(*)
        real(kind=4), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmcarg_omp_offload_ilp64

        subroutine mkl_vm_vzarg_omp_offload_ilp64(n, a, y)
        integer(kind=8), intent(in) :: n
        complex(kind=8), intent(in) :: a(*)
        real(kind=8), intent(inout) :: y(*)
        end subroutine mkl_vm_vzarg_omp_offload_ilp64

        subroutine mkl_vm_vmzarg_omp_offload_ilp64(n, a, y, mode)
        integer(kind=8), intent(in) :: n
        complex(kind=8), intent(in) :: a(*)
        real(kind=8), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmzarg_omp_offload_ilp64


! function: Arg, indexing: strided
        subroutine mkl_vm_vcargi_omp_offload_ilp64(n, a, inca, y, incy)
        integer(kind=8), intent(in) :: n
        complex(kind=4), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        real(kind=4), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        end subroutine mkl_vm_vcargi_omp_offload_ilp64

        subroutine mkl_vm_vmcargi_omp_offload_ilp64(n, a, inca, y, incy, mode)
        integer(kind=8), intent(in) :: n
        complex(kind=4), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        real(kind=4), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmcargi_omp_offload_ilp64

        subroutine mkl_vm_vzargi_omp_offload_ilp64(n, a, inca, y, incy)
        integer(kind=8), intent(in) :: n
        complex(kind=8), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        real(kind=8), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        end subroutine mkl_vm_vzargi_omp_offload_ilp64

        subroutine mkl_vm_vmzargi_omp_offload_ilp64(n, a, inca, y, incy, mode)
        integer(kind=8), intent(in) :: n
        complex(kind=8), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        real(kind=8), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmzargi_omp_offload_ilp64




! function: Asin, indexing: simple
        subroutine mkl_vm_vsasin_omp_offload_ilp64(n, a, y)
        integer(kind=8), intent(in) :: n
        real(kind=4), intent(in) :: a(*)
        real(kind=4), intent(inout) :: y(*)
        end subroutine mkl_vm_vsasin_omp_offload_ilp64

        subroutine mkl_vm_vmsasin_omp_offload_ilp64(n, a, y, mode)
        integer(kind=8), intent(in) :: n
        real(kind=4), intent(in) :: a(*)
        real(kind=4), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmsasin_omp_offload_ilp64

        subroutine mkl_vm_vdasin_omp_offload_ilp64(n, a, y)
        integer(kind=8), intent(in) :: n
        real(kind=8), intent(in) :: a(*)
        real(kind=8), intent(inout) :: y(*)
        end subroutine mkl_vm_vdasin_omp_offload_ilp64

        subroutine mkl_vm_vmdasin_omp_offload_ilp64(n, a, y, mode)
        integer(kind=8), intent(in) :: n
        real(kind=8), intent(in) :: a(*)
        real(kind=8), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmdasin_omp_offload_ilp64

        subroutine mkl_vm_vcasin_omp_offload_ilp64(n, a, y)
        integer(kind=8), intent(in) :: n
        complex(kind=4), intent(in) :: a(*)
        complex(kind=4), intent(inout) :: y(*)
        end subroutine mkl_vm_vcasin_omp_offload_ilp64

        subroutine mkl_vm_vmcasin_omp_offload_ilp64(n, a, y, mode)
        integer(kind=8), intent(in) :: n
        complex(kind=4), intent(in) :: a(*)
        complex(kind=4), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmcasin_omp_offload_ilp64

        subroutine mkl_vm_vzasin_omp_offload_ilp64(n, a, y)
        integer(kind=8), intent(in) :: n
        complex(kind=8), intent(in) :: a(*)
        complex(kind=8), intent(inout) :: y(*)
        end subroutine mkl_vm_vzasin_omp_offload_ilp64

        subroutine mkl_vm_vmzasin_omp_offload_ilp64(n, a, y, mode)
        integer(kind=8), intent(in) :: n
        complex(kind=8), intent(in) :: a(*)
        complex(kind=8), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmzasin_omp_offload_ilp64


! function: Asin, indexing: strided
        subroutine mkl_vm_vsasini_omp_offload_ilp64(n, a, inca, y, incy)
        integer(kind=8), intent(in) :: n
        real(kind=4), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        real(kind=4), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        end subroutine mkl_vm_vsasini_omp_offload_ilp64

        subroutine mkl_vm_vmsasini_omp_offload_ilp64(n, a, inca, y, incy, mode)
        integer(kind=8), intent(in) :: n
        real(kind=4), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        real(kind=4), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmsasini_omp_offload_ilp64

        subroutine mkl_vm_vdasini_omp_offload_ilp64(n, a, inca, y, incy)
        integer(kind=8), intent(in) :: n
        real(kind=8), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        real(kind=8), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        end subroutine mkl_vm_vdasini_omp_offload_ilp64

        subroutine mkl_vm_vmdasini_omp_offload_ilp64(n, a, inca, y, incy, mode)
        integer(kind=8), intent(in) :: n
        real(kind=8), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        real(kind=8), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmdasini_omp_offload_ilp64

        subroutine mkl_vm_vcasini_omp_offload_ilp64(n, a, inca, y, incy)
        integer(kind=8), intent(in) :: n
        complex(kind=4), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        complex(kind=4), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        end subroutine mkl_vm_vcasini_omp_offload_ilp64

        subroutine mkl_vm_vmcasini_omp_offload_ilp64(n, a, inca, y, incy, mode)
        integer(kind=8), intent(in) :: n
        complex(kind=4), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        complex(kind=4), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmcasini_omp_offload_ilp64

        subroutine mkl_vm_vzasini_omp_offload_ilp64(n, a, inca, y, incy)
        integer(kind=8), intent(in) :: n
        complex(kind=8), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        complex(kind=8), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        end subroutine mkl_vm_vzasini_omp_offload_ilp64

        subroutine mkl_vm_vmzasini_omp_offload_ilp64(n, a, inca, y, incy, mode)
        integer(kind=8), intent(in) :: n
        complex(kind=8), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        complex(kind=8), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmzasini_omp_offload_ilp64




! function: Asinh, indexing: simple
        subroutine mkl_vm_vsasinh_omp_offload_ilp64(n, a, y)
        integer(kind=8), intent(in) :: n
        real(kind=4), intent(in) :: a(*)
        real(kind=4), intent(inout) :: y(*)
        end subroutine mkl_vm_vsasinh_omp_offload_ilp64

        subroutine mkl_vm_vmsasinh_omp_offload_ilp64(n, a, y, mode)
        integer(kind=8), intent(in) :: n
        real(kind=4), intent(in) :: a(*)
        real(kind=4), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmsasinh_omp_offload_ilp64

        subroutine mkl_vm_vdasinh_omp_offload_ilp64(n, a, y)
        integer(kind=8), intent(in) :: n
        real(kind=8), intent(in) :: a(*)
        real(kind=8), intent(inout) :: y(*)
        end subroutine mkl_vm_vdasinh_omp_offload_ilp64

        subroutine mkl_vm_vmdasinh_omp_offload_ilp64(n, a, y, mode)
        integer(kind=8), intent(in) :: n
        real(kind=8), intent(in) :: a(*)
        real(kind=8), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmdasinh_omp_offload_ilp64

        subroutine mkl_vm_vcasinh_omp_offload_ilp64(n, a, y)
        integer(kind=8), intent(in) :: n
        complex(kind=4), intent(in) :: a(*)
        complex(kind=4), intent(inout) :: y(*)
        end subroutine mkl_vm_vcasinh_omp_offload_ilp64

        subroutine mkl_vm_vmcasinh_omp_offload_ilp64(n, a, y, mode)
        integer(kind=8), intent(in) :: n
        complex(kind=4), intent(in) :: a(*)
        complex(kind=4), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmcasinh_omp_offload_ilp64

        subroutine mkl_vm_vzasinh_omp_offload_ilp64(n, a, y)
        integer(kind=8), intent(in) :: n
        complex(kind=8), intent(in) :: a(*)
        complex(kind=8), intent(inout) :: y(*)
        end subroutine mkl_vm_vzasinh_omp_offload_ilp64

        subroutine mkl_vm_vmzasinh_omp_offload_ilp64(n, a, y, mode)
        integer(kind=8), intent(in) :: n
        complex(kind=8), intent(in) :: a(*)
        complex(kind=8), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmzasinh_omp_offload_ilp64


! function: Asinh, indexing: strided
        subroutine mkl_vm_vsasinhi_omp_offload_ilp64(n, a, inca, y, incy)
        integer(kind=8), intent(in) :: n
        real(kind=4), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        real(kind=4), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        end subroutine mkl_vm_vsasinhi_omp_offload_ilp64

        subroutine mkl_vm_vmsasinhi_omp_offload_ilp64(n, a, inca, y, incy, mode)
        integer(kind=8), intent(in) :: n
        real(kind=4), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        real(kind=4), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmsasinhi_omp_offload_ilp64

        subroutine mkl_vm_vdasinhi_omp_offload_ilp64(n, a, inca, y, incy)
        integer(kind=8), intent(in) :: n
        real(kind=8), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        real(kind=8), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        end subroutine mkl_vm_vdasinhi_omp_offload_ilp64

        subroutine mkl_vm_vmdasinhi_omp_offload_ilp64(n, a, inca, y, incy, mode)
        integer(kind=8), intent(in) :: n
        real(kind=8), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        real(kind=8), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmdasinhi_omp_offload_ilp64

        subroutine mkl_vm_vcasinhi_omp_offload_ilp64(n, a, inca, y, incy)
        integer(kind=8), intent(in) :: n
        complex(kind=4), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        complex(kind=4), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        end subroutine mkl_vm_vcasinhi_omp_offload_ilp64

        subroutine mkl_vm_vmcasinhi_omp_offload_ilp64(n, a, inca, y, incy, mode)
        integer(kind=8), intent(in) :: n
        complex(kind=4), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        complex(kind=4), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmcasinhi_omp_offload_ilp64

        subroutine mkl_vm_vzasinhi_omp_offload_ilp64(n, a, inca, y, incy)
        integer(kind=8), intent(in) :: n
        complex(kind=8), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        complex(kind=8), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        end subroutine mkl_vm_vzasinhi_omp_offload_ilp64

        subroutine mkl_vm_vmzasinhi_omp_offload_ilp64(n, a, inca, y, incy, mode)
        integer(kind=8), intent(in) :: n
        complex(kind=8), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        complex(kind=8), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmzasinhi_omp_offload_ilp64




! function: Asinpi, indexing: simple
        subroutine mkl_vm_vsasinpi_omp_offload_ilp64(n, a, y)
        integer(kind=8), intent(in) :: n
        real(kind=4), intent(in) :: a(*)
        real(kind=4), intent(inout) :: y(*)
        end subroutine mkl_vm_vsasinpi_omp_offload_ilp64

        subroutine mkl_vm_vmsasinpi_omp_offload_ilp64(n, a, y, mode)
        integer(kind=8), intent(in) :: n
        real(kind=4), intent(in) :: a(*)
        real(kind=4), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmsasinpi_omp_offload_ilp64

        subroutine mkl_vm_vdasinpi_omp_offload_ilp64(n, a, y)
        integer(kind=8), intent(in) :: n
        real(kind=8), intent(in) :: a(*)
        real(kind=8), intent(inout) :: y(*)
        end subroutine mkl_vm_vdasinpi_omp_offload_ilp64

        subroutine mkl_vm_vmdasinpi_omp_offload_ilp64(n, a, y, mode)
        integer(kind=8), intent(in) :: n
        real(kind=8), intent(in) :: a(*)
        real(kind=8), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmdasinpi_omp_offload_ilp64


! function: Asinpi, indexing: strided
        subroutine mkl_vm_vsasinpii_omp_offload_ilp64(n, a, inca, y, incy)
        integer(kind=8), intent(in) :: n
        real(kind=4), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        real(kind=4), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        end subroutine mkl_vm_vsasinpii_omp_offload_ilp64

        subroutine mkl_vm_vmsasinpii_omp_offload_ilp64(n, a, inca, y, incy, mode)
        integer(kind=8), intent(in) :: n
        real(kind=4), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        real(kind=4), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmsasinpii_omp_offload_ilp64

        subroutine mkl_vm_vdasinpii_omp_offload_ilp64(n, a, inca, y, incy)
        integer(kind=8), intent(in) :: n
        real(kind=8), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        real(kind=8), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        end subroutine mkl_vm_vdasinpii_omp_offload_ilp64

        subroutine mkl_vm_vmdasinpii_omp_offload_ilp64(n, a, inca, y, incy, mode)
        integer(kind=8), intent(in) :: n
        real(kind=8), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        real(kind=8), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmdasinpii_omp_offload_ilp64




! function: Atan, indexing: simple
        subroutine mkl_vm_vsatan_omp_offload_ilp64(n, a, y)
        integer(kind=8), intent(in) :: n
        real(kind=4), intent(in) :: a(*)
        real(kind=4), intent(inout) :: y(*)
        end subroutine mkl_vm_vsatan_omp_offload_ilp64

        subroutine mkl_vm_vmsatan_omp_offload_ilp64(n, a, y, mode)
        integer(kind=8), intent(in) :: n
        real(kind=4), intent(in) :: a(*)
        real(kind=4), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmsatan_omp_offload_ilp64

        subroutine mkl_vm_vdatan_omp_offload_ilp64(n, a, y)
        integer(kind=8), intent(in) :: n
        real(kind=8), intent(in) :: a(*)
        real(kind=8), intent(inout) :: y(*)
        end subroutine mkl_vm_vdatan_omp_offload_ilp64

        subroutine mkl_vm_vmdatan_omp_offload_ilp64(n, a, y, mode)
        integer(kind=8), intent(in) :: n
        real(kind=8), intent(in) :: a(*)
        real(kind=8), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmdatan_omp_offload_ilp64

        subroutine mkl_vm_vcatan_omp_offload_ilp64(n, a, y)
        integer(kind=8), intent(in) :: n
        complex(kind=4), intent(in) :: a(*)
        complex(kind=4), intent(inout) :: y(*)
        end subroutine mkl_vm_vcatan_omp_offload_ilp64

        subroutine mkl_vm_vmcatan_omp_offload_ilp64(n, a, y, mode)
        integer(kind=8), intent(in) :: n
        complex(kind=4), intent(in) :: a(*)
        complex(kind=4), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmcatan_omp_offload_ilp64

        subroutine mkl_vm_vzatan_omp_offload_ilp64(n, a, y)
        integer(kind=8), intent(in) :: n
        complex(kind=8), intent(in) :: a(*)
        complex(kind=8), intent(inout) :: y(*)
        end subroutine mkl_vm_vzatan_omp_offload_ilp64

        subroutine mkl_vm_vmzatan_omp_offload_ilp64(n, a, y, mode)
        integer(kind=8), intent(in) :: n
        complex(kind=8), intent(in) :: a(*)
        complex(kind=8), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmzatan_omp_offload_ilp64


! function: Atan, indexing: strided
        subroutine mkl_vm_vsatani_omp_offload_ilp64(n, a, inca, y, incy)
        integer(kind=8), intent(in) :: n
        real(kind=4), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        real(kind=4), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        end subroutine mkl_vm_vsatani_omp_offload_ilp64

        subroutine mkl_vm_vmsatani_omp_offload_ilp64(n, a, inca, y, incy, mode)
        integer(kind=8), intent(in) :: n
        real(kind=4), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        real(kind=4), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmsatani_omp_offload_ilp64

        subroutine mkl_vm_vdatani_omp_offload_ilp64(n, a, inca, y, incy)
        integer(kind=8), intent(in) :: n
        real(kind=8), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        real(kind=8), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        end subroutine mkl_vm_vdatani_omp_offload_ilp64

        subroutine mkl_vm_vmdatani_omp_offload_ilp64(n, a, inca, y, incy, mode)
        integer(kind=8), intent(in) :: n
        real(kind=8), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        real(kind=8), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmdatani_omp_offload_ilp64

        subroutine mkl_vm_vcatani_omp_offload_ilp64(n, a, inca, y, incy)
        integer(kind=8), intent(in) :: n
        complex(kind=4), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        complex(kind=4), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        end subroutine mkl_vm_vcatani_omp_offload_ilp64

        subroutine mkl_vm_vmcatani_omp_offload_ilp64(n, a, inca, y, incy, mode)
        integer(kind=8), intent(in) :: n
        complex(kind=4), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        complex(kind=4), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmcatani_omp_offload_ilp64

        subroutine mkl_vm_vzatani_omp_offload_ilp64(n, a, inca, y, incy)
        integer(kind=8), intent(in) :: n
        complex(kind=8), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        complex(kind=8), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        end subroutine mkl_vm_vzatani_omp_offload_ilp64

        subroutine mkl_vm_vmzatani_omp_offload_ilp64(n, a, inca, y, incy, mode)
        integer(kind=8), intent(in) :: n
        complex(kind=8), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        complex(kind=8), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmzatani_omp_offload_ilp64




! function: Atan2, indexing: simple
        subroutine mkl_vm_vsatan2_omp_offload_ilp64(n, a, b, y)
        integer(kind=8), intent(in) :: n
        real(kind=4), intent(in) :: a(*)
        real(kind=4), intent(in) :: b(*)
        real(kind=4), intent(inout) :: y(*)
        end subroutine mkl_vm_vsatan2_omp_offload_ilp64

        subroutine mkl_vm_vmsatan2_omp_offload_ilp64(n, a, b, y, mode)
        integer(kind=8), intent(in) :: n
        real(kind=4), intent(in) :: a(*)
        real(kind=4), intent(in) :: b(*)
        real(kind=4), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmsatan2_omp_offload_ilp64

        subroutine mkl_vm_vdatan2_omp_offload_ilp64(n, a, b, y)
        integer(kind=8), intent(in) :: n
        real(kind=8), intent(in) :: a(*)
        real(kind=8), intent(in) :: b(*)
        real(kind=8), intent(inout) :: y(*)
        end subroutine mkl_vm_vdatan2_omp_offload_ilp64

        subroutine mkl_vm_vmdatan2_omp_offload_ilp64(n, a, b, y, mode)
        integer(kind=8), intent(in) :: n
        real(kind=8), intent(in) :: a(*)
        real(kind=8), intent(in) :: b(*)
        real(kind=8), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmdatan2_omp_offload_ilp64


! function: Atan2, indexing: strided
        subroutine mkl_vm_vsatan2i_omp_offload_ilp64(n, a, inca, b, incb, y, incy)
        integer(kind=8), intent(in) :: n
        real(kind=4), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        real(kind=4), intent(in) :: b(*)
        integer(kind=8), intent(in) :: incb
        real(kind=4), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        end subroutine mkl_vm_vsatan2i_omp_offload_ilp64

        subroutine mkl_vm_vmsatan2i_omp_offload_ilp64(n, a, inca, b, incb, y, incy, mode)
        integer(kind=8), intent(in) :: n
        real(kind=4), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        real(kind=4), intent(in) :: b(*)
        integer(kind=8), intent(in) :: incb
        real(kind=4), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmsatan2i_omp_offload_ilp64

        subroutine mkl_vm_vdatan2i_omp_offload_ilp64(n, a, inca, b, incb, y, incy)
        integer(kind=8), intent(in) :: n
        real(kind=8), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        real(kind=8), intent(in) :: b(*)
        integer(kind=8), intent(in) :: incb
        real(kind=8), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        end subroutine mkl_vm_vdatan2i_omp_offload_ilp64

        subroutine mkl_vm_vmdatan2i_omp_offload_ilp64(n, a, inca, b, incb, y, incy, mode)
        integer(kind=8), intent(in) :: n
        real(kind=8), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        real(kind=8), intent(in) :: b(*)
        integer(kind=8), intent(in) :: incb
        real(kind=8), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmdatan2i_omp_offload_ilp64




! function: Atan2pi, indexing: simple
        subroutine mkl_vm_vsatan2pi_omp_offload_ilp64(n, a, b, y)
        integer(kind=8), intent(in) :: n
        real(kind=4), intent(in) :: a(*)
        real(kind=4), intent(in) :: b(*)
        real(kind=4), intent(inout) :: y(*)
        end subroutine mkl_vm_vsatan2pi_omp_offload_ilp64

        subroutine mkl_vm_vmsatan2pi_omp_offload_ilp64(n, a, b, y, mode)
        integer(kind=8), intent(in) :: n
        real(kind=4), intent(in) :: a(*)
        real(kind=4), intent(in) :: b(*)
        real(kind=4), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmsatan2pi_omp_offload_ilp64

        subroutine mkl_vm_vdatan2pi_omp_offload_ilp64(n, a, b, y)
        integer(kind=8), intent(in) :: n
        real(kind=8), intent(in) :: a(*)
        real(kind=8), intent(in) :: b(*)
        real(kind=8), intent(inout) :: y(*)
        end subroutine mkl_vm_vdatan2pi_omp_offload_ilp64

        subroutine mkl_vm_vmdatan2pi_omp_offload_ilp64(n, a, b, y, mode)
        integer(kind=8), intent(in) :: n
        real(kind=8), intent(in) :: a(*)
        real(kind=8), intent(in) :: b(*)
        real(kind=8), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmdatan2pi_omp_offload_ilp64


! function: Atan2pi, indexing: strided
        subroutine mkl_vm_vsatan2pii_omp_offload_ilp64(n, a, inca, b, incb, y, incy)
        integer(kind=8), intent(in) :: n
        real(kind=4), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        real(kind=4), intent(in) :: b(*)
        integer(kind=8), intent(in) :: incb
        real(kind=4), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        end subroutine mkl_vm_vsatan2pii_omp_offload_ilp64

        subroutine mkl_vm_vmsatan2pii_omp_offload_ilp64(n, a, inca, b, incb, y, incy, mode)
        integer(kind=8), intent(in) :: n
        real(kind=4), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        real(kind=4), intent(in) :: b(*)
        integer(kind=8), intent(in) :: incb
        real(kind=4), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmsatan2pii_omp_offload_ilp64

        subroutine mkl_vm_vdatan2pii_omp_offload_ilp64(n, a, inca, b, incb, y, incy)
        integer(kind=8), intent(in) :: n
        real(kind=8), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        real(kind=8), intent(in) :: b(*)
        integer(kind=8), intent(in) :: incb
        real(kind=8), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        end subroutine mkl_vm_vdatan2pii_omp_offload_ilp64

        subroutine mkl_vm_vmdatan2pii_omp_offload_ilp64(n, a, inca, b, incb, y, incy, mode)
        integer(kind=8), intent(in) :: n
        real(kind=8), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        real(kind=8), intent(in) :: b(*)
        integer(kind=8), intent(in) :: incb
        real(kind=8), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmdatan2pii_omp_offload_ilp64




! function: Atanh, indexing: simple
        subroutine mkl_vm_vsatanh_omp_offload_ilp64(n, a, y)
        integer(kind=8), intent(in) :: n
        real(kind=4), intent(in) :: a(*)
        real(kind=4), intent(inout) :: y(*)
        end subroutine mkl_vm_vsatanh_omp_offload_ilp64

        subroutine mkl_vm_vmsatanh_omp_offload_ilp64(n, a, y, mode)
        integer(kind=8), intent(in) :: n
        real(kind=4), intent(in) :: a(*)
        real(kind=4), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmsatanh_omp_offload_ilp64

        subroutine mkl_vm_vdatanh_omp_offload_ilp64(n, a, y)
        integer(kind=8), intent(in) :: n
        real(kind=8), intent(in) :: a(*)
        real(kind=8), intent(inout) :: y(*)
        end subroutine mkl_vm_vdatanh_omp_offload_ilp64

        subroutine mkl_vm_vmdatanh_omp_offload_ilp64(n, a, y, mode)
        integer(kind=8), intent(in) :: n
        real(kind=8), intent(in) :: a(*)
        real(kind=8), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmdatanh_omp_offload_ilp64

        subroutine mkl_vm_vcatanh_omp_offload_ilp64(n, a, y)
        integer(kind=8), intent(in) :: n
        complex(kind=4), intent(in) :: a(*)
        complex(kind=4), intent(inout) :: y(*)
        end subroutine mkl_vm_vcatanh_omp_offload_ilp64

        subroutine mkl_vm_vmcatanh_omp_offload_ilp64(n, a, y, mode)
        integer(kind=8), intent(in) :: n
        complex(kind=4), intent(in) :: a(*)
        complex(kind=4), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmcatanh_omp_offload_ilp64

        subroutine mkl_vm_vzatanh_omp_offload_ilp64(n, a, y)
        integer(kind=8), intent(in) :: n
        complex(kind=8), intent(in) :: a(*)
        complex(kind=8), intent(inout) :: y(*)
        end subroutine mkl_vm_vzatanh_omp_offload_ilp64

        subroutine mkl_vm_vmzatanh_omp_offload_ilp64(n, a, y, mode)
        integer(kind=8), intent(in) :: n
        complex(kind=8), intent(in) :: a(*)
        complex(kind=8), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmzatanh_omp_offload_ilp64


! function: Atanh, indexing: strided
        subroutine mkl_vm_vsatanhi_omp_offload_ilp64(n, a, inca, y, incy)
        integer(kind=8), intent(in) :: n
        real(kind=4), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        real(kind=4), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        end subroutine mkl_vm_vsatanhi_omp_offload_ilp64

        subroutine mkl_vm_vmsatanhi_omp_offload_ilp64(n, a, inca, y, incy, mode)
        integer(kind=8), intent(in) :: n
        real(kind=4), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        real(kind=4), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmsatanhi_omp_offload_ilp64

        subroutine mkl_vm_vdatanhi_omp_offload_ilp64(n, a, inca, y, incy)
        integer(kind=8), intent(in) :: n
        real(kind=8), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        real(kind=8), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        end subroutine mkl_vm_vdatanhi_omp_offload_ilp64

        subroutine mkl_vm_vmdatanhi_omp_offload_ilp64(n, a, inca, y, incy, mode)
        integer(kind=8), intent(in) :: n
        real(kind=8), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        real(kind=8), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmdatanhi_omp_offload_ilp64

        subroutine mkl_vm_vcatanhi_omp_offload_ilp64(n, a, inca, y, incy)
        integer(kind=8), intent(in) :: n
        complex(kind=4), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        complex(kind=4), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        end subroutine mkl_vm_vcatanhi_omp_offload_ilp64

        subroutine mkl_vm_vmcatanhi_omp_offload_ilp64(n, a, inca, y, incy, mode)
        integer(kind=8), intent(in) :: n
        complex(kind=4), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        complex(kind=4), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmcatanhi_omp_offload_ilp64

        subroutine mkl_vm_vzatanhi_omp_offload_ilp64(n, a, inca, y, incy)
        integer(kind=8), intent(in) :: n
        complex(kind=8), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        complex(kind=8), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        end subroutine mkl_vm_vzatanhi_omp_offload_ilp64

        subroutine mkl_vm_vmzatanhi_omp_offload_ilp64(n, a, inca, y, incy, mode)
        integer(kind=8), intent(in) :: n
        complex(kind=8), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        complex(kind=8), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmzatanhi_omp_offload_ilp64




! function: Atanpi, indexing: simple
        subroutine mkl_vm_vsatanpi_omp_offload_ilp64(n, a, y)
        integer(kind=8), intent(in) :: n
        real(kind=4), intent(in) :: a(*)
        real(kind=4), intent(inout) :: y(*)
        end subroutine mkl_vm_vsatanpi_omp_offload_ilp64

        subroutine mkl_vm_vmsatanpi_omp_offload_ilp64(n, a, y, mode)
        integer(kind=8), intent(in) :: n
        real(kind=4), intent(in) :: a(*)
        real(kind=4), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmsatanpi_omp_offload_ilp64

        subroutine mkl_vm_vdatanpi_omp_offload_ilp64(n, a, y)
        integer(kind=8), intent(in) :: n
        real(kind=8), intent(in) :: a(*)
        real(kind=8), intent(inout) :: y(*)
        end subroutine mkl_vm_vdatanpi_omp_offload_ilp64

        subroutine mkl_vm_vmdatanpi_omp_offload_ilp64(n, a, y, mode)
        integer(kind=8), intent(in) :: n
        real(kind=8), intent(in) :: a(*)
        real(kind=8), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmdatanpi_omp_offload_ilp64


! function: Atanpi, indexing: strided
        subroutine mkl_vm_vsatanpii_omp_offload_ilp64(n, a, inca, y, incy)
        integer(kind=8), intent(in) :: n
        real(kind=4), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        real(kind=4), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        end subroutine mkl_vm_vsatanpii_omp_offload_ilp64

        subroutine mkl_vm_vmsatanpii_omp_offload_ilp64(n, a, inca, y, incy, mode)
        integer(kind=8), intent(in) :: n
        real(kind=4), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        real(kind=4), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmsatanpii_omp_offload_ilp64

        subroutine mkl_vm_vdatanpii_omp_offload_ilp64(n, a, inca, y, incy)
        integer(kind=8), intent(in) :: n
        real(kind=8), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        real(kind=8), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        end subroutine mkl_vm_vdatanpii_omp_offload_ilp64

        subroutine mkl_vm_vmdatanpii_omp_offload_ilp64(n, a, inca, y, incy, mode)
        integer(kind=8), intent(in) :: n
        real(kind=8), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        real(kind=8), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmdatanpii_omp_offload_ilp64




! function: Cbrt, indexing: simple
        subroutine mkl_vm_vscbrt_omp_offload_ilp64(n, a, y)
        integer(kind=8), intent(in) :: n
        real(kind=4), intent(in) :: a(*)
        real(kind=4), intent(inout) :: y(*)
        end subroutine mkl_vm_vscbrt_omp_offload_ilp64

        subroutine mkl_vm_vmscbrt_omp_offload_ilp64(n, a, y, mode)
        integer(kind=8), intent(in) :: n
        real(kind=4), intent(in) :: a(*)
        real(kind=4), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmscbrt_omp_offload_ilp64

        subroutine mkl_vm_vdcbrt_omp_offload_ilp64(n, a, y)
        integer(kind=8), intent(in) :: n
        real(kind=8), intent(in) :: a(*)
        real(kind=8), intent(inout) :: y(*)
        end subroutine mkl_vm_vdcbrt_omp_offload_ilp64

        subroutine mkl_vm_vmdcbrt_omp_offload_ilp64(n, a, y, mode)
        integer(kind=8), intent(in) :: n
        real(kind=8), intent(in) :: a(*)
        real(kind=8), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmdcbrt_omp_offload_ilp64


! function: Cbrt, indexing: strided
        subroutine mkl_vm_vscbrti_omp_offload_ilp64(n, a, inca, y, incy)
        integer(kind=8), intent(in) :: n
        real(kind=4), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        real(kind=4), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        end subroutine mkl_vm_vscbrti_omp_offload_ilp64

        subroutine mkl_vm_vmscbrti_omp_offload_ilp64(n, a, inca, y, incy, mode)
        integer(kind=8), intent(in) :: n
        real(kind=4), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        real(kind=4), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmscbrti_omp_offload_ilp64

        subroutine mkl_vm_vdcbrti_omp_offload_ilp64(n, a, inca, y, incy)
        integer(kind=8), intent(in) :: n
        real(kind=8), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        real(kind=8), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        end subroutine mkl_vm_vdcbrti_omp_offload_ilp64

        subroutine mkl_vm_vmdcbrti_omp_offload_ilp64(n, a, inca, y, incy, mode)
        integer(kind=8), intent(in) :: n
        real(kind=8), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        real(kind=8), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmdcbrti_omp_offload_ilp64




! function: CdfNorm, indexing: simple
        subroutine mkl_vm_vscdfnorm_omp_offload_ilp64(n, a, y)
        integer(kind=8), intent(in) :: n
        real(kind=4), intent(in) :: a(*)
        real(kind=4), intent(inout) :: y(*)
        end subroutine mkl_vm_vscdfnorm_omp_offload_ilp64

        subroutine mkl_vm_vmscdfnorm_omp_offload_ilp64(n, a, y, mode)
        integer(kind=8), intent(in) :: n
        real(kind=4), intent(in) :: a(*)
        real(kind=4), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmscdfnorm_omp_offload_ilp64

        subroutine mkl_vm_vdcdfnorm_omp_offload_ilp64(n, a, y)
        integer(kind=8), intent(in) :: n
        real(kind=8), intent(in) :: a(*)
        real(kind=8), intent(inout) :: y(*)
        end subroutine mkl_vm_vdcdfnorm_omp_offload_ilp64

        subroutine mkl_vm_vmdcdfnorm_omp_offload_ilp64(n, a, y, mode)
        integer(kind=8), intent(in) :: n
        real(kind=8), intent(in) :: a(*)
        real(kind=8), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmdcdfnorm_omp_offload_ilp64


! function: CdfNorm, indexing: strided
        subroutine mkl_vm_vscdfnormi_omp_offload_ilp64(n, a, inca, y, incy)
        integer(kind=8), intent(in) :: n
        real(kind=4), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        real(kind=4), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        end subroutine mkl_vm_vscdfnormi_omp_offload_ilp64

        subroutine mkl_vm_vmscdfnormi_omp_offload_ilp64(n, a, inca, y, incy, mode)
        integer(kind=8), intent(in) :: n
        real(kind=4), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        real(kind=4), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmscdfnormi_omp_offload_ilp64

        subroutine mkl_vm_vdcdfnormi_omp_offload_ilp64(n, a, inca, y, incy)
        integer(kind=8), intent(in) :: n
        real(kind=8), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        real(kind=8), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        end subroutine mkl_vm_vdcdfnormi_omp_offload_ilp64

        subroutine mkl_vm_vmdcdfnormi_omp_offload_ilp64(n, a, inca, y, incy, mode)
        integer(kind=8), intent(in) :: n
        real(kind=8), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        real(kind=8), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmdcdfnormi_omp_offload_ilp64




! function: CdfNormInv, indexing: simple
        subroutine mkl_vm_vscdfnorminv_omp_offload_ilp64(n, a, y)
        integer(kind=8), intent(in) :: n
        real(kind=4), intent(in) :: a(*)
        real(kind=4), intent(inout) :: y(*)
        end subroutine mkl_vm_vscdfnorminv_omp_offload_ilp64

        subroutine mkl_vm_vmscdfnorminv_omp_offload_ilp64(n, a, y, mode)
        integer(kind=8), intent(in) :: n
        real(kind=4), intent(in) :: a(*)
        real(kind=4), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmscdfnorminv_omp_offload_ilp64

        subroutine mkl_vm_vdcdfnorminv_omp_offload_ilp64(n, a, y)
        integer(kind=8), intent(in) :: n
        real(kind=8), intent(in) :: a(*)
        real(kind=8), intent(inout) :: y(*)
        end subroutine mkl_vm_vdcdfnorminv_omp_offload_ilp64

        subroutine mkl_vm_vmdcdfnorminv_omp_offload_ilp64(n, a, y, mode)
        integer(kind=8), intent(in) :: n
        real(kind=8), intent(in) :: a(*)
        real(kind=8), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmdcdfnorminv_omp_offload_ilp64


! function: CdfNormInv, indexing: strided
        subroutine mkl_vm_vscdfnorminvi_omp_offload_ilp64(n, a, inca, y, incy)
        integer(kind=8), intent(in) :: n
        real(kind=4), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        real(kind=4), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        end subroutine mkl_vm_vscdfnorminvi_omp_offload_ilp64

        subroutine mkl_vm_vmscdfnorminvi_omp_offload_ilp64(n, a, inca, y, incy, mode)
        integer(kind=8), intent(in) :: n
        real(kind=4), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        real(kind=4), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmscdfnorminvi_omp_offload_ilp64

        subroutine mkl_vm_vdcdfnorminvi_omp_offload_ilp64(n, a, inca, y, incy)
        integer(kind=8), intent(in) :: n
        real(kind=8), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        real(kind=8), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        end subroutine mkl_vm_vdcdfnorminvi_omp_offload_ilp64

        subroutine mkl_vm_vmdcdfnorminvi_omp_offload_ilp64(n, a, inca, y, incy, mode)
        integer(kind=8), intent(in) :: n
        real(kind=8), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        real(kind=8), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmdcdfnorminvi_omp_offload_ilp64




! function: Ceil, indexing: simple
        subroutine mkl_vm_vsceil_omp_offload_ilp64(n, a, y)
        integer(kind=8), intent(in) :: n
        real(kind=4), intent(in) :: a(*)
        real(kind=4), intent(inout) :: y(*)
        end subroutine mkl_vm_vsceil_omp_offload_ilp64

        subroutine mkl_vm_vmsceil_omp_offload_ilp64(n, a, y, mode)
        integer(kind=8), intent(in) :: n
        real(kind=4), intent(in) :: a(*)
        real(kind=4), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmsceil_omp_offload_ilp64

        subroutine mkl_vm_vdceil_omp_offload_ilp64(n, a, y)
        integer(kind=8), intent(in) :: n
        real(kind=8), intent(in) :: a(*)
        real(kind=8), intent(inout) :: y(*)
        end subroutine mkl_vm_vdceil_omp_offload_ilp64

        subroutine mkl_vm_vmdceil_omp_offload_ilp64(n, a, y, mode)
        integer(kind=8), intent(in) :: n
        real(kind=8), intent(in) :: a(*)
        real(kind=8), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmdceil_omp_offload_ilp64


! function: Ceil, indexing: strided
        subroutine mkl_vm_vsceili_omp_offload_ilp64(n, a, inca, y, incy)
        integer(kind=8), intent(in) :: n
        real(kind=4), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        real(kind=4), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        end subroutine mkl_vm_vsceili_omp_offload_ilp64

        subroutine mkl_vm_vmsceili_omp_offload_ilp64(n, a, inca, y, incy, mode)
        integer(kind=8), intent(in) :: n
        real(kind=4), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        real(kind=4), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmsceili_omp_offload_ilp64

        subroutine mkl_vm_vdceili_omp_offload_ilp64(n, a, inca, y, incy)
        integer(kind=8), intent(in) :: n
        real(kind=8), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        real(kind=8), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        end subroutine mkl_vm_vdceili_omp_offload_ilp64

        subroutine mkl_vm_vmdceili_omp_offload_ilp64(n, a, inca, y, incy, mode)
        integer(kind=8), intent(in) :: n
        real(kind=8), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        real(kind=8), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmdceili_omp_offload_ilp64




! function: CIS, indexing: simple
        subroutine mkl_vm_vccis_omp_offload_ilp64(n, a, y)
        integer(kind=8), intent(in) :: n
        real(kind=4), intent(in) :: a(*)
        complex(kind=4), intent(inout) :: y(*)
        end subroutine mkl_vm_vccis_omp_offload_ilp64

        subroutine mkl_vm_vmccis_omp_offload_ilp64(n, a, y, mode)
        integer(kind=8), intent(in) :: n
        real(kind=4), intent(in) :: a(*)
        complex(kind=4), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmccis_omp_offload_ilp64

        subroutine mkl_vm_vzcis_omp_offload_ilp64(n, a, y)
        integer(kind=8), intent(in) :: n
        real(kind=8), intent(in) :: a(*)
        complex(kind=8), intent(inout) :: y(*)
        end subroutine mkl_vm_vzcis_omp_offload_ilp64

        subroutine mkl_vm_vmzcis_omp_offload_ilp64(n, a, y, mode)
        integer(kind=8), intent(in) :: n
        real(kind=8), intent(in) :: a(*)
        complex(kind=8), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmzcis_omp_offload_ilp64


! function: CIS, indexing: strided
        subroutine mkl_vm_vccisi_omp_offload_ilp64(n, a, inca, y, incy)
        integer(kind=8), intent(in) :: n
        real(kind=4), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        complex(kind=4), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        end subroutine mkl_vm_vccisi_omp_offload_ilp64

        subroutine mkl_vm_vmccisi_omp_offload_ilp64(n, a, inca, y, incy, mode)
        integer(kind=8), intent(in) :: n
        real(kind=4), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        complex(kind=4), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmccisi_omp_offload_ilp64

        subroutine mkl_vm_vzcisi_omp_offload_ilp64(n, a, inca, y, incy)
        integer(kind=8), intent(in) :: n
        real(kind=8), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        complex(kind=8), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        end subroutine mkl_vm_vzcisi_omp_offload_ilp64

        subroutine mkl_vm_vmzcisi_omp_offload_ilp64(n, a, inca, y, incy, mode)
        integer(kind=8), intent(in) :: n
        real(kind=8), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        complex(kind=8), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmzcisi_omp_offload_ilp64




! function: Conj, indexing: simple
        subroutine mkl_vm_vcconj_omp_offload_ilp64(n, a, y)
        integer(kind=8), intent(in) :: n
        complex(kind=4), intent(in) :: a(*)
        complex(kind=4), intent(inout) :: y(*)
        end subroutine mkl_vm_vcconj_omp_offload_ilp64

        subroutine mkl_vm_vmcconj_omp_offload_ilp64(n, a, y, mode)
        integer(kind=8), intent(in) :: n
        complex(kind=4), intent(in) :: a(*)
        complex(kind=4), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmcconj_omp_offload_ilp64

        subroutine mkl_vm_vzconj_omp_offload_ilp64(n, a, y)
        integer(kind=8), intent(in) :: n
        complex(kind=8), intent(in) :: a(*)
        complex(kind=8), intent(inout) :: y(*)
        end subroutine mkl_vm_vzconj_omp_offload_ilp64

        subroutine mkl_vm_vmzconj_omp_offload_ilp64(n, a, y, mode)
        integer(kind=8), intent(in) :: n
        complex(kind=8), intent(in) :: a(*)
        complex(kind=8), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmzconj_omp_offload_ilp64


! function: Conj, indexing: strided
        subroutine mkl_vm_vcconji_omp_offload_ilp64(n, a, inca, y, incy)
        integer(kind=8), intent(in) :: n
        complex(kind=4), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        complex(kind=4), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        end subroutine mkl_vm_vcconji_omp_offload_ilp64

        subroutine mkl_vm_vmcconji_omp_offload_ilp64(n, a, inca, y, incy, mode)
        integer(kind=8), intent(in) :: n
        complex(kind=4), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        complex(kind=4), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmcconji_omp_offload_ilp64

        subroutine mkl_vm_vzconji_omp_offload_ilp64(n, a, inca, y, incy)
        integer(kind=8), intent(in) :: n
        complex(kind=8), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        complex(kind=8), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        end subroutine mkl_vm_vzconji_omp_offload_ilp64

        subroutine mkl_vm_vmzconji_omp_offload_ilp64(n, a, inca, y, incy, mode)
        integer(kind=8), intent(in) :: n
        complex(kind=8), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        complex(kind=8), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmzconji_omp_offload_ilp64




! function: CopySign, indexing: simple
        subroutine mkl_vm_vscopysign_omp_offload_ilp64(n, a, b, y)
        integer(kind=8), intent(in) :: n
        real(kind=4), intent(in) :: a(*)
        real(kind=4), intent(in) :: b(*)
        real(kind=4), intent(inout) :: y(*)
        end subroutine mkl_vm_vscopysign_omp_offload_ilp64

        subroutine mkl_vm_vmscopysign_omp_offload_ilp64(n, a, b, y, mode)
        integer(kind=8), intent(in) :: n
        real(kind=4), intent(in) :: a(*)
        real(kind=4), intent(in) :: b(*)
        real(kind=4), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmscopysign_omp_offload_ilp64

        subroutine mkl_vm_vdcopysign_omp_offload_ilp64(n, a, b, y)
        integer(kind=8), intent(in) :: n
        real(kind=8), intent(in) :: a(*)
        real(kind=8), intent(in) :: b(*)
        real(kind=8), intent(inout) :: y(*)
        end subroutine mkl_vm_vdcopysign_omp_offload_ilp64

        subroutine mkl_vm_vmdcopysign_omp_offload_ilp64(n, a, b, y, mode)
        integer(kind=8), intent(in) :: n
        real(kind=8), intent(in) :: a(*)
        real(kind=8), intent(in) :: b(*)
        real(kind=8), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmdcopysign_omp_offload_ilp64


! function: CopySign, indexing: strided
        subroutine mkl_vm_vscopysigni_omp_offload_ilp64(n, a, inca, b, incb, y, incy)
        integer(kind=8), intent(in) :: n
        real(kind=4), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        real(kind=4), intent(in) :: b(*)
        integer(kind=8), intent(in) :: incb
        real(kind=4), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        end subroutine mkl_vm_vscopysigni_omp_offload_ilp64

        subroutine mkl_vm_vmscopysigni_omp_offload_ilp64(n, a, inca, b, incb, y, incy, mode)
        integer(kind=8), intent(in) :: n
        real(kind=4), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        real(kind=4), intent(in) :: b(*)
        integer(kind=8), intent(in) :: incb
        real(kind=4), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmscopysigni_omp_offload_ilp64

        subroutine mkl_vm_vdcopysigni_omp_offload_ilp64(n, a, inca, b, incb, y, incy)
        integer(kind=8), intent(in) :: n
        real(kind=8), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        real(kind=8), intent(in) :: b(*)
        integer(kind=8), intent(in) :: incb
        real(kind=8), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        end subroutine mkl_vm_vdcopysigni_omp_offload_ilp64

        subroutine mkl_vm_vmdcopysigni_omp_offload_ilp64(n, a, inca, b, incb, y, incy, mode)
        integer(kind=8), intent(in) :: n
        real(kind=8), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        real(kind=8), intent(in) :: b(*)
        integer(kind=8), intent(in) :: incb
        real(kind=8), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmdcopysigni_omp_offload_ilp64




! function: Cos, indexing: simple
        subroutine mkl_vm_vscos_omp_offload_ilp64(n, a, y)
        integer(kind=8), intent(in) :: n
        real(kind=4), intent(in) :: a(*)
        real(kind=4), intent(inout) :: y(*)
        end subroutine mkl_vm_vscos_omp_offload_ilp64

        subroutine mkl_vm_vmscos_omp_offload_ilp64(n, a, y, mode)
        integer(kind=8), intent(in) :: n
        real(kind=4), intent(in) :: a(*)
        real(kind=4), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmscos_omp_offload_ilp64

        subroutine mkl_vm_vdcos_omp_offload_ilp64(n, a, y)
        integer(kind=8), intent(in) :: n
        real(kind=8), intent(in) :: a(*)
        real(kind=8), intent(inout) :: y(*)
        end subroutine mkl_vm_vdcos_omp_offload_ilp64

        subroutine mkl_vm_vmdcos_omp_offload_ilp64(n, a, y, mode)
        integer(kind=8), intent(in) :: n
        real(kind=8), intent(in) :: a(*)
        real(kind=8), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmdcos_omp_offload_ilp64

        subroutine mkl_vm_vccos_omp_offload_ilp64(n, a, y)
        integer(kind=8), intent(in) :: n
        complex(kind=4), intent(in) :: a(*)
        complex(kind=4), intent(inout) :: y(*)
        end subroutine mkl_vm_vccos_omp_offload_ilp64

        subroutine mkl_vm_vmccos_omp_offload_ilp64(n, a, y, mode)
        integer(kind=8), intent(in) :: n
        complex(kind=4), intent(in) :: a(*)
        complex(kind=4), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmccos_omp_offload_ilp64

        subroutine mkl_vm_vzcos_omp_offload_ilp64(n, a, y)
        integer(kind=8), intent(in) :: n
        complex(kind=8), intent(in) :: a(*)
        complex(kind=8), intent(inout) :: y(*)
        end subroutine mkl_vm_vzcos_omp_offload_ilp64

        subroutine mkl_vm_vmzcos_omp_offload_ilp64(n, a, y, mode)
        integer(kind=8), intent(in) :: n
        complex(kind=8), intent(in) :: a(*)
        complex(kind=8), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmzcos_omp_offload_ilp64


! function: Cos, indexing: strided
        subroutine mkl_vm_vscosi_omp_offload_ilp64(n, a, inca, y, incy)
        integer(kind=8), intent(in) :: n
        real(kind=4), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        real(kind=4), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        end subroutine mkl_vm_vscosi_omp_offload_ilp64

        subroutine mkl_vm_vmscosi_omp_offload_ilp64(n, a, inca, y, incy, mode)
        integer(kind=8), intent(in) :: n
        real(kind=4), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        real(kind=4), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmscosi_omp_offload_ilp64

        subroutine mkl_vm_vdcosi_omp_offload_ilp64(n, a, inca, y, incy)
        integer(kind=8), intent(in) :: n
        real(kind=8), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        real(kind=8), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        end subroutine mkl_vm_vdcosi_omp_offload_ilp64

        subroutine mkl_vm_vmdcosi_omp_offload_ilp64(n, a, inca, y, incy, mode)
        integer(kind=8), intent(in) :: n
        real(kind=8), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        real(kind=8), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmdcosi_omp_offload_ilp64

        subroutine mkl_vm_vccosi_omp_offload_ilp64(n, a, inca, y, incy)
        integer(kind=8), intent(in) :: n
        complex(kind=4), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        complex(kind=4), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        end subroutine mkl_vm_vccosi_omp_offload_ilp64

        subroutine mkl_vm_vmccosi_omp_offload_ilp64(n, a, inca, y, incy, mode)
        integer(kind=8), intent(in) :: n
        complex(kind=4), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        complex(kind=4), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmccosi_omp_offload_ilp64

        subroutine mkl_vm_vzcosi_omp_offload_ilp64(n, a, inca, y, incy)
        integer(kind=8), intent(in) :: n
        complex(kind=8), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        complex(kind=8), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        end subroutine mkl_vm_vzcosi_omp_offload_ilp64

        subroutine mkl_vm_vmzcosi_omp_offload_ilp64(n, a, inca, y, incy, mode)
        integer(kind=8), intent(in) :: n
        complex(kind=8), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        complex(kind=8), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmzcosi_omp_offload_ilp64




! function: Cosd, indexing: simple
        subroutine mkl_vm_vscosd_omp_offload_ilp64(n, a, y)
        integer(kind=8), intent(in) :: n
        real(kind=4), intent(in) :: a(*)
        real(kind=4), intent(inout) :: y(*)
        end subroutine mkl_vm_vscosd_omp_offload_ilp64

        subroutine mkl_vm_vmscosd_omp_offload_ilp64(n, a, y, mode)
        integer(kind=8), intent(in) :: n
        real(kind=4), intent(in) :: a(*)
        real(kind=4), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmscosd_omp_offload_ilp64

        subroutine mkl_vm_vdcosd_omp_offload_ilp64(n, a, y)
        integer(kind=8), intent(in) :: n
        real(kind=8), intent(in) :: a(*)
        real(kind=8), intent(inout) :: y(*)
        end subroutine mkl_vm_vdcosd_omp_offload_ilp64

        subroutine mkl_vm_vmdcosd_omp_offload_ilp64(n, a, y, mode)
        integer(kind=8), intent(in) :: n
        real(kind=8), intent(in) :: a(*)
        real(kind=8), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmdcosd_omp_offload_ilp64


! function: Cosd, indexing: strided
        subroutine mkl_vm_vscosdi_omp_offload_ilp64(n, a, inca, y, incy)
        integer(kind=8), intent(in) :: n
        real(kind=4), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        real(kind=4), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        end subroutine mkl_vm_vscosdi_omp_offload_ilp64

        subroutine mkl_vm_vmscosdi_omp_offload_ilp64(n, a, inca, y, incy, mode)
        integer(kind=8), intent(in) :: n
        real(kind=4), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        real(kind=4), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmscosdi_omp_offload_ilp64

        subroutine mkl_vm_vdcosdi_omp_offload_ilp64(n, a, inca, y, incy)
        integer(kind=8), intent(in) :: n
        real(kind=8), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        real(kind=8), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        end subroutine mkl_vm_vdcosdi_omp_offload_ilp64

        subroutine mkl_vm_vmdcosdi_omp_offload_ilp64(n, a, inca, y, incy, mode)
        integer(kind=8), intent(in) :: n
        real(kind=8), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        real(kind=8), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmdcosdi_omp_offload_ilp64




! function: Cosh, indexing: simple
        subroutine mkl_vm_vscosh_omp_offload_ilp64(n, a, y)
        integer(kind=8), intent(in) :: n
        real(kind=4), intent(in) :: a(*)
        real(kind=4), intent(inout) :: y(*)
        end subroutine mkl_vm_vscosh_omp_offload_ilp64

        subroutine mkl_vm_vmscosh_omp_offload_ilp64(n, a, y, mode)
        integer(kind=8), intent(in) :: n
        real(kind=4), intent(in) :: a(*)
        real(kind=4), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmscosh_omp_offload_ilp64

        subroutine mkl_vm_vdcosh_omp_offload_ilp64(n, a, y)
        integer(kind=8), intent(in) :: n
        real(kind=8), intent(in) :: a(*)
        real(kind=8), intent(inout) :: y(*)
        end subroutine mkl_vm_vdcosh_omp_offload_ilp64

        subroutine mkl_vm_vmdcosh_omp_offload_ilp64(n, a, y, mode)
        integer(kind=8), intent(in) :: n
        real(kind=8), intent(in) :: a(*)
        real(kind=8), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmdcosh_omp_offload_ilp64

        subroutine mkl_vm_vccosh_omp_offload_ilp64(n, a, y)
        integer(kind=8), intent(in) :: n
        complex(kind=4), intent(in) :: a(*)
        complex(kind=4), intent(inout) :: y(*)
        end subroutine mkl_vm_vccosh_omp_offload_ilp64

        subroutine mkl_vm_vmccosh_omp_offload_ilp64(n, a, y, mode)
        integer(kind=8), intent(in) :: n
        complex(kind=4), intent(in) :: a(*)
        complex(kind=4), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmccosh_omp_offload_ilp64

        subroutine mkl_vm_vzcosh_omp_offload_ilp64(n, a, y)
        integer(kind=8), intent(in) :: n
        complex(kind=8), intent(in) :: a(*)
        complex(kind=8), intent(inout) :: y(*)
        end subroutine mkl_vm_vzcosh_omp_offload_ilp64

        subroutine mkl_vm_vmzcosh_omp_offload_ilp64(n, a, y, mode)
        integer(kind=8), intent(in) :: n
        complex(kind=8), intent(in) :: a(*)
        complex(kind=8), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmzcosh_omp_offload_ilp64


! function: Cosh, indexing: strided
        subroutine mkl_vm_vscoshi_omp_offload_ilp64(n, a, inca, y, incy)
        integer(kind=8), intent(in) :: n
        real(kind=4), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        real(kind=4), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        end subroutine mkl_vm_vscoshi_omp_offload_ilp64

        subroutine mkl_vm_vmscoshi_omp_offload_ilp64(n, a, inca, y, incy, mode)
        integer(kind=8), intent(in) :: n
        real(kind=4), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        real(kind=4), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmscoshi_omp_offload_ilp64

        subroutine mkl_vm_vdcoshi_omp_offload_ilp64(n, a, inca, y, incy)
        integer(kind=8), intent(in) :: n
        real(kind=8), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        real(kind=8), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        end subroutine mkl_vm_vdcoshi_omp_offload_ilp64

        subroutine mkl_vm_vmdcoshi_omp_offload_ilp64(n, a, inca, y, incy, mode)
        integer(kind=8), intent(in) :: n
        real(kind=8), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        real(kind=8), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmdcoshi_omp_offload_ilp64

        subroutine mkl_vm_vccoshi_omp_offload_ilp64(n, a, inca, y, incy)
        integer(kind=8), intent(in) :: n
        complex(kind=4), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        complex(kind=4), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        end subroutine mkl_vm_vccoshi_omp_offload_ilp64

        subroutine mkl_vm_vmccoshi_omp_offload_ilp64(n, a, inca, y, incy, mode)
        integer(kind=8), intent(in) :: n
        complex(kind=4), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        complex(kind=4), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmccoshi_omp_offload_ilp64

        subroutine mkl_vm_vzcoshi_omp_offload_ilp64(n, a, inca, y, incy)
        integer(kind=8), intent(in) :: n
        complex(kind=8), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        complex(kind=8), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        end subroutine mkl_vm_vzcoshi_omp_offload_ilp64

        subroutine mkl_vm_vmzcoshi_omp_offload_ilp64(n, a, inca, y, incy, mode)
        integer(kind=8), intent(in) :: n
        complex(kind=8), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        complex(kind=8), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmzcoshi_omp_offload_ilp64




! function: Cospi, indexing: simple
        subroutine mkl_vm_vscospi_omp_offload_ilp64(n, a, y)
        integer(kind=8), intent(in) :: n
        real(kind=4), intent(in) :: a(*)
        real(kind=4), intent(inout) :: y(*)
        end subroutine mkl_vm_vscospi_omp_offload_ilp64

        subroutine mkl_vm_vmscospi_omp_offload_ilp64(n, a, y, mode)
        integer(kind=8), intent(in) :: n
        real(kind=4), intent(in) :: a(*)
        real(kind=4), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmscospi_omp_offload_ilp64

        subroutine mkl_vm_vdcospi_omp_offload_ilp64(n, a, y)
        integer(kind=8), intent(in) :: n
        real(kind=8), intent(in) :: a(*)
        real(kind=8), intent(inout) :: y(*)
        end subroutine mkl_vm_vdcospi_omp_offload_ilp64

        subroutine mkl_vm_vmdcospi_omp_offload_ilp64(n, a, y, mode)
        integer(kind=8), intent(in) :: n
        real(kind=8), intent(in) :: a(*)
        real(kind=8), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmdcospi_omp_offload_ilp64


! function: Cospi, indexing: strided
        subroutine mkl_vm_vscospii_omp_offload_ilp64(n, a, inca, y, incy)
        integer(kind=8), intent(in) :: n
        real(kind=4), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        real(kind=4), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        end subroutine mkl_vm_vscospii_omp_offload_ilp64

        subroutine mkl_vm_vmscospii_omp_offload_ilp64(n, a, inca, y, incy, mode)
        integer(kind=8), intent(in) :: n
        real(kind=4), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        real(kind=4), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmscospii_omp_offload_ilp64

        subroutine mkl_vm_vdcospii_omp_offload_ilp64(n, a, inca, y, incy)
        integer(kind=8), intent(in) :: n
        real(kind=8), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        real(kind=8), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        end subroutine mkl_vm_vdcospii_omp_offload_ilp64

        subroutine mkl_vm_vmdcospii_omp_offload_ilp64(n, a, inca, y, incy, mode)
        integer(kind=8), intent(in) :: n
        real(kind=8), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        real(kind=8), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmdcospii_omp_offload_ilp64




! function: Div, indexing: simple
        subroutine mkl_vm_vsdiv_omp_offload_ilp64(n, a, b, y)
        integer(kind=8), intent(in) :: n
        real(kind=4), intent(in) :: a(*)
        real(kind=4), intent(in) :: b(*)
        real(kind=4), intent(inout) :: y(*)
        end subroutine mkl_vm_vsdiv_omp_offload_ilp64

        subroutine mkl_vm_vmsdiv_omp_offload_ilp64(n, a, b, y, mode)
        integer(kind=8), intent(in) :: n
        real(kind=4), intent(in) :: a(*)
        real(kind=4), intent(in) :: b(*)
        real(kind=4), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmsdiv_omp_offload_ilp64

        subroutine mkl_vm_vddiv_omp_offload_ilp64(n, a, b, y)
        integer(kind=8), intent(in) :: n
        real(kind=8), intent(in) :: a(*)
        real(kind=8), intent(in) :: b(*)
        real(kind=8), intent(inout) :: y(*)
        end subroutine mkl_vm_vddiv_omp_offload_ilp64

        subroutine mkl_vm_vmddiv_omp_offload_ilp64(n, a, b, y, mode)
        integer(kind=8), intent(in) :: n
        real(kind=8), intent(in) :: a(*)
        real(kind=8), intent(in) :: b(*)
        real(kind=8), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmddiv_omp_offload_ilp64

        subroutine mkl_vm_vcdiv_omp_offload_ilp64(n, a, b, y)
        integer(kind=8), intent(in) :: n
        complex(kind=4), intent(in) :: a(*)
        complex(kind=4), intent(in) :: b(*)
        complex(kind=4), intent(inout) :: y(*)
        end subroutine mkl_vm_vcdiv_omp_offload_ilp64

        subroutine mkl_vm_vmcdiv_omp_offload_ilp64(n, a, b, y, mode)
        integer(kind=8), intent(in) :: n
        complex(kind=4), intent(in) :: a(*)
        complex(kind=4), intent(in) :: b(*)
        complex(kind=4), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmcdiv_omp_offload_ilp64

        subroutine mkl_vm_vzdiv_omp_offload_ilp64(n, a, b, y)
        integer(kind=8), intent(in) :: n
        complex(kind=8), intent(in) :: a(*)
        complex(kind=8), intent(in) :: b(*)
        complex(kind=8), intent(inout) :: y(*)
        end subroutine mkl_vm_vzdiv_omp_offload_ilp64

        subroutine mkl_vm_vmzdiv_omp_offload_ilp64(n, a, b, y, mode)
        integer(kind=8), intent(in) :: n
        complex(kind=8), intent(in) :: a(*)
        complex(kind=8), intent(in) :: b(*)
        complex(kind=8), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmzdiv_omp_offload_ilp64


! function: Div, indexing: strided
        subroutine mkl_vm_vsdivi_omp_offload_ilp64(n, a, inca, b, incb, y, incy)
        integer(kind=8), intent(in) :: n
        real(kind=4), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        real(kind=4), intent(in) :: b(*)
        integer(kind=8), intent(in) :: incb
        real(kind=4), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        end subroutine mkl_vm_vsdivi_omp_offload_ilp64

        subroutine mkl_vm_vmsdivi_omp_offload_ilp64(n, a, inca, b, incb, y, incy, mode)
        integer(kind=8), intent(in) :: n
        real(kind=4), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        real(kind=4), intent(in) :: b(*)
        integer(kind=8), intent(in) :: incb
        real(kind=4), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmsdivi_omp_offload_ilp64

        subroutine mkl_vm_vddivi_omp_offload_ilp64(n, a, inca, b, incb, y, incy)
        integer(kind=8), intent(in) :: n
        real(kind=8), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        real(kind=8), intent(in) :: b(*)
        integer(kind=8), intent(in) :: incb
        real(kind=8), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        end subroutine mkl_vm_vddivi_omp_offload_ilp64

        subroutine mkl_vm_vmddivi_omp_offload_ilp64(n, a, inca, b, incb, y, incy, mode)
        integer(kind=8), intent(in) :: n
        real(kind=8), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        real(kind=8), intent(in) :: b(*)
        integer(kind=8), intent(in) :: incb
        real(kind=8), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmddivi_omp_offload_ilp64

        subroutine mkl_vm_vcdivi_omp_offload_ilp64(n, a, inca, b, incb, y, incy)
        integer(kind=8), intent(in) :: n
        complex(kind=4), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        complex(kind=4), intent(in) :: b(*)
        integer(kind=8), intent(in) :: incb
        complex(kind=4), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        end subroutine mkl_vm_vcdivi_omp_offload_ilp64

        subroutine mkl_vm_vmcdivi_omp_offload_ilp64(n, a, inca, b, incb, y, incy, mode)
        integer(kind=8), intent(in) :: n
        complex(kind=4), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        complex(kind=4), intent(in) :: b(*)
        integer(kind=8), intent(in) :: incb
        complex(kind=4), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmcdivi_omp_offload_ilp64

        subroutine mkl_vm_vzdivi_omp_offload_ilp64(n, a, inca, b, incb, y, incy)
        integer(kind=8), intent(in) :: n
        complex(kind=8), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        complex(kind=8), intent(in) :: b(*)
        integer(kind=8), intent(in) :: incb
        complex(kind=8), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        end subroutine mkl_vm_vzdivi_omp_offload_ilp64

        subroutine mkl_vm_vmzdivi_omp_offload_ilp64(n, a, inca, b, incb, y, incy, mode)
        integer(kind=8), intent(in) :: n
        complex(kind=8), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        complex(kind=8), intent(in) :: b(*)
        integer(kind=8), intent(in) :: incb
        complex(kind=8), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmzdivi_omp_offload_ilp64




! function: Erf, indexing: simple
        subroutine mkl_vm_vserf_omp_offload_ilp64(n, a, y)
        integer(kind=8), intent(in) :: n
        real(kind=4), intent(in) :: a(*)
        real(kind=4), intent(inout) :: y(*)
        end subroutine mkl_vm_vserf_omp_offload_ilp64

        subroutine mkl_vm_vmserf_omp_offload_ilp64(n, a, y, mode)
        integer(kind=8), intent(in) :: n
        real(kind=4), intent(in) :: a(*)
        real(kind=4), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmserf_omp_offload_ilp64

        subroutine mkl_vm_vderf_omp_offload_ilp64(n, a, y)
        integer(kind=8), intent(in) :: n
        real(kind=8), intent(in) :: a(*)
        real(kind=8), intent(inout) :: y(*)
        end subroutine mkl_vm_vderf_omp_offload_ilp64

        subroutine mkl_vm_vmderf_omp_offload_ilp64(n, a, y, mode)
        integer(kind=8), intent(in) :: n
        real(kind=8), intent(in) :: a(*)
        real(kind=8), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmderf_omp_offload_ilp64


! function: Erf, indexing: strided
        subroutine mkl_vm_vserfi_omp_offload_ilp64(n, a, inca, y, incy)
        integer(kind=8), intent(in) :: n
        real(kind=4), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        real(kind=4), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        end subroutine mkl_vm_vserfi_omp_offload_ilp64

        subroutine mkl_vm_vmserfi_omp_offload_ilp64(n, a, inca, y, incy, mode)
        integer(kind=8), intent(in) :: n
        real(kind=4), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        real(kind=4), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmserfi_omp_offload_ilp64

        subroutine mkl_vm_vderfi_omp_offload_ilp64(n, a, inca, y, incy)
        integer(kind=8), intent(in) :: n
        real(kind=8), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        real(kind=8), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        end subroutine mkl_vm_vderfi_omp_offload_ilp64

        subroutine mkl_vm_vmderfi_omp_offload_ilp64(n, a, inca, y, incy, mode)
        integer(kind=8), intent(in) :: n
        real(kind=8), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        real(kind=8), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmderfi_omp_offload_ilp64




! function: Erfc, indexing: simple
        subroutine mkl_vm_vserfc_omp_offload_ilp64(n, a, y)
        integer(kind=8), intent(in) :: n
        real(kind=4), intent(in) :: a(*)
        real(kind=4), intent(inout) :: y(*)
        end subroutine mkl_vm_vserfc_omp_offload_ilp64

        subroutine mkl_vm_vmserfc_omp_offload_ilp64(n, a, y, mode)
        integer(kind=8), intent(in) :: n
        real(kind=4), intent(in) :: a(*)
        real(kind=4), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmserfc_omp_offload_ilp64

        subroutine mkl_vm_vderfc_omp_offload_ilp64(n, a, y)
        integer(kind=8), intent(in) :: n
        real(kind=8), intent(in) :: a(*)
        real(kind=8), intent(inout) :: y(*)
        end subroutine mkl_vm_vderfc_omp_offload_ilp64

        subroutine mkl_vm_vmderfc_omp_offload_ilp64(n, a, y, mode)
        integer(kind=8), intent(in) :: n
        real(kind=8), intent(in) :: a(*)
        real(kind=8), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmderfc_omp_offload_ilp64


! function: Erfc, indexing: strided
        subroutine mkl_vm_vserfci_omp_offload_ilp64(n, a, inca, y, incy)
        integer(kind=8), intent(in) :: n
        real(kind=4), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        real(kind=4), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        end subroutine mkl_vm_vserfci_omp_offload_ilp64

        subroutine mkl_vm_vmserfci_omp_offload_ilp64(n, a, inca, y, incy, mode)
        integer(kind=8), intent(in) :: n
        real(kind=4), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        real(kind=4), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmserfci_omp_offload_ilp64

        subroutine mkl_vm_vderfci_omp_offload_ilp64(n, a, inca, y, incy)
        integer(kind=8), intent(in) :: n
        real(kind=8), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        real(kind=8), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        end subroutine mkl_vm_vderfci_omp_offload_ilp64

        subroutine mkl_vm_vmderfci_omp_offload_ilp64(n, a, inca, y, incy, mode)
        integer(kind=8), intent(in) :: n
        real(kind=8), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        real(kind=8), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmderfci_omp_offload_ilp64




! function: ErfcInv, indexing: simple
        subroutine mkl_vm_vserfcinv_omp_offload_ilp64(n, a, y)
        integer(kind=8), intent(in) :: n
        real(kind=4), intent(in) :: a(*)
        real(kind=4), intent(inout) :: y(*)
        end subroutine mkl_vm_vserfcinv_omp_offload_ilp64

        subroutine mkl_vm_vmserfcinv_omp_offload_ilp64(n, a, y, mode)
        integer(kind=8), intent(in) :: n
        real(kind=4), intent(in) :: a(*)
        real(kind=4), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmserfcinv_omp_offload_ilp64

        subroutine mkl_vm_vderfcinv_omp_offload_ilp64(n, a, y)
        integer(kind=8), intent(in) :: n
        real(kind=8), intent(in) :: a(*)
        real(kind=8), intent(inout) :: y(*)
        end subroutine mkl_vm_vderfcinv_omp_offload_ilp64

        subroutine mkl_vm_vmderfcinv_omp_offload_ilp64(n, a, y, mode)
        integer(kind=8), intent(in) :: n
        real(kind=8), intent(in) :: a(*)
        real(kind=8), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmderfcinv_omp_offload_ilp64


! function: ErfcInv, indexing: strided
        subroutine mkl_vm_vserfcinvi_omp_offload_ilp64(n, a, inca, y, incy)
        integer(kind=8), intent(in) :: n
        real(kind=4), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        real(kind=4), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        end subroutine mkl_vm_vserfcinvi_omp_offload_ilp64

        subroutine mkl_vm_vmserfcinvi_omp_offload_ilp64(n, a, inca, y, incy, mode)
        integer(kind=8), intent(in) :: n
        real(kind=4), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        real(kind=4), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmserfcinvi_omp_offload_ilp64

        subroutine mkl_vm_vderfcinvi_omp_offload_ilp64(n, a, inca, y, incy)
        integer(kind=8), intent(in) :: n
        real(kind=8), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        real(kind=8), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        end subroutine mkl_vm_vderfcinvi_omp_offload_ilp64

        subroutine mkl_vm_vmderfcinvi_omp_offload_ilp64(n, a, inca, y, incy, mode)
        integer(kind=8), intent(in) :: n
        real(kind=8), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        real(kind=8), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmderfcinvi_omp_offload_ilp64




! function: ErfInv, indexing: simple
        subroutine mkl_vm_vserfinv_omp_offload_ilp64(n, a, y)
        integer(kind=8), intent(in) :: n
        real(kind=4), intent(in) :: a(*)
        real(kind=4), intent(inout) :: y(*)
        end subroutine mkl_vm_vserfinv_omp_offload_ilp64

        subroutine mkl_vm_vmserfinv_omp_offload_ilp64(n, a, y, mode)
        integer(kind=8), intent(in) :: n
        real(kind=4), intent(in) :: a(*)
        real(kind=4), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmserfinv_omp_offload_ilp64

        subroutine mkl_vm_vderfinv_omp_offload_ilp64(n, a, y)
        integer(kind=8), intent(in) :: n
        real(kind=8), intent(in) :: a(*)
        real(kind=8), intent(inout) :: y(*)
        end subroutine mkl_vm_vderfinv_omp_offload_ilp64

        subroutine mkl_vm_vmderfinv_omp_offload_ilp64(n, a, y, mode)
        integer(kind=8), intent(in) :: n
        real(kind=8), intent(in) :: a(*)
        real(kind=8), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmderfinv_omp_offload_ilp64


! function: ErfInv, indexing: strided
        subroutine mkl_vm_vserfinvi_omp_offload_ilp64(n, a, inca, y, incy)
        integer(kind=8), intent(in) :: n
        real(kind=4), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        real(kind=4), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        end subroutine mkl_vm_vserfinvi_omp_offload_ilp64

        subroutine mkl_vm_vmserfinvi_omp_offload_ilp64(n, a, inca, y, incy, mode)
        integer(kind=8), intent(in) :: n
        real(kind=4), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        real(kind=4), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmserfinvi_omp_offload_ilp64

        subroutine mkl_vm_vderfinvi_omp_offload_ilp64(n, a, inca, y, incy)
        integer(kind=8), intent(in) :: n
        real(kind=8), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        real(kind=8), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        end subroutine mkl_vm_vderfinvi_omp_offload_ilp64

        subroutine mkl_vm_vmderfinvi_omp_offload_ilp64(n, a, inca, y, incy, mode)
        integer(kind=8), intent(in) :: n
        real(kind=8), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        real(kind=8), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmderfinvi_omp_offload_ilp64




! function: Exp, indexing: simple
        subroutine mkl_vm_vsexp_omp_offload_ilp64(n, a, y)
        integer(kind=8), intent(in) :: n
        real(kind=4), intent(in) :: a(*)
        real(kind=4), intent(inout) :: y(*)
        end subroutine mkl_vm_vsexp_omp_offload_ilp64

        subroutine mkl_vm_vmsexp_omp_offload_ilp64(n, a, y, mode)
        integer(kind=8), intent(in) :: n
        real(kind=4), intent(in) :: a(*)
        real(kind=4), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmsexp_omp_offload_ilp64

        subroutine mkl_vm_vdexp_omp_offload_ilp64(n, a, y)
        integer(kind=8), intent(in) :: n
        real(kind=8), intent(in) :: a(*)
        real(kind=8), intent(inout) :: y(*)
        end subroutine mkl_vm_vdexp_omp_offload_ilp64

        subroutine mkl_vm_vmdexp_omp_offload_ilp64(n, a, y, mode)
        integer(kind=8), intent(in) :: n
        real(kind=8), intent(in) :: a(*)
        real(kind=8), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmdexp_omp_offload_ilp64

        subroutine mkl_vm_vcexp_omp_offload_ilp64(n, a, y)
        integer(kind=8), intent(in) :: n
        complex(kind=4), intent(in) :: a(*)
        complex(kind=4), intent(inout) :: y(*)
        end subroutine mkl_vm_vcexp_omp_offload_ilp64

        subroutine mkl_vm_vmcexp_omp_offload_ilp64(n, a, y, mode)
        integer(kind=8), intent(in) :: n
        complex(kind=4), intent(in) :: a(*)
        complex(kind=4), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmcexp_omp_offload_ilp64

        subroutine mkl_vm_vzexp_omp_offload_ilp64(n, a, y)
        integer(kind=8), intent(in) :: n
        complex(kind=8), intent(in) :: a(*)
        complex(kind=8), intent(inout) :: y(*)
        end subroutine mkl_vm_vzexp_omp_offload_ilp64

        subroutine mkl_vm_vmzexp_omp_offload_ilp64(n, a, y, mode)
        integer(kind=8), intent(in) :: n
        complex(kind=8), intent(in) :: a(*)
        complex(kind=8), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmzexp_omp_offload_ilp64


! function: Exp, indexing: strided
        subroutine mkl_vm_vsexpi_omp_offload_ilp64(n, a, inca, y, incy)
        integer(kind=8), intent(in) :: n
        real(kind=4), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        real(kind=4), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        end subroutine mkl_vm_vsexpi_omp_offload_ilp64

        subroutine mkl_vm_vmsexpi_omp_offload_ilp64(n, a, inca, y, incy, mode)
        integer(kind=8), intent(in) :: n
        real(kind=4), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        real(kind=4), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmsexpi_omp_offload_ilp64

        subroutine mkl_vm_vdexpi_omp_offload_ilp64(n, a, inca, y, incy)
        integer(kind=8), intent(in) :: n
        real(kind=8), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        real(kind=8), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        end subroutine mkl_vm_vdexpi_omp_offload_ilp64

        subroutine mkl_vm_vmdexpi_omp_offload_ilp64(n, a, inca, y, incy, mode)
        integer(kind=8), intent(in) :: n
        real(kind=8), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        real(kind=8), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmdexpi_omp_offload_ilp64

        subroutine mkl_vm_vcexpi_omp_offload_ilp64(n, a, inca, y, incy)
        integer(kind=8), intent(in) :: n
        complex(kind=4), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        complex(kind=4), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        end subroutine mkl_vm_vcexpi_omp_offload_ilp64

        subroutine mkl_vm_vmcexpi_omp_offload_ilp64(n, a, inca, y, incy, mode)
        integer(kind=8), intent(in) :: n
        complex(kind=4), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        complex(kind=4), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmcexpi_omp_offload_ilp64

        subroutine mkl_vm_vzexpi_omp_offload_ilp64(n, a, inca, y, incy)
        integer(kind=8), intent(in) :: n
        complex(kind=8), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        complex(kind=8), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        end subroutine mkl_vm_vzexpi_omp_offload_ilp64

        subroutine mkl_vm_vmzexpi_omp_offload_ilp64(n, a, inca, y, incy, mode)
        integer(kind=8), intent(in) :: n
        complex(kind=8), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        complex(kind=8), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmzexpi_omp_offload_ilp64




! function: Exp10, indexing: simple
        subroutine mkl_vm_vsexp10_omp_offload_ilp64(n, a, y)
        integer(kind=8), intent(in) :: n
        real(kind=4), intent(in) :: a(*)
        real(kind=4), intent(inout) :: y(*)
        end subroutine mkl_vm_vsexp10_omp_offload_ilp64

        subroutine mkl_vm_vmsexp10_omp_offload_ilp64(n, a, y, mode)
        integer(kind=8), intent(in) :: n
        real(kind=4), intent(in) :: a(*)
        real(kind=4), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmsexp10_omp_offload_ilp64

        subroutine mkl_vm_vdexp10_omp_offload_ilp64(n, a, y)
        integer(kind=8), intent(in) :: n
        real(kind=8), intent(in) :: a(*)
        real(kind=8), intent(inout) :: y(*)
        end subroutine mkl_vm_vdexp10_omp_offload_ilp64

        subroutine mkl_vm_vmdexp10_omp_offload_ilp64(n, a, y, mode)
        integer(kind=8), intent(in) :: n
        real(kind=8), intent(in) :: a(*)
        real(kind=8), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmdexp10_omp_offload_ilp64


! function: Exp10, indexing: strided
        subroutine mkl_vm_vsexp10i_omp_offload_ilp64(n, a, inca, y, incy)
        integer(kind=8), intent(in) :: n
        real(kind=4), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        real(kind=4), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        end subroutine mkl_vm_vsexp10i_omp_offload_ilp64

        subroutine mkl_vm_vmsexp10i_omp_offload_ilp64(n, a, inca, y, incy, mode)
        integer(kind=8), intent(in) :: n
        real(kind=4), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        real(kind=4), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmsexp10i_omp_offload_ilp64

        subroutine mkl_vm_vdexp10i_omp_offload_ilp64(n, a, inca, y, incy)
        integer(kind=8), intent(in) :: n
        real(kind=8), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        real(kind=8), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        end subroutine mkl_vm_vdexp10i_omp_offload_ilp64

        subroutine mkl_vm_vmdexp10i_omp_offload_ilp64(n, a, inca, y, incy, mode)
        integer(kind=8), intent(in) :: n
        real(kind=8), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        real(kind=8), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmdexp10i_omp_offload_ilp64




! function: Exp2, indexing: simple
        subroutine mkl_vm_vsexp2_omp_offload_ilp64(n, a, y)
        integer(kind=8), intent(in) :: n
        real(kind=4), intent(in) :: a(*)
        real(kind=4), intent(inout) :: y(*)
        end subroutine mkl_vm_vsexp2_omp_offload_ilp64

        subroutine mkl_vm_vmsexp2_omp_offload_ilp64(n, a, y, mode)
        integer(kind=8), intent(in) :: n
        real(kind=4), intent(in) :: a(*)
        real(kind=4), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmsexp2_omp_offload_ilp64

        subroutine mkl_vm_vdexp2_omp_offload_ilp64(n, a, y)
        integer(kind=8), intent(in) :: n
        real(kind=8), intent(in) :: a(*)
        real(kind=8), intent(inout) :: y(*)
        end subroutine mkl_vm_vdexp2_omp_offload_ilp64

        subroutine mkl_vm_vmdexp2_omp_offload_ilp64(n, a, y, mode)
        integer(kind=8), intent(in) :: n
        real(kind=8), intent(in) :: a(*)
        real(kind=8), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmdexp2_omp_offload_ilp64


! function: Exp2, indexing: strided
        subroutine mkl_vm_vsexp2i_omp_offload_ilp64(n, a, inca, y, incy)
        integer(kind=8), intent(in) :: n
        real(kind=4), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        real(kind=4), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        end subroutine mkl_vm_vsexp2i_omp_offload_ilp64

        subroutine mkl_vm_vmsexp2i_omp_offload_ilp64(n, a, inca, y, incy, mode)
        integer(kind=8), intent(in) :: n
        real(kind=4), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        real(kind=4), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmsexp2i_omp_offload_ilp64

        subroutine mkl_vm_vdexp2i_omp_offload_ilp64(n, a, inca, y, incy)
        integer(kind=8), intent(in) :: n
        real(kind=8), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        real(kind=8), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        end subroutine mkl_vm_vdexp2i_omp_offload_ilp64

        subroutine mkl_vm_vmdexp2i_omp_offload_ilp64(n, a, inca, y, incy, mode)
        integer(kind=8), intent(in) :: n
        real(kind=8), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        real(kind=8), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmdexp2i_omp_offload_ilp64




! function: ExpInt1, indexing: simple
        subroutine mkl_vm_vsexpint1_omp_offload_ilp64(n, a, y)
        integer(kind=8), intent(in) :: n
        real(kind=4), intent(in) :: a(*)
        real(kind=4), intent(inout) :: y(*)
        end subroutine mkl_vm_vsexpint1_omp_offload_ilp64

        subroutine mkl_vm_vmsexpint1_omp_offload_ilp64(n, a, y, mode)
        integer(kind=8), intent(in) :: n
        real(kind=4), intent(in) :: a(*)
        real(kind=4), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmsexpint1_omp_offload_ilp64

        subroutine mkl_vm_vdexpint1_omp_offload_ilp64(n, a, y)
        integer(kind=8), intent(in) :: n
        real(kind=8), intent(in) :: a(*)
        real(kind=8), intent(inout) :: y(*)
        end subroutine mkl_vm_vdexpint1_omp_offload_ilp64

        subroutine mkl_vm_vmdexpint1_omp_offload_ilp64(n, a, y, mode)
        integer(kind=8), intent(in) :: n
        real(kind=8), intent(in) :: a(*)
        real(kind=8), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmdexpint1_omp_offload_ilp64


! function: ExpInt1, indexing: strided
        subroutine mkl_vm_vsexpint1i_omp_offload_ilp64(n, a, inca, y, incy)
        integer(kind=8), intent(in) :: n
        real(kind=4), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        real(kind=4), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        end subroutine mkl_vm_vsexpint1i_omp_offload_ilp64

        subroutine mkl_vm_vmsexpint1i_omp_offload_ilp64(n, a, inca, y, incy, mode)
        integer(kind=8), intent(in) :: n
        real(kind=4), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        real(kind=4), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmsexpint1i_omp_offload_ilp64

        subroutine mkl_vm_vdexpint1i_omp_offload_ilp64(n, a, inca, y, incy)
        integer(kind=8), intent(in) :: n
        real(kind=8), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        real(kind=8), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        end subroutine mkl_vm_vdexpint1i_omp_offload_ilp64

        subroutine mkl_vm_vmdexpint1i_omp_offload_ilp64(n, a, inca, y, incy, mode)
        integer(kind=8), intent(in) :: n
        real(kind=8), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        real(kind=8), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmdexpint1i_omp_offload_ilp64




! function: Expm1, indexing: simple
        subroutine mkl_vm_vsexpm1_omp_offload_ilp64(n, a, y)
        integer(kind=8), intent(in) :: n
        real(kind=4), intent(in) :: a(*)
        real(kind=4), intent(inout) :: y(*)
        end subroutine mkl_vm_vsexpm1_omp_offload_ilp64

        subroutine mkl_vm_vmsexpm1_omp_offload_ilp64(n, a, y, mode)
        integer(kind=8), intent(in) :: n
        real(kind=4), intent(in) :: a(*)
        real(kind=4), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmsexpm1_omp_offload_ilp64

        subroutine mkl_vm_vdexpm1_omp_offload_ilp64(n, a, y)
        integer(kind=8), intent(in) :: n
        real(kind=8), intent(in) :: a(*)
        real(kind=8), intent(inout) :: y(*)
        end subroutine mkl_vm_vdexpm1_omp_offload_ilp64

        subroutine mkl_vm_vmdexpm1_omp_offload_ilp64(n, a, y, mode)
        integer(kind=8), intent(in) :: n
        real(kind=8), intent(in) :: a(*)
        real(kind=8), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmdexpm1_omp_offload_ilp64


! function: Expm1, indexing: strided
        subroutine mkl_vm_vsexpm1i_omp_offload_ilp64(n, a, inca, y, incy)
        integer(kind=8), intent(in) :: n
        real(kind=4), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        real(kind=4), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        end subroutine mkl_vm_vsexpm1i_omp_offload_ilp64

        subroutine mkl_vm_vmsexpm1i_omp_offload_ilp64(n, a, inca, y, incy, mode)
        integer(kind=8), intent(in) :: n
        real(kind=4), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        real(kind=4), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmsexpm1i_omp_offload_ilp64

        subroutine mkl_vm_vdexpm1i_omp_offload_ilp64(n, a, inca, y, incy)
        integer(kind=8), intent(in) :: n
        real(kind=8), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        real(kind=8), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        end subroutine mkl_vm_vdexpm1i_omp_offload_ilp64

        subroutine mkl_vm_vmdexpm1i_omp_offload_ilp64(n, a, inca, y, incy, mode)
        integer(kind=8), intent(in) :: n
        real(kind=8), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        real(kind=8), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmdexpm1i_omp_offload_ilp64




! function: Fdim, indexing: simple
        subroutine mkl_vm_vsfdim_omp_offload_ilp64(n, a, b, y)
        integer(kind=8), intent(in) :: n
        real(kind=4), intent(in) :: a(*)
        real(kind=4), intent(in) :: b(*)
        real(kind=4), intent(inout) :: y(*)
        end subroutine mkl_vm_vsfdim_omp_offload_ilp64

        subroutine mkl_vm_vmsfdim_omp_offload_ilp64(n, a, b, y, mode)
        integer(kind=8), intent(in) :: n
        real(kind=4), intent(in) :: a(*)
        real(kind=4), intent(in) :: b(*)
        real(kind=4), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmsfdim_omp_offload_ilp64

        subroutine mkl_vm_vdfdim_omp_offload_ilp64(n, a, b, y)
        integer(kind=8), intent(in) :: n
        real(kind=8), intent(in) :: a(*)
        real(kind=8), intent(in) :: b(*)
        real(kind=8), intent(inout) :: y(*)
        end subroutine mkl_vm_vdfdim_omp_offload_ilp64

        subroutine mkl_vm_vmdfdim_omp_offload_ilp64(n, a, b, y, mode)
        integer(kind=8), intent(in) :: n
        real(kind=8), intent(in) :: a(*)
        real(kind=8), intent(in) :: b(*)
        real(kind=8), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmdfdim_omp_offload_ilp64


! function: Fdim, indexing: strided
        subroutine mkl_vm_vsfdimi_omp_offload_ilp64(n, a, inca, b, incb, y, incy)
        integer(kind=8), intent(in) :: n
        real(kind=4), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        real(kind=4), intent(in) :: b(*)
        integer(kind=8), intent(in) :: incb
        real(kind=4), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        end subroutine mkl_vm_vsfdimi_omp_offload_ilp64

        subroutine mkl_vm_vmsfdimi_omp_offload_ilp64(n, a, inca, b, incb, y, incy, mode)
        integer(kind=8), intent(in) :: n
        real(kind=4), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        real(kind=4), intent(in) :: b(*)
        integer(kind=8), intent(in) :: incb
        real(kind=4), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmsfdimi_omp_offload_ilp64

        subroutine mkl_vm_vdfdimi_omp_offload_ilp64(n, a, inca, b, incb, y, incy)
        integer(kind=8), intent(in) :: n
        real(kind=8), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        real(kind=8), intent(in) :: b(*)
        integer(kind=8), intent(in) :: incb
        real(kind=8), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        end subroutine mkl_vm_vdfdimi_omp_offload_ilp64

        subroutine mkl_vm_vmdfdimi_omp_offload_ilp64(n, a, inca, b, incb, y, incy, mode)
        integer(kind=8), intent(in) :: n
        real(kind=8), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        real(kind=8), intent(in) :: b(*)
        integer(kind=8), intent(in) :: incb
        real(kind=8), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmdfdimi_omp_offload_ilp64




! function: Floor, indexing: simple
        subroutine mkl_vm_vsfloor_omp_offload_ilp64(n, a, y)
        integer(kind=8), intent(in) :: n
        real(kind=4), intent(in) :: a(*)
        real(kind=4), intent(inout) :: y(*)
        end subroutine mkl_vm_vsfloor_omp_offload_ilp64

        subroutine mkl_vm_vmsfloor_omp_offload_ilp64(n, a, y, mode)
        integer(kind=8), intent(in) :: n
        real(kind=4), intent(in) :: a(*)
        real(kind=4), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmsfloor_omp_offload_ilp64

        subroutine mkl_vm_vdfloor_omp_offload_ilp64(n, a, y)
        integer(kind=8), intent(in) :: n
        real(kind=8), intent(in) :: a(*)
        real(kind=8), intent(inout) :: y(*)
        end subroutine mkl_vm_vdfloor_omp_offload_ilp64

        subroutine mkl_vm_vmdfloor_omp_offload_ilp64(n, a, y, mode)
        integer(kind=8), intent(in) :: n
        real(kind=8), intent(in) :: a(*)
        real(kind=8), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmdfloor_omp_offload_ilp64


! function: Floor, indexing: strided
        subroutine mkl_vm_vsfloori_omp_offload_ilp64(n, a, inca, y, incy)
        integer(kind=8), intent(in) :: n
        real(kind=4), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        real(kind=4), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        end subroutine mkl_vm_vsfloori_omp_offload_ilp64

        subroutine mkl_vm_vmsfloori_omp_offload_ilp64(n, a, inca, y, incy, mode)
        integer(kind=8), intent(in) :: n
        real(kind=4), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        real(kind=4), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmsfloori_omp_offload_ilp64

        subroutine mkl_vm_vdfloori_omp_offload_ilp64(n, a, inca, y, incy)
        integer(kind=8), intent(in) :: n
        real(kind=8), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        real(kind=8), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        end subroutine mkl_vm_vdfloori_omp_offload_ilp64

        subroutine mkl_vm_vmdfloori_omp_offload_ilp64(n, a, inca, y, incy, mode)
        integer(kind=8), intent(in) :: n
        real(kind=8), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        real(kind=8), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmdfloori_omp_offload_ilp64




! function: Fmax, indexing: simple
        subroutine mkl_vm_vsfmax_omp_offload_ilp64(n, a, b, y)
        integer(kind=8), intent(in) :: n
        real(kind=4), intent(in) :: a(*)
        real(kind=4), intent(in) :: b(*)
        real(kind=4), intent(inout) :: y(*)
        end subroutine mkl_vm_vsfmax_omp_offload_ilp64

        subroutine mkl_vm_vmsfmax_omp_offload_ilp64(n, a, b, y, mode)
        integer(kind=8), intent(in) :: n
        real(kind=4), intent(in) :: a(*)
        real(kind=4), intent(in) :: b(*)
        real(kind=4), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmsfmax_omp_offload_ilp64

        subroutine mkl_vm_vdfmax_omp_offload_ilp64(n, a, b, y)
        integer(kind=8), intent(in) :: n
        real(kind=8), intent(in) :: a(*)
        real(kind=8), intent(in) :: b(*)
        real(kind=8), intent(inout) :: y(*)
        end subroutine mkl_vm_vdfmax_omp_offload_ilp64

        subroutine mkl_vm_vmdfmax_omp_offload_ilp64(n, a, b, y, mode)
        integer(kind=8), intent(in) :: n
        real(kind=8), intent(in) :: a(*)
        real(kind=8), intent(in) :: b(*)
        real(kind=8), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmdfmax_omp_offload_ilp64


! function: Fmax, indexing: strided
        subroutine mkl_vm_vsfmaxi_omp_offload_ilp64(n, a, inca, b, incb, y, incy)
        integer(kind=8), intent(in) :: n
        real(kind=4), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        real(kind=4), intent(in) :: b(*)
        integer(kind=8), intent(in) :: incb
        real(kind=4), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        end subroutine mkl_vm_vsfmaxi_omp_offload_ilp64

        subroutine mkl_vm_vmsfmaxi_omp_offload_ilp64(n, a, inca, b, incb, y, incy, mode)
        integer(kind=8), intent(in) :: n
        real(kind=4), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        real(kind=4), intent(in) :: b(*)
        integer(kind=8), intent(in) :: incb
        real(kind=4), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmsfmaxi_omp_offload_ilp64

        subroutine mkl_vm_vdfmaxi_omp_offload_ilp64(n, a, inca, b, incb, y, incy)
        integer(kind=8), intent(in) :: n
        real(kind=8), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        real(kind=8), intent(in) :: b(*)
        integer(kind=8), intent(in) :: incb
        real(kind=8), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        end subroutine mkl_vm_vdfmaxi_omp_offload_ilp64

        subroutine mkl_vm_vmdfmaxi_omp_offload_ilp64(n, a, inca, b, incb, y, incy, mode)
        integer(kind=8), intent(in) :: n
        real(kind=8), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        real(kind=8), intent(in) :: b(*)
        integer(kind=8), intent(in) :: incb
        real(kind=8), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmdfmaxi_omp_offload_ilp64




! function: Fmin, indexing: simple
        subroutine mkl_vm_vsfmin_omp_offload_ilp64(n, a, b, y)
        integer(kind=8), intent(in) :: n
        real(kind=4), intent(in) :: a(*)
        real(kind=4), intent(in) :: b(*)
        real(kind=4), intent(inout) :: y(*)
        end subroutine mkl_vm_vsfmin_omp_offload_ilp64

        subroutine mkl_vm_vmsfmin_omp_offload_ilp64(n, a, b, y, mode)
        integer(kind=8), intent(in) :: n
        real(kind=4), intent(in) :: a(*)
        real(kind=4), intent(in) :: b(*)
        real(kind=4), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmsfmin_omp_offload_ilp64

        subroutine mkl_vm_vdfmin_omp_offload_ilp64(n, a, b, y)
        integer(kind=8), intent(in) :: n
        real(kind=8), intent(in) :: a(*)
        real(kind=8), intent(in) :: b(*)
        real(kind=8), intent(inout) :: y(*)
        end subroutine mkl_vm_vdfmin_omp_offload_ilp64

        subroutine mkl_vm_vmdfmin_omp_offload_ilp64(n, a, b, y, mode)
        integer(kind=8), intent(in) :: n
        real(kind=8), intent(in) :: a(*)
        real(kind=8), intent(in) :: b(*)
        real(kind=8), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmdfmin_omp_offload_ilp64


! function: Fmin, indexing: strided
        subroutine mkl_vm_vsfmini_omp_offload_ilp64(n, a, inca, b, incb, y, incy)
        integer(kind=8), intent(in) :: n
        real(kind=4), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        real(kind=4), intent(in) :: b(*)
        integer(kind=8), intent(in) :: incb
        real(kind=4), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        end subroutine mkl_vm_vsfmini_omp_offload_ilp64

        subroutine mkl_vm_vmsfmini_omp_offload_ilp64(n, a, inca, b, incb, y, incy, mode)
        integer(kind=8), intent(in) :: n
        real(kind=4), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        real(kind=4), intent(in) :: b(*)
        integer(kind=8), intent(in) :: incb
        real(kind=4), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmsfmini_omp_offload_ilp64

        subroutine mkl_vm_vdfmini_omp_offload_ilp64(n, a, inca, b, incb, y, incy)
        integer(kind=8), intent(in) :: n
        real(kind=8), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        real(kind=8), intent(in) :: b(*)
        integer(kind=8), intent(in) :: incb
        real(kind=8), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        end subroutine mkl_vm_vdfmini_omp_offload_ilp64

        subroutine mkl_vm_vmdfmini_omp_offload_ilp64(n, a, inca, b, incb, y, incy, mode)
        integer(kind=8), intent(in) :: n
        real(kind=8), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        real(kind=8), intent(in) :: b(*)
        integer(kind=8), intent(in) :: incb
        real(kind=8), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmdfmini_omp_offload_ilp64




! function: Fmod, indexing: simple
        subroutine mkl_vm_vsfmod_omp_offload_ilp64(n, a, b, y)
        integer(kind=8), intent(in) :: n
        real(kind=4), intent(in) :: a(*)
        real(kind=4), intent(in) :: b(*)
        real(kind=4), intent(inout) :: y(*)
        end subroutine mkl_vm_vsfmod_omp_offload_ilp64

        subroutine mkl_vm_vmsfmod_omp_offload_ilp64(n, a, b, y, mode)
        integer(kind=8), intent(in) :: n
        real(kind=4), intent(in) :: a(*)
        real(kind=4), intent(in) :: b(*)
        real(kind=4), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmsfmod_omp_offload_ilp64

        subroutine mkl_vm_vdfmod_omp_offload_ilp64(n, a, b, y)
        integer(kind=8), intent(in) :: n
        real(kind=8), intent(in) :: a(*)
        real(kind=8), intent(in) :: b(*)
        real(kind=8), intent(inout) :: y(*)
        end subroutine mkl_vm_vdfmod_omp_offload_ilp64

        subroutine mkl_vm_vmdfmod_omp_offload_ilp64(n, a, b, y, mode)
        integer(kind=8), intent(in) :: n
        real(kind=8), intent(in) :: a(*)
        real(kind=8), intent(in) :: b(*)
        real(kind=8), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmdfmod_omp_offload_ilp64


! function: Fmod, indexing: strided
        subroutine mkl_vm_vsfmodi_omp_offload_ilp64(n, a, inca, b, incb, y, incy)
        integer(kind=8), intent(in) :: n
        real(kind=4), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        real(kind=4), intent(in) :: b(*)
        integer(kind=8), intent(in) :: incb
        real(kind=4), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        end subroutine mkl_vm_vsfmodi_omp_offload_ilp64

        subroutine mkl_vm_vmsfmodi_omp_offload_ilp64(n, a, inca, b, incb, y, incy, mode)
        integer(kind=8), intent(in) :: n
        real(kind=4), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        real(kind=4), intent(in) :: b(*)
        integer(kind=8), intent(in) :: incb
        real(kind=4), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmsfmodi_omp_offload_ilp64

        subroutine mkl_vm_vdfmodi_omp_offload_ilp64(n, a, inca, b, incb, y, incy)
        integer(kind=8), intent(in) :: n
        real(kind=8), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        real(kind=8), intent(in) :: b(*)
        integer(kind=8), intent(in) :: incb
        real(kind=8), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        end subroutine mkl_vm_vdfmodi_omp_offload_ilp64

        subroutine mkl_vm_vmdfmodi_omp_offload_ilp64(n, a, inca, b, incb, y, incy, mode)
        integer(kind=8), intent(in) :: n
        real(kind=8), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        real(kind=8), intent(in) :: b(*)
        integer(kind=8), intent(in) :: incb
        real(kind=8), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmdfmodi_omp_offload_ilp64




! function: Frac, indexing: simple
        subroutine mkl_vm_vsfrac_omp_offload_ilp64(n, a, y)
        integer(kind=8), intent(in) :: n
        real(kind=4), intent(in) :: a(*)
        real(kind=4), intent(inout) :: y(*)
        end subroutine mkl_vm_vsfrac_omp_offload_ilp64

        subroutine mkl_vm_vmsfrac_omp_offload_ilp64(n, a, y, mode)
        integer(kind=8), intent(in) :: n
        real(kind=4), intent(in) :: a(*)
        real(kind=4), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmsfrac_omp_offload_ilp64

        subroutine mkl_vm_vdfrac_omp_offload_ilp64(n, a, y)
        integer(kind=8), intent(in) :: n
        real(kind=8), intent(in) :: a(*)
        real(kind=8), intent(inout) :: y(*)
        end subroutine mkl_vm_vdfrac_omp_offload_ilp64

        subroutine mkl_vm_vmdfrac_omp_offload_ilp64(n, a, y, mode)
        integer(kind=8), intent(in) :: n
        real(kind=8), intent(in) :: a(*)
        real(kind=8), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmdfrac_omp_offload_ilp64


! function: Frac, indexing: strided
        subroutine mkl_vm_vsfraci_omp_offload_ilp64(n, a, inca, y, incy)
        integer(kind=8), intent(in) :: n
        real(kind=4), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        real(kind=4), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        end subroutine mkl_vm_vsfraci_omp_offload_ilp64

        subroutine mkl_vm_vmsfraci_omp_offload_ilp64(n, a, inca, y, incy, mode)
        integer(kind=8), intent(in) :: n
        real(kind=4), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        real(kind=4), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmsfraci_omp_offload_ilp64

        subroutine mkl_vm_vdfraci_omp_offload_ilp64(n, a, inca, y, incy)
        integer(kind=8), intent(in) :: n
        real(kind=8), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        real(kind=8), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        end subroutine mkl_vm_vdfraci_omp_offload_ilp64

        subroutine mkl_vm_vmdfraci_omp_offload_ilp64(n, a, inca, y, incy, mode)
        integer(kind=8), intent(in) :: n
        real(kind=8), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        real(kind=8), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmdfraci_omp_offload_ilp64




! function: Hypot, indexing: simple
        subroutine mkl_vm_vshypot_omp_offload_ilp64(n, a, b, y)
        integer(kind=8), intent(in) :: n
        real(kind=4), intent(in) :: a(*)
        real(kind=4), intent(in) :: b(*)
        real(kind=4), intent(inout) :: y(*)
        end subroutine mkl_vm_vshypot_omp_offload_ilp64

        subroutine mkl_vm_vmshypot_omp_offload_ilp64(n, a, b, y, mode)
        integer(kind=8), intent(in) :: n
        real(kind=4), intent(in) :: a(*)
        real(kind=4), intent(in) :: b(*)
        real(kind=4), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmshypot_omp_offload_ilp64

        subroutine mkl_vm_vdhypot_omp_offload_ilp64(n, a, b, y)
        integer(kind=8), intent(in) :: n
        real(kind=8), intent(in) :: a(*)
        real(kind=8), intent(in) :: b(*)
        real(kind=8), intent(inout) :: y(*)
        end subroutine mkl_vm_vdhypot_omp_offload_ilp64

        subroutine mkl_vm_vmdhypot_omp_offload_ilp64(n, a, b, y, mode)
        integer(kind=8), intent(in) :: n
        real(kind=8), intent(in) :: a(*)
        real(kind=8), intent(in) :: b(*)
        real(kind=8), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmdhypot_omp_offload_ilp64


! function: Hypot, indexing: strided
        subroutine mkl_vm_vshypoti_omp_offload_ilp64(n, a, inca, b, incb, y, incy)
        integer(kind=8), intent(in) :: n
        real(kind=4), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        real(kind=4), intent(in) :: b(*)
        integer(kind=8), intent(in) :: incb
        real(kind=4), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        end subroutine mkl_vm_vshypoti_omp_offload_ilp64

        subroutine mkl_vm_vmshypoti_omp_offload_ilp64(n, a, inca, b, incb, y, incy, mode)
        integer(kind=8), intent(in) :: n
        real(kind=4), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        real(kind=4), intent(in) :: b(*)
        integer(kind=8), intent(in) :: incb
        real(kind=4), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmshypoti_omp_offload_ilp64

        subroutine mkl_vm_vdhypoti_omp_offload_ilp64(n, a, inca, b, incb, y, incy)
        integer(kind=8), intent(in) :: n
        real(kind=8), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        real(kind=8), intent(in) :: b(*)
        integer(kind=8), intent(in) :: incb
        real(kind=8), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        end subroutine mkl_vm_vdhypoti_omp_offload_ilp64

        subroutine mkl_vm_vmdhypoti_omp_offload_ilp64(n, a, inca, b, incb, y, incy, mode)
        integer(kind=8), intent(in) :: n
        real(kind=8), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        real(kind=8), intent(in) :: b(*)
        integer(kind=8), intent(in) :: incb
        real(kind=8), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmdhypoti_omp_offload_ilp64




! function: Inv, indexing: simple
        subroutine mkl_vm_vsinv_omp_offload_ilp64(n, a, y)
        integer(kind=8), intent(in) :: n
        real(kind=4), intent(in) :: a(*)
        real(kind=4), intent(inout) :: y(*)
        end subroutine mkl_vm_vsinv_omp_offload_ilp64

        subroutine mkl_vm_vmsinv_omp_offload_ilp64(n, a, y, mode)
        integer(kind=8), intent(in) :: n
        real(kind=4), intent(in) :: a(*)
        real(kind=4), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmsinv_omp_offload_ilp64

        subroutine mkl_vm_vdinv_omp_offload_ilp64(n, a, y)
        integer(kind=8), intent(in) :: n
        real(kind=8), intent(in) :: a(*)
        real(kind=8), intent(inout) :: y(*)
        end subroutine mkl_vm_vdinv_omp_offload_ilp64

        subroutine mkl_vm_vmdinv_omp_offload_ilp64(n, a, y, mode)
        integer(kind=8), intent(in) :: n
        real(kind=8), intent(in) :: a(*)
        real(kind=8), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmdinv_omp_offload_ilp64


! function: Inv, indexing: strided
        subroutine mkl_vm_vsinvi_omp_offload_ilp64(n, a, inca, y, incy)
        integer(kind=8), intent(in) :: n
        real(kind=4), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        real(kind=4), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        end subroutine mkl_vm_vsinvi_omp_offload_ilp64

        subroutine mkl_vm_vmsinvi_omp_offload_ilp64(n, a, inca, y, incy, mode)
        integer(kind=8), intent(in) :: n
        real(kind=4), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        real(kind=4), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmsinvi_omp_offload_ilp64

        subroutine mkl_vm_vdinvi_omp_offload_ilp64(n, a, inca, y, incy)
        integer(kind=8), intent(in) :: n
        real(kind=8), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        real(kind=8), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        end subroutine mkl_vm_vdinvi_omp_offload_ilp64

        subroutine mkl_vm_vmdinvi_omp_offload_ilp64(n, a, inca, y, incy, mode)
        integer(kind=8), intent(in) :: n
        real(kind=8), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        real(kind=8), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmdinvi_omp_offload_ilp64




! function: InvCbrt, indexing: simple
        subroutine mkl_vm_vsinvcbrt_omp_offload_ilp64(n, a, y)
        integer(kind=8), intent(in) :: n
        real(kind=4), intent(in) :: a(*)
        real(kind=4), intent(inout) :: y(*)
        end subroutine mkl_vm_vsinvcbrt_omp_offload_ilp64

        subroutine mkl_vm_vmsinvcbrt_omp_offload_ilp64(n, a, y, mode)
        integer(kind=8), intent(in) :: n
        real(kind=4), intent(in) :: a(*)
        real(kind=4), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmsinvcbrt_omp_offload_ilp64

        subroutine mkl_vm_vdinvcbrt_omp_offload_ilp64(n, a, y)
        integer(kind=8), intent(in) :: n
        real(kind=8), intent(in) :: a(*)
        real(kind=8), intent(inout) :: y(*)
        end subroutine mkl_vm_vdinvcbrt_omp_offload_ilp64

        subroutine mkl_vm_vmdinvcbrt_omp_offload_ilp64(n, a, y, mode)
        integer(kind=8), intent(in) :: n
        real(kind=8), intent(in) :: a(*)
        real(kind=8), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmdinvcbrt_omp_offload_ilp64


! function: InvCbrt, indexing: strided
        subroutine mkl_vm_vsinvcbrti_omp_offload_ilp64(n, a, inca, y, incy)
        integer(kind=8), intent(in) :: n
        real(kind=4), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        real(kind=4), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        end subroutine mkl_vm_vsinvcbrti_omp_offload_ilp64

        subroutine mkl_vm_vmsinvcbrti_omp_offload_ilp64(n, a, inca, y, incy, mode)
        integer(kind=8), intent(in) :: n
        real(kind=4), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        real(kind=4), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmsinvcbrti_omp_offload_ilp64

        subroutine mkl_vm_vdinvcbrti_omp_offload_ilp64(n, a, inca, y, incy)
        integer(kind=8), intent(in) :: n
        real(kind=8), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        real(kind=8), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        end subroutine mkl_vm_vdinvcbrti_omp_offload_ilp64

        subroutine mkl_vm_vmdinvcbrti_omp_offload_ilp64(n, a, inca, y, incy, mode)
        integer(kind=8), intent(in) :: n
        real(kind=8), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        real(kind=8), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmdinvcbrti_omp_offload_ilp64




! function: InvSqrt, indexing: simple
        subroutine mkl_vm_vsinvsqrt_omp_offload_ilp64(n, a, y)
        integer(kind=8), intent(in) :: n
        real(kind=4), intent(in) :: a(*)
        real(kind=4), intent(inout) :: y(*)
        end subroutine mkl_vm_vsinvsqrt_omp_offload_ilp64

        subroutine mkl_vm_vmsinvsqrt_omp_offload_ilp64(n, a, y, mode)
        integer(kind=8), intent(in) :: n
        real(kind=4), intent(in) :: a(*)
        real(kind=4), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmsinvsqrt_omp_offload_ilp64

        subroutine mkl_vm_vdinvsqrt_omp_offload_ilp64(n, a, y)
        integer(kind=8), intent(in) :: n
        real(kind=8), intent(in) :: a(*)
        real(kind=8), intent(inout) :: y(*)
        end subroutine mkl_vm_vdinvsqrt_omp_offload_ilp64

        subroutine mkl_vm_vmdinvsqrt_omp_offload_ilp64(n, a, y, mode)
        integer(kind=8), intent(in) :: n
        real(kind=8), intent(in) :: a(*)
        real(kind=8), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmdinvsqrt_omp_offload_ilp64


! function: InvSqrt, indexing: strided
        subroutine mkl_vm_vsinvsqrti_omp_offload_ilp64(n, a, inca, y, incy)
        integer(kind=8), intent(in) :: n
        real(kind=4), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        real(kind=4), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        end subroutine mkl_vm_vsinvsqrti_omp_offload_ilp64

        subroutine mkl_vm_vmsinvsqrti_omp_offload_ilp64(n, a, inca, y, incy, mode)
        integer(kind=8), intent(in) :: n
        real(kind=4), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        real(kind=4), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmsinvsqrti_omp_offload_ilp64

        subroutine mkl_vm_vdinvsqrti_omp_offload_ilp64(n, a, inca, y, incy)
        integer(kind=8), intent(in) :: n
        real(kind=8), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        real(kind=8), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        end subroutine mkl_vm_vdinvsqrti_omp_offload_ilp64

        subroutine mkl_vm_vmdinvsqrti_omp_offload_ilp64(n, a, inca, y, incy, mode)
        integer(kind=8), intent(in) :: n
        real(kind=8), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        real(kind=8), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmdinvsqrti_omp_offload_ilp64




! function: LGamma, indexing: simple
        subroutine mkl_vm_vslgamma_omp_offload_ilp64(n, a, y)
        integer(kind=8), intent(in) :: n
        real(kind=4), intent(in) :: a(*)
        real(kind=4), intent(inout) :: y(*)
        end subroutine mkl_vm_vslgamma_omp_offload_ilp64

        subroutine mkl_vm_vmslgamma_omp_offload_ilp64(n, a, y, mode)
        integer(kind=8), intent(in) :: n
        real(kind=4), intent(in) :: a(*)
        real(kind=4), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmslgamma_omp_offload_ilp64

        subroutine mkl_vm_vdlgamma_omp_offload_ilp64(n, a, y)
        integer(kind=8), intent(in) :: n
        real(kind=8), intent(in) :: a(*)
        real(kind=8), intent(inout) :: y(*)
        end subroutine mkl_vm_vdlgamma_omp_offload_ilp64

        subroutine mkl_vm_vmdlgamma_omp_offload_ilp64(n, a, y, mode)
        integer(kind=8), intent(in) :: n
        real(kind=8), intent(in) :: a(*)
        real(kind=8), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmdlgamma_omp_offload_ilp64


! function: LGamma, indexing: strided
        subroutine mkl_vm_vslgammai_omp_offload_ilp64(n, a, inca, y, incy)
        integer(kind=8), intent(in) :: n
        real(kind=4), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        real(kind=4), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        end subroutine mkl_vm_vslgammai_omp_offload_ilp64

        subroutine mkl_vm_vmslgammai_omp_offload_ilp64(n, a, inca, y, incy, mode)
        integer(kind=8), intent(in) :: n
        real(kind=4), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        real(kind=4), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmslgammai_omp_offload_ilp64

        subroutine mkl_vm_vdlgammai_omp_offload_ilp64(n, a, inca, y, incy)
        integer(kind=8), intent(in) :: n
        real(kind=8), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        real(kind=8), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        end subroutine mkl_vm_vdlgammai_omp_offload_ilp64

        subroutine mkl_vm_vmdlgammai_omp_offload_ilp64(n, a, inca, y, incy, mode)
        integer(kind=8), intent(in) :: n
        real(kind=8), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        real(kind=8), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmdlgammai_omp_offload_ilp64




! function: LinearFrac, indexing: simple
        subroutine mkl_vm_vslinearfrac_omp_offload_ilp64(n, a, b, c, d, e, f, y)
        integer(kind=8), intent(in) :: n
        real(kind=4), intent(in) :: a(*)
        real(kind=4), intent(in) :: b(*)
        real(kind=4), intent(in) :: c
        real(kind=4), intent(in) :: d
        real(kind=4), intent(in) :: e
        real(kind=4), intent(in) :: f
        real(kind=4), intent(inout) :: y(*)
        end subroutine mkl_vm_vslinearfrac_omp_offload_ilp64

        subroutine mkl_vm_vmslinearfrac_omp_offload_ilp64(n, a, b, c, d, e, f, y, mode)
        integer(kind=8), intent(in) :: n
        real(kind=4), intent(in) :: a(*)
        real(kind=4), intent(in) :: b(*)
        real(kind=4), intent(in) :: c
        real(kind=4), intent(in) :: d
        real(kind=4), intent(in) :: e
        real(kind=4), intent(in) :: f
        real(kind=4), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmslinearfrac_omp_offload_ilp64

        subroutine mkl_vm_vdlinearfrac_omp_offload_ilp64(n, a, b, c, d, e, f, y)
        integer(kind=8), intent(in) :: n
        real(kind=8), intent(in) :: a(*)
        real(kind=8), intent(in) :: b(*)
        real(kind=8), intent(in) :: c
        real(kind=8), intent(in) :: d
        real(kind=8), intent(in) :: e
        real(kind=8), intent(in) :: f
        real(kind=8), intent(inout) :: y(*)
        end subroutine mkl_vm_vdlinearfrac_omp_offload_ilp64

        subroutine mkl_vm_vmdlinearfrac_omp_offload_ilp64(n, a, b, c, d, e, f, y, mode)
        integer(kind=8), intent(in) :: n
        real(kind=8), intent(in) :: a(*)
        real(kind=8), intent(in) :: b(*)
        real(kind=8), intent(in) :: c
        real(kind=8), intent(in) :: d
        real(kind=8), intent(in) :: e
        real(kind=8), intent(in) :: f
        real(kind=8), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmdlinearfrac_omp_offload_ilp64


! function: LinearFrac, indexing: strided
        subroutine mkl_vm_vslinearfraci_omp_offload_ilp64(n, a, inca, b, incb, c, d, e, f, y, incy)
        integer(kind=8), intent(in) :: n
        real(kind=4), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        real(kind=4), intent(in) :: b(*)
        integer(kind=8), intent(in) :: incb
        real(kind=4), intent(in) :: c
        real(kind=4), intent(in) :: d
        real(kind=4), intent(in) :: e
        real(kind=4), intent(in) :: f
        real(kind=4), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        end subroutine mkl_vm_vslinearfraci_omp_offload_ilp64

        subroutine mkl_vm_vmslinearfraci_omp_offload_ilp64(n, a, inca, b, incb, c, d, e, f, y, incy, mode)
        integer(kind=8), intent(in) :: n
        real(kind=4), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        real(kind=4), intent(in) :: b(*)
        integer(kind=8), intent(in) :: incb
        real(kind=4), intent(in) :: c
        real(kind=4), intent(in) :: d
        real(kind=4), intent(in) :: e
        real(kind=4), intent(in) :: f
        real(kind=4), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmslinearfraci_omp_offload_ilp64

        subroutine mkl_vm_vdlinearfraci_omp_offload_ilp64(n, a, inca, b, incb, c, d, e, f, y, incy)
        integer(kind=8), intent(in) :: n
        real(kind=8), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        real(kind=8), intent(in) :: b(*)
        integer(kind=8), intent(in) :: incb
        real(kind=8), intent(in) :: c
        real(kind=8), intent(in) :: d
        real(kind=8), intent(in) :: e
        real(kind=8), intent(in) :: f
        real(kind=8), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        end subroutine mkl_vm_vdlinearfraci_omp_offload_ilp64

        subroutine mkl_vm_vmdlinearfraci_omp_offload_ilp64(n, a, inca, b, incb, c, d, e, f, y, incy, mode)
        integer(kind=8), intent(in) :: n
        real(kind=8), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        real(kind=8), intent(in) :: b(*)
        integer(kind=8), intent(in) :: incb
        real(kind=8), intent(in) :: c
        real(kind=8), intent(in) :: d
        real(kind=8), intent(in) :: e
        real(kind=8), intent(in) :: f
        real(kind=8), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmdlinearfraci_omp_offload_ilp64




! function: Ln, indexing: simple
        subroutine mkl_vm_vsln_omp_offload_ilp64(n, a, y)
        integer(kind=8), intent(in) :: n
        real(kind=4), intent(in) :: a(*)
        real(kind=4), intent(inout) :: y(*)
        end subroutine mkl_vm_vsln_omp_offload_ilp64

        subroutine mkl_vm_vmsln_omp_offload_ilp64(n, a, y, mode)
        integer(kind=8), intent(in) :: n
        real(kind=4), intent(in) :: a(*)
        real(kind=4), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmsln_omp_offload_ilp64

        subroutine mkl_vm_vdln_omp_offload_ilp64(n, a, y)
        integer(kind=8), intent(in) :: n
        real(kind=8), intent(in) :: a(*)
        real(kind=8), intent(inout) :: y(*)
        end subroutine mkl_vm_vdln_omp_offload_ilp64

        subroutine mkl_vm_vmdln_omp_offload_ilp64(n, a, y, mode)
        integer(kind=8), intent(in) :: n
        real(kind=8), intent(in) :: a(*)
        real(kind=8), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmdln_omp_offload_ilp64

        subroutine mkl_vm_vcln_omp_offload_ilp64(n, a, y)
        integer(kind=8), intent(in) :: n
        complex(kind=4), intent(in) :: a(*)
        complex(kind=4), intent(inout) :: y(*)
        end subroutine mkl_vm_vcln_omp_offload_ilp64

        subroutine mkl_vm_vmcln_omp_offload_ilp64(n, a, y, mode)
        integer(kind=8), intent(in) :: n
        complex(kind=4), intent(in) :: a(*)
        complex(kind=4), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmcln_omp_offload_ilp64

        subroutine mkl_vm_vzln_omp_offload_ilp64(n, a, y)
        integer(kind=8), intent(in) :: n
        complex(kind=8), intent(in) :: a(*)
        complex(kind=8), intent(inout) :: y(*)
        end subroutine mkl_vm_vzln_omp_offload_ilp64

        subroutine mkl_vm_vmzln_omp_offload_ilp64(n, a, y, mode)
        integer(kind=8), intent(in) :: n
        complex(kind=8), intent(in) :: a(*)
        complex(kind=8), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmzln_omp_offload_ilp64


! function: Ln, indexing: strided
        subroutine mkl_vm_vslni_omp_offload_ilp64(n, a, inca, y, incy)
        integer(kind=8), intent(in) :: n
        real(kind=4), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        real(kind=4), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        end subroutine mkl_vm_vslni_omp_offload_ilp64

        subroutine mkl_vm_vmslni_omp_offload_ilp64(n, a, inca, y, incy, mode)
        integer(kind=8), intent(in) :: n
        real(kind=4), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        real(kind=4), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmslni_omp_offload_ilp64

        subroutine mkl_vm_vdlni_omp_offload_ilp64(n, a, inca, y, incy)
        integer(kind=8), intent(in) :: n
        real(kind=8), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        real(kind=8), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        end subroutine mkl_vm_vdlni_omp_offload_ilp64

        subroutine mkl_vm_vmdlni_omp_offload_ilp64(n, a, inca, y, incy, mode)
        integer(kind=8), intent(in) :: n
        real(kind=8), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        real(kind=8), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmdlni_omp_offload_ilp64

        subroutine mkl_vm_vclni_omp_offload_ilp64(n, a, inca, y, incy)
        integer(kind=8), intent(in) :: n
        complex(kind=4), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        complex(kind=4), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        end subroutine mkl_vm_vclni_omp_offload_ilp64

        subroutine mkl_vm_vmclni_omp_offload_ilp64(n, a, inca, y, incy, mode)
        integer(kind=8), intent(in) :: n
        complex(kind=4), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        complex(kind=4), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmclni_omp_offload_ilp64

        subroutine mkl_vm_vzlni_omp_offload_ilp64(n, a, inca, y, incy)
        integer(kind=8), intent(in) :: n
        complex(kind=8), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        complex(kind=8), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        end subroutine mkl_vm_vzlni_omp_offload_ilp64

        subroutine mkl_vm_vmzlni_omp_offload_ilp64(n, a, inca, y, incy, mode)
        integer(kind=8), intent(in) :: n
        complex(kind=8), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        complex(kind=8), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmzlni_omp_offload_ilp64




! function: Log10, indexing: simple
        subroutine mkl_vm_vslog10_omp_offload_ilp64(n, a, y)
        integer(kind=8), intent(in) :: n
        real(kind=4), intent(in) :: a(*)
        real(kind=4), intent(inout) :: y(*)
        end subroutine mkl_vm_vslog10_omp_offload_ilp64

        subroutine mkl_vm_vmslog10_omp_offload_ilp64(n, a, y, mode)
        integer(kind=8), intent(in) :: n
        real(kind=4), intent(in) :: a(*)
        real(kind=4), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmslog10_omp_offload_ilp64

        subroutine mkl_vm_vdlog10_omp_offload_ilp64(n, a, y)
        integer(kind=8), intent(in) :: n
        real(kind=8), intent(in) :: a(*)
        real(kind=8), intent(inout) :: y(*)
        end subroutine mkl_vm_vdlog10_omp_offload_ilp64

        subroutine mkl_vm_vmdlog10_omp_offload_ilp64(n, a, y, mode)
        integer(kind=8), intent(in) :: n
        real(kind=8), intent(in) :: a(*)
        real(kind=8), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmdlog10_omp_offload_ilp64

        subroutine mkl_vm_vclog10_omp_offload_ilp64(n, a, y)
        integer(kind=8), intent(in) :: n
        complex(kind=4), intent(in) :: a(*)
        complex(kind=4), intent(inout) :: y(*)
        end subroutine mkl_vm_vclog10_omp_offload_ilp64

        subroutine mkl_vm_vmclog10_omp_offload_ilp64(n, a, y, mode)
        integer(kind=8), intent(in) :: n
        complex(kind=4), intent(in) :: a(*)
        complex(kind=4), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmclog10_omp_offload_ilp64

        subroutine mkl_vm_vzlog10_omp_offload_ilp64(n, a, y)
        integer(kind=8), intent(in) :: n
        complex(kind=8), intent(in) :: a(*)
        complex(kind=8), intent(inout) :: y(*)
        end subroutine mkl_vm_vzlog10_omp_offload_ilp64

        subroutine mkl_vm_vmzlog10_omp_offload_ilp64(n, a, y, mode)
        integer(kind=8), intent(in) :: n
        complex(kind=8), intent(in) :: a(*)
        complex(kind=8), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmzlog10_omp_offload_ilp64


! function: Log10, indexing: strided
        subroutine mkl_vm_vslog10i_omp_offload_ilp64(n, a, inca, y, incy)
        integer(kind=8), intent(in) :: n
        real(kind=4), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        real(kind=4), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        end subroutine mkl_vm_vslog10i_omp_offload_ilp64

        subroutine mkl_vm_vmslog10i_omp_offload_ilp64(n, a, inca, y, incy, mode)
        integer(kind=8), intent(in) :: n
        real(kind=4), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        real(kind=4), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmslog10i_omp_offload_ilp64

        subroutine mkl_vm_vdlog10i_omp_offload_ilp64(n, a, inca, y, incy)
        integer(kind=8), intent(in) :: n
        real(kind=8), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        real(kind=8), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        end subroutine mkl_vm_vdlog10i_omp_offload_ilp64

        subroutine mkl_vm_vmdlog10i_omp_offload_ilp64(n, a, inca, y, incy, mode)
        integer(kind=8), intent(in) :: n
        real(kind=8), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        real(kind=8), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmdlog10i_omp_offload_ilp64

        subroutine mkl_vm_vclog10i_omp_offload_ilp64(n, a, inca, y, incy)
        integer(kind=8), intent(in) :: n
        complex(kind=4), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        complex(kind=4), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        end subroutine mkl_vm_vclog10i_omp_offload_ilp64

        subroutine mkl_vm_vmclog10i_omp_offload_ilp64(n, a, inca, y, incy, mode)
        integer(kind=8), intent(in) :: n
        complex(kind=4), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        complex(kind=4), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmclog10i_omp_offload_ilp64

        subroutine mkl_vm_vzlog10i_omp_offload_ilp64(n, a, inca, y, incy)
        integer(kind=8), intent(in) :: n
        complex(kind=8), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        complex(kind=8), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        end subroutine mkl_vm_vzlog10i_omp_offload_ilp64

        subroutine mkl_vm_vmzlog10i_omp_offload_ilp64(n, a, inca, y, incy, mode)
        integer(kind=8), intent(in) :: n
        complex(kind=8), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        complex(kind=8), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmzlog10i_omp_offload_ilp64




! function: Log1p, indexing: simple
        subroutine mkl_vm_vslog1p_omp_offload_ilp64(n, a, y)
        integer(kind=8), intent(in) :: n
        real(kind=4), intent(in) :: a(*)
        real(kind=4), intent(inout) :: y(*)
        end subroutine mkl_vm_vslog1p_omp_offload_ilp64

        subroutine mkl_vm_vmslog1p_omp_offload_ilp64(n, a, y, mode)
        integer(kind=8), intent(in) :: n
        real(kind=4), intent(in) :: a(*)
        real(kind=4), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmslog1p_omp_offload_ilp64

        subroutine mkl_vm_vdlog1p_omp_offload_ilp64(n, a, y)
        integer(kind=8), intent(in) :: n
        real(kind=8), intent(in) :: a(*)
        real(kind=8), intent(inout) :: y(*)
        end subroutine mkl_vm_vdlog1p_omp_offload_ilp64

        subroutine mkl_vm_vmdlog1p_omp_offload_ilp64(n, a, y, mode)
        integer(kind=8), intent(in) :: n
        real(kind=8), intent(in) :: a(*)
        real(kind=8), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmdlog1p_omp_offload_ilp64


! function: Log1p, indexing: strided
        subroutine mkl_vm_vslog1pi_omp_offload_ilp64(n, a, inca, y, incy)
        integer(kind=8), intent(in) :: n
        real(kind=4), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        real(kind=4), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        end subroutine mkl_vm_vslog1pi_omp_offload_ilp64

        subroutine mkl_vm_vmslog1pi_omp_offload_ilp64(n, a, inca, y, incy, mode)
        integer(kind=8), intent(in) :: n
        real(kind=4), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        real(kind=4), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmslog1pi_omp_offload_ilp64

        subroutine mkl_vm_vdlog1pi_omp_offload_ilp64(n, a, inca, y, incy)
        integer(kind=8), intent(in) :: n
        real(kind=8), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        real(kind=8), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        end subroutine mkl_vm_vdlog1pi_omp_offload_ilp64

        subroutine mkl_vm_vmdlog1pi_omp_offload_ilp64(n, a, inca, y, incy, mode)
        integer(kind=8), intent(in) :: n
        real(kind=8), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        real(kind=8), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmdlog1pi_omp_offload_ilp64




! function: Log2, indexing: simple
        subroutine mkl_vm_vslog2_omp_offload_ilp64(n, a, y)
        integer(kind=8), intent(in) :: n
        real(kind=4), intent(in) :: a(*)
        real(kind=4), intent(inout) :: y(*)
        end subroutine mkl_vm_vslog2_omp_offload_ilp64

        subroutine mkl_vm_vmslog2_omp_offload_ilp64(n, a, y, mode)
        integer(kind=8), intent(in) :: n
        real(kind=4), intent(in) :: a(*)
        real(kind=4), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmslog2_omp_offload_ilp64

        subroutine mkl_vm_vdlog2_omp_offload_ilp64(n, a, y)
        integer(kind=8), intent(in) :: n
        real(kind=8), intent(in) :: a(*)
        real(kind=8), intent(inout) :: y(*)
        end subroutine mkl_vm_vdlog2_omp_offload_ilp64

        subroutine mkl_vm_vmdlog2_omp_offload_ilp64(n, a, y, mode)
        integer(kind=8), intent(in) :: n
        real(kind=8), intent(in) :: a(*)
        real(kind=8), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmdlog2_omp_offload_ilp64


! function: Log2, indexing: strided
        subroutine mkl_vm_vslog2i_omp_offload_ilp64(n, a, inca, y, incy)
        integer(kind=8), intent(in) :: n
        real(kind=4), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        real(kind=4), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        end subroutine mkl_vm_vslog2i_omp_offload_ilp64

        subroutine mkl_vm_vmslog2i_omp_offload_ilp64(n, a, inca, y, incy, mode)
        integer(kind=8), intent(in) :: n
        real(kind=4), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        real(kind=4), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmslog2i_omp_offload_ilp64

        subroutine mkl_vm_vdlog2i_omp_offload_ilp64(n, a, inca, y, incy)
        integer(kind=8), intent(in) :: n
        real(kind=8), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        real(kind=8), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        end subroutine mkl_vm_vdlog2i_omp_offload_ilp64

        subroutine mkl_vm_vmdlog2i_omp_offload_ilp64(n, a, inca, y, incy, mode)
        integer(kind=8), intent(in) :: n
        real(kind=8), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        real(kind=8), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmdlog2i_omp_offload_ilp64




! function: Logb, indexing: simple
        subroutine mkl_vm_vslogb_omp_offload_ilp64(n, a, y)
        integer(kind=8), intent(in) :: n
        real(kind=4), intent(in) :: a(*)
        real(kind=4), intent(inout) :: y(*)
        end subroutine mkl_vm_vslogb_omp_offload_ilp64

        subroutine mkl_vm_vmslogb_omp_offload_ilp64(n, a, y, mode)
        integer(kind=8), intent(in) :: n
        real(kind=4), intent(in) :: a(*)
        real(kind=4), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmslogb_omp_offload_ilp64

        subroutine mkl_vm_vdlogb_omp_offload_ilp64(n, a, y)
        integer(kind=8), intent(in) :: n
        real(kind=8), intent(in) :: a(*)
        real(kind=8), intent(inout) :: y(*)
        end subroutine mkl_vm_vdlogb_omp_offload_ilp64

        subroutine mkl_vm_vmdlogb_omp_offload_ilp64(n, a, y, mode)
        integer(kind=8), intent(in) :: n
        real(kind=8), intent(in) :: a(*)
        real(kind=8), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmdlogb_omp_offload_ilp64


! function: Logb, indexing: strided
        subroutine mkl_vm_vslogbi_omp_offload_ilp64(n, a, inca, y, incy)
        integer(kind=8), intent(in) :: n
        real(kind=4), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        real(kind=4), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        end subroutine mkl_vm_vslogbi_omp_offload_ilp64

        subroutine mkl_vm_vmslogbi_omp_offload_ilp64(n, a, inca, y, incy, mode)
        integer(kind=8), intent(in) :: n
        real(kind=4), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        real(kind=4), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmslogbi_omp_offload_ilp64

        subroutine mkl_vm_vdlogbi_omp_offload_ilp64(n, a, inca, y, incy)
        integer(kind=8), intent(in) :: n
        real(kind=8), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        real(kind=8), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        end subroutine mkl_vm_vdlogbi_omp_offload_ilp64

        subroutine mkl_vm_vmdlogbi_omp_offload_ilp64(n, a, inca, y, incy, mode)
        integer(kind=8), intent(in) :: n
        real(kind=8), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        real(kind=8), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmdlogbi_omp_offload_ilp64




! function: MaxMag, indexing: simple
        subroutine mkl_vm_vsmaxmag_omp_offload_ilp64(n, a, b, y)
        integer(kind=8), intent(in) :: n
        real(kind=4), intent(in) :: a(*)
        real(kind=4), intent(in) :: b(*)
        real(kind=4), intent(inout) :: y(*)
        end subroutine mkl_vm_vsmaxmag_omp_offload_ilp64

        subroutine mkl_vm_vmsmaxmag_omp_offload_ilp64(n, a, b, y, mode)
        integer(kind=8), intent(in) :: n
        real(kind=4), intent(in) :: a(*)
        real(kind=4), intent(in) :: b(*)
        real(kind=4), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmsmaxmag_omp_offload_ilp64

        subroutine mkl_vm_vdmaxmag_omp_offload_ilp64(n, a, b, y)
        integer(kind=8), intent(in) :: n
        real(kind=8), intent(in) :: a(*)
        real(kind=8), intent(in) :: b(*)
        real(kind=8), intent(inout) :: y(*)
        end subroutine mkl_vm_vdmaxmag_omp_offload_ilp64

        subroutine mkl_vm_vmdmaxmag_omp_offload_ilp64(n, a, b, y, mode)
        integer(kind=8), intent(in) :: n
        real(kind=8), intent(in) :: a(*)
        real(kind=8), intent(in) :: b(*)
        real(kind=8), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmdmaxmag_omp_offload_ilp64


! function: MaxMag, indexing: strided
        subroutine mkl_vm_vsmaxmagi_omp_offload_ilp64(n, a, inca, b, incb, y, incy)
        integer(kind=8), intent(in) :: n
        real(kind=4), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        real(kind=4), intent(in) :: b(*)
        integer(kind=8), intent(in) :: incb
        real(kind=4), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        end subroutine mkl_vm_vsmaxmagi_omp_offload_ilp64

        subroutine mkl_vm_vmsmaxmagi_omp_offload_ilp64(n, a, inca, b, incb, y, incy, mode)
        integer(kind=8), intent(in) :: n
        real(kind=4), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        real(kind=4), intent(in) :: b(*)
        integer(kind=8), intent(in) :: incb
        real(kind=4), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmsmaxmagi_omp_offload_ilp64

        subroutine mkl_vm_vdmaxmagi_omp_offload_ilp64(n, a, inca, b, incb, y, incy)
        integer(kind=8), intent(in) :: n
        real(kind=8), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        real(kind=8), intent(in) :: b(*)
        integer(kind=8), intent(in) :: incb
        real(kind=8), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        end subroutine mkl_vm_vdmaxmagi_omp_offload_ilp64

        subroutine mkl_vm_vmdmaxmagi_omp_offload_ilp64(n, a, inca, b, incb, y, incy, mode)
        integer(kind=8), intent(in) :: n
        real(kind=8), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        real(kind=8), intent(in) :: b(*)
        integer(kind=8), intent(in) :: incb
        real(kind=8), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmdmaxmagi_omp_offload_ilp64




! function: MinMag, indexing: simple
        subroutine mkl_vm_vsminmag_omp_offload_ilp64(n, a, b, y)
        integer(kind=8), intent(in) :: n
        real(kind=4), intent(in) :: a(*)
        real(kind=4), intent(in) :: b(*)
        real(kind=4), intent(inout) :: y(*)
        end subroutine mkl_vm_vsminmag_omp_offload_ilp64

        subroutine mkl_vm_vmsminmag_omp_offload_ilp64(n, a, b, y, mode)
        integer(kind=8), intent(in) :: n
        real(kind=4), intent(in) :: a(*)
        real(kind=4), intent(in) :: b(*)
        real(kind=4), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmsminmag_omp_offload_ilp64

        subroutine mkl_vm_vdminmag_omp_offload_ilp64(n, a, b, y)
        integer(kind=8), intent(in) :: n
        real(kind=8), intent(in) :: a(*)
        real(kind=8), intent(in) :: b(*)
        real(kind=8), intent(inout) :: y(*)
        end subroutine mkl_vm_vdminmag_omp_offload_ilp64

        subroutine mkl_vm_vmdminmag_omp_offload_ilp64(n, a, b, y, mode)
        integer(kind=8), intent(in) :: n
        real(kind=8), intent(in) :: a(*)
        real(kind=8), intent(in) :: b(*)
        real(kind=8), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmdminmag_omp_offload_ilp64


! function: MinMag, indexing: strided
        subroutine mkl_vm_vsminmagi_omp_offload_ilp64(n, a, inca, b, incb, y, incy)
        integer(kind=8), intent(in) :: n
        real(kind=4), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        real(kind=4), intent(in) :: b(*)
        integer(kind=8), intent(in) :: incb
        real(kind=4), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        end subroutine mkl_vm_vsminmagi_omp_offload_ilp64

        subroutine mkl_vm_vmsminmagi_omp_offload_ilp64(n, a, inca, b, incb, y, incy, mode)
        integer(kind=8), intent(in) :: n
        real(kind=4), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        real(kind=4), intent(in) :: b(*)
        integer(kind=8), intent(in) :: incb
        real(kind=4), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmsminmagi_omp_offload_ilp64

        subroutine mkl_vm_vdminmagi_omp_offload_ilp64(n, a, inca, b, incb, y, incy)
        integer(kind=8), intent(in) :: n
        real(kind=8), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        real(kind=8), intent(in) :: b(*)
        integer(kind=8), intent(in) :: incb
        real(kind=8), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        end subroutine mkl_vm_vdminmagi_omp_offload_ilp64

        subroutine mkl_vm_vmdminmagi_omp_offload_ilp64(n, a, inca, b, incb, y, incy, mode)
        integer(kind=8), intent(in) :: n
        real(kind=8), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        real(kind=8), intent(in) :: b(*)
        integer(kind=8), intent(in) :: incb
        real(kind=8), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmdminmagi_omp_offload_ilp64




! function: Modf, indexing: simple
        subroutine mkl_vm_vsmodf_omp_offload_ilp64(n, a, y, z)
        integer(kind=8), intent(in) :: n
        real(kind=4), intent(in) :: a(*)
        real(kind=4), intent(inout) :: y(*)
        real(kind=4), intent(inout) :: z(*)
        end subroutine mkl_vm_vsmodf_omp_offload_ilp64

        subroutine mkl_vm_vmsmodf_omp_offload_ilp64(n, a, y, z, mode)
        integer(kind=8), intent(in) :: n
        real(kind=4), intent(in) :: a(*)
        real(kind=4), intent(inout) :: y(*)
        real(kind=4), intent(inout) :: z(*)
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmsmodf_omp_offload_ilp64

        subroutine mkl_vm_vdmodf_omp_offload_ilp64(n, a, y, z)
        integer(kind=8), intent(in) :: n
        real(kind=8), intent(in) :: a(*)
        real(kind=8), intent(inout) :: y(*)
        real(kind=8), intent(inout) :: z(*)
        end subroutine mkl_vm_vdmodf_omp_offload_ilp64

        subroutine mkl_vm_vmdmodf_omp_offload_ilp64(n, a, y, z, mode)
        integer(kind=8), intent(in) :: n
        real(kind=8), intent(in) :: a(*)
        real(kind=8), intent(inout) :: y(*)
        real(kind=8), intent(inout) :: z(*)
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmdmodf_omp_offload_ilp64


! function: Modf, indexing: strided
        subroutine mkl_vm_vsmodfi_omp_offload_ilp64(n, a, inca, y, incy, z, incz)
        integer(kind=8), intent(in) :: n
        real(kind=4), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        real(kind=4), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        real(kind=4), intent(inout) :: z(*)
        integer(kind=8), intent(in) :: incz
        end subroutine mkl_vm_vsmodfi_omp_offload_ilp64

        subroutine mkl_vm_vmsmodfi_omp_offload_ilp64(n, a, inca, y, incy, z, incz, mode)
        integer(kind=8), intent(in) :: n
        real(kind=4), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        real(kind=4), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        real(kind=4), intent(inout) :: z(*)
        integer(kind=8), intent(in) :: incz
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmsmodfi_omp_offload_ilp64

        subroutine mkl_vm_vdmodfi_omp_offload_ilp64(n, a, inca, y, incy, z, incz)
        integer(kind=8), intent(in) :: n
        real(kind=8), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        real(kind=8), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        real(kind=8), intent(inout) :: z(*)
        integer(kind=8), intent(in) :: incz
        end subroutine mkl_vm_vdmodfi_omp_offload_ilp64

        subroutine mkl_vm_vmdmodfi_omp_offload_ilp64(n, a, inca, y, incy, z, incz, mode)
        integer(kind=8), intent(in) :: n
        real(kind=8), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        real(kind=8), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        real(kind=8), intent(inout) :: z(*)
        integer(kind=8), intent(in) :: incz
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmdmodfi_omp_offload_ilp64




! function: Mul, indexing: simple
        subroutine mkl_vm_vsmul_omp_offload_ilp64(n, a, b, y)
        integer(kind=8), intent(in) :: n
        real(kind=4), intent(in) :: a(*)
        real(kind=4), intent(in) :: b(*)
        real(kind=4), intent(inout) :: y(*)
        end subroutine mkl_vm_vsmul_omp_offload_ilp64

        subroutine mkl_vm_vmsmul_omp_offload_ilp64(n, a, b, y, mode)
        integer(kind=8), intent(in) :: n
        real(kind=4), intent(in) :: a(*)
        real(kind=4), intent(in) :: b(*)
        real(kind=4), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmsmul_omp_offload_ilp64

        subroutine mkl_vm_vdmul_omp_offload_ilp64(n, a, b, y)
        integer(kind=8), intent(in) :: n
        real(kind=8), intent(in) :: a(*)
        real(kind=8), intent(in) :: b(*)
        real(kind=8), intent(inout) :: y(*)
        end subroutine mkl_vm_vdmul_omp_offload_ilp64

        subroutine mkl_vm_vmdmul_omp_offload_ilp64(n, a, b, y, mode)
        integer(kind=8), intent(in) :: n
        real(kind=8), intent(in) :: a(*)
        real(kind=8), intent(in) :: b(*)
        real(kind=8), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmdmul_omp_offload_ilp64

        subroutine mkl_vm_vcmul_omp_offload_ilp64(n, a, b, y)
        integer(kind=8), intent(in) :: n
        complex(kind=4), intent(in) :: a(*)
        complex(kind=4), intent(in) :: b(*)
        complex(kind=4), intent(inout) :: y(*)
        end subroutine mkl_vm_vcmul_omp_offload_ilp64

        subroutine mkl_vm_vmcmul_omp_offload_ilp64(n, a, b, y, mode)
        integer(kind=8), intent(in) :: n
        complex(kind=4), intent(in) :: a(*)
        complex(kind=4), intent(in) :: b(*)
        complex(kind=4), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmcmul_omp_offload_ilp64

        subroutine mkl_vm_vzmul_omp_offload_ilp64(n, a, b, y)
        integer(kind=8), intent(in) :: n
        complex(kind=8), intent(in) :: a(*)
        complex(kind=8), intent(in) :: b(*)
        complex(kind=8), intent(inout) :: y(*)
        end subroutine mkl_vm_vzmul_omp_offload_ilp64

        subroutine mkl_vm_vmzmul_omp_offload_ilp64(n, a, b, y, mode)
        integer(kind=8), intent(in) :: n
        complex(kind=8), intent(in) :: a(*)
        complex(kind=8), intent(in) :: b(*)
        complex(kind=8), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmzmul_omp_offload_ilp64


! function: Mul, indexing: strided
        subroutine mkl_vm_vsmuli_omp_offload_ilp64(n, a, inca, b, incb, y, incy)
        integer(kind=8), intent(in) :: n
        real(kind=4), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        real(kind=4), intent(in) :: b(*)
        integer(kind=8), intent(in) :: incb
        real(kind=4), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        end subroutine mkl_vm_vsmuli_omp_offload_ilp64

        subroutine mkl_vm_vmsmuli_omp_offload_ilp64(n, a, inca, b, incb, y, incy, mode)
        integer(kind=8), intent(in) :: n
        real(kind=4), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        real(kind=4), intent(in) :: b(*)
        integer(kind=8), intent(in) :: incb
        real(kind=4), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmsmuli_omp_offload_ilp64

        subroutine mkl_vm_vdmuli_omp_offload_ilp64(n, a, inca, b, incb, y, incy)
        integer(kind=8), intent(in) :: n
        real(kind=8), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        real(kind=8), intent(in) :: b(*)
        integer(kind=8), intent(in) :: incb
        real(kind=8), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        end subroutine mkl_vm_vdmuli_omp_offload_ilp64

        subroutine mkl_vm_vmdmuli_omp_offload_ilp64(n, a, inca, b, incb, y, incy, mode)
        integer(kind=8), intent(in) :: n
        real(kind=8), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        real(kind=8), intent(in) :: b(*)
        integer(kind=8), intent(in) :: incb
        real(kind=8), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmdmuli_omp_offload_ilp64

        subroutine mkl_vm_vcmuli_omp_offload_ilp64(n, a, inca, b, incb, y, incy)
        integer(kind=8), intent(in) :: n
        complex(kind=4), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        complex(kind=4), intent(in) :: b(*)
        integer(kind=8), intent(in) :: incb
        complex(kind=4), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        end subroutine mkl_vm_vcmuli_omp_offload_ilp64

        subroutine mkl_vm_vmcmuli_omp_offload_ilp64(n, a, inca, b, incb, y, incy, mode)
        integer(kind=8), intent(in) :: n
        complex(kind=4), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        complex(kind=4), intent(in) :: b(*)
        integer(kind=8), intent(in) :: incb
        complex(kind=4), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmcmuli_omp_offload_ilp64

        subroutine mkl_vm_vzmuli_omp_offload_ilp64(n, a, inca, b, incb, y, incy)
        integer(kind=8), intent(in) :: n
        complex(kind=8), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        complex(kind=8), intent(in) :: b(*)
        integer(kind=8), intent(in) :: incb
        complex(kind=8), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        end subroutine mkl_vm_vzmuli_omp_offload_ilp64

        subroutine mkl_vm_vmzmuli_omp_offload_ilp64(n, a, inca, b, incb, y, incy, mode)
        integer(kind=8), intent(in) :: n
        complex(kind=8), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        complex(kind=8), intent(in) :: b(*)
        integer(kind=8), intent(in) :: incb
        complex(kind=8), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmzmuli_omp_offload_ilp64




! function: MulByConj, indexing: simple
        subroutine mkl_vm_vcmulbyconj_omp_offload_ilp64(n, a, b, y)
        integer(kind=8), intent(in) :: n
        complex(kind=4), intent(in) :: a(*)
        complex(kind=4), intent(in) :: b(*)
        complex(kind=4), intent(inout) :: y(*)
        end subroutine mkl_vm_vcmulbyconj_omp_offload_ilp64

        subroutine mkl_vm_vmcmulbyconj_omp_offload_ilp64(n, a, b, y, mode)
        integer(kind=8), intent(in) :: n
        complex(kind=4), intent(in) :: a(*)
        complex(kind=4), intent(in) :: b(*)
        complex(kind=4), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmcmulbyconj_omp_offload_ilp64

        subroutine mkl_vm_vzmulbyconj_omp_offload_ilp64(n, a, b, y)
        integer(kind=8), intent(in) :: n
        complex(kind=8), intent(in) :: a(*)
        complex(kind=8), intent(in) :: b(*)
        complex(kind=8), intent(inout) :: y(*)
        end subroutine mkl_vm_vzmulbyconj_omp_offload_ilp64

        subroutine mkl_vm_vmzmulbyconj_omp_offload_ilp64(n, a, b, y, mode)
        integer(kind=8), intent(in) :: n
        complex(kind=8), intent(in) :: a(*)
        complex(kind=8), intent(in) :: b(*)
        complex(kind=8), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmzmulbyconj_omp_offload_ilp64


! function: MulByConj, indexing: strided
        subroutine mkl_vm_vcmulbyconji_omp_offload_ilp64(n, a, inca, b, incb, y, incy)
        integer(kind=8), intent(in) :: n
        complex(kind=4), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        complex(kind=4), intent(in) :: b(*)
        integer(kind=8), intent(in) :: incb
        complex(kind=4), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        end subroutine mkl_vm_vcmulbyconji_omp_offload_ilp64

        subroutine mkl_vm_vmcmulbyconji_omp_offload_ilp64(n, a, inca, b, incb, y, incy, mode)
        integer(kind=8), intent(in) :: n
        complex(kind=4), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        complex(kind=4), intent(in) :: b(*)
        integer(kind=8), intent(in) :: incb
        complex(kind=4), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmcmulbyconji_omp_offload_ilp64

        subroutine mkl_vm_vzmulbyconji_omp_offload_ilp64(n, a, inca, b, incb, y, incy)
        integer(kind=8), intent(in) :: n
        complex(kind=8), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        complex(kind=8), intent(in) :: b(*)
        integer(kind=8), intent(in) :: incb
        complex(kind=8), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        end subroutine mkl_vm_vzmulbyconji_omp_offload_ilp64

        subroutine mkl_vm_vmzmulbyconji_omp_offload_ilp64(n, a, inca, b, incb, y, incy, mode)
        integer(kind=8), intent(in) :: n
        complex(kind=8), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        complex(kind=8), intent(in) :: b(*)
        integer(kind=8), intent(in) :: incb
        complex(kind=8), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmzmulbyconji_omp_offload_ilp64




! function: NearbyInt, indexing: simple
        subroutine mkl_vm_vsnearbyint_omp_offload_ilp64(n, a, y)
        integer(kind=8), intent(in) :: n
        real(kind=4), intent(in) :: a(*)
        real(kind=4), intent(inout) :: y(*)
        end subroutine mkl_vm_vsnearbyint_omp_offload_ilp64

        subroutine mkl_vm_vmsnearbyint_omp_offload_ilp64(n, a, y, mode)
        integer(kind=8), intent(in) :: n
        real(kind=4), intent(in) :: a(*)
        real(kind=4), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmsnearbyint_omp_offload_ilp64

        subroutine mkl_vm_vdnearbyint_omp_offload_ilp64(n, a, y)
        integer(kind=8), intent(in) :: n
        real(kind=8), intent(in) :: a(*)
        real(kind=8), intent(inout) :: y(*)
        end subroutine mkl_vm_vdnearbyint_omp_offload_ilp64

        subroutine mkl_vm_vmdnearbyint_omp_offload_ilp64(n, a, y, mode)
        integer(kind=8), intent(in) :: n
        real(kind=8), intent(in) :: a(*)
        real(kind=8), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmdnearbyint_omp_offload_ilp64


! function: NearbyInt, indexing: strided
        subroutine mkl_vm_vsnearbyinti_omp_offload_ilp64(n, a, inca, y, incy)
        integer(kind=8), intent(in) :: n
        real(kind=4), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        real(kind=4), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        end subroutine mkl_vm_vsnearbyinti_omp_offload_ilp64

        subroutine mkl_vm_vmsnearbyinti_omp_offload_ilp64(n, a, inca, y, incy, mode)
        integer(kind=8), intent(in) :: n
        real(kind=4), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        real(kind=4), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmsnearbyinti_omp_offload_ilp64

        subroutine mkl_vm_vdnearbyinti_omp_offload_ilp64(n, a, inca, y, incy)
        integer(kind=8), intent(in) :: n
        real(kind=8), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        real(kind=8), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        end subroutine mkl_vm_vdnearbyinti_omp_offload_ilp64

        subroutine mkl_vm_vmdnearbyinti_omp_offload_ilp64(n, a, inca, y, incy, mode)
        integer(kind=8), intent(in) :: n
        real(kind=8), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        real(kind=8), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmdnearbyinti_omp_offload_ilp64




! function: NextAfter, indexing: simple
        subroutine mkl_vm_vsnextafter_omp_offload_ilp64(n, a, b, y)
        integer(kind=8), intent(in) :: n
        real(kind=4), intent(in) :: a(*)
        real(kind=4), intent(in) :: b(*)
        real(kind=4), intent(inout) :: y(*)
        end subroutine mkl_vm_vsnextafter_omp_offload_ilp64

        subroutine mkl_vm_vmsnextafter_omp_offload_ilp64(n, a, b, y, mode)
        integer(kind=8), intent(in) :: n
        real(kind=4), intent(in) :: a(*)
        real(kind=4), intent(in) :: b(*)
        real(kind=4), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmsnextafter_omp_offload_ilp64

        subroutine mkl_vm_vdnextafter_omp_offload_ilp64(n, a, b, y)
        integer(kind=8), intent(in) :: n
        real(kind=8), intent(in) :: a(*)
        real(kind=8), intent(in) :: b(*)
        real(kind=8), intent(inout) :: y(*)
        end subroutine mkl_vm_vdnextafter_omp_offload_ilp64

        subroutine mkl_vm_vmdnextafter_omp_offload_ilp64(n, a, b, y, mode)
        integer(kind=8), intent(in) :: n
        real(kind=8), intent(in) :: a(*)
        real(kind=8), intent(in) :: b(*)
        real(kind=8), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmdnextafter_omp_offload_ilp64


! function: NextAfter, indexing: strided
        subroutine mkl_vm_vsnextafteri_omp_offload_ilp64(n, a, inca, b, incb, y, incy)
        integer(kind=8), intent(in) :: n
        real(kind=4), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        real(kind=4), intent(in) :: b(*)
        integer(kind=8), intent(in) :: incb
        real(kind=4), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        end subroutine mkl_vm_vsnextafteri_omp_offload_ilp64

        subroutine mkl_vm_vmsnextafteri_omp_offload_ilp64(n, a, inca, b, incb, y, incy, mode)
        integer(kind=8), intent(in) :: n
        real(kind=4), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        real(kind=4), intent(in) :: b(*)
        integer(kind=8), intent(in) :: incb
        real(kind=4), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmsnextafteri_omp_offload_ilp64

        subroutine mkl_vm_vdnextafteri_omp_offload_ilp64(n, a, inca, b, incb, y, incy)
        integer(kind=8), intent(in) :: n
        real(kind=8), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        real(kind=8), intent(in) :: b(*)
        integer(kind=8), intent(in) :: incb
        real(kind=8), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        end subroutine mkl_vm_vdnextafteri_omp_offload_ilp64

        subroutine mkl_vm_vmdnextafteri_omp_offload_ilp64(n, a, inca, b, incb, y, incy, mode)
        integer(kind=8), intent(in) :: n
        real(kind=8), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        real(kind=8), intent(in) :: b(*)
        integer(kind=8), intent(in) :: incb
        real(kind=8), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmdnextafteri_omp_offload_ilp64




! function: Pow, indexing: simple
        subroutine mkl_vm_vspow_omp_offload_ilp64(n, a, b, y)
        integer(kind=8), intent(in) :: n
        real(kind=4), intent(in) :: a(*)
        real(kind=4), intent(in) :: b(*)
        real(kind=4), intent(inout) :: y(*)
        end subroutine mkl_vm_vspow_omp_offload_ilp64

        subroutine mkl_vm_vmspow_omp_offload_ilp64(n, a, b, y, mode)
        integer(kind=8), intent(in) :: n
        real(kind=4), intent(in) :: a(*)
        real(kind=4), intent(in) :: b(*)
        real(kind=4), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmspow_omp_offload_ilp64

        subroutine mkl_vm_vdpow_omp_offload_ilp64(n, a, b, y)
        integer(kind=8), intent(in) :: n
        real(kind=8), intent(in) :: a(*)
        real(kind=8), intent(in) :: b(*)
        real(kind=8), intent(inout) :: y(*)
        end subroutine mkl_vm_vdpow_omp_offload_ilp64

        subroutine mkl_vm_vmdpow_omp_offload_ilp64(n, a, b, y, mode)
        integer(kind=8), intent(in) :: n
        real(kind=8), intent(in) :: a(*)
        real(kind=8), intent(in) :: b(*)
        real(kind=8), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmdpow_omp_offload_ilp64

        subroutine mkl_vm_vcpow_omp_offload_ilp64(n, a, b, y)
        integer(kind=8), intent(in) :: n
        complex(kind=4), intent(in) :: a(*)
        complex(kind=4), intent(in) :: b(*)
        complex(kind=4), intent(inout) :: y(*)
        end subroutine mkl_vm_vcpow_omp_offload_ilp64

        subroutine mkl_vm_vmcpow_omp_offload_ilp64(n, a, b, y, mode)
        integer(kind=8), intent(in) :: n
        complex(kind=4), intent(in) :: a(*)
        complex(kind=4), intent(in) :: b(*)
        complex(kind=4), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmcpow_omp_offload_ilp64

        subroutine mkl_vm_vzpow_omp_offload_ilp64(n, a, b, y)
        integer(kind=8), intent(in) :: n
        complex(kind=8), intent(in) :: a(*)
        complex(kind=8), intent(in) :: b(*)
        complex(kind=8), intent(inout) :: y(*)
        end subroutine mkl_vm_vzpow_omp_offload_ilp64

        subroutine mkl_vm_vmzpow_omp_offload_ilp64(n, a, b, y, mode)
        integer(kind=8), intent(in) :: n
        complex(kind=8), intent(in) :: a(*)
        complex(kind=8), intent(in) :: b(*)
        complex(kind=8), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmzpow_omp_offload_ilp64


! function: Pow, indexing: strided
        subroutine mkl_vm_vspowi_omp_offload_ilp64(n, a, inca, b, incb, y, incy)
        integer(kind=8), intent(in) :: n
        real(kind=4), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        real(kind=4), intent(in) :: b(*)
        integer(kind=8), intent(in) :: incb
        real(kind=4), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        end subroutine mkl_vm_vspowi_omp_offload_ilp64

        subroutine mkl_vm_vmspowi_omp_offload_ilp64(n, a, inca, b, incb, y, incy, mode)
        integer(kind=8), intent(in) :: n
        real(kind=4), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        real(kind=4), intent(in) :: b(*)
        integer(kind=8), intent(in) :: incb
        real(kind=4), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmspowi_omp_offload_ilp64

        subroutine mkl_vm_vdpowi_omp_offload_ilp64(n, a, inca, b, incb, y, incy)
        integer(kind=8), intent(in) :: n
        real(kind=8), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        real(kind=8), intent(in) :: b(*)
        integer(kind=8), intent(in) :: incb
        real(kind=8), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        end subroutine mkl_vm_vdpowi_omp_offload_ilp64

        subroutine mkl_vm_vmdpowi_omp_offload_ilp64(n, a, inca, b, incb, y, incy, mode)
        integer(kind=8), intent(in) :: n
        real(kind=8), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        real(kind=8), intent(in) :: b(*)
        integer(kind=8), intent(in) :: incb
        real(kind=8), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmdpowi_omp_offload_ilp64

        subroutine mkl_vm_vcpowi_omp_offload_ilp64(n, a, inca, b, incb, y, incy)
        integer(kind=8), intent(in) :: n
        complex(kind=4), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        complex(kind=4), intent(in) :: b(*)
        integer(kind=8), intent(in) :: incb
        complex(kind=4), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        end subroutine mkl_vm_vcpowi_omp_offload_ilp64

        subroutine mkl_vm_vmcpowi_omp_offload_ilp64(n, a, inca, b, incb, y, incy, mode)
        integer(kind=8), intent(in) :: n
        complex(kind=4), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        complex(kind=4), intent(in) :: b(*)
        integer(kind=8), intent(in) :: incb
        complex(kind=4), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmcpowi_omp_offload_ilp64

        subroutine mkl_vm_vzpowi_omp_offload_ilp64(n, a, inca, b, incb, y, incy)
        integer(kind=8), intent(in) :: n
        complex(kind=8), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        complex(kind=8), intent(in) :: b(*)
        integer(kind=8), intent(in) :: incb
        complex(kind=8), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        end subroutine mkl_vm_vzpowi_omp_offload_ilp64

        subroutine mkl_vm_vmzpowi_omp_offload_ilp64(n, a, inca, b, incb, y, incy, mode)
        integer(kind=8), intent(in) :: n
        complex(kind=8), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        complex(kind=8), intent(in) :: b(*)
        integer(kind=8), intent(in) :: incb
        complex(kind=8), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmzpowi_omp_offload_ilp64




! function: Pow2o3, indexing: simple
        subroutine mkl_vm_vspow2o3_omp_offload_ilp64(n, a, y)
        integer(kind=8), intent(in) :: n
        real(kind=4), intent(in) :: a(*)
        real(kind=4), intent(inout) :: y(*)
        end subroutine mkl_vm_vspow2o3_omp_offload_ilp64

        subroutine mkl_vm_vmspow2o3_omp_offload_ilp64(n, a, y, mode)
        integer(kind=8), intent(in) :: n
        real(kind=4), intent(in) :: a(*)
        real(kind=4), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmspow2o3_omp_offload_ilp64

        subroutine mkl_vm_vdpow2o3_omp_offload_ilp64(n, a, y)
        integer(kind=8), intent(in) :: n
        real(kind=8), intent(in) :: a(*)
        real(kind=8), intent(inout) :: y(*)
        end subroutine mkl_vm_vdpow2o3_omp_offload_ilp64

        subroutine mkl_vm_vmdpow2o3_omp_offload_ilp64(n, a, y, mode)
        integer(kind=8), intent(in) :: n
        real(kind=8), intent(in) :: a(*)
        real(kind=8), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmdpow2o3_omp_offload_ilp64


! function: Pow2o3, indexing: strided
        subroutine mkl_vm_vspow2o3i_omp_offload_ilp64(n, a, inca, y, incy)
        integer(kind=8), intent(in) :: n
        real(kind=4), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        real(kind=4), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        end subroutine mkl_vm_vspow2o3i_omp_offload_ilp64

        subroutine mkl_vm_vmspow2o3i_omp_offload_ilp64(n, a, inca, y, incy, mode)
        integer(kind=8), intent(in) :: n
        real(kind=4), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        real(kind=4), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmspow2o3i_omp_offload_ilp64

        subroutine mkl_vm_vdpow2o3i_omp_offload_ilp64(n, a, inca, y, incy)
        integer(kind=8), intent(in) :: n
        real(kind=8), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        real(kind=8), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        end subroutine mkl_vm_vdpow2o3i_omp_offload_ilp64

        subroutine mkl_vm_vmdpow2o3i_omp_offload_ilp64(n, a, inca, y, incy, mode)
        integer(kind=8), intent(in) :: n
        real(kind=8), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        real(kind=8), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmdpow2o3i_omp_offload_ilp64




! function: Pow3o2, indexing: simple
        subroutine mkl_vm_vspow3o2_omp_offload_ilp64(n, a, y)
        integer(kind=8), intent(in) :: n
        real(kind=4), intent(in) :: a(*)
        real(kind=4), intent(inout) :: y(*)
        end subroutine mkl_vm_vspow3o2_omp_offload_ilp64

        subroutine mkl_vm_vmspow3o2_omp_offload_ilp64(n, a, y, mode)
        integer(kind=8), intent(in) :: n
        real(kind=4), intent(in) :: a(*)
        real(kind=4), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmspow3o2_omp_offload_ilp64

        subroutine mkl_vm_vdpow3o2_omp_offload_ilp64(n, a, y)
        integer(kind=8), intent(in) :: n
        real(kind=8), intent(in) :: a(*)
        real(kind=8), intent(inout) :: y(*)
        end subroutine mkl_vm_vdpow3o2_omp_offload_ilp64

        subroutine mkl_vm_vmdpow3o2_omp_offload_ilp64(n, a, y, mode)
        integer(kind=8), intent(in) :: n
        real(kind=8), intent(in) :: a(*)
        real(kind=8), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmdpow3o2_omp_offload_ilp64


! function: Pow3o2, indexing: strided
        subroutine mkl_vm_vspow3o2i_omp_offload_ilp64(n, a, inca, y, incy)
        integer(kind=8), intent(in) :: n
        real(kind=4), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        real(kind=4), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        end subroutine mkl_vm_vspow3o2i_omp_offload_ilp64

        subroutine mkl_vm_vmspow3o2i_omp_offload_ilp64(n, a, inca, y, incy, mode)
        integer(kind=8), intent(in) :: n
        real(kind=4), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        real(kind=4), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmspow3o2i_omp_offload_ilp64

        subroutine mkl_vm_vdpow3o2i_omp_offload_ilp64(n, a, inca, y, incy)
        integer(kind=8), intent(in) :: n
        real(kind=8), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        real(kind=8), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        end subroutine mkl_vm_vdpow3o2i_omp_offload_ilp64

        subroutine mkl_vm_vmdpow3o2i_omp_offload_ilp64(n, a, inca, y, incy, mode)
        integer(kind=8), intent(in) :: n
        real(kind=8), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        real(kind=8), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmdpow3o2i_omp_offload_ilp64




! function: Powr, indexing: simple
        subroutine mkl_vm_vspowr_omp_offload_ilp64(n, a, b, y)
        integer(kind=8), intent(in) :: n
        real(kind=4), intent(in) :: a(*)
        real(kind=4), intent(in) :: b(*)
        real(kind=4), intent(inout) :: y(*)
        end subroutine mkl_vm_vspowr_omp_offload_ilp64

        subroutine mkl_vm_vmspowr_omp_offload_ilp64(n, a, b, y, mode)
        integer(kind=8), intent(in) :: n
        real(kind=4), intent(in) :: a(*)
        real(kind=4), intent(in) :: b(*)
        real(kind=4), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmspowr_omp_offload_ilp64

        subroutine mkl_vm_vdpowr_omp_offload_ilp64(n, a, b, y)
        integer(kind=8), intent(in) :: n
        real(kind=8), intent(in) :: a(*)
        real(kind=8), intent(in) :: b(*)
        real(kind=8), intent(inout) :: y(*)
        end subroutine mkl_vm_vdpowr_omp_offload_ilp64

        subroutine mkl_vm_vmdpowr_omp_offload_ilp64(n, a, b, y, mode)
        integer(kind=8), intent(in) :: n
        real(kind=8), intent(in) :: a(*)
        real(kind=8), intent(in) :: b(*)
        real(kind=8), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmdpowr_omp_offload_ilp64


! function: Powr, indexing: strided
        subroutine mkl_vm_vspowri_omp_offload_ilp64(n, a, inca, b, incb, y, incy)
        integer(kind=8), intent(in) :: n
        real(kind=4), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        real(kind=4), intent(in) :: b(*)
        integer(kind=8), intent(in) :: incb
        real(kind=4), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        end subroutine mkl_vm_vspowri_omp_offload_ilp64

        subroutine mkl_vm_vmspowri_omp_offload_ilp64(n, a, inca, b, incb, y, incy, mode)
        integer(kind=8), intent(in) :: n
        real(kind=4), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        real(kind=4), intent(in) :: b(*)
        integer(kind=8), intent(in) :: incb
        real(kind=4), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmspowri_omp_offload_ilp64

        subroutine mkl_vm_vdpowri_omp_offload_ilp64(n, a, inca, b, incb, y, incy)
        integer(kind=8), intent(in) :: n
        real(kind=8), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        real(kind=8), intent(in) :: b(*)
        integer(kind=8), intent(in) :: incb
        real(kind=8), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        end subroutine mkl_vm_vdpowri_omp_offload_ilp64

        subroutine mkl_vm_vmdpowri_omp_offload_ilp64(n, a, inca, b, incb, y, incy, mode)
        integer(kind=8), intent(in) :: n
        real(kind=8), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        real(kind=8), intent(in) :: b(*)
        integer(kind=8), intent(in) :: incb
        real(kind=8), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmdpowri_omp_offload_ilp64




! function: Powx, indexing: simple
        subroutine mkl_vm_vspowx_omp_offload_ilp64(n, a, b, y)
        integer(kind=8), intent(in) :: n
        real(kind=4), intent(in) :: a(*)
        real(kind=4), intent(in) :: b
        real(kind=4), intent(inout) :: y(*)
        end subroutine mkl_vm_vspowx_omp_offload_ilp64

        subroutine mkl_vm_vmspowx_omp_offload_ilp64(n, a, b, y, mode)
        integer(kind=8), intent(in) :: n
        real(kind=4), intent(in) :: a(*)
        real(kind=4), intent(in) :: b
        real(kind=4), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmspowx_omp_offload_ilp64

        subroutine mkl_vm_vdpowx_omp_offload_ilp64(n, a, b, y)
        integer(kind=8), intent(in) :: n
        real(kind=8), intent(in) :: a(*)
        real(kind=8), intent(in) :: b
        real(kind=8), intent(inout) :: y(*)
        end subroutine mkl_vm_vdpowx_omp_offload_ilp64

        subroutine mkl_vm_vmdpowx_omp_offload_ilp64(n, a, b, y, mode)
        integer(kind=8), intent(in) :: n
        real(kind=8), intent(in) :: a(*)
        real(kind=8), intent(in) :: b
        real(kind=8), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmdpowx_omp_offload_ilp64

        subroutine mkl_vm_vcpowx_omp_offload_ilp64(n, a, b, y)
        integer(kind=8), intent(in) :: n
        complex(kind=4), intent(in) :: a(*)
        complex(kind=4), intent(in) :: b
        complex(kind=4), intent(inout) :: y(*)
        end subroutine mkl_vm_vcpowx_omp_offload_ilp64

        subroutine mkl_vm_vmcpowx_omp_offload_ilp64(n, a, b, y, mode)
        integer(kind=8), intent(in) :: n
        complex(kind=4), intent(in) :: a(*)
        complex(kind=4), intent(in) :: b
        complex(kind=4), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmcpowx_omp_offload_ilp64

        subroutine mkl_vm_vzpowx_omp_offload_ilp64(n, a, b, y)
        integer(kind=8), intent(in) :: n
        complex(kind=8), intent(in) :: a(*)
        complex(kind=8), intent(in) :: b
        complex(kind=8), intent(inout) :: y(*)
        end subroutine mkl_vm_vzpowx_omp_offload_ilp64

        subroutine mkl_vm_vmzpowx_omp_offload_ilp64(n, a, b, y, mode)
        integer(kind=8), intent(in) :: n
        complex(kind=8), intent(in) :: a(*)
        complex(kind=8), intent(in) :: b
        complex(kind=8), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmzpowx_omp_offload_ilp64


! function: Powx, indexing: strided
        subroutine mkl_vm_vspowxi_omp_offload_ilp64(n, a, inca, b, y, incy)
        integer(kind=8), intent(in) :: n
        real(kind=4), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        real(kind=4), intent(in) :: b
        real(kind=4), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        end subroutine mkl_vm_vspowxi_omp_offload_ilp64

        subroutine mkl_vm_vmspowxi_omp_offload_ilp64(n, a, inca, b, y, incy, mode)
        integer(kind=8), intent(in) :: n
        real(kind=4), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        real(kind=4), intent(in) :: b
        real(kind=4), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmspowxi_omp_offload_ilp64

        subroutine mkl_vm_vdpowxi_omp_offload_ilp64(n, a, inca, b, y, incy)
        integer(kind=8), intent(in) :: n
        real(kind=8), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        real(kind=8), intent(in) :: b
        real(kind=8), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        end subroutine mkl_vm_vdpowxi_omp_offload_ilp64

        subroutine mkl_vm_vmdpowxi_omp_offload_ilp64(n, a, inca, b, y, incy, mode)
        integer(kind=8), intent(in) :: n
        real(kind=8), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        real(kind=8), intent(in) :: b
        real(kind=8), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmdpowxi_omp_offload_ilp64

        subroutine mkl_vm_vcpowxi_omp_offload_ilp64(n, a, inca, b, y, incy)
        integer(kind=8), intent(in) :: n
        complex(kind=4), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        complex(kind=4), intent(in) :: b
        complex(kind=4), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        end subroutine mkl_vm_vcpowxi_omp_offload_ilp64

        subroutine mkl_vm_vmcpowxi_omp_offload_ilp64(n, a, inca, b, y, incy, mode)
        integer(kind=8), intent(in) :: n
        complex(kind=4), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        complex(kind=4), intent(in) :: b
        complex(kind=4), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmcpowxi_omp_offload_ilp64

        subroutine mkl_vm_vzpowxi_omp_offload_ilp64(n, a, inca, b, y, incy)
        integer(kind=8), intent(in) :: n
        complex(kind=8), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        complex(kind=8), intent(in) :: b
        complex(kind=8), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        end subroutine mkl_vm_vzpowxi_omp_offload_ilp64

        subroutine mkl_vm_vmzpowxi_omp_offload_ilp64(n, a, inca, b, y, incy, mode)
        integer(kind=8), intent(in) :: n
        complex(kind=8), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        complex(kind=8), intent(in) :: b
        complex(kind=8), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmzpowxi_omp_offload_ilp64




! function: Remainder, indexing: simple
        subroutine mkl_vm_vsremainder_omp_offload_ilp64(n, a, b, y)
        integer(kind=8), intent(in) :: n
        real(kind=4), intent(in) :: a(*)
        real(kind=4), intent(in) :: b(*)
        real(kind=4), intent(inout) :: y(*)
        end subroutine mkl_vm_vsremainder_omp_offload_ilp64

        subroutine mkl_vm_vmsremainder_omp_offload_ilp64(n, a, b, y, mode)
        integer(kind=8), intent(in) :: n
        real(kind=4), intent(in) :: a(*)
        real(kind=4), intent(in) :: b(*)
        real(kind=4), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmsremainder_omp_offload_ilp64

        subroutine mkl_vm_vdremainder_omp_offload_ilp64(n, a, b, y)
        integer(kind=8), intent(in) :: n
        real(kind=8), intent(in) :: a(*)
        real(kind=8), intent(in) :: b(*)
        real(kind=8), intent(inout) :: y(*)
        end subroutine mkl_vm_vdremainder_omp_offload_ilp64

        subroutine mkl_vm_vmdremainder_omp_offload_ilp64(n, a, b, y, mode)
        integer(kind=8), intent(in) :: n
        real(kind=8), intent(in) :: a(*)
        real(kind=8), intent(in) :: b(*)
        real(kind=8), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmdremainder_omp_offload_ilp64


! function: Remainder, indexing: strided
        subroutine mkl_vm_vsremainderi_omp_offload_ilp64(n, a, inca, b, incb, y, incy)
        integer(kind=8), intent(in) :: n
        real(kind=4), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        real(kind=4), intent(in) :: b(*)
        integer(kind=8), intent(in) :: incb
        real(kind=4), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        end subroutine mkl_vm_vsremainderi_omp_offload_ilp64

        subroutine mkl_vm_vmsremainderi_omp_offload_ilp64(n, a, inca, b, incb, y, incy, mode)
        integer(kind=8), intent(in) :: n
        real(kind=4), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        real(kind=4), intent(in) :: b(*)
        integer(kind=8), intent(in) :: incb
        real(kind=4), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmsremainderi_omp_offload_ilp64

        subroutine mkl_vm_vdremainderi_omp_offload_ilp64(n, a, inca, b, incb, y, incy)
        integer(kind=8), intent(in) :: n
        real(kind=8), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        real(kind=8), intent(in) :: b(*)
        integer(kind=8), intent(in) :: incb
        real(kind=8), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        end subroutine mkl_vm_vdremainderi_omp_offload_ilp64

        subroutine mkl_vm_vmdremainderi_omp_offload_ilp64(n, a, inca, b, incb, y, incy, mode)
        integer(kind=8), intent(in) :: n
        real(kind=8), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        real(kind=8), intent(in) :: b(*)
        integer(kind=8), intent(in) :: incb
        real(kind=8), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmdremainderi_omp_offload_ilp64




! function: Rint, indexing: simple
        subroutine mkl_vm_vsrint_omp_offload_ilp64(n, a, y)
        integer(kind=8), intent(in) :: n
        real(kind=4), intent(in) :: a(*)
        real(kind=4), intent(inout) :: y(*)
        end subroutine mkl_vm_vsrint_omp_offload_ilp64

        subroutine mkl_vm_vmsrint_omp_offload_ilp64(n, a, y, mode)
        integer(kind=8), intent(in) :: n
        real(kind=4), intent(in) :: a(*)
        real(kind=4), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmsrint_omp_offload_ilp64

        subroutine mkl_vm_vdrint_omp_offload_ilp64(n, a, y)
        integer(kind=8), intent(in) :: n
        real(kind=8), intent(in) :: a(*)
        real(kind=8), intent(inout) :: y(*)
        end subroutine mkl_vm_vdrint_omp_offload_ilp64

        subroutine mkl_vm_vmdrint_omp_offload_ilp64(n, a, y, mode)
        integer(kind=8), intent(in) :: n
        real(kind=8), intent(in) :: a(*)
        real(kind=8), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmdrint_omp_offload_ilp64


! function: Rint, indexing: strided
        subroutine mkl_vm_vsrinti_omp_offload_ilp64(n, a, inca, y, incy)
        integer(kind=8), intent(in) :: n
        real(kind=4), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        real(kind=4), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        end subroutine mkl_vm_vsrinti_omp_offload_ilp64

        subroutine mkl_vm_vmsrinti_omp_offload_ilp64(n, a, inca, y, incy, mode)
        integer(kind=8), intent(in) :: n
        real(kind=4), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        real(kind=4), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmsrinti_omp_offload_ilp64

        subroutine mkl_vm_vdrinti_omp_offload_ilp64(n, a, inca, y, incy)
        integer(kind=8), intent(in) :: n
        real(kind=8), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        real(kind=8), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        end subroutine mkl_vm_vdrinti_omp_offload_ilp64

        subroutine mkl_vm_vmdrinti_omp_offload_ilp64(n, a, inca, y, incy, mode)
        integer(kind=8), intent(in) :: n
        real(kind=8), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        real(kind=8), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmdrinti_omp_offload_ilp64




! function: Round, indexing: simple
        subroutine mkl_vm_vsround_omp_offload_ilp64(n, a, y)
        integer(kind=8), intent(in) :: n
        real(kind=4), intent(in) :: a(*)
        real(kind=4), intent(inout) :: y(*)
        end subroutine mkl_vm_vsround_omp_offload_ilp64

        subroutine mkl_vm_vmsround_omp_offload_ilp64(n, a, y, mode)
        integer(kind=8), intent(in) :: n
        real(kind=4), intent(in) :: a(*)
        real(kind=4), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmsround_omp_offload_ilp64

        subroutine mkl_vm_vdround_omp_offload_ilp64(n, a, y)
        integer(kind=8), intent(in) :: n
        real(kind=8), intent(in) :: a(*)
        real(kind=8), intent(inout) :: y(*)
        end subroutine mkl_vm_vdround_omp_offload_ilp64

        subroutine mkl_vm_vmdround_omp_offload_ilp64(n, a, y, mode)
        integer(kind=8), intent(in) :: n
        real(kind=8), intent(in) :: a(*)
        real(kind=8), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmdround_omp_offload_ilp64


! function: Round, indexing: strided
        subroutine mkl_vm_vsroundi_omp_offload_ilp64(n, a, inca, y, incy)
        integer(kind=8), intent(in) :: n
        real(kind=4), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        real(kind=4), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        end subroutine mkl_vm_vsroundi_omp_offload_ilp64

        subroutine mkl_vm_vmsroundi_omp_offload_ilp64(n, a, inca, y, incy, mode)
        integer(kind=8), intent(in) :: n
        real(kind=4), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        real(kind=4), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmsroundi_omp_offload_ilp64

        subroutine mkl_vm_vdroundi_omp_offload_ilp64(n, a, inca, y, incy)
        integer(kind=8), intent(in) :: n
        real(kind=8), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        real(kind=8), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        end subroutine mkl_vm_vdroundi_omp_offload_ilp64

        subroutine mkl_vm_vmdroundi_omp_offload_ilp64(n, a, inca, y, incy, mode)
        integer(kind=8), intent(in) :: n
        real(kind=8), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        real(kind=8), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmdroundi_omp_offload_ilp64




! function: Sin, indexing: simple
        subroutine mkl_vm_vssin_omp_offload_ilp64(n, a, y)
        integer(kind=8), intent(in) :: n
        real(kind=4), intent(in) :: a(*)
        real(kind=4), intent(inout) :: y(*)
        end subroutine mkl_vm_vssin_omp_offload_ilp64

        subroutine mkl_vm_vmssin_omp_offload_ilp64(n, a, y, mode)
        integer(kind=8), intent(in) :: n
        real(kind=4), intent(in) :: a(*)
        real(kind=4), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmssin_omp_offload_ilp64

        subroutine mkl_vm_vdsin_omp_offload_ilp64(n, a, y)
        integer(kind=8), intent(in) :: n
        real(kind=8), intent(in) :: a(*)
        real(kind=8), intent(inout) :: y(*)
        end subroutine mkl_vm_vdsin_omp_offload_ilp64

        subroutine mkl_vm_vmdsin_omp_offload_ilp64(n, a, y, mode)
        integer(kind=8), intent(in) :: n
        real(kind=8), intent(in) :: a(*)
        real(kind=8), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmdsin_omp_offload_ilp64

        subroutine mkl_vm_vcsin_omp_offload_ilp64(n, a, y)
        integer(kind=8), intent(in) :: n
        complex(kind=4), intent(in) :: a(*)
        complex(kind=4), intent(inout) :: y(*)
        end subroutine mkl_vm_vcsin_omp_offload_ilp64

        subroutine mkl_vm_vmcsin_omp_offload_ilp64(n, a, y, mode)
        integer(kind=8), intent(in) :: n
        complex(kind=4), intent(in) :: a(*)
        complex(kind=4), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmcsin_omp_offload_ilp64

        subroutine mkl_vm_vzsin_omp_offload_ilp64(n, a, y)
        integer(kind=8), intent(in) :: n
        complex(kind=8), intent(in) :: a(*)
        complex(kind=8), intent(inout) :: y(*)
        end subroutine mkl_vm_vzsin_omp_offload_ilp64

        subroutine mkl_vm_vmzsin_omp_offload_ilp64(n, a, y, mode)
        integer(kind=8), intent(in) :: n
        complex(kind=8), intent(in) :: a(*)
        complex(kind=8), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmzsin_omp_offload_ilp64


! function: Sin, indexing: strided
        subroutine mkl_vm_vssini_omp_offload_ilp64(n, a, inca, y, incy)
        integer(kind=8), intent(in) :: n
        real(kind=4), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        real(kind=4), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        end subroutine mkl_vm_vssini_omp_offload_ilp64

        subroutine mkl_vm_vmssini_omp_offload_ilp64(n, a, inca, y, incy, mode)
        integer(kind=8), intent(in) :: n
        real(kind=4), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        real(kind=4), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmssini_omp_offload_ilp64

        subroutine mkl_vm_vdsini_omp_offload_ilp64(n, a, inca, y, incy)
        integer(kind=8), intent(in) :: n
        real(kind=8), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        real(kind=8), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        end subroutine mkl_vm_vdsini_omp_offload_ilp64

        subroutine mkl_vm_vmdsini_omp_offload_ilp64(n, a, inca, y, incy, mode)
        integer(kind=8), intent(in) :: n
        real(kind=8), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        real(kind=8), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmdsini_omp_offload_ilp64

        subroutine mkl_vm_vcsini_omp_offload_ilp64(n, a, inca, y, incy)
        integer(kind=8), intent(in) :: n
        complex(kind=4), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        complex(kind=4), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        end subroutine mkl_vm_vcsini_omp_offload_ilp64

        subroutine mkl_vm_vmcsini_omp_offload_ilp64(n, a, inca, y, incy, mode)
        integer(kind=8), intent(in) :: n
        complex(kind=4), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        complex(kind=4), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmcsini_omp_offload_ilp64

        subroutine mkl_vm_vzsini_omp_offload_ilp64(n, a, inca, y, incy)
        integer(kind=8), intent(in) :: n
        complex(kind=8), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        complex(kind=8), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        end subroutine mkl_vm_vzsini_omp_offload_ilp64

        subroutine mkl_vm_vmzsini_omp_offload_ilp64(n, a, inca, y, incy, mode)
        integer(kind=8), intent(in) :: n
        complex(kind=8), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        complex(kind=8), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmzsini_omp_offload_ilp64




! function: SinCos, indexing: simple
        subroutine mkl_vm_vssincos_omp_offload_ilp64(n, a, y, z)
        integer(kind=8), intent(in) :: n
        real(kind=4), intent(in) :: a(*)
        real(kind=4), intent(inout) :: y(*)
        real(kind=4), intent(inout) :: z(*)
        end subroutine mkl_vm_vssincos_omp_offload_ilp64

        subroutine mkl_vm_vmssincos_omp_offload_ilp64(n, a, y, z, mode)
        integer(kind=8), intent(in) :: n
        real(kind=4), intent(in) :: a(*)
        real(kind=4), intent(inout) :: y(*)
        real(kind=4), intent(inout) :: z(*)
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmssincos_omp_offload_ilp64

        subroutine mkl_vm_vdsincos_omp_offload_ilp64(n, a, y, z)
        integer(kind=8), intent(in) :: n
        real(kind=8), intent(in) :: a(*)
        real(kind=8), intent(inout) :: y(*)
        real(kind=8), intent(inout) :: z(*)
        end subroutine mkl_vm_vdsincos_omp_offload_ilp64

        subroutine mkl_vm_vmdsincos_omp_offload_ilp64(n, a, y, z, mode)
        integer(kind=8), intent(in) :: n
        real(kind=8), intent(in) :: a(*)
        real(kind=8), intent(inout) :: y(*)
        real(kind=8), intent(inout) :: z(*)
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmdsincos_omp_offload_ilp64


! function: SinCos, indexing: strided
        subroutine mkl_vm_vssincosi_omp_offload_ilp64(n, a, inca, y, incy, z, incz)
        integer(kind=8), intent(in) :: n
        real(kind=4), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        real(kind=4), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        real(kind=4), intent(inout) :: z(*)
        integer(kind=8), intent(in) :: incz
        end subroutine mkl_vm_vssincosi_omp_offload_ilp64

        subroutine mkl_vm_vmssincosi_omp_offload_ilp64(n, a, inca, y, incy, z, incz, mode)
        integer(kind=8), intent(in) :: n
        real(kind=4), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        real(kind=4), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        real(kind=4), intent(inout) :: z(*)
        integer(kind=8), intent(in) :: incz
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmssincosi_omp_offload_ilp64

        subroutine mkl_vm_vdsincosi_omp_offload_ilp64(n, a, inca, y, incy, z, incz)
        integer(kind=8), intent(in) :: n
        real(kind=8), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        real(kind=8), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        real(kind=8), intent(inout) :: z(*)
        integer(kind=8), intent(in) :: incz
        end subroutine mkl_vm_vdsincosi_omp_offload_ilp64

        subroutine mkl_vm_vmdsincosi_omp_offload_ilp64(n, a, inca, y, incy, z, incz, mode)
        integer(kind=8), intent(in) :: n
        real(kind=8), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        real(kind=8), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        real(kind=8), intent(inout) :: z(*)
        integer(kind=8), intent(in) :: incz
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmdsincosi_omp_offload_ilp64




! function: Sind, indexing: simple
        subroutine mkl_vm_vssind_omp_offload_ilp64(n, a, y)
        integer(kind=8), intent(in) :: n
        real(kind=4), intent(in) :: a(*)
        real(kind=4), intent(inout) :: y(*)
        end subroutine mkl_vm_vssind_omp_offload_ilp64

        subroutine mkl_vm_vmssind_omp_offload_ilp64(n, a, y, mode)
        integer(kind=8), intent(in) :: n
        real(kind=4), intent(in) :: a(*)
        real(kind=4), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmssind_omp_offload_ilp64

        subroutine mkl_vm_vdsind_omp_offload_ilp64(n, a, y)
        integer(kind=8), intent(in) :: n
        real(kind=8), intent(in) :: a(*)
        real(kind=8), intent(inout) :: y(*)
        end subroutine mkl_vm_vdsind_omp_offload_ilp64

        subroutine mkl_vm_vmdsind_omp_offload_ilp64(n, a, y, mode)
        integer(kind=8), intent(in) :: n
        real(kind=8), intent(in) :: a(*)
        real(kind=8), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmdsind_omp_offload_ilp64


! function: Sind, indexing: strided
        subroutine mkl_vm_vssindi_omp_offload_ilp64(n, a, inca, y, incy)
        integer(kind=8), intent(in) :: n
        real(kind=4), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        real(kind=4), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        end subroutine mkl_vm_vssindi_omp_offload_ilp64

        subroutine mkl_vm_vmssindi_omp_offload_ilp64(n, a, inca, y, incy, mode)
        integer(kind=8), intent(in) :: n
        real(kind=4), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        real(kind=4), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmssindi_omp_offload_ilp64

        subroutine mkl_vm_vdsindi_omp_offload_ilp64(n, a, inca, y, incy)
        integer(kind=8), intent(in) :: n
        real(kind=8), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        real(kind=8), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        end subroutine mkl_vm_vdsindi_omp_offload_ilp64

        subroutine mkl_vm_vmdsindi_omp_offload_ilp64(n, a, inca, y, incy, mode)
        integer(kind=8), intent(in) :: n
        real(kind=8), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        real(kind=8), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmdsindi_omp_offload_ilp64




! function: Sinh, indexing: simple
        subroutine mkl_vm_vssinh_omp_offload_ilp64(n, a, y)
        integer(kind=8), intent(in) :: n
        real(kind=4), intent(in) :: a(*)
        real(kind=4), intent(inout) :: y(*)
        end subroutine mkl_vm_vssinh_omp_offload_ilp64

        subroutine mkl_vm_vmssinh_omp_offload_ilp64(n, a, y, mode)
        integer(kind=8), intent(in) :: n
        real(kind=4), intent(in) :: a(*)
        real(kind=4), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmssinh_omp_offload_ilp64

        subroutine mkl_vm_vdsinh_omp_offload_ilp64(n, a, y)
        integer(kind=8), intent(in) :: n
        real(kind=8), intent(in) :: a(*)
        real(kind=8), intent(inout) :: y(*)
        end subroutine mkl_vm_vdsinh_omp_offload_ilp64

        subroutine mkl_vm_vmdsinh_omp_offload_ilp64(n, a, y, mode)
        integer(kind=8), intent(in) :: n
        real(kind=8), intent(in) :: a(*)
        real(kind=8), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmdsinh_omp_offload_ilp64

        subroutine mkl_vm_vcsinh_omp_offload_ilp64(n, a, y)
        integer(kind=8), intent(in) :: n
        complex(kind=4), intent(in) :: a(*)
        complex(kind=4), intent(inout) :: y(*)
        end subroutine mkl_vm_vcsinh_omp_offload_ilp64

        subroutine mkl_vm_vmcsinh_omp_offload_ilp64(n, a, y, mode)
        integer(kind=8), intent(in) :: n
        complex(kind=4), intent(in) :: a(*)
        complex(kind=4), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmcsinh_omp_offload_ilp64

        subroutine mkl_vm_vzsinh_omp_offload_ilp64(n, a, y)
        integer(kind=8), intent(in) :: n
        complex(kind=8), intent(in) :: a(*)
        complex(kind=8), intent(inout) :: y(*)
        end subroutine mkl_vm_vzsinh_omp_offload_ilp64

        subroutine mkl_vm_vmzsinh_omp_offload_ilp64(n, a, y, mode)
        integer(kind=8), intent(in) :: n
        complex(kind=8), intent(in) :: a(*)
        complex(kind=8), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmzsinh_omp_offload_ilp64


! function: Sinh, indexing: strided
        subroutine mkl_vm_vssinhi_omp_offload_ilp64(n, a, inca, y, incy)
        integer(kind=8), intent(in) :: n
        real(kind=4), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        real(kind=4), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        end subroutine mkl_vm_vssinhi_omp_offload_ilp64

        subroutine mkl_vm_vmssinhi_omp_offload_ilp64(n, a, inca, y, incy, mode)
        integer(kind=8), intent(in) :: n
        real(kind=4), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        real(kind=4), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmssinhi_omp_offload_ilp64

        subroutine mkl_vm_vdsinhi_omp_offload_ilp64(n, a, inca, y, incy)
        integer(kind=8), intent(in) :: n
        real(kind=8), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        real(kind=8), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        end subroutine mkl_vm_vdsinhi_omp_offload_ilp64

        subroutine mkl_vm_vmdsinhi_omp_offload_ilp64(n, a, inca, y, incy, mode)
        integer(kind=8), intent(in) :: n
        real(kind=8), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        real(kind=8), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmdsinhi_omp_offload_ilp64

        subroutine mkl_vm_vcsinhi_omp_offload_ilp64(n, a, inca, y, incy)
        integer(kind=8), intent(in) :: n
        complex(kind=4), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        complex(kind=4), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        end subroutine mkl_vm_vcsinhi_omp_offload_ilp64

        subroutine mkl_vm_vmcsinhi_omp_offload_ilp64(n, a, inca, y, incy, mode)
        integer(kind=8), intent(in) :: n
        complex(kind=4), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        complex(kind=4), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmcsinhi_omp_offload_ilp64

        subroutine mkl_vm_vzsinhi_omp_offload_ilp64(n, a, inca, y, incy)
        integer(kind=8), intent(in) :: n
        complex(kind=8), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        complex(kind=8), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        end subroutine mkl_vm_vzsinhi_omp_offload_ilp64

        subroutine mkl_vm_vmzsinhi_omp_offload_ilp64(n, a, inca, y, incy, mode)
        integer(kind=8), intent(in) :: n
        complex(kind=8), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        complex(kind=8), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmzsinhi_omp_offload_ilp64




! function: Sinpi, indexing: simple
        subroutine mkl_vm_vssinpi_omp_offload_ilp64(n, a, y)
        integer(kind=8), intent(in) :: n
        real(kind=4), intent(in) :: a(*)
        real(kind=4), intent(inout) :: y(*)
        end subroutine mkl_vm_vssinpi_omp_offload_ilp64

        subroutine mkl_vm_vmssinpi_omp_offload_ilp64(n, a, y, mode)
        integer(kind=8), intent(in) :: n
        real(kind=4), intent(in) :: a(*)
        real(kind=4), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmssinpi_omp_offload_ilp64

        subroutine mkl_vm_vdsinpi_omp_offload_ilp64(n, a, y)
        integer(kind=8), intent(in) :: n
        real(kind=8), intent(in) :: a(*)
        real(kind=8), intent(inout) :: y(*)
        end subroutine mkl_vm_vdsinpi_omp_offload_ilp64

        subroutine mkl_vm_vmdsinpi_omp_offload_ilp64(n, a, y, mode)
        integer(kind=8), intent(in) :: n
        real(kind=8), intent(in) :: a(*)
        real(kind=8), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmdsinpi_omp_offload_ilp64


! function: Sinpi, indexing: strided
        subroutine mkl_vm_vssinpii_omp_offload_ilp64(n, a, inca, y, incy)
        integer(kind=8), intent(in) :: n
        real(kind=4), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        real(kind=4), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        end subroutine mkl_vm_vssinpii_omp_offload_ilp64

        subroutine mkl_vm_vmssinpii_omp_offload_ilp64(n, a, inca, y, incy, mode)
        integer(kind=8), intent(in) :: n
        real(kind=4), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        real(kind=4), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmssinpii_omp_offload_ilp64

        subroutine mkl_vm_vdsinpii_omp_offload_ilp64(n, a, inca, y, incy)
        integer(kind=8), intent(in) :: n
        real(kind=8), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        real(kind=8), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        end subroutine mkl_vm_vdsinpii_omp_offload_ilp64

        subroutine mkl_vm_vmdsinpii_omp_offload_ilp64(n, a, inca, y, incy, mode)
        integer(kind=8), intent(in) :: n
        real(kind=8), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        real(kind=8), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmdsinpii_omp_offload_ilp64




! function: Sqr, indexing: simple
        subroutine mkl_vm_vssqr_omp_offload_ilp64(n, a, y)
        integer(kind=8), intent(in) :: n
        real(kind=4), intent(in) :: a(*)
        real(kind=4), intent(inout) :: y(*)
        end subroutine mkl_vm_vssqr_omp_offload_ilp64

        subroutine mkl_vm_vmssqr_omp_offload_ilp64(n, a, y, mode)
        integer(kind=8), intent(in) :: n
        real(kind=4), intent(in) :: a(*)
        real(kind=4), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmssqr_omp_offload_ilp64

        subroutine mkl_vm_vdsqr_omp_offload_ilp64(n, a, y)
        integer(kind=8), intent(in) :: n
        real(kind=8), intent(in) :: a(*)
        real(kind=8), intent(inout) :: y(*)
        end subroutine mkl_vm_vdsqr_omp_offload_ilp64

        subroutine mkl_vm_vmdsqr_omp_offload_ilp64(n, a, y, mode)
        integer(kind=8), intent(in) :: n
        real(kind=8), intent(in) :: a(*)
        real(kind=8), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmdsqr_omp_offload_ilp64


! function: Sqr, indexing: strided
        subroutine mkl_vm_vssqri_omp_offload_ilp64(n, a, inca, y, incy)
        integer(kind=8), intent(in) :: n
        real(kind=4), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        real(kind=4), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        end subroutine mkl_vm_vssqri_omp_offload_ilp64

        subroutine mkl_vm_vmssqri_omp_offload_ilp64(n, a, inca, y, incy, mode)
        integer(kind=8), intent(in) :: n
        real(kind=4), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        real(kind=4), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmssqri_omp_offload_ilp64

        subroutine mkl_vm_vdsqri_omp_offload_ilp64(n, a, inca, y, incy)
        integer(kind=8), intent(in) :: n
        real(kind=8), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        real(kind=8), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        end subroutine mkl_vm_vdsqri_omp_offload_ilp64

        subroutine mkl_vm_vmdsqri_omp_offload_ilp64(n, a, inca, y, incy, mode)
        integer(kind=8), intent(in) :: n
        real(kind=8), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        real(kind=8), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmdsqri_omp_offload_ilp64




! function: Sqrt, indexing: simple
        subroutine mkl_vm_vssqrt_omp_offload_ilp64(n, a, y)
        integer(kind=8), intent(in) :: n
        real(kind=4), intent(in) :: a(*)
        real(kind=4), intent(inout) :: y(*)
        end subroutine mkl_vm_vssqrt_omp_offload_ilp64

        subroutine mkl_vm_vmssqrt_omp_offload_ilp64(n, a, y, mode)
        integer(kind=8), intent(in) :: n
        real(kind=4), intent(in) :: a(*)
        real(kind=4), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmssqrt_omp_offload_ilp64

        subroutine mkl_vm_vdsqrt_omp_offload_ilp64(n, a, y)
        integer(kind=8), intent(in) :: n
        real(kind=8), intent(in) :: a(*)
        real(kind=8), intent(inout) :: y(*)
        end subroutine mkl_vm_vdsqrt_omp_offload_ilp64

        subroutine mkl_vm_vmdsqrt_omp_offload_ilp64(n, a, y, mode)
        integer(kind=8), intent(in) :: n
        real(kind=8), intent(in) :: a(*)
        real(kind=8), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmdsqrt_omp_offload_ilp64

        subroutine mkl_vm_vcsqrt_omp_offload_ilp64(n, a, y)
        integer(kind=8), intent(in) :: n
        complex(kind=4), intent(in) :: a(*)
        complex(kind=4), intent(inout) :: y(*)
        end subroutine mkl_vm_vcsqrt_omp_offload_ilp64

        subroutine mkl_vm_vmcsqrt_omp_offload_ilp64(n, a, y, mode)
        integer(kind=8), intent(in) :: n
        complex(kind=4), intent(in) :: a(*)
        complex(kind=4), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmcsqrt_omp_offload_ilp64

        subroutine mkl_vm_vzsqrt_omp_offload_ilp64(n, a, y)
        integer(kind=8), intent(in) :: n
        complex(kind=8), intent(in) :: a(*)
        complex(kind=8), intent(inout) :: y(*)
        end subroutine mkl_vm_vzsqrt_omp_offload_ilp64

        subroutine mkl_vm_vmzsqrt_omp_offload_ilp64(n, a, y, mode)
        integer(kind=8), intent(in) :: n
        complex(kind=8), intent(in) :: a(*)
        complex(kind=8), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmzsqrt_omp_offload_ilp64


! function: Sqrt, indexing: strided
        subroutine mkl_vm_vssqrti_omp_offload_ilp64(n, a, inca, y, incy)
        integer(kind=8), intent(in) :: n
        real(kind=4), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        real(kind=4), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        end subroutine mkl_vm_vssqrti_omp_offload_ilp64

        subroutine mkl_vm_vmssqrti_omp_offload_ilp64(n, a, inca, y, incy, mode)
        integer(kind=8), intent(in) :: n
        real(kind=4), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        real(kind=4), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmssqrti_omp_offload_ilp64

        subroutine mkl_vm_vdsqrti_omp_offload_ilp64(n, a, inca, y, incy)
        integer(kind=8), intent(in) :: n
        real(kind=8), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        real(kind=8), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        end subroutine mkl_vm_vdsqrti_omp_offload_ilp64

        subroutine mkl_vm_vmdsqrti_omp_offload_ilp64(n, a, inca, y, incy, mode)
        integer(kind=8), intent(in) :: n
        real(kind=8), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        real(kind=8), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmdsqrti_omp_offload_ilp64

        subroutine mkl_vm_vcsqrti_omp_offload_ilp64(n, a, inca, y, incy)
        integer(kind=8), intent(in) :: n
        complex(kind=4), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        complex(kind=4), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        end subroutine mkl_vm_vcsqrti_omp_offload_ilp64

        subroutine mkl_vm_vmcsqrti_omp_offload_ilp64(n, a, inca, y, incy, mode)
        integer(kind=8), intent(in) :: n
        complex(kind=4), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        complex(kind=4), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmcsqrti_omp_offload_ilp64

        subroutine mkl_vm_vzsqrti_omp_offload_ilp64(n, a, inca, y, incy)
        integer(kind=8), intent(in) :: n
        complex(kind=8), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        complex(kind=8), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        end subroutine mkl_vm_vzsqrti_omp_offload_ilp64

        subroutine mkl_vm_vmzsqrti_omp_offload_ilp64(n, a, inca, y, incy, mode)
        integer(kind=8), intent(in) :: n
        complex(kind=8), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        complex(kind=8), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmzsqrti_omp_offload_ilp64




! function: Sub, indexing: simple
        subroutine mkl_vm_vssub_omp_offload_ilp64(n, a, b, y)
        integer(kind=8), intent(in) :: n
        real(kind=4), intent(in) :: a(*)
        real(kind=4), intent(in) :: b(*)
        real(kind=4), intent(inout) :: y(*)
        end subroutine mkl_vm_vssub_omp_offload_ilp64

        subroutine mkl_vm_vmssub_omp_offload_ilp64(n, a, b, y, mode)
        integer(kind=8), intent(in) :: n
        real(kind=4), intent(in) :: a(*)
        real(kind=4), intent(in) :: b(*)
        real(kind=4), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmssub_omp_offload_ilp64

        subroutine mkl_vm_vdsub_omp_offload_ilp64(n, a, b, y)
        integer(kind=8), intent(in) :: n
        real(kind=8), intent(in) :: a(*)
        real(kind=8), intent(in) :: b(*)
        real(kind=8), intent(inout) :: y(*)
        end subroutine mkl_vm_vdsub_omp_offload_ilp64

        subroutine mkl_vm_vmdsub_omp_offload_ilp64(n, a, b, y, mode)
        integer(kind=8), intent(in) :: n
        real(kind=8), intent(in) :: a(*)
        real(kind=8), intent(in) :: b(*)
        real(kind=8), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmdsub_omp_offload_ilp64

        subroutine mkl_vm_vcsub_omp_offload_ilp64(n, a, b, y)
        integer(kind=8), intent(in) :: n
        complex(kind=4), intent(in) :: a(*)
        complex(kind=4), intent(in) :: b(*)
        complex(kind=4), intent(inout) :: y(*)
        end subroutine mkl_vm_vcsub_omp_offload_ilp64

        subroutine mkl_vm_vmcsub_omp_offload_ilp64(n, a, b, y, mode)
        integer(kind=8), intent(in) :: n
        complex(kind=4), intent(in) :: a(*)
        complex(kind=4), intent(in) :: b(*)
        complex(kind=4), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmcsub_omp_offload_ilp64

        subroutine mkl_vm_vzsub_omp_offload_ilp64(n, a, b, y)
        integer(kind=8), intent(in) :: n
        complex(kind=8), intent(in) :: a(*)
        complex(kind=8), intent(in) :: b(*)
        complex(kind=8), intent(inout) :: y(*)
        end subroutine mkl_vm_vzsub_omp_offload_ilp64

        subroutine mkl_vm_vmzsub_omp_offload_ilp64(n, a, b, y, mode)
        integer(kind=8), intent(in) :: n
        complex(kind=8), intent(in) :: a(*)
        complex(kind=8), intent(in) :: b(*)
        complex(kind=8), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmzsub_omp_offload_ilp64


! function: Sub, indexing: strided
        subroutine mkl_vm_vssubi_omp_offload_ilp64(n, a, inca, b, incb, y, incy)
        integer(kind=8), intent(in) :: n
        real(kind=4), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        real(kind=4), intent(in) :: b(*)
        integer(kind=8), intent(in) :: incb
        real(kind=4), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        end subroutine mkl_vm_vssubi_omp_offload_ilp64

        subroutine mkl_vm_vmssubi_omp_offload_ilp64(n, a, inca, b, incb, y, incy, mode)
        integer(kind=8), intent(in) :: n
        real(kind=4), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        real(kind=4), intent(in) :: b(*)
        integer(kind=8), intent(in) :: incb
        real(kind=4), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmssubi_omp_offload_ilp64

        subroutine mkl_vm_vdsubi_omp_offload_ilp64(n, a, inca, b, incb, y, incy)
        integer(kind=8), intent(in) :: n
        real(kind=8), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        real(kind=8), intent(in) :: b(*)
        integer(kind=8), intent(in) :: incb
        real(kind=8), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        end subroutine mkl_vm_vdsubi_omp_offload_ilp64

        subroutine mkl_vm_vmdsubi_omp_offload_ilp64(n, a, inca, b, incb, y, incy, mode)
        integer(kind=8), intent(in) :: n
        real(kind=8), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        real(kind=8), intent(in) :: b(*)
        integer(kind=8), intent(in) :: incb
        real(kind=8), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmdsubi_omp_offload_ilp64

        subroutine mkl_vm_vcsubi_omp_offload_ilp64(n, a, inca, b, incb, y, incy)
        integer(kind=8), intent(in) :: n
        complex(kind=4), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        complex(kind=4), intent(in) :: b(*)
        integer(kind=8), intent(in) :: incb
        complex(kind=4), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        end subroutine mkl_vm_vcsubi_omp_offload_ilp64

        subroutine mkl_vm_vmcsubi_omp_offload_ilp64(n, a, inca, b, incb, y, incy, mode)
        integer(kind=8), intent(in) :: n
        complex(kind=4), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        complex(kind=4), intent(in) :: b(*)
        integer(kind=8), intent(in) :: incb
        complex(kind=4), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmcsubi_omp_offload_ilp64

        subroutine mkl_vm_vzsubi_omp_offload_ilp64(n, a, inca, b, incb, y, incy)
        integer(kind=8), intent(in) :: n
        complex(kind=8), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        complex(kind=8), intent(in) :: b(*)
        integer(kind=8), intent(in) :: incb
        complex(kind=8), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        end subroutine mkl_vm_vzsubi_omp_offload_ilp64

        subroutine mkl_vm_vmzsubi_omp_offload_ilp64(n, a, inca, b, incb, y, incy, mode)
        integer(kind=8), intent(in) :: n
        complex(kind=8), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        complex(kind=8), intent(in) :: b(*)
        integer(kind=8), intent(in) :: incb
        complex(kind=8), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmzsubi_omp_offload_ilp64




! function: Tan, indexing: simple
        subroutine mkl_vm_vstan_omp_offload_ilp64(n, a, y)
        integer(kind=8), intent(in) :: n
        real(kind=4), intent(in) :: a(*)
        real(kind=4), intent(inout) :: y(*)
        end subroutine mkl_vm_vstan_omp_offload_ilp64

        subroutine mkl_vm_vmstan_omp_offload_ilp64(n, a, y, mode)
        integer(kind=8), intent(in) :: n
        real(kind=4), intent(in) :: a(*)
        real(kind=4), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmstan_omp_offload_ilp64

        subroutine mkl_vm_vdtan_omp_offload_ilp64(n, a, y)
        integer(kind=8), intent(in) :: n
        real(kind=8), intent(in) :: a(*)
        real(kind=8), intent(inout) :: y(*)
        end subroutine mkl_vm_vdtan_omp_offload_ilp64

        subroutine mkl_vm_vmdtan_omp_offload_ilp64(n, a, y, mode)
        integer(kind=8), intent(in) :: n
        real(kind=8), intent(in) :: a(*)
        real(kind=8), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmdtan_omp_offload_ilp64

        subroutine mkl_vm_vctan_omp_offload_ilp64(n, a, y)
        integer(kind=8), intent(in) :: n
        complex(kind=4), intent(in) :: a(*)
        complex(kind=4), intent(inout) :: y(*)
        end subroutine mkl_vm_vctan_omp_offload_ilp64

        subroutine mkl_vm_vmctan_omp_offload_ilp64(n, a, y, mode)
        integer(kind=8), intent(in) :: n
        complex(kind=4), intent(in) :: a(*)
        complex(kind=4), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmctan_omp_offload_ilp64

        subroutine mkl_vm_vztan_omp_offload_ilp64(n, a, y)
        integer(kind=8), intent(in) :: n
        complex(kind=8), intent(in) :: a(*)
        complex(kind=8), intent(inout) :: y(*)
        end subroutine mkl_vm_vztan_omp_offload_ilp64

        subroutine mkl_vm_vmztan_omp_offload_ilp64(n, a, y, mode)
        integer(kind=8), intent(in) :: n
        complex(kind=8), intent(in) :: a(*)
        complex(kind=8), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmztan_omp_offload_ilp64


! function: Tan, indexing: strided
        subroutine mkl_vm_vstani_omp_offload_ilp64(n, a, inca, y, incy)
        integer(kind=8), intent(in) :: n
        real(kind=4), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        real(kind=4), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        end subroutine mkl_vm_vstani_omp_offload_ilp64

        subroutine mkl_vm_vmstani_omp_offload_ilp64(n, a, inca, y, incy, mode)
        integer(kind=8), intent(in) :: n
        real(kind=4), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        real(kind=4), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmstani_omp_offload_ilp64

        subroutine mkl_vm_vdtani_omp_offload_ilp64(n, a, inca, y, incy)
        integer(kind=8), intent(in) :: n
        real(kind=8), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        real(kind=8), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        end subroutine mkl_vm_vdtani_omp_offload_ilp64

        subroutine mkl_vm_vmdtani_omp_offload_ilp64(n, a, inca, y, incy, mode)
        integer(kind=8), intent(in) :: n
        real(kind=8), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        real(kind=8), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmdtani_omp_offload_ilp64

        subroutine mkl_vm_vctani_omp_offload_ilp64(n, a, inca, y, incy)
        integer(kind=8), intent(in) :: n
        complex(kind=4), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        complex(kind=4), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        end subroutine mkl_vm_vctani_omp_offload_ilp64

        subroutine mkl_vm_vmctani_omp_offload_ilp64(n, a, inca, y, incy, mode)
        integer(kind=8), intent(in) :: n
        complex(kind=4), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        complex(kind=4), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmctani_omp_offload_ilp64

        subroutine mkl_vm_vztani_omp_offload_ilp64(n, a, inca, y, incy)
        integer(kind=8), intent(in) :: n
        complex(kind=8), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        complex(kind=8), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        end subroutine mkl_vm_vztani_omp_offload_ilp64

        subroutine mkl_vm_vmztani_omp_offload_ilp64(n, a, inca, y, incy, mode)
        integer(kind=8), intent(in) :: n
        complex(kind=8), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        complex(kind=8), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmztani_omp_offload_ilp64




! function: Tand, indexing: simple
        subroutine mkl_vm_vstand_omp_offload_ilp64(n, a, y)
        integer(kind=8), intent(in) :: n
        real(kind=4), intent(in) :: a(*)
        real(kind=4), intent(inout) :: y(*)
        end subroutine mkl_vm_vstand_omp_offload_ilp64

        subroutine mkl_vm_vmstand_omp_offload_ilp64(n, a, y, mode)
        integer(kind=8), intent(in) :: n
        real(kind=4), intent(in) :: a(*)
        real(kind=4), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmstand_omp_offload_ilp64

        subroutine mkl_vm_vdtand_omp_offload_ilp64(n, a, y)
        integer(kind=8), intent(in) :: n
        real(kind=8), intent(in) :: a(*)
        real(kind=8), intent(inout) :: y(*)
        end subroutine mkl_vm_vdtand_omp_offload_ilp64

        subroutine mkl_vm_vmdtand_omp_offload_ilp64(n, a, y, mode)
        integer(kind=8), intent(in) :: n
        real(kind=8), intent(in) :: a(*)
        real(kind=8), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmdtand_omp_offload_ilp64


! function: Tand, indexing: strided
        subroutine mkl_vm_vstandi_omp_offload_ilp64(n, a, inca, y, incy)
        integer(kind=8), intent(in) :: n
        real(kind=4), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        real(kind=4), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        end subroutine mkl_vm_vstandi_omp_offload_ilp64

        subroutine mkl_vm_vmstandi_omp_offload_ilp64(n, a, inca, y, incy, mode)
        integer(kind=8), intent(in) :: n
        real(kind=4), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        real(kind=4), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmstandi_omp_offload_ilp64

        subroutine mkl_vm_vdtandi_omp_offload_ilp64(n, a, inca, y, incy)
        integer(kind=8), intent(in) :: n
        real(kind=8), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        real(kind=8), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        end subroutine mkl_vm_vdtandi_omp_offload_ilp64

        subroutine mkl_vm_vmdtandi_omp_offload_ilp64(n, a, inca, y, incy, mode)
        integer(kind=8), intent(in) :: n
        real(kind=8), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        real(kind=8), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmdtandi_omp_offload_ilp64




! function: Tanh, indexing: simple
        subroutine mkl_vm_vstanh_omp_offload_ilp64(n, a, y)
        integer(kind=8), intent(in) :: n
        real(kind=4), intent(in) :: a(*)
        real(kind=4), intent(inout) :: y(*)
        end subroutine mkl_vm_vstanh_omp_offload_ilp64

        subroutine mkl_vm_vmstanh_omp_offload_ilp64(n, a, y, mode)
        integer(kind=8), intent(in) :: n
        real(kind=4), intent(in) :: a(*)
        real(kind=4), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmstanh_omp_offload_ilp64

        subroutine mkl_vm_vdtanh_omp_offload_ilp64(n, a, y)
        integer(kind=8), intent(in) :: n
        real(kind=8), intent(in) :: a(*)
        real(kind=8), intent(inout) :: y(*)
        end subroutine mkl_vm_vdtanh_omp_offload_ilp64

        subroutine mkl_vm_vmdtanh_omp_offload_ilp64(n, a, y, mode)
        integer(kind=8), intent(in) :: n
        real(kind=8), intent(in) :: a(*)
        real(kind=8), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmdtanh_omp_offload_ilp64

        subroutine mkl_vm_vctanh_omp_offload_ilp64(n, a, y)
        integer(kind=8), intent(in) :: n
        complex(kind=4), intent(in) :: a(*)
        complex(kind=4), intent(inout) :: y(*)
        end subroutine mkl_vm_vctanh_omp_offload_ilp64

        subroutine mkl_vm_vmctanh_omp_offload_ilp64(n, a, y, mode)
        integer(kind=8), intent(in) :: n
        complex(kind=4), intent(in) :: a(*)
        complex(kind=4), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmctanh_omp_offload_ilp64

        subroutine mkl_vm_vztanh_omp_offload_ilp64(n, a, y)
        integer(kind=8), intent(in) :: n
        complex(kind=8), intent(in) :: a(*)
        complex(kind=8), intent(inout) :: y(*)
        end subroutine mkl_vm_vztanh_omp_offload_ilp64

        subroutine mkl_vm_vmztanh_omp_offload_ilp64(n, a, y, mode)
        integer(kind=8), intent(in) :: n
        complex(kind=8), intent(in) :: a(*)
        complex(kind=8), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmztanh_omp_offload_ilp64


! function: Tanh, indexing: strided
        subroutine mkl_vm_vstanhi_omp_offload_ilp64(n, a, inca, y, incy)
        integer(kind=8), intent(in) :: n
        real(kind=4), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        real(kind=4), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        end subroutine mkl_vm_vstanhi_omp_offload_ilp64

        subroutine mkl_vm_vmstanhi_omp_offload_ilp64(n, a, inca, y, incy, mode)
        integer(kind=8), intent(in) :: n
        real(kind=4), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        real(kind=4), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmstanhi_omp_offload_ilp64

        subroutine mkl_vm_vdtanhi_omp_offload_ilp64(n, a, inca, y, incy)
        integer(kind=8), intent(in) :: n
        real(kind=8), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        real(kind=8), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        end subroutine mkl_vm_vdtanhi_omp_offload_ilp64

        subroutine mkl_vm_vmdtanhi_omp_offload_ilp64(n, a, inca, y, incy, mode)
        integer(kind=8), intent(in) :: n
        real(kind=8), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        real(kind=8), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmdtanhi_omp_offload_ilp64

        subroutine mkl_vm_vctanhi_omp_offload_ilp64(n, a, inca, y, incy)
        integer(kind=8), intent(in) :: n
        complex(kind=4), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        complex(kind=4), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        end subroutine mkl_vm_vctanhi_omp_offload_ilp64

        subroutine mkl_vm_vmctanhi_omp_offload_ilp64(n, a, inca, y, incy, mode)
        integer(kind=8), intent(in) :: n
        complex(kind=4), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        complex(kind=4), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmctanhi_omp_offload_ilp64

        subroutine mkl_vm_vztanhi_omp_offload_ilp64(n, a, inca, y, incy)
        integer(kind=8), intent(in) :: n
        complex(kind=8), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        complex(kind=8), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        end subroutine mkl_vm_vztanhi_omp_offload_ilp64

        subroutine mkl_vm_vmztanhi_omp_offload_ilp64(n, a, inca, y, incy, mode)
        integer(kind=8), intent(in) :: n
        complex(kind=8), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        complex(kind=8), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmztanhi_omp_offload_ilp64




! function: Tanpi, indexing: simple
        subroutine mkl_vm_vstanpi_omp_offload_ilp64(n, a, y)
        integer(kind=8), intent(in) :: n
        real(kind=4), intent(in) :: a(*)
        real(kind=4), intent(inout) :: y(*)
        end subroutine mkl_vm_vstanpi_omp_offload_ilp64

        subroutine mkl_vm_vmstanpi_omp_offload_ilp64(n, a, y, mode)
        integer(kind=8), intent(in) :: n
        real(kind=4), intent(in) :: a(*)
        real(kind=4), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmstanpi_omp_offload_ilp64

        subroutine mkl_vm_vdtanpi_omp_offload_ilp64(n, a, y)
        integer(kind=8), intent(in) :: n
        real(kind=8), intent(in) :: a(*)
        real(kind=8), intent(inout) :: y(*)
        end subroutine mkl_vm_vdtanpi_omp_offload_ilp64

        subroutine mkl_vm_vmdtanpi_omp_offload_ilp64(n, a, y, mode)
        integer(kind=8), intent(in) :: n
        real(kind=8), intent(in) :: a(*)
        real(kind=8), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmdtanpi_omp_offload_ilp64


! function: Tanpi, indexing: strided
        subroutine mkl_vm_vstanpii_omp_offload_ilp64(n, a, inca, y, incy)
        integer(kind=8), intent(in) :: n
        real(kind=4), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        real(kind=4), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        end subroutine mkl_vm_vstanpii_omp_offload_ilp64

        subroutine mkl_vm_vmstanpii_omp_offload_ilp64(n, a, inca, y, incy, mode)
        integer(kind=8), intent(in) :: n
        real(kind=4), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        real(kind=4), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmstanpii_omp_offload_ilp64

        subroutine mkl_vm_vdtanpii_omp_offload_ilp64(n, a, inca, y, incy)
        integer(kind=8), intent(in) :: n
        real(kind=8), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        real(kind=8), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        end subroutine mkl_vm_vdtanpii_omp_offload_ilp64

        subroutine mkl_vm_vmdtanpii_omp_offload_ilp64(n, a, inca, y, incy, mode)
        integer(kind=8), intent(in) :: n
        real(kind=8), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        real(kind=8), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmdtanpii_omp_offload_ilp64




! function: TGamma, indexing: simple
        subroutine mkl_vm_vstgamma_omp_offload_ilp64(n, a, y)
        integer(kind=8), intent(in) :: n
        real(kind=4), intent(in) :: a(*)
        real(kind=4), intent(inout) :: y(*)
        end subroutine mkl_vm_vstgamma_omp_offload_ilp64

        subroutine mkl_vm_vmstgamma_omp_offload_ilp64(n, a, y, mode)
        integer(kind=8), intent(in) :: n
        real(kind=4), intent(in) :: a(*)
        real(kind=4), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmstgamma_omp_offload_ilp64

        subroutine mkl_vm_vdtgamma_omp_offload_ilp64(n, a, y)
        integer(kind=8), intent(in) :: n
        real(kind=8), intent(in) :: a(*)
        real(kind=8), intent(inout) :: y(*)
        end subroutine mkl_vm_vdtgamma_omp_offload_ilp64

        subroutine mkl_vm_vmdtgamma_omp_offload_ilp64(n, a, y, mode)
        integer(kind=8), intent(in) :: n
        real(kind=8), intent(in) :: a(*)
        real(kind=8), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmdtgamma_omp_offload_ilp64


! function: TGamma, indexing: strided
        subroutine mkl_vm_vstgammai_omp_offload_ilp64(n, a, inca, y, incy)
        integer(kind=8), intent(in) :: n
        real(kind=4), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        real(kind=4), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        end subroutine mkl_vm_vstgammai_omp_offload_ilp64

        subroutine mkl_vm_vmstgammai_omp_offload_ilp64(n, a, inca, y, incy, mode)
        integer(kind=8), intent(in) :: n
        real(kind=4), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        real(kind=4), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmstgammai_omp_offload_ilp64

        subroutine mkl_vm_vdtgammai_omp_offload_ilp64(n, a, inca, y, incy)
        integer(kind=8), intent(in) :: n
        real(kind=8), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        real(kind=8), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        end subroutine mkl_vm_vdtgammai_omp_offload_ilp64

        subroutine mkl_vm_vmdtgammai_omp_offload_ilp64(n, a, inca, y, incy, mode)
        integer(kind=8), intent(in) :: n
        real(kind=8), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        real(kind=8), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmdtgammai_omp_offload_ilp64




! function: Trunc, indexing: simple
        subroutine mkl_vm_vstrunc_omp_offload_ilp64(n, a, y)
        integer(kind=8), intent(in) :: n
        real(kind=4), intent(in) :: a(*)
        real(kind=4), intent(inout) :: y(*)
        end subroutine mkl_vm_vstrunc_omp_offload_ilp64

        subroutine mkl_vm_vmstrunc_omp_offload_ilp64(n, a, y, mode)
        integer(kind=8), intent(in) :: n
        real(kind=4), intent(in) :: a(*)
        real(kind=4), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmstrunc_omp_offload_ilp64

        subroutine mkl_vm_vdtrunc_omp_offload_ilp64(n, a, y)
        integer(kind=8), intent(in) :: n
        real(kind=8), intent(in) :: a(*)
        real(kind=8), intent(inout) :: y(*)
        end subroutine mkl_vm_vdtrunc_omp_offload_ilp64

        subroutine mkl_vm_vmdtrunc_omp_offload_ilp64(n, a, y, mode)
        integer(kind=8), intent(in) :: n
        real(kind=8), intent(in) :: a(*)
        real(kind=8), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmdtrunc_omp_offload_ilp64


! function: Trunc, indexing: strided
        subroutine mkl_vm_vstrunci_omp_offload_ilp64(n, a, inca, y, incy)
        integer(kind=8), intent(in) :: n
        real(kind=4), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        real(kind=4), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        end subroutine mkl_vm_vstrunci_omp_offload_ilp64

        subroutine mkl_vm_vmstrunci_omp_offload_ilp64(n, a, inca, y, incy, mode)
        integer(kind=8), intent(in) :: n
        real(kind=4), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        real(kind=4), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmstrunci_omp_offload_ilp64

        subroutine mkl_vm_vdtrunci_omp_offload_ilp64(n, a, inca, y, incy)
        integer(kind=8), intent(in) :: n
        real(kind=8), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        real(kind=8), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        end subroutine mkl_vm_vdtrunci_omp_offload_ilp64

        subroutine mkl_vm_vmdtrunci_omp_offload_ilp64(n, a, inca, y, incy, mode)
        integer(kind=8), intent(in) :: n
        real(kind=8), intent(in) :: a(*)
        integer(kind=8), intent(in) :: inca
        real(kind=8), intent(inout) :: y(*)
        integer(kind=8), intent(in) :: incy
        integer(kind=8), intent(in) :: mode
        end subroutine mkl_vm_vmdtrunci_omp_offload_ilp64


    end interface
