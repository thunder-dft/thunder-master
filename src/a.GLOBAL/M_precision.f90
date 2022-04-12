!! NAME
!! M_precision
!!
!! FUNCTION
!! This module contains definitions for a number of constants:
!! - Mathematical constants
!! - Computational constants
!! - Physical constants
!!
!!
!! NOTES
!! The content of this file is derived from 'Numerical Recipes in Fortran 90'
!! W.H. Press et al., volume 2 of 'Fortran Numerical Recipes', Cambridge
!! University Press, Second Edition (1996), p. 937 and 1361
!!
!! SOURCE

        module M_precision
        implicit none

!Keyword 'integer' stands for default integer type
!and may be used whenever integer are presumed to be small

!nb of bytes related to an integer subtype n such as -10^(argument) < n < 10^(argument) (this is standard F90)
        integer, parameter :: i1b=selected_int_kind(2)
        integer, parameter :: i2b=selected_int_kind(4)
        integer, parameter :: i4b=selected_int_kind(9)
        integer, parameter :: i8b=selected_int_kind(18)

!nb of bytes related to default simple-precision real/complex subtypes
!(= 4 for many machine architectures, = 8 for e.g. Cray)
        integer, parameter :: sp=kind(1.0)          ! Single precision should not be used
        integer, parameter :: spc=kind((1.0,1.0))

!nb of bytes related to default double-precision real/complex subtypes
!(= 8 for many machine architectures)
        integer, parameter :: dp=kind(1.0d0)
        integer, parameter :: dpc=kind((1.0_dp,1.0_dp))  ! Complex should not be used presently
                                                         ! except for use of libraries

!Identifiers for important files
        integer, parameter :: ilogfile=21                ! output.log unit number

!Mathematical constants
        double precision, parameter :: pi = 3.141592653589793238462643d0  ! pi = 4.0d0*atan(1.0d0)
        double precision, parameter :: pisq3 = 29.6088132032680740d0
        double precision, parameter :: epsilon = 1.0d-15

        end module M_precision
