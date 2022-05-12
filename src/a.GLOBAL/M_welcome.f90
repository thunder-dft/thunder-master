! copyright info:
!
!                             @Copyright 2022
!                           Fireball Committee
! Hong Kong Quantum AI Laboratory, Ltd. - James P. Lewis, Chair
! Universidad de Madrid - Jose Ortega
! Academy of Sciences of the Czech Republic - Pavel Jelinek
! Arizona State University - Otto F. Sankey

! Previous and/or current contributors:
! Auburn University - Jian Jun Dong
! California Institute of Technology - Brandon Keith
! Czech Institute of Physics - Prokop Hapala
! Czech Institute of Physics - Vladimír Zobač
! Dublin Institute of Technology - Barry Haycock
! Pacific Northwest National Laboratory - Kurt Glaesemann
! University of Texas at Austin - Alex Demkov
! Ohio University - Dave Drabold
! Synfuels China Technology Co., Ltd. - Pengju Ren
! Washington University - Pete Fedders
! West Virginia University - Ning Ma and Hao Wang
! also Gary Adams, Juergen Frisch, John Tomfohr, Kevin Schmidt,
!      and Spencer Shellman

!
! RESTRICTED RIGHTS LEGEND
! Use, duplication, or disclosure of this software and its documentation
! by the Government is subject to restrictions as set forth in subdivision
! { (b) (3) (ii) } of the Rights in Technical Data and Computer Software
! clause at 52.227-7013.

! M_welcome
! Program Description
! ===========================================================================
!      This is a module containing subroutines that make the various
! welcome banners.
! ===========================================================================
! Code written by:
! James P. Lewis
! Unit 909 of Buidling 17W
! 17 Science Park West Avenue
! Pak Shek Kok, New Territories 999077
! Hong Kong
!
! Phone: +852 6612 9539 (mobile)
! ===========================================================================
! Module Declaration
! ===========================================================================
        module M_welcome
        use M_precision

! Type Declaration
! ===========================================================================
! None

! module procedures
        contains

! ===========================================================================
! banner.f90
! ===========================================================================
! Program Description
! ===========================================================================
!       Writes the banner used in all routines displaying contributors.
!
! ===========================================================================
! Code written by:
!> @author James P. Lewis
! Unit 909 of Buidling 17W
! 17 Science Park West Avenue
! Pak Shek Kok, New Territories 999077
! Hong Kong
!
! Phone: +852 6612 9539 (mobile)
! ===========================================================================
!
! Program Declaration
! ===========================================================================
        subroutine banner

! /GLOBAL
        use M_precision

        implicit none

! Argument Declaration and Description
! ===========================================================================
! None

! Local Parameters and Data Declaration
! ===========================================================================
! None

! Local Variable Declaration and Description
! ===========================================================================
! None

! Procedure
! ===========================================================================
        write (ilogfile,'(4x,A)')
        write (ilogfile,'(4x,A)') '                   James P. Lewis '
        write (ilogfile,'(4x,A)') '                Hong Kong University '
        write (ilogfile,'(4x,A)') '        Hong Kong Quantum AI Laboratory, Ltd. '
        write (ilogfile,'(4x,A)')
        write (ilogfile,'(4x,A)') '                    Jose Ortega '
        write (ilogfile,'(4x,A)') '               Departmento de Fisica '
        write (ilogfile,'(4x,A)') '          Teorica de la Materia Condensada '
        write (ilogfile,'(4x,A)') '               Universidad de Madrid '
        write (ilogfile,'(4x,A)')
        write (ilogfile,'(4x,A)') '                   Pavel Jelinek '
        write (ilogfile,'(4x,A)') '                Institute of Physics '
        write (ilogfile,'(4x,A)') '               Prague, Czech Republic '
        write (ilogfile,'(4x,A)')
        write (ilogfile,'(4x,A)') '                     Pengju Ren '
        write (ilogfile,'(4x,A)') '          Synfuels China Technology Co., Ltd. '
        write (ilogfile,'(4x,A)') '                 Beijing, P.R. China '
        write (ilogfile,'(4x,A)')
        write (ilogfile,'(4x,A)') '                   Otto F. Sankey '
        write (ilogfile,'(4x,A)') '         Department of Physics and Astronomy '
        write (ilogfile,'(4x,A)') '              Arizona State University '
        write (ilogfile,*)

        write (ilogfile,*)
        write (ilogfile,'(4x,A)') ' Contributions from: '
        write (ilogfile,'(4x,A)') ' - Alex A. Demkov (University of Texas - Austin) '
        write (ilogfile,'(4x,A)') ' - Jian Jun Dong (Auburn University) '
        write (ilogfile,'(4x,A)') ' - David A. Drabold (Ohio University) '
        write (ilogfile,'(4x,A)') ' - Peter A. Fedders (Washington University) '
        write (ilogfile,'(4x,A)') ' - Kurt R. Glaesemann (Pacific Northwest National Laboratory) '
        write (ilogfile,'(4x,A)') ' - Prokop Hapala (Czech Institute of Physics) '
        write (ilogfile,'(4x,A)') ' - Barry Haycock (Dublin Institute of Technology) '
        write (ilogfile,'(4x,A)') ' - Brandon Keith (Caltech) '
        write (ilogfile,'(4x,A)') ' - Hao Wang (West Virginia University) '
        write (ilogfile,'(4x,A)') ' - Ning Ma (West Virginia University) '
        write (ilogfile,'(4x,A)') ' - Vladimír Zobač (Czech Institute of Physics) '
        write (ilogfile,'(4x,A)')
        write (ilogfile,'(4x,A)') ' also Gary Adams, John Tomfohr, Juergen Frisch, '
        write (ilogfile,'(4x,A)') '      Kevin Schmidt, and Spencer Shellman '

        write (ilogfile,*)
        write (ilogfile,'(4x,A)') ' Copyright Information: '
        write (ilogfile,'(4x,A)')
        write (ilogfile,'(4x,A)') ' Usable only with permission from the FIREBALL executive '
        write (ilogfile,'(4x,A)') ' committee. This program is NOT, under any circumstances, '
        write (ilogfile,'(4x,A)') ' to be transfered to an unauthorized user. '
        write (ilogfile,*)
        write (ilogfile,'(4x,A)') ' ==================================================== '
        write (ilogfile,*)

! End Subroutine
! ===========================================================================
        return
        end subroutine banner

! ===========================================================================
! welcome_fireball.f90
! ===========================================================================
! Program Description
! ===========================================================================
!       Writes the banner for FIREBALL.
!
! ===========================================================================
! Code written by:
!> @author James P. Lewis
! Unit 909 of Buidling 17W
! 17 Science Park West Avenue
! Pak Shek Kok, New Territories 999077
! Hong Kong
!
! Phone: +852 6612 9539 (mobile)
! ===========================================================================
!
! Program Declaration
! ===========================================================================
        subroutine welcome_fireball

        use M_precision
        implicit none

! Argument Declaration and Description
! ===========================================================================
! None

! Local Parameters and Data Declaration
! ===========================================================================
! None

! Local Variable Declaration and Description
! ===========================================================================
! None

! Procedure
! ===========================================================================
        write (ilogfile,*)
        write (ilogfile,'(4x,A)') '       ______  _             _             _  _  '
        write (ilogfile,'(4x,A)') '      |  ____|(_)           | |           | || | '
        write (ilogfile,'(4x,A)') '      | |__    _  _ __  ___ | |__    __ _ | || | '
        write (ilogfile,'(4x,A)') '      |  __|  | || `__|/ _ \| `_ \  / _` || || | '
        write (ilogfile,'(4x,A)') '      | |     | || |  |  __/| |_) || (_| || || | '
        write (ilogfile,'(4x,A)') '      |_|     |_||_|   \___||_.__/  \__,_||_||_| '
        write (ilogfile,*)
        write (ilogfile,'(4x,A)') ' ==================================================== '
        write (ilogfile,*)
        write (ilogfile,'(4x,A)') '                    Version 2022 '
        write (ilogfile,'(4x,A)') '          A fast local orbital QMD Package '

! End Subroutine
! ===========================================================================
        return
        end subroutine welcome_fireball

! ===========================================================================
! welcome_begin.f90
! ===========================================================================
! Program Description
! ===========================================================================
!       Writes the banner for BEGIN.
!
! ===========================================================================
! Code written by:
!> @author James P. Lewis
! Unit 909 of Buidling 17W
! 17 Science Park West Avenue
! Pak Shek Kok, New Territories 999077
! Hong Kong
!
! Phone: +852 6612 9539 (mobile)
! ===========================================================================
!
! Program Declaration
! ===========================================================================
        subroutine welcome_begin

        use M_precision
        implicit none

! Argument Declaration and Description
! ===========================================================================
! None

! Local Parameters and Data Declaration
! ===========================================================================
! None

! Local Variable Declaration and Description
! ===========================================================================
! None

! Procedure
! ===========================================================================
        write (ilogfile,*)
        write (ilogfile,'(4x,A)') '            ____             _       '
        write (ilogfile,'(4x,A)') '           |  _ \           (_)      '
        write (ilogfile,'(4x,A)') '           | |_) | ___  __ _ _ _ __  '
        write (ilogfile,'(4x,A)') '           |  _ < / _ \/ _` | | `_ \ '
        write (ilogfile,'(4x,A)') '           | |_) |  __/ (_| | | | | |'
        write (ilogfile,'(4x,A)') '           |____/ \___|\__, |_|_| |_|'
        write (ilogfile,'(4x,A)') '                        __/ |        '
        write (ilogfile,'(4x,A)') '                       |___/         '

        write (ilogfile,*)
        write (ilogfile,'(4x,A)') ' ==================================================== '
        write (ilogfile,*)
        write (ilogfile,'(4x,A)') '                    Version 2022 '
        write (ilogfile,'(4x,A)') '      Sankey-Niklewski wave-functions generator '

! End Subroutine
! ===========================================================================
        return
        end subroutine welcome_begin

! ===========================================================================
! welcome_create.f90
! ===========================================================================
! Program Description
! ===========================================================================
!       Writes the banner for CREATE.
!
! ===========================================================================
! Code written by:
!> @author James P. Lewis
! Unit 909 of Buidling 17W
! 17 Science Park West Avenue
! Pak Shek Kok, New Territories 999077
! Hong Kong
!
! Phone: +852 6612 9539 (mobile)
! ===========================================================================
!
! Program Declaration
! ===========================================================================
        subroutine welcome_create

        use M_precision
        implicit none

! Argument Declaration and Description
! ===========================================================================
! None

! Local Parameters and Data Declaration
! ===========================================================================
! None

! Local Variable Declaration and Description
! ===========================================================================
! None

! Procedure
! ===========================================================================
        write (ilogfile,*)
        write (ilogfile,'(4x,A)') '             _____                _        '
        write (ilogfile,'(4x,A)') '            / ____|              | |       '
        write (ilogfile,'(4x,A)') '           | |     _ __ ___  __ _| |_ ___  '
        write (ilogfile,'(4x,A)') '           | |    | `__/ _ \/ _` | __/ _ \ '
        write (ilogfile,'(4x,A)') '           | |____| | |  __/ (_| | ||  __/ '
        write (ilogfile,'(4x,A)') '            \_____|_|  \___|\__,_|\__\___| '
        write (ilogfile,*)
        write (ilogfile,'(4x,A)') ' ==================================================== '
        write (ilogfile,*)
        write (ilogfile,'(4x,A)') '                    Version 2022 '
        write (ilogfile,'(4x,A)') '      Sankey-Niklewski wave-functions generator '

! End Subroutine
! ===========================================================================
        return
        end subroutine welcome_create

! ===========================================================================
! welcome_lightning.f90
! ===========================================================================
! Program Description
! ===========================================================================
!       Writes the banner for LIGHTNING.
!
! ===========================================================================
! Code written by:
!> @author James P. Lewis
! Unit 909 of Buidling 17W
! 17 Science Park West Avenue
! Pak Shek Kok, New Territories 999077
! Hong Kong
!
! Phone: +852 6612 9539 (mobile)
! ===========================================================================
!
! Program Declaration
! ===========================================================================
        subroutine welcome_lightning

        use M_precision
        implicit none

! Argument Declaration and Description
! ===========================================================================
! None

! Local Parameters and Data Declaration
! ===========================================================================
! None

! Local Variable Declaration and Description
! ===========================================================================
! None

! Procedure
! ===========================================================================
        write (ilogfile,*)
        write (ilogfile,'(4x,A)') '     _      _       _     _         _             '
        write (ilogfile,'(4x,A)') '    | |    (_)     | |   | |       (_)            '
        write (ilogfile,'(4x,A)') '    | |     _  __ _| |__ | |_ _ __  _ _ __   __ _ '
        write (ilogfile,'(4x,A)') '    | |    | |/ _` | `_ \| __| `_ \| | `_ \ / _` |'
        write (ilogfile,'(4x,A)') '    | |____| | (_| | | | | |_| | | | | | | | (_| |'
        write (ilogfile,'(4x,A)') '    |______|_|\__, |_| |_|\__|_| |_|_|_| |_|\__, |'
        write (ilogfile,'(4x,A)') '               __/ |                         __/ |'
        write (ilogfile,'(4x,A)') '              |___/                         |___/ '
        write (ilogfile,*)
        write (ilogfile,'(4x,A)') ' ==================================================== '
        write (ilogfile,*)
        write (ilogfile,'(4x,A)') '                    Version 2022 '
        write (ilogfile,'(4x,A)') '          A fast local orbital QMD Package '

! End Subroutine
! ===========================================================================
        return
        end subroutine welcome_lightning

! ===========================================================================
! welcome_thunder.f90
! ===========================================================================
! Program Description
! ===========================================================================
!       Writes the banner for THUNDER.
!
! ===========================================================================
! Code written by:
!> @author James P. Lewis
! Unit 909 of Buidling 17W
! 17 Science Park West Avenue
! Pak Shek Kok, New Territories 999077
! Hong Kong
!
! Phone: +852 6612 9539 (mobile)
! ===========================================================================
!
! Program Declaration
! ===========================================================================
        subroutine welcome_thunder

        use M_precision
        implicit none

! Argument Declaration and Description
! ===========================================================================
! None

! Local Parameters and Data Declaration
! ===========================================================================
! None

! Local Variable Declaration and Description
! ===========================================================================
! None

! Procedure
! ===========================================================================
        write (ilogfile,*)
        write (ilogfile,'(4x,A)') '      _______ _                     _           '
        write (ilogfile,'(4x,A)') '     |__   __| |                   | |          '
        write (ilogfile,'(4x,A)') '        | |  | |__  _   _ _ __   __| | ___ _ __ '
        write (ilogfile,'(4x,A)') '        | |  | `_ \| | | | `_ \ / _` |/ _ \ `__|'
        write (ilogfile,'(4x,A)') '        | |  | | | | |_| | | | | (_| |  __/ |   '
        write (ilogfile,'(4x,A)') '        |_|  |_| |_|\__,_|_| |_|\__,_|\___|_|   '
        write (ilogfile,*)
        write (ilogfile,'(4x,A)') ' ==================================================== '
        write (ilogfile,*)
        write (ilogfile,'(4x,A)') '                    Version 2022 '
        write (ilogfile,'(4x,A)') '          A fast local orbital QMD Package '

! End Subroutine
! ===========================================================================
        return
        end subroutine welcome_thunder


! End Module
! =============================================================================
        end module
