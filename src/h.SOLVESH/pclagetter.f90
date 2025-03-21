        subroutine pclagetter( a, desca, alpha, norbitals)
        implicit none
!
!     written by kurt r. glaesemann
!
!     extension of, but not part of
!  -- scalapack tools routine (version 1.5) --
!
        double precision, intent (in), dimension (*) :: a
        integer, intent (in), dimension (*) :: desca
        integer, intent (in) :: norbitals
        double precision, intent (out), dimension (norbitals,norbitals) :: alpha

!     .. things that used to be scalar arguments but we fix ..
        integer :: ia, icprnt, irprnt, ja, m, n
!     ..
!
!  purpose
!  =======
!
!  this stores values of a in local array alpha (*) on
!  the the process of coordinates (irprnt, icprnt)
!
!  notes
!  =====
!
!  each global data object is described by an associated description
!  vector.  this vector stores the information required to establish
!  the mapping between an object element and its corresponding process
!  and memory location.
!
!  arguments
!  =========
!
!  m       (global input) integer
!          the number of rows to be operated on i.e the number of rows
!          of the distributed submatrix sub( a ). m >= 0.
!
!  n       (global input) integer
!          the number of columns to be operated on i.e the number of
!          columns of the distributed submatrix sub( a ). n >= 0.
!
!  a       (local input) complex pointer into the local memory to a
!          local array of dimension (lld_a, locc(ja+n-1) ) containing
!          the local pieces of the distributed matrix sub( a ).
!
!  ia      (global input) integer
!          the row index in the global array a indicating the first
!          row of sub( a ).
!
!  ja      (global input) integer
!          the column index in the global array a indicating the
!          first column of sub( a ).
!
!  desca   (global and local input) integer array
!          the array descriptor for the distributed matrix a.
!
!  irprnt  (global input) integer
!          the row index of the printing process.
!
!  icprnt  (global input) integer
!          the column index of the printing process.
!
!  norbitals (local input) integer size of alpha matrix
!  =====================================================================
!
!     .. parameters ..
        integer, parameter ::   ctxt_=2, mb_=5, nb_=6, lld_=9
!     .. local scalars ..
        integer            kk, i, iacol, iarow, ib, ictxt, icurcol,  &
        &                   icurrow, ii, iia, in, j, jb, jj, jja, jn, k,  &
        &                   lda, nb_a, mb_a, mycol, myrow, npcol, nprow
!     ..
!     .. external functions ..
        integer            iceil
        external           iceil
!     ..
!     .. executable statements ..
!
!     these are usually inputs, but we will require square matrices
!     and that a begins at (0,0), and that the master get the data
!     this simplyfies the subroutine CALL.
!
        ia=1
        ja=1
        irprnt=0
        icprnt=0
        n=norbitals
        m=norbitals
!
!     get grid parameters
!
        ictxt = desca( ctxt_ )
        call blacs_gridinfo( ictxt, nprow, npcol, myrow, mycol )
!
        call infog2l( ia, ja, desca, nprow, npcol, myrow, mycol,  &
        &              iia, jja, iarow, iacol )
        icurrow = iarow
        icurcol = iacol
        ii = iia
        jj = jja
        lda = desca( lld_ )
        nb_a = desca( nb_ )
        mb_a = desca( mb_ )
!
!     handle the first block of column separately
!
        jn = min( iceil( ja, nb_a ) * nb_a, ja+n-1 )
        in = min( iceil( ia, mb_a ) * mb_a, ia+m-1 )
        jb = jn-ja+1
        ib = in-ia+1
        j = ja
        i = ia
        if( icurrow.eq.irprnt .and. icurcol.eq.icprnt ) then
          if( myrow.eq.irprnt .and. mycol.eq.icprnt ) then
            do kk = 0, jb-1
              do k = 0, ib-1
                alpha (i+k, j+kk)=a(ii+k+(jj+kk-1)*lda)
              end do
            end do
          end if
        else
          if( myrow.eq.icurrow .and. mycol.eq.icurcol ) then
            call dgesd2d( ictxt, ib, jb, a( ii+(jj-1)*lda ), lda,  &
            &                    irprnt, icprnt )
          else if( myrow.eq.irprnt .and. mycol.eq.icprnt ) then
            call dgerv2d( ictxt, ib, jb, alpha (i, j), norbitals,  &
            &                    icurrow, icurcol )
          end if
        end if
        if( myrow.eq.icurrow ) ii = ii + ib
        icurrow = mod( icurrow+1, nprow )
        call blacs_barrier( ictxt, 'all' )
!
!     loop over remaining block of rows (in first column)
!
        do 50 i = in+1, ia+m-1, mb_a
          ib = min( mb_a, ia+m-i )
          if( icurrow.eq.irprnt .and. icurcol.eq.icprnt ) then
            if( myrow.eq.irprnt .and. mycol.eq.icprnt ) then
              do kk = 0, jb-1
                do k = 0, ib-1
                  alpha(i+k, j+kk)=a( ii+k+(jj+kk-1)*lda )
                end do
              end do
            end if
          else
            if( myrow.eq.icurrow .and. mycol.eq.icurcol ) then
              call dgesd2d( ictxt, ib, jb, a( ii+(jj-1)*lda ),  &
              &                       lda, irprnt, icprnt )
            else if( myrow.eq.irprnt .and. mycol.eq.icprnt ) then
              call dgerv2d( ictxt, ib, jb, alpha(i, j), norbitals,  &
              &                       icurrow, icurcol )
            end if
          end if
          if( myrow.eq.icurrow ) ii = ii + ib
          icurrow = mod( icurrow+1, nprow )
          call blacs_barrier( ictxt, 'all' )
50      continue
!
        ii = iia
        icurrow = iarow
!
        if( mycol.eq.icurcol ) jj = jj + jb
        icurcol = mod( icurcol+1, npcol )
        call blacs_barrier( ictxt, 'all' )
!
!     loop over remaining column blocks
!
        do 130 j = jn+1, ja+n-1, nb_a
          jb = min(  nb_a, ja+n-j )
          in = min( iceil( ia, mb_a ) * mb_a, ia+m-1 )
          ib = in-ia+1
          i = ia
          if( icurrow.eq.irprnt .and. icurcol.eq.icprnt ) then
            if( myrow.eq.irprnt .and. mycol.eq.icprnt ) then
              do kk = 0, jb-1
                do k = 0, ib-1
                  alpha(i+k, j+kk)=a( ii+k+(jj+kk-1)*lda )
                end do
              end do
            end if
          else
            if( myrow.eq.icurrow .and. mycol.eq.icurcol ) then
              call dgesd2d( ictxt, ib, jb, a( ii+(jj-1)*lda ),  &
              &                       lda, irprnt, icprnt )
            else if( myrow.eq.irprnt .and. mycol.eq.icprnt ) then
              call dgerv2d( ictxt, ib, jb, alpha(i, j), norbitals,  &
              &                       icurrow, icurcol )
            end if
          end if
          if( myrow.eq.icurrow ) ii = ii + ib
          icurrow = mod( icurrow+1, nprow )
          call blacs_barrier( ictxt, 'all' )
!
!        loop over remaining block of rows
!
          do 110 i = in+1, ia+m-1, mb_a
            ib = min( mb_a, ia+m-i )
            if( icurrow.eq.irprnt .and. icurcol.eq.icprnt ) then
              if( myrow.eq.irprnt .and. mycol.eq.icprnt ) then
                do kk = 0, jb-1
                  do k = 0, ib-1
                    alpha(i+k, j+kk)=a( ii+k+(jj+kk-1)*lda )
                  end do
                end do
              end if
            else
              if( myrow.eq.icurrow .and. mycol.eq.icurcol ) then
                call dgesd2d( ictxt, ib, jb, a( ii+(jj-1)*lda ),  &
                &                          lda, irprnt, icprnt )
              else if( myrow.eq.irprnt .and. mycol.eq.icprnt ) then
                call dgerv2d( ictxt, ib, jb, alpha(i, j),  &
                &                          norbitals, icurrow, icurcol )
              end if
            end if
            if( myrow.eq.icurrow ) ii = ii + ib
            icurrow = mod( icurrow+1, nprow )
            call blacs_barrier( ictxt, 'all' )
110       continue
!
          ii = iia
          icurrow = iarow
!
          if( mycol.eq.icurcol ) jj = jj + jb
          icurcol = mod( icurcol+1, npcol )
          call blacs_barrier( ictxt, 'all' )
!
130     continue
!
        return
        end subroutine pclagetter
