!=======================================================================
!> @brief Special purpose routine for the output of selected parameters
!> @details All user-specified expressions must be written in standard
!! Fortran 90 or in the Fortran version available on your machine. <br>
!! To activate the statements, please remove the 'c' from the
!! first column. <br>
!! This routine will be executed at the start of calculation. <br>
!! See the list of variables for the use of this routine.
      SUBROUTINE useout(idum)
!=======================================================================
!.....contact fire@avl.com
!-----
!-----------------------------------------------------------------------
!-----
      USE comm0
      USE comm1
      USE prec_mod, ONLY: PREC
      USE fire_mpi_var

      IMPLICIT NONE
!-----PARAMETER
      INTEGER, INTENT(IN) :: IDUM               !< not used parameter
      integer :: i,num,j , isel, ityp, nelem, k, nc, mat
      real*8 :: ri
      real*8 :: temp,wgt, pos_min
      logical :: lop
      real*8, allocatable, dimension(:) :: r_ave
      real*8 :: avesot12, ms, actsr
      CHARACTER(len=128) :: sname
      character(len=16),dimension(:),allocatable :: cname
!      data cname /'sen','1'/
!      data cname /'vol'/
!ms      include 'comdp.inc'
!ms      include 'com0.inc'
      include 'SwiftIO_FortranFunctions.inc'
      
      double precision, dimension (3) :: v1
      REAL(PREC) :: pos_p
      REAL(PREC) :: pos_y
      INTEGER :: i_fnd_p
      REAL(prec), ALLOCATABLE, DIMENSION(:) :: dconv      !< cell array to be plot
      CHARACTER(len=256) :: txt1,txt2
      integer ::  ip1, ip2

!-----
!-----------------------------------------------------------------------
!-----
      IF(I_USEOUT == 1) THEN
!-----
!-----
      !        Check  Continuity equation

        ALLOCATE (dconv(ncell)); dconv = zero
        mat = 1

!----- Initialise mass fluxes
!         u(1,1:ncell) = 3.0d0 !2.0d+1
        ! CALL initfl(mat)
        do j = 1, nface
          ip1 = lf(1,j)
          ip2 = lf(2,j)
          dconv(ip1) = dconv(ip1) + F(j)
          dconv(ip2) = dconv(ip2) - F(j)
        enddo
        do j= 1, nbfac
          ip1 = lb(j)
          dconv(ip1) = dconv(ip1) + FB(j)
        enddo

        write(txt1,'(a,1es13.5)')" R.M.S. convection ",SUM(dconv(nsp(mat):nep(mat)))/(nep(mat)-nsp(mat)+1)
        txt2=''
        ! IF(iampro<2) CALL UCMESS('usedef','F',txt1,txt2)
      
        CALL FIO_UUMESS ('USEOUT','I',txt1,txt2)

        DEALLOCATE (dconv)

      ELSE IF(I_USEOUT == 2) THEN
!-----
      END IF
!-----
      RETURN
      END SUBROUTINE useout
