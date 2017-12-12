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
      REAL(PREC) :: pos_p, dtmp
      REAL(PREC) :: pos_y
      INTEGER :: i_fnd_p
      REAL(prec), ALLOCATABLE, DIMENSION(:) :: dconv      !< cell array to be plot
      CHARACTER(len=256) :: txt1,txt2
      integer ::  ip1, ip2, ir, ib, jb
      REAL(prec), DIMENSION(1:3) :: dtmp3      !< cell array to be plot


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
        ! UPDATE INternal cells
        CALL initfl(mat)

        ! update bnd faces

        DO ir = 0,nreg
          IF (ibc(2,ir) .NE. mat) CYCLE
!           IF (ibc(1,ir) .EQ. 12 .AND. ibc(5,ir) .GT. 1) THEN
          IF (ibc(1,ir) .EQ. 12) THEN
            DO ib = nsr(ir),ner(ir)
              jb = lbj(ib)
              ip1 = lb(ib)
              ip2 = lb(jb)
              fb(ib)=denb(ib)*(ub(1,ib)*sb(1,ib)+ub(2,ib)*sb(2,ib)+ub(3,ib)*sb(3,ib))
              dtmp3(1:3) = (u(1:3,ip1) + u(1:3,ip2)) * 5.0d-1
              dtmp = 5.0d-1 * (den(ip1)+den(ip2))
              write(txt1,'(a,i5,5es13.5)')" Periodic Vels ",ip1,fb(ib),ub(1,ib),ub(1,jb),u(1,ip1),u(1,ip2)
              fb(ib)=dtmp*(dtmp3(1)*sb(1,ib)+dtmp3(2)*sb(2,ib)+dtmp3(3)*sb(3,ib))
              write(txt2,'(a,i5,4es13.5)')" Corrected Vels",ip1,fb(ib),dtmp3(1:3)
              ! txt2=''
              CALL FIO_UUMESS ('USEOUT','I',txt1,txt2)
            ENDDO
          ENDIF
        ENDDO

      !        Check  Continuity equation
        do j = 1, nface
          ip1 = lf(1,j)
          ip2 = lf(2,j)
          dconv(ip1) = dconv(ip1) - F(j)
          ! dconv(ip2) = dconv(ip2) + F(j)
        enddo
        do j= 1, nbfac
          ip1 = lb(j)
          dconv(ip1) = dconv(ip1) + FB(j)
        enddo

        DO ir = 0,nreg
          IF (ibc(2,ir) .NE. mat) CYCLE
!           IF (ibc(1,ir) .EQ. 12 .AND. ibc(5,ir) .GT. 1) THEN
          IF (ibc(1,ir) .EQ. 12) THEN
            DO ib = nsr(ir),ner(ir)
              jb = lbj(ib)
              ip1 = lb(ib)
              ip2 = lb(jb)
              write(txt1,'(a,i5,3es13.5)')" Periodic Cell  ",jb,fb(ib),dconv(ip1),dconv(ip2)
              txt2=''
              CALL FIO_UUMESS ('USEOUT','I',txt1,txt2)
            ENDDO
          ENDIF
        ENDDO

        write(txt1,'(a,1es13.5)')" R.M.S. conversion ",SUM(dconv(nsp(mat):nep(mat)))/(nep(mat)-nsp(mat)+1)
        txt2=''
        CALL FIO_UUMESS ('USEOUT','I',txt1,txt2)

        DEALLOCATE (dconv)

      ELSE IF(I_USEOUT == 2) THEN
!-----
      END IF
!-----
      RETURN
      END SUBROUTINE useout
