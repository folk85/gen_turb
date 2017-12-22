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
      REAL(PREC) :: pos_p, dtmp, dcell
      REAL(PREC) :: pos_y
      INTEGER :: i_fnd_p
      REAL(prec), ALLOCATABLE, DIMENSION(:) :: dconv      !< cell array to be plot
      CHARACTER(len=256) :: txt1,txt2
      integer ::  ip1, ip2, ir, ib, jb
      REAL(prec), DIMENSION(1:3) :: dtmp3      !< cell array to be plot
      REAL(prec), DIMENSION(1:3) :: dmean      !< mean velocity 
      REAL(prec), DIMENSION(1:3) :: drms      !< rms velocity 
      REAL(prec) :: dlocmax      !< coordinate of MAX locations of the flame
      REAL(prec) :: dlocmin      !< coordinate of MIN locations of the flame
      REAL(prec) :: dlymin, dlymax
      REAL(prec), DIMENSION(:), ALLOCATABLE :: dr1,dr2,dr3,nr      !< Correlation vector
!-----
      INTEGER :: err_out,cfd_io_unit, ifile
!-----
      !DATA cname /  /

        mat = 1

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
          dconv(ip2) = dconv(ip2) + F(j)
        enddo
        do j= 1, nbfac
          ip1 = lb(j)
          dconv(ip1) = dconv(ip1) - FB(j)
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

        write(txt1,'(a)')" Income "
        txt2=''
        CALL FIO_UUMESS ('USEOUT','I',txt1,txt2)
  
        !  Define mean values over the domain
        CALL FIO_UUMESS ('USEOUT','I','Define mean values over the domain',txt2)
        dtmp = 0.0d0
        do i=nsp(mat), nep(mat)
          if (t(i) > 1.20d+3) then
            dtmp = dtmp + 1.0d0
          endif
        enddo
        CALL FIO_UUMESS ('USEOUT','I','set calc',txt2)
        dcell = DBLE(nep(mat) - nsp(mat)+ 1)
        CALL FIO_UUMESS ('USEOUT','I','Before DGLSUM(dtmp)',txt2)
        CALL DGLSUM(dtmp)
        CALL FIO_UUMESS ('USEOUT','I','Before DGLSUM(dcell)',txt2)
        CALL DGLSUM(dcell)
        dtmp = dtmp / dcell

        CALL FIO_UUMESS ('USEOUT','I','Calc product fraction',txt2)

        ! Calc size of Y segment
        dlymax = MAXVAL(xb(2,1:nbfac))
        dlymin = MINVAL(xb(2,1:nbfac))
        CALL DGLMAX(dlymax)
        CALL DGLMIN(dlymin)
        CALL FIO_UUMESS ('USEOUT','I','Calc Max and min Y positions',txt2)

        ! r.m.s. and mean values
        do i = 1, 3
          dmean(i) = SUM(u(i,nsp(mat):nep(mat)))
        end do
        CALL DGLSUMVEC(dmean(1),3)
        dmean(1:3) = dmean(1:3) / dcell
        CALL FIO_UUMESS ('USEOUT','I','Calc Mean velocities values',txt2)

        do i = 1, 3
          drms(i) = SUM((dmean(i) - u(i,nsp(mat):nep(mat)))**2)
        end do
        CALL DGLSUMVEC(drms(1),3)
        drms(1:3) = SQRT(drms(1:3) / dcell)
        CALL FIO_UUMESS ('USEOUT','I','Calc SDT of velocities',txt2)


        dlocmax = MAXVAL(xp(2,nsp(1):nep(1)),MASK=t(nsp(1):nep(1)).GT.1.20d+3)
        dlocmin = MINVAL(xp(2,nsp(1):nep(1)),MASK=t(nsp(1):nep(1)).LE.1.20d+3)
        CALL DGLMAX(dlocmax)
        CALL DGLMIN(dlocmin)
        CALL FIO_UUMESS ('USEOUT','I','Calc MAX/MIN flame locations',txt2)


        if (mpi_master) then
          ifile= cfd_io_unit(err_out)
          CALL FIO_UUMESS ('USEOUT','I','cfd_io_unit',txt2)
          write(txt1,'(a,i5)')" Work with file num :",ifile
          txt2=''
          CALL FIO_UUMESS ('USEOUT','I',txt1,txt2)
          ifile  = 212
          INQUIRE(ifile,OPENED= lop)
          IF(.NOT. lop) THEN
            WRITE(sname,'(2a)') 'user_io_stat','.dat'
            OPEN(ifile,FILE=TRIM(sname))
            WRITE(ifile,'(2a)') 'time Product_Fraction Mean_Flame_position Min_flame_loc ', &
            'Max_flame_loc mean.X mean.Y mean.Z std.X std.Y std.Z'
          END IF ! itst.eq.1
          write(ifile,'(11es13.5)')time, dtmp, dtmp*(dlymax - dlymin)+dlymin,dlocmin, dlocmax,dmean(1:3), drms(1:3)
        endif
!-----
      ELSE IF(I_USEOUT == 99) THEN
         CALL USE_FORMULA('USEOUT', 1, idum, 0, 0)
      END IF
!-----
      RETURN
      END SUBROUTINE useout
