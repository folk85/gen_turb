!=======================================================================
!> @brief Special purpose routine for changing the run parameters
!> @details All user-specified expressions must be written in standard
!! Fortran 90 or in the Fortran version available on your machine. <br>
!! This routine will be executed after each iteration. <br>
!! See the list of variables for the use of this routine.
      SUBROUTINE useite(itout)
!=======================================================================
!.....contact fire@avl.com
!-----
!-----------------------------------------------------------------------
!-----
      USE comm0
      USE comm1
      USE comm2
      USE constants
!
      IMPLICIT NONE
      integer :: i, j, ib,jb, ip1,ip2, ir, mat

!
!-----parameters
      INTEGER, INTENT(IN) :: itout  !< current outer iteration index
!-----
!-----------------------------------------------------------------------
!-----
      IF(I_USEITE == 1) THEN

        mat = 1

!----- Initialise mass fluxes
!         u(1,1:ncell) = 3.0d0 !2.0d+1
        CALL initfl(mat)
!-----        

!         do i=1, nface
!           WRITE(*,'(a,i10,es10.2)') "USEite: Internal faces ",i, F(i)
!         enddo
        DO ir = 0,nreg
!           WRITE(*,'(a,4i)')"useite: ir:",ir, ibc(1,ir),
!      x                      ibc(2,ir),ibc(5,ir)
          IF (ibc(2,ir) .NE. mat) CYCLE
!           IF (ibc(1,ir) .EQ. 12 .AND. ibc(5,ir) .GT. 1) THEN
          IF (ibc(1,ir) .EQ. 12) THEN
            DO ib = nsr(ir),ner(ir)
              jb = lbj(ib)
              ip1 = lb(ib)
              ip2 = lb(jb)
              fb(ib)=denb(ib)*(ub(1,ib)*sb(1,ib)+ub(2,ib)*sb(2,ib)+ub(3,ib)*sb(3,ib))
!               fb(jb) = - fb(ib)
!               WRITE(*,'(a,4i,7es10.2)') "useini bond: ",ib, jb, ip1,ip2
!      x  , xp(1,ip1),xp(1,ip2), FB(ib), FB(jb),ub(1,ib),u(1,ip1),denb(ib)
!            CALL rotate(db(:,jb),dbrot(:))
!            dbrot(:) = db(:,ib)-dbrot(:)
!            wf1 = SQRT(DOT_PRODUCT(db(:,jb),db(:,jb)) / DOT_PRODUCT(dbrot,dbrot))
            ENDDO
          ENDIF
        ENDDO

!-----
      ELSE IF(I_USEITE == 2) THEN
!-----
      ELSE IF(I_USEITE == 99) THEN
         CALL USE_FORMULA('USEITE', 1, itout, 0, 0)
      END IF
!-----
      RETURN
      END SUBROUTINE useite
