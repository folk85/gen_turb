c=======================================================================
!> @brief Special purpose routine for generating user defined output
!> @details All user-specified expressions must be written in standard
!! Fortran 90 or in the Fortran version available on your machine. <br>
!! To activate the statements, please remove the 'c' from the
!! first column. <br>
!! This routine will be executed at the end of calculation. <br>
!! See the list of variables for the use of this routine.
!! <br>
!! New way to plot user-defined data into the fl3-file:
!! <ol>
!! <li>allocate 2 dynamic arrays: one with dimension ncell
!!     the other with dimension nbfac</li>
!! <li>store your user-defined quantity in these two arrays</li>
!! <li>call the routine which will perform the work for you:
!!      WRITE_USER_FL3 (mat,mph,ifile,Sname,Sunit,ncell,nbfac,
!!                      floatcell,floatbd)
!!     - the first 3 arguments are the same as USEPLO
!!     - Sname is a string containing the name of your data,
!!          as it will be displayed in IMPRESS
!!     - Sunit is a string containing the unit of your data,
!!          displayed in IMPRESS as well
!!     - ncell and nbfac are the dimensions
!!     - floatcell and floatbd are the local arrays
!! <li>perform operations 1 to 3 for all your user data</li>
!! <li>at the end deallocate the dynamic arrays</li>
!! </ol>
!======================================================================
      SUBROUTINE useplo(mat,mph1,ifile)
!=======================================================================
!
!     USEPLO IS A SPECIAL PURPOSE ROUTINE FOR GENERATING USER
!     DEFINED OUTPUT
!     ALL USER-SPECIFIED EXPRESSIONS MUST BE WRITTEN IN STANDARD
!     FORTRAN 90 OR IN THE FORTRAN VERSION AVAILABLE ON YOUR MACHINE
!
!     THIS ROUTINE WILL BE EXECUTED AT THE END OF CALCULATION
!
!     See the list of variables for the use of this routine
!
!.....contact fire@avl.com
!
!-----------------------------------------------------------------------
!     TO ACTIVATE THE EXAMPLE CODE SET VARIABLE "IS_EXAMPLE_CODE" TO
!     "FALSE".
!
!     Note: "I_USEPLO" is the integer value specified in the GUI at
!           "User-functions / Activation".
!
!-----------------------------------------------------------------------
!
      USE comm0
      USE comm1
      USE comm2
      USE prec_mod, ONLY: prec     !...precision in FIRE
      USE cfd_mpi_mod
!
      IMPLICIT NONE
!
      LOGICAL, PARAMETER :: IS_EXAMPLE_CODE = .TRUE.
!
!     arguments
      INTEGER, INTENT(IN) :: mat   !< material index
      INTEGER, INTENT(IN) :: mph1  !< multiphase flag
      INTEGER, INTENT(IN) :: ifile !< file number
!
!     local variables
      CHARACTER(len=1)   :: c_mph                         !< string for multiphase index
      INTEGER            :: ib                            !< boundary face index
      INTEGER            :: ib1                           !< boundary face index for multiphase array
      INTEGER            :: ip                            !< cell index
      INTEGER            :: ip1                           !< cell index for multiphase array
      INTEGER            :: ir                            !< boundary region index
      INTEGER            :: mph                           !< phase index
      REAL(prec), ALLOCATABLE, DIMENSION(:) :: floatcell  !< cell array to be plot
      REAL(prec), ALLOCATABLE, DIMENSION(:) :: floatbd    !< boundary array to be plot
      REAL(prec), ALLOCATABLE, DIMENSION(:) :: dconv      !< cell array to be plot
      integer             :: i,j, ip2  !
!
!-----------------------------------------------------------------------
      IF(I_USEPLO == 1) THEN
!
!       allocate plot arrays with dimension "number of cells" and "number of boundary faces"
        ALLOCATE (floatcell(ncell)); floatcell = zero
        ALLOCATE (dconv(ncell)); dconv = zero
        ALLOCATE (floatbd(nbfac));  floatbd = zero
      !        Check  Continuity equation

!----- Initialise mass fluxes
!         u(1,1:ncell) = 3.0d0 !2.0d+1
        ! CALL initfl(mat)
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

!       loop over all cells in the current domain (material)
        DO ip = nsp(mat), nep(mat)
!
!         store specific heat capacity on plot array
          floatcell(ip) = dconv(ip)
!
        END DO
!
!       loop over all boundary regions
        DO ir = 0, nreg
!
!         check for boundary region in current domain (material)
          IF (ibc(2,ir) /= mat) CYCLE
!
!         loop over boundary faces of the current region
          DO ib = nsr(ir), ner(ir)
            floatbd(ib) = dconv(lb(ib))
          END DO
!
        END DO
!
!       write 3D result file
        CALL Write_User_Fl3 (mat,mph1,ifile,'Conversation',
     &       'kg/s',ncell,nbfac,floatcell,floatbd)
!
!       deallocate plot arrays
        DEALLOCATE (dconv)

        floatcell(:) = zero
        floatbd(:) = zero
!
!       loop over all cells in the current domain (material)
        DO ip = nsp(mat), nep(mat)
!
!         store specific heat capacity on plot array
          floatcell(ip) = xp(2,ip)
!
        END DO
!
!       loop over all boundary regions
        DO ir = 0, nreg
!
!         check for boundary region in current domain (material)
          IF (ibc(2,ir) /= mat) CYCLE
!
!         loop over boundary faces of the current region
          DO ib = nsr(ir), ner(ir)
            floatbd(ib) = xb(2,ib)
          END DO
!
        END DO
!
!       write 3D result file
        CALL Write_User_Fl3 (mat,mph1,ifile,'Position.Y',
     &       'm',ncell,nbfac,floatcell,floatbd)
!
!       deallocate plot arrays
        DEALLOCATE (floatcell)
        DEALLOCATE (floatbd)

!
!-----------------------------------------------------------------------
!     Example 1: Plot single phase 3D results of specific heat capacity
!-----------------------------------------------------------------------
      ELSE IF(I_USEPLO == 2) THEN
!
       IF (.NOT. IS_EXAMPLE_CODE) THEN
!
!       allocate plot arrays with dimension "number of cells" and "number of boundary faces"
        ALLOCATE (floatcell(ncell)); floatcell = zero
        ALLOCATE (floatbd(nbfac));  floatbd = zero
!
!       loop over all cells in the current domain (material)
        DO ip = nsp(mat), nep(mat)
!
!         store specific heat capacity on plot array
          floatcell(ip) = cpcof(ip)
!
        END DO
!
!       loop over all boundary regions
        DO ir = 0, nreg
!
!         check for boundary region in current domain (material)
          IF (ibc(2,ir) /= mat) CYCLE
!
!         loop over boundary faces of the current region
          DO ib = nsr(ir), ner(ir)
            floatbd(ib) = cpcofb(ib)
          END DO
!
        END DO
!
!       write 3D result file
        CALL Write_User_Fl3 (mat,mph1,ifile,'SpecificHeat',
     &       'J/kgK',ncell,nbfac,floatcell,floatbd)
!
!       deallocate plot arrays
        DEALLOCATE (floatcell)
        DEALLOCATE (floatbd)
!
       END IF
!
!--------------------------------------------------------------------------------
!     Example 2: Plot temperature in mulitphase in degree Celsius for all phases
!--------------------------------------------------------------------------------
      ELSE IF(I_USEPLO == 3) THEN
!
       IF (.NOT. IS_EXAMPLE_CODE) THEN
!
!       allocate plot arrays with dimension "number of cells" and "number of boundary faces"
        ALLOCATE (floatcell(ncell)); floatcell = zero
        ALLOCATE (floatbd(nbfac));  floatbd = zero
!
!       check for multiphase
        IF (lphase) then
!
!        loop over all phases
         DO mph = 1, nph
!
!         convert phase index into character
          WRITE(c_mph,'(I1)') mph
!
!         loop over all cells
          DO ip = nsp(mat), nep(mat)
!
!           index for multiphase variables
            ip1 = (mph - 1) * ncell + ip
!
!           store cell temperature on array
            floatcell(ip) = tm(ip1)-273.15
!
          END DO
!
!         loop over all boundary regions
          DO ir = 0, nreg
!
!          check for boundary region in current domain (material)
           IF (ibc(2,ir) /= mat) CYCLE
!
!          loop over boundary faces of the current region
           DO ib = nsr(ir), ner(ir)
!
!            index for multiphase variables
             ib1 = (mph - 1) * nbfac + ib
!
!            store boundary temperature of phase on array
             floatbd(ib) = tmb(ib1)-273.15
!
           END DO
!
          END DO
!
!         write 3D results for every phase
          CALL Write_User_Fl3 (mat,mph1,ifile,'TempC'//'_Ph'//c_mph,
     &         'degC',ncell,nbfac,floatcell,floatbd)
!
!        END phase loop
         END DO
!
        END IF
!
!       deallocate plot arrays
        DEALLOCATE (floatcell)
        DEALLOCATE (floatbd)
!
       END IF
!
      END IF
!
      RETURN
      END SUBROUTINE useplo
!

      SUBROUTINE CVS_USEPLO_F()
        CHARACTER(LEN=150) :: CVS_USEPLO_F_ID
        DATA CVS_USEPLO_F_ID
     1/'@(#)
     1  $Revision: 1.16.uf $
     1  $ASTFile: 14eb5b33fb04e6d1 $
     1  '/
        WRITE(*,*) CVS_USEPLO_F_ID
      END SUBROUTINE CVS_USEPLO_F
