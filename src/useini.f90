!=======================================================================
!> @brief Special purpose routine for initialisation based on user coding
!> @details All user-specified expressions must be written in standard
!! Fortran 90 or in the Fortran version available on your machine. <br>
!! To activate the statements, please remove the 'c' from the
!! first column. <br>
!! See the list of variables for the use of this routine.
      SUBROUTINE useini(k,mph)
!=======================================================================
!.....contact fire@avl.com
!-----
!-----------------------------------------------------------------------
!-----
      USE comm0
      USE comm1
      
      USE fire_mpi_var
      
      USE tmp_mod

      IMPLICIT NONE
!-----PARAMETER
      INTEGER K
      INTEGER MPH

      real(prec) :: dlt !< Time duration
      real(prec) :: dlx !< Length by X-corrdinate
      real(prec) :: dly !< Length by Y-corrdinate
      real(prec) :: dlz !< Length by Z-corrdinate
      integer :: nt !< number of elements by time
      ! integer :: nx !< number of elements by X-corrdinate
      ! integer :: ny !< number of elements by Y-corrdinate
      ! integer :: nz !<@var number of elements by Z-corrdinate
      integer :: mat !< Number of material
      real(prec) :: dsigma    !< define Deviation of velocities
      real(prec) :: dlength   !< define Integral Length Scale 
      real(prec) :: dtau      !< define Integral Time Scale 
      real(prec) :: dstd         !<
      real(prec) :: dmean         !<
      real(prec) :: dcell         !<
      real(prec) :: dtmp         !<
      integer, dimension(3)  :: nels            !< 
      real(prec), dimension(3)  :: dels         !<
      real(prec), dimension(3)  :: std_i         !<
      real(prec), dimension(3)  :: dmean_i         !<
      real(prec), dimension(3)  :: dxx         !<
      real(prec), dimension(3)  :: dvels         !<
      CHARACTER(len=256) :: txt1,txt2
      integer :: i,j
      logical :: lread2restart
!-----
!-----------------------------------------------------------------------
!-----
      mat = 1

      IF(I_USEINI == 1) THEN

        if (itst .LT. 2) then
          CALL CFD_PROFILE('PUSH','INIT_TURB')

          call mkl_set_num_threads(1)
          CALL FIO_UUMESS ('USEINI','I','Set number of OMP to 1','')

          ! set  time duration
          dlt = 1.0d-7
          ! set space Length
          dlx = 4.0d-3
          dly = dlx*2.
          dlz = dlx
          ! set number of nodes
          nx = 80
          ny = nx*2
          nz = nx

          !set number of timesteps
          nt = 1
          ntimes = nt

          !set number of Modes
          nmodes = 10000
          tcell  = nx * ny * nz
          CALL tmp_alloc()
          

          !set the Integral values
          dlength = 1.0d-2
          dsigma = 1.0d0*12.0d+0
          dtau = dlength / dsigma

          !generate arrays 
          dels(1) = dlx
          dels(2) = dly
          dels(3) = dlz

          nels(1) = nx
          nels(2) = ny 
          nels(3) = nz

          write(*,*) "work in 3D-space + Time"

          if (mpi_master) then

            !----------------------------------------------------------------
            lread2restart = .TRUE.
            lread2restart = .FALSE.
            !----------------------------------------------------------------
            IF (lread2restart) THEN
              CALL read_coef()
            ELSE
              do i = 1, ntimes
                in_time = nmodes * (i-1) + 1
                ! generate fields
                CALL gen_flow_3d(dels, nels, dsigma, dlength, dtau)
                write(*,'(2(a,i))') "USEINI: gen_flow_3d:",i,'/',ntimes
              end do
              CALL write_coef()
            END IF !(lread2restart) THEN
          endif

          i = ntimes * nmodes
          if (mode_mpi) THEN
            if (mpi_master) write(*,*) "CALL dbcast ac_m"
            CALL dbcast(0, ac_m(1,1), 3*i)
            if (mpi_master) write(*,*) "CALL dbcast as_m"
            CALL dbcast(0, as_m(1,1), 3*i)
            if (mpi_master) write(*,*) "CALL dbcast b_m"
            CALL dbcast(0, b_m(1,1), 3*i)
            if (mpi_master) write(*,*) "CALL dbcast c_m"
            CALL dbcast(0, c_m(1), i)
          END IF

          ! if (mpi_master) write(*,*) "CALL dglscatter ac_m"
          ! CALL dglscatter(ac_m,ac_m,3*i)
          ! if (mpi_master) write(*,*) "CALL dglscatter as_m"
          ! CALL dglscatter(as_m,as_m,3*i)
          ! if (mpi_master) write(*,*) "CALL dglscatter b_m"
          ! CALL dglscatter(b_m,b_m,3*i)
          ! if (mpi_master) write(*,*) "CALL dglscatter c_m"
          ! CALL dglscatter(c_m,c_m,i)
          ! if (mpi_master) write(*,*) "CALL dglscatter cs_m"
          ! CALL dglscatter(cs_m,cs_m,i)
          ! if (mpi_master) write(*,*) "CALL dglscatter dphi_m"
          ! CALL dglscatter(dphi_m,dphi_m,3*i)
          if (mpi_master) write(*,*) "Check "

          dxx(1) = 5.0d-5 ; dxx(2) = 5.0d-5 ; dxx(3) = 5.0d-5 
          dvels(:) = 0.0d0
          CALL set_vels_at_space_time(time,dxx,dvels)
          ! write(txt1,'(i,s,3es10.2)') iampro,": vels ",dvels(1),dvels(2),dvels(3)
          write(txt1,'(a,4es13.5)')" return First  ",dxx(1), dvels(1:3)
          txt2=''
          ! IF(iampro<2) CALL UCMESS('usedef','F',txt1,txt2)
        
          CALL FIO_UUMESS ('USEINI','I',txt1,txt2)

          dxx(1) = 5.0d-5+dlx ; dxx(2) = 5.0d-5 ; dxx(3) = 5.0d-5 
          dvels(:) = 0.0d0
          CALL set_vels_at_space_time(time,dxx,dvels)
          ! write(txt1,'(i,s,3es10.2)') iampro,": vels ",dvels(1),dvels(2),dvels(3)
          write(txt1,'(a,4es13.5)')" return Second ",dxx(1), dvels(1:3)
          txt2=''
          ! IF(iampro<2) CALL UCMESS('usedef','F',txt1,txt2)
          ! write(txt1,'(s)') "vels "!,dxx(1) !, dvels(1)
          ! write(txt2,'(s,i,l)') "Use MPI mode",iampro, mpi_slave
          CALL FIO_UUMESS ('USEINI','I',txt1,txt2)
          dxx(1) = 5.0d-5 ; dxx(2) = 5.0d-5 ; dxx(3) = 5.0d-5 
          dvels(:) = 0.0d0
          CALL set_vels_at_space_time(time+dlt,dxx,dvels)
          ! write(txt1,'(i,s,3es10.2)') iampro,": vels ",dvels(1),dvels(2),dvels(3)
          write(txt1,'(a,4es13.5)')" return Third  ",dxx(1), dvels(1:3)
          txt2=''
          ! IF(iampro<2) CALL UCMESS('usedef','F',txt1,txt2)
        
          CALL FIO_UUMESS ('USEINI','I',txt1,txt2)
          ! write(*,*)txt1
          ! CALL CFD_STOP()
          CALL CFD_PROFILE('POP','INIT_TURB')
        endif

        ! Define velocities
        CALL CFD_PROFILE('PUSH','GEN_VELS')

        ! CALL set_vels_at_time(time)

        do i=nsp(mat), nep(mat)
          CALL set_vels_at_space_time(time,xp(1:3,i),u(1:3,i))
        enddo
        ! u(:,:) =  0.0d0
        CALL exchng(u,3,1)
        CALL exchng(u,3,2)
        CALL exchng(u,3,3)
        ! call exchng3(u,3,3)

        CALL CFD_PROFILE('POP','GEN_VELS')

        ! write(*,*) "Calcs mean and deviation"
        do i=1, 3
          dmean = SUM(u(i,nsp(mat):nep(mat))) !/ DBLE(nep(mat)-nsp(mat)+1)
          dcell = DBLE(nep(mat)-nsp(mat)+1)
          CALL DGLSUM(dmean)
          CALL DGLSUM(dcell)
          dstd  = SUM((u(i,nsp(mat):nep(mat))-dmean)**2) !/ DBLE(nep(mat)-nsp(mat)+1)
          CALL DGLSUM(dstd)
          if (mpi_master) THEN
            dtmp = SQRT(dstd)/dcell
            write(txt1,'(a,i13,2es13.5)') "Mean and std vel: ",i, dmean/dcell, dtmp
            txt2 = ''
            CALL FIO_UUMESS ('USEINI','I',txt1,txt2)
          END IF ! (mpi_master) THEN
        enddo

!-----
      ELSE IF(I_USEINI == 2) THEN
!-----
      ELSE IF(I_USEINI == 99) THEN
         CALL USE_FORMULA('USEINI', 2, K, MPH, 0)
      END IF
!-----
      RETURN
      END
