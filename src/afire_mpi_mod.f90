!=======================================================================
!> Module with global mpi variables
MODULE fire_mpi_var
!=======================================================================
 IMPLICIT NONE
 SAVE

!> \addtogroup fire_mpi_var Global mpi variables
!> Global mpi variables
!> \{
 INTEGER :: iampro = 0   !< \ro_var current process number [-]
                         !! - Serial: 0
                         !! - MPI 1 process: 0
                         !! - MPI n processes (n > 1): 1 .. n
 INTEGER :: numpro = 0   !< \ro_var number of processes [-]
                         !! - Serial: 0
                         !! - MPI 1 process: 0
                         !! - MPI n processes (n > 1): n
 LOGICAL :: mpi_master  = .false. !< \ro_var process is the MPI master process (or the serial process)
 LOGICAL :: mpi_slave   = .false. !< \ro_var process is a MPI slave process
 LOGICAL :: mode_mpi    = .false. !< \ro_var calculation is a calculation with more than 1 process
 LOGICAL :: mode_serial = .false. !< \ro_var calculation is a calculation with 1 process
!> \}

END MODULE
