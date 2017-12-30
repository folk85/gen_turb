program mk_array_seq
  implicit none
  real*8, dimension(2,3) :: arr
  real*8 :: et,tt(2), start, dend, dmid
  integer :: i, j 
  integer :: m = 100000
  integer :: n = 10000
  real*8, dimension(:),allocatable :: c,d,t
  real*8, dimension(:,:),allocatable :: a,b,e,f,x,u, y,v

  do i = 1, 2
    do j = 1, 3
      arr(i,j) = 3*(i-1) + j
    enddo
  enddo
  print *, arr(:,:)
  !
  do i = 1, 3
    do j = 1, 2
      arr(j,i) = 2*(i-1) + j
    enddo
  enddo
  print *, arr(:,:)

  allocate(c(m)); c(:) = 0.0d0
  allocate(d(m)); d(:) = 0.0d0
  allocate(t(m)); t(:) = 0.0d0

  allocate(a(1:3,1:m)); a(:,:) = 0.0d0
  allocate(b(1:3,1:m)); b(:,:) = 0.0d0

  allocate(x(1:3,1:n)); x(:,:) = 0.0d0
  allocate(u(1:3,1:n)); u(:,:) = 0.0d0

  allocate(e(1:m,1:3)); e(:,:) = 0.0d0
  allocate(f(1:m,1:3)); f(:,:) = 0.0d0
  allocate(y(1:n,1:3)); y(:,:) = 0.0d0
  allocate(v(1:n,1:3)); v(:,:) = 0.0d0

  ! call random_number()

  call cpu_time(start)
  do  i = 1, n
    ! c(1:m) = COS(b(1,1:m) * x(1,i) + b(2,1:m) * x(2,i) + b(3,1:m) * x(3,i) + d(1:m))
    c(1:m) = b(1,1:m) * x(1,i) + b(2,1:m) * x(2,i) + b(3,1:m) * x(3,i) + d(1:m)
    ! d(1:m) = COS(c(1:m))
    ! t(1:m) = SIN(c(1:m))
    CALL vssincos(m,c,d,t)

    do j=1,3
      ! u(j,i) = 2.0d0 * SUM(a(j,1:m) * COS(c(1:m)))
      u(j,i) = 2.0d0 * SUM(a(j,1:m) * d(1:m)+a(j,1:m) * t(1:m))
    enddo
  enddo

  call cpu_time(dmid)
  do i = 1, n
    ! c(1:m) = COS(e(1:m,1) * y(i,1) + e(1:m,2) * y(i,2) + e(1:m,3) * y(i,3) + d(1:m))
    c(1:m) = e(1:m,1) * y(i,1) + e(1:m,2) * y(i,2) + e(1:m,3) * y(i,3) + d(1:m)
    ! d(1:m) = COS(c(1:m))
    ! t(1:m) = SIN(c(1:m))
    CALL vssincos(m,c,d,t)
    do j = 1, 3
      ! v(i,j) = 2.0d0 * SUM(f(1:m,j) * COS(c(1:m)))
      ! v(i,j) = 2.0d0 * SUM(f(1:m,j) * c(1:m))
      v(i,j) = 2.0d0 * SUM(f(1:m,j) * d(1:m)+f(1:m,j) * t(1:m))
    end do      
  end do      

  call cpu_time(dend)
  call dtime(tt,et)
  print *, start, dmid, dend
  print *, et, tt(1:2)
end program mk_array_seq
