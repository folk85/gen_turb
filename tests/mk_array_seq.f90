program mk_array_seq
  implicit none
  real, dimension(2,3) :: arr
  integer :: i, j 

  do i = 1, 2
    do j = 1, 3
      arr(i,j) = 3*(i-1) + j
    enddo
  enddo
  print *, arr(:,:)
end program mk_array_seq