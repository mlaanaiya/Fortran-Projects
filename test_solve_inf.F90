program test_solve

  implicit none

  integer :: i, j, ierr, n
  double precision, dimension (:,:), allocatable :: L
  double precision, dimension (:), allocatable :: x, b
  double precision :: backward_error
  real :: start, finish

  write(*,*) 'n?'
  99 read(*,*) n
  if (n<=0) then
    print*, 'Entrez un entier strictement positif'
    goto 99
  endif

  ! Initialization: L is lower triangular
  write(*,*) 'Initialization...'
  write(*,*)
  
  allocate(L(n,n), stat=ierr)
  if(ierr.ne.0) then
    write(*,*)'Could not allocate L(n,n) with n=',n
    goto 999
  end if

  allocate(x(n), stat=ierr)
  if(ierr.ne.0) then
    write(*,*)'Could not allocate x(n) with n=',n
    goto 999
  end if

  allocate(b(n), stat=ierr)
  if(ierr.ne.0) then
    write(*,*)'Could not allocate b(n) with n=',n
    goto 999
  end if

  L = 0.D0
  do i = 1, n  
    L(i,i) = n + 1.D0
    do j = 1, i-1
      L(i,j) = 1.D0
    end do
  end do
  b = 1.D0

  ! Left-looking triangular solve of Lx=b
  write(*,*) 'Solving with a left-looking method...'

  call cpu_time(start)
    call left_looking_solve(L,x,b,n)
  call cpu_time(finish)
  print '("Time = ",f6.3," seconds.")',finish-start

  print*,'l´erreur est de : ', backward_error(L,x,b,n)


 
  ! Right-looking triangular solve of Lx=b
  write(*,*) 'Solving with a right-looking method...'

  call cpu_time(start)
    call right_looking_solve(L,x,b,n)
  call cpu_time(finish)
  print '("Time = ",f6.3," seconds.")',finish-start

  print*,'l´erreur est de : ', backward_error(L,x,b,n)



999 if(allocated(L)) deallocate(L)
    if(allocated(x)) deallocate(x)
    if(allocated(b)) deallocate(b)

end program test_solve


! Implement sub-programs left_looking_solve
  subroutine left_looking_solve(L,x,b,n)
  implicit none
  integer, intent(in) :: n
  double precision, intent(in), dimension(n) :: b
  double precision, intent(in), dimension(n,n) :: L
  double precision, intent(inout), dimension(n) :: x
  integer :: i
  integer :: j
  x = b
  b1 : do j = 1,n
         b2 : do i = 1, j-1
                x(j) = x(j) - L(j,i)*x(i)
              end do b2
         x(j) = x(j)/L(j,j)
       end do b1
       if (n<=10) then
           write(*,*) 'Le vecteur x : ',x
       end if
  end subroutine left_looking_solve

! Implement sub-programs right_looking_solve
  subroutine right_looking_solve(L,x,b,n)
  implicit none
  integer, intent(in) :: n
  double precision, intent(in), dimension(n) :: b
  double precision, intent(in), dimension(n,n) :: L
  double precision, intent(out), dimension(n) :: x
  integer :: i
  integer :: j
  x = b
  b1 : do j = 1,n
         x(j) = x(j)/L(j,j)
         b2 : do i = j+1, n
                x(i) = x(i) - L(i,j)*x(j)
              end do b2
       end do b1
       if (n<=10) then
           write(*,*)'Le vecteur x : ',x
       end if
  end subroutine right_looking_solve

! Implement sub-programs backward_error
  double precision function backward_error(L,x,b,n)
  implicit none
  integer, intent(in) :: n
  double precision, intent(in), dimension(n) :: b
  double precision, intent(in), dimension(n,n) :: L
  double precision, intent(in), dimension(n) :: x
  double precision :: first_vector
  double precision :: second_vector
  first_vector = dot_product(matmul(L, x) - b, matmul(L, x) - b)
  second_vector = dot_product(b,b)
  backward_error = SQRT(first_vector)/SQRT(second_vector)
  return
  end function backward_error 

!Comparaison:
!Le temps d'exécution avec la right method est inférieur à celui 
!de la left method. L'erreur est le même dans les deux méthodes.



