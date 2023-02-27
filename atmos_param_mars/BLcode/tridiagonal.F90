module tridiagonal_mod
!   Solves the tridiagonal system of equations.

!     The following schematic represents the system of equations solved,
!     where X is the solution.
!
!     | B(1)  A(1)   0     0                .......            0    |  |X(1)|   |D(1)|
!     | C(2)  B(2)  A(2)   0                .......            0    |  |X(2)|   |D(2)|
!     |  0    C(3)  B(3)  A(3)  0           .......            0    |  | .. |   | .. |
!     |  ..........................................                 |  | .. | = | .. |
!     |  ..........................................                 |  | .. |   | .. |
!     |                                  C(N-2) B(N-2) A(N-2)  0    |  | .. |   | .. |
!     |                                    0    C(N-1) B(N-1) A(N-1)|  | .. |   | .. |
!     |                                    0      0    C(N)   B(N)  |  |X(N)|   |D(N)|
! 
! For simplicity (?), A and C are assumed to be dimensioned the same size 
! as B, D, and X, although any input values for A(N) and C(1) are ignored.
! (some checks are needed here)
!
! If A is not present, it is assumed that the matrix (A,B.C) has not been changed 
! since the last call to tri_invert.
!
! To release memory, 


!--------------------------------------------------------------------------
real,    private, allocatable, dimension(:,:,:) :: e,g,cc
real,    private, allocatable, dimension(:,:)   :: bb
logical, private :: init_tridiagonal = .false.
!--------------------------------------------------------------------------

contains

!--------------------------------------------------------------------------

subroutine tri_invert(x,d,a,b,c)
!     Sets up and solves the tridiagonal system of equations.

implicit none

real, intent(out), dimension(:,:,:) :: x
real, intent(in),  dimension(:,:,:) :: d
real, optional,    dimension(:,:,:) :: a,b,c

real, dimension(size(x,1),size(x,2),size(x,3)) :: f
integer :: k

if(present(a)) then
    init_tridiagonal = .true.

    if(allocated(e))     deallocate(e)
    if(allocated(g))     deallocate(g)
    if(allocated(bb))    deallocate(bb)
    if(allocated(cc))    deallocate(cc)
    allocate(e (size(x,1),size(x,2),size(x,3)))
    allocate(g (size(x,1),size(x,2),size(x,3)))
    allocate(bb(size(x,1),size(x,2)))
    allocate(cc(size(x,1),size(x,2),size(x,3)))

    e(:,:,1) = - a(:,:,1)/b(:,:,1)
    a(:,:,size(x,3)) = 0.0

    do  k= 2,size(x,3)
        g(:,:,k) = 1/(b(:,:,k)+c(:,:,k)*e(:,:,k-1))
        e(:,:,k) = - a(:,:,k)*g(:,:,k)
    end do
    cc = c
    bb = 1.0/b(:,:,1)

end if

f(:,:,1) =  d(:,:,1)*bb
do k= 2, size(x,3)
    f(:,:,k) = (d(:,:,k) - cc(:,:,k)*f(:,:,k-1))*g(:,:,k)
end do

x(:,:,size(x,3)) = f(:,:,size(x,3))
do k = size(x,3)-1,1,-1
    x(:,:,k) = e(:,:,k)*x(:,:,k+1)+f(:,:,k)
end do

return
end subroutine tri_invert

!-----------------------------------------------------------------

subroutine close_tridiagonal
!     Releases memory used by the solver.
implicit none

deallocate(e)
deallocate(g)
deallocate(bb)
deallocate(cc)

return
end subroutine close_tridiagonal

!----------------------------------------------------------------



end module tridiagonal_mod
