program newtonMethod
    implicit none
    integer :: i, n
    real(8)    :: y(30), x(30), r(30), J(30,4)
    real(8)    :: beta(4), betaOld(4),sol(4), inverz(4,4)
    real(8)    :: eta0, etainf, tau, rr, eta, S
    open(10,file = "podatki2.dat",status = "old")
    do i = 1, 30
        read(10,*) y(i), x(i)
    end do
    close(10)

    !zacetni priblizek
    eta0   = 5.0
    etainf = 6.0
    tau    = 0.005
    rr     = 0.6

    beta = (/eta0,etainf,tau,rr/)
    print*, "beta", beta

    !zacetna vrednost rešitve
    sol(1:4)=(/1.0e0,1.0e0,1.0e0,1.0e0/)

    n = 0
    do while (any(sol>1.0d-6) .and. n < 100)

        !vrednosti funkcij, morajo biti enake 0, tu primerjave izmerjene z izracunano
        do i = 1, 30
            r(i) = y(i) - eta(x(i),beta(1),beta(2),beta(3),beta(4))
        end do
        S = sum(r(:)**2)

        print*,"S",S

        !izracunan Jacobijeve matrike
        do i = 1, 30
            J(i,1) = exp(-(x(i)/beta(3))**beta(4))
            J(i,2) = 1.0 - exp(-(x(i)/beta(3))**beta(4))
            J(i,3) = (beta(4)*(beta(1)-beta(2))*exp(-(x(i)/beta(3))**beta(4))*(x(i)/beta(3))**beta(4))/beta(3)
            J(i,4) = (beta(1)-beta(2))*log(x(i)/beta(3))*exp(-(x(i)/beta(3))**beta(4))*(x(i)/beta(3))**beta(4)
        end do

        betaOld = beta

        call inverse(matmul(transpose(J),J),inverz,4)

        beta = beta - matmul(matmul(inverz,transpose(J)),r)

        print*, "beta", beta
        sol=abs((beta - betaOld)/beta)

        N=N+1
    end do

end program newtonMethod

function eta(x,eta0,etainf,tau,r)
    implicit none
    real(8):: eta,eta0,etainf,tau,r,x

    eta = eta0-(eta0-etainf)*(1.0e0-exp(-((x/tau)**r)))

end function eta

subroutine inverse(a,c,n)
!============================================================
! Inverse matrix
! Method: Based on Doolittle LU factorization for Ax=b
! Alex G. December 2009
!-----------------------------------------------------------
! input ...
! a(n,n) - array of coefficients for matrix A
! n      - dimension
! output ...
! c(n,n) - inverse matrix of A
! comments ...
! the original matrix a(n,n) will be destroyed
! during the calculation
!===========================================================
implicit none
integer, intent(in):: n
real(8) :: a(n,n)
real(8),intent(out):: c(n,n)
real(8):: L(n,n), U(n,n), b(n), d(n), x(n)
real(8):: coeff
integer:: i, j, k

! step 0: initialization for matrices L and U and b
! Fortran 90/95 allows such operations on matrices
L=0.0d0
U=0.0d0
b=0.0d0

! step 1: forward elimination
do k=1, n-1
   do i=k+1,n
      coeff=a(i,k)/a(k,k)
      L(i,k) = coeff
      do j=k+1,n
         a(i,j) = a(i,j)-coeff*a(k,j)
      end do
   end do
end do

! Step 2: prepare L and U matrices
! L matrix is a matrix of the elimination coefficient
! + the diagonal elements are 1.0
do i=1,n
  L(i,i) = 1.0d0
end do
! U matrix is the upper triangular part of A
do j=1,n
  do i=1,j
    U(i,j) = a(i,j)
  end do
end do

! Step 3: compute columns of the inverse matrix C
do k=1,n
  b(k)=1.0
  d(1) = b(1)
! Step 3a: Solve Ld=b using the forward substitution
  do i=2,n
    d(i)=b(i)
    do j=1,i-1
      d(i) = d(i) - L(i,j)*d(j)
    end do
  end do
! Step 3b: Solve Ux=d using the back substitution
  x(n)=d(n)/U(n,n)
  do i = n-1,1,-1
    x(i) = d(i)
    do j=n,i+1,-1
      x(i)=x(i)-U(i,j)*x(j)
    end do
    x(i) = x(i)/u(i,i)
  end do
! Step 3c: fill the solutions x(n) into column k of C
  do i=1,n
    c(i,k) = x(i)
  end do
  b(k)=0.0d0
end do
end subroutine inverse
