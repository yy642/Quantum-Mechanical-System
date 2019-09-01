!ifort -O3 -o calculate_chi_para-ifort.bin -qopenmp response_function_threaded.f90
!export OMP_NUM_THREADS=3
PROGRAM response_function
!
IMPLICIT NONE
!
! Global Parameters (consistent with Quantum Espresso)
!
INTEGER, PARAMETER  :: DP = selected_real_kind(14,200)  !double-precision kind
REAL(DP), PARAMETER :: pi = 3.14159265358979323846_DP   !pi to double-precision
!
! Local Parameters
!
INTEGER, PARAMETER  :: N = 5                            !number of grid points along each direction
INTEGER, PARAMETER  :: Ngrid = N*N*N                    !total number of grid points
INTEGER, PARAMETER  :: Nvirt = 5                        !highest virtual included
!REAL(DP), PARAMETER :: omega = 12.7062_DP                   !frequency (real part of i*omega)
!
! Local Variables
!
INTEGER :: i,j,n1,n2,n3,ix,iy,iz,jx,jy,jz
INTEGER, DIMENSION(:,:), ALLOCATABLE :: imat
REAL(DP), DIMENSION(:,:), ALLOCATABLE :: hmat,chi
REAL(DP) :: tmp1,tmp2,tmp3,tmp4,tmp5,x1,x2,x,t0,t1,omega
DOUBLE COMPLEX :: fp,fm 
!
CALL CPU_TIME(t0)
!
! Open files used during program execution...
!
OPEN(UNIT=10,FILE="imat.in",ACTION='READ',FORM='FORMATTED')
OPEN(UNIT=20,FILE="hmat.in",ACTION='READ',FORM='FORMATTED')
OPEN(UNIT=30,FILE="chi.out",ACTION='WRITE',FORM='FORMATTED')
OPEN(UNIT=40,FILE="input.in",ACTION='READ',FORM='FORMATTED')
!
! Allocate and read in grid from imat.in...
!
READ(40,*) omega 
!
ALLOCATE(imat(Ngrid,3)); imat=0
!
DO i=1,Ngrid
  !
  READ(10,*) imat(i,1),imat(i,2),imat(i,3)
  !
END DO !i
!
! Allocate and read in QHO wavefunctions from hmat.in...
!
ALLOCATE(hmat(Nvirt+1,N)); hmat=0.0_DP
!
DO i=1,Nvirt+1
  !
  READ(20,*) hmat(i,1),hmat(i,2),hmat(i,3),hmat(i,4),hmat(i,5)
  !
END DO !i
!
! Allocate and form chi...
!
ALLOCATE(chi(Ngrid,Ngrid)); chi=0.0_DP
!
!$omp parallel private(i,ix,iy,iz,tmp1,j,jx,jy,jz,tmp2,n1,x1,tmp3,n2,x2,tmp4,n3,x,tmp5)
!$omp do
DO i=1,Ngrid
  !
  ix = imat(i,1); iy = imat(i,2); iz = imat(i,3)
  !
  tmp1 = hmat(1,ix+3) * hmat(1,iy+3) * hmat(1,iz+3)
  !
  DO j=i,Ngrid
    !
    jx = imat(j,1); jy = imat(j,2); jz = imat(j,3)
    !
    tmp2 = tmp1 * hmat(1,jx+3) * hmat(1,jy+3) * hmat(1,jz+3)
    !
    DO n1=1,Nvirt+1
      !
      x1 = DBLE(n1) - 1.0_DP
      !
      tmp3 = tmp2 * hmat(n1,ix+3) * hmat(n1,jx+3)
      !
      DO n2=1,Nvirt+1
        !
        x2 = x1 + DBLE(n2) - 1.0_DP
        !
        tmp4 = tmp3 * hmat(n2,iy+3) * hmat(n2,jy+3)
        !
        DO n3=1,Nvirt+1
          !
          x = x2 + DBLE(n3) - 1.0_DP
          !
          fp = CMPLX(+1.0_DP*x,+1.0_DP*omega)
          fm = CMPLX(+1.0_DP*x,-1.0_DP*omega)
          !
          tmp5 = tmp4 * hmat(n3,iz+3) * hmat(n3,jz+3) * (1.0_DP/fp + 1.0_DP/fm)
          !
          chi(i,j) = chi(i,j) - tmp5     !negative sign from frequency here...
          !
          !
        END DO !n3
        !
      END DO !n2
      !
    END DO !n1
    !
  END DO !j
  !
END DO !i
!$omp end do
!$omp end parallel
!
! Symmeterize chi...
!
DO i=1,Ngrid-1
  !
  DO j=i+1,Ngrid
    !
    chi(j,i) = chi(i,j)
    !
  END DO !j
  !
END DO !i
!
! Write final chi to chi.out...
!
DO i=1,Ngrid
  !
  DO j=1,Ngrid
    !
    WRITE(30,'(ES25.15)') chi(i,j)
    !
  END DO !j
  !
END DO !i
!
! Deallocate arrays used during program execution...
!
DEALLOCATE(imat)
DEALLOCATE(hmat)
DEALLOCATE(chi)
!
! Close files used during program execution...
!
CLOSE(10)
CLOSE(20)
CLOSE(30)
CLOSE(40)
!
CALL CPU_TIME(t1)
!
PRINT '("Time = ",F12.3," seconds.")', t1-t0
!
END PROGRAM response_function
