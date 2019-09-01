PROGRAM integral 
!
!Using casimir polder integral to compute dispersion energy using frequency-dependent polarzability
! input: eigenmat.in, the eigenvalue of the target system.
! input: norder.in, target order of the energyexpansion. n=2 to infinty.
! output: the energy
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
! Local Variables
!
INTEGER :: LD,Nvirt,i,j,n1,n2,n3,ix,iy,iz,jx,jy,jz,iunit,INFO,windex,iorder,norder
INTEGER, DIMENSION(:,:), ALLOCATABLE :: imat
REAL(DP), DIMENSION(:,:), ALLOCATABLE :: hmat,chiring,chiatom,chi,Vmat,chi2,chiV,chiVnew
REAL(DP), DIMENSION(:), ALLOCATABLE :: eigenmat,WORK
REAL(DP) :: tmp1,tmp2,tmp3,tmp4,tmp5,x1,x2,x,t0,t1,tmpd,tmpd2,s2,trace,total_energy,xr,yr,zr,xr2,yr2,omega,d,energyn
DOUBLE COMPLEX :: fp,fm 
REAL(DP) :: casimir_omega_weight(1:20)
 casimir_omega_weight(0) = 0.0000000
 casimir_omega_weight(1) = 0.0786611
 casimir_omega_weight(2) = 0.0796400
 casimir_omega_weight(3) = 0.0816475
 casimir_omega_weight(4) = 0.0847872
 casimir_omega_weight(5) = 0.0892294
 casimir_omega_weight(6) = 0.0952317
 casimir_omega_weight(7) = 0.1031720
 casimir_omega_weight(8) = 0.1136050
 casimir_omega_weight(9) = 0.1273500
casimir_omega_weight(10) = 0.1456520
casimir_omega_weight(11) = 0.1704530
casimir_omega_weight(12) = 0.2049170
casimir_omega_weight(13) = 0.2544560
casimir_omega_weight(14) = 0.3289620
casimir_omega_weight(15) = 0.4480920
casimir_omega_weight(16) = 0.6556060
casimir_omega_weight(17) = 1.0659600
casimir_omega_weight(18) = 2.0635700
casimir_omega_weight(19) = 5.6851000
casimir_omega_weight(20) = 50.955800
!
!
CALL CPU_TIME(t0)
!
! Open files used during program execution...
!
OPEN(UNIT=10,FILE="eigenmat.in",ACTION='READ',FORM='FORMATTED')
OPEN(UNIT=20,FILE="norder.in",ACTION='READ',FORM='FORMATTED')
OPEN(UNIT=60,FILE="energy_nthorder.out",ACTION='WRITE',FORM='FORMATTED')
!
READ(10,*) LD 
ALLOCATE(eigenmat(20*LD)); eigenmat=0.0_DP
DO i=1,20*LD
  !
  READ(10,*) (eigenmat(i))
  !
END DO !i
READ(20,*) norder
!!
!! loop over all w...
!!
total_energy=0.0_DP
!
DO iunit=1,20 ! loop w
  !
  trace=0.0_DP
  !
  windex = (iunit-1)*LD
  !
  DO iorder=2,norder
    !
    DO i=1, LD
      !
      energyn = eigenmat(i+windex)**iorder
      trace = trace + energyn/iorder 
      !
    END DO !i
  END DO !i
  !
  total_energy = total_energy + trace * casimir_omega_weight(iunit)
  !
END DO !w
!
total_energy = total_energy / 40.0_DP *(-1.0_DP)
!
WRITE(60,'("norder=", I10, " total_energy= " ES25.15)') norder, total_energy 
!
! Deallocate arrays used during program execution...
!
DEALLOCATE(eigenmat)
!
! Close files used during program execution...
!
CLOSE(10)
CLOSE(20)
CLOSE(60)
!
CALL CPU_TIME(t1)
!
PRINT '("Time = ",F12.3," seconds.")', t1-t0
!
END PROGRAM integral
