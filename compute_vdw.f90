PROGRAM compute_vdw
!
IMPLICIT NONE
!
! Global Parameters (consistent with Quantum Espresso)
!
INTEGER, PARAMETER :: DP = selected_real_kind(14,200)  !double-precision kind
REAL(DP), PARAMETER :: pi = 3.14159265358979323846_DP  !pi to double-precision
REAL(DP), PARAMETER :: bohr2Ang = 0.52917720859_DP     !bohr (au) to Angstrom conversion factor
REAL(DP), PARAMETER :: vdw_econv_thr = 1E-3_DP         !threshold for periodic convergence
!
! Local Parameters
!
REAL(DP), PARAMETER :: vHfree = 10.384749018802095_DP  !free atom volume (Hydrogen)
REAL(DP), PARAMETER :: vBfree = 50.835592620550344_DP  !free atom volume (Boron)
REAL(DP), PARAMETER :: vCfree = 38.560088325018114_DP  !free atom volume (Carbon)
REAL(DP), PARAMETER :: vOfree = 23.522305161381485_DP  !free atom volume (Oxygen)
!
REAL(DP), PARAMETER :: alphaHfree =  4.50_DP            !free atom static dipole polarizability (Hydrogen)
REAL(DP), PARAMETER :: alphaBfree = 21.00_DP            !free atom static dipole polarizability (Boron)
REAL(DP), PARAMETER :: alphaCfree = 12.00_DP            !free atom static dipole polarizability (Carbon)
REAL(DP), PARAMETER :: alphaOfree =  5.40_DP            !free atom static dipole polarizability (Oxygen)
!
REAL(DP), PARAMETER :: R0Hfree = 3.100_DP               !free atom radius (Hydrogen)
REAL(DP), PARAMETER :: R0Bfree = 3.890_DP               !free atom radius (Boron)
REAL(DP), PARAMETER :: R0Cfree = 3.590_DP               !free atom radius (Carbon)
REAL(DP), PARAMETER :: R0Ofree = 3.190_DP               !free atom radius (Oxygen)
!
REAL(DP), PARAMETER :: C6Hfree =  6.50_DP               !free atom C_6 coefficients (Hydrogen)
REAL(DP), PARAMETER :: C6Bfree = 99.50_DP               !free atom C_6 coefficients (Boron)
REAL(DP), PARAMETER :: C6Cfree = 46.60_DP              !free atom C_6 coefficients (Carbon)
REAL(DP), PARAMETER :: C6Ofree = 15.60_DP              !free atom C_6 coefficients (Oxygen)
!
REAL(DP), PARAMETER :: sR = 0.94_DP                    !damping function parameter #1 (PBE)
REAL(DP), PARAMETER :: d = 20.0_DP                     !damping function parameter #2 (All XC functionals)
!
!
! Local Variables
!
INTEGER :: ia,ib,nat,n1,n2,n_period,Dori,Dstart,iD
CHARACTER(LEN=2), DIMENSION(:), ALLOCATABLE :: atsym
REAL(DP), DIMENSION(:), ALLOCATABLE :: vhirsh,ratio,R0eff,alphaeff,C6eff
REAL(DP), DIMENSION(:,:), ALLOCATABLE :: atxyz
REAL(DP) :: RAB,R0effAB,feffAB,C6effAB,EvdW_period,Edft,s3,lx,x2,y2,z2,b2,deltax,DD,etsvdw0,etsvdwm1,etsvdwp1
LOGICAL :: periodic_converged
!
! Open files used during program execution...
!
OPEN(UNIT=10,FILE="input.in",ACTION='READ',FORM='FORMATTED')
OPEN(UNIT=20,FILE="output.out",ACTION='WRITE',FORM='FORMATTED')
!
! Read in number of atoms from input.in...
!
READ(10,*) nat
READ(10,*) lx
READ(10,*) Dstart
READ(10,*) id
lx = lx/bohr2Ang
!
! Allocate and Initialize atomic symbol, Cartesian coordinate, and
! Hirshfeld volume/ratio arrays...
!
ALLOCATE(atsym(nat)); atsym=''
ALLOCATE(atxyz(nat,3)); atxyz=0.0_DP 
ALLOCATE(vhirsh(nat)); vhirsh=0.0_DP 
ALLOCATE(ratio(nat)); vhirsh=0.0_DP 
ALLOCATE(R0eff(nat)); R0eff=0.0_DP 
ALLOCATE(alphaeff(nat)); alphaeff=0.0_DP 
ALLOCATE(C6eff(nat)); C6eff=0.0_DP 
!
! Read in Cartesian coordinates and Hirshfeld volume/ratio arrays...
!
DO ia=1,nat
  !
  READ(10,*) atsym(ia),atxyz(ia,1),atxyz(ia,2),atxyz(ia,3),vhirsh(ia)
  !
  ! Convert atomic coordinates from Angstroms into Bohr...
  !
  atxyz(ia,1)=atxyz(ia,1)/bohr2Ang
  atxyz(ia,2)=atxyz(ia,2)/bohr2Ang
  atxyz(ia,3)=atxyz(ia,3)/bohr2Ang
  !
END DO !ia
!
! Construct Hirshfeld ratio, effective vdW radius, effective alpha, effective C6AA... 
!
DO ia=1,nat
  ! 
  IF (atsym(ia).EQ.'H') THEN
    !
    ratio(ia)=vhirsh(ia)/vHfree
    R0eff(ia)=ratio(ia)**(1.0_DP/3.0_DP)*R0Hfree
    alphaeff(ia)=ratio(ia)*alphaHfree
    C6eff(ia)=ratio(ia)**(2.0_DP)*C6Hfree
    !
  ELSE IF (atsym(ia).EQ.'B') THEN
    !
    ratio(ia)=vhirsh(ia)/vBfree
    R0eff(ia)=ratio(ia)**(1.0_DP/3.0_DP)*R0Bfree
    alphaeff(ia)=ratio(ia)*alphaBfree
    C6eff(ia)=ratio(ia)**(2.0_DP)*C6Bfree
    !
  ELSE IF (atsym(ia).EQ.'C') THEN
    !
    ratio(ia)=vhirsh(ia)/vCfree
    R0eff(ia)=ratio(ia)**(1.0_DP/3.0_DP)*R0Cfree
    alphaeff(ia)=ratio(ia)*alphaCfree
    C6eff(ia)=ratio(ia)**(2.0_DP)*C6Cfree
    !
  ELSE IF (atsym(ia).EQ.'O') THEN
    !
    ratio(ia)=vhirsh(ia)/vOfree
    R0eff(ia)=ratio(ia)**(1.0_DP/3.0_DP)*R0Ofree
    alphaeff(ia)=ratio(ia)*alphaOfree
    C6eff(ia)=ratio(ia)**(2.0_DP)*C6Ofree
    !
  END IF
  !
END DO  !ia
!
! Compute distance, effective vdW radius, damping function, and energy on the fly...
! Using finite-difference to compute the gradient
!
deltax=0.001_DP

DO Dori=Dstart,40,id
   !
   DD=Dori/bohr2Ang
   CALL EvdW( 0.0_DP       ,etsvdw0 )
   CALL EvdW( 1.0_DP*deltax,etsvdwp1)
   CALL EvdW(-1.0_DP*deltax,etsvdwm1)
   b2=DD/etsvdw0*(etsvdwm1-etsvdwp1)/(2.0_DP*deltax)
   Write(20,'(4F25.15)') DD*bohr2Ang,b2
END DO
!
! Deallocate arrays used during program execution...
!
DEALLOCATE(atsym)
DEALLOCATE(atxyz)
DEALLOCATE(vhirsh)
DEALLOCATE(ratio)
DEALLOCATE(R0eff)
DEALLOCATE(alphaeff)
DEALLOCATE(C6eff)
!
! Close files used during program execution...
!
CLOSE(10)
CLOSE(20)
!
    CONTAINS
      !
      SUBROUTINE EvdW(deltad,etsvdw)
          IMPLICIT NONE
          REAL(DP),INTENT(IN)  ::deltad
          REAL(DP),INTENT(OUT) :: etsvdw
          REAL(DP) :: EvdW_period
          etsvdw=0.0_DP

          s3=1.73205080757_DP
          n_period = 0
          periodic_converged = .FALSE.
          !
              DO WHILE (.NOT.periodic_converged)
              
                EvdW_period = 0.0_DP
                  
                !$omp parallel private(ia,ib,n1,n2,RAB,C6effAB,R0effAB,feffAB,x2,y2),reduction(+:EvdW_period)
                !$omp do

                DO ia=1,nat
                  !
                  DO ib=1,nat
                    !
                    DO n1 = -n_period, n_period
                     !
                      DO n2 = -n_period, n_period
                        !
                        IF((ABS(n1).EQ.n_period) .OR. (ABS(n2).EQ.n_period)) THEN 
                          !
                          x2=atxyz(ia,1)-(DBLE(n2)*0.5_DP+DBLE(n1))*lx
                          y2=atxyz(ia,2)-DBLE(n2)*0.5_DP*s3*lx
                          z2=0
                          RAB=DSQRT((x2-atxyz(ib,1))**(2.0_DP)+(y2-atxyz(ib,2))**(2.0_DP)+(DD+deltad)**(2.0_DP))
                          R0effAB=R0eff(ia)+R0eff(ib)
                          feffAB=1.0_DP/(1.0_DP+DEXP(-d*(RAB/(sR*R0effAB)-1.0_DP)))
                          C6effAB=(2.0_DP*C6eff(ia)*C6eff(ib)) & 
                          /((alphaeff(ib)/alphaeff(ia))*C6eff(ia)+(alphaeff(ia)/alphaeff(ib))*C6eff(ib))
                          !
                          EvdW_period=EvdW_period+(C6effAB * feffAB / RAB**(6.0_DP))
                          !
                        END IF !end if(n1=n_period)
                        !
                      END DO !n2
                      !
                    END DO !end perido n1
                    !
                  END DO !ib 
                  ! 
                END DO  !ia
                !$omp end do 
                !$omp end parallel

                !
                IF (n_period.EQ.0) THEN
                  etsvdw = EvdW_period
                  n_period = n_period + 1
                !
                ELSE IF (ABS(EvdW_period)/etsvdw < vdw_econv_thr) THEN
                  periodic_converged=.TRUE.
                  etsvdw = etsvdw + EvdW_period
                !
                ELSE
                  etsvdw = etsvdw + EvdW_period
                  n_period = n_period + 1
                !
                END IF
                !
              END DO
              !
            etsvdw=-1.0_DP*etsvdw 
        !
      END SUBROUTINE EvdW

END PROGRAM compute_vdw
