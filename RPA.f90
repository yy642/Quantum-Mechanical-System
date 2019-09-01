!ifort -O3 -o ring-atom.o ring-atom.f90 -lblas -llapack
!export OMP_NUM_THREADS=3
PROGRAM RPA 
!
IMPLICIT NONE
!
! Global Parameters (consistent with Quantum Espresso)
!
INTEGER, PARAMETER  :: DP = selected_real_kind(14,200)  !double-precision kind
REAL(DP), PARAMETER :: pi = 3.14159265358979323846_DP   !pi to double-precision
REAL(DP), PARAMETER :: NcoverNg               = 0.5d0 !YY: twice as many grid points than carbon atoms...
REAL(DP), PARAMETER :: beta                    = 0.83_DP 
REAL(DP), PARAMETER :: bohr = 0.52917721d0
!
! Local Parameters
!
! Local Variables
!
INTEGER :: n_atoms,io_line,i_row,n_atoma,i_col,i_index,j_index,i_freq,i_myatom,INFO,LWORK
REAL(DP), DIMENSION(:,:), ALLOCATABLE :: alpha0_matrix,coords,T,alpha0T,relay_matrix,coupled_atom_pol,alpha1_matrix,alpha1T 
REAL(DP), dimension(3,3) ::TPP,mol_pol_tensor
REAL(DP), dimension(3) :: dxyz
REAL(DP), DIMENSION(:), ALLOCATABLE :: hirshfeld_volume,R_p,Rvdw_iso,alpha_omega,alpha_eff,C6_eff,eigenmat1,eigenmat0,WORK
REAL(DP) :: C6_free
REAL(DP) :: alpha_free
REAL(DP) :: R_vdw_free 
REAL(DP) :: r_ij 
REAL(DP) :: r_pp
REAL(DP) :: Rvdw12 
REAL(DP) :: energy0 
REAL(DP) :: energy1 
REAL(DP) :: freq_energy0 
REAL(DP) :: freq_energy1 

character*2,allocatable,dimension(:) :: atom_name
real*8 :: casimir_omega(0:20)
real*8 :: casimir_omega_weight(0:20)      
!
!
! Open files used during program execution...
!
OPEN(UNIT=999,FILE="input.in",ACTION='READ',FORM='FORMATTED')
OPEN(UNIT=20,FILE="eigen0.out",ACTION='WRITE',FORM='FORMATTED')
OPEN(UNIT=30,FILE="eigen1.out",ACTION='WRITE',FORM='FORMATTED')
!OPEN(UNIT=40,FILE="alpha0T.out",ACTION='WRITE',FORM='FORMATTED')
!OPEN(UNIT=50,FILE="eigen.out",ACTION='WRITE',FORM='FORMATTED')
!
READ(999,*)n_atoms
READ(999,*)
!---------------------------------------------------------------------
IF(.NOT. ALLOCATED(coords))                 ALLOCATE(coords(3,n_atoms))
IF(.NOT. ALLOCATED(atom_name))              ALLOCATE(atom_name(n_atoms))
IF(.NOT. ALLOCATED(hirshfeld_volume))       ALLOCATE(hirshfeld_volume(n_atoms))
  DO i_index=1,n_atoms,1
    READ(999,*,IOSTAT=io_line) atom_name(i_index), coords(:,i_index),hirshfeld_volume(i_index)
  END DO
CLOSE(999)
coords = coords/bohr!
LWORK = n_atoms*9-1
!print '(I10)',LWORK
! Allocate and read in grid from imat.in...
!
allocate(alpha0_matrix(3*n_atoms,3*n_atoms))
allocate(alpha1_matrix(3*n_atoms,3*n_atoms))
allocate(relay_matrix(3*n_atoms,3*n_atoms))
allocate(alpha0T(3*n_atoms,3*n_atoms))
allocate(alpha1T(3*n_atoms,3*n_atoms))
allocate(T(3*n_atoms,3*n_atoms))
allocate(R_p(n_atoms))
allocate(Rvdw_iso(n_atoms))
allocate(alpha_omega(n_atoms))
allocate(alpha_eff(n_atoms))
allocate(C6_eff(n_atoms))
allocate(WORK(LWORK)); WORK=0.0_DP
allocate(eigenmat0(3*n_atoms)); eigenmat0=0.0_DP
allocate(eigenmat1(3*n_atoms)); eigenmat1=0.0_DP
allocate(coupled_atom_pol(20,n_atoms)); coupled_atom_pol=0.0_DP
!
casimir_omega        = 0.d0
casimir_omega_weight = 0.d0 
call gauss_legendre_grid()
energy0=0.d0
energy1=0.d0
!
WRITE(20,'(I5)') n_atoms*3
WRITE(30,'(I5)') n_atoms*3

do i_freq=1,20,1
  freq_energy0=0.d0
  freq_energy1=0.d0

  !! zeroout array before getting actual frequency dependent
  !parameters
  relay_matrix= 0.0d0
  R_p = 0.0d0
  alpha_omega = 0.0d0
  alpha0_matrix= 0.0d0
  alpha0T=0.0d0
  T=0.0d0
  mol_pol_tensor= 0.0d0
  Rvdw_iso=0.d0
  ! loop over atoms
  do i_myatom=1,n_atoms,1
    !
    call get_vdw_param(atom_name(i_myatom),C6_free,alpha_free,R_vdw_free)
    call get_alpha_omega_and_Rp(hirshfeld_volume(i_myatom),C6_free,&
    alpha_free,casimir_omega(i_freq),alpha_omega(i_myatom),R_p(i_myatom))
    Rvdw_iso(i_myatom)= (hirshfeld_volume(i_myatom)**0.333333333333333333333333333d0)*R_vdw_free !YY:when is this quantity used?
    !
  enddo ! end loop over atoms
!
! calculate alpha1
  call calculate_scs_matrix()
!
! calculate alpha0
!
  do i_row=1,n_atoms,1 !#1
    !
    do i_col=i_row,n_atoms,1 !#2
      !
      if(i_row.eq.i_col) then  !$1
        !
        do i_index=1,3,1
          !
          do j_index=1,3,1
            !
            if(i_index.eq.j_index) then
              alpha0_matrix(3*i_row-3+i_index,3*i_col-3+j_index)=-1.d0*alpha_omega(i_row)
            else
              alpha0_matrix(3*i_row-3+i_index,3*i_col-3+j_index)=0.d0
            endif
            !
          enddo !j_index
          !
        enddo !i_index
        !
      end if !$1
      !
    end do !i_col
    !
  end do !i_row

  ! compute T 
  do i_row=1,n_atoms,1 !#1
    !
    do i_col=1,n_atoms,1 !#2
      !
      if (i_row .ne. i_col) then
        !
        TPP=0.d0
        dxyz(:) = coords(:,i_col)-coords(:,i_row)
        r_ij = dsqrt((dxyz(1))**2.0d0 +(dxyz(2))**2.0d0 + (dxyz(3))**2.0d0 )
        r_pp = dSqrt(R_p(i_row)**2 + R_p(i_col)**2)
        Rvdw12 = Rvdw_iso(i_row) +  Rvdw_iso(i_col)
        !
        call bare_dipole(dxyz,r_ij,Rvdw12,beta,TPP)
        !
        do i_index=1,3,1
          !
          do j_index=1,3,1
            !
            T(3*i_row-3+i_index,3*i_col-3+j_index)=TPP(i_index,j_index)
            T(3*i_col-3+j_index,3*i_row-3+i_index)=TPP(i_index,j_index)
            !
          enddo
          !
        enddo
        !
      end if
      !
    enddo   !#2
    !
  enddo  !#1
  CALL dgemm('N','N',3*n_atoms,3*n_atoms,3*n_atoms,1.0_DP,alpha0_matrix,3*n_atoms,T,3*n_atoms,0.0_DP,alpha0T,3*n_atoms)
  CALL DSYEV('N','U',3*n_atoms, alpha0T,3*n_atoms, eigenmat0, WORK, LWORK, INFO)

  CALL dgemm('N','N',3*n_atoms,3*n_atoms,3*n_atoms,1.0_DP,relay_matrix,3*n_atoms,T,3*n_atoms,0.0_DP,alpha1T,3*n_atoms)
  CALL DSYEV('N','U',3*n_atoms, alpha1T,3*n_atoms, eigenmat1, WORK, LWORK, INFO)


  do i_index=1, 3*n_atoms,1
    freq_energy0 = freq_energy0 + eigenmat0(i_index)*eigenmat0(i_index)
    freq_energy1 = freq_energy1 + eigenmat1(i_index)*eigenmat1(i_index)
  end do
  !
  energy0 = energy0 + freq_energy0*casimir_omega_weight(i_freq)
  energy1 = energy1 + freq_energy1*casimir_omega_weight(i_freq)

  !
  do i_index=1, 3*n_atoms,1
    !
    WRITE(20,'(F25.15)') eigenmat0(i_index)
    WRITE(30,'(F25.15)') eigenmat1(i_index)
    !
  end do
  !
end do !i_freq
!
print '(2F25.15)', energy0, energy1 
!
! Deallocate arrays used during program execution...
!
deallocate(alpha0_matrix)
deallocate(alpha1_matrix)
deallocate(alpha0T)
deallocate(alpha1T)
deallocate(T)
deallocate(R_p)
deallocate(Rvdw_iso)
deallocate(alpha_omega)
deallocate(alpha_eff)
deallocate(C6_eff)
deallocate(coords)
deallocate(atom_name)
deallocate(hirshfeld_volume)
deallocate(coupled_atom_pol)
deallocate(eigenmat0)
deallocate(eigenmat1)
!
! Close files used during program execution...
!
CLOSE(999)
CLOSE(20)
CLOSE(30)
!
CONTAINS
 
subroutine get_alpha_omega_and_Rp(HF,C6_free,alpha_free,w,a_w,Rpw)
       implicit none
       real*8,intent(in)::HF,C6_free,alpha_free,w
       real*8,intent(out)::a_w,Rpw
       real*8::eff_freq,w_eff_freq_ratio
       ! alpha(iomega)
       eff_freq=0.0d0
       w_eff_freq_ratio=0.0d0
       eff_freq= ((4.d0/3.d0)*C6_free/(alpha_free**2.0))**2 !eta**2
       w_eff_freq_ratio=(w*w)/eff_freq
       a_w=(HF*alpha_free)/(1.0+w_eff_freq_ratio)

       !Rp(iomega) ! Rp_eff definition from   A. Mayer, Phys. Rev. B,
       !75,045407(2007)
       Rpw= ((a_w/3.d0)*dsqrt(2.0d0/pi))**0.333333333333333333333333333d0 !sigma

       return
endsubroutine get_alpha_omega_and_Rp
!
subroutine gauss_legendre_grid()
         casimir_omega(0) = 0.0000000 ;  casimir_omega_weight(0) = 0.0000000
         casimir_omega(1) = 0.0392901 ;  casimir_omega_weight(1) = 0.0786611
         casimir_omega(2) = 0.1183580 ;  casimir_omega_weight(2) = 0.0796400
         casimir_omega(3) = 0.1989120 ;  casimir_omega_weight(3) = 0.0816475
         casimir_omega(4) = 0.2820290 ;  casimir_omega_weight(4) = 0.0847872
         casimir_omega(5) = 0.3689190 ;  casimir_omega_weight(5) = 0.0892294
         casimir_omega(6) = 0.4610060 ;  casimir_omega_weight(6) = 0.0952317
         casimir_omega(7) = 0.5600270 ;  casimir_omega_weight(7) = 0.1031720
         casimir_omega(8) = 0.6681790 ;  casimir_omega_weight(8) = 0.1136050
         casimir_omega(9) = 0.7883360 ;  casimir_omega_weight(9) = 0.1273500
        casimir_omega(10) = 0.9243900 ; casimir_omega_weight(10) = 0.1456520
        casimir_omega(11) = 1.0817900 ; casimir_omega_weight(11) = 0.1704530
        casimir_omega(12) = 1.2684900 ; casimir_omega_weight(12) = 0.2049170
        casimir_omega(13) = 1.4966100 ; casimir_omega_weight(13) = 0.2544560
        casimir_omega(14) = 1.7856300 ; casimir_omega_weight(14) = 0.3289620
        casimir_omega(15) = 2.1691700 ; casimir_omega_weight(15) = 0.4480920
        casimir_omega(16) = 2.7106200 ; casimir_omega_weight(16) = 0.6556060
        casimir_omega(17) = 3.5457300 ; casimir_omega_weight(17) = 1.0659600
        casimir_omega(18) = 5.0273400 ; casimir_omega_weight(18) = 2.0635700
        casimir_omega(19) = 8.4489600 ; casimir_omega_weight(19) = 5.6851000
        casimir_omega(20) = 25.451700 ; casimir_omega_weight(20) = 50.955800
return
endsubroutine gauss_legendre_grid

subroutine bare_dipole(dxyz,r_ij,Rvdw12,beta,Tij)
      real*8,dimension(3),intent(in) :: dxyz
      real*8,intent(in) :: r_ij
      real*8,intent(in) :: Rvdw12
      real*8,intent(in) :: beta
      real*8,dimension(3,3),intent(out) :: Tij
      real*8 :: erf
      real*8 :: d_param,fermi_param
      ! local vars
      integer :: i_index, j_index
             do i_index=1,3,1
               do j_index=1,3,1
                Tij(i_index,j_index)=3.d0*dxyz(i_index)*dxyz(j_index)
               enddo
             enddo

             do i_index=1,3,1
               Tij(i_index,i_index) = Tij(i_index,i_index) - r_ij**2
             enddo
             Tij=Tij/r_ij**5
             Tij=-1.d0*Tij
      return
endsubroutine bare_dipole 


subroutine calculate_scs_matrix()
    integer ::i_index,j_index
    real*8,dimension(3,3)::alpha_0
    real*8,dimension(3,3)::TPP
    real*8,dimension(3) :: dxyz
    real*8,dimension(3) :: coord_curr
    real*8 :: r_ij
    real*8 :: r_pp
    real*8 :: Rvdw12
    real*8 :: beta
    integer :: i_row, i_col
    integer :: errorflag

    !For LAPACK
    integer,dimension(3*n_atoms):: IPIV
    real*8,dimension(3*n_atoms):: WORK

    ! initio values
    relay_matrix=0.0d0

    beta=0.83

    !bare alpha matrix:

    ! compute relay matrix of  cluster or unit cell
    do i_row=1,n_atoms,1 !#1
      !
      do i_col=i_row,n_atoms,1 !#2
        !
        TPP=0.d0
        !
        if(i_row.eq.i_col) then  !$1
          !
          do i_index=1,3,1
            !
            do j_index=1,3,1
              !
              if(i_index.eq.j_index) then
                !
                relay_matrix(3*i_row-3+i_index,3*i_col-3+j_index)=1.d0/alpha_omega(i_row)
                !
              else
                !
                relay_matrix(3*i_row-3+i_index,3*i_col-3+j_index)=0.d0
                !
              endif !i_index=j_index
              !
            enddo!j_index
            !
          enddo!i_index
          !
        else!i_row=i_col
          !
          dxyz(:) = coords(:,i_col)-coords(:,i_row)
          r_ij = dsqrt((dxyz(1))**2.0d0 +(dxyz(2))**2.0d0 + (dxyz(3))**2.0d0 )
          r_pp = dSqrt(R_p(i_row)**2 + R_p(i_col)**2)
          Rvdw12 = Rvdw_iso(i_row) +  Rvdw_iso(i_col)
          call bare_dipole(dxyz,r_ij,Rvdw12,beta,TPP)
          do i_index=1,3,1
            !
            do j_index=1,3,1
              !
              relay_matrix(3*i_row-3+i_index,3*i_col-3+j_index)=TPP(i_index,j_index)
              relay_matrix(3*i_col-3+j_index,3*i_row-3+i_index)=TPP(i_index,j_index)
              !
            enddo
            !
          enddo
          !
        endif !i_row=i_col
        !
      enddo !#2
      !
    enddo  !#1
    call DGETRF(3*n_atoms, 3*n_atoms, relay_matrix, 3*n_atoms, IPIV, errorflag)
    if(errorflag.ne.0) then
      print '(4A)', "Error** Matrix inversion failed"
    endif

    call DGETRI(3*n_atoms, relay_matrix, 3*n_atoms, IPIV, WORK,3*n_atoms,errorflag )
    if(errorflag.ne.0) then
      print '(4A)', "Error** Matrix inversion failed"
    endif
  return
endsubroutine  calculate_scs_matrix

subroutine contract_matrix(tensor)
  implicit none

  integer::ir,ic,i,j,i_row,i_col
  real*8,dimension(3,3)::tensor
  tensor(:,:)=0.0d0

  do ir=1,n_atoms,1
       do ic=1,n_atoms,1
           i_row=0
              do i=3*ir-2,3*ir,1
                 i_row=i_row+1
                 i_col=0
                    do j=3*ic-2,3*ic,1
                    i_col=i_col+1
           tensor(i_row,i_col)=tensor(i_row,i_col) + relay_matrix(i,j)
                    enddo
              enddo
       enddo
  enddo

  return
endsubroutine contract_matrix

subroutine calculate_screened_polarizability(iso_polar_coupled)
  implicit none

  integer::i_row,i_col,i_index,j_index
  real*8,dimension(3,3)::matrix
  real*8,dimension(n_atoms),intent(OUT)::iso_polar_coupled
  real*8 :: WORK(9),eigen(3)
  integer ::   LWORK,errorflag
  iso_polar_coupled=0.d0

  ! o polar_coupled is formed by summing either rows 'blocks' or column
  ! 'blocks' of relay_matrix_screened and contains full anisotropic atom
  ! resolved polarizability 

  do i_row=1,n_atoms,1
  matrix=0.0
    do i_col=1,n_atoms,1
      do i_index=1,3,1
        do j_index=1,3,1
          matrix(i_index,j_index) = matrix(i_index,j_index) + relay_matrix(3*i_row-3+i_index,3*i_col-3+j_index)
        enddo
      enddo
    enddo
 !o iso_polar_coupled contains average of trances of atom resolved
 !polarizabilities in polar_coupled
  LWORK=9
  call DSYEV('V','U',3,matrix,3,eigen,WORK,LWORK,errorflag)
  iso_polar_coupled(i_row) =(sum(eigen))/3.d0
  enddo
  return
endsubroutine calculate_screened_polarizability
!
! Free-atom C6, polarizability and vdW radius
subroutine get_vdw_param(atom, C6, alpha, R0)
  implicit none
  
  ! local variables
  character*2 :: atom
  real*8 :: C6
  real*8 :: alpha
  real*8 :: R0
  
  select case (atom)
  
  case('H')
  alpha=4.500000
  C6=6.500000
  R0=3.100000
  
  case('C')
  alpha=12.000000
  C6=46.600000
  R0=3.590000
  
  case('G')
  alpha=6.000000*NcoverNg
  C6=11.6500000*NcoverNg*NcoverNg
  R0=2.8500000*NcoverNg !YY: is this ok for ring AND disk?
  
  case default
     C6=0.0 
     alpha=1.0
     R0=1.0
     print '(1X,4A)', '*** WARNING: VdW parameters not defined for atom: ', atom
     stop
  
  end select

end subroutine get_vdw_param


END PROGRAM RPA 
