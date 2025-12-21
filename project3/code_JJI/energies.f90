! This module contains 3 functions that calculate energies

module energies

  implicit none

contains

  ! 1 - Calculate total kinetic energy T

  double precision function T(Natoms, velocity, mass)
    implicit none
    integer, intent(in) :: Natoms
    double precision, intent(in) :: velocity(Natoms,3)
    double precision, intent(in) :: mass(Natoms)
    integer :: i
    double precision :: v2

    T = 0.0d0 ! initialize variable

    ! Loop over all atoms
    do i = 1,Natoms
      ! First compute magnitude of velocity
      v2 = velocity(i,1)**2 + velocity(i,2)**2 + velocity(i,3)**2
      ! Then compute kinetic energy T
      T = T + mass(i) * v2
    end do

    ! Finally multiply by the factor 1/2
    T = 0.5d0*T

  end function T


  ! 2 - Calculate total potential energy V using Lennard-Jones potential

  double precision function V(epsilon, sigma, Natoms, distance)
    implicit none
    double precision, intent(in) :: epsilon, sigma
    integer, intent(in) :: Natoms
    double precision, intent(in) :: distance(Natoms,Natoms)
    double precision :: rij
    integer :: i, j
    
    V = 0.0d0 ! initialize variable

    ! same nested loop as in compute_distances
    ! without a repetition of pairs (only upper diagonal part of the matrix is filled and used)
    do i = 1,Natoms-1
      do j = i+1,Natoms
        ! Compute total potential energy using Lennard Jones potential
        rij = distance(i,j)
        V = V + 4.0d0*epsilon * ( (sigma/rij)**12 - (sigma/rij)**6 )
      end do
    end do

  end function V


  ! 3 - Calculate Total energy E = T + V

  double precision function E(V_total,T_total)
    double precision, intent(in) :: V_total,T_total

    E = T_total + V_total

  end function E

end module energies
