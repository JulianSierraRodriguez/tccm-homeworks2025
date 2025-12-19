! This module contains 1 function that calculates total potential energy using Lennard-Jones potential

module potential_energy

  implicit none

contains

  double precision function VLJ(epsilon, sigma, Natoms, distance)
    implicit none
    double precision, intent(in) :: epsilon, sigma
    integer, intent(in) :: Natoms
    double precision, intent(in) :: distance(Natoms,Natoms)
    double precision :: rij
    integer :: i, j
    
    VLJ = 0.0d0 ! initialize variable

    ! same nested loop as in compute_distances
    ! without a repetition of pairs (only upper diagonal part of the matrix is filled and used)
    do i = 1,Natoms-1
      do j = i+1,Natoms
        ! Compute total potential energy using Lennard Jones potential
        rij = distance(i,j)
        VLJ = VLJ + 4.0d0*epsilon * ( (sigma/rij)**12 - (sigma/rij)**6 )
      end do
    end do

  end function VLJ

end module potential_energy
