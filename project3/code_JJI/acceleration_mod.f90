! This module contains 1 subroutine that computes acceleration

module acceleration_mod

  implicit none

contains

  subroutine compute_acc(Natoms, coord, mass, distance, epsilon, sigma, acceleration)
    implicit none
    integer, intent(in) :: Natoms
    double precision, intent(in) :: coord(Natoms,3)
    double precision, intent(in) :: mass(Natoms)
    double precision, intent(in) :: distance(Natoms,Natoms)
    double precision, intent(in) :: epsilon, sigma
    double precision, intent(out) :: acceleration(Natoms,3)

    double precision :: rij, U_r, F_x, F_y, F_z, force(Natoms,3)
    integer :: i, j
    
    acceleration = 0.0d0 ! initialize variables
    force = 0.0d0

    ! 1 - Compute all different forces (upper triangular part)
    do i = 1,Natoms-1
      do j = i+1,Natoms
        
        ! IDEA: añadir un if con treshold para rij (evitar dividir entre 0)
        
        rij = distance(i,j)

        ! Compute U(r_ij)
        U_r = 24.0d0*(epsilon/rij) * ( (sigma/rij)**6 - 2.0d0*(sigma/rij)**12 )

        ! Compute force components F_ij for x, y, z
        F_x = U_r * (coord(i,1)-coord(j,1)) / rij
        F_y = U_r * (coord(i,2)-coord(j,2)) / rij
        F_x = U_r * (coord(i,3)-coord(j,3)) / rij

        ! Finally compute total forces for each atom F_i, for each coordinate x, y, z
        ! Add all components contributing to each total force
        force(i,1) = force(i,1) + F_x
        force(i,2) = force(i,2) + F_y
        force(i,3) = force(i,3) + F_z

        ! Calculate the rest of force components:
        ! Components equivalence: e.g. Force_x(21) = - Force_x(12)
        force(j,1) = force(j,1) - F_x
        force(j,2) = force(j,2) - F_y
        force(j,3) = force(j,3) - F_z

      end do
    end do

    ! 2 - Compute acceleration for each atom
      
    do i = 1,Natoms
      acceleration(i,1) = - force(i,1) / mass(i)
      acceleration(i,2) = - force(i,2) / mass(i)
      acceleration(i,3) = - force(i,3) / mass(i)
    end do

  end subroutine compute_acc

end module acceleration_mod


