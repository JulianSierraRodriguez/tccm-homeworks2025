! This module contains 1 subroutine that produces Verlet Algorithm.

module verlet_al

  implicit none

contains

  subroutine verlet_algorithm(Natoms,coord,velocity,acceleration, time_step,epsilon,sigma,mass,T_total, V_total, &
  E_total,u_xyz,u_output,iter,M)
    use treat_input
    use treat_output
    use acceleration_mod
    use energies

    implicit none
    integer, intent(in) :: Natoms, u_xyz, u_output, iter, M
    double precision, intent(in) :: time_step, epsilon, sigma
    double precision :: coord(Natoms,3)
    double precision :: velocity(Natoms,3)
    double precision :: acceleration(Natoms,3)
    double precision :: mass(Natoms)
    double precision :: T_total, V_total, E_total

    integer :: i,j
    double precision :: coord_new(Natoms,3)
    double precision :: distance_new(Natoms,Natoms)
    double precision :: velocity_new(Natoms,3)
    double precision :: acceleration_new(Natoms,3)

    ! We use equation 7 to compute the new positions.
    ! We use nested loops to run over the three coordinates of all the atoms.
    do i = 1,Natoms
      do j = 1,3
        coord_new(i,j) = coord(i,j) + velocity(i,j)*time_step + acceleration(i,j)*((time_step**2)/2)
      end do
    end do

    ! We compute the new distances (module treat_input).
    call compute_distances(Natoms, coord_new, distance_new)

    ! We compute the new accelerations (module acceleration_mod).
    call compute_acc(Natoms, coord, mass, distance_new, epsilon, sigma, acceleration_new)

    ! We use equation 8 to compute the new velocities.
    ! We use nested loops to run over the three velocities of all the atoms.
    do i = 1,Natoms
      do j = 1,3
        velocity_new(i,j) = velocity(i,j) + 0.5*(acceleration(i,j)+ acceleration_new(i,j))*time_step
      end do
    end do

    ! We compute the new set of Enery terms and the total energy (the three functions are in the energies module).
    V_total = V(epsilon, sigma, Natoms, distance_new) ! kJ/mol
    T_total = T(Natoms, velocity_new, mass) ! kJ/mol
    E_total = E(V_total,T_total) ! kJ/mol

    ! We output the new positions every M steps on the trajectory file.
    
    ! Explanation about the formating, after the comma inside write():
    ! A -> string of text
    ! 1X -> one space between what we print, we can change the number for different number of spaces.
    ! I0 -> to write an integer.
    ! FY.Z -> To write a float using a total of Y digits (including the minus sign) with Z digits after the decimal point. 

    if (mod(iter,M) .eq. 0) then
      ! We only output the coordinates every M steps in the u_xyz file (trajectory file).
      call write_xyz(u_xyz, Natoms, coord_new, E_total, V_total, T_total)

      if (mod(iter,M*10) .eq. 0) then
        ! We only output the energies every 10*M steps in the u_output file (output file), this is only to have a quick overview 
        ! to see if the simulation performed correctly.
        write(u_output, '(A,1X,I0)') "Step =",iter
        write(u_output,'(A,1X,F8.5,1X,A,1X,F8.5,1X,A,1X,F8.5,1X,A)') &
        "E_total = ", E_total,"| V_total = ", V_total,"| T_total = ", T_total, "(kJ/mol)"
        write(u_output, '(A)') "------------------------------------"
      end if
    end if

    ! We upload our set of variables with their data in the new step to continue the simulation.
    coord = coord_new
    velocity = velocity_new
    acceleration = acceleration_new

    ! The internal arrays are not allocatable, so they free the space automatically when the subroutine ends.

  end subroutine verlet_algorithm

end module verlet_al

