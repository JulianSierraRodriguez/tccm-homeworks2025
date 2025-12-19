program main
  use read_input ! calls the module to be able to use its functions and subroutines.
  use potential_energy
  use kinetic_energy
  implicit none

  integer :: Natoms, i_stat, steps, M, i, j
  double precision :: epsilon, sigma, V_total, T_total, E_total, time_step
  double precision, external :: E
  character(100) :: input_file ! these are more than enough characters for this input.
  double precision, allocatable :: coord(:,:), mass(:), distance(:,:), velocity(:,:)

  input_file = "inp.txt"
  epsilon   = 0.997 ! kJ/mol
  sigma     = 3.405 ! angstrom
  sigma = sigma / 10 ! nm
  time_step = 0.002
  steps     = 1000
  M         = 10

  ! We call the read_Natoms function from the read_input module to obtain the number of atoms from the file.
  Natoms = read_Natoms(input_file)

  write(*,*) "Natoms = ", Natoms

  ! We allocate the arrays for the initial coordinates and the masses, and check if they were correctly allocated.
  allocate (coord(Natoms,3), stat=i_stat)
  call error_allocate(i_stat)

  allocate (mass(Natoms),stat=i_stat)
  call error_allocate(i_stat)

  ! We call the subroutine read_molecule, from the read_input module to obtain the coordinates and masses.
  call read_molecule(input_file, Natoms, coord,mass)

  write(*,*) "mass = ", mass  
  write(*,*) "coord = "
  do i = 1,Natoms
    write(*,*) coord(i,:)
  end do  

  ! We allocate the array for the distances and check if it was correctly allocated
  allocate (distance(Natoms,Natoms),stat=i_stat)
  call error_allocate(i_stat)

  ! We call the subroutine compute_distances, from the read_input module to compute the distances between the atoms.
  call compute_distances(Natoms, coord, distance)

  write(*,*) "distance = "
  do i = 1,Natoms
    write(*,*) distance(i,:)
  end do

  !!!!!! Second module potential_energy

  V_total = V(epsilon, sigma, Natoms, distance) ! kJ/mol
  write(*,*) "V_total = ", V_total, "kJ/mol"

  ! allocate the velocities array, check it, and initialize it to 0.
  allocate (velocity(Natoms,Natoms),stat=i_stat)
  call error_allocate(i_stat)

  do i = 1,Natoms
    do j = 1,Natoms
      velocity(i,j) = 0.0d0
    end do
  end do

  write(*,*) "velocity = "
  do i = 1,Natoms
    write(*,*) velocity(i,:)
  end do

  T_total = T(Natoms, velocity, mass)

  write(*,*) "T_total = ", T_total, "kJ/mol"

  E_total = E(V_total,T_total)
  write(*,*) "E_total = ", E_total, "kJ/mol"


end program main

subroutine error_allocate(i_stat)
  ! This subroutine checks if the allocation failed or not.
  implicit none
  integer, intent(in) :: i_stat
  if (i_stat /= 0) then
    print *, "Memory allocation failed!"
  end if
end subroutine error_allocate

double precision function E(V_total,T_total)
  double precision, intent(in) :: V_total,T_total

  E = T_total + V_total
  
end function E