program main
  use read_input ! calls the module to be able to use its functions and subroutines.
  use potential_energy
  implicit none

  integer :: Natoms, i_stat, i
  double precision :: epsilon, sigma
  character(100) :: input_file ! these are more than enough characters for this input.
  double precision, allocatable :: coord(:,:), mass(:), distance(:,:), V(:,:)

  ! We call the read_Natoms function from the read_input module to obtain the number of atoms from the file.
  input_file = "inp.txt"
  epsilon = 0.997
  sigma   = 3.405

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

  ! We allocate the array for the potential energy and check if it was correctly allocated
  allocate (V(Natoms,Natoms),stat=i_stat)
  call error_allocate(i_stat)

  
  V = VLJ(epsilon, sigma, Natoms, distance)
  write(*,*) "V = "
  do i = 1,Natoms
    write(*,*) V(i,:)
  end do

end program main

subroutine error_allocate(i_stat)
  ! This subroutine checks if the allocation failed or not.
  implicit none
  integer, intent(in) :: i_stat
  if (i_stat /= 0) then
    print *, "Memory allocation failed!"
  end if
end subroutine error_allocate