program main
  use read_input ! calls the module to be able to use its functions and subroutines.
  use potential_energy
  use kinetic_energy
  use acceleration_mod

  implicit none

  integer :: Natoms, i_stat, steps, M, i, j
  integer :: u_output, u_xyz
  double precision :: epsilon, sigma, V_total, T_total, E_total, time_step
  double precision, external :: E
  character(100) :: input_file, output_file, xyz_file ! these are more than enough characters for this input.
  double precision, allocatable :: coord(:,:), mass(:), distance(:,:), velocity(:,:), acceleration(:,:)

  !       MD parameters 

  input_file = "inp.txt"
  epsilon   = 0.997 ! kJ/mol
  sigma     = 3.405 ! angstrom
  sigma = sigma / 10 ! nm
  time_step = 0.002
  steps     = 1000
  M         = 10

  !    Outputs parameters

  u_output = 10
  u_xyz    = 20
  output_file = "MD.out"
  xyz_file    = "trajectory.xyz"

  call create_output(u_output,output_file)
  call create_output(u_xyz,xyz_file)

  ! We call the read_Natoms function from the read_input module to obtain the number of atoms from the file.
  Natoms = read_Natoms(input_file)

  call write_output_head(u_output, input_file, epsilon, sigma, time_step, steps, M, output_file, u_xyz, xyz_file,Natoms)

  ! write(*,*) "Natoms = ", Natoms

  ! We allocate the arrays for the initial coordinates and the masses, and check if they were correctly allocated.
  allocate (coord(Natoms,3), stat=i_stat)
  call error_allocate(i_stat)

  allocate (mass(Natoms),stat=i_stat)
  call error_allocate(i_stat)

  ! We call the subroutine read_molecule, from the read_input module to obtain the coordinates and masses.
  call read_molecule(input_file, Natoms, coord,mass)

  ! write(*,*) "mass = ", mass  
  ! write(*,*) "coord = "
  ! do i = 1,Natoms
  !   write(*,*) coord(i,:)
  ! end do  

  ! We allocate the array for the distances and check if it was correctly allocated
  allocate (distance(Natoms,Natoms),stat=i_stat)
  call error_allocate(i_stat)

  ! We call the subroutine compute_distances, from the read_input module to compute the distances between the atoms.
  call compute_distances(Natoms, coord, distance)

  ! write(*,*) "distance = "
  ! do i = 1,Natoms
  !   write(*,*) distance(i,:)
  ! end do

  !!!!!! Second module -> potential_energy

  ! We use the V function to obtain the initial potential energy
  V_total = V(epsilon, sigma, Natoms, distance) ! kJ/mol
  ! write(*,*) "V_total = ", V_total, "kJ/mol"

  !!!!!! Third module -> kinetic_energy

  ! allocate the velocities array, check it, and initialize it to 0.
  allocate (velocity(Natoms,3),stat=i_stat)
  call error_allocate(i_stat)

  do i = 1,Natoms
    do j = 1,Natoms
      velocity(i,j) = 0.0d0
    end do
  end do

  ! write(*,*) "velocity = "
  ! do i = 1,Natoms
  !   write(*,*) velocity(i,:)
  ! end do

  ! We use the T function to obtain the initial kinetic energy
  T_total = T(Natoms, velocity, mass)

  ! write(*,*) "T_total = ", T_total, "kJ/mol"

  ! E_total = E(V_total,T_total)
  ! write(*,*) "E_total = ", E_total, "kJ/mol"

  !!!!!! Fourth module -> acceleration

  allocate (acceleration(Natoms,3),stat=i_stat)
  call error_allocate(i_stat)

  call compute_acc(Natoms, coord, mass, distance, epsilon, sigma, acceleration)

  ! write(*,*) "acceleration = "
  ! do i = 1,Natoms
  !   write(*,*) acceleration(i,:)
  ! end do

  call write_xyz(u_xyz, Natoms, coord, E_total, V_total, T_total)

  close(u_output)
  close(u_xyz)
end program main

subroutine error_allocate(i_stat)
  ! This subroutine checks if the allocation failed or not.
  implicit none
  integer, intent(in) :: i_stat
  if (i_stat /= 0) then
    write(10,'(A)') "Memory allocation failed!"
    stop 1
  end if
end subroutine error_allocate

double precision function E(V_total,T_total)
  double precision, intent(in) :: V_total,T_total

  E = T_total + V_total

end function E

subroutine create_output(u,name_file)
  ! This subroutine creates or replaces (if it already exists) the output.
  implicit none
  integer, intent(in) :: u
  integer :: ios
  character(100), intent(in) :: name_file
  logical :: exists

  inquire(file=name_file, exist=exists)
  if (exists) then
    open(unit=u, file=name_file, status='REPLACE', action='WRITE', iostat=ios)
  else
    open(unit=u, file=name_file, status='NEW', action='WRITE', iostat=ios)
  end if

  if (ios /= 0) then
     write(10,'(A,I0)') 'Error opening', name_file,' IOSTAT = ', ios
     stop 1
  end if

end subroutine create_output


subroutine write_output_head(u_output, input_file, epsilon, sigma, time_step, steps, M, output_file, u_xyz, xyz_file, Natoms)

integer, intent(in) :: u_output, Natoms, steps, u_xyz, M
character(100), intent(in) :: input_file, output_file, xyz_file
double precision :: epsilon, sigma, time_step

write(u_output,'(A)') "|---------------|"
write(u_output,'(A)') "|   MD OUTPUT   |"
write(u_output,'(A)') "|---------------|"
write(u_output,'()')
write(u_output,'(A)')"### MD parameters ###"
write(u_output,'()')
write(u_output,'(A,1X,A)') "Input Filename =",trim(input_file)
write(u_output,'(A,1X,I0,1X,A,1X,F5.3,1X,A,1X,I0,1X,A)') &
"Steps of the simulation =", steps, &
"| Time step =", time_step, "ps | Save every",M, "steps"
write(u_output,'(A,1X,I0,1X,A)') "There are", Natoms, "atoms"
write(u_output,'()')
write(u_output,'(A)')"### Output parameters ###"
write(u_output,'()')
write(u_output,'(A,1X,A,1X,A,1X,I0)')"Output Filename =", trim(output_file), "| u =", u_output
write(u_output,'(A,1X,A,1X,A,1X,I0)')"Trajectory Filename =", trim(xyz_file), "| u =", u_xyz
write(u_output,'()')
write(u_output,'(A)')"### Lennard-Jones parameters ###"
write(u_output,'()')
write(u_output,'(A,1X,F7.4,1X,A,1X,F6.3,1X,A)')"sigma =", sigma, "nm | epsilon =", epsilon, "kJ/mol"
write(u_output,'()')

end subroutine write_output_head

subroutine write_xyz(u_xyz, Natoms, coord, E_total, V_total, T_total)
  integer, intent(in) :: u_xyz, Natoms
  double precision, intent(in) :: coord(Natoms,3)
  double precision, intent(in) :: E_total, V_total, T_total
  integer :: i
  write(u_xyz,'(I0)') Natoms
  write(u_xyz,'(A,1X,F8.5,1X,A,1X,F8.5,1X,A,1X,F8.5,1X,A)') &
  "E_total = ", E_total,"V_total = ", V_total,"T_total = ", T_total, "(kJ/mol)"
  do i = 1,Natoms
    write(u_xyz,'(A,1X,F10.7,1X,F10.7,1X,F10.7)')'Ar',coord(i,:)
  end do
end subroutine write_xyz