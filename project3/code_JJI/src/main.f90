program main
  use treat_input ! calls the module to be able to use its functions and subroutines.
  use treat_output
  use energies
  use acceleration_mod
  use error_treat
  use verlet_al

  implicit none

  integer :: Natoms, i_stat, steps, M, i, j, iter
  integer :: u_output, u_xyz
  double precision :: epsilon, sigma, V_total, T_total, E_total, time_step
  ! double precision, external :: E
  character(100) :: input_file, output_file, xyz_file ! these are more than enough characters for these namefiles.
  double precision, allocatable :: coord(:,:), mass(:), distance(:,:), velocity(:,:), acceleration(:,:)


  !------------------------------!
  ! PREPARING FOR THE SIMULATION !
  !------------------------------!


  ! Simulation parameters 

  input_file = "inp.txt"
  epsilon   = 0.997 ! kJ/mol
  sigma     = 3.405 ! angstrom
  sigma = sigma / 10 ! nm
  time_step = 0.002; steps = 1000;  M = 10

  ! Outputs parameters

  u_output = 10; u_xyz = 20
  output_file = "MD.out"; xyz_file = "trajectory.xyz"

  ! We call the create output (treat_output module) to generate or replace the output and trajectory files.
  call create_output(u_output,output_file)
  call create_output(u_xyz,xyz_file)

  ! We call the read_Natoms function (treat_input module) to obtain the number of atoms from the file.
  Natoms = read_Natoms(input_file)

  ! We call the write_output_head function (treat_output module) to output the initial parameters to the output file.
  call write_output_head(u_output, input_file, epsilon, sigma, time_step, steps, M, output_file, u_xyz, xyz_file, Natoms)

  ! We allocate the arrays for the initial coordinates and the masses, and check if they were correctly allocated using the error_allocate subroutine (error_treat module).
  allocate (coord(Natoms,3), stat=i_stat)
  call error_allocate(i_stat)

  allocate (mass(Natoms), stat=i_stat)
  call error_allocate(i_stat)

  ! We call the subroutine read_molecule (treat_input module) to obtain the coordinates and masses from the input file.
  call read_molecule(input_file, Natoms, coord, mass)

  ! We allocate the array for the distances and check if it was correctly allocated.
  allocate (distance(Natoms,Natoms), stat=i_stat)
  call error_allocate(i_stat)

  ! We call the subroutine compute_distances (treat_input module) to compute the distances between the atoms.
  call compute_distances(Natoms, coord, distance)


  !----------------------!
  ! INITIAL ENERGY TERMS !
  !----------------------!


  ! We use the V function (energies module) to obtain the initial potential energy
  V_total = V(epsilon, sigma, Natoms, distance) ! kJ/mol

  ! We allocate the velocities array, check it, and then we initialize it to 0.
  allocate (velocity(Natoms,3), stat=i_stat)
  call error_allocate(i_stat)

  velocity = 0.0d0

  ! We use the T function (energies module) to obtain the initial kinetic energy.
  ! It should be 0, we don't have any velocities.
  T_total = T(Natoms, velocity, mass) ! kJ/mol

  ! We use the E function (energies module) to obtain the initial total energy.
  E_total = E(V_total, T_total) ! kJ/mol


  !----------------------------!
  ! Step "0" of the simulation !
  !----------------------------!


  ! We allocate the acceleration array, and check it.
  allocate (acceleration(Natoms,3),stat=i_stat)
  call error_allocate(i_stat)

  ! We call compute_acc (acceleration_mod module) to obtain the initial acceleration from the forces from the initial configuration of atoms.
  call compute_acc(Natoms, coord, mass, distance, epsilon, sigma, acceleration)

  ! We write the initial step into the trajectory file.
  call write_xyz(u_xyz, Natoms, coord, E_total, V_total, T_total)

  
  !------------!
  ! Simulation !
  !------------!


  write(u_output,'(A)') "|---------------------------|"
  write(u_output,'(A)') "|     Starting the MD       |"
  write(u_output,'(A)') "|---------------------------|"

  ! We do an iterative loop for all the steps of the simulation, where we call the verlet_algorithm subroutine (verlet_al module), that does the verlet algorithm, every M steps it prints the coordinates into the trajectory file, with the correct format, and every 10*M it
  ! uploads the energy to the output file, this way it gives the user an overview of the simulation process.
  
  do iter = 1,steps
    call verlet_algorithm(Natoms,coord,velocity,acceleration, time_step,epsilon,sigma,mass,T_total, V_total, &
    E_total,u_xyz,u_output,iter,M)
  end do

  write(u_output,'(A)') "|---------------------------|"
  write(u_output,'(A)') "| The MD ended succesfully! |"
  write(u_output,'(A)') "|---------------------------|"


  ! We close the files.
  close(u_output)
  close(u_xyz)

  ! We deallocate the arrays.
  deallocate(mass)
  deallocate(coord)
  deallocate(velocity)
  deallocate(acceleration)
  deallocate(distance)
end program main



