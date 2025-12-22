! This program tests the verlet_al module, where we have our verlet algorithm procedure.

program test_verlet

  ! We call the needed modules.
  use verlet_al
  use treat_input
  use treat_output
  use acceleration_mod

  implicit none
  integer, parameter :: Natoms = 2
  double precision :: coord(Natoms,3), velocity(Natoms,3), acc(Natoms,3), mass(Natoms), dist(Natoms,Natoms)
  double precision :: dt = 0.1d0, epsilon = 1.0d0, sigma = 1.0d0
  double precision :: E, V_p, T
  integer :: u_out = 10, u_xyz = 20, iter = 1, M = 1
  character(len=100) :: name_out, name_xyz
  logical :: velo_test = .false., acc_test = .false., coord_test = .false.
  
  write(*,*) "=== TEST MODULE: VERLET ALGORITHM ==="

  ! We prepare some dummy files so we can test verlet and write some results.

  name_out = 'dummy_verlet.out'; name_xyz = 'dummy_verlet.xyz'
  call create_output(u_out, name_out)
  call create_output(u_xyz, name_xyz)

  ! We give initial values to some dummy arrays. We give the first atom x = -0.5 and the second x = 0.5.
  coord = 0.0d0
  coord(1,1) = -0.5d0
  coord(2,1) =  0.5d0
  velocity   = 0.0d0
  mass  = 1.0d0
  
  ! We give the acceleration and distance arrays zeros to all values.
  acc = 0.0d0
  dist = 0.0d0
  ! We compute the initial distance, needed to compute the initial acceleration.
  call compute_distances(Natoms, coord, dist)
  
  ! We compute the initial acceleration, it is needed to start the Verlet algorithm.
  call compute_acc(Natoms, coord, mass, dist, epsilon, sigma, acc)

  ! We check that we obtained the correct values, this is not yet part of the verlet algorithm.
  ! The correct values are a_x1 = -24.0 and a_x2 = 24.0, the rest are all 0.0.
  if ((abs(acc(1,1) + 24.0d0) < 1e-6) .and. (abs(acc(2,1) - 24.0d0) < 1e-6)) then
     write(*,*) "[OK] Acceleration mod: Initial acceleration is correct."
  else
     write(*,*) "[FAIL] Acceleration mod: Incorrect acceleration.", "Value acc_x1 :", acc(1,1), "Expected: -24.0"
     write(*,*) "                                                ", "Value acc_x2 :", acc(2,1), "Expected:  24.0"
  end if

  ! Call Verlet for 1 step, for our two-atom system.
  call verlet_algorithm(Natoms, coord, velocity, acc, dt, epsilon, sigma, mass, T, V_p, E, u_xyz, u_out, iter, M)


  ! TEST 1: Check  coordinates. The correct result is x1 = -0.62 and x2 = 0.62.
  ! If it is correct we update the logical for the coordinates to true.

  if ((abs(coord(1,1) + 0.62d0) < 1e-6) .and. (abs(coord(2,1) - 0.62d0) < 1e-6)) then
    write(*,*) "[OK] Verlet: Updated positions are correct."
    coord_test = .true.
  else
    write(*,*) "[FAIL] Verlet: Incorrect position.", "Value x1 :", coord(1,1), "Expected: -0.62"
    write(*,*) "                                  ", "Value x2 :", coord(2,1), "Expected:  0.62"
  end if


  ! TEST 2: Check accelerations. The correct result is a_x1 = 1.931445... and a_x2 = -1.931445....
  ! If it is correct we update the logical for the acceleration to true.

  if ((abs(acc(1,1) - 1.931445d0) < 1e-6) .and. (abs(acc(2,1) + 1.931445d0) < 1e-6)) then
    write(*,*) "[OK] Verlet: Updated accelerations are correct."
    acc_test = .true.
  else
    write(*,*) "[FAIL] Verlet: Incorrect acceleration.", "Value a_x1 :", acc(1,1), "Expected:  1.931445..."
    write(*,*) "                                      ", "Value a_x2 :", acc(2,1), "Expected: -1.931445..."
  end if


  ! TEST 3: Check velocities. The correct result is v_x1 = -1.103428... and v_x2 = 1.103428....
  ! If it is correct we update the logical for the velocities to true.

  if ((abs(velocity(1,1) + 1.103428d0) < 1e-6) .and. (abs(velocity(2,1) - 1.103428d0) < 1e-6)) then
    write(*,*) "[OK] Verlet: Updated velocities are correct."
    velo_test = .true.
  else
    write(*,*) "[FAIL] Verlet: Incorrect velocities.", "Value v_x1 :", velocity(1,1), "Expected: -1.103428..."
    write(*,*) "                                    ", "Value v_x2 :", velocity(2,1), "Expected:  1.103428..."
  end if

  write(*,*)

  ! RESUME TESTS: If all the tests passed we obtain the pass.
  if (coord_test .and. velo_test .and. acc_test) then
    write(*,*) "[OK] Verlet: The verlet algorithm works correctly."
  else
    write(*,*) "[FAIL] Verlet: One or more tests did not pass. Check above to identify the component or components that failed."
  end if

  write(*,*) ! Separation line

  ! Cleanup: Close and delete dummy files
  close(u_out, status='delete')
  close(u_xyz, status='delete')

end program test_verlet
