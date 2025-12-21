! This test checks that the error handler does not crash the program when the status is successful


program test_error_treat
  use error_treat
  implicit none
  integer :: stat

  write(*,*) "=== TEST MODULE: ERROR_TREAT ==="

  ! --- TEST 1: error_allocate (Success Case) ---
  stat = 0
  call error_allocate(stat)
  write(*,*) "[OK] error_allocate: Program continued successfully with stat=0."

  ! Note for the user
  write(*,*) "NOTE: To test actual failure, change stat=1 in the code."
  write(*,*) "      The program should stop and write 'Memory allocation failed'."

end program test_error_treat
