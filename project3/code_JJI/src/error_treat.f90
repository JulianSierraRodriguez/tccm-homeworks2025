! This module contains 1 subroutine that controls if the allocation failed or not

module error_treat

  implicit none

contains

  ! 1 - Checks that the array was allocated correctly
  subroutine error_allocate(i_stat)
    implicit none
    integer, intent(in) :: i_stat

    ! Requirement: i_stat has to be 0 or the allocation failed.
    if (i_stat /= 0) then
      write(10,'(A)') "Memory allocation failed!"
      ! Stops the program if allocation failed.
      stop 1
    end if

  end subroutine error_allocate

end module error_treat

