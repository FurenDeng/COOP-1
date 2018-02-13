program Test
#include "constants.h"    
  use coop_wrapper_utils
  implicit none

  COOP_INT::i
  COOP_REAL::homework(90), midterm(90)
  type(coop_file)::fp
  type(coop_asy)::fig
  call fp%open("tmp.txt", "r")
  call fig%open("scorestats.txt")
  call fig%init(ylabel="homework", xlabel="midterm")
  do i=1, 90
     read(fp%unit, *) midterm(i), homework(i)
     print*, midterm(i), homework(i)
  enddo
  call fig%dots(midterm, homework)
  call fig%close()
  call fp%close()
end program Test



