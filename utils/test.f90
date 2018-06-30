program Test
#include "constants.h"    
  use coop_wrapper_utils
  implicit none
  COOP_REAL,dimension(:,:),allocatable::mat
  COOP_INT::n
  COOP_UNKNOWN_STRING,parameter::filename = "test.covmat"
  type(coop_file)::fp
  logical::success
  n = coop_file_numColumns(filename)
  allocate(mat(n,n))
  call fp%open(filename,"r")
  call coop_read_matrix(fp%unit, mat, n, n, success)
  if(.not. success) stop "broken covmat"
  call fp%close()
  call coop_matrix_inverse(mat)
  print*, mat(:,1)
  print*, mat(:,2)
  
end program Test



