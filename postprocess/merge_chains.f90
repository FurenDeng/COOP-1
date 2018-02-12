program getdist
  use coop_wrapper_utils
  use coop_inifile_mod
  implicit none
#include "constants.h"
  COOP_STRING root1, root2, tmpstr
  COOP_LONG_STRING:: line
  type(coop_file)::fp1, fp2
  COOP_INT::row1, col1, row2, col2, ind, i,perc, j
  root1 = coop_InputArgs(1)
  if(trim(root1).eq.'')stop "./MergeChains chainname1 chainname2 [num_chains] [discard_percent]"
  root2 = coop_InputArgs(2)
  if(trim(root2).eq.'' .or. trim(root2).eq.trim(root1))stop "./MergeChains chainname1 chainname2 [numchains]"  
  tmpstr = coop_InputArgs(3)
  if(trim(tmpstr).eq.'')then
     ind = 8
  else
     read(tmpstr, *) ind
  endif
  tmpstr = coop_InputArgs(4)
  if(trim(tmpstr).eq.'')then
     perc = 30
  else
     read(tmpstr, *) perc
  endif
  
  do i=1, ind
     row2 = coop_file_numlines(trim(root2)//"_"//COOP_STR_OF(i)//".txt")
     call fp1%open(trim(root1)//"_"//COOP_STR_OF(i)//".txt", "a")
     call fp2%open(trim(root2)//"_"//COOP_STR_OF(i)//".txt")
     do j=1, row2*perc/100
        read(fp2%unit, "(a)")line
     enddo
     do j= row2*perc/100+1, row2
        read(fp2%unit, "(a)")line
        write(fp1%unit, "(a)") trim(line)
     enddo
     call fp1%close()
     call fp2%close()
  enddo
end program getdist


