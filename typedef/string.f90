module coop_string_mod
  use coop_constants_mod
  use coop_basicutils_mod
  implicit none
#include "constants.h"

  private

  integer,parameter::sp = kind(1.)
  integer,parameter::dl = kind(1.d0)

  public::coop_num2str,  coop_ndigits, coop_str2int, coop_str2real, coop_str2logical, coop_substr, coop_str_replace, coop_str_numalpha, coop_str2lower, coop_str2upper, coop_case_insensitive_eq, coop_file_path_of, coop_file_name_of, coop_file_add_postfix, coop_convert_to_C_string, coop_convert_to_Fortran_String, coop_data_type, coop_string_contain_numbers, coop_numstr2goodstr, coop_num2goodstr, coop_string_strip_quotes, coop_str_numUpperAlpha, coop_str_numLowerAlpha

  Interface coop_num2str
     module procedure coop_int2str, coop_real2str, coop_logical2str, coop_double2str
  end Interface coop_num2str


contains

  function coop_num2goodstr(x, minus, point) result(str)
    COOP_REAL x
    COOP_UNKNOWN_STRING,optional::minus, point
    COOP_SHORT_STRING str, str_int, str_float, str_zero, str_sign, str_point
    COOP_REAL absx
    COOP_INT  n
    absx = abs(x)
    if(absx .ge. 1.d5)then
       write(str,*) x
       str = trim(adjustl(str))
       return
    endif
    if(x.lt.0.d0)then
       if(present(minus))then
          str_sign = trim(adjustl(minus))
       else
          str_sign = "-"
       endif
    else
       str_sign = ""
    endif
    if(abs(absx - nint(absx)) .lt. max(1.d-5, absx*1.d-5))then
       write(str_int, *) nint(absx)
       str = trim(str_sign)//trim(adjustl(str_int))
       return
    endif
    n = floor(absx)
    write(str_int, *) n
    if(present(point))then
       str_point = trim(adjustl(point))
    else
       str_point = "."
    endif
    absx = absx - n
    str_zero = ""
    do while(absx .lt. 0.09999d0)
       absx = absx*10.d0
       str_zero = trim(str_zero)//"0"
    enddo
    absx = absx*10000.d0
    n = nint(absx)
    do while(mod(n, 10).eq.0)
       n = n/10
    enddo
    write(str_float,*) n
    str = trim(str_sign)//trim(adjustl(str_int))//trim(str_point)//trim(str_zero)//trim(adjustl(str_float))
  end function coop_num2goodstr

  function coop_numstr2goodstr(str) result(str_out)
    COOP_STRING str, str_out
    COOP_INT idot, iminus
    idot = index(str, ".")
    if(idot.eq.0)then
       str_out = adjustl(trim(str))
    else
       str_out = adjustl(trim(str(1:idot-1)//"pt"//str(idot+1:len_trim(str))))
    endif
  end function coop_numstr2goodstr


  function coop_Ndigits(i, ndigits, base) result(str_ndigits)
    COOP_INT, optional::base
    COOP_INT i, ndigits, j 
    COOP_STRING Str_Ndigits, strzero
    if(present(base))then
       str_ndigits = coop_int2str(i, base)
    else
       str_ndigits = coop_int2str(i)
    endif
    strzero = ""
    do j = len_trim(str_ndigits)+1, ndigits
       strzero = trim(strzero)//"0"
    enddo
    str_Ndigits = trim(strzero)//trim(Str_Ndigits)
  end function coop_Ndigits

  function coop_int2str(i, base)
    COOP_INT,INTENT(IN)::i
    COOP_INT, optional::base
    COOP_STRING coop_int2str
    COOP_INT j, k
    character,dimension(36),parameter::dig = (/ "0", "1", "2", "3", "4", "5", "6", "7", "8", "9", &
         "A", "B", "C", "D", "E", "F", "G", "H", "I", "J", &
         "K", "L", "M", "N", "O", "P", "Q", "R", "S", "T", &
         "U", "V", "W", "X", "Y", "Z" /)
    coop_int2str = ""
    if(present(base))then
       if(base .gt. 36 .or. base .lt. 2) stop "coop_int2str: base overflow"
       j = i
       do while(j.gt. 0)
          k = j/base
          coop_int2str = dig(j-k*base+1)//trim(coop_int2str)
          j = k
       enddo
    else
       WRITE (coop_int2str,*) i
       coop_int2str= Trim(ADJUSTL(coop_int2str))
    endif
  End function coop_int2str

  function coop_double2str(x, fmt) result(str)
    real(dl) x
    COOP_STRING str
    COOP_UNKNOWN_STRING ,optional::fmt
    if(present(fmt))then
       str = coop_real2str(real(x, sp), fmt)
    else
       str = coop_real2str(real(x, sp))
    endif
  end function coop_double2str

  Function Coop_real2str(x,fmt) !!This is a smater algorithm that convert number to string
    COOP_INT, parameter::ndigits = 5
    COOP_STRING Coop_real2str, str
    real(sp) x, ax
    COOP_INT ind,  n, i
    COOP_INT ix
    COOP_UNKNOWN_STRING ,optional::fmt
    if(present(fmt))then
       if(fmt .eq. "none")then
          coop_real2str = ""
          return
       endif
       if(trim(fmt).eq. "d")then
          write(str,"(G22.15)") x
          if(scan(str, "e").eq.0 .and. scan(str, "E").eq.0)then
             str=trim(adjustl(str))
             n=scan(str,".")
             if(n>0)then
                n=max(n-1,verify(trim(str),"0",.true.))
             else
                n=Len_Trim(str)
             endif
             Coop_real2str =  Str(1:n)//"d0"
          else
             Coop_real2str = trim(adjustl(str))
             n = scan(Coop_real2str, "e")
             if(n.eq.0) n = scan(Coop_real2str,"E")
             if(n.le.1) stop "error in Coop_real2str: does not look like a valid number"
             ind = n-1
             if(scan(Coop_real2str(1:ind),".") .gt. 0)then
                do while(Coop_real2str(ind:ind).eq."0")
                   ind = ind -1
                enddo
             endif
             if(Coop_real2str(n+1:n+1) == "+")then
                Coop_real2str = Coop_real2str(1:ind)//"d"//trim(Coop_real2str(n+2:))
             else
                Coop_real2str = Coop_real2str(1:ind)//"d"//trim(Coop_real2str(n+1:))
             endif
          endif
          return
       endif
       if(trim(fmt).ne."f" .and. trim(fmt).ne."")then
          write(str,"("//trim(fmt)//")") x
          Coop_real2str = trim(adjustl(str))
          return
       endif
    endif
    ax = abs(x)
    if(ax.eq.0.)then
       Coop_real2str="0"
       return
    endif
    ind=0
    do while(ax.lt.1.d0)
       ind = ind -1
       ax = ax*10.d0
    enddo
    do while(ax.ge.10.)
       ind = ind + 1
       ax = ax/10.d0
    enddo
    if(ind .gt. ndigits .or. ind .lt. -ndigits)then
       write(Str,*) nint(ax * 10.d0**ndigits)
       Str=trim(adjustl(Str))
       n=verify(trim(str), "0", .true.)
       if(n.gt.1)then
          Coop_real2str=str(1:1)//"."//str(2:n)
       else
          Coop_real2str=Str(1:1)
       endif
       write(str,*) ind
       coop_real2str=Trim(Coop_real2str)//"e"//Trim(adjustl(Str))
    else
       ax=abs(x)
       write(Str,'(F16.'//trim(coop_int2str(ndigits-ind+1))//')') ax
       str=trim(adjustl(str))
       n=scan(str,".")
       i  = verify(trim(str),"0",.true.)
       if(n>0)then
          if(n.eq.i)then
             n = n -1
          else
             n = i
          endif
       else
          n=Len_Trim(str)
       endif
       coop_real2str = str(1:n)
    endif
    if(x.lt.0.d0) Coop_real2str='-'//trim(Coop_real2str)
  End Function coop_real2str
  


  function Coop_logical2str(x)
    logical x
    COOP_STRING coop_logical2str
    if(x)then
       Coop_logical2str="T"
    else
       Coop_logical2str="F"
    endif
  end function Coop_logical2str


  Function coop_str2int(str)
    COOP_UNKNOWN_STRING  str
    COOP_INT coop_str2int
    read(str,*)coop_str2int
  End Function coop_str2int

  function coop_str2real(str)
    COOP_UNKNOWN_STRING  str
    COOP_REAL coop_str2real
    read(str,*)coop_str2real
  end function coop_str2real

  function coop_str2logical(str)
    COOP_UNKNOWN_STRING str
    logical coop_str2logical
    read(str,*) coop_str2logical
  end function coop_str2logical

!!================= Now the really cool stuff ====================
  function coop_substr(str, pos_start, pos_end, terminator)
    COOP_LONG_STRING coop_substr
    COOP_UNKNOWN_STRING  str
    COOP_INT pos_start, pos_end, i
    character,optional::terminator
    if(present(terminator))then
       do i = pos_start, pos_end
          if(str(i:i) .eq. terminator) exit
       enddo
       coop_substr = str(pos_start:i-1)
    else
       coop_substr = str(pos_start:pos_end)
    endif
  end function coop_substr

  function coop_str_replace(str, needle, replacement) result(repl)
    COOP_UNKNOWN_STRING str, needle, replacement
    COOP_LONG_STRING repl
    COOP_INT i, istart, k, lenr, lenn, lens
    lenr = len(replacement)
    lenn = len(needle)
    lens = len(str)
    repl = ""
    k = 1
    istart = 1
    i = index(str(istart:), needle)
    do while(i.ne.0)
       repl(k:k+i-2+lenr) = str(istart:istart+i-2)//replacement
       k = k + i + lenr-1
       istart = istart + i + lenn-1
       if(istart .gt. lens)exit
       i = index(str(istart:), needle)
    enddo
    repl(k:k+lens - istart) = str(istart:lens)
  end function coop_str_replace


  function coop_str_numalpha(str) result(s)
    COOP_UNKNOWN_STRING str
    COOP_STRING s
    COOP_INT i, lens, j, k
    lens = len_trim(str)
    s = ""
    k = 1
    do i=1, lens
       j = ichar(str(i:i))
       if(j.ge.48 .and.j.le.57 .or. j.ge.65 .and. j.le.90 .or. j.ge.97 .and. j.le. 122 )then
          s(k:k) = str(i:i)
          k  = k  + 1
       end if
    enddo
  end function coop_str_numalpha


  function coop_str_numUpperalpha(str) result(s)
    COOP_UNKNOWN_STRING str
    COOP_STRING s
    COOP_INT i, lens, j, k
    lens = len_trim(str)
    s = ""
    k = 1
    do i=1, lens
       j = ichar(str(i:i))
       if(j.ge.48 .and.j.le.57 .or. j.ge.65 .and. j.le.90 .or. j.ge.97 .and. j.le. 122 )then
          if(j .ge. 97)then
             s(k:k) = char(j-32)
          else
             s(k:k) = str(i:i)
          endif
          k  = k  + 1
       end if
    enddo
  end function coop_str_numUpperalpha


  function coop_str_numLoweralpha(str) result(s)
    COOP_UNKNOWN_STRING str
    COOP_STRING s
    COOP_INT i, lens, j, k
    lens = len_trim(str)
    s = ""
    k = 1
    do i=1, lens
       j = ichar(str(i:i))
       if(j.ge.48 .and.j.le.57 .or. j.ge.65 .and. j.le.90 .or. j.ge.97 .and. j.le. 122 )then
          if(j .ge. 65 .and. j .le. 90)then
             s(k:k) = char(j+ 32)
          else
             s(k:k) = str(i:i)
          endif
          k  = k  + 1
       end if
    enddo
  end function coop_str_numLoweralpha

  

  subroutine coop_Str2Lower(Str)
    COOP_UNKNOWN_STRING  Str
    COOP_INT i, n, k
    n = len_trim(Str)
    do i=1, n
       k = ichar(str(i:i))
       if(k.ge.65 .and. k .le. 90)then
          str(i:i) = char(k+32)
       endif
    enddo
  End subroutine Coop_Str2Lower

  subroutine coop_Str2Upper(Str)
    COOP_UNKNOWN_STRING  Str
    COOP_INT n, i, k
    n = len_trim(Str)
    do i = 1, n
       k = ichar(str(i:i))
       if(k.ge.97 .and. k.le. 122)then
          str(i:i) = char(k-32)
       end if
    enddo
  End subroutine Coop_Str2Upper

  function Coop_case_insensitive_eq(Str1, Str2)
    COOP_UNKNOWN_STRING Str1,Str2
    logical Coop_case_insensitive_eq
    COOP_INT n1, n2, i, i1, i2
    n1=len_trim(Str1)
    n2=len_trim(Str2)
    if(n1 .ne. n2)then
       Coop_case_insensitive_eq=.false.
       return
    end if
    do i = 1, n1
       If(Str1(i:i) .ne. Str2(i:i))then
          i1 = ichar(str1(i:i))
          i2 = ichar(str2(i:i))
          if(i1 .lt. 65 .or. i1 .gt.122 .or. (i1 .gt. 90 .and. i1 .lt. 97) .or. &
               i2 .lt. 65 .or. i2 .gt.122 .or. (i2 .gt. 90 .and. i2 .lt. 97))then
             
             Coop_case_insensitive_eq=.false.
             return
          else
             if(i1 .ne. i2 .and. abs(i1-i2) .ne. 32)then
                coop_case_insensitive_eq = .false.
                return
             endif
          endif
       endif
    enddo
    Coop_case_insensitive_eq=.true.
  End function Coop_case_insensitive_eq


  function coop_file_path_of(fstr) result(path)
    COOP_UNKNOWN_STRING  fstr
    COOP_LONG_STRING path
    COOP_INT n
    n = scan(fstr, '\/', .true.)
    if(n.eq.0)then
       path = ""
    else
       path = fstr(1:n)
    endif
    path = trim(adjustl(path))
  end function coop_file_path_of

  function coop_file_name_of(fstr, want_ext) result(fname)
    COOP_UNKNOWN_STRING  fstr
    COOP_STRING fname
    COOP_INT n
    logical,optional:: want_ext
    n = scan(fstr, '\/', .true.)
    if(n.eq.0)then
       fname = trim(fstr)
    else
       fname = trim(fstr(n+1:))
    endif
    fname = trim(adjustl(fname))
    if(present(want_ext))then
       if(.not. want_ext)then
          n = scan(fname, ".")
          if(n.ne.0)then
             fname = trim(fname(1:n-1))
          endif
       endif
    endif
  end function coop_file_name_of

  function coop_file_add_postfix(fstr, postfix) result(fname)
    COOP_UNKNOWN_STRING fstr, postfix
    COOP_LONG_STRING fname
    COOP_INT n
    n = scan(fstr, ".", .true.)
    if(n.eq.0)then
       fname = trim(fstr)//trim(postfix)
    else
       fname = fstr(1:n-1)//trim(postfix)//trim(fstr(n:))
    endif
    fname = trim(adjustl(fname))
  end function coop_file_add_postfix

  subroutine coop_convert_to_C_string(str)
    COOP_UNKNOWN_STRING str
    COOP_INT i
    do i = len_trim(str)+1, len(str)
       str(i:i) = char(0)
    enddo
  end subroutine coop_convert_to_C_string

  subroutine coop_convert_to_Fortran_String(str)
    COOP_UNKNOWN_STRING str
    COOP_INT i
    do i=1, len(str)
       if(str(i:i).eq.char(0))then
          str(i:i) = " "
       endif
    enddo
  end subroutine coop_convert_to_Fortran_String


  function coop_string_contain_numbers(str) result(nums)
    COOP_UNKNOWN_STRING::str
    COOP_INT nums
    COOP_INT length, i, j
    nums = 0
    Length=scan(Str,"1234567890+-.eE",BACK=.true.)
    if(length .eq. 0) return
    i = scan(str(i:Length), "1234567890-+.eE")
    do while(i.ne. 0)
       nums = nums + 1
       j = verify(Str(i:Length), "1234567890-+.eE")
       if(j.eq.0)return
       i = i + j
       if(i .ge. length)return
       j=scan(Str(i:Length),"1234567890-+.eE")
       if(j.eq.0)return
       i = i + j - 1
    enddo
  end function coop_string_contain_numbers

  !! 0 string; 1 integer; 2 float
  function coop_data_type(str) result(dtype)
    COOP_UNKNOWN_STRING str
    COOP_STRING s
    integer dtype, i, l, idot, ieE, ipm
    s = trim(adjustl(str))
    if(trim(s).eq."T" .or. trim(s).eq."F")then
       dtype = COOP_FORMAT_LOGICAL
       return
    endif
    l = len_trim(s)
    i = verify(s(1:l), "0123456789.+-eE")
    if(i.ne.0)then
       dtype = COOP_FORMAT_STRING
       return
    endif
    i = scan(s(1:l), "0123456789")
    if(i.eq.0)then
       dtype = COOP_FORMAT_STRING
       return
    endif
    i = scan(s(1:l), ".eE")
    if(i.eq.0)then
       i = scan(s(1:l), "+-", .true.)
       if(i.gt.1)then
          dtype = COOP_FORMAT_STRING
       else
          dtype = COOP_FORMAT_INTEGER
       endif
       return
    endif
    ieE = scan(s(1:l), "eE")
    idot = scan(s(1:l), ".")
    ipm = scan(s(1:l), "+-")
    if(idot .ne. 0 .and. ieE.ne.0 .and. idot .gt. ieE)then
       dtype = COOP_FORMAT_STRING
       return
    endif
    if(scan(s(idot+1:l), ".").ne.0 .or. scan(s(ieE+1:l), "eE").ne.0)then
       dtype = COOP_FORMAT_STRING
       return
    endif
    if(ipm .ne. 1 .and. ipm .ne. 0 .and. ipm .ne. ieE+1)then
       dtype = COOP_FORMAT_STRING
       return
    endif
    if(ipm .eq. 1)then
       ipm = scan(s(2:l), "+-")
       if(ipm.ne.0)then
          if(ipm .ne. ieE)then
             dtype = COOP_FORMAT_STRING
             return
          endif
       endif
    endif
    dtype = COOP_FORMAT_FLOAT
    return
  end function coop_data_type

  function coop_string_strip_quotes(str) result(restr)
    COOP_UNKNOWN_STRING::str
    COOP_STRING::restr
    COOP_INT::istart, iend, ic
    ic = scan(str, "/")
    if(ic .eq. 0)then
       ic  = len_trim(str)
    else
       ic = ic - 1
       if(ic.eq.0)then
          restr = ""
          return
       endif
    endif
    istart = verify(str(1:ic), " '""")
    iend = verify(str(1:ic)," '""", .true.)
    if(istart .le.iend)then
       restr  = str(istart:iend)
    else
       restr = ""
    endif
  end function coop_string_strip_quotes

end module coop_string_mod
