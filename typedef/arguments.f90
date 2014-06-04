module coop_arguments_mod
  use coop_constants_mod
  use coop_basicutils_mod
  implicit none
#include "constants.h"

private

public::coop_arguments

  type coop_arguments
     COOP_INT n_int, n_real, n_logical
     COOP_INT,dimension(:),allocatable::i
     COOP_REAL,dimension(:),allocatable::r
     logical,dimension(:),allocatable::l
   contains
     procedure::init => coop_arguments_initialize
     procedure::free => coop_arguments_free
     procedure::clone => coop_arguments_clone
  end type coop_arguments

  interface coop_arguments
     module procedure coop_arguments_constructor
  end interface coop_arguments

contains


  function coop_arguments_constructor(i, r , l) result(this)
    type(coop_arguments) this
    COOP_INT, dimension(:),optional::i
    COOP_REAL, dimension(:),optional::r
    logical, dimension(:),optional::l
    call this%free
    if(present(i))then
       this%n_int = size(i)
       allocate(this%i(this%n_int))
       this%i = i
    endif
    if(present(r))then
       this%n_real = size(r)
       allocate(this%r(this%n_real))
       this%r = r
    endif
    if(present(l))then
       this%n_logical = size(l)
       allocate(this%l(this%n_logical))
       this%l = l
    endif
  end function coop_arguments_constructor

  subroutine coop_arguments_initialize(this,i,r,l)
    class(coop_arguments) this
    COOP_INT, dimension(:),optional::i
    COOP_REAL, dimension(:),optional::r
    logical, dimension(:),optional::l
    call this%free
    if(present(i))then
       this%n_int = size(i)
       allocate(this%i(this%n_int))
       this%i = i
    endif
    if(present(r))then
       this%n_real = size(r)
       allocate(this%r(this%n_real))
       this%r = r
    endif
    if(present(l))then
       this%n_logical = size(l)
       allocate(this%l(this%n_logical))
       this%l = l
    endif
  end subroutine coop_arguments_initialize
    
  subroutine coop_arguments_free(this)
    class(coop_arguments) this
    if(allocated(this%i))deallocate(this%i)
    if(allocated(this%r))deallocate(this%r)
    if(allocated(this%l))deallocate(this%l)
    this%n_int = 0
    this%n_real = 0
    this%n_logical = 0
  end subroutine coop_arguments_free

  subroutine coop_arguments_clone(this, args)
    class(coop_arguments) this
    type(coop_arguments) args
    call this%free()
    if(allocated(args%r))then
       allocate(this%r(args%n_real))
       this%r = args%r
       this%n_real = args%n_real
    endif
    if(allocated(args%i))then
       allocate(this%i(args%n_int))
       this%i = args%i
       this%n_int = args%n_int
    endif
    if(allocated(args%l))then
       allocate(this%l(args%n_logical))
       this%l = args%l
       this%n_logical = args%n_logical
    endif
  end subroutine coop_arguments_clone

end module coop_arguments_mod