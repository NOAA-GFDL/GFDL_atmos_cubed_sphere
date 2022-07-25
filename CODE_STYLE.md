# Coding Style

## General

* Trim all trailing whitespace from every line (some editors can do this
  automatically).
* No tab characters.
* A copy of the FV3 Gnu Lesser General Public License Header
  must be included at the top of each file.
* Supply an author block for each file with a description of the file and the author(s)
  name or GitHub ID.
* Documentation may be written so that it can be parsed by [Doxygen](http://www.doxygen.nl/).
* All variables should be defined, and include units. Unit-less variables should be marked `unitless`
* Provide detailed descriptions of modules, interfaces, functions, and subroutines
* Define all function/subroutine arguments, and function results (see below)
* Follow coding style of the current file, as much as possible.

## Fortran

### General

* Use Fortran 95 standard or newer
* Two space indentation
* Never use implicit variables (i.e., always specify `IMPLICIT NONE`)
* Lines must be <= 120 characters long (including comments)
* logical, compound logical, and relational if statements may be one line,
  using “&” for line continuation if necessary:
  ```Fortran
  if(file_exists(fileName)) call open_file(fileObj,fileName, is_restart=.false)
  ```
* Avoid the use of `GOTO` statements
* Avoid the use of Fortran reserved words as variables (e.g. `DATA`, `NAME`)
* `COMMON` blocks should never be used

### Derived types

* Type names must be in CapitalWord format.
* Description on the line before the type definition.
* Inline doxygen descriptions for all member variables.

## Functions
* Functions should include a result variable on its own line, that does not have
  a specific intent.
* Inline descriptions for all arguments, except the result variable.
* Description on the line(s) before the function definition.  Specify what the function is returning (with the `@return` doxygen keyword if using doxygen).

## Blocks
* terminate `do` loops with `enddo`
* terminate block `if`, `then` statements with `endif`

## OpenMP

* All openMP directives should specify default(none), and then explicitly list
  all shared and private variables.
* All critical sections must have a unique name.

## Fortran Example

```Fortran

!***********************************************************************
!*                   GNU Lesser General Public License
!*
!* This file is part of the FV3 dynamical core.
!*
!* The FV3 dynamical core is free software: you can redistribute it
!* and/or modify it under the terms of the
!* GNU Lesser General Public License as published by the
!* Free Software Foundation, either version 3 of the License, or
!* (at your option) any later version.
!*
!* The FV3 dynamical core is distributed in the hope that it will be
!* useful, but WITHOUT ANY WARRANTY; without even the implied warranty
!* of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
!* See the GNU General Public License for more details.
!*
!* You should have received a copy of the GNU Lesser General Public
!* License along with the FV3 dynamical core.
!* If not, see <http://www.gnu.org/licenses/>.

!***********************************************************************

!> @file
!! @brief Example code
!! @author <developer>
!! @email <email>

module example_mod
  use util_mod, only: util_func1
  implicit none
  private

  public :: sub1
  public :: func1

  !> @brief Doxygen description of type.
  type,public :: CustomType
    integer(kind=<KIND>) :: a_var !< Inline doxygen description.
    real(kind=<KIND>),dimension(:),allocatable :: b_arr !< long description
                                                        !! continued on
                                                        !! multiple lines.
  endtype CustomType

  contains

  !> @brief Doxygen description.
  subroutine sub1(arg1, &
    & arg2, &
    & arg3)
    real(kind=<KIND>),intent(in) :: arg1 !< Inline doxygen description.
    integer(kind=<KIND>),intent(inout) :: arg2 !< Inline doxygen description.
    character(len=*),intent(inout) :: arg3 !< Long inline doxygen
                                           !! description.
  end subroutine sub1

  !> @brief Doxygen description
  !! @return Function return value.
  function func1(arg1, &
    & arg2) &
    & result(res)
    integer(kind=<KIND>),intent(in) :: arg1 !< Inline doxygen description
    integer(kind=<KIND>),intent(in) :: arg2 !< Inline doxygen description
    integer(kind=<KIND>) :: res
  end function func1

end module example_mod
```
