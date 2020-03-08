! Willem is een vette baas

PROGRAM GET_3D_CUBE
USE ISO_C_BINDING
IMPLICIT NONE

INTERFACE
  ! Get the extensions of the data cube
  SUBROUTINE getdims(nx,ny,nz) BIND(C)
    IMPORT :: C_INT
    INTEGER(C_INT) :: nx,ny,nz
  END SUBROUTINE getdims

  ! Get the actual data
  SUBROUTINE getarr(ptr) BIND(C)
    IMPORT :: C_DOUBLE
    REAL(C_DOUBLE), DIMENSION(*), INTENT(INOUT) :: ptr
  END SUBROUTINE getarr

  ! Load the primordial field
  SUBROUTINE load_primordial_field(fname) BIND(C)
    IMPORT :: C_CHAR
    CHARACTER(KIND=C_CHAR) :: fname(*)
  END SUBROUTINE load_primordial_field

  ! Run CLASS
  SUBROUTINE run_class(fname_ini, fname_pre) BIND(C)
    IMPORT :: C_CHAR
    CHARACTER(KIND=C_CHAR) :: fname_ini(*)
    CHARACTER(KIND=C_CHAR) :: fname_pre(*)
  END SUBROUTINE run_class
END INTERFACE

INTEGER :: nx, ny, nz                               ! Extensions of the data cube
REAL (C_DOUBLE), ALLOCATABLE, TARGET :: cube(:,:,:) ! The data cube

CHARACTER(LEN = :), ALLOCATABLE :: fname_primordial_field   ! primordial Gaussian field
CHARACTER(LEN = :), ALLOCATABLE :: fname_class_ini          ! CLASS parameter file
CHARACTER(LEN = :), ALLOCATABLE :: fname_class_pre          ! optional precision file

! Retrieve the dimensions and allocate the array
CALL getdims(nx,ny,nz)
WRITE(*,*) "The dimensions are", nx, ny, nz
ALLOCATE(cube(nx,ny,nz))

! Set everything to 1
cube = 10
cube(1,1,1) = 2

! Get the array from the C code
WRITE(*,*) cube(1,1,1)
CALL getarr(cube)
WRITE(*,*) cube(1,1,1)

! Free up memory
DEALLOCATE(cube)

fname_primordial_field = "gaussian_pure.hdf5"   ! primordial Gaussian field
fname_class_ini = "neutrino_example.ini"        ! CLASS parameter file
fname_class_pre = ""                            ! optional precision file

CALL load_primordial_field(fname_primordial_field)
CALL run_class(fname_class_ini, fname_class_pre)


ENDPROGRAM GET_3D_CUBE
