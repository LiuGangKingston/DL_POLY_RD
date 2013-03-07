      module atomicvirial_module

c----------------------------------------------------------
c     
c     Module for the computation of the atomic stress
c     virials.
c     
c     Author  R. Darkins - University College London
c     Date    Mar 2013
c     
c     NB. (xx, xy, xz, yy, yz, zz)
c     
c----------------------------------------------------------

      implicit none
      private

c----------------------------------------------------------
c     Private routines
c----------------------------------------------------------

      public :: init_atomicvirial
      public :: reset_atomicvirial
      public :: update_atomicvirial
      public :: units_atomicvirial
      public :: cleanup_atomicvirial

c----------------------------------------------------------
c     Public variables
c----------------------------------------------------------

      public :: lvirial
      public :: lvirupdate
      public :: atmvir
      public :: thisvir
      public :: ONE_THIRD
      public :: dbgv
      
      logical,save :: lvirial = .false.
      logical,save :: lvirupdate = .false.
      real(8), dimension(:,:), allocatable, save :: atmvir
      real(8), dimension(6), save :: thisvir
      real(8), save :: ONE_THIRD = 1.0d0 / 3.0d0
      real(8), save :: dbgv = 0.0d0

c----------------------------------------------------------
c     Internal miscellany
c----------------------------------------------------------

      integer :: ierr = 0

      contains
      
c----------------------------------------------------------
c     Subroutines
c----------------------------------------------------------

      subroutine init_atomicvirial(natms)

      implicit none
      integer, intent(in) :: natms
      integer :: i,j

      allocate(atmvir(6,natms), stat=ierr)
      do i=1,natms
         do j=1,6
            atmvir(j,i)=0.0d0
         enddo
      enddo

      return

      end subroutine init_atomicvirial

      subroutine reset_atomicvirial(natms)

      implicit none
      integer, intent(in) :: natms
      integer :: i,j

      do i=1,natms
         do j=1,6
            atmvir(j,i)=0.0d0
         enddo
      enddo

      return

      end subroutine reset_atomicvirial

c----------------------------------------------------------
c     
c     cid (for debugging purposes):
c     
c      0 - angles
c      1 - bonds
c      2 - coulomb (including ewald)
c      3 - dihedral
c      4 - four-body
c      5 - inversion
c      6 - three-body
c      7 - vdw
c      8 - neu_coul
c      9 - kinetic
c     10 - core shell
c
c----------------------------------------------------------
      
      subroutine update_atomicvirial(atmid,cid)
      
      implicit none
      integer, intent(in) :: atmid
      integer, intent(in) :: cid
      integer :: j
      
      do j=1,6
         atmvir(j,atmid)=atmvir(j,atmid)+thisvir(j)
      enddo

      return
      
      end subroutine update_atomicvirial

c     Convert units to k-atm
      subroutine units_atomicvirial(natms)
      
      implicit none
      integer, intent(in) :: natms
      integer :: i,j

      do i=1,natms
         do j=1,6
            atmvir(j,i)=atmvir(j,i)*0.163882576d0
         enddo
      enddo

      return
      
      end subroutine units_atomicvirial

      subroutine cleanup_atomicvirial
      
      implicit none

      deallocate(atmvir)

      return
      
      end subroutine cleanup_atomicvirial

      end module atomicvirial_module

