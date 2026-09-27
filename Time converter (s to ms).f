      program scale_first_column
      implicit none

      real(8) :: t, re, im
      integer :: ios

      open(unit=10, file="fid.dat", status="old", action="read")
      open(unit=20, file="fid_new.dat", status="unknown")


      do
        read(10, *, iostat=ios) t, re, im
        if (ios /= 0) exit

        t = t * 1000.0d0

        write(20, '(F12.6, 1X, ES15.6, 1X, ES15.6)') t, re, im
      end do

      close(10)
      close(20)

      end program scale_first_column
