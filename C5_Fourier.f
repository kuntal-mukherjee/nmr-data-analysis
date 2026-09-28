      IMPLICIT NONE

      INTEGER, PARAMETER :: NP = 100

      REAL, PARAMETER :: PI = 4.0 * ATAN(1.0)
      !The integrals with labels 1 and 2 refer to the I+ operator and
      !integrals with labels 3 and 4 refer to the I- operator

      COMPLEX A1(NP,NP), A2(NP,NP) ! A1 & B1, integrals from 0 to tp (1st element)

      COMPLEX B1(NP,NP), B2(NP,NP) ! A2 & B2, integrals from tp to 2tp (2nd element)

      COMPLEX A3(NP,NP), A4(NP,NP) ! A3 & B3, integrals from 2tp to 3tp (3rd element)

      COMPLEX B3(NP,NP), B4(NP,NP) ! A4 & B4, integrals from 3tp to 4tp (4th element)

      COMPLEX A5(NP,NP), A6(NP,NP) ! A5 & B5, integrals from 4tp to 5tp (5th element)

      COMPLEX B5(NP,NP), B6(NP,NP) ! A6 & B6, integrals from 5tp to 6tp (6th element)

      COMPLEX A7(NP,NP), A8(NP,NP) ! A7 & B7, integrals from 6tp to 7tp (7th element)

      COMPLEX B7(NP,NP), B8(NP,NP) ! A8 & B8, integrals from 7tp to 8tp (8th element)

      COMPLEX A9(NP,NP), A10(NP,NP) ! A9 & B9, integrals from 8tp to 9tp (9th element)

      COMPLEX B9(NP,NP), B10(NP,NP) ! A10 & B10, integrals from 9tp to 10tp (10th element)

      COMPLEX    SUM1(NP,NP), SUM2(NP,NP)


      REAL    AN1, AN2, tp, WC

      REAL    phi1, phi2, phi3, phi4, phi5, phi6, phi7, phi8
      REAL    phi9, phi10

      REAL    wc_t1, wc_t2, wc_t3, wc_t4, wc_t5, wc_t6, wc_t7, wc_t8
      REAL    wc_t9, wc_t10

      COMPLEX AM

      INTEGER :: K, J, I

      ! Complex unit (imaginary no. defined)
      AM = (0.0, 1.0)

      ! Initialize
      DO J = 1,1   !This is just to form array
         A1(J,J) = (0.0,0.0)
         A2(J,J) = (0.0,0.0)
         A3(J,J) = (0.0,0.0)
         A4(J,J) = (0.0,0.0)
         A5(J,J) = (0.0,0.0)
         A6(J,J) = (0.0,0.0)
         A7(J,J) = (0.0,0.0)
         A8(J,J) = (0.0,0.0)
         B1(J,J) = (0.0,0.0)
         B2(J,J) = (0.0,0.0)
         B3(J,J) = (0.0,0.0)
         B4(J,J) = (0.0,0.0)
         B5(J,J) = (0.0,0.0)
         B6(J,J) = (0.0,0.0)
         B7(J,J) = (0.0,0.0)
         B8(J,J) = (0.0,0.0)

         SUM1(J,J)=(0.0,0.0)
         SUM2(J,J)=(0.0,0.0)
      ENDDO

      K = 5  !Type of Sequence
      AN1 = -2.0              ! n-value for I+ operator
      AN2 = -2.0              ! n-value for I- operator

      phi1 = 0.0             ! phase of 1st element
      phi2= PI               ! phase of 2nd element
      phi3= 2.0*PI/5.0       ! phase of 3rd element
      phi4= PI+2.0*PI/5.0    ! phase of 4th element
      phi5= 4.0*PI/5.0       ! phase of 5th element
      phi6= PI+4.0*PI/5.0    ! phase of 6th element
      phi7= 6.0*PI/5.0       ! phase of 7th element
      phi8= PI+6.0*PI/5.0    ! phase of 8th element
      phi9= 8.0*PI/5.0       ! phase of 9th element
      phi10= PI+8.0*PI/5.0   ! phase of 10th element

      wc_t1 = PI/5.0            ! corresponds to width of one element, tp
      wc_t2 = 2.0*PI/5.0        ! corresponds to 2tp
      wc_t3 = 3.0*PI/5.0        ! corresponds to 3tp
      wc_t4 = 4.0*PI/5.0        ! corresponds to 4tp
      wc_t5 = PI                ! corresponds to 5tp
      wc_t6 = 6.0*PI/5.0        ! corresponds to 6tp
      wc_t7 = 7.0*PI/5.0        ! corresponds to 7tp
      wc_t8 = 8.0*PI/5.0        ! corresponds to 8tp
      wc_t9 = 9.0*PI/5.0        ! corresponds to 9tp
      wc_t10 = 2.0*PI           ! corresponds to 10tp



    ! Calculating integrals associated with I+ operator

      DO J = 1,1

         A1(J,J) = EXP(-AM*phi1) *
     1   ( EXP(-AM*AN1*wc_t1) - EXP(-AM*AN1*0.0) ) + A1(J,J)    ! 1st integral

          A2(J,J) = EXP(-AM*phi2) *
     1   ( EXP(-AM*AN1*wc_t2) - EXP(-AM*AN1*wc_t1) ) + A2(J,J)     ! 2nd integral

         A3(J,J) = EXP(-AM*phi3) *
     1   ( EXP(-AM*AN1*wc_t3) - EXP(-AM*AN1*wc_t2) ) + A3(J,J)    ! 3rd integral

          A4(J,J) = EXP(-AM*phi4) *
     1   ( EXP(-AM*AN1*wc_t4) - EXP(-AM*AN1*wc_t3) ) + A4(J,J)     ! 4th integral

         A5(J,J) = EXP(-AM*phi5) *
     1   ( EXP(-AM*AN1*wc_t5) - EXP(-AM*AN1*wc_t4) ) + A5(J,J)    ! 5th integral

          A6(J,J) = EXP(-AM*phi6) *
     1   ( EXP(-AM*AN1*wc_t6) - EXP(-AM*AN1*wc_t5) ) + A6(J,J)     ! 6th integral

          A7(J,J) = EXP(-AM*phi7) *
     1   ( EXP(-AM*AN1*wc_t7) - EXP(-AM*AN1*wc_t6) ) + A7(J,J)    ! 7th integral

          A8(J,J) = EXP(-AM*phi8) *
     1   ( EXP(-AM*AN1*wc_t8) - EXP(-AM*AN1*wc_t7) ) + A8(J,J)     ! 8th integral

          A9(J,J) = EXP(-AM*phi9) *
     1   ( EXP(-AM*AN1*wc_t9) - EXP(-AM*AN1*wc_t8) ) + A9(J,J)    ! 9th integral

          A10(J,J) = EXP(-AM*phi10) *
     1   ( EXP(-AM*AN1*wc_t10) - EXP(-AM*AN1*wc_t9) ) + A10(J,J)    ! 10th integral



      SUM1(J,J)=AM*0.5*(A1(J,J)+A2(J,J)+A3(J,J)+A4(J,J)+A5(J,J)+A6(J,J)
     1  +A7(J,J)+A8(J,J)+A9(J,J)+A10(J,J))
      ENDDO

    ! calculating integrals associated with I- operator
      DO J = 1,1

         B1(J,J) = EXP(AM*phi1) *
     1   ( EXP(-AM*AN2*wc_t1) - EXP(-AM*AN2*0.0) ) + B1(J,J)    ! 1st integral

          B2(J,J) = EXP(AM*phi2) *
     1   ( EXP(-AM*AN2*wc_t2) - EXP(-AM*AN2*wc_t1) ) + B2(J,J)     ! 2nd integral

          B3(J,J) = EXP(AM*phi3) *
     1   ( EXP(-AM*AN2*wc_t3) - EXP(-AM*AN2*wc_t2) ) + B3(J,J)    ! 3rd integral

          B4(J,J) = EXP(AM*phi4) *
     1   ( EXP(-AM*AN2*wc_t4) - EXP(-AM*AN2*wc_t3) ) + B4(J,J)     ! 4th integral

         B5(J,J) = EXP(AM*phi5) *
     1   ( EXP(-AM*AN2*wc_t5) - EXP(-AM*AN2*wc_t4) ) + B5(J,J)    ! 5th integral

          B6(J,J) = EXP(AM*phi6) *
     1   ( EXP(-AM*AN2*wc_t6) - EXP(-AM*AN2*wc_t5) ) + B6(J,J)     ! 6th integral

         B7(J,J) = EXP(AM*phi7) *
     1   ( EXP(-AM*AN2*wc_t7) - EXP(-AM*AN2*wc_t6) ) + B7(J,J)    ! 7th integral

          B8(J,J) = EXP(AM*phi8) *
     1   ( EXP(-AM*AN2*wc_t8) - EXP(-AM*AN2*wc_t7) ) + B8(J,J)     ! 8th integral

          B9(J,J) = EXP(AM*phi9) *
     1   ( EXP(-AM*AN2*wc_t9) - EXP(-AM*AN2*wc_t8) ) + B9(J,J)    ! 9th integral

          B10(J,J) = EXP(AM*phi10) *
     1   ( EXP(-AM*AN2*wc_t10) - EXP(-AM*AN2*wc_t9) ) + B10(J,J)   ! 10th integral

      SUM2(J,J)=AM*0.5*(B1(J,J)+B2(J,J)+B3(J,J)+B4(J,J)+B5(J,J)+B6(J,J)
     1   +B7(J,J)+B8(J,J)+B9(J,J)+B10(J,J))
      ENDDO

      PRINT "(A,I0,A)", "C", K, "-SEQUENCE"
c      WRITE(*,*) A1(1,1), A2(1,1), A3(1,1), A4(1,1), A5(1,1), A6(1,1)

      PRINT "(A,F6.2,A)", "For n = ", AN1, ", the coefficient of I+ is"
      WRITE(*,'(A,2F12.6)') "SUM1 = ", REAL(SUM1(1,1)), AIMAG(SUM1(1,1))

c      WRITE(*,*) B1(1,1), B2(1,1), B3(1,1), B4(1,1), B5(1,1), B6(1,1)

      PRINT "(A,F6.2,A)", "For n = ", AN2, ", the coefficient of I- is"
      WRITE(*,'(A,2F12.6)') "SUM2 = ", REAL(SUM2(1,1)), AIMAG(SUM2(1,1))
      END
