C*****************************************************************
C    PROGRAM TO CALCULATE THE MAGNETIZATION EXCHANGE IN ROTATIONAL
C           RESONANCE EXPERIMENTS (TWO-SPIN SYSTEM, C-C)
C     SPINNING FREQUENCY DEPENDENT EXPERIMENTS
C             N=1 ROTATIONAL RESONANCE CONDITION
C*****************************************************************

        IMPLICIT REAL*8 (A-H,O-Z)
        PARAMETER (NS=2,IP=6044,IP1=1000,MF=4)
C*****************************************************************
C------DECLARATION OF VARIABLES  USED IN THE PROGRAM--------------
C------------------------ PART-1 --------------------------------

        COMPLEX AM
        DIMENSION A(IP),B(IP),G(IP),AREA(IP)
        DIMENSION RES1(IP1),RES2(IP1)
        COMPLEX ARES1(IP1),ARES2(IP1)

        COMPLEX SM3,SM3A,SM3B
        COMPLEX SM4,SM4A,SM4B
        COMPLEX SM5,SM5A,SM5B

        COMPLEX Zpp12,Zmm12,Epp12,Emm12

        DIMENSION JS(NS)

C*****************************************************************
C-----------------------------------------------------------------
C-------DECLARATION OF VARIABLES FOR CSA  ------------------------
C------------------------PART-2 ----------------------------------

        DIMENSION RNU(NS),DELTA(NS)
        DIMENSION RHO20(NS),RHO2P2(NS),RHO2M2(NS)
        DIMENSION AL_CSA(NS),BET_CSA(NS),GAM_CSA(NS)

        COMPLEX GPCS(NS,MF),GMCS(NS,MF),CPCS(NS,MF),CMCS(NS,MF)

        COMPLEX GCP1_20(NS),GCP2_20(NS),GCM1_20(NS),GCM2_20(NS)
        COMPLEX GCP1_2P2(NS),GCP2_2P2(NS),GCM1_2P2(NS),GCM2_2P2(NS)
        COMPLEX GCP1_2M2(NS),GCP2_2M2(NS),GCM1_2M2(NS),GCM2_2M2(NS)
C-----------------------------------------------------------------
          COMPLEX GP_P(NS,MF),GM_P(NS,MF),GZ_P(NS,MF)
        COMPLEX CP_P(NS,MF),CM_P(NS,MF),CZ_P(NS,MF)

        COMPLEX GP_M(NS,MF),GM_M(NS,MF),GZ_M(NS,MF)
        COMPLEX CP_M(NS,MF),CM_M(NS,MF),CZ_M(NS,MF)

        COMPLEX GX_P(NS,MF),GX_M(NS,MF),GY_P(NS,MF),GY_M(NS,MF)
        COMPLEX CX_P(NS,MF),CX_M(NS,MF),CY_P(NS,MF),CY_M(NS,MF)

C----------------------------------------------------------------
        COMPLEX RX_P(NS,MF),RX_M(NS,MF),RY_P(NS,MF),RY_M(NS,MF)
        COMPLEX temp1, temp2, temp11,temp21

C*****************************************************************
C----- DECLARATION OF VARIABLES FOR DIPOLAR INTERACTIONS ---------
C----------------------- PART-3 ----------------------------------

         DIMENSION R(NS,NS),SF(NS,NS)

         COMPLEX wig_mat(5,5),wig_pm(5,5),wig_mr(5,5)

         COMPLEX W_D1(NS,NS),W_D2(NS,NS),W_D11(NS,NS),W_D21(NS,NS)
         COMPLEX G1D(NS,NS),G2D(NS,NS),G11D(NS,NS),G21D(NS,NS)
         COMPLEX GP1(NS,NS),GP2(NS,NS),GM1(NS,NS),GM2(NS,NS)

C------------ PAIR 1-2--------------------------------------------

         COMPLEX wig_pm12(5,5)
          COMPLEX Gzz12_P(MF),Gzz12_M(MF),Gpm12_P(MF),Gpm12_M(MF)
         COMPLEX Gmp12_P(MF),Gmp12_M(MF),Gpp12_P(MF),Gpp12_M(MF)
         COMPLEX Czz12_P(MF),Czz12_M(MF),Cpm12_P(MF),Cpm12_M(MF)
         COMPLEX Cmp12_P(MF),Cmp12_M(MF),Cpp12_P(MF),Cpp12_M(MF)
         COMPLEX Gmm12_P(MF),Gmm12_M(MF),Cmm12_P(MF),Cmm12_M(MF)


C*****************************************************************
C----------END OF DECLARATION OF VARIABLES------------------------
C-----------------------------------------------------------------


C--------------------END INPUT FILES------------------------------

C------------SOME INTERNAL DEFINITIONS----------------------------
         PI=4.0*ATAN(1.0)
c        PI=22.0/7.0
         AM=(0.0,1.0)

         THETAH=54.7356*pi/180.0
         CTH=COS(THETAH)
      	 STH=SIN(THETAH)
     	 CTH2=(CTH)**2
      	 STH2=1.0-CTH2
      	 S2TH=2.0*CTH*STH

C----------------------------------------------------------------
C--- BEGIN INITIALIZATION OF ARRAYS -- PART-1--------------------------

            DO I=1,IP1
               ARES1(I)=cmplx(0.0,0.0)
               ARES2(I)=cmplx(0.0,0.0)

               RES1(I)=0.0
               RES2(I)=0.0
	      ENDDO

C----- END OF INITIALIZATION OF ARRAYS -- PART-1-------------------
C******************************************************************

C--- BEGIN INITIALIZATION OF ARRAYS -- PART-2--------------------------

           DO I=1,NS
               DELTA(I)=0.0
               RNU(I)=0.0
               RHO20(I)=0.0
               RHO2M2(I)=0.0
               RHO2P2(I)=0.0

               JS(I)=0
            enddo
C-----------------------
C--------------------------------------------------------------
             DO I=1,NS
        	  DO J=1,MF

              GZ_P(I,J)=CMPLX(0.0,0.0)
            GP_P(I,J)=CMPLX(0.0,0.0)
            GM_P(I,J)=CMPLX(0.0,0.0)

            CZ_P(I,J)=CMPLX(0.0,0.0)
	        CP_P(I,J)=CMPLX(0.0,0.0)
	        CM_P(I,J)=CMPLX(0.0,0.0)

            GZ_M(I,J)=CMPLX(0.0,0.0)
	        GP_M(I,J)=CMPLX(0.0,0.0)
	        GM_M(I,J)=CMPLX(0.0,0.0)

            CZ_M(I,J)=CMPLX(0.0,0.0)
	        CP_M(I,J)=CMPLX(0.0,0.0)
	        CM_M(I,J)=CMPLX(0.0,0.0)


            GX_P(I,J)=CMPLX(0.0,0.0)
            GX_M(I,J)=CMPLX(0.0,0.0)
            CX_P(I,J)=CMPLX(0.0,0.0)
            CX_M(I,J)=CMPLX(0.0,0.0)

            GY_P(I,J)=CMPLX(0.0,0.0)
            GY_M(I,J)=CMPLX(0.0,0.0)
            CY_P(I,J)=CMPLX(0.0,0.0)
            CY_M(I,J)=CMPLX(0.0,0.0)

            RX_P(I,J)=CMPLX(0.0,0.0)
            RX_M(I,J)=CMPLX(0.0,0.0)
            RY_P(I,J)=CMPLX(0.0,0.0)
            RY_M(I,J)=CMPLX(0.0,0.0)

            GPCS(I,J)=CMPLX(0.0,0.0)
	        GMCS(I,J)=CMPLX(0.0,0.0)
	        CPCS(I,J)=CMPLX(0.0,0.0)
	        CMCS(I,J)=CMPLX(0.0,0.0)

                    ENDDO
		    ENDDO
C------------------------------------------------------------------

            DO I=1,NS

             GCP1_20(I)=CMPLX(0.0,0.0)
             GCP2_20(I)=CMPLX(0.0,0.0)
             GCM1_20(I)=CMPLX(0.0,0.0)
             GCM2_20(I)=CMPLX(0.0,0.0)

             GCP1_2P2(I)=CMPLX(0.0,0.0)
             GCP2_2P2(I)=CMPLX(0.0,0.0)
             GCM1_2P2(I)=CMPLX(0.0,0.0)
             GCM2_2P2(I)=CMPLX(0.0,0.0)

             GCP1_2M2(I)=CMPLX(0.0,0.0)
             GCP2_2M2(I)=CMPLX(0.0,0.0)
             GCM1_2M2(I)=CMPLX(0.0,0.0)
             GCM2_2M2(I)=CMPLX(0.0,0.0)

            ENDDO
C----- END OF INITIALIZATION OF ARRAYS -- PART-2-------------------
C******************************************************************

C--- BEGIN INITIALIZATION OF ARRAYS -- PART-3----------------------
C---------INITIALIZE THE WIGNER ROTATION MATRICES-----------------
            DO i=1,5
                DO j=1,5
             wig_mat(i,j)=cmplx(0.0,0.0)
             wig_pm(i,j)=cmplx(0.0,0.0)
             wig_mr(i,j)=cmplx(0.0,0.0)
                ENDDO
            ENDDO
C-----------------------------------------------------------------
             DO I=1,NS
              DO J=1,NS
                  R(I,J)=0.0
                  SF(I,J)=0.0
              ENDDO
             ENDDO
C-----------------------------------------------------------------
            DO I=1,NS
             DO J=1,NS

               W_D1(I,J)=CMPLX(0.0,0.0)
               W_D11(I,J)=CMPLX(0.0,0.0)
               W_D2(I,J)=CMPLX(0.0,0.0)
               W_D21(I,J)=CMPLX(0.0,0.0)

               G1D(I,J)=CMPLX(0.0,0.0)
               G2D(I,J)=CMPLX(0.0,0.0)
               G11D(I,J)=CMPLX(0.0,0.0)
               G21D(I,J)=CMPLX(0.0,0.0)

               GP1(I,J)=CMPLX(0.0,0.0)
               GP2(I,J)=CMPLX(0.0,0.0)
               GM1(I,J)=CMPLX(0.0,0.0)
               GM2(I,J)=CMPLX(0.0,0.0)

             ENDDO
            ENDDO
C----------------------------------------------------------------
           DO I=1,MF

            Gzz12_P(I)=CMPLX(0.0,0.0)
             Gzz12_M(I)=CMPLX(0.0,0.0)
             Gpm12_P(I)=CMPLX(0.0,0.0)
             Gpm12_M(I)=CMPLX(0.0,0.0)
             Gmp12_P(I)=CMPLX(0.0,0.0)
             Gmp12_M(I)=CMPLX(0.0,0.0)
             Gpp12_P(I)=CMPLX(0.0,0.0)
             Gpp12_M(I)=CMPLX(0.0,0.0)
             Gmm12_P(I)=CMPLX(0.0,0.0)
             Gmm12_M(I)=CMPLX(0.0,0.0)

             Czz12_P(I)=CMPLX(0.0,0.0)
             Czz12_M(I)=CMPLX(0.0,0.0)
             Cpm12_P(I)=CMPLX(0.0,0.0)
             Cpm12_M(I)=CMPLX(0.0,0.0)
             Cmp12_P(I)=CMPLX(0.0,0.0)
             Cmp12_M(I)=CMPLX(0.0,0.0)
             Cpp12_P(I)=CMPLX(0.0,0.0)
             Cpp12_M(I)=CMPLX(0.0,0.0)
             Cmm12_P(I)=CMPLX(0.0,0.0)
             Cmm12_M(I)=CMPLX(0.0,0.0)


            ENDDO
C---------------------------------------------------------------
          DO i=1,5
            DO j=1,5
             wig_pm12(i,j)=cmplx(0.0,0.0)

            ENDDO
          ENDDO
C---------------------------------------------------------------
C----- END OF INITIALIZATION OF ARRAYS -- PART-3----------------
C**************************************************************
C-------------END OF INITILIZATION------------------------------

C--------------------BEGIN INPUT FILES---------------------------

C         OPEN(11,FILE='angle50.dat')
C         OPEN(11,FILE='angle300.dat')
         OPEN(11,FILE='angle6044.dat')
         OPEN(12,FILE='mHORROR_N1inp.dat')
         OPEN(15,FILE='mHORROR_N1out.dat')

C-------------END OF INPUT FILES--------------------------------

C------------BEGIN INPUT PARAMETERS-----------------------------

         DO I=1,NS
         READ(12,*)DELTA(I),RNU(I),JS(I)
         READ(12,*)AL_CSA(I),BET_CSA(I),GAM_CSA(I)
         ENDDO
         READ(12,*)AD12,BD12,GD12,R(1,2),SF(1,2)
         READ(12,*)W10,W20
         READ(12,*)WRF1,WRF2
C         READ(12,*)TMIX
C         READ(12,*)NPOINTS,WR_INT,WR_INC

         READ(12,*)WR
         READ(12,*)NPOINTS,TMIX_INT,TMIX_INC

C-------------------END OF INPUT PARAMETERS--------------------

C************** BEGIN  READING ORIENTATIONS *******************
          DO I=1,IP
           READ(11,*)A(I),B(I),G(I),AREA(I)
          ENDDO
C-------------------END OF READING ORIENTATIONS ----------------

C**************************************************************
C--------- CALCULATIONS BEGIN FROM HERE ***********************
C**************************************************************

C------------CALCULATION OF  CSA PARAMETERS--------------------
             DO I=1,NS
             RHO20(I)=DELTA(I)
             RHO2M2(I)=-0.4082*DELTA(I)*RNU(I)
             RHO2P2(I)=RHO2M2(I)
             ENDDO
C**************************************************************

C----CALCULATE DIPOLAR COEFFICIENTS FROM PAS TO MOL-AXIS------

C********************* DIPOLAR 1-2  *****************************

             AD12=AD12*pi/180.0
             BD12=BD12*pi/180.0
             GD12=GD12*pi/180.0

         CALL WIGNER(AD12,BD12,GD12,Wig_mat)

              DO k=1,5
               DO m=1,5
             wig_pm12(k,m)=wig_mat(k,m)

               ENDDO
              ENDDO
C------END OF DIPOLAR PAS TO MOL AXIS CALCULATIONS--------------

C***************************************************************
C---------------BEGIN POWDER AVERRAGING--------------------------

 	    DO iangle=1,IP

C        WRITE(*,*)IANGLE
      	  alpha=a(iangle)
      	  beta=b(iangle)
      	  gamma=g(iangle)

            alfa=alpha
            beita=beta
            gm_pm=gamma
c----------------- for single orientation ----------
c            alfa=alpha*0.0
c            beita=beta*0.0
c            gm_pm=gamma*0.0


         CALL WIGNER(alfa,beita,gm_pm,Wig_mat)

              DO k=1,5
               DO m=1,5

             wig_mr(k,m)=Wig_mat(k,m)

               ENDDO
              ENDDO
C---------------------------------------------------

            DO I=1,NS
             DO J=I+1,NS

              gp1(I,J)=cmplx(0.0,0.0)
              gp2(I,J)=cmplx(0.0,0.0)
	      gm1(I,J)=cmplx(0.0,0.0)
              gm2(I,J)=cmplx(0.0,0.0)

		    ENDDO
		     ENDDO

c     calculate dipolar g1 coefficients for all the pairs

            DO k=1,5
              gp1(1,2)=wig_pm12(3,k)*wig_mr(k,2)+gp1(1,2)
              gp2(1,2)=wig_pm12(3,k)*wig_mr(k,1)+gp2(1,2)
              gm1(1,2)=wig_pm12(3,k)*wig_mr(k,4)+gm1(1,2)
              gm2(1,2)=wig_pm12(3,k)*wig_mr(k,5)+gm2(1,2)
            ENDDO

c-----------------------------------------------------------
          DO I=1,NS
	     DO J=I+1,NS

          g1d(I,J)=sqrt(3.0/8.0)*s2th*gp1(I,J)
          g2d(I,J)=sqrt(3.0/8.0)*sth2*gp2(I,J)

          g11d(I,J)=-sqrt(3.0/8.0)*s2th*gm1(I,J)
          g21d(I,J)=sqrt(3.0/8.0)*sth2*gm2(I,J)

	      ENDDO
		    ENDDO
C-------------------------------------------------------------

C----------------DIPOLAR CONSTANTS---------------------------------

C           b=mu*gamma*gamma*hbar/(4*Pi*r^3)  (rad/s)

C           mu=4*Pi  *e-7  T^2 m^3 J^-1
C           gamma(proton)      = 26.752    *e7 rad/ s T
c           gamma(carbon-13)   =  6.72828  *e7 rad/ s T
c           gamma(nitrogen-15) = -2.71262  *e7 rad/ s T
C           hbar=1.054e-34 J s
C           h=6.626e-34 Js
c
c     In this program the dipolar coupling is given in kHz (i.e b/2Pi)
c
c        wd1=119.93 kHz  (for two protons separated by 1 Angstorm)
c        wd1=7.5859 kHz  (for two carbons separated by 1 Angstorm)
c        wd1=3.058  kHz  (for carbon-nitrogen separated by 1 Angstorm)
c        wd1=30.163 kHz  (for carbon-hydrogen separated by 1 Angstorm)
c        wd1=12.161 kHz  (for nitrogen-hydrogen separated by 1 Angstorm)
C-------------------------------------------------------------------

              RS=-1.0
              wd12=7.5859

C---------------DIPOLAR COEFF FOR PAIR - 1 & 2  ------------------

              w_d1(1,2)=(RS*wd12*g1d(1,2))/(r(1,2)**3)
              w_d2(1,2)=(RS*wd12*g2d(1,2))/(r(1,2)**3)
              w_d11(1,2)=(RS*wd12*g11d(1,2))/(r(1,2)**3)
              w_d21(1,2)=(RS*wd12*g21d(1,2))/(r(1,2)**3)

C----------------------------------------------------------------
C************** END OF DIPOLAR CALCULATIONS *********************
C----------------------------------------------------------------
C****** CALCULATES G1,G2 COEFFICIENTS FOR CSA INTERACTION ********

          DO i=1,NS

              alfa=al_csa(i)*pi/180.0
              beita=bet_csa(i)*pi/180.0
              gm_pm=gam_csa(i)*pi/180.0

         call Wigner(alfa,beita,gm_pm,Wig_mat)

              do k=1,5
               do m=1,5
             wig_pm(k,m)=wig_mat(k,m)
               enddo
              enddo
c-----------------------------------------------

              gcp1_20(i)=cmplx(0.0,0.0)
              gcp1_2p2(i)=cmplx(0.0,0.0)
              gcp1_2m2(i)=cmplx(0.0,0.0)
              gcp2_20(i)=cmplx(0.0,0.0)
              gcp2_2p2(i)=cmplx(0.0,0.0)
              gcp2_2m2(i)=cmplx(0.0,0.0)

              gcm1_20(i)=cmplx(0.0,0.0)
              gcm1_2p2(i)=cmplx(0.0,0.0)
              gcm1_2m2(i)=cmplx(0.0,0.0)
              gcm2_20(i)=cmplx(0.0,0.0)
              gcm2_2p2(i)=cmplx(0.0,0.0)
              gcm2_2m2(i)=cmplx(0.0,0.0)

c  calculate  g1 coefficients

            do k=1,5
              gcp1_20(i)=wig_pm(3,k)*wig_mr(k,2)+gcp1_20(i)
              gcp2_20(i)=wig_pm(3,k)*wig_mr(k,1)+gcp2_20(i)
              gcm1_20(i)=wig_pm(3,k)*wig_mr(k,4)+gcm1_20(i)
              gcm2_20(i)=wig_pm(3,k)*wig_mr(k,5)+gcm2_20(i)
            enddo

            do k=1,5
              gcp1_2m2(i)=wig_pm(1,k)*wig_mr(k,2)+gcp1_2m2(i)
              gcp2_2m2(i)=wig_pm(1,k)*wig_mr(k,1)+gcp2_2m2(i)
              gcm1_2m2(i)=wig_pm(1,k)*wig_mr(k,4)+gcm1_2m2(i)
              gcm2_2m2(i)=wig_pm(1,k)*wig_mr(k,5)+gcm2_2m2(i)
            enddo

            do k=1,5
              gcp1_2p2(i)=wig_pm(5,k)*wig_mr(k,2)+gcp1_2p2(i)
              gcp2_2p2(i)=wig_pm(5,k)*wig_mr(k,1)+gcp2_2p2(i)
              gcm1_2p2(i)=wig_pm(5,k)*wig_mr(k,4)+gcm1_2p2(i)
              gcm2_2p2(i)=wig_pm(5,k)*wig_mr(k,5)+gcm2_2p2(i)
            enddo

c------------------------------------------------------------

         temp1=rho20(i)*gcp1_20(i)+rho2m2(i)*gcp1_2m2(i)
     2    +rho2p2(i)*gcp1_2p2(i)

         temp2=rho20(i)*gcp2_20(i)+rho2m2(i)*gcp2_2m2(i)
     2    +rho2p2(i)*gcp2_2p2(i)

         temp11=rho20(i)*gcm1_20(i)+rho2m2(i)*gcm1_2m2(i)
     2    +rho2p2(i)*gcm1_2p2(i)

         temp21=rho20(i)*gcm2_20(i)+rho2m2(i)*gcm2_2m2(i)
     2    +rho2p2(i)*gcm2_2p2(i)


         gpcs(i,1)=sqrt(3.0/8.0)*s2th*temp1
         gpcs(i,2)=sqrt(3.0/8.0)*sth2*temp2

         gmcs(i,1)=-sqrt(3.0/8.0)*s2th*temp11
         gmcs(i,2)=sqrt(3.0/8.0)*sth2*temp21


          ENDDO
C----------------------END OF CSA CALCULATIONS ---------------

C******* ZERO-ORDER EFFECTIVE HAMILTONIAN *********************

         Zpp12=W_D11(1,2)
         Zmm12=W_D1(1,2)


C******************************************************************
C---BEGIN TIME EVOLUTION AS A FUNCTION OF SPINNING FREQUENCY--------

          DO II=1,NPOINTS

              TMIX=TMIX_INT+TMIX_INC*(II-1)



C**************************************************************
C    G and C -- COEFFICIENTS INVOLVING SINGLE SPIN OPERATORS
C--------------------------------------------------------------

          DO I=1,NS
           DO J=1,MF
	        GZ_P(I,J)=GPCS(I,J)
            GZ_M(I,J)=GMCS(I,J)
             GP_P(I,J)=GPCS(I,J)
            GP_M(I,J)=GMCS(I,J)
            GM_P(I,J)=GPCS(I,J)
            GM_M(I,J)=GMCS(I,J)



            ENDDO
          ENDDO
C--------------------------------------------------------------
          DO I=1,NS
           DO J=1,MF
            CZ_P(I,J)=GZ_P(I,J)/(J*WR)
            CZ_M(I,J)=-GZ_M(I,J)/(J*WR)
            CP_P(I,J)=GP_P(I,J)/((4.0*J*J-1.0)*WR)
            CP_M(I,J)=GP_M(I,J)/((4.0*J*J-1.0)*WR)
            CM_P(I,J)=GM_P(I,J)/((4.0*J*J-1.0)*WR)
            CM_M(I,J)=GM_M(I,J)/((4.0*J*J-1.0)*WR)



           ENDDO
          ENDDO

C****************************************************************
C------G AND C  COEFFICIENTS FOR DIPOLAR INTERACTIONS -----------

C----------- PAIR 1-2 ------------------------------------------

C        Gzz12_P(1)=2.0*W_D1(1,2)
C        Gzz12_P(2)=2.0*W_D2(1,2)
C        Gzz12_M(1)=2.0*W_D11(1,2)
C        Gzz12_M(2)=2.0*W_D21(1,2)

C        Gpm12_P(2)=-0.5*W_D1(1,2)
C        Gpm12_P(3)=-0.5*W_D2(1,2)
C        Gpm12_M(1)=-0.5*W_D21(1,2)

C        Gmp12_M(2)=-0.5*W_D11(1,2)
C        Gmp12_M(3)=-0.5*W_D21(1,2)
C        Gmp12_P(1)=-0.5*W_D2(1,2)


        Gpp12_P(2)=W_D1(1,2)
        Gpp12_P(3)=W_D2(1,2)
        Gpp12_M(1)=W_D21(1,2)

        Gmm12_M(2)=W_D11(1,2)
        Gmm12_M(3)=W_D21(1,2)
        Gmm12_P(1)=W_D2(1,2)



C        DO I=1,MF
C         Czz12_P(I)=Gzz12_P(I)/(I*WR)
C         Czz12_M(I)=-Gzz12_M(I)/(I*WR)

C        Cpm12_P(I)=Gpm12_P(I)/(-I*WR-WRF1+WRF2)
C         Cpm12_M(I)=Gpm12_M(I)/(I*WR-WRF1+WRF2)

C         Cmp12_P(I)=Gmp12_P(I)/(-I*WR+WRF1-WRF2)
C         Cmp12_M(I)=Gmp12_M(I)/(I*WR+WRF1-WRF2)

C         ENDDO

         DO I=1,MF

         Cpp12_P(I)=Gpp12_P(I)/(I*WR)
         Cpp12_M(I)=Gpp12_M(I)/(-I*WR)

         Cmm12_P(I)=Gmm12_P(I)/(I*WR)
         Cmm12_M(I)=Gmm12_M(I)/(-I*WR)

        ENDDO


C***** CALCULATIONS OF G AND C COEFF COMPLETE **********************

C        WRITE(*,*)r(1,2),r(1,3),r(2,3)

C--------CALCULATION OF SECOND ORDER TERMS----------------------------
C*********************************************************************

C--- CROSS TERMS BETWEEN CSA (OF SPINS 1 AND 2) WITH DIPOLAR (1-2) ----

        SM4A=CMPLX(0.0,0.0)
        SM4B=CMPLX(0.0,0.0)


        DO I=2,MF
            DO J=1,MF
      SM4A=((Cpp12_P(I)*Gmm12_M(I)+Cpp12_M(1)*Gmm12_P(1))-(W_D11(1,2)*
     1 (Cmm12_M(I)+Cmm12_P(1)))-(Cmm12_P(1)*Gpp12_M(1)+Cmm12_M(I)*
     2  Gpp12_P(I))-(W_D1(1,2)*(Cpp12_P(I)+Cpp12_M(1))))+SM4A

        SM4B=((Cpp12_P(I)*Gmm12_M(I)+Cpp12_M(1)*Gmm12_P(1))-(W_D11(1,2)*
     1 (Cmm12_M(I)+Cmm12_P(1)))-(Cmm12_P(1)*Gpp12_M(1)+Cmm12_M(I)*
     2  Gpp12_P(I))-(W_D1(1,2)*(Cpp12_P(I)+Cpp12_M(1))))+SM4B

         ENDDO
         ENDDO


        SM3A=CMPLX(0.0,0.0)
        SM3B=CMPLX(0.0,0.0)


        DO I=1,MF
      SM3A=CP_P(1,I)*GP_M(1,I)+CP_M(1,I)*GP_P(1,I)+SM3A

        SM3B=CP_P(2,I)*GP_M(2,I)+CP_M(2,I)*GP_P(2,I)+SM3B

         ENDDO
C----------- END OF SECOND ORDER CORRECTIONS --------------------------


        Epp12=(3.0/4.0)*Zpp12
        Emm12=(3.0/4.0)*Zmm12

C           W1=WRF1+SM4A-AN*WR
C           W2=WRF2+SM4B-AN*WR

          W1=(9.0/16.0)*SM4A-SM3A
          W2=(9.0/16.0)*SM4B-SM3B

        E12=ABS(Epp12)
        E21=ABS(Emm12)
        E12R=E12*E21
        Z12=(W1+W2)*0.5
        Z12R=Z12*Z12
        X12=SQRT(E12R+Z12R)
        Y12=SQRT(E12R+Z12R)
        X12R=X12*X12
        Y12R=E12R/X12R









C---BEGIN TIME EVOLUTION AS A FUNCTION OF SPINNING FREQUENCY--------

C          DO II=1,NPOINTS
C           TMIX=TMIX_INT+TMIX_INC*(II-1)

C-----------------CALCULATES THE OBSERVABLES-----------------------

        ARES1(II)=(1.0-Y12R*SIN(2.0*PI*Y12*TMIX)*
     1       SIN(2.0*PI*Y12*TMIX))*AREA(IANGLE)+ARES1(II)

        ARES2(II)=(-Y12R*SIN(2.0*PI*Y12*TMIX)*
     1       SIN(2.0*PI*Y12*TMIX))*AREA(IANGLE)+ARES2(II)



C-------------------------------------------------------------------

                 ENDDO

C-------------------END OF LOOP FOR MIXING TIME--------------------

                 ENDDO

C-------------------END OF LOOP FOR POWDER AVERAGING---------------


C-----------------END OF CALCULATIONS-----------------------------


                DO IK=1,NPOINTS

                RES1(IK)=REAL(ARES1(IK))
	     	RES2(IK)=REAL(ARES2(IK))

                TMIX=TMIX_INT+TMIX_INC*(IK-1)
c                WR=WR_INT+WR_INC*(IK-1)
                WRITE(15,*)TMIX,RES1(IK),RES2(IK)


                ENDDO

C-----------------------------------------------------------

         CLOSE(11)
         CLOSE(12)
         CLOSE(15)

          END

C-----------------------------------------------------------
C***********************************************************

C***********************************************************


         subroutine Wigner(alfa,beita,gm_pm,Wig_mat)

         complex Wig_mat(5,5)
         real*8 d_red(5,5),bt,alfa,beita,gm_pm

               do i=1,5
                 do j=1,5
                   d_red(i,j)=0.0
                  enddo
                enddo

             do i=1,5
               do j=1,5
                wig_mat(i,j)=cmplx(0.0,0.0)
               enddo
              enddo


        bt=beita

         d_red(1,1)=(cos(bt/2.0))**4
         d_red(5,5)=d_red(1,1)
         d_red(2,1)=-0.50*(1.0+cos(bt))*sin(bt)
         d_red(1,2)=-d_red(2,1)
         d_red(5,4)=d_red(2,1)
         d_red(4,5)=-d_red(2,1)
         d_red(3,1)=sqrt(3.0/8.0)*(sin(bt))**2
         d_red(1,3)=d_red(3,1)
         d_red(3,5)=d_red(3,1)
         d_red(5,3)=d_red(3,1)
         d_red(4,1)=0.50*(cos(bt)-1.0)*sin(bt)
         d_red(1,4)=-d_red(4,1)
         d_red(2,5)=d_red(1,4)
         d_red(5,2)=-d_red(1,4)
         d_red(5,1)=(sin(bt/2.0))**4
         d_red(1,5)=d_red(5,1)
         d_red(2,2)=0.50*(2.0*cos(bt)-1.0)
     1   *(1.0+cos(bt))
         d_red(4,4)=d_red(2,2)
         d_red(3,2)=-sqrt(3.0/8.0)*(sin(2.0*bt))
         d_red(2,3)=-d_red(3,2)
         d_red(3,4)=-d_red(3,2)
         d_red(4,3)=d_red(3,2)
         d_red(3,3)=0.50*(3.0*(cos(bt))**2-1.0)
         d_red(4,2)=0.50*(2.0*cos(bt)+1.0)
     1   *(1.0-cos(bt))
         d_red(2,4)=d_red(4,2)


         do ii=1,5
          do jj=1,5
         Wig_mat(ii,jj)=d_red(ii,jj)
     1   *exp((0.0,1.0)*(jj-3.0)*gm_pm*(-1.0))
     2   *exp((0.0,1.0)*(ii-3.0)*alfa*(-1.0))

           enddo
         enddo

         end
C****************************************************************

