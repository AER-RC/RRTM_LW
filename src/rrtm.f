****************************************************************************
*                                                                          *
*                               RRTM                                       *
*                                                                          *
*                                                                          *
*                                                                          *
*                   RAPID RADIATIVE TRANSFER MODEL                         *
*                                                                          *
*                                                                          *
*            ATMOSPHERIC AND ENVIRONMENTAL RESEARCH, INC.                  *
*                        840 MEMORIAL DRIVE                                *
*                        CAMBRIDGE, MA 02139                               *
*                                                                          *
*                                                                          *
*                           ELI J. MLAWER                                  *
*                         STEVEN J. TAUBMAN~                               *
*                         SHEPARD A. CLOUGH                                *
*                                                                          *
*                                                                          *
*                         ~currently at GFDL                               *
*                                                                          *
*                                                                          *
*                                                                          *
*                       email:  mlawer@aer.com                             *
*                                                                          *
*        The authors wish to acknowledge the contributions of the          *
*        following people:  Patrick D. Brown, Michael J. Iacono,           *
*        Ronald E. Farren, Luke Chen, Robert Bergstrom.                    *
*                                                                          *
****************************************************************************

      PROGRAM RRTM
                    
C *** This program is the driver for RRTM, the AER rapid model.  Before
C     the first atmosphere is processed, this routine calls INIT,
C     whick sets up any general information needed to calculate the   
C     fluxes and heating rates for any atmosphere.  Then 
C     for each atmosphere the user wishes to analyze, this routine
C     a) calls SETCOEF to calculate various quantities needed for 
C        the radiative transfer algorithm
C     b) calls RT or RTREG (depending on angular quadrature
C         method) to do the radiative transfer calculation
C     c) writes out the upward, downward, and net flux for each
C        level and the heating rate for each layer
C     Keep in mind that AER_RRTM is still a work in progress.

      PARAMETER (MXLAY=200)
      PARAMETER (MG = 16)
      PARAMETER (NBANDS = 16)

      COMMON /FEATURES/  NG(NBANDS),NSPA(MG),NSPB(MG)
      COMMON /CONTROL/   NUMANGS
      COMMON /PROFILE/   NLAYERS,PAVEL(MXLAY),TAVEL(MXLAY),
     &                   PZ(0:MXLAY),TZ(0:MXLAY)
      COMMON /OUTPUT/    TOTUFLUX(0:MXLAY), TOTDFLUX(0:MXLAY),
     &                   FNET(0:MXLAY), HTR(0:MXLAY)
      COMMON /CONSTANTS/ PI,FLUXFAC,HEATFAC

      DATA NG /16,16,16,16,16,16,16,16,16,16,16,16,16,16,16,16/
      DATA NSPA /1,1,7,7,7,1,0,0,0,0,0,0,0,0,0,0/
      DATA NSPB /1,1,3,3,3,0,0,0,0,0,0,0,0,0,0,0/

C     HEATFAC is the factor by which one must multiply delta-flux/ 
C     delta-pressure, with flux in w/m-2 and pressure in mbar, to get 
C     the heating rate in units of degrees/day.  It is equal to 
C           (g)x(#sec/day)x(1e-5)/(specific heat of air at const. p)
C        =  (9.8066)(3600)(1e-5)/(1.004)
      DATA HEATFAC /8.4391/

      PI = 2.*ASIN(1.)
      FLUXFAC = PI * 2.D4  

      ISETUP = 1
      IANGFLAG = 1
      NUMANGS = 1
      IWR = 10

      IF (ISETUP .EQ. 1) CALL INIT

      NUMATMOS = 1

C *** The input for AER_RRTM_RTREG is as follows:
C     NUMANGS - the number of angles over which the Gaussian quadrature
C               is performed.  This may be 1,2,3, or 4 for standard
C               quadrature and must be 3 for numerical quadrature.
C     NLAYERS - the number of atmospheric layers over which the
C               recursive radiative transfer algorithm is to be
C               performed.  This is one less than the number of output
C               values of upward flux, downward flux, and net flux,
C               and is equal to the number of output values of
C               heating rate.
C     Only NLAYERS is needed by RT.

C     Multiple atmospheres not yet implemented. 
      DO 4000 IATMOS = 1, NUMATMOS

C ***    Input layer information from INPUT_RRTM
         CALL READPROF

C ***    Calculate information needed by the radiative transfer routine
C        that is specific to this atmosphere, especially the coef-
C        ficients and indices needed to compute the optical depths
C        by interpolating data from stored reference atmospheres. 
         CALL SETCOEF

C ***    Call the radiative transfer routine.
         CALL RT 

C ***    Process output for this atmosphere.
         OPEN (IWR,FILE='OUTPUT_RRTM',FORM='FORMATTED')
         WRITE(IWR,9900)
         WRITE(IWR,9901)
C
         DO 3000 I = NLAYERS, 0, -1
            WRITE(IWR,9902) I, PZ(I), TOTUFLUX(I), TOTDFLUX(I),
     &           FNET(I), HTR(I)
 3000    CONTINUE

         CLOSE(IWR)

 4000 CONTINUE

 9900 FORMAT(1X,'LEVEL    PRESSURE   UPWARD FLUX   DOWNWARD FLUX    NET    
     &FLUX        HEATING RATE')
 9901 FORMAT(1X,'            mb          W/m2           W/m2          W/
     &m2           degree/day')
 9902 FORMAT(1X,I3,3X,F11.6,4X,1P,3(G12.6,3X),G16.9,0P)
      STOP
      END

C************************  SUBROUTINE READPROF  *****************************C

      SUBROUTINE READPROF                                                     
                                                                         
      IMPLICIT DOUBLE PRECISION (V)                                      
                                                                         
      PARAMETER (MXLAY=200)
      DIMENSION WBRODL(MXLAY), WKL(35,MXLAY), ALTZ(0:MXLAY)

      COMMON /PROFILE/ NLAYERS,PAVEL(MXLAY),TAVEL(MXLAY),
     &               PZ(0:MXLAY),TZ(0:MXLAY)
      COMMON /SPECIES/COLDRY(MXLAY),WATER(MXLAY),CO2(MXLAY),OZ(MXLAY)
 
   10 FORMAT (////)
   11 FORMAT (E10.3,E10.3)
   12 FORMAT (1X,I1,I3,I5)                                     
   13 FORMAT (G15.7,G10.4,G10.4,A3,I2,1X,2(G7.2,G8.3,G7.2))
   14 FORMAT (G15.7,G10.4,G10.4,A3,I2,23X,(G7.2,G8.3,G7.2))
   15 FORMAT (8G15.7)
                                                                         
      IRD = 9
      OPEN (IRD,FILE='INPUT_RRTM',FORM='FORMATTED')
      READ (IRD,10)
      READ (IRD,11) V1,V2
      READ (IRD,12) IFORM,NLAYERS,NMOL

      IF (NMOL.EQ.0) NMOL = 7                                          
                                                                         
      READ (IRD,13) PAVEL(1),TAVEL(1),SECNTK,CINP,IPTHAK,
     &              ALTZ(0),PZ(0),TZ(0),ALTZ(1),PZ(1),TZ(1)
      READ (IRD,15) (WKL(M,1),M=1,7), WBRODL(1)
      WATER(1) = WKL(1,1)
      CO2(1) = WKL(2,1)
      OZ(1) = WKL(3,1)
      COLDRY(1) = WBRODL(1)/(1.-CO2(1)-WKL(3,1)-WKL(4,1)-WKL(5,1)
     &     -WKL(6,1)-WKL(7,1))
      IF(NMOL .GT. 7) READ (IRD,15) (WKL(M,1),M=8,NMOL)

      DO 100 L = 2, NLAYERS                                                 
           READ (IRD,14) PAVEL(L),TAVEL(L),SECNTK,CINP,IPTHRK,
     &                   ALTZ(L),PZ(L),TZ(L)
           READ (IRD,15) (WKL(M,L),M=1,7), WBRODL(L)
           WATER(L) = WKL(1,L)
           CO2(L) = WKL(2,L)
           OZ(L) = WKL(3,L)
           COLDRY(L) = WBRODL(L)/(1.-CO2(L)-WKL(3,L)-WKL(4,L)
     &          -WKL(5,L)-WKL(6,L)-WKL(7,L))
C           WRITE(*,*)L,CO2(L),COLCO2(L),COLDRY
           IF(NMOL .GT. 7) READ (IRD,15) (WKL(M,L),M=8,NMOL)
  100 CONTINUE                                                            

      CLOSE(IRD)

      RETURN
      END 

