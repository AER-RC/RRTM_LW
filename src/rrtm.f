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
                    
C *** This program is the driver for RRTM, the AER rapid model.  
C     For each atmosphere the user wishes to analyze, this routine
C     a) calls READPROF to read in the atmospheric profile
C     b) calls SETCOEF to calculate various quantities needed for 
C        the radiative transfer algorithm
C     c) calls RT or RTREG (depending on angular quadrature
C         method) to do the radiative transfer calculation
C     d) writes out the upward, downward, and net flux for each
C        level and the heating rate for each layer

      PARAMETER (MXLAY=203)
      PARAMETER (MG = 16)
      PARAMETER (NBANDS = 16)

      COMMON /CONSTANTS/ PI,FLUXFAC,HEATFAC
      COMMON /FEATURES/  NG(NBANDS),NSPA(MG),NSPB(MG)
      COMMON /BANDS/     WAVENUM1(NBANDS),WAVENUM2(NBANDS),
     &                   DELWAVE(NBANDS)
      COMMON /CONTROL/   NUMANGS, IOUT, ISTART, IEND
      COMMON /PROFILE/   NLAYERS,PAVEL(MXLAY),TAVEL(MXLAY),
     &                   PZ(0:MXLAY),TZ(0:MXLAY),TBOUND
      COMMON /OUTPUT/    TOTUFLUX(0:MXLAY), TOTDFLUX(0:MXLAY),
     &                   FNET(0:MXLAY), HTR(0:MXLAY)

      DATA WAVENUM1(1) /10./, WAVENUM2(1) /250./, DELWAVE(1) /240./
      DATA WAVENUM1(2) /250./, WAVENUM2(2) /500./, DELWAVE(2) /250./
      DATA WAVENUM1(3) /500./, WAVENUM2(3) /630./, DELWAVE(3) /130./
      DATA WAVENUM1(4) /630./, WAVENUM2(4) /700./, DELWAVE(4) /70./
      DATA WAVENUM1(5) /700./, WAVENUM2(5) /820./, DELWAVE(5) /120./
      DATA WAVENUM1(6) /820./, WAVENUM2(6) /980./, DELWAVE(6) /160./
      DATA WAVENUM1(7) /980./, WAVENUM2(7) /1080./, DELWAVE(7) /100./
      DATA WAVENUM1(8) /1080./, WAVENUM2(8) /1180./, DELWAVE(8) /100./
      DATA WAVENUM1(9) /1180./, WAVENUM2(9) /1390./, DELWAVE(9) /210./
      DATA WAVENUM1(10) /1390./,WAVENUM2(10) /1480./,DELWAVE(10) /90./
      DATA WAVENUM1(11) /1480./,WAVENUM2(11) /1800./,DELWAVE(11) /320./
      DATA WAVENUM1(12) /1800./,WAVENUM2(12) /2080./,DELWAVE(12) /280./
      DATA WAVENUM1(13) /2080./,WAVENUM2(13) /2250./,DELWAVE(13) /170./
      DATA WAVENUM1(14) /2250./,WAVENUM2(14) /2380./,DELWAVE(14) /130./
      DATA WAVENUM1(15) /2380./,WAVENUM2(15) /2600./,DELWAVE(15) /220./
      DATA WAVENUM1(16) /2600./,WAVENUM2(16) /3000./,DELWAVE(16) /400./

      DATA NG /16,16,16,16,16,16,16,16,16,16,16,16,16,16,16,16/
      DATA NSPA /1,1,7,7,7,1,9,1,11,1,1,9,9,1,9,9/
      DATA NSPB /1,1,3,3,3,0,1,1, 1,1,1,0,0,1,0,0/

C     HEATFAC is the factor by which one must multiply delta-flux/ 
C     delta-pressure, with flux in w/m-2 and pressure in mbar, to get 
C     the heating rate in units of degrees/day.  It is equal to 
C           (g)x(#sec/day)x(1e-5)/(specific heat of air at const. p)
C        =  (9.8066)(3600)(1e-5)/(1.004)
      DATA HEATFAC /8.4391/
      CHARACTER PAGE

      PI = 2.*ASIN(1.)
      FLUXFAC = PI * 2.D4  

      IWR = 10
      PAGE = CHAR(12)

C     Multiple atmospheres not yet implemented. 
      NUMATMOS = 1
      DO 4000 IATMOS = 1, NUMATMOS

C ***    Input atmospheric profile from INPUT_RRTM.
         CALL READPROF

         ISTART = 1
         IEND = 16
         IFLAG = IOUT

 1000    CONTINUE
         IF (IFLAG .GT. 0 .AND. IFLAG .LE. 16) THEN
            ISTART = IFLAG
            IEND = IFLAG
         ENDIF

C ***    Calculate information needed by the radiative transfer routine
C        that is specific to this atmosphere, especially some of the 
C        coefficients and indices needed to compute the optical depths
C        by interpolating data from stored reference atmospheres. 

         CALL SETCOEF

C ***    Call the radiative transfer routine.
         IF (NUMANGS .EQ. 0) THEN
            CALL RTR
         ELSE
            CALL RTREG
         ENDIF

         IF (IOUT .LT. 0) GO TO 4000

C ***    Process output for this atmosphere.
         OPEN (IWR,FILE='OUTPUT_RRTM',FORM='FORMATTED')
         WRITE(IWR,9899)WAVENUM1(ISTART),WAVENUM2(IEND)
         WRITE(IWR,9900)
         WRITE(IWR,9901)
C
         DO 3000 I = NLAYERS, 0, -1
            WRITE(IWR,9902) I, PZ(I), TOTUFLUX(I), TOTDFLUX(I),
     &           FNET(I), HTR(I)
 3000    CONTINUE
         WRITE(IWR,9903)PAGE
         
         IF (IOUT .LE. 16 .OR. IFLAG .EQ. 16) GO TO 3500
         IF (IFLAG .EQ. 99) THEN
            IFLAG = 1
         ELSEIF (IOUT .EQ. 99) THEN
            IFLAG = IFLAG + 1
         ENDIF
         GO TO 1000

 3500    CONTINUE
         CLOSE(IWR)

 4000 CONTINUE

 9899 FORMAT(1X,'Wavenumbers: ',F6.1,' - ',F6.1,' cm-1')
 9900 FORMAT(1X,'LEVEL    PRESSURE   UPWARD FLUX   DOWNWARD FLUX    NET    
     &FLUX       HEATING RATE')
 9901 FORMAT(1X,'            mb          W/m2          W/m2           W/
     &m2          degree/day')
 9902 FORMAT(1X,I3,3X,F11.6,4X,1P,2(G12.6,2X),G13.6,3X,G16.9,0P)
 9903 FORMAT(A)

      STOP
      END

C************************  SUBROUTINE READPROF  *****************************C

      SUBROUTINE READPROF                                                     
                                                                         
C     Read in atmospheric profile.

      IMPLICIT DOUBLE PRECISION (V)                                      
                                                                         
      PARAMETER (MXLAY=203)
      DIMENSION ALTZ(0:MXLAY), SUMMOL(MXLAY)


      COMMON /CONTROL/ NUMANGS, IOUT, ISTART, IEND
      COMMON /PROFILE/ NLAYERS,PAVEL(MXLAY),TAVEL(MXLAY),
     &                 PZ(0:MXLAY),TZ(0:MXLAY),TBOUND
      COMMON /SPECIES/ COLDRY(MXLAY),WKL(35,MXLAY),WBRODL(MXLAY),
     &                 NMOL

      CHARACTER*80 FORM1(0:1),FORM2(0:1),FORM3(0:1)
      CHARACTER*1 CTEST, CDOLLAR

      DATA CDOLLAR /'$'/

      FORM1(0) = '(3F10.4,A3,I2,1X,2(F7.2,F8.3,F7.2))'
      FORM2(0) = '(3F10.4,A3,I2,23X,(F7.2,F8.3,F7.2))'
      FORM3(0) = '(8E10.3)'
      FORM1(1) = '(G15.7,G10.4,G10.4,A3,I2,1X,2(G7.2,G8.3,G7.2))'
      FORM2(1) = '(G15.7,G10.4,G10.4,A3,I2,23X,(G7.2,G8.3,G7.2))'
      FORM3(1) = '(8G15.7)'

      IRD = 9
      OPEN (IRD,FILE='INPUT_RRTM',FORM='FORMATTED')

 1000 CONTINUE
      READ (IRD,9010,END=8800) CTEST
      IF (CTEST .NE. CDOLLAR) GO TO 1000

      READ(IRD,9011) IATM, NUMANGS, IOUT
      READ (IRD,9012) TBOUND

      IF (IATM .EQ. 0) THEN
         READ (IRD,9013) IFORM,NLAYERS,NMOL
         IF (NMOL.EQ.0) NMOL = 7                                    
         READ (IRD,FORM1(IFORM)) PAVEL(1),TAVEL(1),SECNTK,CINP,
     &        IPTHAK,ALTZ(0),PZ(0),TZ(0),ALTZ(1),PZ(1),TZ(1)
         READ (IRD,FORM3(IFORM)) (WKL(M,1),M=1,7), WBRODL(1)
         IF(NMOL .GT. 7) READ (IRD,FORM3(IFORM)) (WKL(M,1),M=8,NMOL)

         DO 4000 L = 2, NLAYERS
            READ (IRD,FORM2(IFORM)) PAVEL(L),TAVEL(L),SECNTK,CINP,
     &           IPTHRK,ALTZ(L),PZ(L),TZ(L)
            READ (IRD,FORM3(IFORM)) (WKL(M,L),M=1,7), WBRODL(L)
            IF(NMOL .GT. 7) READ (IRD,FORM3(IFORM)) (WKL(M,L),M=8,NMOL)
 4000    CONTINUE                                                            
           
      ELSE
         CALL RRTATM
      ENDIF

C     Test for mixing ratio input.
      IMIX = 0
      IF (WKL(2,1) .LE. 1.0) IMIX = 1

      DO 5000 L = 1, NLAYERS
         SUMMOL(L) = 0.0
         DO 4100 IMOL = 2, NMOL
            SUMMOL(L) = SUMMOL(L) + WKL(IMOL,L)
 4100    CONTINUE
         IF (IMIX .EQ. 1) THEN
            COLDRY(L) = WBRODL(L) / (1. - SUMMOL(L))
            DO 4200 IMOL = 1, NMOL
               WKL(IMOL,L) = COLDRY(L) * WKL(IMOL,L)
 4200       CONTINUE
         ELSE
            COLDRY(L) = WBRODL(L) + SUMMOL(L)
         ENDIF
 5000 CONTINUE

      CLOSE(IRD)
      GO TO 9000

 8800 CONTINUE
      STOP ' INVALID INPUT_RRTM '

 9000 CONTINUE

 9010 FORMAT (A1)
 9011 FORMAT (49X,I1,34X,I1,2X,I3)
 9012 FORMAT (E10.3)
 9013 FORMAT (1X,I1,I3,I5)                                     

      RETURN
      END 

