C     path:      %P%
C     revision:  $Revision$
C     created:   $Date$  
C     presently: %H%  %T%
      SUBROUTINE RTR

C *** This program calculates the upward fluxes, downward fluxes,
C     and heating rates for an arbitrary atmosphere.  The input to
C     this program is the atmospheric profile and all Planck function
C     information.  First-order "numerical" quadrature is used for the 
C     angle integration, i.e. only one exponential is computed per layer
C     per g-value per band.

      PARAMETER (MG=16)
      PARAMETER (MXLAY=203)
      PARAMETER (MXANG = 4)
      PARAMETER (NBANDS = 16)

      IMPLICIT DOUBLE PRECISION (V)                                     

      COMMON /CONSTANTS/ PI,FLUXFAC,HEATFAC
      COMMON /FEATURES/  NG(NBANDS),NSPA(MG),NSPB(MG)
      COMMON /BANDS/     WAVENUM1(NBANDS),WAVENUM2(NBANDS),
     &                   DELWAVE(NBANDS)
      COMMON /CONTROL/   NUMANGS, IOUT, ISTART, IEND
      COMMON /PROFILE/   NLAYERS,PAVEL(MXLAY),TAVEL(MXLAY),
     &                   PZ(0:MXLAY),TZ(0:MXLAY),TBOUND
      COMMON /PLNKDAT/   PLANKLAY(MXLAY,NBANDS),
     &                   PLANKLEV(0:MXLAY,NBANDS),PLANKBND(NBANDS)
      COMMON /PLANKG/    FRACS(MXLAY,MG)
      COMMON /TAUGCOM/   TAUG(MXLAY,MG)
      COMMON /OUTPUT/    TOTUFLUX(0:MXLAY), TOTDFLUX(0:MXLAY),
     &                   FNET(0:MXLAY), HTR(0:MXLAY)
      COMMON /HVERSN/    HVRRTM,HVRREG,HVRRTR,HVRATM,HVRSET,HVRTAU,
     *                   HVDUM1(4),HVRUTL,HVREXT
      COMMON /HVRSNB/    HVRKG(NBANDS)

      CHARACTER*8 HVRRTM,HVRREG,HVRRTR,HVRATM,HVRSET,HVRTAU,
     *            HVDUM1,HVRUTL,HVREXT
      CHARACTER*8 HVRKG
                                       
      DIMENSION BBD1(MXLAY),BBD2(MXLAY),BBD3(MXLAY)
      DIMENSION ATRANS1(MXLAY),ATRANS2(MXLAY),ATRANS3(MXLAY)
      DIMENSION UFLUX(0:MXLAY),DFLUX(0:MXLAY)
      DIMENSION DRAD1(0:MXLAY),URAD1(0:MXLAY)
      DIMENSION DRAD2(0:MXLAY),URAD2(0:MXLAY)
      DIMENSION DRAD3(0:MXLAY),URAD3(0:MXLAY)
      DIMENSION WTNUM(MXANG)

C *** These weights correspond to angles of 34.9, 65.8, and 74.1
C     degrees, resp.  They are used when "numerical" Gaussian
C     quadrature is chosen.
      DATA WTNUM(3) /0.1084176674/, WTNUM(2) /0.0424353369/
      DATA WTNUM(1) /0.3491473794/

      HVRRTR = '$Revision$'

C *** SECANG is equal to the secant of the first angle.
      SECANG = 1.219512195

      DO 200 LAY = 0, NLAYERS
         URAD1(LAY) = 0.0
         URAD2(LAY) = 0.0
         URAD3(LAY) = 0.0
         DRAD1(LAY) = 0.0
         DRAD2(LAY) = 0.0
         DRAD3(LAY) = 0.0
         TOTUFLUX(LAY) = 0.0
         TOTDFLUX(LAY) = 0.0
 200  CONTINUE

C *** Loop over frequency bands.
      DO 6000 IBAND = ISTART, IEND
         IF (IBAND .EQ. 1) THEN
            CALL TAUGB1
         ELSEIF (IBAND .EQ. 2) THEN
            CALL TAUGB2
         ELSEIF (IBAND .EQ. 3) THEN
            CALL TAUGB3
         ELSEIF (IBAND .EQ. 4) THEN
            CALL TAUGB4
         ELSEIF (IBAND .EQ. 5) THEN
            CALL TAUGB5
         ELSEIF (IBAND .EQ. 6) THEN
            CALL TAUGB6
         ELSEIF (IBAND .EQ. 7) THEN
            CALL TAUGB7
         ELSEIF (IBAND .EQ. 8) THEN
            CALL TAUGB8
         ELSEIF (IBAND .EQ. 9) THEN
            CALL TAUGB9
         ELSEIF (IBAND .EQ. 10) THEN
            CALL TAUGB10
         ELSEIF (IBAND .EQ. 11) THEN
            CALL TAUGB11
         ELSEIF (IBAND .EQ. 12) THEN
            CALL TAUGB12
         ELSEIF (IBAND .EQ. 13) THEN
            CALL TAUGB13
         ELSEIF (IBAND .EQ. 14) THEN
            CALL TAUGB14
         ELSEIF (IBAND .EQ. 15) THEN
            CALL TAUGB15
         ELSE
            CALL TAUGB16
         ENDIF
        
C ***    Loop over g-channels.
         IG = 1
 1000    CONTINUE

C ***    Radiative transfer starts here.
         RADLU1 = FRACS(1,IG) * PLANKBND(IBAND)
         RADLU2 = RADLU1
         RADLU3 = RADLU1
         URAD1(0) = URAD1(0) + RADLU1
         RADLD1 = 0.
         RADLD2 = 0.
         RADLD3 = 0.
C ***    Upward radiative transfer.  Due to the simple form taken by
C        certain equations when the optical depth is small, this 
C        condition is tested for.  In either case, the labels 1, 2, and
C        3 refer to the three angles.
         BGLEV = FRACS(1,IG) * PLANKLEV(0,IBAND)
         DO 2500 LEV = 1, NLAYERS
            BLAY = PLANKLAY(LEV,IBAND)
            PLFRAC = FRACS(LEV,IG)
            BGLAY = PLFRAC * BLAY
            ODEPTH = SECANG * TAUG(LEV,IG)
            IF (ODEPTH .LE. 1.E-2) THEN
               ATRANS1(LEV) = ODEPTH
               ATRANS2(LEV) = ODEPTH + ODEPTH
               ATRANS3(LEV) = ATRANS2(LEV) + ODEPTH
               BBD1(LEV) = BGLAY
               BBD2(LEV) = BGLAY
               BBD3(LEV) = BGLAY
               RADLU1 = RADLU1 + (BGLAY-RADLU1)*ATRANS1(LEV)
               RADLU2 = RADLU2 + (BGLAY-RADLU2)*ATRANS2(LEV)
               RADLU3 = RADLU3 + (BGLAY-RADLU3)*ATRANS3(LEV)
               BGLEV = PLFRAC * PLANKLEV(LEV,IBAND)
            ELSEIF (ODEPTH .GE. 10.) THEN
               ATRANS1(LEV) = 1.
               ATRANS2(LEV) = 1.
               ATRANS3(LEV) = 1.
               DELBGDN = BGLEV - BGLAY
               BGLEV = PLFRAC * PLANKLEV(LEV,IBAND)
               DELBGUP = BGLEV - BGLAY
               TAUSFAC = ODEPTH/(5.+ODEPTH)
               BBD1(LEV) = BGLAY + TAUSFAC * DELBGDN
               RADLU1 = BGLAY + TAUSFAC * DELBGUP
               TAUSFAC = ODEPTH/(2.5+ODEPTH)
               BBD2(LEV) = BGLAY + TAUSFAC * DELBGDN
               RADLU2 = BGLAY + TAUSFAC * DELBGUP
               TAUSFAC = ODEPTH/(1.666667+ODEPTH)
               BBD3(LEV) = BGLAY + TAUSFAC * DELBGDN
               RADLU3 = BGLAY + TAUSFAC * DELBGUP
            ELSE
               TRANS1 = EXP(-ODEPTH)
               ATRANS1(LEV) = 1. - TRANS1
               ATRANS2(LEV) = 1. - TRANS1*TRANS1
               ATRANS3(LEV) = ATRANS1(LEV) + TRANS1*ATRANS2(LEV)
               DELBGDN = BGLEV - BGLAY
               BGLEV = PLFRAC * PLANKLEV(LEV,IBAND)
               DELBGUP = BGLEV - BGLAY

C ***          TAUSFAC is needed for the Pade approximation for
C              the Planck function.
               TAUSFAC = ODEPTH/(5.+ODEPTH)
               BBU = BGLAY + TAUSFAC * DELBGUP
               RADLU1 = RADLU1 + (BBU-RADLU1)*ATRANS1(LEV)
               BBD1(LEV) = BGLAY + TAUSFAC * DELBGDN

               TAUSFAC = ODEPTH/(2.5+ODEPTH)
               BBU = BGLAY + TAUSFAC * DELBGUP
               RADLU2 = RADLU2 + (BBU-RADLU2)*ATRANS2(LEV)
               BBD2(LEV) = BGLAY + TAUSFAC * DELBGDN

               TAUSFAC = ODEPTH/(1.6666667+ODEPTH)
               BBU = BGLAY + TAUSFAC * DELBGUP
               RADLU3 = RADLU3 + (BBU-RADLU3)*ATRANS3(LEV)
               BBD3(LEV) = BGLAY + TAUSFAC * DELBGDN
            ENDIF

            URAD1(LEV) = URAD1(LEV) + RADLU1
            URAD2(LEV) = URAD2(LEV) + RADLU2
            URAD3(LEV) = URAD3(LEV) + RADLU3
 2500    CONTINUE

C ***    Downward radiative transfer.
         DO 2600 LEV = NLAYERS-1, 0, -1
            RADLD1 = RADLD1 + (BBD1(LEV+1)-RADLD1)*ATRANS1(LEV+1)
            RADLD2 = RADLD2 + (BBD2(LEV+1)-RADLD2)*ATRANS2(LEV+1)
            RADLD3 = RADLD3 + (BBD3(LEV+1)-RADLD3)*ATRANS3(LEV+1)
            DRAD1(LEV) = DRAD1(LEV) + RADLD1
            DRAD2(LEV) = DRAD2(LEV) + RADLD2
            DRAD3(LEV) = DRAD3(LEV) + RADLD3
 2600    CONTINUE

         IG = IG + 1
         IF (IG .LE. 16) GO TO 1000

C ***    First level radiance is independent of angle.
         URAD2(0) = URAD1(0)
         URAD3(0) = URAD1(0)
C ***    Calculate upward, downward, and net flux.
         DO 5000 LEV = NLAYERS, 0, -1
            UFLUX(LEV) = URAD1(LEV)*WTNUM(1) + URAD2(LEV)*WTNUM(2) + 
     &           URAD3(LEV)*WTNUM(3) 
            DFLUX(LEV) = DRAD1(LEV)*WTNUM(1) + DRAD2(LEV)*WTNUM(2) + 
     &           DRAD3(LEV)*WTNUM(3)
            URAD1(LEV) = 0.0
            URAD2(LEV) = 0.0
            URAD3(LEV) = 0.0
            DRAD1(LEV) = 0.0
            DRAD2(LEV) = 0.0
            DRAD3(LEV) = 0.0
            TOTUFLUX(LEV) = TOTUFLUX(LEV) + UFLUX(LEV) * DELWAVE(IBAND)
            TOTDFLUX(LEV) = TOTDFLUX(LEV) + DFLUX(LEV) * DELWAVE(IBAND)
 5000    CONTINUE
 6000 CONTINUE

      TOTUFLUX(0) = TOTUFLUX(0) * FLUXFAC
      TOTDFLUX(0) = TOTDFLUX(0) * FLUXFAC
      FNET(0) = TOTUFLUX(0) - TOTDFLUX(0)
      DO 7000 LEV = 1, NLAYERS
         TOTUFLUX(LEV) = TOTUFLUX(LEV) * FLUXFAC
         TOTDFLUX(LEV) = TOTDFLUX(LEV) * FLUXFAC
         FNET(LEV) = TOTUFLUX(LEV) - TOTDFLUX(LEV)
         L = LEV - 1

C        Calculate Heating Rates.
         HTR(L)=HEATFAC*(FNET(L)-FNET(LEV))/(PZ(L)-PZ(LEV)) 
 7000 CONTINUE

      HTR(NLAYERS) = 0.0

      RETURN
      END   
