C     path:      %P%
C     revision:  $Revision$
C     created:   $Date$  
C     presently: %H%  %T%
      SUBROUTINE RT

C *** This program calculates the upward fluxes, downward fluxes,
C     and heating rates for an arbitrary atmosphere.  The input to
C     this program is the atmospheric profile and interpolation
C     coefficients that are used to calculate the optical depths
C     from values stored for reference atmospheric profiles. 
C     First-order "numerical" quadrature is used for the angle
C     integration, i.e. only one exponential is computed per layer
C     per g-value per band.

      PARAMETER (MG=16)
      PARAMETER (MXLAY=203)
      PARAMETER (MXANG = 4)
      PARAMETER (NBANDS = 16)

      IMPLICIT DOUBLE PRECISION (V)                                     

      COMMON /FEATURES/  NG(NBANDS),NSPA(MG),NSPB(MG)
      COMMON /CONTROL/   NUMANGS
      COMMON /PROFILE/   NLAYERS,PAVEL(MXLAY),TAVEL(MXLAY),
     &                   PZ(0:MXLAY),TZ(0:MXLAY),TSFC
      COMMON /BANDS/     WAVENUM1(NBANDS),WAVENUM2(NBANDS),
     &                   DELWAVE(NBANDS)
      COMMON /OUTPUT/    TOTUFLUX(0:MXLAY), TOTDFLUX(0:MXLAY),
     &                   FNET(0:MXLAY), HTR(0:MXLAY)
      COMMON /FRACS/     TLEVFRAC(0:MXLAY),TLAYFRAC(MXLAY)
      COMMON /MISC/      INDLEV(0:MXLAY),INDLAY(MXLAY)
      COMMON /CONSTANTS/ PI,FLUXFAC,HEATFAC
      COMMON /TAUGCOM/   TAUG(MXLAY,MG)
      COMMON /SOURCE/    AVGPLANK(16,141,16)

      COMMON /HVERSN/ HVRRTM,HVRINI,HVRRT0,HVRATM,HVRSET,HVRTAU,
     *                HVDUM1(4),HVRUTL,HVREXT

      EQUIVALENCE (AVGPLANK,PLANK)

      DIMENSION BAVE(MXLAY,16), PLANK(2256,16)
      DIMENSION BBD1(MXLAY),BBD2(MXLAY),BBD3(MXLAY)
      DIMENSION DELTUP(MXLAY),DELTDN(MXLAY),ILAY0(MXLAY)
      DIMENSION DIFFD(MXLAY,MG),DIFFU(MXLAY,MG), BBND0(MG)
      DIMENSION ATRANS1(MXLAY),ATRANS2(MXLAY),ATRANS3(MXLAY)
      DIMENSION UFLUX(0:MXLAY),DFLUX(0:MXLAY)
      DIMENSION DRAD1(0:MXLAY-1),URAD1(0:MXLAY)
      DIMENSION DRAD2(0:MXLAY-1),URAD2(0:MXLAY)
      DIMENSION DRAD3(0:MXLAY-1),URAD3(0:MXLAY)
      DIMENSION WTNUM(MXANG)

      CHARACTER*8 HVRRTM,HVRINI,HVRRT0,HVRATM,HVRSET,HVRTAU,
     *            HVDUM1,HVRUTL,HVREXT

C *** These weights correspond to angles of 34.9, 65.8, and 74.1
C     degrees, resp.  They are used when "numerical" Gaussian
C     quadrature is chosen.
      DATA WTNUM(3) /0.1084176674/, WTNUM(2) /0.0424353369/
      DATA WTNUM(1) /0.3491473794/

      HVRRT0 = '$Revision$'

C *** SECANG is equal to the secant of the first angle.
      SECANG = 1.219512195
      
      TOTUFLUX(0) = 0.0
      TOTDFLUX(0) = 0.0
      DO 200 LAY = 1, NLAYERS
         TOTUFLUX(LAY) = 0.0
         TOTDFLUX(LAY) = 0.0
         ILAY0(LAY) = (INDLAY(LAY) - 1) * 16
         DELTUP(LAY) = TZ(LAY) - TAVEL(LAY)
         DELTDN(LAY) = TZ(LAY-1) - TAVEL(LAY)
 200  CONTINUE
      INDEX0 = (INDLEV(0) - 1) * 16

C *** Loop over frequency bands.
 500  CONTINUE
      DO 6000 IBAND = 1, 2

C *** Compute the Planck function at the surface for the frequency
C     chosen as the representative frequency for each IG.
         INDEX = INDEX0
         DO 600 IG = 1, 16
            INDEX = INDEX + 1
            BBND0(IG) = PLANK(INDEX,IBAND) + TLEVFRAC(0) *
     &           (PLANK(INDEX+16,IBAND) - PLANK(INDEX,IBAND))
 600     CONTINUE
         
         DO 800 LEV = 0, NLAYERS-1
            LAY = LEV + 1
C ***       Compute the Planck function for each layer and the upward
C           and downward change in Planck function between layer "center"
C           and respective layer boundary.  These computations are done
C           at the frequency chosen as the representative frequency for 
C           each IG.
            ILAY = ILAY0(LAY)
            DO 700 IG = 1, 16
               ILAY = ILAY + 1
               DPLNKDT = PLANK(ILAY+16,IBAND) - PLANK(ILAY,IBAND)
               BAVE(LAY,IG) = PLANK(ILAY,IBAND) + TLAYFRAC(LAY)*DPLNKDT
               DIFFU(LAY,IG) = DELTUP(LAY)*DPLNKDT 
               DIFFD(LAY,IG) = DELTDN(LAY)*DPLNKDT
 700        CONTINUE
 800     CONTINUE

         DRAD1(NLAYERS) = 0.
         DRAD2(NLAYERS) = 0.
         DRAD3(NLAYERS) = 0.

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
         ENDIF
        
C ***    Loop over g-channels.
         IG = 1
 1000    CONTINUE

C ***    Radiative transfer starts here.
         RADLU1 = BBND0(IG)
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
         DO 2500 LEV = 1, NLAYERS
            ODEPTH = SECANG * TAUG(LEV,IG)
            IF (ODEPTH .LE. 1.E-2) THEN
               ATRANS1(LEV) = ODEPTH
               ATRANS2(LEV) = ODEPTH+ODEPTH
               ATRANS3(LEV) = ATRANS2(LEV)+ODEPTH
               BBU = BAVE(LEV,IG)
               BBD1(LEV) = BBU
               BBD2(LEV) = BBU
               BBD3(LEV) = BBU
               RADLU1 = RADLU1 + (BBU-RADLU1)*ATRANS1(LEV)
               RADLU2 = RADLU2 + (BBU-RADLU2)*ATRANS2(LEV)
               RADLU3 = RADLU3 + (BBU-RADLU3)*ATRANS3(LEV)
            ELSE
               TRANS1 = EXP(-ODEPTH)
               ATRANS1(LEV) = 1.-TRANS1
               ATRANS2(LEV) = 1.-TRANS1*TRANS1
               ATRANS3(LEV) = ATRANS1(LEV) + TRANS1*ATRANS2(LEV)

C ***          TAUSFAC is needed for the Pade approximation for
C              the Planck function.
               TAUSFAC = ODEPTH/(5.+ODEPTH)
               BBU = BAVE(LEV,IG) + TAUSFAC*DIFFU(LEV,IG)
               RADLU1 = RADLU1 + (BBU-RADLU1)*ATRANS1(LEV)
               BBD1(LEV) = BAVE(LEV,IG) + TAUSFAC*DIFFD(LEV,IG)

               TAUSFAC = ODEPTH/(2.5+ODEPTH)
               BBU = BAVE(LEV,IG) + TAUSFAC*DIFFU(LEV,IG)
               RADLU2 = RADLU2 + (BBU-RADLU2)*ATRANS2(LEV)
               BBD2(LEV) = BAVE(LEV,IG) + TAUSFAC*DIFFD(LEV,IG)

               TAUSFAC = ODEPTH/(1.6666667+ODEPTH)
               BBU = BAVE(LEV,IG) + TAUSFAC*DIFFU(LEV,IG)
               RADLU3 = RADLU3 + (BBU-RADLU3)*ATRANS3(LEV)
               BBD3(LEV) = BAVE(LEV,IG) + TAUSFAC*DIFFD(LEV,IG)
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
      DO 7000 LEV = 0, NLAYERS
         TOTUFLUX(LEV) = TOTUFLUX(LEV) * FLUXFAC
         TOTDFLUX(LEV) = TOTDFLUX(LEV) * FLUXFAC
         FNET(LEV) = TOTUFLUX(LEV) - TOTDFLUX(LEV)
 7000 CONTINUE

C *** Calculate Heating Rates.

      HTR(NLAYERS) = 0.0
      DO 8000 LEV  = NLAYERS-1, 0, -1
         LAY = LEV + 1
         HTR(LEV)=HEATFAC*(FNET(LEV)-FNET(LAY))/(PZ(LEV)-PZ(LAY)) 
 8000 CONTINUE

      RETURN
      END   
