C     path:      %P%
C     revision:  $Revision$
C     created:   $Date$  
C     presently: %H%  %T%
      SUBROUTINE RTRCLD

C *** This program calculates the upward fluxes, downward fluxes,
C     and heating rates for an arbitrary cloudy atmosphere.  The input 
C     to this program is the atmospheric profile, including cloud 
C     properties, and all Planck function information.  First-order 
C     "numerical" quadrature is used for the angle integration, i.e.
C     only one exponential is computed per layer per g-value per band.

      PARAMETER (MG=16)
      PARAMETER (MXLAY=203)
      PARAMETER (NBANDS = 16)
      PARAMETER (MXCBANDS = 5)

      IMPLICIT DOUBLE PRECISION (V)                                     

      COMMON /CONSTANTS/ PI,FLUXFAC,HEATFAC
      COMMON /FEATURES/  NG(NBANDS),NSPA(NBANDS),NSPB(NBANDS)
      COMMON /BANDS/     WAVENUM1(NBANDS),WAVENUM2(NBANDS),
     &                   DELWAVE(NBANDS)
      COMMON /CONTROL/   NUMANGS, IOUT, ISTART, IEND
      COMMON /PROFILE/   NLAYERS,PAVEL(MXLAY),TAVEL(MXLAY),
     &                   PZ(0:MXLAY),TZ(0:MXLAY),TBOUND
      COMMON /CLOUDDAT/  NCBANDS,CLDFRAC(MXLAY),TAUCLOUD(MXLAY,MXCBANDS)
      COMMON /PLNKDAT/   PLANKLAY(MXLAY,NBANDS),
     &                   PLANKLEV(0:MXLAY,NBANDS),PLANKBND(NBANDS)
      COMMON /PLANKG/    FRACS(MXLAY,MG)
      COMMON /TAUGCOM/   TAUG(MXLAY,MG)
      COMMON /OUTPUT/    TOTUFLUX(0:MXLAY), TOTDFLUX(0:MXLAY),
     &                   FNET(0:MXLAY), HTR(0:MXLAY)
      COMMON /HVERSN/    HVRRTM,HVRREG,HVRRTR,HVRATM,HVRSET,HVRTAU,
     *                   HVRRGC,HVRRTC,HVRCLD,HVRDUM,HVRUTL,HVREXT

      CHARACTER*8 HVRRTM,HVRREG,HVRRTR,HVRATM,HVRSET,HVRTAU,
     *            HVRRGC,HVRRTC,HVRCLD,HVRDUM,HVRUTL,HVREXT
                                       
      DIMENSION BBDGAS1(MXLAY),BBDGAS2(MXLAY),BBDGAS3(MXLAY)
      DIMENSION BBDTOT1(MXLAY),BBDTOT2(MXLAY),BBDTOT3(MXLAY)
      DIMENSION ATRANS1(MXLAY),ATRANS2(MXLAY),ATRANS3(MXLAY)
      DIMENSION UFLUX(0:MXLAY),DFLUX(0:MXLAY)
      DIMENSION DRAD1(0:MXLAY),URAD1(0:MXLAY)
      DIMENSION DRAD2(0:MXLAY),URAD2(0:MXLAY)
      DIMENSION DRAD3(0:MXLAY),URAD3(0:MXLAY)
      DIMENSION WTNUM(3), ODCLD(MXLAY,MXCBANDS)
      DIMENSION ATOT1(MXLAY),ATOT2(MXLAY),ATOT3(MXLAY)
      DIMENSION ABSCLD1(MXLAY,MXCBANDS),ABSCLD2(MXLAY,MXCBANDS)
      DIMENSION EFCLFRAC1(MXLAY,MXCBANDS),ABSCLD3(MXLAY,MXCBANDS)
      DIMENSION EFCLFRAC2(MXLAY,MXCBANDS),EFCLFRAC3(MXLAY,MXCBANDS)
      DIMENSION ICEBNDA(NBANDS),ICEBNDB(NBANDS)

C     These arrays indicate the spectral 'region' (used in the 
C     calculation of ice cloud optical depths) corresponding
C     to each spectral band.  See cldprop.f for more details.
      DATA ICEBNDA /1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1/
      DATA ICEBNDB /1,2,3,3,3,4,4,4,5,5,5,5,5,5,5,5/

C *** These weights correspond to angles of 34.9, 65.8, and 74.1
C     degrees, resp.  (They are used when "numerical" Gaussian
C     quadrature is chosen.) Note that, using these angles,
C     Trans_2 = (Trans_1) ** 2 and Trans_3 = (Trans_1) ** 3
      DATA WTNUM(3) /0.1084176674/, WTNUM(2) /0.0424353369/
      DATA WTNUM(1) /0.3491473794/

      HVRRTC = '$Revision$'

C *** SECANG is equal to the secant of the first angle.
      SECANG = 1.219512195

      URAD1(0) = 0.0
      URAD2(0) = 0.0
      URAD3(0) = 0.0
      DRAD1(0) = 0.0
      DRAD2(0) = 0.0
      DRAD3(0) = 0.0
      TOTUFLUX(0) = 0.0
      TOTDFLUX(0) = 0.0
      DO 200 LAY = 1, NLAYERS
         URAD1(LAY) = 0.0
         URAD2(LAY) = 0.0
         URAD3(LAY) = 0.0
         DRAD1(LAY) = 0.0
         DRAD2(LAY) = 0.0
         DRAD3(LAY) = 0.0
         TOTUFLUX(LAY) = 0.0
         TOTDFLUX(LAY) = 0.0
         DO 100 IB = 1, NCBANDS
            IF (CLDFRAC(LAY) .GE. 1.E-6) THEN
               ODCLD(LAY,IB) = SECANG * TAUCLOUD(LAY,IB)
               TRANSCLD = EXP(-ODCLD(LAY,IB))
               ABSCLD1(LAY,IB) = 1. - TRANSCLD
               ABSCLD2(LAY,IB) = 1. - TRANSCLD*TRANSCLD
               ABSCLD3(LAY,IB) = ABSCLD1(LAY,IB) + 
     &              TRANSCLD * ABSCLD2(LAY,IB)
               EFCLFRAC1(LAY,IB) = ABSCLD1(LAY,IB) * CLDFRAC(LAY)
               EFCLFRAC2(LAY,IB) = ABSCLD2(LAY,IB) * CLDFRAC(LAY)
               EFCLFRAC3(LAY,IB) = ABSCLD3(LAY,IB) * CLDFRAC(LAY)
            ELSE
               ODCLD(LAY,IB) = 0.0
               ABSCLD1(LAY,IB) = 0.0
               ABSCLD2(LAY,IB) = 0.0
               ABSCLD3(LAY,IB) = 0.0
               EFCLFRAC1(LAY,IB) = 0.0
               EFCLFRAC2(LAY,IB) = 0.0
               EFCLFRAC3(LAY,IB) = 0.0
            ENDIF
 100     CONTINUE
 200  CONTINUE

C *** Loop over frequency bands.
      DO 6000 IBAND = ISTART, IEND
         IF (NCBANDS .EQ. 1) THEN
            IB = ICEBNDA(IBAND)
         ELSE
            IB = ICEBNDB(IBAND)
         ENDIF
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
         ELSEIF (IBAND .EQ. 16) THEN
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
         URAD2(0) = URAD2(0) + RADLU1
         URAD3(0) = URAD3(0) + RADLU1
         RADLD1 = 0.
         RADLD2 = 0.
         RADLD3 = 0.
C ***    Upward radiative transfer loop.  Due to the simple form taken
C        by certain equations when the optical depth is small, this 
C        condition is tested for.  In all cases, the labels 1, 2, and
C        3 refer to the three angles.
         BGLEV = FRACS(1,IG) * PLANKLEV(0,IBAND)
         DO 2500 LEV = 1, NLAYERS
            BLAY = PLANKLAY(LEV,IBAND)
            PLFRAC = FRACS(LEV,IG)
            BGLAY = PLFRAC * BLAY
            ODEPTH = SECANG * TAUG(LEV,IG) 
            DELBGDN = BGLEV - BGLAY

C           Here are some variable definitions:
C             ODTOT      optical depth of gas and cloud
C             ATRANS     absorptivity for only gas
C             ATOT       absorptivity for gas and cloud
C             TFACGAS    gas-only Pade factor, used for Planck fn
C             TFACTOT    gas and cloud Pade factor, used for Planck fn
C             BBDGAS     gas-only Planck function for downward rt
C             BBDTOT     gas and cloud Planck function for downward rt
C             BBUTOT     gas and cloud Planck function for upward calc.
C             GASSRC     source radiance due to gas only

            IF (ODEPTH .LE. 1.E-2 .AND. ODCLD(LEV,IB) .LE. 1.E-2) THEN
               ATRANS1(LEV) = ODEPTH
               ATRANS2(LEV) = ODEPTH + ODEPTH
               ATRANS3(LEV) = ATRANS2(LEV) + ODEPTH
               BBDGAS1(LEV) = BGLAY
               BBDGAS2(LEV) = BGLAY
               BBDGAS3(LEV) = BGLAY
               BBDTOT1(LEV) = BGLAY
               BBDTOT2(LEV) = BGLAY
               BBDTOT3(LEV) = BGLAY
               BGLEV = PLFRAC * PLANKLEV(LEV,IBAND)
               BBUTOT1 = BGLAY
               BBUTOT2 = BGLAY
               BBUTOT3 = BGLAY
               ATOT1(LEV) = ATRANS1(LEV) + ABSCLD1(LEV,IB)
               ATOT2(LEV) = ATRANS2(LEV) + ABSCLD2(LEV,IB)
               ATOT3(LEV) = ATRANS3(LEV) + ABSCLD3(LEV,IB)
               GASSRC1 = BGLAY * ATRANS1(LEV)
               GASSRC2 = BGLAY * ATRANS2(LEV)
               GASSRC3 = BGLAY * ATRANS3(LEV)
            ELSEIF (ODEPTH .LE. 1.E-2) THEN
               ODTOT = ODEPTH + ODCLD(LEV,IB)
               ATRANS1(LEV) = ODEPTH
               ATRANS2(LEV) = ODEPTH + ODEPTH
               ATRANS3(LEV) = ATRANS2(LEV) + ODEPTH
               BBDGAS1(LEV) = BGLAY
               BBDGAS2(LEV) = BGLAY
               BBDGAS3(LEV) = BGLAY
               TFACTOT1 = ODTOT/(5.+ODTOT)
               BBDTOT1(LEV) = BGLAY + TFACTOT1*DELBGDN
               TFACTOT2 = ODTOT/(2.5+ODTOT)
               BBDTOT2(LEV) = BGLAY + TFACTOT2*DELBGDN
               TFACTOT3 = ODTOT/(1.666667+ODTOT)
               BBDTOT3(LEV) = BGLAY + TFACTOT3*DELBGDN
               BGLEV = PLFRAC * PLANKLEV(LEV,IBAND)
               DELBGUP = BGLEV - BGLAY
               BBUTOT1 = BGLAY + TFACTOT1*DELBGUP
               BBUTOT2 = BGLAY + TFACTOT2*DELBGUP
               BBUTOT3 = BGLAY + TFACTOT3*DELBGUP
               ATOT1(LEV) = ATRANS1(LEV) + ABSCLD1(LEV,IB)
               ATOT2(LEV) = ATRANS2(LEV) + ABSCLD2(LEV,IB)
               ATOT3(LEV) = ATRANS3(LEV) + ABSCLD3(LEV,IB)
               GASSRC1 = BGLAY * ATRANS1(LEV)
               GASSRC2 = BGLAY * ATRANS2(LEV)
               GASSRC3 = BGLAY * ATRANS3(LEV)
            ELSEIF (ODEPTH .GE. 10.) THEN
               ODTOT = ODEPTH + ODCLD(LEV,IB)
               ATRANS1(LEV) = 1.
               ATRANS2(LEV) = 1.
               ATRANS3(LEV) = 1.
               TFACGAS1 = ODEPTH/(5.+ODEPTH)
               BBDGAS1(LEV) = BGLAY + TFACGAS1*DELBGDN
               TFACGAS2 = ODEPTH/(2.5+ODEPTH)
               BBDGAS2(LEV) = BGLAY + TFACGAS2*DELBGDN
               TFACGAS3 = ODEPTH/(1.666667+ODEPTH)
               BBDGAS3(LEV) = BGLAY + TFACGAS3*DELBGDN
               TFACTOT1 = ODTOT/(5.+ODTOT)
               BBDTOT1(LEV) = BGLAY + TFACTOT1*DELBGDN
               TFACTOT2 = ODTOT/(2.5+ODTOT)
               BBDTOT2(LEV) = BGLAY + TFACTOT2*DELBGDN
               TFACTOT3 = ODTOT/(1.666667+ODTOT)
               BBDTOT3(LEV) = BGLAY + TFACTOT3*DELBGDN
               BGLEV = PLFRAC * PLANKLEV(LEV,IBAND)
               DELBGUP = BGLEV - BGLAY
               BBUTOT1 = BGLAY + TFACTOT1*DELBGUP
               BBUTOT2 = BGLAY + TFACTOT2*DELBGUP
               BBUTOT3 = BGLAY + TFACTOT3*DELBGUP
               ATOT1(LEV) = 1.
               ATOT2(LEV) = 1.
               ATOT3(LEV) = 1.
               GASSRC1 = BGLAY + TFACGAS1*DELBGUP
               GASSRC2 = BGLAY + TFACGAS2*DELBGUP
               GASSRC3 = BGLAY + TFACGAS3*DELBGUP
            ELSEIF (ODCLD(LEV,IB) .GE. 10.) THEN
               ODTOT = ODEPTH + ODCLD(LEV,IB)
               TRANS1 = EXP(-ODEPTH)
               ATRANS1(LEV) = 1. - TRANS1
               ATRANS2(LEV) = 1. - TRANS1*TRANS1
               ATRANS3(LEV) = ATRANS1(LEV) + TRANS1 * ATRANS2(LEV)
               TFACGAS1 = ODEPTH/(5.+ODEPTH)
               BBDGAS1(LEV) = BGLAY + TFACGAS1*DELBGDN
               TFACGAS2 = ODEPTH/(2.5+ODEPTH)
               BBDGAS2(LEV) = BGLAY + TFACGAS2*DELBGDN
               TFACGAS3 = ODEPTH/(1.666667+ODEPTH)
               BBDGAS3(LEV) = BGLAY + TFACGAS3*DELBGDN
               TFACTOT1 = ODTOT/(5.+ODTOT)
               BBDTOT1(LEV) = BGLAY + TFACTOT1*DELBGDN
               TFACTOT2 = ODTOT/(2.5+ODTOT)
               BBDTOT2(LEV) = BGLAY + TFACTOT2*DELBGDN
               TFACTOT3 = ODTOT/(1.666667+ODTOT)
               BBDTOT3(LEV) = BGLAY + TFACTOT3*DELBGDN
               BGLEV = PLFRAC * PLANKLEV(LEV,IBAND)
               DELBGUP = BGLEV - BGLAY
               BBUTOT1 = BGLAY + TFACTOT1*DELBGUP
               BBUTOT2 = BGLAY + TFACTOT2*DELBGUP
               BBUTOT3 = BGLAY + TFACTOT3*DELBGUP
               ATOT1(LEV) = 1.
               ATOT2(LEV) = 1.
               ATOT3(LEV) = 1.
               GASSRC1 = (BGLAY + TFACGAS1*DELBGUP) * ATRANS1(LEV)
               GASSRC2 = (BGLAY + TFACGAS2*DELBGUP) * ATRANS2(LEV)
               GASSRC3 = (BGLAY + TFACGAS3*DELBGUP) * ATRANS3(LEV)
            ELSE
               ODTOT = ODEPTH + ODCLD(LEV,IB)
               TRANS1 = EXP(-ODEPTH)
               ATRANS1(LEV) = 1. - TRANS1
               ATRANS2(LEV) = 1. - TRANS1*TRANS1
               ATRANS3(LEV) = ATRANS1(LEV) + TRANS1 * ATRANS2(LEV)
               TFACGAS1 = ODEPTH/(5.+ODEPTH)
               BBDGAS1(LEV) = BGLAY + TFACGAS1*DELBGDN
               TFACGAS2 = ODEPTH/(2.5+ODEPTH)
               BBDGAS2(LEV) = BGLAY + TFACGAS2*DELBGDN
               TFACGAS3 = ODEPTH/(1.666667+ODEPTH)
               BBDGAS3(LEV) = BGLAY + TFACGAS3*DELBGDN
               TFACTOT1 = ODTOT/(5.+ODTOT)
               BBDTOT1(LEV) = BGLAY + TFACTOT1*DELBGDN
               TFACTOT2 = ODTOT/(2.5+ODTOT)
               BBDTOT2(LEV) = BGLAY + TFACTOT2*DELBGDN
               TFACTOT3 = ODTOT/(1.666667+ODTOT)
               BBDTOT3(LEV) = BGLAY + TFACTOT3*DELBGDN
               BGLEV = PLFRAC * PLANKLEV(LEV,IBAND)
               DELBGUP = BGLEV - BGLAY
               BBUTOT1 = BGLAY + TFACTOT1*DELBGUP
               BBUTOT2 = BGLAY + TFACTOT2*DELBGUP
               BBUTOT3 = BGLAY + TFACTOT3*DELBGUP
               ATOT1(LEV) = ATRANS1(LEV) + ABSCLD1(LEV,IB) -
     &              ATRANS1(LEV) * ABSCLD1(LEV,IB)
               ATOT2(LEV) = ATRANS2(LEV) + ABSCLD2(LEV,IB) -
     &              ATRANS1(LEV) * ABSCLD1(LEV,IB)
               ATOT3(LEV) = ATRANS3(LEV) + ABSCLD3(LEV,IB) -
     &              ATRANS1(LEV) * ABSCLD1(LEV,IB)
               GASSRC1 = (BGLAY + TFACGAS1*DELBGUP) * ATRANS1(LEV)
               GASSRC2 = (BGLAY + TFACGAS2*DELBGUP) * ATRANS2(LEV)
               GASSRC3 = (BGLAY + TFACGAS3*DELBGUP) * ATRANS3(LEV)
            ENDIF
C           Upward radiative transfer occurs here.
            RADLU1 = RADLU1 - RADLU1 * (ATRANS1(LEV) +
     &           EFCLFRAC1(LEV,IB) * (1. - ATRANS1(LEV))) +
     &           GASSRC1 + CLDFRAC(LEV) * 
     &           (BBUTOT1 * ATOT1(LEV) - GASSRC1)
            RADLU2 = RADLU2 - RADLU2 * (ATRANS2(LEV) +
     &           EFCLFRAC2(LEV,IB) * (1. - ATRANS2(LEV))) +
     &           GASSRC2 + CLDFRAC(LEV) * 
     &           (BBUTOT2 * ATOT2(LEV) - GASSRC2)
            RADLU3 = RADLU3 - RADLU3 * (ATRANS3(LEV) +
     &           EFCLFRAC3(LEV,IB) * (1. - ATRANS3(LEV))) +
     &           GASSRC3 + CLDFRAC(LEV) * 
     &           (BBUTOT3 * ATOT3(LEV) - GASSRC3)
C           Keep running total (over all IG) for each level.
            URAD1(LEV) = URAD1(LEV) + RADLU1
            URAD2(LEV) = URAD2(LEV) + RADLU2
            URAD3(LEV) = URAD3(LEV) + RADLU3
 2500    CONTINUE

C ***    Downward radiative transfer.
         DO 2600 LEV = NLAYERS-1, 0, -1
            GASSRC1 = BBDGAS1(LEV+1) * ATRANS1(LEV+1)
            RADLD1 = RADLD1 - RADLD1 * (ATRANS1(LEV+1) +
     &           EFCLFRAC1(LEV+1,IB) * (1. - ATRANS1(LEV+1))) +
     &           GASSRC1 + CLDFRAC(LEV+1) * 
     &           (BBDTOT1(LEV+1) * ATOT1(LEV+1) - GASSRC1)
            GASSRC2 = BBDGAS2(LEV+1) * ATRANS2(LEV+1)
            RADLD2 = RADLD2 - RADLD2 * (ATRANS2(LEV+1) +
     &           EFCLFRAC2(LEV+1,IB) * (1. - ATRANS2(LEV+1))) +
     &           GASSRC2 + CLDFRAC(LEV+1) * 
     &           (BBDTOT2(LEV+1) * ATOT2(LEV+1) - GASSRC2)
            GASSRC3 = BBDGAS3(LEV+1) * ATRANS3(LEV+1)
            RADLD3 = RADLD3 - RADLD3 * (ATRANS3(LEV+1) +
     &           EFCLFRAC3(LEV+1,IB) * (1. - ATRANS3(LEV+1))) +
     &           GASSRC3 + CLDFRAC(LEV+1) * 
     &           (BBDTOT3(LEV+1) * ATOT3(LEV+1) - GASSRC3)
            DRAD1(LEV) = DRAD1(LEV) + RADLD1
            DRAD2(LEV) = DRAD2(LEV) + RADLD2
            DRAD3(LEV) = DRAD3(LEV) + RADLD3
 2600    CONTINUE

         IG = IG + 1
         IF (IG .LE. NG(IBAND)) GO TO 1000
         
C ***    Process longwave output from band.
C ***    Calculate upward, downward, and net flux.
         DO 5000 LEV = NLAYERS, 0, -1
            UFLUX(LEV) = URAD1(LEV)*WTNUM(1) +
     &           URAD2(LEV)*WTNUM(2) + URAD3(LEV)*WTNUM(3)
            DFLUX(LEV) = DRAD1(LEV)*WTNUM(1) +
     &           DRAD2(LEV)*WTNUM(2) + DRAD3(LEV)*WTNUM(3)
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

      TOTUCLFL = TOTUCLFL * FLUXFAC
      TOTDCLFL = TOTDCLFL * FLUXFAC

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

 9000 CONTINUE

      RETURN
      END   