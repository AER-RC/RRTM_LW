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

      IMPLICIT DOUBLE PRECISION (V)                                     

      COMMON /CONSTANTS/ PI,FLUXFAC,HEATFAC
      COMMON /FEATURES/  NG(NBANDS),NSPA(NBANDS),NSPB(NBANDS)
      COMMON /BANDS/     WAVENUM1(NBANDS),WAVENUM2(NBANDS),
     &                   DELWAVE(NBANDS)
      COMMON /CONTROL/   NUMANGS, IOUT, ISTART, IEND
      COMMON /PROFILE/   NLAYERS,PAVEL(MXLAY),TAVEL(MXLAY),
     &                   PZ(0:MXLAY),TZ(0:MXLAY)
      COMMON /SURFACE/   TBOUND,IREFLECT,SEMISS(NBANDS)
      COMMON /CLOUDDAT/  NCBANDS,CLDFRAC(MXLAY),TAUCLOUD(MXLAY,NBANDS)
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
                                       
      DIMENSION BBUGAS1(MXLAY),BBUGAS2(MXLAY),BBUGAS3(MXLAY)
      DIMENSION BBUTOT1(MXLAY),BBUTOT2(MXLAY),BBUTOT3(MXLAY)
      DIMENSION ATRANS1(MXLAY),ATRANS2(MXLAY),ATRANS3(MXLAY)
      DIMENSION UFLUX(0:MXLAY),DFLUX(0:MXLAY)
      DIMENSION DRAD1(0:MXLAY),URAD1(0:MXLAY)
      DIMENSION DRAD2(0:MXLAY),URAD2(0:MXLAY)
      DIMENSION DRAD3(0:MXLAY),URAD3(0:MXLAY)
      DIMENSION WTNUM(3), ODCLD(MXLAY,NBANDS)
      DIMENSION ATOT1(MXLAY),ATOT2(MXLAY),ATOT3(MXLAY)
      DIMENSION ABSCLD1(MXLAY,NBANDS),ABSCLD2(MXLAY,NBANDS)
      DIMENSION EFCLFRAC1(MXLAY,NBANDS),ABSCLD3(MXLAY,NBANDS)
      DIMENSION EFCLFRAC2(MXLAY,NBANDS),EFCLFRAC3(MXLAY,NBANDS)
      DIMENSION IPAT(16,0:2)

C     These arrays indicate the spectral 'region' (used in the 
C     calculation of ice cloud optical depths) corresponding
C     to each spectral band.  See cldprop.f for more details.
      DATA IPAT /1,1,1,1,1,1,1,1,1, 1, 1, 1, 1, 1, 1, 1,
     &           1,2,3,3,3,4,4,4,5, 5, 5, 5, 5, 5, 5, 5,
     &           1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16/

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
            IB = IPAT(IBAND,0)
         ELSEIF (NCBANDS .EQ.  5) THEN
            IB = IPAT(IBAND,1)
         ELSEIF (NCBANDS .EQ. 16) THEN
            IB = IPAT(IBAND,2)
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
         RADLD1 = 0.
         RADLD2 = 0.
         RADLD3 = 0.
C ***    Downward radiative transfer loop.  Due to the simple form taken
C        by certain equations when the optical depth is small, this 
C        condition is tested for.  In all cases, the labels 1, 2, and
C        3 refer to the three angles.
         BGLEV = FRACS(NLAYERS,IG) * PLANKLEV(NLAYERS,IBAND)
         DO 2500 LEV = NLAYERS, 1, -1
            BLAY = PLANKLAY(LEV,IBAND)
            PLFRAC = FRACS(LEV,IG)
            BGLAY = PLFRAC * BLAY
            ODEPTH = SECANG * TAUG(LEV,IG) 
            DELBGUP = BGLEV - BGLAY

C           Here are some variable definitions:
C             ODTOT      optical depth of gas and cloud
C             ATRANS     absorptivity for gas only
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
               BBUGAS1(LEV) = BGLAY
               BBUGAS2(LEV) = BGLAY
               BBUGAS3(LEV) = BGLAY
               BBUTOT1(LEV) = BGLAY
               BBUTOT2(LEV) = BGLAY
               BBUTOT3(LEV) = BGLAY
               BGLEV = PLFRAC * PLANKLEV(LEV-1,IBAND)
               BBDTOT1 = BGLAY
               BBDTOT2 = BGLAY
               BBDTOT3 = BGLAY
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
               BBUGAS1(LEV) = BGLAY
               BBUGAS2(LEV) = BGLAY
               BBUGAS3(LEV) = BGLAY
               TFACTOT1 = ODTOT/(5.+ODTOT)
               BBUTOT1(LEV) = BGLAY + TFACTOT1*DELBGUP
               TFACTOT2 = ODTOT/(2.5+ODTOT)
               BBUTOT2(LEV) = BGLAY + TFACTOT2*DELBGUP
               TFACTOT3 = ODTOT/(1.666667+ODTOT)
               BBUTOT3(LEV) = BGLAY + TFACTOT3*DELBGUP
               BGLEV = PLFRAC * PLANKLEV(LEV-1,IBAND)
               DELBGDN = BGLEV - BGLAY
               BBDTOT1 = BGLAY + TFACTOT1*DELBGDN
               BBDTOT2 = BGLAY + TFACTOT2*DELBGDN
               BBDTOT3 = BGLAY + TFACTOT3*DELBGDN
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
               BBUGAS1(LEV) = BGLAY + TFACGAS1*DELBGUP
               TFACGAS2 = ODEPTH/(2.5+ODEPTH)
               BBUGAS2(LEV) = BGLAY + TFACGAS2*DELBGUP
               TFACGAS3 = ODEPTH/(1.666667+ODEPTH)
               BBUGAS3(LEV) = BGLAY + TFACGAS3*DELBGUP
               TFACTOT1 = ODTOT/(5.+ODTOT)
               BBUTOT1(LEV) = BGLAY + TFACTOT1*DELBGUP
               TFACTOT2 = ODTOT/(2.5+ODTOT)
               BBUTOT2(LEV) = BGLAY + TFACTOT2*DELBGUP
               TFACTOT3 = ODTOT/(1.666667+ODTOT)
               BBUTOT3(LEV) = BGLAY + TFACTOT3*DELBGUP
               BGLEV = PLFRAC * PLANKLEV(LEV-1,IBAND)
               DELBGDN = BGLEV - BGLAY
               BBDTOT1 = BGLAY + TFACTOT1*DELBGDN
               BBDTOT2 = BGLAY + TFACTOT2*DELBGDN
               BBDTOT3 = BGLAY + TFACTOT3*DELBGDN
               ATOT1(LEV) = 1.
               ATOT2(LEV) = 1.
               ATOT3(LEV) = 1.
               GASSRC1 = BGLAY + TFACGAS1*DELBGDN
               GASSRC2 = BGLAY + TFACGAS2*DELBGDN
               GASSRC3 = BGLAY + TFACGAS3*DELBGDN
            ELSEIF (ODCLD(LEV,IB) .GE. 10.) THEN
               ODTOT = ODEPTH + ODCLD(LEV,IB)
               TRANS1 = EXP(-ODEPTH)
               ATRANS1(LEV) = 1. - TRANS1
               ATRANS2(LEV) = 1. - TRANS1*TRANS1
               ATRANS3(LEV) = ATRANS1(LEV) + TRANS1 * ATRANS2(LEV)
               TFACGAS1 = ODEPTH/(5.+ODEPTH)
               BBUGAS1(LEV) = BGLAY + TFACGAS1*DELBGUP
               TFACGAS2 = ODEPTH/(2.5+ODEPTH)
               BBUGAS2(LEV) = BGLAY + TFACGAS2*DELBGUP
               TFACGAS3 = ODEPTH/(1.666667+ODEPTH)
               BBUGAS3(LEV) = BGLAY + TFACGAS3*DELBGUP
               TFACTOT1 = ODTOT/(5.+ODTOT)
               BBUTOT1(LEV) = BGLAY + TFACTOT1*DELBGUP
               TFACTOT2 = ODTOT/(2.5+ODTOT)
               BBUTOT2(LEV) = BGLAY + TFACTOT2*DELBGUP
               TFACTOT3 = ODTOT/(1.666667+ODTOT)
               BBUTOT3(LEV) = BGLAY + TFACTOT3*DELBGUP
               BGLEV = PLFRAC * PLANKLEV(LEV-1,IBAND)
               DELBGDN = BGLEV - BGLAY
               BBDTOT1 = BGLAY + TFACTOT1*DELBGDN
               BBDTOT2 = BGLAY + TFACTOT2*DELBGDN
               BBDTOT3 = BGLAY + TFACTOT3*DELBGDN
               ATOT1(LEV) = 1.
               ATOT2(LEV) = 1.
               ATOT3(LEV) = 1.
               GASSRC1 = (BGLAY + TFACGAS1*DELBGDN) * ATRANS1(LEV)
               GASSRC2 = (BGLAY + TFACGAS2*DELBGDN) * ATRANS2(LEV)
               GASSRC3 = (BGLAY + TFACGAS3*DELBGDN) * ATRANS3(LEV)
            ELSE
               ODTOT = ODEPTH + ODCLD(LEV,IB)
               TRANS1 = EXP(-ODEPTH)
               ATRANS1(LEV) = 1. - TRANS1
               ATRANS2(LEV) = 1. - TRANS1*TRANS1
               ATRANS3(LEV) = ATRANS1(LEV) + TRANS1 * ATRANS2(LEV)
               TFACGAS1 = ODEPTH/(5.+ODEPTH)
               BBUGAS1(LEV) = BGLAY + TFACGAS1*DELBGUP
               TFACGAS2 = ODEPTH/(2.5+ODEPTH)
               BBUGAS2(LEV) = BGLAY + TFACGAS2*DELBGUP
               TFACGAS3 = ODEPTH/(1.666667+ODEPTH)
               BBUGAS3(LEV) = BGLAY + TFACGAS3*DELBGUP
               TFACTOT1 = ODTOT/(5.+ODTOT)
               BBUTOT1(LEV) = BGLAY + TFACTOT1*DELBGUP
               TFACTOT2 = ODTOT/(2.5+ODTOT)
               BBUTOT2(LEV) = BGLAY + TFACTOT2*DELBGUP
               TFACTOT3 = ODTOT/(1.666667+ODTOT)
               BBUTOT3(LEV) = BGLAY + TFACTOT3*DELBGUP
               BGLEV = PLFRAC * PLANKLEV(LEV-1,IBAND)
               DELBGDN = BGLEV - BGLAY
               BBDTOT1 = BGLAY + TFACTOT1*DELBGDN
               BBDTOT2 = BGLAY + TFACTOT2*DELBGDN
               BBDTOT3 = BGLAY + TFACTOT3*DELBGDN
               ATOT1(LEV) = ATRANS1(LEV) + ABSCLD1(LEV,IB) -
     &              ATRANS1(LEV) * ABSCLD1(LEV,IB)
               ATOT2(LEV) = ATRANS2(LEV) + ABSCLD2(LEV,IB) -
     &              ATRANS2(LEV) * ABSCLD2(LEV,IB)
               ATOT3(LEV) = ATRANS3(LEV) + ABSCLD3(LEV,IB) -
     &              ATRANS3(LEV) * ABSCLD3(LEV,IB)
               GASSRC1 = (BGLAY + TFACGAS1*DELBGDN) * ATRANS1(LEV)
               GASSRC2 = (BGLAY + TFACGAS2*DELBGDN) * ATRANS2(LEV)
               GASSRC3 = (BGLAY + TFACGAS3*DELBGDN) * ATRANS3(LEV)
            ENDIF
C           Downward radiative transfer occurs here.
            RADLD1 = RADLD1 - RADLD1 * (ATRANS1(LEV) +
     &           EFCLFRAC1(LEV,IB) * (1. - ATRANS1(LEV))) +
     &           GASSRC1 + CLDFRAC(LEV) * 
     &           (BBDTOT1 * ATOT1(LEV) - GASSRC1)
            RADLD2 = RADLD2 - RADLD2 * (ATRANS2(LEV) +
     &           EFCLFRAC2(LEV,IB) * (1. - ATRANS2(LEV))) +
     &           GASSRC2 + CLDFRAC(LEV) * 
     &           (BBDTOT2 * ATOT2(LEV) - GASSRC2)
            RADLD3 = RADLD3 - RADLD3 * (ATRANS3(LEV) +
     &           EFCLFRAC3(LEV,IB) * (1. - ATRANS3(LEV))) +
     &           GASSRC3 + CLDFRAC(LEV) * 
     &           (BBDTOT3 * ATOT3(LEV) - GASSRC3)
C           Keep running total (over all IG) for each level.
            DRAD1(LEV-1) = DRAD1(LEV-1) + RADLD1
            DRAD2(LEV-1) = DRAD2(LEV-1) + RADLD2
            DRAD3(LEV-1) = DRAD3(LEV-1) + RADLD3
 2500    CONTINUE

C ***    Upward radiative transfer.
         RAD0 = FRACS(1,IG) * PLANKBND(IBAND)
C        Add in reflection of surface downward radiance.
         REFLECT = 1. - SEMISS(IBAND)
         IF (IREFLECT .EQ. 1) THEN
C           Specular reflection.
            RADLU1 = RAD0 + REFLECT * RADLD1
            RADLU2 = RAD0 + REFLECT * RADLD2
            RADLU3 = RAD0 + REFLECT * RADLD3
         ELSE
C           Lambertian reflection.
            RAD = 2. * (RADLD1*WTNUM(1) + RADLD2*WTNUM(2) + 
     &           RADLD3*WTNUM(3))
            RADLU1 = RAD0 + REFLECT * RAD
            RADLU2 = RADLU1
            RADLU3 = RADLU1
         ENDIF
         URAD1(0) = URAD1(0) + RADLU1
         URAD2(0) = URAD2(0) + RADLU2
         URAD3(0) = URAD3(0) + RADLU3
         DO 2600 LEV = 1, NLAYERS
            GASSRC1 = BBUGAS1(LEV) * ATRANS1(LEV)
            RADLU1 = RADLU1 - RADLU1 * (ATRANS1(LEV) +
     &           EFCLFRAC1(LEV,IB) * (1. - ATRANS1(LEV))) +
     &           GASSRC1 + CLDFRAC(LEV) * 
     &           (BBUTOT1(LEV) * ATOT1(LEV) - GASSRC1)
            GASSRC2 = BBUGAS2(LEV) * ATRANS2(LEV)
            RADLU2 = RADLU2 - RADLU2 * (ATRANS2(LEV) +
     &           EFCLFRAC2(LEV,IB) * (1. - ATRANS2(LEV))) +
     &           GASSRC2 + CLDFRAC(LEV) * 
     &           (BBUTOT2(LEV) * ATOT2(LEV) - GASSRC2)
            GASSRC3 = BBUGAS3(LEV) * ATRANS3(LEV)
            RADLU3 = RADLU3 - RADLU3 * (ATRANS3(LEV) +
     &           EFCLFRAC3(LEV,IB) * (1. - ATRANS3(LEV))) +
     &           GASSRC3 + CLDFRAC(LEV) * 
     &           (BBUTOT3(LEV) * ATOT3(LEV) - GASSRC3)
            URAD1(LEV) = URAD1(LEV) + RADLU1
            URAD2(LEV) = URAD2(LEV) + RADLU2
            URAD3(LEV) = URAD3(LEV) + RADLU3
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
