      SUBROUTINE RTREG

C *** This program calculates the upward fluxes, downward fluxes,
C     and heating rates for an arbitrary atmosphere.  The input to
C     this program is the atmospheric profile and all needed Planck
C     function information.  First-order standard Gaussian quadrature 
C     is used for the angle integration.

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
     &                   FNET(0:MXLAY), HTR(MXLAY)

      DIMENSION ATRANS(MXLAY),BBD(MXLAY)
      DIMENSION UFLUX(0:MXLAY),DFLUX(0:MXLAY)
      DIMENSION DRAD(0:MXLAY-1,MXANG),URAD(0:MXLAY,MXANG)
      DIMENSION SECANG(MXANG),ANGWEIGH(MXANG)
      DIMENSION SECREG(MXANG,MXANG),WTREG(MXANG,MXANG)

C *** When standard first-order Gaussian quadrature is chosen as
C     the method to approximate the integral over angles that yields
C     flux from radiances, then SECREG(I,J) is the secant of the Ith  
C     (out of a total of J angles) and WTREG(I,J) is the corresponding
C     weight.
      DATA SECREG(1,1) / 1.5/
      DATA SECREG(2,2) / 2.81649655/, SECREG(1,2) / 1.18350343/
      DATA SECREG(3,3) / 4.70941630/, SECREG(2,3) / 1.69338507/
      DATA SECREG(1,3) / 1.09719858/
      DATA SECREG(4,4) / 7.15513024/, SECREG(3,4) / 2.40148179/
      DATA SECREG(2,4) / 1.38282560/, SECREG(1,4) / 1.06056257/
      DATA WTREG(1,1) / 0.50000000/
      DATA WTREG(2,2) /0.1819586183/, WTREG(1,2) /0.3180413817/
      DATA WTREG(3,3) /0.0698269799/, WTREG(2,3) /0.2292411064/
      DATA WTREG(1,3) /0.2009319137/
      DATA WTREG(4,4) /0.0311809710/, WTREG(3,4) /0.1298475476/
      DATA WTREG(2,4) /0.2034645680/, WTREG(1,4) /0.1355069134/

C *** Load angle data in arrays depending on angular quadrature scheme.
      DO 100 IANG = 1, NUMANGS
         SECANG(IANG) = SECREG(IANG,NUMANGS)
         ANGWEIGH(IANG) = WTREG(IANG,NUMANGS)
 100  CONTINUE

      TOTUFLUX(0) = 0.0
      TOTDFLUX(0) = 0.0
      DO 200 LAY = 1, NLAYERS
         TOTUFLUX(LAY) = 0.0
         TOTDFLUX(LAY) = 0.0
 200  CONTINUE

C *** Loop over frequency bands.
      DO 6000 IBAND = ISTART, IEND
         DO 1100 IANG = 1, NUMANGS
            DRAD(NLAYERS,IANG) = 0.
 1100    CONTINUE

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
         
C ***    Loop over each angle for which the radiance is to be computed.
         DO 3000 IANG = 1, NUMANGS
C ***       Radiative transfer starts here.
            RADLU = FRACS(1,IG) * PLANKBND(IBAND)
            URAD(0,IANG) = URAD(0,IANG) + RADLU
            RADLD = 0.

C ***       Upward radiative transfer.  
            DO 2500 LEV = 1, NLAYERS
               BLAY = PLANKLAY(LEV,IBAND)
               PLFRAC = FRACS(LEV,IG)
               ODEPTH = SECANG(IANG) * TAUG(LEV,IG)
               ATRANS(LEV) = 1. - EXP(-ODEPTH)
               DPLANKUP = PLANKLEV(LEV,IBAND) - BLAY
               DPLANKDN = PLANKLEV(LEV-1,IBAND) - BLAY
               TAUSFAC = ODEPTH/(5.+ODEPTH)
               BBU = PLFRAC * (BLAY + TAUSFAC * DPLANKUP)
               RADLU = RADLU + (BBU-RADLU)*ATRANS(LEV)
               URAD(LEV,IANG) = URAD(LEV,IANG) + RADLU
               BBD(LEV) = PLFRAC * (BLAY + TAUSFAC * DPLANKDN)
 2500       CONTINUE

C ***       Downward radiative transfer.
            DO 2600 LEV = NLAYERS-1, 0, -1
               RADLD = RADLD + (BBD(LEV+1)-RADLD)*ATRANS(LEV+1)
               DRAD(LEV,IANG) = DRAD(LEV,IANG) + RADLD
 2600       CONTINUE
 3000    CONTINUE 

         IG = IG + 1
         IF (IG .LE. 16) GO TO 1000

C ***    Calculate upward, downward, and net flux.
         DO 5000 LEV = NLAYERS, 0, -1
            UFLUX(LEV) = 0.0
            DFLUX(LEV) = 0.0
            DO 4500 IANG = 1, NUMANGS
               UFLUX(LEV) = UFLUX(LEV) + URAD(LEV,IANG)*ANGWEIGH(IANG)
               DFLUX(LEV) = DFLUX(LEV) + DRAD(LEV,IANG)*ANGWEIGH(IANG)
               URAD(LEV,IANG) = 0.0
               DRAD(LEV,IANG) = 0.0
 4500       CONTINUE
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


