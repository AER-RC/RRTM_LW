C     path:      %P%
C     revision:  $Revision$
C     created:   $Date$  
C     presently: %H%  %T%
      SUBROUTINE CLDPROP(ICLDATM)

C     Purpose:  Compute the cloud optical depth(s) for each cloudy
C               layer.

      PARAMETER (MXLAY=203)
      PARAMETER (MXCBANDS = 5)

      COMMON /CONTROL/   NUMANGS, IOUT, ISTART, IEND, ICLD
      COMMON /PROFILE/   NLAYERS,PAVEL(MXLAY),TAVEL(MXLAY),
     &                   PZ(0:MXLAY),TZ(0:MXLAY)
      COMMON /CLOUDIN/   INFLAG,CLDDAT1(MXLAY),CLDDAT2(MXLAY),
     &                   CLDDAT3(MXLAY)
      COMMON /CLOUDDAT/  NCBANDS,CLDFRAC(MXLAY),TAUCLOUD(MXLAY,MXCBANDS)
      COMMON /HVERSN/    HVRRTM,HVRREG,HVRRTR,HVRATM,HVRSET,HVRTAU,
     *                   HVRRGC,HVRRTC,HVRCLD,HVRDUM,HVRUTL,HVREXT

      CHARACTER*8 HVRRTM,HVRREG,HVRRTR,HVRATM,HVRSET,HVRTAU,
     *            HVRRGC,HVRRTC,HVRCLD,HVRDUM,HVRUTL,HVREXT
      DIMENSION ABSICE3(2), ABSICE4(2,MXCBANDS)

      DATA EPS /1.E-6/

C     Explanation of the method for each value of INFLAG.
C     INFLAG = 1:  For each cloudy layer, the cloud fraction and (gray)
C                  optical depth are input.  
C     INFLAG = 2:  For each cloudy layer, the cloud fraction and liquid
C                  water path (g/m2) are input.  The (gray) cloud optical 
C                  depth is computed as in CCM2.
C     INFLAG = 3:  For each cloudy layer, the cloud fraction, cloud 
C                  water path (g/m2), cloud ice fraction, and ice effective
C                  radius (microns) are input.  The (gray) cloud optical 
C                  depth is computed as in CCM3.
C     INFLAG = 4:  For each cloudy layer, the cloud fraction, cloud 
C                  water path (g/m2), cloud ice fraction, and ice effective
C                  radius (microns) are input.  The cloud optical depth 
C                  is computed as in CCM3, except that the parameters in the
C                  formula for the optical depth due to clouds with ice 
C                  depend on the spectral region, as in Ebert and Curry, JGR
C                  97, 3831-3836 (1992).  The spectral regions in this work
C                  have been matched with the spectral bands in RRTM to as
C                  great an extent as possible:  
C                  E&C 1      IB = 5      RRTM bands 9-16
C                  E&C 2      IB = 4      RRTM bands 6-8
C                  E&C 3      IB = 3      RRTM bands 3-5
C                  E&C 4      IB = 2      RRTM band 2
C                  E&C 5      IB = 1      RRTM band 1
C  
C     ABSLIQn is the liquid water absorption coefficient (m2/g) 
C     for INFLAG=n.
      DATA ABSLIQ2 /0.0602410/
      DATA ABSLIQ3 /0.0903614/
      DATA ABSLIQ4 /0.0903614/

C     ABSICEn(J,IB) are the parameters needed to compute the liquid water 
C     absorption coefficient in spectral region IB for INFLAG=n.  The 
C     units of ABSICEn(1,IB) are m2/g and ABSICEn(2,IB) has units 
C     (microns (m2/g)).
      DATA ABSICE3 /0.005,  1.0/
      DATA ABSICE4 /0.0036, 1.136,
     &              0.0068, 0.600,
     &              0.0003, 1.338,
     &              0.0016, 1.166,
     &              0.0020, 1.118/

      HVRCLD = '$Revision$'

      ICLDATM = 0
      NCBANDS = 1
      DO 3000 LAY = 1, NLAYERS
         IF (CLDFRAC(LAY) .GE. EPS) THEN
            ICLDATM = 1
            IF (INFLAG .EQ. 1) THEN
               TAUCLOUD(LAY,1) = CLDDAT1(LAY)
            ELSEIF(INFLAG .EQ. 2) THEN
               CWP = CLDDAT1(LAY)
               TAUCLOUD(LAY,1) = ABSLIQ2 * CWP
            ELSEIF(INFLAG .EQ. 3) THEN
               CWP = CLDDAT1(LAY)
               FICE = CLDDAT2(LAY)
               RADICE = CLDDAT3(LAY)
               IF (RADICE .EQ. 0.0) THEN
                  IF (FICE .NE. 0.0) STOP 'ICE RADIUS MUST NOT BE 0.0'
                  RADICE = 1.
               ENDIF
               TAUCLOUD(LAY,1) = CWP * ((1.-FICE) * ABSLIQ3 +
     &              FICE * (ABSICE3(1) + ABSICE3(2)/RADICE))
            ELSE
               NCBANDS = 5
               CWP = CLDDAT1(LAY)
               FICE = CLDDAT2(LAY)
               RADICE = CLDDAT3(LAY)
               IF (RADICE .EQ. 0.0) THEN
                  IF (FICE .NE. 0.0) STOP 'ICE RADIUS MUST NOT BE 0.0'
                  RADICE = 1.
               ENDIF
C              Loop over spectral regions.  See note above to see which
C              spectral bands in RRTM correspond to each spectral region.
               DO 2000 IB = 1, NCBANDS
                  TAUCLOUD(LAY,IB) = CWP * ((1.-FICE) * ABSLIQ4 +
     &                 FICE * (ABSICE4(1,IB) + ABSICE4(2,IB)/RADICE))
 2000          CONTINUE
            ENDIF
         ENDIF
 3000 CONTINUE

      RETURN
      END
