*******************************************************************************
*                                                                             *
*                  Optical depths developed for the                           *
*                                                                             *
*                RAPID RADIATIVE TRANSFER MODEL (RRTM)                        *
*                                                                             *
*                                                                             *
*            ATMOSPHERIC AND ENVIRONMENTAL RESEARCH, INC.                     *
*                        840 MEMORIAL DRIVE                                   *
*                        CAMBRIDGE, MA 02139                                  *
*                                                                             *
*                                                                             *
*                           ELI J. MLAWER                                     *
*                         STEVEN J. TAUBMAN~                                  *
*                         SHEPARD A. CLOUGH                                   *
*                                                                             *
*                                                                             *
*                         ~currently at GFDL                                  *
*                                                                             *
*                                                                             *
*                                                                             *
*                       email:  mlawer@aer.com                                *
*                                                                             *
*        The authors wish to acknowledge the contributions of the             *
*        following people:  Patrick D. Brown, Michael J. Iacono,              *
*        Ronald E. Farren, Luke Chen, Robert Bergstrom.                       *
*                                                                             *
*******************************************************************************
*     TAUMOL                                                                  *
*                                                                             *
*     This file contains the subroutines TAUGBn (where n goes from            *
*     1 to 16).  TAUGBn calculates the optical depths per g-value             *
*     and layer for band n.                                                   *
*                                                                             *
*  Output:  optical depths (unitless)                                         *
*                                                                             *
*     COMMON /TAUGCOM/  TAUG(MXLAY,MG)                                        *
*                                                                             *
*  Input                                                                      *
*                                                                             *
*     COMMON /FEATURES/ NG(NBANDS),NSPA(NBANDS),NSPB(NBANDS)                  *
*     COMMON /PROFILE/  NLAYERS,PAVEL(MXLAY),TAVEL(MXLAY),                    *
*    &                  PZ(0:MXLAY),TZ(0:MXLAY)                               *
*     COMMON /PROFDATA/ LAYTROP,COLH2O(MXLAY),COLCO2(MXLAY),                  *
*                       COLO3(MXLAY)                                          *
*     COMMON /INTFAC/   FAC000(MXLAY,NBANDS),FAC001(MXLAY,NBANDS),            *
*    &                  FAC010(MXLAY,NBANDS),FAC011(MXLAY,NBANDS),            *
*    &                  FAC100(MXLAY,NBANDS),FAC101(MXLAY,NBANDS),            *
*    &                  FAC110(MXLAY,NBANDS),FAC111(MXLAY,NBANDS)             *
*     COMMON /INTIND/   JP(MXLAY),JT(MXLAY),JT1(MXLAY),                       *
*    &                  JS(MXLAY,NBANDS),JS1(MXLAY,NBANDS)                    *
*     COMMON /SELF/     SELFFAC(MXLAY), SELFFRAC(MXLAY), INDSELF(MXLAY)       *
*                                                                             *
*     Description:                                                            *
*     NG(IBAND) - number of g-values in band IBAND                            *
*     NSPA(IBAND) - for the lower atmosphere, the number of reference         *
*                   atmospheres that are stored for band IBAND per            *
*                   pressure level and temperature.  Each of these            *
*                   atmospheres has different relative amounts of the         *
*                   key species for the band (i.e. different binary           *
*                   species parameters).                                      *
*     NSPB(IBAND) - same for upper atmosphere                                 *
*     PAVEL - layer pressures (mb)                                            *
*     TAVEL - layer temperatures (degrees K)                                  *
*     PZ - level pressures (mb)                                               *
*     TZ - level tmeperatures (degrees K)                                     *
*     LAYTROP - layer at which switch is made from one combination of         *
*               key species to another                                        *
*     COLH2O, COLCO2,COLO3 - column amounts of water vapor,carbon dioxide,    *
*                            and ozone, respectively (molecules/cm**2)        *
*     FACijk(LAY,IBAND) - for band IBAND and layer LAY, these are the         *
*                         interpolation factors that multiply the             *
*                         appropriate reference k-values.  These factors      *
*                         sum to 1.  A value of 0 (1) for i,j,k indicates     *
*                         that the corresponding factor multiplies the        *
*                         reference k-value for the lower (higher) of the     *
*                         two appropriate binary species parameters,          *
*                         temperatures, and altitudes, respectively.          *
*     JP - the index of the lower (in altitude) of the two appropriate        *
*          reference pressure levels needed for interpolation                 *
*     JT, JT1 - the indices of the lower of the two appropriate reference     *
*               temperatures needed for interpolation (for pressure           *
*               levels JP and JP+1, respectively)                             *
*     JS, JS1 - the indices of the lower of the two appropriate reference     *
*               binary species parameters needed for interpolation for each   *
*               band (for pressure levels JP and JP+1, respectively)          *
*     SELFFAC - scale factor needed to water vapor self-continuum, equals     *
*               (water vapor density)/(atmospheric density at 296K and        *
*               1013 mb)                                                      *
*     SELFFRAC - factor needed for temperature interpolation of reference     *
*                water vapor self-continuum data                              *
*     INDSELF - index of the lower of the two appropriate reference           *
*               temperatures needed for the self-continuum interpolation      *
*                                                                             *
*  Data input                                                                 *
*     COMMON /Kn/ KA(NSPA(n),5,13,MG), KB(NSPB(n),5,13:59,MG), SELFREF(10,MG) *
*        (note:  n is the band number)                                        *
*                                                                             *
*     Description:                                                            *
*     KA - k-values for low reference atmospheres (no water vapor             *
*          self-continuum) (units: cm**2/molecule)                            *
*     KB - k-values for high reference atmospheres (all sources)              *
*          (units: cm**2/molecule)                                            *
*     SELFREF - k-values for water vapor self-continuum for reference         *
*               atmospheres (used below LAYTROP)                              *
*               (units: cm**2/molecule)                                       *
*                                                                             *
*     DIMENSION ABSA(65*NSPA(n),MG), ABSB(235*NSPB(n),MG)                     *
*     EQUIVALENCE (KA,ABSA),(KB,ABSB)                                         *
*                                                                             *
*******************************************************************************


      SUBROUTINE TAUGB1

C     Written by Eli J. Mlawer, Atmospheric & Environmental Research.

C     BAND 1:  10-250 cm-1 (low - H2O; high - H2O)

      PARAMETER (MG=16, MXLAY=200, NBANDS=16)

C  Output

      COMMON /TAUGCOM/  TAUG(MXLAY,MG)

C  Input

      COMMON /FEATURES/ NG(NBANDS),NSPA(NBANDS),NSPB(NBANDS)
      COMMON /PROFILE/  NLAYERS,PAVEL(MXLAY),TAVEL(MXLAY),
     &                  PZ(0:MXLAY),TZ(0:MXLAY)
      COMMON /PROFDATA/ LAYTROP,COLH2O(MXLAY),COLCO2(MXLAY),
     &                  COLO3(MXLAY)
      COMMON /INTFAC/   FAC000(MXLAY,NBANDS),FAC001(MXLAY,NBANDS),
     &                  FAC010(MXLAY,NBANDS),FAC011(MXLAY,NBANDS),
     &                  FAC100(MXLAY,NBANDS),FAC101(MXLAY,NBANDS),
     &                  FAC110(MXLAY,NBANDS),FAC111(MXLAY,NBANDS)
      COMMON /INTIND/   JP(MXLAY),JT(MXLAY),JT1(MXLAY),
     &                  JS(MXLAY,NBANDS),JS1(MXLAY,NBANDS)
      COMMON /SELF/     SELFFAC(MXLAY), SELFFRAC(MXLAY), INDSELF(MXLAY)
      COMMON /K1/   KA(5,13,MG), KB(5,13:59,MG) , SELFREF(10,MG)

      EQUIVALENCE (KA,ABSA),(KB,ABSB)

      DIMENSION ABSA(65,MG),ABSB(235,MG)
      REAL KA,KB

C     Compute the optical depth by interpolating in ln(pressure) and 
C     temperature.  Below LAYTROP, the water vapor self-continuum 
C     is interpolated (in temperature) separately.  
      DO 2500 LAY = 1, LAYTROP
         IND0 = ((JP(LAY)-1)*5+(JT(LAY)-1))*NSPA(1) + JS(LAY,1)
         IND1 = (JP(LAY)*5+(JT1(LAY)-1))*NSPA(1) + JS1(LAY,1)
         INDS = INDSELF(LAY)
         DO 2000 IG = 1, NG(1)
            TAUG(LAY,IG) = COLH2O(LAY) *
     &          (FAC000(LAY,1) * ABSA(IND0,IG) +
     &           FAC010(LAY,1) * ABSA(IND0+1,IG) +
     &           FAC001(LAY,1) * ABSA(IND1,IG) + 
     &           FAC011(LAY,1) * ABSA(IND1+1,IG) +
     &           SELFFAC(LAY) * (SELFREF(INDS,IG) + 
     &           SELFFRAC(LAY) *
     &           (SELFREF(INDS+1,IG) - SELFREF(INDS,IG))))
 2000    CONTINUE
 2500 CONTINUE
      DO 3500 LAY = LAYTROP+1, NLAYERS
         IND0 = ((JP(LAY)-13)*5+(JT(LAY)-1))*NSPB(1) + JS(LAY,1)
         IND1 = ((JP(LAY)-12)*5+(JT1(LAY)-1))*NSPB(1) + JS1(LAY,1)
         DO 3000 IG = 1, NG(1)
            TAUG(LAY,IG) = COLH2O(LAY) * 
     &          (FAC000(LAY,1) * ABSB(IND0,IG) +
     &           FAC010(LAY,1) * ABSB(IND0+1,IG) +
     &           FAC001(LAY,1) * ABSB(IND1,IG) + 
     &           FAC011(LAY,1) * ABSB(IND1+1,IG)) 
 3000    CONTINUE
 3500 CONTINUE
      RETURN
      END

C----------------------------------------------------------------------------

      SUBROUTINE TAUGB2

C     BAND 2:  250-500 cm-1 (low - H2O; high - H2O)

      PARAMETER (MG=16, MXLAY=200, NBANDS=16)

C  Output

      COMMON /TAUGCOM/  TAUG(MXLAY,MG)

C  Input

      COMMON /FEATURES/ NG(NBANDS),NSPA(NBANDS),NSPB(NBANDS)
      COMMON /PROFILE/  NLAYERS,PAVEL(MXLAY),TAVEL(MXLAY),
     &                  PZ(0:MXLAY),TZ(0:MXLAY)
      COMMON /PROFDATA/ LAYTROP,COLH2O(MXLAY),COLCO2(MXLAY),
     &                  COLO3(MXLAY)
      COMMON /INTFAC/   FAC000(MXLAY,NBANDS),FAC001(MXLAY,NBANDS),
     &                  FAC010(MXLAY,NBANDS),FAC011(MXLAY,NBANDS),
     &                  FAC100(MXLAY,NBANDS),FAC101(MXLAY,NBANDS),
     &                  FAC110(MXLAY,NBANDS),FAC111(MXLAY,NBANDS)
      COMMON /INTIND/   JP(MXLAY),JT(MXLAY),JT1(MXLAY),
     &                  JS(MXLAY,NBANDS),JS1(MXLAY,NBANDS)
      COMMON /SELF/     SELFFAC(MXLAY), SELFFRAC(MXLAY), INDSELF(MXLAY)
      COMMON /K2/       KA(5,13,MG), KB(5,13:59,MG) , SELFREF(10,MG)

      EQUIVALENCE (KA,ABSA),(KB,ABSB)

      DIMENSION ABSA(65,MG),ABSB(235,MG)
      REAL KA,KB

C     Compute the optical depth by interpolating in ln(pressure) and 
C     temperature.  Below LAYTROP, the water vapor self-continuum is 
C     interpolated (in temperature) separately.
      DO 2500 LAY = 1, LAYTROP
         IND0 = ((JP(LAY)-1)*5+(JT(LAY)-1))*NSPA(2) + JS(LAY,2)
         IND1 = (JP(LAY)*5+(JT1(LAY)-1))*NSPA(2) + JS1(LAY,2)
         INDS = INDSELF(LAY)
         DO 2000 IG = 1, NG(2)
            TAUG(LAY,IG) = COLH2O(LAY) * 
     &          (FAC000(LAY,2) * ABSA(IND0,IG) +
     &           FAC010(LAY,2) * ABSA(IND0+1,IG) +
     &           FAC001(LAY,2) * ABSA(IND1,IG) + 
     &           FAC011(LAY,2) * ABSA(IND1+1,IG) + 
     &           SELFFAC(LAY) * (SELFREF(INDS,IG) + 
     &           SELFFRAC(LAY) *
     &           (SELFREF(INDS+1,IG) - SELFREF(INDS,IG))))
 2000    CONTINUE
 2500 CONTINUE
      DO 3500 LAY = LAYTROP+1, NLAYERS
         IND0 = ((JP(LAY)-13)*5+(JT(LAY)-1))*NSPB(2) + JS(LAY,2)
         IND1 = ((JP(LAY)-12)*5+(JT1(LAY)-1))*NSPB(2) + JS1(LAY,2)
         DO 3000 IG = 1, NG(2)
            TAUG(LAY,IG) = COLH2O(LAY) * 
     &          (FAC000(LAY,2) * ABSB(IND0,IG) +
     &           FAC010(LAY,2) * ABSB(IND0+1,IG) +
     &           FAC001(LAY,2) * ABSB(IND1,IG) + 
     &           FAC011(LAY,2) * ABSB(IND1+1,IG))  
 3000    CONTINUE
 3500 CONTINUE
      RETURN
      END

C----------------------------------------------------------------------------

      SUBROUTINE TAUGB3

C     BAND 3:  500-630 cm-1 (low - H2O,CO2; high - H2O,CO2)

      PARAMETER (MG=16, MXLAY=200, NBANDS=16)

C  Output

      COMMON /TAUGCOM/  TAUG(MXLAY,MG)

C  Input

      COMMON /FEATURES/ NG(NBANDS),NSPA(NBANDS),NSPB(NBANDS)
      COMMON /PROFILE/  NLAYERS,PAVEL(MXLAY),TAVEL(MXLAY),
     &                  PZ(0:MXLAY),TZ(0:MXLAY)
      COMMON /PROFDATA/ LAYTROP,COLH2O(MXLAY),COLCO2(MXLAY),
     &                  COLO3(MXLAY)
      COMMON /INTFAC/   FAC000(MXLAY,NBANDS),FAC001(MXLAY,NBANDS),
     &                  FAC010(MXLAY,NBANDS),FAC011(MXLAY,NBANDS),
     &                  FAC100(MXLAY,NBANDS),FAC101(MXLAY,NBANDS),
     &                  FAC110(MXLAY,NBANDS),FAC111(MXLAY,NBANDS)
      COMMON /INTIND/   JP(MXLAY),JT(MXLAY),JT1(MXLAY),
     &                  JS(MXLAY,NBANDS),JS1(MXLAY,NBANDS)
      COMMON /SELF/     SELFFAC(MXLAY), SELFFRAC(MXLAY), INDSELF(MXLAY)
      COMMON /K3/       KA(7,5,13,MG), KB(3,5,13:59,MG) , SELFREF(10,MG)

      EQUIVALENCE (KA,ABSA),(KB,ABSB)

      DIMENSION ABSA(455,MG),ABSB(705,MG)
      REAL KA,KB
      STRRAT31 = 1.403829 
      STRRAT32 = 1.403829

C     Compute the optical depth by interpolating in ln(pressure), 
C     temperature, and appropriate species.  Below LAYTROP, the water
C     vapor self-continuum is interpolated (in temperature) separately.  
      DO 2500 LAY = 1, LAYTROP
         IND0 = ((JP(LAY)-1)*5+(JT(LAY)-1))*NSPA(3) + JS(LAY,3)
         IND1 = (JP(LAY)*5+(JT1(LAY)-1))*NSPA(3) + JS1(LAY,3)
         INDS = INDSELF(LAY)
         DO 2000 IG = 1, NG(3)
            TAUG(LAY,IG) = (COLH2O(LAY) + STRRAT31*COLCO2(LAY)) * 
     &          (FAC000(LAY,3) * ABSA(IND0,IG) +
     &           FAC100(LAY,3) * ABSA(IND0+1,IG) +
     &           FAC010(LAY,3) * ABSA(IND0+7,IG) +
     &           FAC110(LAY,3) * ABSA(IND0+8,IG) +
     &           FAC001(LAY,3) * ABSA(IND1,IG) + 
     &           FAC101(LAY,3) * ABSA(IND1+1,IG) +
     &           FAC011(LAY,3) * ABSA(IND1+7,IG) +
     &           FAC111(LAY,3) * ABSA(IND1+8,IG)) +
     &           COLH2O(LAY) * 
     &           SELFFAC(LAY) * (SELFREF(INDS,IG) + 
     &           SELFFRAC(LAY) *
     &           (SELFREF(INDS+1,IG) - SELFREF(INDS,IG)))
 2000    CONTINUE
 2500 CONTINUE
      DO 3500 LAY = LAYTROP+1, NLAYERS
         IND0 = ((JP(LAY)-13)*5+(JT(LAY)-1))*NSPB(3) + JS(LAY,3)
         IND1 = ((JP(LAY)-12)*5+(JT1(LAY)-1))*NSPB(3) + JS1(LAY,3)
         DO 3000 IG = 1, NG(3)
            TAUG(LAY,IG) = (COLH2O(LAY) + STRRAT32*COLCO2(LAY)) * 
     &          (FAC000(LAY,3) * ABSB(IND0,IG) +
     &           FAC100(LAY,3) * ABSB(IND0+1,IG) +
     &           FAC010(LAY,3) * ABSB(IND0+3,IG) +
     &           FAC110(LAY,3) * ABSB(IND0+4,IG) +
     &           FAC001(LAY,3) * ABSB(IND1,IG) + 
     &           FAC101(LAY,3) * ABSB(IND1+1,IG) +
     &           FAC011(LAY,3) * ABSB(IND1+3,IG) +
     &           FAC111(LAY,3) * ABSB(IND1+4,IG))
 3000    CONTINUE
 3500 CONTINUE
      RETURN
      END

C----------------------------------------------------------------------------

      SUBROUTINE TAUGB4

C     BAND 4:  630-700 cm-1 (low - H2O,CO2; high - O3,CO2)

      PARAMETER (MG=16, MXLAY=200, NBANDS=16)

C  Output

      COMMON /TAUGCOM/  TAUG(MXLAY,MG)

C  Input

      COMMON /FEATURES/ NG(NBANDS),NSPA(NBANDS),NSPB(NBANDS)
      COMMON /PROFILE/  NLAYERS,PAVEL(MXLAY),TAVEL(MXLAY),
     &                  PZ(0:MXLAY),TZ(0:MXLAY)
      COMMON /PROFDATA/ LAYTROP,COLH2O(MXLAY),COLCO2(MXLAY),
     &                  COLO3(MXLAY)
      COMMON /INTFAC/   FAC000(MXLAY,NBANDS),FAC001(MXLAY,NBANDS),
     &                  FAC010(MXLAY,NBANDS),FAC011(MXLAY,NBANDS),
     &                  FAC100(MXLAY,NBANDS),FAC101(MXLAY,NBANDS),
     &                  FAC110(MXLAY,NBANDS),FAC111(MXLAY,NBANDS)
      COMMON /INTIND/   JP(MXLAY),JT(MXLAY),JT1(MXLAY),
     &                  JS(MXLAY,NBANDS),JS1(MXLAY,NBANDS)
      COMMON /SELF/     SELFFAC(MXLAY), SELFFRAC(MXLAY), INDSELF(MXLAY)
      COMMON /K4/       KA(7,5,13,MG), KB(3,5,13:59,MG) , SELFREF(10,MG)

      EQUIVALENCE (KA,ABSA),(KB,ABSB)

      DIMENSION ABSA(455,MG),ABSB(705,MG)
      REAL KA,KB
      STRRAT41 = 866.9486
      STRRAT42 = 32.959

C     Compute the optical depth by interpolating in ln(pressure),
C     temperature, and appropriate species.  Below LAYTROP, the water 
C     vapor self-continuum is interpolated (in temperature) separately.
      DO 2500 LAY = 1, LAYTROP
         IND0 = ((JP(LAY)-1)*5+(JT(LAY)-1))*NSPA(4) + JS(LAY,4)
         IND1 = (JP(LAY)*5+(JT1(LAY)-1))*NSPA(4) + JS1(LAY,4)
         INDS = INDSELF(LAY)
         DO 2000 IG = 1, NG(4)
            TAUG(LAY,IG) = (COLH2O(LAY) + STRRAT41*COLCO2(LAY)) * 
     &          (FAC000(LAY,4) * ABSA(IND0,IG) +
     &           FAC100(LAY,4) * ABSA(IND0+1,IG) +
     &           FAC010(LAY,4) * ABSA(IND0+7,IG) +
     &           FAC110(LAY,4) * ABSA(IND0+8,IG) +
     &           FAC001(LAY,4) * ABSA(IND1,IG) + 
     &           FAC101(LAY,4) * ABSA(IND1+1,IG) +
     &           FAC011(LAY,4) * ABSA(IND1+7,IG) +
     &           FAC111(LAY,4) * ABSA(IND1+8,IG)) +
     &           COLH2O(LAY) *
     &           SELFFAC(LAY) * (SELFREF(INDS,IG) + 
     &           SELFFRAC(LAY) *
     &           (SELFREF(INDS+1,IG) - SELFREF(INDS,IG)))
 2000    CONTINUE
 2500 CONTINUE
      DO 3500 LAY = LAYTROP+1, NLAYERS
         IND0 = ((JP(LAY)-13)*5+(JT(LAY)-1))*NSPB(4) + JS(LAY,4)
         IND1 = ((JP(LAY)-12)*5+(JT1(LAY)-1))*NSPB(4) + JS1(LAY,4)
         DO 3000 IG = 1, NG(4)
            TAUG(LAY,IG) = (COLO3(LAY) + STRRAT42*COLCO2(LAY)) * 
     &          (FAC000(LAY,4) * ABSB(IND0,IG) +
     &           FAC100(LAY,4) * ABSB(IND0+1,IG) +
     &           FAC010(LAY,4) * ABSB(IND0+3,IG) +
     &           FAC110(LAY,4) * ABSB(IND0+4,IG) +
     &           FAC001(LAY,4) * ABSB(IND1,IG) + 
     &           FAC101(LAY,4) * ABSB(IND1+1,IG) +
     &           FAC011(LAY,4) * ABSB(IND1+3,IG) +
     &           FAC111(LAY,4) * ABSB(IND1+4,IG))
 3000    CONTINUE 
 3500 CONTINUE
      RETURN
      END

C----------------------------------------------------------------------------

      SUBROUTINE TAUGB5

C     BAND 5:  700-820 cm-1 (low - H2O,CO2; high - O3,CO2)

      PARAMETER (MG=16, MXLAY=200, NBANDS=16)

C  Output

      COMMON /TAUGCOM/  TAUG(MXLAY,MG)

C  Input

      COMMON /FEATURES/ NG(NBANDS),NSPA(NBANDS),NSPB(NBANDS)
      COMMON /PROFILE/  NLAYERS,PAVEL(MXLAY),TAVEL(MXLAY),
     &                  PZ(0:MXLAY),TZ(0:MXLAY)
      COMMON /PROFDATA/ LAYTROP,COLH2O(MXLAY),COLCO2(MXLAY),
     &                  COLO3(MXLAY)
      COMMON /INTFAC/   FAC000(MXLAY,NBANDS),FAC001(MXLAY,NBANDS),
     &                  FAC010(MXLAY,NBANDS),FAC011(MXLAY,NBANDS),
     &                  FAC100(MXLAY,NBANDS),FAC101(MXLAY,NBANDS),
     &                  FAC110(MXLAY,NBANDS),FAC111(MXLAY,NBANDS)
      COMMON /INTIND/   JP(MXLAY),JT(MXLAY),JT1(MXLAY),
     &                  JS(MXLAY,NBANDS),JS1(MXLAY,NBANDS)
      COMMON /SELF/     SELFFAC(MXLAY), SELFFRAC(MXLAY), INDSELF(MXLAY)
      COMMON /K5/       KA(7,5,13,MG), KB(3,5,13:59,MG) , SELFREF(10,MG)

      EQUIVALENCE (KA,ABSA),(KB,ABSB)

      DIMENSION ABSA(455,MG),ABSB(705,MG)
      REAL KA,KB
      STRRAT51 = 86.9028
      STRRAT52 = 0.820692

C     Compute the optical depth by interpolating in ln(pressure), 
C     temperature, and appropriate species.  Below LAYTROP, the water 
C     vapor self-continuum is interpolated (in temperature) separately.
      DO 2500 LAY = 1, LAYTROP
         IND0 = ((JP(LAY)-1)*5+(JT(LAY)-1))*NSPA(5) + JS(LAY,5)
         IND1 = (JP(LAY)*5+(JT1(LAY)-1))*NSPA(5) + JS1(LAY,5)
         INDS = INDSELF(LAY)
         DO 2000 IG = 1, NG(5)
            TAUG(LAY,IG) = (COLH2O(LAY) + STRRAT51*COLCO2(LAY)) * 
     &          (FAC000(LAY,5) * ABSA(IND0,IG) +
     &           FAC100(LAY,5) * ABSA(IND0+1,IG) +
     &           FAC010(LAY,5) * ABSA(IND0+7,IG) +
     &           FAC110(LAY,5) * ABSA(IND0+8,IG) +
     &           FAC001(LAY,5) * ABSA(IND1,IG) + 
     &           FAC101(LAY,5) * ABSA(IND1+1,IG) +
     &           FAC011(LAY,5) * ABSA(IND1+7,IG) +
     &           FAC111(LAY,5) * ABSA(IND1+8,IG)) +
     &           COLH2O(LAY) *
     &           SELFFAC(LAY) * (SELFREF(INDS,IG) + 
     &           SELFFRAC(LAY) *
     &           (SELFREF(INDS+1,IG) - SELFREF(INDS,IG)))
 2000    CONTINUE
 2500 CONTINUE
      DO 3500 LAY = LAYTROP+1, NLAYERS
         IND0 = ((JP(LAY)-13)*5+(JT(LAY)-1))*NSPB(5) + JS(LAY,5)
         IND1 = ((JP(LAY)-12)*5+(JT1(LAY)-1))*NSPB(5) + JS1(LAY,5)
         DO 3000 IG = 1, NG(5)
            TAUG(LAY,IG) = (COLO3(LAY) + STRRAT52*COLCO2(LAY)) * 
     &          (FAC000(LAY,5) * ABSB(IND0,IG) +
     &           FAC100(LAY,5) * ABSB(IND0+1,IG) +
     &           FAC010(LAY,5) * ABSB(IND0+3,IG) +
     &           FAC110(LAY,5) * ABSB(IND0+4,IG) +
     &           FAC001(LAY,5) * ABSB(IND1,IG) + 
     &           FAC101(LAY,5) * ABSB(IND1+1,IG) +
     &           FAC011(LAY,5) * ABSB(IND1+3,IG) +
     &           FAC111(LAY,5) * ABSB(IND1+4,IG))
 3000    CONTINUE
 3500 CONTINUE
      RETURN
      END

C----------------------------------------------------------------------------

      SUBROUTINE TAUGB6

C     BAND 6:  820-980 cm-1 (low - H2O; high - nothing)

      PARAMETER (MG=16, MXLAY=200, NBANDS=16)

C  Output

      COMMON /TAUGCOM/  TAUG(MXLAY,MG)

C  Input

      COMMON /FEATURES/ NG(NBANDS),NSPA(NBANDS),NSPB(NBANDS)
      COMMON /PROFILE/  NLAYERS,PAVEL(MXLAY),TAVEL(MXLAY),
     &                  PZ(0:MXLAY),TZ(0:MXLAY)
      COMMON /PROFDATA/ LAYTROP,COLH2O(MXLAY),COLCO2(MXLAY),
     &                  COLO3(MXLAY)
      COMMON /INTFAC/   FAC000(MXLAY,NBANDS),FAC001(MXLAY,NBANDS),
     &                  FAC010(MXLAY,NBANDS),FAC011(MXLAY,NBANDS),
     &                  FAC100(MXLAY,NBANDS),FAC101(MXLAY,NBANDS),
     &                  FAC110(MXLAY,NBANDS),FAC111(MXLAY,NBANDS)
      COMMON /INTIND/   JP(MXLAY),JT(MXLAY),JT1(MXLAY),
     &                  JS(MXLAY,NBANDS),JS1(MXLAY,NBANDS)
      COMMON /SELF/     SELFFAC(MXLAY), SELFFRAC(MXLAY), INDSELF(MXLAY)
      COMMON /K6/       KA(5,13,MG), SELFREF(10,MG)

      EQUIVALENCE (KA,ABSA),(KB,ABSB)

      DIMENSION ABSA(65,MG)
      REAL KA,KB

C     Compute the optical depth by interpolating in ln(pressure) and
C     temperature. The water vapor self-continuum is interpolated
C     (in temperature) separately.  
      DO 2500 LAY = 1, LAYTROP
         IND0 = ((JP(LAY)-1)*5+(JT(LAY)-1))*NSPA(6) + JS(LAY,6)
         IND1 = (JP(LAY)*5+(JT1(LAY)-1))*NSPA(6) + JS1(LAY,6)
         INDS = INDSELF(LAY)
         DO 2000 IG = 1, NG(6)
            TAUG(LAY,IG) = COLH2O(LAY) * 
     &          (FAC000(LAY,6) * ABSA(IND0,IG) +
     &           FAC010(LAY,6) * ABSA(IND0+1,IG) +
     &           FAC001(LAY,6) * ABSA(IND1,IG) + 
     &           FAC011(LAY,6) * ABSA(IND1+1,IG) + 
     &           SELFFAC(LAY) * (SELFREF(INDS,IG) + 
     &           SELFFRAC(LAY)*
     &           (SELFREF(INDS+1,IG)-SELFREF(INDS,IG))))
 2000    CONTINUE
 2500 CONTINUE
C     Nothing important goes on above LAYTROP in this band.
      DO 3500 LAY = LAYTROP+1, NLAYERS
         DO 3000 IG = 1, NG(6)
            TAUG(LAY,IG) = 0.0
 3000    CONTINUE
 3500 CONTINUE
      RETURN
      END

