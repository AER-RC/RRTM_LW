C     path:      %P%
C     revision:  $Revision$
C     created:   $Date$  
C     presently: %H%  %T%
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

      PARAMETER (MXLAY=203)
      PARAMETER (MG = 16)
      PARAMETER (NBANDS = 16)

      CHARACTER*8 HVRRTM,HVRINI,HVRRT0,HVRATM,HVRSET,HVRTAU,
     *            HVDUM1,HVRUTL,HVREXT
      CHARACTER*8 HVRKG

      COMMON /FEATURES/  NG(NBANDS),NSPA(MG),NSPB(MG)
      COMMON /CONTROL/   NUMANGS
      COMMON /PROFILE/   NLAYERS,PAVEL(MXLAY),TAVEL(MXLAY),
     &                   PZ(0:MXLAY),TZ(0:MXLAY),TSFC
      COMMON /OUTPUT/    TOTUFLUX(0:MXLAY), TOTDFLUX(0:MXLAY),
     &                   FNET(0:MXLAY), HTR(0:MXLAY)
      COMMON /CONSTANTS/ PI,FLUXFAC,HEATFAC
      COMMON /HVERSN/ HVRRTM,HVRINI,HVRRT0,HVRATM,HVRSET,HVRTAU,
     *                HVDUM1(4),HVRUTL,HVREXT
      COMMON /HVRSNB/ HVRKG(NBANDS)


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
         OPEN (IWR,FILE='OUTPUT_RRTM',FORM='FORMATTED')

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
         WRITE(IWR,9900)
         WRITE(IWR,9901)
C
         DO 3000 I = NLAYERS, 0, -1
            WRITE(IWR,9902) I, PZ(I), TOTUFLUX(I), TOTDFLUX(I),
     &           FNET(I), HTR(I)
 3000    CONTINUE

C
C ***    Output module version numbers
C
         WRITE(IWR,9910) HVRRTM,HVRINI,HVRRT0,HVRATM,HVRSET,HVRTAU,
     *                   HVRUTL,HVREXT,(HVRKG(NB),NB=1,2)

         CLOSE(IWR)


 9900 FORMAT(1X,'LEVEL    PRESSURE   UPWARD FLUX   DOWNWARD FLUX    NET    
     &FLUX       HEATING RATE')
 9901 FORMAT(1X,'            mb          W/m2          W/m2           W/
     &m2          degree/day')
 9902 FORMAT(1X,I3,3X,F11.6,4X,1P,2(G12.6,2X),G13.6,3X,G16.9,0P)
 9910 FORMAT('  Modules and versions used in this calculation:',/,/,5X,
     *        '    rrtm.f: ',6X,A8,10X, '  init.f: ',6X,A8,/,5X,
     *        '      rt.f: ',6X,A8,10X, 'rrtatm.f: ',6X,A8,/,5X,
     *        ' setcoef.f: ',6X,A8,10X, 'taumol.f: ',6X,A8,/,5X,
     *        'util_xxx.f: ',6X,A8,10X, ' extra.f: ',6X,A8,/,5X,
     *        '   k_gB1.f: ',6X,A8,10X, ' k_gB2.f: ',6X,A8,/,5X)
c    *        '   k_gB3.f: ',6X,A8,10X, ' k_gB4.f: ',6X,A8,/,5X,
c    *        '   k_gB5.f: ',6X,A8,10X, ' k_gB6.f: ',6X,A8,/,5X,
c    *        '   k_gB7.f: ',6X,A8,10X, ' k_gB8.f: ',6X,A8,/,5X,
c    *        '   k_gB9.f: ',6X,A8,10X, 'k_gB10.f: ',6X,A8,/,5X,
c    *        '  k_gB11.f: ',6X,A8,10X, 'k_gB12.f: ',6X,A8,/,5X,
c    *        '  k_gB13.f: ',6X,A8,10X, 'k_gB14.f: ',6X,A8,/,5X,
c    *        '  k_gB15.f: ',6X,A8,10X, 'k_gB16.f: ',6X,A8,/)


      STOP
      END

C************************  SUBROUTINE READPROF  *****************************C

      SUBROUTINE READPROF                                                     
                                                                         
      IMPLICIT DOUBLE PRECISION (V)                                      
                                                                         
      PARAMETER (MXLAY=203)
      DIMENSION ALTZ(0:MXLAY)


      COMMON /IFIL/ IRD,IPR,IPU,DUM(15)

      COMMON /PROFILE/ NLAYERS,PAVEL(MXLAY),TAVEL(MXLAY),
     &                 PZ(0:MXLAY),TZ(0:MXLAY),TSFC
      COMMON /SPECIES/ COLDRY(MXLAY),WKL(35,MXLAY),WBRODL(MXLAY),NMOL
      COMMON /MOLAMNT/DRYAIR(MXLAY),WH2O(MXLAY),WCO2(MXLAY),WO3(MXLAY)

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
      OPEN (IRD,FILE='TAPE5',FORM='FORMATTED')

 1000 CONTINUE
      READ (IRD,9010,END=8800) CTEST
      IF (CTEST .NE. CDOLLAR) GO TO 1000

      READ(IRD,9011) IATM
      READ (IRD,9012) TSFC

      IF (IATM .EQ. 0) THEN
         READ (IRD,9013) IFORM,NLAYERS,NMOL
         IF (NMOL.EQ.0) NMOL = 7                                    
         READ (IRD,FORM1(IFORM)) PAVEL(1),TAVEL(1),SECNTK,CINP,
     &        IPTHAK,ALTZ(0),PZ(0),TZ(0),ALTZ(1),PZ(1),TZ(1)
         READ (IRD,FORM3(IFORM)) (WKL(M,1),M=1,7), WBRODL(1)
         IF(NMOL .GT. 7) READ (IRD,FORM3(IFORM)) (WKL(M,1),M=8,NMOL)

         SUMMOL = 0.0
         DO 2000 IMOL = 2, NMOL
            SUMMOL = SUMMOL + WKL(IMOL,1)
 2000    CONTINUE

C        Test for mixing ratio input.
         IMIX = 0
         IF (WKL(2,1) .LE. 1.0) IMIX = 1
         IF (IMIX .EQ. 1) THEN
            DRAIR = WBRODL(1) / (1. - SUMMOL)
            DO 2500 IMOL = 1, NMOL
               WKL(IMOL,1) = DRAIR * WKL(IMOL,1)
 2500       CONTINUE
         ENDIF

      DO 4000 L = 2, NLAYERS
           READ (IRD,FORM2(IFORM)) PAVEL(L),TAVEL(L),SECNTK,CINP,
     &        IPTHRK,ALTZ(L),PZ(L),TZ(L)
           READ (IRD,FORM3(IFORM)) (WKL(M,L),M=1,7), WBRODL(L)
           IF(NMOL .GT. 7) READ (IRD,FORM3(IFORM)) (WKL(M,L),M=8,NMOL)

           SUMMOL = 0.0
           DO 3000 IMOL = 2, NMOL
              SUMMOL = SUMMOL + WKL(IMOL,L)
 3000      CONTINUE
           
           IF (IMIX .EQ. 1) THEN
              DRAIR = WBRODL(L) / (1. - SUMMOL)
              DO 3500 IMOL = 1, NMOL
                 WKL(IMOL,L) = DRAIR * WKL(IMOL,L)
 3500         CONTINUE
           ENDIF
 4000   CONTINUE                                                            

      ELSE
         IPU = 7
         IPR = 66
         OPEN(UNIT=IPR,FILE='TAPE6',STATUS='UNKNOWN')
         CALL RRTATM
      ENDIF

      CLOSE(IRD)
      GO TO 9000

 8800 CONTINUE
      STOP ' INVALID INPUT_RRTM '

 9000 CONTINUE
c
c
      DO 100 L = 1, NLAYERS                                                 
         COLDRY(L) = WBRODL(L)
         DO 90 M = 2,NMOL
            COLDRY(L) = COLDRY(L) + WKL(M,L)
 90      CONTINUE
         DRYAIR(L) = COLDRY(L)
         WH2O(L) = WKL(1,L)/COLDRY(L)
         WCO2(L) = WKL(2,L)/COLDRY(L)
         WO3(L)  = WKL(3,L)/COLDRY(L)
 100  CONTINUE                                                            
c
 9010 FORMAT (A1)
 9011 FORMAT (49X,I1)
 9012 FORMAT (E10.3)
 9013 FORMAT (1X,I1,I3,I5)                                     

      RETURN
      END 

      BLOCK DATA


      COMMON /HVERSN/ HVRRTM,HVRINI,HVRRT0,HVRATM,HVRSET,HVRTAU,
     *                HVDUM1(4),HVRUTL,HVREXT

      CHARACTER*8 HVRRTM,HVRINI,HVRRT0,HVRATM,HVRSET,HVRTAU,
     *            HVDUM1,HVRUTL,HVREXT

      DATA HVRRTM / '$Revision$' /,        HVRINI / 'NOT USED' /,
     *     HVRRT0 / 'NOT USED' /,   HVRATM / 'NOT USED' /,
     *     HVRSET / 'NOT USED' /,   HVRTAU / 'NOT USED' /,
     *     HVDUM1 / 4*'NOT USED' /, HVRUTL / 'NOT USED' /,
     *     HVREXT / 'NOT USED' /

      END



