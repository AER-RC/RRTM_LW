C     path:      %P%
C     revision:  $Revision$
C     created:   $Date$  
C     presently: %H%  %T%
      SUBROUTINE SETCOEF

C     Purpose:  For a given atmosphere, calculate the indices and
C     fractions needed to interpolate between absorption coefficients 
C     stored for a variety of reference atmospheres.

      PARAMETER (MXLAY = 203)
      PARAMETER (NBANDS = 16)
      PARAMETER (MG =16)

C  Input      
      COMMON /PROFILE/  NLAYERS,PAVEL(MXLAY),TAVEL(MXLAY),
     &                  PZ(0:MXLAY),TZ(0:MXLAY),TSFC
      COMMON /MOLAMNT/  COLDRY(MXLAY),WATER(MXLAY),CO2(MXLAY),
     &                  OZ(MXLAY)

C  Output
      COMMON /PROFDATA/ LAYTROP,COLH2O(MXLAY),COLCO2(MXLAY),
     &                  COLO3(MXLAY)
      COMMON /FRACS/    TLEVFRAC(0:MXLAY),TLAYFRAC(MXLAY)
      COMMON /MISC/     INDLEV(0:MXLAY),INDLAY(MXLAY)
      COMMON /INTFAC/   FAC000(MXLAY,NBANDS),FAC001(MXLAY,NBANDS),
     &                  FAC010(MXLAY,NBANDS),FAC011(MXLAY,NBANDS),
     &                  FAC100(MXLAY,NBANDS),FAC101(MXLAY,NBANDS),
     &                  FAC110(MXLAY,NBANDS),FAC111(MXLAY,NBANDS)
      COMMON /INTIND/   JP(MXLAY),JT(MXLAY),JT1(MXLAY),
     &                  JS(MXLAY,NBANDS),JS1(MXLAY,NBANDS)
      COMMON /SELF/     SELFFAC, SELFFRAC, INDSELF

      CHARACTER*30 A1
      CHARACTER*10 A2

      DIMENSION SELFFAC(MXLAY),SELFFRAC(MXLAY),INDSELF(MXLAY)
      DIMENSION PREF(59),PREFLOG(59),TREF(59)
      DIMENSION SPARM31(7,13),SPARM32(3,13:59)
      DIMENSION SPARM41(7,13),SPARM42(3,13:59)
      DIMENSION SPARM51(7,13),SPARM52(3,13:59)

C     These pressures are chosen such that the ln of the first pressure
C     has only a few non-zero digits (i.e. ln(PREF(1)) = 6.96000) and
C     each subsequent ln(pressure) differs from the previous one by 0.2.
      DATA PREF /
     &    1.05363E+03,8.62642E+02,7.06272E+02,5.78246E+02,4.73428E+02,
     &    3.87610E+02,3.17348E+02,2.59823E+02,2.12725E+02,1.74164E+02,
     &    1.42594E+02,1.16746E+02,9.55835E+01,7.82571E+01,6.40715E+01,
     &    5.24573E+01,4.29484E+01,3.51632E+01,2.87892E+01,2.35706E+01,
     &    1.92980E+01,1.57998E+01,1.29358E+01,1.05910E+01,8.67114E+00,
     &    7.09933E+00,5.81244E+00,4.75882E+00,3.89619E+00,3.18993E+00,
     &    2.61170E+00,2.13828E+00,1.75067E+00,1.43333E+00,1.17351E+00,
     &    9.60789E-01,7.86628E-01,6.44036E-01,5.27292E-01,4.31710E-01,
     &    3.53455E-01,2.89384E-01,2.36928E-01,1.93980E-01,1.58817E-01,
     &    1.30029E-01,1.06458E-01,8.71608E-02,7.13612E-02,5.84256E-02,
     &    4.78349E-02,3.91639E-02,3.20647E-02,2.62523E-02,2.14936E-02,
     &    1.75975E-02,1.44076E-02,1.17959E-02,9.65769E-03/
      DATA PREFLOG /
     &     6.9600E+00, 6.7600E+00, 6.5600E+00, 6.3600E+00, 6.1600E+00,
     &     5.9600E+00, 5.7600E+00, 5.5600E+00, 5.3600E+00, 5.1600E+00,
     &     4.9600E+00, 4.7600E+00, 4.5600E+00, 4.3600E+00, 4.1600E+00,
     &     3.9600E+00, 3.7600E+00, 3.5600E+00, 3.3600E+00, 3.1600E+00,
     &     2.9600E+00, 2.7600E+00, 2.5600E+00, 2.3600E+00, 2.1600E+00,
     &     1.9600E+00, 1.7600E+00, 1.5600E+00, 1.3600E+00, 1.1600E+00,
     &     9.6000E-01, 7.6000E-01, 5.6000E-01, 3.6000E-01, 1.6000E-01,
     &    -4.0000E-02,-2.4000E-01,-4.4000E-01,-6.4000E-01,-8.4000E-01,
     &    -1.0400E+00,-1.2400E+00,-1.4400E+00,-1.6400E+00,-1.8400E+00,
     &    -2.0400E+00,-2.2400E+00,-2.4400E+00,-2.6400E+00,-2.8400E+00,
     &    -3.0400E+00,-3.2400E+00,-3.4400E+00,-3.6400E+00,-3.8400E+00,
     &    -4.0400E+00,-4.2400E+00,-4.4400E+00,-4.6400E+00/
C     These are the temperatures associated with the respective 
C     pressures for the MLS standard atmosphere. 
      DATA TREF /
     &     2.9420E+02, 2.8799E+02, 2.7894E+02, 2.6925E+02, 2.5983E+02,
     &     2.5017E+02, 2.4077E+02, 2.3179E+02, 2.2306E+02, 2.1578E+02,
     &     2.1570E+02, 2.1570E+02, 2.1570E+02, 2.1706E+02, 2.1858E+02,
     &     2.2018E+02, 2.2174E+02, 2.2328E+02, 2.2479E+02, 2.2655E+02,
     &     2.2834E+02, 2.3113E+02, 2.3401E+02, 2.3703E+02, 2.4022E+02,
     &     2.4371E+02, 2.4726E+02, 2.5085E+02, 2.5457E+02, 2.5832E+02,
     &     2.6216E+02, 2.6606E+02, 2.6999E+02, 2.7340E+02, 2.7536E+02,
     &     2.7568E+02, 2.7372E+02, 2.7163E+02, 2.6955E+02, 2.6593E+02,
     &     2.6211E+02, 2.5828E+02, 2.5360E+02, 2.4854E+02, 2.4348E+02,
     &     2.3809E+02, 2.3206E+02, 2.2603E+02, 2.2000E+02, 2.1435E+02,
     &     2.0887E+02, 2.0340E+02, 1.9792E+02, 1.9290E+02, 1.8809E+02,
     &     1.8329E+02, 1.7849E+02, 1.7394E+02, 1.7212E+02/

C     BAND 1 (10-250 cm-1):
C     Since there is only one key specie in this band (H2O), there is 
C     no need for a binary species parameter as in other bands.

C     BAND 2 (250-500 cm-1):
C     Since there is only one key specie in this band (H2O), there is 
C     no need for a binary species parameter as in other bands.

C     BAND 3 (500-630 cm-1):
C     These are the 7 values of the binary species parameter for each
C     reference layer in the lower atmosphere.  The parameter is 
C     defined by H2O/(H2O + STRRAT31*CO2).  
      DATA SPARM31 / 
     &0.0,.92399E+0,.96050E+0,.97331E+0,.97985E+0,.98381E+0,.98648E+0,
     &0.0,.90909E+0,.95238E+0,.96774E+0,.97561E+0,.98039E+0,.98361E+0,
     &0.0,.86888E+0,.92984E+0,.95211E+0,.96364E+0,.97070E+0,.97547E+0,
     &0.0,.80014E+0,.88898E+0,.92314E+0,.94122E+0,.95242E+0,.96003E+0,
     &0.0,.70018E+0,.82365E+0,.87509E+0,.90330E+0,.92111E+0,.93339E+0,
     &0.0,.55649E+0,.71505E+0,.79010E+0,.83386E+0,.86252E+0,.88274E+0,
     &0.0,.39084E+0,.56202E+0,.65810E+0,.71961E+0,.76236E+0,.79380E+0,
     &0.0,.24065E+0,.38794E+0,.48737E+0,.55901E+0,.61309E+0,.65535E+0,
     &0.0,.12898E+0,.22849E+0,.30760E+0,.37199E+0,.42542E+0,.47048E+0,
     &0.0,.54244E-1,.10291E+0,.14681E+0,.18661E+0,.22286E+0,.25603E+0,
     &0.0,.34384E-1,.66482E-1,.96514E-1,.12467E+0,.15113E+0,.17604E+0,
     &0.0,.87717E-2,.17391E-1,.25861E-1,.34187E-1,.42372E-1,.50419E-1,
     &0.0,.63879E-2,.12695E-1,.18922E-1,.25071E-1,.31144E-1,.37141E-1/
C     &0.0,.37805E+0,.54867E+0,.64583E+0,.70857E+0,.75242E+0,.78481E+0,
C     &0.0,.33334E+0,.50001E+0,.60001E+0,.66668E+0,.71429E+0,.75001E+0,
C     &0.0,.24887E+0,.39855E+0,.49849E+0,.56995E+0,.62358E+0,.66532E+0,
C     &0.0,.16679E+0,.28589E+0,.37521E+0,.44466E+0,.50022E+0,.54567E+0,
C     &0.0,.10456E+0,.18932E+0,.25942E+0,.31837E+0,.36862E+0,.41197E+0,
C     &0.0,.59032E-1,.11148E+0,.15840E+0,.20060E+0,.23878E+0,.27348E+0,
C     &0.0,.31083E-1,.60292E-1,.87791E-1,.11373E+0,.13823E+0,.16141E+0,
C     &0.0,.15598E-1,.30718E-1,.45379E-1,.59604E-1,.73411E-1,.86819E-1,
C     &0.0,.73497E-2,.14592E-1,.21730E-1,.28764E-1,.35699E-1,.42535E-1,
C     &0.0,.28596E-2,.57028E-2,.85299E-2,.11341E-1,.14136E-1,.16916E-1,
C     &0.0,.17772E-2,.35482E-2,.53129E-2,.70713E-2,.88235E-2,.10570E-1,
C     &0.0,.44227E-3,.88415E-3,.13256E-2,.17667E-2,.22074E-2,.26478E-2,
C     &0.0,.32135E-3,.64248E-3,.96342E-3,.12841E-2,.16047E-2,.19250E-2/

C     These are the 3 values of the binary species parameter for each
C     reference layer in the upper atmosphere.  The parameter is defined 
C     by H2O/(H2O + STRRAT32*CO2).  
      DATA SPARM32 / 
     &0.0,.63879E-02,.12695E-01,0.0,.63045E-02,.12530E-01,
     &0.0,.64838E-02,.12884E-01,0.0,.68209E-02,.13549E-01,
     &0.0,.72289E-02,.14354E-01,0.0,.77942E-02,.15468E-01,
     &0.0,.82561E-02,.16377E-01,0.0,.85705E-02,.16995E-01,
     &0.0,.88346E-02,.17514E-01,0.0,.91022E-02,.18040E-01,
     &0.0,.93599E-02,.18546E-01,0.0,.95285E-02,.18877E-01,
     &0.0,.96769E-02,.19168E-01,0.0,.97877E-02,.19386E-01,
     &0.0,.98681E-02,.19543E-01,0.0,.99259E-02,.19657E-01,
     &0.0,.10037E-01,.19874E-01,0.0,.10182E-01,.20159E-01,
     &0.0,.10425E-01,.20635E-01,0.0,.10635E-01,.21046E-01,
     &0.0,.10819E-01,.21407E-01,0.0,.10882E-01,.21531E-01,
     &0.0,.10916E-01,.21596E-01,0.0,.10916E-01,.21596E-01,
     &0.0,.10825E-01,.21417E-01,0.0,.10729E-01,.21229E-01,
     &0.0,.10633E-01,.21041E-01,0.0,.10431E-01,.20647E-01,
     &0.0,.10216E-01,.20225E-01,0.0,.10000E-01,.19802E-01,
     &0.0,.96906E-02,.19195E-01,0.0,.93386E-02,.18504E-01,
     &0.0,.89865E-02,.17813E-01,0.0,.86255E-02,.17104E-01,
     &0.0,.82470E-02,.16359E-01,0.0,.78683E-02,.15614E-01,
     &0.0,.74892E-02,.14867E-01,0.0,.71169E-02,.14133E-01,
     &0.0,.67474E-02,.13404E-01,0.0,.63778E-02,.12675E-01,
     &0.0,.60078E-02,.11944E-01,0.0,.56437E-02,.11224E-01,
     &0.0,.52816E-02,.10508E-01,0.0,.49184E-02,.97886E-02,
     &0.0,.45540E-02,.90667E-02,0.0,.41951E-02,.83551E-02,
     &0.0,.39044E-02,.77785E-02/

C     BAND 4 (630-700 cm-1):
C     These are the 7 values of the binary species parameter for each
C     reference layer in the lower atmosphere.  The parameter is defined 
C     by H2O/(H2O + STRRAT41*CO2).  
      DATA SPARM41 / 
     &0.0,.19305E-1,.37879E-1,.55762E-1,.72993E-1,.89606E-1,.10563E+0,
     &0.0,.15935E-1,.31371E-1,.46330E-1,.60833E-1,.74903E-1,.88557E-1,
     &0.0,.10616E-1,.21010E-1,.31187E-1,.41155E-1,.50919E-1,.60487E-1,
     &0.0,.64410E-2,.12800E-1,.19077E-1,.25276E-1,.31396E-1,.37440E-1,
     &0.0,.37673E-2,.75062E-2,.11217E-1,.14901E-1,.18557E-1,.22186E-1,
     &0.0,.20276E-2,.40470E-2,.60583E-2,.80614E-2,.10057E-1,.12044E-1,
     &0.0,.10379E-2,.20735E-2,.31071E-2,.41385E-2,.51678E-2,.61950E-2,
     &0.0,.51290E-3,.10253E-2,.15371E-2,.20485E-2,.25593E-2,.30695E-2,
     &0.0,.23973E-3,.47934E-3,.71884E-3,.95822E-3,.11975E-2,.14366E-2,
     &0.0,.92865E-4,.18571E-3,.27854E-3,.37136E-3,.46415E-3,.55693E-3,
     &0.0,.57656E-4,.11531E-3,.17295E-3,.23058E-3,.28821E-3,.34584E-3,
     &0.0,.14329E-4,.28658E-4,.42987E-4,.57315E-4,.71642E-4,.85969E-4,
     &0.0,.10410E-4,.20820E-4,.31230E-4,.41639E-4,.52049E-4,.62458E-4/

C     These are the 3 values of the binary species parameter for each
C     reference layer in the upper atmosphere.  The parameter is defined 
C     by O3/(O3 + STRRAT42*CO2).  
      DATA SPARM42 / 
     &0.0,.59487E-04,.11897E-03,0.0,.95597E-04,.19118E-03,
     &0.0,.15055E-03,.30105E-03,0.0,.19883E-03,.39759E-03,
     &0.0,.25272E-03,.50531E-03,0.0,.31265E-03,.62511E-03,
     &0.0,.39257E-03,.78483E-03,0.0,.45439E-03,.90836E-03,
     &0.0,.50928E-03,.10180E-02,0.0,.55619E-03,.11118E-02,
     &0.0,.60333E-03,.12059E-02,0.0,.65695E-03,.13130E-02,
     &0.0,.70526E-03,.14095E-02,0.0,.74371E-03,.14863E-02,
     &0.0,.75431E-03,.15075E-02,0.0,.74428E-03,.14874E-02,
     &0.0,.69132E-03,.13817E-02,0.0,.62614E-03,.12515E-02,
     &0.0,.53902E-03,.10775E-02,0.0,.45851E-03,.91659E-03,
     &0.0,.38299E-03,.76569E-03,0.0,.32801E-03,.65581E-03,
     &0.0,.28048E-03,.56081E-03,0.0,.24126E-03,.48240E-03,
     &0.0,.21282E-03,.42555E-03,0.0,.18496E-03,.36986E-03,
     &0.0,.15710E-03,.31416E-03,0.0,.14201E-03,.28398E-03,
     &0.0,.12863E-03,.25722E-03,0.0,.11524E-03,.23046E-03,
     &0.0,.10231E-03,.20460E-03,0.0,.89581E-04,.17915E-03,
     &0.0,.76851E-04,.15369E-03,0.0,.65244E-04,.13048E-03,
     &0.0,.55876E-04,.11175E-03,0.0,.46508E-04,.93011E-04,
     &0.0,.37140E-04,.74276E-04,0.0,.31127E-04,.62252E-04,
     &0.0,.26660E-04,.53318E-04,0.0,.22192E-04,.44384E-04,
     &0.0,.17725E-04,.35450E-04,0.0,.16376E-04,.32752E-04,
     &0.0,.16583E-04,.33166E-04,0.0,.16791E-04,.33582E-04,
     &0.0,.16999E-04,.33998E-04,0.0,.17752E-04,.35504E-04,
     &0.0,.24351E-04,.48700E-04/

C     BAND 5 (700-820cm-1):
C     These are the 7 values of the binary species parameter for each
C     reference layer in the lower atmosphere.  The parameter is defined 
C     by H2O/(H2O + STRRAT51*CO2).  
      DATA SPARM51 / 
     &0.0,.16414E+0,.28200E+0,.37073E+0,.43994E+0,.49543E+0,.54092E+0,
     &0.0,.13908E+0,.24420E+0,.32644E+0,.39254E+0,.44682E+0,.49220E+0,
     &0.0,.96695E-1,.17634E+0,.24308E+0,.29981E+0,.34863E+0,.39109E+0,
     &0.0,.60744E-1,.11453E+0,.16249E+0,.20552E+0,.24435E+0,.27956E+0,
     &0.0,.36353E-1,.70156E-1,.10167E+0,.13111E+0,.15869E+0,.18457E+0,
     &0.0,.19866E-1,.38958E-1,.57321E-1,.74995E-1,.92018E-1,.10843E+0,
     &0.0,.10258E-1,.20308E-1,.30156E-1,.39807E-1,.49269E-1,.58546E-1,
     &0.0,.50933E-2,.10135E-1,.15126E-1,.20067E-1,.24958E-1,.29801E-1,
     &0.0,.23864E-2,.47614E-2,.71252E-2,.94777E-2,.11819E-1,.14150E-1,
     &0.0,.92566E-3,.18496E-2,.27718E-2,.36924E-2,.46112E-2,.55284E-2,
     &0.0,.57488E-3,.11491E-2,.17227E-2,.22956E-2,.28678E-2,.34394E-2,
     &0.0,.14293E-3,.28582E-3,.42867E-3,.57148E-3,.71425E-3,.85697E-3,
     &0.0,.10384E-3,.20766E-3,.31146E-3,.41524E-3,.51900E-3,.62273E-3/

C     These are the 3 values of the binary species parameter for each
C     reference layer in the upper atmosphere.  The parameter is defined 
C     by O3/(O3 + STRRAT52*CO2).  
      DATA SPARM52 / 
     &0.0,.23834E-02,.47556E-02,0.0,.38249E-02,.76206E-02,
     &0.0,.60106E-02,.11949E-01,0.0,.79234E-02,.15722E-01,
     &0.0,.10050E-01,.19900E-01,0.0,.12404E-01,.24504E-01,
     &0.0,.15527E-01,.30579E-01,0.0,.17929E-01,.35227E-01,
     &0.0,.20053E-01,.39317E-01,0.0,.21861E-01,.42786E-01,
     &0.0,.23671E-01,.46246E-01,0.0,.25721E-01,.50153E-01,
     &0.0,.27562E-01,.53646E-01,0.0,.29022E-01,.56408E-01,
     &0.0,.29424E-01,.57166E-01,0.0,.29044E-01,.56448E-01,
     &0.0,.27032E-01,.52640E-01,0.0,.24544E-01,.47912E-01,
     &0.0,.21200E-01,.41519E-01,0.0,.18089E-01,.35535E-01,
     &0.0,.15154E-01,.29855E-01,0.0,.13006E-01,.25678E-01,
     &0.0,.11142E-01,.22038E-01,0.0,.95983E-02,.19014E-01,
     &0.0,.84762E-02,.16810E-01,0.0,.73747E-02,.14641E-01,
     &0.0,.62707E-02,.12463E-01,0.0,.56716E-02,.11279E-01,
     &0.0,.51397E-02,.10227E-01,0.0,.46073E-02,.91723E-02,
     &0.0,.40923E-02,.81513E-02,0.0,.35850E-02,.71443E-02,
     &0.0,.30771E-02,.61353E-02,0.0,.26135E-02,.52134E-02,
     &0.0,.22391E-02,.44681E-02,0.0,.18644E-02,.37218E-02,
     &0.0,.14894E-02,.29743E-02,0.0,.12485E-02,.24939E-02,
     &0.0,.10695E-02,.21368E-02,0.0,.89047E-03,.17794E-02,
     &0.0,.71135E-03,.14217E-02,0.0,.65725E-03,.13136E-02,
     &0.0,.66556E-03,.13302E-02,0.0,.67389E-03,.13469E-02,
     &0.0,.68223E-03,.13635E-02,0.0,.71244E-03,.14239E-02,
     &0.0,.97699E-03,.19521E-02/

C     BAND 6 (820-980 cm-1):
C     Since there is at most one key specie in this band, there is no
C     need for for a ratio between two key species as in other bands.

C ****************** START OF EXECUTABLE CODE ***************************

      STPFAC = 296./1013.

C     The following are the ratios of the line strengths in the 
C     respective bands of that band's two key species.
      STRRAT31 = 1.403829 
      STRRAT32 = 1.403829
      STRRAT41 = 866.9486
      STRRAT42 = 32.959
      STRRAT51 = 86.9028
      STRRAT52 = 0.820692

C     Since the output absorption coefficient might not be too accurate
C     if it resulted from extrapolation (in any of the variables) rather
C     than interpolation, print out a message warning the user that this
C     has occurred.
      IERR = 13
      OPEN(IERR,FILE='EXTRAP.WARNINGS',FORM='FORMATTED')

      LAYTROP = 0
      DO 6000 LAY = 1, NLAYERS
C        Find the two reference pressures on either side of the
C        layer pressure.  Store them in JP and JP1.  Store in FP the
C        fraction of the difference (in ln(pressure)) between these
C        two values that the layer pressure lies.
         PLOG = ALOG(PAVEL(LAY))
         JP(LAY) = INT(36. - 5*(PLOG+0.04))
         IF (JP(LAY) .LT. 1) THEN
            JP(LAY) = 1
            A1 = 'PRESSURE'
            A2 = 'HIGH'
            WRITE(IERR,9810)A1,PAVEL(LAY)
            WRITE(IERR,9811)LAY,A2
            WRITE(IERR,9812)
         ELSEIF (JP(LAY) .GT. 58) THEN
            JP(LAY) = 58
            A1 = 'PRESSURE'
            A2 = 'LOW'
            WRITE(IERR,9810)A1,PAVEL(LAY)
            WRITE(IERR,9811)LAY,A2
            WRITE(IERR,9812)
         ENDIF
         JP1 = JP(LAY) + 1
         FP = 5. * (PREFLOG(JP(LAY)) - PLOG)
C        Determine, for each reference pressure (JP and JP1), which
C        reference temperature (these are different for each  
C        reference pressure) is nearest the layer temperature but does
C        not exceed it.  Store these indices in JT and JT1, resp.
C        Store in FT (resp. FT1) the fraction of the way between JT
C        (JT1) and the next highest reference temperature that the 
C        layer temperature falls.
         JT(LAY) = INT(3. + (TAVEL(LAY)-TREF(JP(LAY)))/15.)
         IF (JT(LAY) .LT. 1) THEN
            JT(LAY) = 1
            A1 = 'TEMPERATURE'
            A2 = 'LOW'
            WRITE(IERR,9810)A1,TAVEL(LAY)
            WRITE(IERR,9811)LAY,A2
            WRITE(IERR,9812)
         ELSEIF(JT(LAY) .GT. 4) THEN
            JT(LAY) = 4
         ENDIF
         FT = ((TAVEL(LAY)-TREF(JP(LAY)))/15.) - FLOAT(JT(LAY)-3)
         JT1(LAY) = INT(3. + (TAVEL(LAY)-TREF(JP1))/15.)
         IF (JT1(LAY) .LT. 1) THEN
            JT1(LAY) = 1
         ELSEIF(JT1(LAY) .GT. 4) THEN
            JT1(LAY) = 4
            A1 = 'TEMPERATURE'
            A2 = 'HIGH'
            WRITE(IERR,9810)A1,TAVEL(LAY)
            WRITE(IERR,9811)LAY,A2
            WRITE(IERR,9812)
         ENDIF
         FT1 = ((TAVEL(LAY)-TREF(JP1))/15.) - FLOAT(JT1(LAY)-3)

C        Calculate needed column amounts.
         COLH2O(LAY) = 1.E-20 * COLDRY(LAY) * WATER(LAY)
         COLCO2(LAY) = 1.E-20 * COLDRY(LAY) * CO2(LAY)
         COLO3(LAY) = 1.E-20 * COLDRY(LAY) * OZ(LAY)

C        If the pressure is less than ~100mb, perform a different
C        set of species interpolations.
         IF (PLOG .LE. 4.56) GO TO 5300
         LAYTROP =  LAYTROP + 1

C        Set up factors needed to separately include the water vapor
C        self-continuum in the calculation of absorption coefficient.
         SCALEFAC = PAVEL(LAY) * STPFAC / TAVEL(LAY)
         SELFFAC(LAY) = SCALEFAC * WATER(LAY)/(1.+WATER(LAY))
         FACTOR = (TAVEL(LAY)-188.0)/7.2
         INDSELF(LAY) = MIN(9, MAX(1, INT(FACTOR)-7))
         SELFFRAC(LAY) = FACTOR - FLOAT(INDSELF(LAY) + 7)

C        Determine, for reference pressures JP and JP1, which
C        reference species ratio is just smaller than the actual
C        species ratio.  Store these in JS and JS1, respectively.
C        If there is no binary species parameter for a band, then
C        JS is set to 1.  Also, store the respective fractions 
C        in F*S and F*S1.

         JS(LAY,1) = 1
         JS1(LAY,1) = 1

         JS(LAY,2) = 1
         JS1(LAY,2) = 1

         SPECPARM = WATER(LAY)/(WATER(LAY) + STRRAT31*CO2(LAY))
         DO 3310 IREF = 2, 6
            IF(SPECPARM .LE. SPARM31(IREF,JP(LAY))) GO TO 3320
 3310    CONTINUE
 3320    CONTINUE
         JS(LAY,3) = IREF - 1
	 F3S = (SPECPARM - SPARM31(IREF-1,JP(LAY))) /
     &        (SPARM31(IREF,JP(LAY)) - SPARM31(IREF-1,JP(LAY)))
         DO 3350 IREF = 2, 6
            IF(SPECPARM .LE. SPARM31(IREF,JP1)) GO TO 3360
 3350    CONTINUE
 3360    CONTINUE
         JS1(LAY,3) = IREF - 1
	 F3S1 = (SPECPARM - SPARM31(IREF-1,JP1)) /
     &        (SPARM31(IREF,JP1) - SPARM31(IREF-1,JP1))

         SPECPARM = WATER(LAY)/(WATER(LAY) + STRRAT41*CO2(LAY))
         DO 3410 IREF = 2, 6
            IF(SPECPARM .LE. SPARM41(IREF,JP(LAY))) GO TO 3420
 3410    CONTINUE
 3420    CONTINUE
         JS(LAY,4) = IREF - 1
	 F4S = (SPECPARM - SPARM41(IREF-1,JP(LAY))) /
     &        (SPARM41(IREF,JP(LAY)) - SPARM41(IREF-1,JP(LAY)))
         DO 3450 IREF = 2, 6
            IF(SPECPARM .LE. SPARM41(IREF,JP1)) GO TO 3460
 3450    CONTINUE
 3460    CONTINUE
         JS1(LAY,4) = IREF - 1
	 F4S1 = (SPECPARM - SPARM41(IREF-1,JP1)) /
     &        (SPARM41(IREF,JP1) - SPARM41(IREF-1,JP1))

         SPECPARM = WATER(LAY)/(WATER(LAY) + STRRAT51*CO2(LAY))
         DO 3510 IREF = 2, 6
            IF(SPECPARM .LE. SPARM51(IREF,JP(LAY))) GO TO 3520
 3510    CONTINUE
 3520    CONTINUE
         JS(LAY,5) = IREF - 1
	 F5S = (SPECPARM - SPARM51(IREF-1,JP(LAY))) /
     &        (SPARM51(IREF,JP(LAY)) - SPARM51(IREF-1,JP(LAY)))
         DO 3550 IREF = 2, 6
            IF(SPECPARM .LE. SPARM51(IREF,JP1)) GO TO 3560
 3550    CONTINUE
 3560    CONTINUE
         JS1(LAY,5) = IREF - 1
	 F5S1 = (SPECPARM - SPARM51(IREF-1,JP1)) /
     &        (SPARM51(IREF,JP1) - SPARM51(IREF-1,JP1))

         JS(LAY,6) = 1
         JS1(LAY,6) = 1

         GO TO 5400

C        Above LAYTROP.
 5300    CONTINUE

         JS(LAY,1) = 1
         JS1(LAY,1) = 1

         JS(LAY,2) = 1
         JS1(LAY,2) = 1

         SPECPARM = WATER(LAY)/(WATER(LAY) + STRRAT32*CO2(LAY))
         IREF = 3
         IF(SPECPARM .LE. SPARM32(2,JP(LAY))) IREF = 2
         JS(LAY,3) = IREF - 1
	 F3S = (SPECPARM - SPARM32(IREF-1,JP(LAY))) /
     &        (SPARM32(IREF,JP(LAY)) - SPARM32(IREF-1,JP(LAY)))
         IREF = 3
         IF(SPECPARM .LE. SPARM32(2,JP1)) IREF = 2
         JS1(LAY,3) = IREF - 1
	 F3S1 = (SPECPARM - SPARM32(IREF-1,JP1)) /
     &        (SPARM32(IREF,JP1) - SPARM32(IREF-1,JP1))

         SPECPARM = OZ(LAY)/(OZ(LAY) + STRRAT42*CO2(LAY))
         IREF = 3
         IF(SPECPARM .LE. SPARM42(2,JP(LAY))) IREF = 2
         JS(LAY,4) = IREF - 1
	 F4S = (SPECPARM - SPARM42(IREF-1,JP(LAY))) /
     &        (SPARM42(IREF,JP(LAY)) - SPARM42(IREF-1,JP(LAY)))
         IREF = 3
         IF(SPECPARM .LE. SPARM42(2,JP1)) IREF = 2
         JS1(LAY,4) = IREF - 1
	 F4S1 = (SPECPARM - SPARM42(IREF-1,JP1)) /
     &        (SPARM42(IREF,JP1) - SPARM42(IREF-1,JP1))

         SPECPARM = OZ(LAY)/(OZ(LAY) + STRRAT52*CO2(LAY))
         IREF = 3
         IF(SPECPARM .LE. SPARM52(2,JP(LAY))) IREF = 2
         JS(LAY,5) = IREF - 1
	 F5S = (SPECPARM - SPARM52(IREF-1,JP(LAY))) /
     &        (SPARM52(IREF,JP(LAY)) - SPARM52(IREF-1,JP(LAY)))
         IREF = 3
         IF(SPECPARM .LE. SPARM52(2,JP1)) IREF = 2
         JS1(LAY,5) = IREF - 1
	 F5S1 = (SPECPARM - SPARM52(IREF-1,JP1)) /
     &        (SPARM52(IREF,JP1) - SPARM52(IREF-1,JP1))

 5400    CONTINUE
C
C        We have now isolated the layer ln pressure, temperature,
C        and binary species parameter between two reference pressures, 
C        two reference temperatures (for each reference pressure), and
C        two reference binary species parameters (for each reference 
C        pressure).  We multiply the pressure fraction FP (and 1 - FP) 
C        with the appropriate temperature and species fractions to get 
C        the combined factors needed for the interpolation that yields
C        the optical depths (performed in routines TAUGBn for band n).
C        Depending on the band (due to differing pressure dependence
C        of the reference data), the ln(pressure) interpolation can 
C        either be linear in FP or go as the square root of FP. 

         RTFP = SQRT(FP)
         COMPFP = 1. - FP

         FAC010(LAY,1) = (1.-RTFP) * FT
         FAC000(LAY,1) = (1.-RTFP) * (1. - FT)
         FAC011(LAY,1) = RTFP * FT1
         FAC001(LAY,1) = RTFP * (1. - FT1)
         FAC010(LAY,2) = (1.-RTFP) * FT
         FAC000(LAY,2) = (1.-RTFP) * (1. - FT)
         FAC011(LAY,2) = RTFP * FT1
         FAC001(LAY,2) = RTFP * (1. - FT1)
         FAC = FT * F3S
         FAC110(LAY,3) = COMPFP * FAC
         FAC100(LAY,3) = COMPFP * (F3S - FAC)
         FAC010(LAY,3) = COMPFP * (FT - FAC)
         FAC000(LAY,3) = COMPFP * (1. - FT - F3S + FAC)
         FAC1 = FT1 * F3S1
         FAC111(LAY,3) = FP * FAC1
         FAC101(LAY,3) = FP * (F3S1 - FAC1)
         FAC011(LAY,3) = FP * (FT1 - FAC1)
         FAC001(LAY,3) = FP * (1. - FT1 - F3S1 + FAC1)
         FAC = FT * F4S
         FAC110(LAY,4) = COMPFP * FAC
         FAC100(LAY,4) = COMPFP * (F4S - FAC)
         FAC010(LAY,4) = COMPFP * (FT - FAC)
         FAC000(LAY,4) = COMPFP * (1. - FT - F4S + FAC)
         FAC1 = FT1 * F4S1
         FAC111(LAY,4) = FP * FAC1
         FAC101(LAY,4) = FP * (F4S1 - FAC1)
         FAC011(LAY,4) = FP * (FT1 - FAC1)
         FAC001(LAY,4) = FP * (1. - FT1 - F4S1 + FAC1)
         FAC = FT * F5S
         FAC110(LAY,5) = COMPFP * FAC
         FAC100(LAY,5) = COMPFP * (F5S - FAC)
         FAC010(LAY,5) = COMPFP * (FT - FAC)
         FAC000(LAY,5) = COMPFP * (1. - FT - F5S + FAC)
         FAC1 = FT1 * F5S1
         FAC111(LAY,5) = FP * FAC1
         FAC101(LAY,5) = FP * (F5S1 - FAC1)
         FAC011(LAY,5) = FP * (FT1 - FAC1)
         FAC001(LAY,5) = FP * (1. - FT1 - F5S1 + FAC1)
         FAC010(LAY,6) = COMPFP * FT
         FAC000(LAY,6) = COMPFP * (1. - FT)
         FAC011(LAY,6) = FP * FT1
         FAC001(LAY,6) = FP * (1. - FT1)
 6000 CONTINUE
      CLOSE(IERR)

C     Find out between what two degrees the level and layer 
C     temperatures fall.  This is needed for the interpolation
C     that calculates the Planck functions.
      TLEVFRAC(0) = TSFC - INT(TSFC)
      INDLEV(0) =  TSFC - 179.
      DO 7000 LEV = 1, NLAYERS
         TLEVFRAC(LEV) = TZ(LEV) - INT(TZ(LEV))
         TLAYFRAC(LEV) = TAVEL(LEV) - INT(TAVEL(LEV))
         INDLEV(LEV) = TZ(LEV) - 179.
         INDLAY(LEV) = TAVEL(LEV) - 179.
 7000 CONTINUE

 9810 FORMAT(' EXTRAPOLATION WARNING:  YOUR ',A11,' OF ',1P,G13.7,0P)
 9811 FORMAT(' FROM LAYER ',I3,' IS TOO ',A4,' FOR INTERPOLATION.')
 9812 FORMAT(' THE OPTICAL DEPTHS CALCULATED MAY NOT BE ACCURATE.')

      RETURN
      END





