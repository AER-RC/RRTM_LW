C     path:      %P%
C     revision:  $Revision$
C     created:   $Date$  
C     presently: %H%  %T%
      SUBROUTINE CLDPROP(ICLDATM)

C     Purpose:  Compute the cloud optical depth(s) for each cloudy
C               layer.

      PARAMETER (MXLAY=203)
      PARAMETER (NBANDS = 16)

      COMMON /CONTROL/   NUMANGS, IOUT, ISTART, IEND, ICLD
      COMMON /PROFILE/   NLAYERS,PAVEL(MXLAY),TAVEL(MXLAY),
     &                   PZ(0:MXLAY),TZ(0:MXLAY)
      COMMON /CLOUDIN/   INFLAG,CLDDAT1(MXLAY),CLDDAT2(MXLAY),
     &                   ICEFLAG,LIQFLAG,CLDDAT3(MXLAY),CLDDAT4(MXLAY)
      COMMON /CLOUDDAT/  NCBANDS,CLDFRAC(MXLAY),TAUCLOUD(MXLAY,NBANDS)
      COMMON /HVERSN/    HVRRTM,HVRREG,HVRRTR,HVRATM,HVRSET,HVRTAU,
     *                   HVRRGC,HVRRTC,HVRCLD,HVRDUM,HVRUTL,HVREXT

      CHARACTER*8 HVRRTM,HVRREG,HVRRTR,HVRATM,HVRSET,HVRTAU,
     *            HVRRGC,HVRRTC,HVRCLD,HVRDUM,HVRUTL,HVREXT
      DIMENSION ABSICE0(2), ABSICE1(2,5), ABSICE2(40,16)
      DIMENSION ABSLIQ1(58,16)
      DIMENSION ABSCOICE(NBANDS), ABSCOLIQ(NBANDS)
      DIMENSION IPAT(16,0:2)

      DATA EPS /1.E-6/

      DATA IPAT /1,1,1,1,1,1,1,1,1, 1, 1, 1, 1, 1, 1, 1,
     &           1,2,3,3,3,4,4,4,5, 5, 5, 5, 5, 5, 5, 5,
     &           1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16/

C     Explanation of the method for each value of INFLAG.  Values of
C     0 or 1 for INFLAG do not distingish being liquid and ice clouds.
C     INFLAG = 2 does distinguish between liquid and ice clouds, and
C     requires further user input to specify the method to be used to 
C     compute the aborption due to each.
C     INFLAG = 0:  For each cloudy layer, the cloud fraction and (gray)
C                  optical depth are input.  
C     INFLAG = 1:  For each cloudy layer, the cloud fraction and cloud
C                  water path (g/m2) are input.  The (gray) cloud optical 
C                  depth is computed as in CCM2.
C     INFLAG = 2:  For each cloudy layer, the cloud fraction, cloud 
C                  water path (g/m2), and cloud ice fraction are input.
C       ICEFLAG = 0:  The ice effective radius (microns) is input and the
C                     optical depths due to ice clouds are computed as in CCM3.
C       ICEFLAG = 1:  The ice effective radius (microns) is input and the
C                     optical depths due to ice clouds are computed as in 
C                     Ebert and Curry, JGR, 97, 3831-3836 (1992).  The 
C                     spectral regions in this work have been matched with
C                     the spectral bands in RRTM to as great an extent 
C                     as possible:  
C                     E&C 1      IB = 5      RRTM bands 9-16
C                     E&C 2      IB = 4      RRTM bands 6-8
C                     E&C 3      IB = 3      RRTM bands 3-5
C                     E&C 4      IB = 2      RRTM band 2
C                     E&C 5      IB = 1      RRTM band 1
C       ICEFLAG = 2:  The ice effective radius (microns) is input and the
C                     optical depths due to ice clouds are computed as in 
C                     Streamer (reference:  J. Key, Streamer User's Guide, 
C                     Technical Report 96-01, Department of Geography,
C                     Boston University, 85 pp. (1996)).  The values of 
C                     absorption coefficients appropriate for the spectral 
C                     bands of RRTM were obtained by an averaging procedure 
C                     based on the work of J. Pinto (private communication).
C       LIQFLAG = 0:  The optical depths due to water clouds are computed as
C                     in CCM3.
C       LIQFLAG = 1:  The water droplet effective radius (microns) is input 
C                     and the optical depths due to water clouds are computed 
C                     as in Hu and Stamnes, J., Clim., 6, 728-742, (1993).
C                     The values for absorption coefficients appropriate for
C                     the spectral bands in RRTM have been obtained for a 
C                     range of effective radii by an averaging procedure 
C                     based on the work of J. Pinto (private communication).
C                     Linear interpolation is used to get the absorption 
C                     coefficients for the input effective radius.

C     ABSCLDn is the liquid water absorption coefficient (m2/g). 
C     For INFLAG = 1.
      DATA ABSCLD1 /0.0602410/
C  
C     Everything below is for INFLAG = 2.

C     ABSICEn(J,IB) are the parameters needed to compute the liquid water 
C     absorption coefficient in spectral region IB for ICEFLAG=n.  The 
C     units of ABSICEn(1,IB) are m2/g and ABSICEn(2,IB) has units 
C     (microns (m2/g)).
C     For ICEFLAG = 0.
      DATA ABSICE0 /0.005,  1.0/
C     For ICEFLAG = 1.
      DATA ABSICE1 /0.0036, 1.136,
     &              0.0068, 0.600,
     &              0.0003, 1.338,
     &              0.0016, 1.166,
     &              0.0020, 1.118/

C     For ICEFLAG = 2.  In each band, the absorption
C     coefficients are listed for a range of effective radii from 13.0
C     to 130.0 microns in increments of 3.0 micronS.
      DATA (ABSICE2(I,1),I=1,40)/
     &5.308137E-02,4.610615E-02,4.128355E-02,3.767141E-02,3.481773E-02,
     &3.247569E-02,3.049760E-02,2.878894E-02,2.728603E-02,2.594420E-02,
     &2.473109E-02,2.362260E-02,2.260035E-02,2.165006E-02,2.076046E-02,
     &1.992250E-02,1.912880E-02,1.837331E-02,1.765101E-02,1.695766E-02,
     &1.628969E-02,1.564405E-02,1.501811E-02,1.440959E-02,1.381653E-02,
     &1.323720E-02,1.267006E-02,1.211378E-02,1.156715E-02,1.102911E-02,
     &1.049869E-02,9.975041E-03,9.457376E-03,8.944992E-03,8.437249E-03,
     &7.933561E-03,7.433394E-03,6.936257E-03,6.441701E-03,5.949308E-03/
      DATA (ABSICE2(I,2),I=1,40)/
     &3.808872E-02,3.502664E-02,3.231765E-02,2.998940E-02,2.798500E-02,
     &2.624287E-02,2.471175E-02,2.335146E-02,2.213092E-02,2.102607E-02,
     &2.001804E-02,1.909193E-02,1.823581E-02,1.744003E-02,1.669669E-02,
     &1.599927E-02,1.534234E-02,1.472132E-02,1.413234E-02,1.357209E-02,
     &1.303774E-02,1.252681E-02,1.203718E-02,1.156698E-02,1.111455E-02,
     &1.067846E-02,1.025741E-02,9.850260E-03,9.455976E-03,9.073645E-03,
     &8.702438E-03,8.341606E-03,7.990473E-03,7.648424E-03,7.314896E-03,
     &6.989379E-03,6.671402E-03,6.360537E-03,6.056387E-03,5.758587E-03/
      DATA (ABSICE2(I,3),I=1,40)/
     &7.402665E-02,6.209102E-02,5.216790E-02,4.453142E-02,3.864415E-02,
     &3.402727E-02,3.033647E-02,2.733239E-02,2.484767E-02,2.276330E-02,
     &2.099304E-02,1.947312E-02,1.815559E-02,1.700375E-02,1.598911E-02,
     &1.508925E-02,1.428628E-02,1.356580E-02,1.291609E-02,1.232749E-02,
     &1.179201E-02,1.130298E-02,1.085478E-02,1.044263E-02,1.006249E-02,
     &9.710860E-03,9.384746E-03,9.081546E-03,8.798995E-03,8.535113E-03,
     &8.288163E-03,8.056611E-03,7.839102E-03,7.634433E-03,7.441532E-03,
     &7.259442E-03,7.087304E-03,6.924348E-03,6.769881E-03,6.623274E-03/
      DATA (ABSICE2(I,4),I=1,40)/
     &8.094820E-02,6.403291E-02,5.250755E-02,4.428943E-02,3.818391E-02,
     &3.349259E-02,2.978790E-02,2.679591E-02,2.433399E-02,2.227614E-02,
     &2.053280E-02,1.903874E-02,1.774534E-02,1.661571E-02,1.562134E-02,
     &1.473990E-02,1.395366E-02,1.324836E-02,1.261241E-02,1.203630E-02,
     &1.151218E-02,1.103347E-02,1.059465E-02,1.019105E-02,9.818699E-03,
     &9.474168E-03,9.154525E-03,8.857223E-03,8.580047E-03,8.321061E-03,
     &8.078565E-03,7.851061E-03,7.637226E-03,7.435882E-03,7.245985E-03,
     &7.066598E-03,6.896884E-03,6.736091E-03,6.583542E-03,6.438624E-03/
      DATA (ABSICE2(I,5),I=1,40)/
     &8.107771E-02,6.356966E-02,5.202664E-02,4.388361E-02,3.785256E-02,
     &3.321862E-02,2.955444E-02,2.658946E-02,2.414439E-02,2.209589E-02,
     &2.035644E-02,1.886226E-02,1.756586E-02,1.643112E-02,1.543016E-02,
     &1.454108E-02,1.374647E-02,1.303233E-02,1.238726E-02,1.180189E-02,
     &1.126846E-02,1.078049E-02,1.033252E-02,9.919924E-03,9.538748E-03,
     &9.185608E-03,8.857581E-03,8.552134E-03,8.267059E-03,8.000423E-03,
     &7.750529E-03,7.515880E-03,7.295149E-03,7.087158E-03,6.890857E-03,
     &6.705306E-03,6.529664E-03,6.363173E-03,6.205150E-03,6.054977E-03/
      DATA (ABSICE2(I,6),I=1,40)/
     &6.617399E-02,5.382356E-02,4.548024E-02,3.943650E-02,3.483962E-02,
     &3.121482E-02,2.827613E-02,2.584066E-02,2.378586E-02,2.202637E-02,
     &2.050087E-02,1.916411E-02,1.798196E-02,1.692814E-02,1.598211E-02,
     &1.512754E-02,1.435129E-02,1.364266E-02,1.299284E-02,1.239453E-02,
     &1.184158E-02,1.132883E-02,1.085186E-02,1.040690E-02,9.990698E-03,
     &9.600436E-03,9.233660E-03,8.888220E-03,8.562226E-03,8.254011E-03,
     &7.962094E-03,7.685159E-03,7.422029E-03,7.171651E-03,6.933078E-03,
     &6.705456E-03,6.488015E-03,6.280053E-03,6.080937E-03,5.890087E-03/
      DATA (ABSICE2(I,7),I=1,40)/
     &5.525622E-02,4.882834E-02,4.307237E-02,3.828016E-02,3.432924E-02,
     &3.105137E-02,2.830273E-02,2.597146E-02,2.397252E-02,2.224126E-02,
     &2.072819E-02,1.939498E-02,1.821160E-02,1.715425E-02,1.620385E-02,
     &1.534496E-02,1.456495E-02,1.385341E-02,1.320165E-02,1.260242E-02,
     &1.204956E-02,1.153786E-02,1.106285E-02,1.062069E-02,1.020806E-02,
     &9.822069E-03,9.460187E-03,9.120204E-03,8.800168E-03,8.498356E-03,
     &8.213235E-03,7.943442E-03,7.687759E-03,7.445092E-03,7.214460E-03,
     &6.994979E-03,6.785846E-03,6.586339E-03,6.395797E-03,6.213622E-03/
      DATA (ABSICE2(I,8),I=1,40)/
     &6.204709E-02,5.447577E-02,4.693652E-02,4.068410E-02,3.566129E-02,
     &3.162111E-02,2.833719E-02,2.563384E-02,2.338009E-02,2.147883E-02,
     &1.985755E-02,1.846151E-02,1.724884E-02,1.618712E-02,1.525092E-02,
     &1.442006E-02,1.367838E-02,1.301277E-02,1.241252E-02,1.186877E-02,
     &1.137421E-02,1.092266E-02,1.050894E-02,1.012866E-02,9.778052E-03,
     &9.453889E-03,9.153392E-03,8.874146E-03,8.614049E-03,8.371264E-03,
     &8.144174E-03,7.931357E-03,7.731552E-03,7.543641E-03,7.366629E-03,
     &7.199625E-03,7.041832E-03,6.892534E-03,6.751085E-03,6.616904E-03/
      DATA (ABSICE2(I,9),I=1,40)/
     &6.886530E-02,5.796039E-02,4.871664E-02,4.157422E-02,3.606747E-02,
     &3.175486E-02,2.831396E-02,2.551931E-02,2.321294E-02,2.128248E-02,
     &1.964647E-02,1.824476E-02,1.703216E-02,1.597411E-02,1.504384E-02,
     &1.422028E-02,1.348667E-02,1.282952E-02,1.223786E-02,1.170268E-02,
     &1.121653E-02,1.077318E-02,1.036740E-02,9.994772E-03,9.651518E-03,
     &9.334412E-03,9.040674E-03,8.767897E-03,8.513989E-03,8.277122E-03,
     &8.055693E-03,7.848290E-03,7.653665E-03,7.470712E-03,7.298447E-03,
     &7.135991E-03,6.982557E-03,6.837440E-03,6.700003E-03,6.569673E-03/
      DATA (ABSICE2(I,10),I=1,40)/
     &7.124487E-02,5.846725E-02,4.869304E-02,4.140966E-02,3.588534E-02,
     &3.159352E-02,2.818254E-02,2.541683E-02,2.313524E-02,2.122481E-02,
     &1.960443E-02,1.821456E-02,1.701063E-02,1.595869E-02,1.503244E-02,
     &1.421124E-02,1.347864E-02,1.282143E-02,1.222886E-02,1.169209E-02,
     &1.120381E-02,1.075791E-02,1.034925E-02,9.973491E-03,9.626911E-03,
     &9.306335E-03,9.009024E-03,8.732604E-03,8.475009E-03,8.234434E-03,
     &8.009292E-03,7.798186E-03,7.599882E-03,7.413280E-03,7.237405E-03,
     &7.071384E-03,6.914434E-03,6.765853E-03,6.625009E-03,6.491330E-03/
      DATA (ABSICE2(I,11),I=1,40)/
     &6.946157E-02,5.694109E-02,4.754750E-02,4.056198E-02,3.525241E-02,
     &3.111421E-02,2.781424E-02,2.512985E-02,2.290863E-02,2.104354E-02,
     &1.945751E-02,1.809386E-02,1.691003E-02,1.587352E-02,1.495910E-02,
     &1.414692E-02,1.342116E-02,1.276904E-02,1.218016E-02,1.164598E-02,
     &1.115939E-02,1.071445E-02,1.030618E-02,9.930324E-03,9.583261E-03,
     &9.261888E-03,8.963524E-03,8.685843E-03,8.426821E-03,8.184683E-03,
     &7.957872E-03,7.745014E-03,7.544891E-03,7.356422E-03,7.178643E-03,
     &7.010693E-03,6.851799E-03,6.701265E-03,6.558467E-03,6.422836E-03/
      DATA (ABSICE2(I,12),I=1,40)/
     &4.162544E-02,3.774206E-02,3.387559E-02,3.048216E-02,2.759398E-02,
     &2.514770E-02,2.306767E-02,2.128679E-02,1.975015E-02,1.841389E-02,
     &1.724326E-02,1.621062E-02,1.529388E-02,1.447525E-02,1.374029E-02,
     &1.307717E-02,1.247616E-02,1.192916E-02,1.142938E-02,1.097111E-02,
     &1.054952E-02,1.016048E-02,9.800441E-03,9.466351E-03,9.155566E-03,
     &8.865782E-03,8.594987E-03,8.341414E-03,8.103506E-03,7.879884E-03,
     &7.669326E-03,7.470744E-03,7.283166E-03,7.105719E-03,6.937623E-03,
     &6.778171E-03,6.626727E-03,6.482715E-03,6.345612E-03,6.214942E-03/
      DATA (ABSICE2(I,13),I=1,40)/
     &5.989806E-02,5.056817E-02,4.311112E-02,3.735098E-02,3.285293E-02,
     &2.927343E-02,2.637042E-02,2.397535E-02,2.196938E-02,2.026708E-02,
     &1.880583E-02,1.753882E-02,1.643044E-02,1.545319E-02,1.458549E-02,
     &1.381020E-02,1.311354E-02,1.248431E-02,1.191334E-02,1.139300E-02,
     &1.091697E-02,1.047989E-02,1.007724E-02,9.705182E-03,9.360395E-03,
     &9.040037E-03,8.741639E-03,8.463052E-03,8.202396E-03,7.958020E-03,
     &7.728467E-03,7.512449E-03,7.308820E-03,7.116561E-03,6.934759E-03,
     &6.762596E-03,6.599337E-03,6.444320E-03,6.296944E-03,6.156667E-03/
      DATA (ABSICE2(I,14),I=1,40)/
     &5.427438E-02,4.682405E-02,4.056412E-02,3.555101E-02,3.153408E-02,
     &2.827650E-02,2.559668E-02,2.336126E-02,2.147261E-02,1.985856E-02,
     &1.846503E-02,1.725091E-02,1.618448E-02,1.524093E-02,1.440064E-02,
     &1.364787E-02,1.296990E-02,1.235630E-02,1.179850E-02,1.128934E-02,
     &1.082285E-02,1.039396E-02,9.998390E-03,9.632453E-03,9.292999E-03,
     &8.977300E-03,8.682986E-03,8.407991E-03,8.150503E-03,7.908927E-03,
     &7.681856E-03,7.468042E-03,7.266375E-03,7.075864E-03,6.895622E-03,
     &6.724854E-03,6.562843E-03,6.408943E-03,6.262568E-03,6.123188E-03/
      DATA (ABSICE2(I,15),I=1,40)/
     &3.734649E-02,3.412194E-02,3.079772E-02,2.783212E-02,2.528118E-02,
     &2.310335E-02,2.123972E-02,1.963551E-02,1.824481E-02,1.703044E-02,
     &1.596258E-02,1.501737E-02,1.417560E-02,1.342171E-02,1.274301E-02,
     &1.212909E-02,1.157132E-02,1.106249E-02,1.059658E-02,1.016849E-02,
     &9.773873E-03,9.409024E-03,9.070757E-03,8.756320E-03,8.463322E-03,
     &8.189678E-03,7.933561E-03,7.693367E-03,7.467678E-03,7.255238E-03,
     &7.054931E-03,6.865762E-03,6.686840E-03,6.517366E-03,6.356622E-03,
     &6.203960E-03,6.058793E-03,5.920590E-03,5.788869E-03,5.663191E-03/
      DATA (ABSICE2(I,16),I=1,40)/
     &4.226924E-02,3.719702E-02,3.276247E-02,2.909934E-02,2.609178E-02,
     &2.360471E-02,2.152589E-02,1.976868E-02,1.826735E-02,1.697197E-02,
     &1.584428E-02,1.485463E-02,1.397979E-02,1.320136E-02,1.250459E-02,
     &1.187752E-02,1.131042E-02,1.079523E-02,1.032526E-02,9.894930E-03,
     &9.499496E-03,9.134956E-03,8.797885E-03,8.485341E-03,8.194783E-03,
     &7.924007E-03,7.671090E-03,7.434349E-03,7.212302E-03,7.003645E-03,
     &6.807219E-03,6.621995E-03,6.447055E-03,6.281580E-03,6.124830E-03,
     &5.976146E-03,5.834927E-03,5.700635E-03,5.572778E-03,5.450912E-03/

C     For LIQFLAG = 0.
      DATA ABSLIQ0 /0.0903614/

C     For LIQFLAG = 1.  In each band, the absorption
C     coefficients are listed for a range of effective radii from 2.5
C     to 59.5 microns in increments of 1.0 micron.
      DATA (ABSLIQ1(I,1),I=1,58)/
C     Band 1.
     &2.317311E-02,6.109513E-02,6.850127E-02,7.025095E-02,7.045720E-02,
     &7.014763E-02,6.964456E-02,6.906029E-02,6.843251E-02,6.776883E-02,
     &6.740428E-02,6.471351E-02,6.215508E-02,5.977013E-02,5.755809E-02,
     &5.550328E-02,5.358573E-02,5.178553E-02,5.008436E-02,4.846591E-02,
     &4.691585E-02,4.542158E-02,4.397195E-02,4.255702E-02,4.116779E-02,
     &3.979603E-02,3.843408E-02,3.707472E-02,3.593146E-02,3.477353E-02,
     &3.368067E-02,3.264653E-02,3.166568E-02,3.073339E-02,2.984557E-02,
     &2.899858E-02,2.818924E-02,2.741469E-02,2.667239E-02,2.596005E-02,
     &2.527561E-02,2.461722E-02,2.398318E-02,2.337194E-02,2.278211E-02,
     &2.221240E-02,2.166162E-02,2.112870E-02,2.061262E-02,2.011246E-02,
     &1.962737E-02,1.915655E-02,1.869927E-02,1.825485E-02,1.782263E-02,
     &1.740204E-02,1.699251E-02,1.659352E-02/
      DATA (ABSLIQ1(I,2),I=1,58)/
C     Band 2.
     &9.103496E-02,0.1480995,0.1402309,0.1300658,0.1214233,
     &0.1140277,0.1073661,1.010236E-01,9.466848E-02,8.801306E-02,
     &8.222715E-02,7.610165E-02,7.087689E-02,6.633009E-02,6.231412E-02,
     &5.872623E-02,5.549144E-02,5.255308E-02,4.986711E-02,4.739856E-02,
     &4.511911E-02,4.300551E-02,4.103841E-02,3.920156E-02,3.748114E-02,
     &3.586532E-02,3.434392E-02,3.290809E-02,3.176277E-02,3.060455E-02,
     &2.951870E-02,2.849868E-02,2.753870E-02,2.663363E-02,2.577891E-02,
     &2.497048E-02,2.420469E-02,2.347828E-02,2.278828E-02,2.213206E-02,
     &2.150719E-02,2.091151E-02,2.034301E-02,1.979989E-02,1.928050E-02,
     &1.878332E-02,1.830698E-02,1.785018E-02,1.741177E-02,1.699066E-02,
     &1.658586E-02,1.619643E-02,1.582154E-02,1.546038E-02,1.511222E-02,
     &1.477637E-02,1.445221E-02,1.413913E-02/
      DATA (ABSLIQ1(I,3),I=1,58)/
C     Band 3.
     &0.2951735,0.2347648,0.1980376,0.1721142,0.1520826,
     &0.1356541,0.1216134,0.1092517,9.812634E-02,8.794481E-02,
     &8.125663E-02,7.445630E-02,6.863743E-02,6.360421E-02,5.920939E-02,
     &5.534017E-02,5.190868E-02,4.884551E-02,4.609508E-02,4.361243E-02,
     &4.136074E-02,3.930963E-02,3.743377E-02,3.571191E-02,3.412609E-02,
     &3.266100E-02,3.130355E-02,3.004245E-02,2.884966E-02,2.780765E-02,
     &2.683173E-02,2.591580E-02,2.505449E-02,2.424303E-02,2.347722E-02,
     &2.275331E-02,2.206794E-02,2.141810E-02,2.080110E-02,2.021450E-02,
     &1.965609E-02,1.912390E-02,1.861611E-02,1.813107E-02,1.766729E-02,
     &1.722339E-02,1.679811E-02,1.639031E-02,1.599892E-02,1.562296E-02,
     &1.526155E-02,1.491384E-02,1.457907E-02,1.425652E-02,1.394554E-02,
     &1.364551E-02,1.335586E-02,1.307606E-02/
      DATA (ABSLIQ1(I,4),I=1,58)/
C     Band 4.
     &0.3009248,0.2369488,0.1969469,0.1686924,0.1471901,
     &0.1299859,0.1157190,1.035677E-01,9.300276E-02,8.366584E-02,
     &7.710747E-02,7.070016E-02,6.522841E-02,6.050237E-02,5.638011E-02,
     &5.275342E-02,4.953838E-02,4.666896E-02,4.409246E-02,4.176635E-02,
     &3.965594E-02,3.773263E-02,3.597266E-02,3.435610E-02,3.286616E-02,
     &3.148853E-02,3.021100E-02,2.902305E-02,2.789483E-02,2.691087E-02,
     &2.598836E-02,2.512166E-02,2.430581E-02,2.353642E-02,2.280958E-02,
     &2.212182E-02,2.147004E-02,2.085146E-02,2.026356E-02,1.970411E-02,
     &1.917105E-02,1.866254E-02,1.817690E-02,1.771260E-02,1.726825E-02,
     &1.684257E-02,1.643440E-02,1.604265E-02,1.566636E-02,1.530459E-02,
     &1.495653E-02,1.462140E-02,1.429847E-02,1.398708E-02,1.368662E-02,
     &1.339651E-02,1.311622E-02,1.284525E-02/
      DATA (ABSLIQ1(I,5),I=1,58)/
C     Band 5.
     &0.2646912,0.2120182,0.1780086,0.1535393,0.1347209,
     &0.1195802,0.1069960,9.627716E-02,8.697099E-02,7.876703E-02,
     &7.292724E-02,6.709201E-02,6.209771E-02,5.777322E-02,5.399095E-02,
     &5.065378E-02,4.768655E-02,4.503014E-02,4.263744E-02,4.047042E-02,
     &3.849806E-02,3.669482E-02,3.503944E-02,3.351411E-02,3.210378E-02,
     &3.079567E-02,2.957880E-02,2.844375E-02,2.737899E-02,2.643902E-02,
     &2.555652E-02,2.472627E-02,2.394368E-02,2.320467E-02,2.250563E-02,
     &2.184333E-02,2.121488E-02,2.061771E-02,2.004948E-02,1.950810E-02,
     &1.899165E-02,1.849842E-02,1.802685E-02,1.757549E-02,1.714305E-02,
     &1.672834E-02,1.633026E-02,1.594781E-02,1.558006E-02,1.522615E-02,
     &1.488531E-02,1.455681E-02,1.423996E-02,1.393415E-02,1.363879E-02,
     &1.335334E-02,1.307730E-02,1.281020E-02/
      DATA (ABSLIQ1(I,6),I=1,58)/
C     Band 6.
     &8.811824E-02,0.1067453,9.797529E-02,8.996253E-02,8.352002E-02,
     &7.818990E-02,7.359390E-02,6.946955E-02,6.562659E-02,6.191484E-02,
     &5.833548E-02,5.493056E-02,5.196417E-02,4.933254E-02,4.696588E-02,
     &4.481484E-02,4.284310E-02,4.102307E-02,3.933320E-02,3.775632E-02,
     &3.627847E-02,3.488817E-02,3.357580E-02,3.233326E-02,3.115363E-02,
     &3.003095E-02,2.896007E-02,2.793646E-02,2.705024E-02,2.626182E-02,
     &2.550246E-02,2.477275E-02,2.407256E-02,2.340134E-02,2.275825E-02,
     &2.214224E-02,2.155220E-02,2.098694E-02,2.044526E-02,1.992599E-02,
     &1.942798E-02,1.895011E-02,1.849132E-02,1.805060E-02,1.762699E-02,
     &1.721959E-02,1.682755E-02,1.645004E-02,1.608633E-02,1.573569E-02,
     &1.539746E-02,1.507101E-02,1.475575E-02,1.445112E-02,1.415662E-02,
     &1.387173E-02,1.359601E-02,1.332903E-02/
      DATA (ABSLIQ1(I,7),I=1,58)/
C     Band 7.
     &4.321742E-02,7.360776E-02,6.983400E-02,6.652314E-02,6.419483E-02,
     &6.235512E-02,6.066381E-02,5.886798E-02,5.671241E-02,5.386294E-02,
     &4.995791E-02,4.862892E-02,4.701198E-02,4.528539E-02,4.354655E-02,
     &4.184799E-02,4.021685E-02,3.866581E-02,3.719923E-02,3.581678E-02,
     &3.451553E-02,3.329122E-02,3.213899E-02,3.105381E-02,3.003074E-02,
     &2.906507E-02,2.815236E-02,2.728853E-02,2.628212E-02,2.557439E-02,
     &2.487989E-02,2.420288E-02,2.354602E-02,2.291084E-02,2.229808E-02,
     &2.170793E-02,2.114022E-02,2.059448E-02,2.007010E-02,1.956633E-02,
     &1.908238E-02,1.861742E-02,1.817060E-02,1.774108E-02,1.732806E-02,
     &1.693073E-02,1.654832E-02,1.618011E-02,1.582540E-02,1.548352E-02,
     &1.515384E-02,1.483577E-02,1.452875E-02,1.423223E-02,1.394571E-02,
     &1.366873E-02,1.340083E-02,1.314159E-02/
      DATA (ABSLIQ1(I,8),I=1,58)/
C     Band 8.
     &0.1418807,7.154191E-02,6.303346E-02,6.111324E-02,6.019305E-02,
     &5.924194E-02,5.789680E-02,5.588757E-02,5.289232E-02,4.844620E-02,
     &4.608385E-02,4.560133E-02,4.454102E-02,4.318657E-02,4.170261E-02,
     &4.018504E-02,3.868917E-02,3.724607E-02,3.587214E-02,3.457487E-02,
     &3.335635E-02,3.221551E-02,3.114941E-02,3.015414E-02,2.922530E-02,
     &2.835836E-02,2.754884E-02,2.679246E-02,2.576920E-02,2.507042E-02,
     &2.439184E-02,2.373495E-02,2.310052E-02,2.248877E-02,2.189956E-02,
     &2.133248E-02,2.078695E-02,2.026227E-02,1.975769E-02,1.927241E-02,
     &1.880560E-02,1.835644E-02,1.792414E-02,1.750791E-02,1.710699E-02,
     &1.672066E-02,1.634822E-02,1.598900E-02,1.564237E-02,1.530773E-02,
     &1.498452E-02,1.467219E-02,1.437023E-02,1.407816E-02,1.379552E-02,
     &1.352189E-02,1.325685E-02,1.300002E-02/
      DATA (ABSLIQ1(I,9),I=1,58)/
C     Band 9.
     &6.727262E-02,6.610126E-02,6.478659E-02,6.337803E-02,6.189851E-02,
     &6.033352E-02,5.861357E-02,5.658760E-02,5.398387E-02,5.035361E-02,
     &4.716084E-02,4.636301E-02,4.503130E-02,4.345260E-02,4.178762E-02,
     &4.012609E-02,3.851708E-02,3.698595E-02,3.554424E-02,3.419543E-02,
     &3.293837E-02,3.176932E-02,3.068324E-02,2.967450E-02,2.873737E-02,
     &2.786621E-02,2.705567E-02,2.630077E-02,2.524495E-02,2.454240E-02,
     &2.386557E-02,2.321441E-02,2.258852E-02,2.198727E-02,2.140986E-02,
     &2.085540E-02,2.032297E-02,1.981160E-02,1.932034E-02,1.884823E-02,
     &1.839436E-02,1.795783E-02,1.753780E-02,1.713345E-02,1.674400E-02,
     &1.636871E-02,1.600688E-02,1.565786E-02,1.532101E-02,1.499576E-02,
     &1.468153E-02,1.437781E-02,1.408409E-02,1.379992E-02,1.352485E-02,
     &1.325846E-02,1.300036E-02,1.275018E-02/
      DATA (ABSLIQ1(I,10),I=1,58)/
C     Band 10.
     &7.970399E-02,7.638436E-02,7.364989E-02,7.135247E-02,6.930427E-02,
     &6.728066E-02,6.502274E-02,6.223945E-02,5.860926E-02,5.378152E-02,
     &5.146821E-02,4.972139E-02,4.773920E-02,4.569609E-02,4.368575E-02,
     &4.175686E-02,3.993275E-02,3.822235E-02,3.662650E-02,3.514157E-02,
     &3.376166E-02,3.247981E-02,3.128873E-02,3.018124E-02,2.915047E-02,
     &2.819001E-02,2.729394E-02,2.645681E-02,2.541650E-02,2.468318E-02,
     &2.397828E-02,2.330174E-02,2.265305E-02,2.203141E-02,2.143587E-02,
     &2.086533E-02,2.031867E-02,1.979474E-02,1.929243E-02,1.881062E-02,
     &1.834825E-02,1.790431E-02,1.747781E-02,1.706784E-02,1.667353E-02,
     &1.629405E-02,1.592863E-02,1.557655E-02,1.523712E-02,1.490969E-02,
     &1.459367E-02,1.428848E-02,1.399359E-02,1.370850E-02,1.343274E-02,
     &1.316586E-02,1.290746E-02,1.265714E-02/
      DATA (ABSLIQ1(I,11),I=1,58)/
C     Band 11.
     &0.1494379,0.1335347,0.1215419,0.1117433,1.032629E-01,
     &9.557744E-02,8.833824E-02,8.129427E-02,7.425327E-02,6.706086E-02,
     &6.387608E-02,5.977878E-02,5.598406E-02,5.253178E-02,4.941317E-02,
     &4.660136E-02,4.406436E-02,4.177059E-02,3.969101E-02,3.779976E-02,
     &3.607421E-02,3.449469E-02,3.304416E-02,3.170785E-02,3.047295E-02,
     &2.932829E-02,2.826416E-02,2.727201E-02,2.617887E-02,2.532765E-02,
     &2.452371E-02,2.376351E-02,2.304378E-02,2.236152E-02,2.171398E-02,
     &2.109866E-02,2.051327E-02,1.995571E-02,1.942406E-02,1.891659E-02,
     &1.843168E-02,1.796787E-02,1.752382E-02,1.709828E-02,1.669013E-02,
     &1.629830E-02,1.592185E-02,1.555986E-02,1.521153E-02,1.487608E-02,
     &1.455281E-02,1.424106E-02,1.394021E-02,1.364970E-02,1.336900E-02,
     &1.309761E-02,1.283507E-02,1.258095E-02/
      DATA (ABSLIQ1(I,12),I=1,58)/
C     Band 12.
     &3.719852E-02,3.885858E-02,3.990704E-02,4.043514E-02,4.046101E-02,
     &3.998338E-02,3.899533E-02,3.748863E-02,3.545507E-02,3.288698E-02,
     &3.325759E-02,3.224436E-02,3.123840E-02,3.025843E-02,2.931462E-02,
     &2.841198E-02,2.755245E-02,2.673607E-02,2.596181E-02,2.522801E-02,
     &2.453266E-02,2.387362E-02,2.324871E-02,2.265581E-02,2.209286E-02,
     &2.155791E-02,2.104914E-02,2.056482E-02,1.997489E-02,1.957035E-02,
     &1.917308E-02,1.878385E-02,1.840320E-02,1.803148E-02,1.766887E-02,
     &1.731546E-02,1.697124E-02,1.663615E-02,1.631005E-02,1.599278E-02,
     &1.568416E-02,1.538398E-02,1.509202E-02,1.480804E-02,1.453181E-02,
     &1.426309E-02,1.400164E-02,1.374723E-02,1.349962E-02,1.325858E-02,
     &1.302390E-02,1.279536E-02,1.257275E-02,1.235586E-02,1.214450E-02,
     &1.193848E-02,1.173761E-02,1.154172E-02/
      DATA (ABSLIQ1(I,13),I=1,58)/
C     Band 13.
     &3.118678E-02,4.483573E-02,4.902244E-02,4.964059E-02,4.868059E-02,
     &4.696101E-02,4.486300E-02,4.257948E-02,4.021383E-02,3.782356E-02,
     &3.742657E-02,3.603843E-02,3.470737E-02,3.344343E-02,3.224992E-02,
     &3.112640E-02,3.007038E-02,2.907835E-02,2.814632E-02,2.727021E-02,
     &2.644598E-02,2.566981E-02,2.493808E-02,2.424745E-02,2.359483E-02,
     &2.297739E-02,2.239252E-02,2.183785E-02,2.117925E-02,2.070756E-02,
     &2.024699E-02,1.979812E-02,1.936131E-02,1.893669E-02,1.852428E-02,
     &1.812396E-02,1.773555E-02,1.735881E-02,1.699345E-02,1.663916E-02,
     &1.629560E-02,1.596244E-02,1.563933E-02,1.532591E-02,1.502185E-02,
     &1.472680E-02,1.444043E-02,1.416242E-02,1.389246E-02,1.363023E-02,
     &1.337545E-02,1.312783E-02,1.288711E-02,1.265302E-02,1.242532E-02,
     &1.220376E-02,1.198811E-02,1.177816E-02/
      DATA (ABSLIQ1(I,14),I=1,58)/
C     Band 14.
     &1.589879E-02,3.506523E-02,4.008515E-02,4.072696E-02,3.981012E-02,
     &3.833056E-02,3.668287E-02,3.503265E-02,3.344965E-02,3.196089E-02,
     &3.137121E-02,3.033481E-02,2.934152E-02,2.839725E-02,2.750374E-02,
     &2.666039E-02,2.586536E-02,2.511614E-02,2.440996E-02,2.374397E-02,
     &2.311538E-02,2.252152E-02,2.195989E-02,2.142815E-02,2.092415E-02,
     &2.044588E-02,1.999154E-02,1.955943E-02,1.902540E-02,1.865979E-02,
     &1.829955E-02,1.794554E-02,1.759834E-02,1.725841E-02,1.692600E-02,
     &1.660128E-02,1.628433E-02,1.597515E-02,1.567370E-02,1.537987E-02,
     &1.509356E-02,1.481461E-02,1.454286E-02,1.427815E-02,1.402027E-02,
     &1.376906E-02,1.352431E-02,1.328583E-02,1.305344E-02,1.282695E-02,
     &1.260616E-02,1.239091E-02,1.218100E-02,1.197627E-02,1.177655E-02,
     &1.158167E-02,1.139148E-02,1.120583E-02/
      DATA (ABSLIQ1(I,15),I=1,58)/
C     Band 15.
     &5.020792E-03,2.176149E-02,2.554494E-02,2.594839E-02,2.536500E-02,
     &2.452808E-02,2.368427E-02,2.291585E-02,2.224506E-02,2.167161E-02,
     &2.114512E-02,2.058169E-02,2.004540E-02,1.953718E-02,1.905667E-02,
     &1.860280E-02,1.817419E-02,1.776930E-02,1.738656E-02,1.702444E-02,
     &1.668148E-02,1.635633E-02,1.604770E-02,1.575444E-02,1.547545E-02,
     &1.520974E-02,1.495642E-02,1.471463E-02,1.436841E-02,1.417275E-02,
     &1.397622E-02,1.377965E-02,1.358376E-02,1.338908E-02,1.319607E-02,
     &1.300508E-02,1.281637E-02,1.263017E-02,1.244664E-02,1.226591E-02,
     &1.208805E-02,1.191314E-02,1.174120E-02,1.157227E-02,1.140633E-02,
     &1.124338E-02,1.108341E-02,1.092637E-02,1.077224E-02,1.062098E-02,
     &1.047253E-02,1.032685E-02,1.018389E-02,1.004360E-02,9.905925E-03,
     &9.770801E-03,9.638176E-03,9.507998E-03/
      DATA (ABSLIQ1(I,16),I=1,58)/
C     Band 16.
     &1.398301E-02,3.704987E-02,3.928249E-02,3.759325E-02,3.534122E-02,
     &3.328979E-02,3.158282E-02,3.020512E-02,2.910579E-02,2.823262E-02,
     &2.653908E-02,2.559978E-02,2.471139E-02,2.387688E-02,2.309574E-02,
     &2.236567E-02,2.168353E-02,2.104590E-02,2.044933E-02,1.989053E-02,
     &1.936639E-02,1.887404E-02,1.841087E-02,1.797450E-02,1.756277E-02,
     &1.717372E-02,1.680557E-02,1.645673E-02,1.600733E-02,1.571652E-02,
     &1.543078E-02,1.515061E-02,1.487636E-02,1.460824E-02,1.434636E-02,
     &1.409077E-02,1.384144E-02,1.359831E-02,1.336129E-02,1.313026E-02,
     &1.290509E-02,1.268563E-02,1.247173E-02,1.226323E-02,1.205997E-02,
     &1.186180E-02,1.166854E-02,1.148005E-02,1.129616E-02,1.111674E-02,
     &1.094163E-02,1.077068E-02,1.060377E-02,1.044075E-02,1.028151E-02,
     &1.012590E-02,9.973826E-03,9.825151E-03/
      
      HVRCLD = '$Revision$'

      ICLDATM = 0
      NCBANDS = 1
      DO 3000 LAY = 1, NLAYERS
         IF (CLDFRAC(LAY) .GE. EPS) THEN
            ICLDATM = 1

C           Ice clouds and water clouds combined.
            IF (INFLAG .EQ. 0) THEN
               TAUCLOUD(LAY,1) = CLDDAT1(LAY)
            ELSEIF(INFLAG .EQ. 1) THEN
               CWP = CLDDAT1(LAY)
               TAUCLOUD(LAY,1) = ABSCLD1 * CWP

C           Separate treatement of ice clouds and water clouds.
            ELSEIF(INFLAG .EQ. 2) THEN
               CWP = CLDDAT1(LAY)
               FICE = CLDDAT2(LAY)
               RADICE = CLDDAT3(LAY)

C              Calculation of absorption coefficients due to ice clouds.
               IF (FICE .EQ. 0.0) THEN
                  ABSCOICE(1) = 0.0
                  ICEPAT = 0
               ELSEIF (ICEFLAG .EQ. 0) THEN
                  IF (RADICE .LT. 10.0) STOP 'ICE RADIUS TOO SMALL'
                  ABSCOICE(1) = FICE * (ABSICE0(1) + 
     &                 ABSICE0(2)/RADICE)
                  ICEPAT = 0
               ELSEIF (ICEFLAG .EQ. 1) THEN
                  IF (RADICE .LT. 13.0 .OR. RADICE .GT. 130.) STOP
     &                 'ICE RADIUS OUT OF BOUNDS'
                  NCBANDS = 5
                  DO 2000 IB = 1, NCBANDS
                     ABSCOICE(IB) = FICE * (ABSICE1(1,IB) + 
     &                    ABSICE1(2,IB)/RADICE)
 2000             CONTINUE
                  ICEPAT = 1
               ELSEIF (ICEFLAG .EQ. 2) THEN
                  IF (RADICE .LT. 13.0 .OR. RADICE .GT. 130.) STOP
     &                 'ICE RADIUS OUT OF BOUNDS'
                  NCBANDS = 16
                  FACTOR = (RADICE - 10.)/3.
                  INDEX = INT(FACTOR)
                  IF (INDEX .EQ. 40) INDEX = 39
                  FINT = FACTOR - FLOAT(INDEX)
                  DO 2200 IB = 1, NCBANDS
                     ABSCOICE(IB) = FICE * (ABSICE2(INDEX,IB) + FINT *
     &                    (ABSICE2(INDEX+1,IB) - (ABSICE2(INDEX,IB))))
 2200             CONTINUE
                  ICEPAT = 2
               ENDIF
                  
C              Calculation of absorption coefficients due to water clouds.
               FLIQ = 1. - FICE
               IF (FLIQ .EQ. 0.0) THEN
                  ABSCOLIQ(1) = 0.0
                  LIQPAT = 0
                  IF (ICEPAT .EQ. 1) ICEPAT = 2
               ELSEIF (LIQFLAG .EQ. 0) THEN
                  ABSCOLIQ(1) = FLIQ * ABSLIQ0
                  LIQPAT = 0
                  IF (ICEPAT .EQ. 1) ICEPAT = 2
               ELSEIF (LIQFLAG .EQ. 1) THEN
                  RADLIQ = CLDDAT4(LAY)
                  IF (RADLIQ .LT. 2.5 .OR. RADLIQ .GT. 60.) STOP
     &                 'LIQUID EFFECTIVE RADIUS OUT OF BOUNDS'
                  INDEX = RADLIQ - 1.5
                  IF (INDEX .EQ. 58) INDEX = 57
                  FINT = RADLIQ - 1.5 - INDEX
                  NCBANDS = 16
                  DO 2400 IB = 1, NCBANDS
                     ABSCOLIQ(IB) = FLIQ * (ABSLIQ1(INDEX,IB) + FINT *
     &                    (ABSLIQ1(INDEX+1,IB) - (ABSLIQ1(INDEX,IB))))
 2400             CONTINUE
                  LIQPAT = 2
               ENDIF

               DO 2800 IB = 1, NCBANDS
                  TAUCLOUD(LAY,IB) = CWP * (ABSCOICE(IPAT(IB,ICEPAT)) +
     &                 ABSCOLIQ(IPAT(IB,LIQPAT)))
 2800          CONTINUE
            ENDIF
         ENDIF
 3000 CONTINUE

      RETURN
      END
