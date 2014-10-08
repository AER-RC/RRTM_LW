C     path:      $Source$
C     author:    $Author: jdelamer $
C     revision:  $Revision: 11468 $
C     created:   $Date: 2004-07-29 08:49:51 -0800 (Thu, 29 Jul 2004) $
C
C  --------------------------------------------------------------------------
C |                                                                          |
C |  Copyright 2002, 2003, Atmospheric & Environmental Research, Inc. (AER). |
C |  This software may be used, copied, or redistributed as long as it is    |
C |  not sold and this copyright notice is reproduced on each copy made.     |
C |  This model is provided as is without any express or implied warranties. |
C |                       (http://www.rtweb.aer.com/)                        |
C |                                                                          |
C  --------------------------------------------------------------------------

      PARAMETER (MG=16)
      REAL KA(5,13,MG), KB(5,13:59,MG)
      DIMENSION SELFREF(10,MG), FORREF(4,MG)
      REAL KA_MN2(19,MG),KB_MN2(19,MG)
      DIMENSION FRACREFA(MG),FRACREFB(MG)

      COMMON /CVRSN1/ HNAMKG1,HVRKG1

      COMMON /K1/ KA ,KB, FORREF, SELFREF, KA_MN2, KB_MN2, 
     &   FRACREFA, FRACREFB

      CHARACTER*18 HVRKG1

      CHARACTER*18 HNAMKG1

      DATA HVRKG1 /'$Revision: 11468 $'/
      DATA HNAMKG1 / '         k_gB01.f:' /

      DATA (KA_MN2(JT, 1),JT=1,19)  /
     & 5.12042e-08, 5.51239e-08, 5.93436e-08, 6.38863e-08, 6.87767e-08,
     & 7.40415e-08, 7.97093e-08, 8.58110e-08, 9.23797e-08, 9.94513e-08,
     & 1.07064e-07, 1.15260e-07, 1.24083e-07, 1.33581e-07, 1.43807e-07,
     & 1.54815e-07, 1.66666e-07, 1.79424e-07, 1.93159e-07/
      DATA (KA_MN2(JT, 2),JT=1,19)  /
     & 2.30938e-07, 2.41696e-07, 2.52955e-07, 2.64738e-07, 2.77071e-07,
     & 2.89978e-07, 3.03486e-07, 3.17623e-07, 3.32419e-07, 3.47904e-07,
     & 3.64111e-07, 3.81072e-07, 3.98824e-07, 4.17402e-07, 4.36846e-07,
     & 4.57196e-07, 4.78494e-07, 5.00784e-07, 5.24112e-07/
      DATA (KA_MN2(JT, 3),JT=1,19)  /
     & 6.70458e-07, 7.04274e-07, 7.39795e-07, 7.77109e-07, 8.16304e-07,
     & 8.57476e-07, 9.00724e-07, 9.46154e-07, 9.93876e-07, 1.04400e-06,
     & 1.09666e-06, 1.15197e-06, 1.21008e-06, 1.27111e-06, 1.33522e-06,
     & 1.40256e-06, 1.47331e-06, 1.54761e-06, 1.62567e-06/
      DATA (KA_MN2(JT, 4),JT=1,19)  /
     & 1.84182e-06, 1.89203e-06, 1.94360e-06, 1.99658e-06, 2.05101e-06,
     & 2.10692e-06, 2.16435e-06, 2.22335e-06, 2.28396e-06, 2.34622e-06,
     & 2.41017e-06, 2.47587e-06, 2.54337e-06, 2.61270e-06, 2.68392e-06,
     & 2.75708e-06, 2.83224e-06, 2.90944e-06, 2.98875e-06/
      DATA (KA_MN2(JT, 5),JT=1,19)  /
     & 3.41996e-06, 3.32758e-06, 3.23770e-06, 3.15024e-06, 3.06515e-06,
     & 2.98235e-06, 2.90180e-06, 2.82341e-06, 2.74715e-06, 2.67294e-06,
     & 2.60074e-06, 2.53049e-06, 2.46214e-06, 2.39563e-06, 2.33092e-06,
     & 2.26796e-06, 2.20670e-06, 2.14709e-06, 2.08910e-06/
      DATA (KA_MN2(JT, 6),JT=1,19)  /
     & 3.38746e-06, 3.25966e-06, 3.13669e-06, 3.01836e-06, 2.90449e-06,
     & 2.79491e-06, 2.68947e-06, 2.58801e-06, 2.49037e-06, 2.39642e-06,
     & 2.30601e-06, 2.21902e-06, 2.13530e-06, 2.05475e-06, 1.97723e-06,
     & 1.90264e-06, 1.83086e-06, 1.76179e-06, 1.69532e-06/
      DATA (KA_MN2(JT, 7),JT=1,19)  /
     & 3.17530e-06, 3.07196e-06, 2.97199e-06, 2.87527e-06, 2.78170e-06,
     & 2.69118e-06, 2.60360e-06, 2.51887e-06, 2.43690e-06, 2.35759e-06,
     & 2.28087e-06, 2.20664e-06, 2.13483e-06, 2.06536e-06, 1.99814e-06,
     & 1.93312e-06, 1.87021e-06, 1.80934e-06, 1.75046e-06/
      DATA (KA_MN2(JT, 8),JT=1,19)  /
     & 2.84701e-06, 2.77007e-06, 2.69521e-06, 2.62237e-06, 2.55150e-06,
     & 2.48254e-06, 2.41545e-06, 2.35017e-06, 2.28666e-06, 2.22486e-06,
     & 2.16473e-06, 2.10623e-06, 2.04930e-06, 1.99392e-06, 1.94003e-06,
     & 1.88760e-06, 1.83659e-06, 1.78695e-06, 1.73866e-06/
      DATA (KA_MN2(JT, 9),JT=1,19)  /
     & 2.79917e-06, 2.73207e-06, 2.66658e-06, 2.60266e-06, 2.54027e-06,
     & 2.47937e-06, 2.41994e-06, 2.36192e-06, 2.30530e-06, 2.25004e-06,
     & 2.19610e-06, 2.14346e-06, 2.09208e-06, 2.04193e-06, 1.99298e-06,
     & 1.94520e-06, 1.89857e-06, 1.85306e-06, 1.80864e-06/
      DATA (KA_MN2(JT,10),JT=1,19)  /
     & 2.74910e-06, 2.64462e-06, 2.54412e-06, 2.44743e-06, 2.35442e-06,
     & 2.26495e-06, 2.17887e-06, 2.09606e-06, 2.01641e-06, 1.93978e-06,
     & 1.86606e-06, 1.79514e-06, 1.72692e-06, 1.66129e-06, 1.59815e-06,
     & 1.53742e-06, 1.47899e-06, 1.42278e-06, 1.36871e-06/
      DATA (KA_MN2(JT,11),JT=1,19)  /
     & 2.63952e-06, 2.60263e-06, 2.56626e-06, 2.53039e-06, 2.49503e-06,
     & 2.46016e-06, 2.42578e-06, 2.39188e-06, 2.35845e-06, 2.32549e-06,
     & 2.29299e-06, 2.26094e-06, 2.22934e-06, 2.19819e-06, 2.16747e-06,
     & 2.13717e-06, 2.10731e-06, 2.07786e-06, 2.04882e-06/
      DATA (KA_MN2(JT,12),JT=1,19)  /
     & 2.94106e-06, 2.82819e-06, 2.71966e-06, 2.61528e-06, 2.51492e-06,
     & 2.41841e-06, 2.32560e-06, 2.23635e-06, 2.15053e-06, 2.06800e-06,
     & 1.98863e-06, 1.91232e-06, 1.83893e-06, 1.76836e-06, 1.70049e-06,
     & 1.63524e-06, 1.57248e-06, 1.51214e-06, 1.45411e-06/
      DATA (KA_MN2(JT,13),JT=1,19)  /
     & 2.94607e-06, 2.87369e-06, 2.80309e-06, 2.73422e-06, 2.66705e-06,
     & 2.60152e-06, 2.53760e-06, 2.47526e-06, 2.41445e-06, 2.35513e-06,
     & 2.29726e-06, 2.24082e-06, 2.18577e-06, 2.13207e-06, 2.07969e-06,
     & 2.02859e-06, 1.97875e-06, 1.93014e-06, 1.88272e-06/
      DATA (KA_MN2(JT,14),JT=1,19)  /
     & 2.58051e-06, 2.48749e-06, 2.39782e-06, 2.31139e-06, 2.22807e-06,
     & 2.14775e-06, 2.07033e-06, 1.99570e-06, 1.92376e-06, 1.85441e-06,
     & 1.78756e-06, 1.72313e-06, 1.66101e-06, 1.60114e-06, 1.54342e-06,
     & 1.48778e-06, 1.43415e-06, 1.38245e-06, 1.33262e-06/
      DATA (KA_MN2(JT,15),JT=1,19)  /
     & 3.03447e-06, 2.88559e-06, 2.74401e-06, 2.60938e-06, 2.48135e-06,
     & 2.35961e-06, 2.24384e-06, 2.13375e-06, 2.02906e-06, 1.92951e-06,
     & 1.83484e-06, 1.74481e-06, 1.65921e-06, 1.57780e-06, 1.50039e-06,
     & 1.42677e-06, 1.35677e-06, 1.29020e-06, 1.22690e-06/
      DATA (KA_MN2(JT,16),JT=1,19)  /
     & 1.48655e-06, 1.48283e-06, 1.47913e-06, 1.47543e-06, 1.47174e-06,
     & 1.46806e-06, 1.46439e-06, 1.46072e-06, 1.45707e-06, 1.45343e-06,
     & 1.44979e-06, 1.44617e-06, 1.44255e-06, 1.43894e-06, 1.43534e-06,
     & 1.43176e-06, 1.42817e-06, 1.42460e-06, 1.42104e-06/
      DATA (KB_MN2(JT, 1),JT=1,19)  /
     & 5.12042e-08, 5.51239e-08, 5.93436e-08, 6.38863e-08, 6.87767e-08,
     & 7.40415e-08, 7.97093e-08, 8.58110e-08, 9.23797e-08, 9.94513e-08,
     & 1.07064e-07, 1.15260e-07, 1.24083e-07, 1.33581e-07, 1.43807e-07,
     & 1.54815e-07, 1.66666e-07, 1.79424e-07, 1.93159e-07/
      DATA (KB_MN2(JT, 2),JT=1,19)  /
     & 2.30938e-07, 2.41696e-07, 2.52955e-07, 2.64738e-07, 2.77071e-07,
     & 2.89978e-07, 3.03486e-07, 3.17623e-07, 3.32419e-07, 3.47904e-07,
     & 3.64111e-07, 3.81072e-07, 3.98824e-07, 4.17402e-07, 4.36846e-07,
     & 4.57196e-07, 4.78494e-07, 5.00784e-07, 5.24112e-07/
      DATA (KB_MN2(JT, 3),JT=1,19)  /
     & 6.70458e-07, 7.04274e-07, 7.39795e-07, 7.77109e-07, 8.16304e-07,
     & 8.57476e-07, 9.00724e-07, 9.46154e-07, 9.93876e-07, 1.04400e-06,
     & 1.09666e-06, 1.15197e-06, 1.21008e-06, 1.27111e-06, 1.33522e-06,
     & 1.40256e-06, 1.47331e-06, 1.54761e-06, 1.62567e-06/
      DATA (KB_MN2(JT, 4),JT=1,19)  /
     & 1.84182e-06, 1.89203e-06, 1.94360e-06, 1.99658e-06, 2.05101e-06,
     & 2.10692e-06, 2.16435e-06, 2.22335e-06, 2.28396e-06, 2.34622e-06,
     & 2.41017e-06, 2.47587e-06, 2.54337e-06, 2.61270e-06, 2.68392e-06,
     & 2.75708e-06, 2.83224e-06, 2.90944e-06, 2.98875e-06/
      DATA (KB_MN2(JT, 5),JT=1,19)  /
     & 3.41996e-06, 3.32758e-06, 3.23770e-06, 3.15024e-06, 3.06515e-06,
     & 2.98235e-06, 2.90180e-06, 2.82341e-06, 2.74715e-06, 2.67294e-06,
     & 2.60074e-06, 2.53049e-06, 2.46214e-06, 2.39563e-06, 2.33092e-06,
     & 2.26796e-06, 2.20670e-06, 2.14709e-06, 2.08910e-06/
      DATA (KB_MN2(JT, 6),JT=1,19)  /
     & 3.38746e-06, 3.25966e-06, 3.13669e-06, 3.01836e-06, 2.90449e-06,
     & 2.79491e-06, 2.68947e-06, 2.58801e-06, 2.49037e-06, 2.39642e-06,
     & 2.30601e-06, 2.21902e-06, 2.13530e-06, 2.05475e-06, 1.97723e-06,
     & 1.90264e-06, 1.83086e-06, 1.76179e-06, 1.69532e-06/
      DATA (KB_MN2(JT, 7),JT=1,19)  /
     & 3.17530e-06, 3.07196e-06, 2.97199e-06, 2.87527e-06, 2.78170e-06,
     & 2.69118e-06, 2.60360e-06, 2.51887e-06, 2.43690e-06, 2.35759e-06,
     & 2.28087e-06, 2.20664e-06, 2.13483e-06, 2.06536e-06, 1.99814e-06,
     & 1.93312e-06, 1.87021e-06, 1.80934e-06, 1.75046e-06/
      DATA (KB_MN2(JT, 8),JT=1,19)  /
     & 2.84701e-06, 2.77007e-06, 2.69521e-06, 2.62237e-06, 2.55150e-06,
     & 2.48254e-06, 2.41545e-06, 2.35017e-06, 2.28666e-06, 2.22486e-06,
     & 2.16473e-06, 2.10623e-06, 2.04930e-06, 1.99392e-06, 1.94003e-06,
     & 1.88760e-06, 1.83659e-06, 1.78695e-06, 1.73866e-06/
      DATA (KB_MN2(JT, 9),JT=1,19)  /
     & 2.79917e-06, 2.73207e-06, 2.66658e-06, 2.60266e-06, 2.54027e-06,
     & 2.47937e-06, 2.41994e-06, 2.36192e-06, 2.30530e-06, 2.25004e-06,
     & 2.19610e-06, 2.14346e-06, 2.09208e-06, 2.04193e-06, 1.99298e-06,
     & 1.94520e-06, 1.89857e-06, 1.85306e-06, 1.80864e-06/
      DATA (KB_MN2(JT,10),JT=1,19)  /
     & 2.74910e-06, 2.64462e-06, 2.54412e-06, 2.44743e-06, 2.35442e-06,
     & 2.26495e-06, 2.17887e-06, 2.09606e-06, 2.01641e-06, 1.93978e-06,
     & 1.86606e-06, 1.79514e-06, 1.72692e-06, 1.66129e-06, 1.59815e-06,
     & 1.53742e-06, 1.47899e-06, 1.42278e-06, 1.36871e-06/
      DATA (KB_MN2(JT,11),JT=1,19)  /
     & 2.63952e-06, 2.60263e-06, 2.56626e-06, 2.53039e-06, 2.49503e-06,
     & 2.46016e-06, 2.42578e-06, 2.39188e-06, 2.35845e-06, 2.32549e-06,
     & 2.29299e-06, 2.26094e-06, 2.22934e-06, 2.19819e-06, 2.16747e-06,
     & 2.13717e-06, 2.10731e-06, 2.07786e-06, 2.04882e-06/
      DATA (KB_MN2(JT,12),JT=1,19)  /
     & 2.94106e-06, 2.82819e-06, 2.71966e-06, 2.61528e-06, 2.51492e-06,
     & 2.41841e-06, 2.32560e-06, 2.23635e-06, 2.15053e-06, 2.06800e-06,
     & 1.98863e-06, 1.91232e-06, 1.83893e-06, 1.76836e-06, 1.70049e-06,
     & 1.63524e-06, 1.57248e-06, 1.51214e-06, 1.45411e-06/
      DATA (KB_MN2(JT,13),JT=1,19)  /
     & 2.94607e-06, 2.87369e-06, 2.80309e-06, 2.73422e-06, 2.66705e-06,
     & 2.60152e-06, 2.53760e-06, 2.47526e-06, 2.41445e-06, 2.35513e-06,
     & 2.29726e-06, 2.24082e-06, 2.18577e-06, 2.13207e-06, 2.07969e-06,
     & 2.02859e-06, 1.97875e-06, 1.93014e-06, 1.88272e-06/
      DATA (KB_MN2(JT,14),JT=1,19)  /
     & 2.58051e-06, 2.48749e-06, 2.39782e-06, 2.31139e-06, 2.22807e-06,
     & 2.14775e-06, 2.07033e-06, 1.99570e-06, 1.92376e-06, 1.85441e-06,
     & 1.78756e-06, 1.72313e-06, 1.66101e-06, 1.60114e-06, 1.54342e-06,
     & 1.48778e-06, 1.43415e-06, 1.38245e-06, 1.33262e-06/
      DATA (KB_MN2(JT,15),JT=1,19)  /
     & 3.03447e-06, 2.88559e-06, 2.74401e-06, 2.60938e-06, 2.48135e-06,
     & 2.35961e-06, 2.24384e-06, 2.13375e-06, 2.02906e-06, 1.92951e-06,
     & 1.83484e-06, 1.74481e-06, 1.65921e-06, 1.57780e-06, 1.50039e-06,
     & 1.42677e-06, 1.35677e-06, 1.29020e-06, 1.22690e-06/
      DATA (KB_MN2(JT,16),JT=1,19)  /
     & 1.48655e-06, 1.48283e-06, 1.47913e-06, 1.47543e-06, 1.47174e-06,
     & 1.46806e-06, 1.46439e-06, 1.46072e-06, 1.45707e-06, 1.45343e-06,
     & 1.44979e-06, 1.44617e-06, 1.44255e-06, 1.43894e-06, 1.43534e-06,
     & 1.43176e-06, 1.42817e-06, 1.42460e-06, 1.42104e-06/

C     The array FORREF contains the coefficient of the water vapor
C     foreign-continuum (including the energy term).  The first 
C     index refers to reference temperature (296,260,224,260) and 
C     pressure (970,475,219,3 mbar) levels.  The second index 
C     runs over the g-channel (1 to 16).

