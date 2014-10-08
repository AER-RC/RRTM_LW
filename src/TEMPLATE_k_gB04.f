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

      REAL KA(9,5,13,MG), KB(5,5,13:59,MG)
      DIMENSION SELFREF(10,MG), FORREF(4,MG), 
      DIMENSION FRACREFA(MG,9), FRACREFB(MG,5)

      COMMON /CVRSN4/ HNAMKG4,HVRKG4
      COMMON /K4/ KA, KB, FORREF, SELFREF, 
     &  FRACREFA, FRACREFB

      CHARACTER*18 HVRKG4

      CHARACTER*18 HNAMKG4

      DATA HVRKG4 /'$Revision: 11468 $'/
      DATA HNAMKG4 / '         k_gB04.f:' /
