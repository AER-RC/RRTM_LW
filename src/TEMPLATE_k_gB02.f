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
      DIMENSION FRACREFA(MG), FRACREFB(MG)

      COMMON /CVRSN2/ HNAMKG2,HVRKG2

      COMMON /K2/ KA ,KB, FORREF, SELFREF, 
     &  FRACREFA, FRACREFB

      CHARACTER*18 HVRKG2

      CHARACTER*18 HNAMKG2

      DATA HVRKG2  / '$Revision: 11468 $' /
      DATA HNAMKG2 / '         k_gB02.f:' /




