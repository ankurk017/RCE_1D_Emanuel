C A1.0 Fundamental constants
      REAL RPI,RCLUM,RHPLA,RKBOL,RNAVO
C A1.1 Astronomical constants
      REAL RDAY,REA,REPSM,RSIYEA,RSIDAY,ROMEGA
c A1.1.bis Constantes concernant l'orbite de la Terre:
      REAL R_ecc, R_peri, R_incl
C A1.2 Geoide
      REAL RA,RG,R1SA
C A1.3 Radiation
      REAL RSIGMA,RI0
C A1.4 Thermodynamic gas phase
c sb      REAL R,RMD,RMV,RD,RV,RCPD,RCPV,RCVD,RCVV
      REAL RAP,RMD,RMV,RD,RV,RCPD,RCPV,RCVD,RCVV
      REAL RKAPPA,RETV
C A1.5,6 Thermodynamic liquid,solid phases
      REAL RCW,RCS
C A1.7 Thermodynamic transition of phase
      REAL RLVTT,RLSTT,RLMLT,RTT,RATM
C A1.8 Curve of saturation
      REAL RESTT,RALPW,RBETW,RGAMW,RALPS,RBETS,RGAMS
      REAL RALPD,RBETD,RGAMD
C
      COMMON/YOMCST/RPI   ,RCLUM ,RHPLA ,RKBOL ,RNAVO
     S      ,RDAY  ,REA   ,REPSM ,RSIYEA,RSIDAY,ROMEGA
     s      ,R_ecc, R_peri, R_incl
     S      ,RA    ,RG    ,R1SA
     S      ,RSIGMA,RI0
c sb     S      ,R     ,RMD   ,RMV   ,RD    ,RV    ,RCPD
     S      ,RAP   ,RMD   ,RMV   ,RD    ,RV    ,RCPD
     S      ,RCPV  ,RCVD  ,RCVV  ,RKAPPA,RETV
     S      ,RCW   ,RCS
     S      ,RLVTT ,RLSTT ,RLMLT ,RTT   ,RATM
     S      ,RESTT ,RALPW ,RBETW ,RGAMW ,RALPS ,RBETS ,RGAMS
     S      ,RALPD ,RBETD ,RGAMD
C    ------------------------------------------------------------------
