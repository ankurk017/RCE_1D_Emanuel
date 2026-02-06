        PROGRAM CONRAD 
C******************************************************************************
C
C     This program calculates single column evolution under the influence
C of surface fluxes, radiation and convection. Radiative cooling may be 
C specified or calculated using a code which accounts for variability
C of water vapor and clouds. One has the ability to specify a vertical
C velocity profile, which may vary periodically in time. Input is
C from the file params.in, which gives the parameters used by the model
C and may be created using the separate routine parameter.f, and the file
C sounding.in, which may be created using the separate routine sounding.f.
C Output is directed to the files sounding.out, which is of the same format
C as sounding.in and may be used to restart the integration; error.out,
C which contains any error messages, and NCAR graphics metacode, which the
C user may replace with his own graphics if so desired.
C
C The main program and subroutine CONVECT were created by Kerry Emanuel,
C while the radiation code was developed by Chou et al. (J. Climate, 4,
C 1991, p. 424) and modified by Ridgway and N. Renno. The radiation subroutine
C may be deleted if radiative cooling is specified; this will result in faster
C compilation and running. 
C
C*****************************************************************************
C
C Note: The parameter NA must be greater than or equal to the total number 
C       of model levels at which the state variables are defined.  This will
C       equal the number of tropospheric levels plus 9. NOTE THAT NA MUST BE
C       SET TO THE SAME VALUE IN SUBROUTINE CONVECT. NTIME is at least the
C       number of time steps that certain quanitities like precipitation are
C	output to the graphics arrays.
C
      PARAMETER (NA=120, NAP=NA+1, NTIME=5000, NCL=1)
c	PARAMETER(NK=500,NI=42,NKRAD=2*NK)
C 
C  Dimension various arrays
C
C  Temperatures:
C
        REAL T1(NA), T2(NA), T3(NA), TTEMP(NA)
        REAL TRAD(NA), TVP(NA), TP(NA),TMEAN(NA)
	  REAL TINIT(NA), RINIT(NA)
C
C  Specific humidities:
C
        REAL R1(NA), R2(NA), R3(NA), RS(NA), RRAD(NA), RTEMP(NA) 
	  REAL RT(NA), RV2(NA) ! total water and water vapor (sb)
	  REAL RMEAN(NA)
C	  REAL RHTEMP(NA)
C
C  Other quantities in convection scheme:
C
        REAL U1(NA),U2(NA),U3(NA),V1(NA),V2(NA),V3(NA)
        REAL TRA1(NA,1),TRA2(NA,1),TRA3(NA,1)
C
C  Pressures:
C
        REAL P(NA), PH(NA), PHRAD(NAP)
C
C  Ozone used in radiative calculation:
C
        REAL O3(NA)
C
C  Specified mean pressure velocity:
C
        REAL OMEGA(NA)
        CHARACTER*1 VCUBE
C
C  Forcings:
C
        REAL FT(NA), FR(NA), FTNET(NA), FRNET(NA)
        REAL FTRA(NA,1),FU(NA),FV(NA),FUNET(NA),FVNET(NA),FTRANET(NA)
        REAL FTADJ(NA), FRADJ(NA) ! addendum sb
C
        LOGICAL ok, correc, correc2, rad1st, radcomp
     :        ,newobs,nodomega,newqlw,evap_prec
     :        ,correctd
C
C  Radiative arrays:
C
	REAL RCOOL(NA)
	REAL COOLR(NA), HEATR(NA), HOUR
	REAL HOUR0, HOURE
	REAL QCONDT(NA)
	CHARACTER*1 RADINT, TSINT, TDEP, DDEP, ANRAD, DARAD, CALB
	CHARACTER*1 FLUXINT,RESTART,ICLDS
	REAL RLAT(1), TS2(1)
C
      REAL MEAN,QSAT,CLDF(NA),CLDQ(NA),SIGSUB(NA),QSUBGRID
     :      ,gsol,minimum,xo,ho,pi,qmin,qmax,pas,signe,x
     :      ,delta(NA),TIMEL,essai
     :      ,RSNEWLS(NA), RNEWLS(NA), TNEWLS(NA)
     :      ,CLDFRAD(NA), CLDQRAD(NA)
     :      ,CLDFB(NA),CLDQB(NA)
C
C  Clouds and radiation 
C
        integer julien
        real CLDEMI(NA), CLDTAU(NA), CLDFICE(NA), CLDT, CLDWP
     :      ,zlongi, dist, RLON(1), gmtime, albsol(1), RCOTN(NA)
     :      ,pplay(NA), paprs(NAP), RCOTNOLD(NA), RCOTU(NA)
     :      ,RCLW(NA),RCSW(NA),RCLWOLD(NA),RCSWOLD(NA)
     :      ,RCLWU(NA),RCSWU(NA)
     :      ,fract(1),rmu0(1)
     :      ,heat(NA),heat0(NA),cool(NA),cool0(NA)
     :      ,radsol(1),albpla(1),topsw(1),toplw(1),solsw(1),sollw(1)
     :      ,topsw0(1),toplw0(1),solsw0(1),sollw0(1)
     :      ,radsolold,radsoln
c
c sb --
C
C  Other convective quantities used for diagnostics:
C
        REAL MP(NA), RP(NA), SIJ(NA,NA), WATER(NA),WT(NA)
        REAL FRSUB(NA)
        REAL MUP(NA), MDOWN(NA), ENT(NA), DET(NA), FRDET(NA)
C
C  Time-dependent quantities for graphics output:
C
        REAL RAIN_EVAP(NTIME,2), SST_TS(NTIME,2), GTIME(NTIME)
        REAL ENTHALPY(NTIME),  TQ(NTIME,2), GTIMEM(NTIME)
        REAL BTIME(NTIME,NA),  WTIME(NTIME,NA)
	  REAL SWTOPTIME(NTIME), LWTOPTIME(NTIME)
	  REAL VPOT(NTIME), PPOT(NTIME), TOUT(NTIME)
	  REAL CLDHOV(NTIME,NA),CTEMPBAR(NA),BTEMPBAR(NA)
	  REAL THOV(NTIME,NA),TTEMPBAR(NA)
	  REAL QHOV(NTIME,NA), QTEMPBAR(NA)
	  REAL RHHOV(NTIME,NA), RHTEMPBAR(NA)
	  REAL MHOV(NTIME,NA), MTEMPBAR(NA)
	  REAL MDHOV(NTIME,NA), MDTEMPBAR(NA)
	  REAL MPHOV(NTIME,NA), MPTEMPBAR(NA)
C
C  Mean quantities for output:
C
        REAL FTBAR(NA,6), FRBAR(NA,6), BUOYBAR(NA,2), RHBAR(NA)
        REAL HBAR(NA,2), MBAR(NA,3), MDET(NA,2), WATERBAR(NA,2)
        REAL TRAGRAPH(NA), GZ(NA)
	  REAL CFRACBAR(NA),CLQBAR(NA)
C
C  Latent heat at T=0 C
C
        REAL LV0
C
C  Common block with subroutine CONVECT:
C
        COMMON / CVT / TVP,MP,RP,SIJ,MUP,MDOWN,ENT,DET,FRDET,FRSUB,
     1     WATER, INB, TP,WT
C
C  Common block with subroutine CLOUDS_SUB_LS:
C
        REAL CLDL(NA),CLDA(NA),CLDK(NA),CLDS(NA)
        COMMON / cldpara / CLDL, CLDA, CLDK, CLDS
c -------------------------------------------------------
c OPTIONS OF THE SIMULATION:
c -------------------------------------------------------
c if radcomp=T: the rad tendency added to net tendencies is computed:
        radcomp = .TRUE. 
c maximum precipitation efficiency:
        EPMAX = 0.999
c minimum precipitation efficiency:
        EPMIN = 0.0  
C
        rad1st = .FALSE.
C
C  Open input and output files:
C
      OPEN(UNIT=8,FILE='params.in',STATUS='OLD')
      OPEN(UNIT=11,FILE='output/error.out',
     1   STATUS='UNKNOWN')
C
C  Read in the parameters of the model from the file params.in
C
      READ(8,10)RESTART,RADINT,ICLDS,TSINT,SCON,RLAT(1),MONTH,IDAY,
     1    HOUR,TDEP,DDEP,DARAD,ANRAD,CALB,ALB,CO2,DEPTH,VS0,
     2    VCUBE,PERIOD,VMAX,PVZERO,PVMAX,PFIX,DT,ENDTIME,AVTIME,
     3    RADFREQ,PRINFREQ
C
   10	FORMAT(//,4(T37,A1,/),2(T32,F12.6,/),2(T32,I6,/),T32,F12.6,
     1    /,5(T37,A1,/),4(T32,F12.6,/),T37,A1,10(/,T32,F12.6))
C
	IF(RESTART.EQ.'Y'.OR.RESTART.EQ.'y')THEN
	 OPEN(UNIT=9,FILE='output/sounding.out',STATUS='UNKNOWN')
	ELSE
       OPEN(UNIT=9,FILE='sounding.in',STATUS='OLD')
	END IF
c	RDAMP=4./(RADFREQ*3600.)
	RDAMP=0.0
C
C  Renormalize certain parameters
C
        IF(PERIOD.LE.1000.0)AVTIME=PERIOD
        PERIOD=PERIOD*24.*3600.
        VFREQ=2.*3.14159/MAX(PERIOD,1.0)
        IF(PERIOD.GT.(1000.0*24.*3600.))VFREQ=0.0
        VMAX=VMAX/3600.0
        PVZERO=MAX(PVZERO,100.0)
        DT=DT*60.0
        ENDTIME=ENDTIME*3600.0*24.0
        AVTIME=AVTIME*3600.0*24.0
        RADFREQ=RADFREQ*3600.0
        PRINFREQ=PRINFREQ*3600.0
	  IF(RADINT.EQ.'Y')RADINT='y'
	  IF(RADINT.EQ.'N')RADINT='n'
	  IF(RADINT.EQ.'n')radcomp=.FALSE.
	  IF(TSINT.EQ.'Y')TSINT='y'
	  IF(TSINT.EQ.'N')TSINT='n'
	  IF(ANRAD.EQ.'Y')ANRAD='y'
	  IF(ANRAD.EQ.'N')ANRAD='n'
	  IF(TDEP.EQ.'Y')TDEP='y'
	  IF(TDEP.EQ.'N')TDEP='n'
	  IF(DDEP.EQ.'Y')DDEP='y'
	  IF(DDEP.EQ.'N')DDEP='n'
	  IF(CALB.EQ.'Y')CALB='y'
	  IF(CALB.EQ.'N')CALB='n'
        IF(ANRAD.EQ.'y')THEN
         TDEP='n'
         DDEP='n'
        END IF
        IF(DARAD.EQ.'y')THEN
         TDEP='n'
        END IF
        MONTH0=MONTH
        IDAY0=IDAY
        HOUR0=HOUR
	IF(MONTH.EQ.1)THEN
	 JULIENDAY=IDAY
	ELSE IF(MONTH.EQ.2)THEN
	 JULIENDAY=31+IDAY
	ELSE IF(MONTH.EQ.3)THEN
	 JULIENDAY=31+28+IDAY
	ELSE IF(MONTH.EQ.4)THEN
	 JULIENDAY=31+28+31+IDAY
	ELSE IF(MONTH.EQ.5)THEN
	 JULIENDAY=31+28+31+30+IDAY
	ELSE IF(MONTH.EQ.6)THEN
	 JULIENDAY=31+28+31+30+31+IDAY
	ELSE IF(MONTH.EQ.7)THEN
	 JULIENDAY=31+28+31+30+31+30+IDAY
	ELSE IF(MONTH.EQ.8)THEN
	 JULIENDAY=31+28+31+30+31+30+31+IDAY
	ELSE IF(MONTH.EQ.9)THEN
	 JULIENDAY=31+28+31+30+31+30+31+31 +IDAY
	ELSE IF(MONTH.EQ.10)THEN
	 JULIENDAY=31+28+31+30+31+30+31+31+30+IDAY
	ELSE IF(MONTH.EQ.11)THEN
	 JULIENDAY=31+28+31+30+31+30+31+31+30+31+IDAY
	ELSE IF(MONTH.EQ.12)THEN
	 JULIENDAY=31+28+31+30+31+30+31+31+30+31+30+IDAY
	END IF
C
C  Calculate annual average earth-sun distance, solar zenith angle and albedo
C
	IF(ANRAD.EQ.'y')THEN
	 distm=0.0
	 rmumean=0.0
	 fracmean=0.0
	 albmean=0.0
	 DO JDAY=1,365
        CALL orbite(FLOAT(JDAY),zlongi,dist)
        CALL angle(zlongi,RLAT,fract,rmu0)
        CALL alboc(FLOAT(JDAY),RLAT,albsol)
	  distm=distm+dist
	  rmumean=rmumean+rmu0(1)
	  fracmean=fracmean+fract(1)
	  albmean=albmean+albsol(1)
	 END DO
	 distm=distm/365.
	 rmumean=rmumean/365.
	 fracmean=fracmean/365.
	 albmean=albmean/365.
	END IF
C
C  Read in the initial sounding from the file sounding.in
C
      READ(9,20)NP,CBMF,TS1
   20 FORMAT(T15,I3,T35,F8.3,T54,F6.2,////)
	CBMF=CBMF*0.001
      N=NP-5
      TS1=TS1+273.15
      TS2(1)=TS1
      PG=1012.5
      DO 40 I=1,NP
        READ(9,30)P(I),T1(I),R1(I),O3(I),OMEGA(I),RCOOL(I)
   30 FORMAT(3X,F6.1,T12,F8.3,T25,F10.6,T39,F6.3,T47,F6.2,
     1    T58,F7.3)
c
c	r1(i)=0.4*r1(i)*(p(i)*0.001)**2
c	r1(i)=0.1*r1(i)*(p(i)*0.001)
c	r1(i)=0.01*r1(i)
C
C  Reset mixing ratios to 0.002 gm/kg above 100 mb.   
C       This prevents spurious high level moisture when starting at 
C       colder temperatures.
C
C	 IF(P(I).LT.100.0)R1(I)=0.002
C
C  Renormalize certain quantities and produce reverse order arrays for
C       radiative subroutine
C
         O3(I)=1.0E-6*O3(I)
c
c	   TC=T1(I)
c         TK=T1(I)+273.15
c	   DENOM=MAX(0.01,(243.5+TC))
c	   EPS=287.04/461.5
c         IF(TC.GE.0.0)THEN
c          ES=6.112*EXP(17.67*TC/DENOM)
c         ELSE
c          ES=EXP(23.33086-6111.72784/TK+0.15215*LOG(TK))
c         END IF
c         ES=MIN(ES,P(I))
c         RSTEMP=EPS*ES/(P(I)-(1.-EPS)*ES)
c	   RHTEMP(I)=0.001*R1(I)/MAX(RSTEMP,1.0E-7)
c	   rhtemp(i)=min(rhtemp(i),1.0)
c
         R1(I)=0.001*R1(I)
         T1(I)=T1(I)+273.15
         T2(I)=T1(I)
         R2(I)=R1(I)
	   TINIT(I)=T2(I)
	   RINIT(I)=R2(I)
	   TMEAN(I)=T1(I)
	   RMEAN(I)=R1(I)
         OMEGA(I)=OMEGA(I)/3600.0
         RCOOL(I)=RCOOL(I)/(24.*3600.)
         IF(I.GT.1)THEN
          PH(I)=0.5*(P(I)+P(I-1))
          PHRAD(NP+2-I)=PH(I)
         END IF
         U1(I)=30.*(P(1)-P(I))/P(1)
         U2(I)=U1(I)
         V1(I)=20.*(P(1)-P(I))/P(1)
         V2(I)=V1(I)
         TRA1(I,1)=0.0
         TRA2(I,1)=TRA1(I,1)
   40	CONTINUE 
        PH(1)=2.*P(1)-PH(2)
        PHRAD(NP+1)=PH(1)
        PH(NP+1)=2.*PH(NP)-PH(NP-1)
        PH(NP+1)=MAX(0.1,PH(NP+1))
        PHRAD(1)=PH(NP+1)
C
        CLOSE(8)
        CLOSE(9)
c
c	 T1(21)=T1(21)-1.
c	 T1(22)=T1(22)-2.
c	 T1(23)=T1(23)-3.
c	 T1(24)=T1(24)-2.
c	 T1(25)=T1(25)-1.
c	DO I=21,25
c	 T2(I)=T1(I)
c	END DO
C
C  If non-interactive clouds, read in time-average cloud fraction and water from last run
C
	IF(ICLDS.EQ.'n'.OR.ICLDS.EQ.'N')THEN
C
	 OPEN(UNIT=9,FILE='output/profile.out',STATUS='OLD')
	 DO I=1,NP
	  READ(9,43)CLDFB(I),CLDQB(I)
	  CLDQB(I)=0.001*CLDQB(I)
	 END DO
	 CLOSE(9)
C
	END IF
   43 FORMAT(151X,F9.3,1X,E10.2)
C
C  If so specified, calculate cubic omega profile
C
        IF(VCUBE.EQ.'y'.OR.VCUBE.EQ.'Y')THEN
         XM=PH(1)-PVMAX
         XT=PH(1)-PVZERO
         A=-VMAX*(3.*XM-2.*XT)/(XM*(XT-XM)**2)
         B=-VMAX*(XT-2.*XM)/(XM*XM*(XT-XM)**2)
         DO 45 I=1,N
          X=PH(1)-PH(I)
          OMEGA(I)=(XT-X)*(A*X+B*X*X)
          IF(P(I).LT.PVZERO)OMEGA(I)=0.0
   45	 CONTINUE
        END IF
C
C  Assign thermodynamic constants
C
        CPD=1005.7
        CPV=1870.0
        CL0=2500.0
        CL=CL0
        RD=287.04
        RV=461.5
        LV0=2.501E6
        G=9.8
        ROWL=1000.0
        EPS=RD/RV
        EPSI=1./EPS
        CPVMCL=CL-CPV
	  SIGD=0.05
C
C  Asselin filter constant
C
        AFILT=0.3
c
c   *** Stratospheric Dobson circulation   ***
c
c	DO 47 I=1,NP
c	 IF(P(I).LT.200.0.AND.P(I).GT.8.0)THEN
c	  OMEGA(I)=-5.0E-4*P(I)*G/(RD*T2(I))
c	 END IF
c   47	CONTINUE
C
C  Initialize certain quantities and arrays
C
	  ICB=2                                     
        FTS=0.0
C
        DO 60 I=1,NP
         FTNET(I)=0.0
         FRNET(I)=0.0
         FUNET(I)=0.0
         FVNET(I)=0.0
         FTRANET(I)=0.0
         RHBAR(I)=0.0
         COOLR(I)=0.0
         HEATR(I)=0.0
	   CFRACBAR(I)=0.0
	   CLQBAR(I)=0.0
	   QCONDT(I)=0.0
	   CLDF(I)=0.0
	   CLDQ(I)=0.0
	   SIGSUB(I)=0.0
	   CTEMPBAR(I)=0.0
	   BTEMPBAR(I)=0.0
	   TTEMPBAR(I)=0.0
	   QTEMPBAR(I)=0.0
	   RHTEMPBAR(I)=0.0
	   MTEMPBAR(I)=0.0
	   MDTEMPBAR(I)=0.0
	   MPTEMPBAR(I)=0.0
         DO 52 J=1,2
          HBAR(I,J)=0.0
          MDET(I,J)=0.0
          WATERBAR(I,J)=0.0
          BUOYBAR(I,J)=0.0
   52	 CONTINUE
         DO 54 J=1,3
          MBAR(I,J)=0.0
   54	 CONTINUE
         DO 56 J=1,4
          FTBAR(I,J)=0.0
          FRBAR(I,J)=0.0
   56	 CONTINUE
   60	CONTINUE
      FS=0.0
	FL=0.0
C
C  Initialize counters
C
        TIME=HOUR*3600.
        RADTIME=0.0
        PRINTIME=0.0
        NT=0
        NPR=0
        NAV=0	
        PREBAR=0.0
        NPBAR=0
	  TPERR=0.0
        OLRERR=0.
        ASRERR=0.
        OLRBIA=0.
        ASRBIA=0.
        OLRNB=0.
        ASRNB=0.
c -- sb
        NRHI=0
        NTPI=0
        RHERRI=0.
        TPERRI=0.
c
	  nradc=0
	  topbar=0.0
	  botbar=0.0
	  vmax=64.0
C
C  Begin time loop at next statement
C
  100	CONTINUE
C
        TIME=TIME+DT
        RADTIME=RADTIME+DT
        PRINTIME=PRINTIME+DT
        NT=NT+1
        NPBAR=NPBAR+1
c        
c        if(nt.eq.3)then
c          do i=1,np
c            if(p(i).le.860.0.and.p(i).ge.390.0)then
c               r2(i)=0.5*r2(i)
c               r1(i)=0.5*r1(i)
c            end if
c          end do
c        end if        
C
C  End integration when TIME = ENDTIME
C
        IF(TIME.GT.ENDTIME)GOTO 200
C
C  Calculate saturation specific humidity and some temporary arrays
C
        TC=TS2(1)-273.15
        ES=6.112*EXP(17.67*TC/(243.5+TC))
        RSS=0.98*EPS*ES/(PG-(1.-EPS)*ES)
        PRADJ=0.0
        FRADJC=0.0
	  RSATS=RSS
C
        DO 110 I=1,NP
         FTNET(I)=0.0
         FRNET(I)=0.0
         FUNET(I)=0.0
         FVNET(I)=0.0
         FTRANET(I)=0.0
         TC=T2(I)-273.15
         DENOM=243.5+TC
         DENOM=MAX(DENOM,1.0)
C
C  Calculate saturation specific humidity, assuming water saturation
C    for T> 0C and ice saturation otherwise
C
         IF(TC.GE.0.0)THEN
          ES=6.112*EXP(17.67*TC/DENOM)
         ELSE
          ES=EXP(23.33086-6111.72784/T2(I)+0.15215*LOG(T2(I)))
         END IF
         ES=MIN(ES,P(I))
         RS(I)=EPS*ES/(P(I)-(1.-EPS)*ES)
         TTEMP(I)=T2(I)
         RTEMP(I)=R2(I)
C
  110	CONTINUE
C
C  Call convection subroutine
C
        TG=TS2(1)
        QG=RSS
        CALL CONVECT(T2,R2,RS,FS,FL,U2,V2,TRA2,P,PH,NA,N,0,DT,IFLAG,FT,
     1 FR,FU,FV,FTRA,PRECIP,WD,TPRIME,QPRIME,QCONDT,CBMF,ENTPRO,ENTPRO2)
C
C  If necessary, write to error file
C
        IF(IFLAG.EQ.4)THEN
         WRITE(11,120)NT
  120	 FORMAT(5X,'CFL condition violated at time step ',I8)
        END IF
C
C  Add convective tendencies to net tendencies
C
        INBMAX=N
        SUM1=0.0
        DO 125 I=1,N
	   ALV=LV0-CPVMCL*(T2(I)-273.15)
	   ALV1=LV0-CPVMCL*(TTEMP(I)-273.15) 
	   SUM1=SUM1+((CPD*(1.-R2(I))+CL*R2(I))*T2(I)+ALV*R2(I)-
     1    (CPD*(1.-RTEMP(I))+CL*RTEMP(I))*TTEMP(I)-ALV1*RTEMP(I))*
     2    100.*(PH(I)-PH(I+1))/(2.*DT*G)
         FTNET(I)=FTNET(I)+FT(I)
         FRNET(I)=FRNET(I)+FR(I)
         FUNET(I)=FUNET(I)+FU(I)
         FVNET(I)=FVNET(I)+FV(I)
         FTRANET(I)=FTRANET(I)+FTRA(I,1)
         IF(P(I).GT.248.0)INBMAX=I
  125	CONTINUE
        INB=MAX(INB,INBMAX)
        FRCADJ=SUM1
C
C  *** Relax water vapor above tropopause to tropopause value  ***
C
        DO 126 I=INB+1,NP
c	 RTINB=MIN(R1(INB),RS(INB),RS(I))
         RTINB=1.0E-6
         FRNET(I)=FRNET(I)+2.0E-6*(RTINB*P(I)/P(INB)-R1(I))
  126   CONTINUE
C
C -- sb:
C -------------------------------------------------------------------
C
C  Compute the cloud fraction and the cloud water content:
C
C    1) make supersaturation adjustment
C       (precipitate only a fraction of the supersaturation)
C
C    2) compute the cloud fraction using the in-cloud water
C       content derived from the convection scheme (the subgrid-scale
C       part of the in-cloud water) and that derived from the
C       large-scale condensation scheme (or supersaturation sdjustment).
C
C -------------------------------------------------------------------
C
         CALL CLOUDS_SUB_LS_40(NP,R2,RS,T2,P,PH,DT,qcondt
     :              ,CLDF,CLDQ,PRADJ,FTADJ,FRADJ,SIGSUB)
C
C  If non-interactive clouds, set cloud fraction and water to mean values from last run
C
	IF(ICLDS.EQ.'n'.OR.ICLDS.EQ.'N')THEN
C
	   DO I=1,NP
	    CLDF(I)=CLDFB(I)
	    CLDQ(I)=CLDQB(I)
	   END DO
C
	END IF	 
C
C  Add tendencies due to supersaturation adjustment:
C
         DO 111 I = 1, NP
C
          FRNET(I)=FRNET(I)+FRADJ(I)
          FTNET(I)=FTNET(I)+FTADJ(I)
	    CTEMPBAR(I)=CTEMPBAR(I)+CLDF(I)
	    BTEMPBAR(I)=BTEMPBAR(I)+TVP(I)-T2(I)*(1.+R2(I)*EPSI-R2(I))
C
          IF(TIME.GT.(ENDTIME-AVTIME))THEN
           FTBAR(I,1)=FTBAR(I,1)+FTADJ(I)
           FRBAR(I,2)=FRBAR(I,2)+FRADJ(I)
           FTBAR(I,6)=FTBAR(I,6)+FTADJ(I) 
	     CFRACBAR(I)=CFRACBAR(I)+CLDF(I)
	     CLQBAR(I)=CLQBAR(I)+CLDQ(I)
          END IF
C
  111   CONTINUE
C
C  Add precip from supersaturation adjustment to convective precip:
C
        PRECIP=PRECIP+PRADJ
C
C Accumulate precipitation over time interval PRINTIME for time
C      series graphical output
C
        PREBAR=PREBAR+PRECIP
C
C -------------------------------------------------------------------
C
C  Compute cloud optical properties:
C
C -------------------------------------------------------------------
c pressure: mb -> Pa:
        DO I = 1, NP
         paprs(I) = PH(I)*100.
         pplay(I) = P(I)*100.
        ENDDO
        paprs(NP+1) = PH(NP+1)*100.
C
        DO I = 1, NP
          CLDFRAD(I) = CLDF(I)
          CLDQRAD(I) = CLDQ(I)
        ENDDO

       CALL OPTICAL(NP,T2,paprs,pplay,CLDFRAD,CLDQRAD
     :             ,CLDEMI,CLDTAU,CLDFICE,CLDT,CLDWP)
	  	
C -------------------------------------------------------------------
C
C  Prepare and call radiation:
C
C -------------------------------------------------------------------
C
        IF(RADINT.EQ.'y')THEN
          IF(TIME.LT.(1.5*DT).OR.RADTIME.GT.(RADFREQ-10.0))THEN
          RADTIME=0.0
c -- memory of previous radiative cooling rates:
          if (rad1st) then
            radsolold=radsoln
           do I=1,NP
            RCOTNOLD(I)=RCOTN(I)
            RCLWOLD(I)=RCLW(I)
            RCSWOLD(I)=RCSW(I)
           enddo
          endif

c -- local hour (RLON=longitude):
c
C        RLON(1) = 156.
	  RLON(1)=0.0
        TIMEL = TIME + RLON(1)/360.*24.*3600.
C
C --  Use annual means, if specified:
C
	IF(ANRAD.EQ.'y')THEN
	 dist=distm
	 rmu0(1)=rmumean
	 fract(1)=fracmean
	 albsol(1)=albmean
	ELSE
c
c -- Earth-Sun distance:
c
	IF(DDEP.EQ.'y')THEN
        julien = JULIENDAY + INT(TIMEL)/(3600*24)
	ELSE
	  julien = JULIENDAY
	END IF
c
        CALL orbite(FLOAT(julien),zlongi,dist)

c -- solar zenith angle and surface albedo:
c
        IF(DARAD.EQ.'y')THEN ! no diurnal cycle
          CALL angle(zlongi,RLAT,fract,rmu0)
          CALL alboc(FLOAT(julien),RLAT,albsol)
        ELSE
	   IF(TDEP.EQ.'y')THEN
          gmtime = MOD(TIME,3600.*24.)/86400. - 0.5*RADFREQ/86400.
          CALL zenang(zlongi,gmtime,RADFREQ,RLAT,RLON,rmu0,fract)
	   ELSE
	    gmtime=HOUR/24.-0.05
          CALL zenang(zlongi,gmtime,0.1,RLAT,RLON,rmu0,fract)
	   END IF
         CALL alboc_cd(rmu0,albsol)
        ENDIF
C
	END IF
c
c -- call radiation:

	  IF(CALB.EQ.'n')albsol(1)=ALB 
        DO I = 1, NP
          RV2(I) = R2(I) - CLDQ(I)*CLDF(I)
          RV2(I) = MAX( MIN(RV2(I),R2(I)) , 0.0 )
c	   if(i.ge.1.and.i.le.5)rv2(i)=0.8*rv2(i)
	    TMEAN(I)=T2(I)
	    RMEAN(I)=R2(I)
        ENDDO
c
c	 CO2=360.0*(VMAX/64.0)**2
c
      CALL radlwsw 
     e            (dist, rmu0, fract, CO2, SCON,
     e             paprs, pplay,TS2,albsol,T2,RV2,O3,
     e             CLDFRAD, CLDEMI, CLDTAU,
     s             heat,heat0,cool,cool0,radsol,albpla,
     s             topsw,toplw,solsw,sollw,
     s             topsw0,toplw0,solsw0,sollw0)
c
	topbar=topbar+toplw(1)
	botbar=botbar+solsw(1)+sollw(1)
	nradc=nradc+1
C
C  Calculate net radiative cooling rate
C
          radsoln=radsol(1)
          DO 130 I=1,NP
           RCOTN(I)=(cool(I)-heat(I))/(24.*3600.)
           RCLW(I)=cool0(I)-cool(I) ! diagnostic only
           RCSW(I)=heat(I)-heat0(I) ! diagnostic only
c thermo constants used in Morcrette's radiation are different from those
c used in the calling program, so we renormalize the cooling rates:
           RCOTN(I)=RCOTN(I)*1004.709/(CPD*(1.-R2(I))+CPV*R2(I))
           RCOTN(I)=RCOTN(I)*G/9.80665
  130     CONTINUE

          if (radcomp) rad1st = .TRUE.

         END IF
        END IF

C sb --
C
C  Add radiative cooling to net heating rate
C
       radsolu=solsw(1)+sollw(1)
      DO I = 1, NP
       RCOTU(I)=RCOTN(I)
       RCLWU(I)=RCLW(I)
       RCSWU(I)=RCSW(I)
      ENDDO

	DO 135 I=1,NP   
        if (radcomp) then
	  FTNET(I)=FTNET(I)-RCOTU(I)
	  RCOOL(I)=RCOTU(I)
	  FTNET(I)=FTNET(I)-RDAMP*(T1(I)-TMEAN(I))
	  FRNET(I)=FRNET(I)-RDAMP*(R1(I)-RMEAN(I))
        else
	 IF(P(I).GT.95.0)THEN
c -- sb:
cc	  FTNET(I)=FTNET(I)-RCOT   
cc	  RCOOL(I)=RCOT
c	  FTNET(I)=FTNET(I)-RCOTU(I)
	  FTNET(I)=FTNET(I)-RCOOL(I)
c	  RCOOL(I)=RCOTU(I)
        FTRANET(I)=FTRANET(I)-TRA1(I,1)/(36.*3600.)
c sb --
         ELSE
	  RCOOL(I)=0.0
         END IF
        endif ! radcomp
  135	CONTINUE
C
C  Calculate surface flux forcing
C
        CD=1.2E-3
        TSA=T2(1)*(PG/P(1))**(RD/CPD)
        ROWS=PG/(RD*TSA*(1.+R2(1)*(1.*EPSI-1.)))
        VS1=SQRT(VS0*VS0+WD*WD)
C
C  Calculate v' of convective downdrafts for use in surface fluxes
C
        VSPRIME=(VS1-VS0)
        VSURF=VS1
        FTSURF=G*ROWS*CD*(VSURF*(TS1-TSA)-VSPRIME*TPRIME)
     1    /(PH(1)-PH(2))
	  FS=100.*ROWS*CD*CPD*(VSURF*(TS1-TSA)-VSPRIME*TPRIME)
        FTNET(1)=FTNET(1)+FTSURF
        FRSURF=G*ROWS*CD*(VSURF*(RSS-R1(1))-VSPRIME*QPRIME)
     1    /(PH(1)-PH(2))
	  FL=100.*2.5E6*ROWS*CD*(VSURF*(RSS-R1(1))-VSPRIME*QPRIME)
        FRNET(1)=FRNET(1)+FRSURF
        FTRANET(1)=FTRANET(1)+G*ROWS*CD*VSURF*(1.-TRA1(1,1))/
     1   (PH(1)-PH(2))
C
C  Calculate forcing of ocean temperature
C
        IF(TSINT.EQ.'y')THEN
         TC=T2(1)-273.15
         ALV=LV0-CPVMCL*TC
         SFLUX=(CPD*(1.-R2(1))+CPV*R2(1))*100.0*ROWS*CD*(VSURF*
     1      (TS1-TSA)-VSPRIME*TPRIME)+
     2      ALV*100.0*ROWS*CD*(VSURF*(RSS-R1(1))-VSPRIME*QPRIME)
C
C  Input interactive net surface radiative flux to ocean
C
         IF(RADINT.EQ.'y')THEN
          RADFLUX=radsolu
c          XSH = (CPD*(1.-R1(1))+CPV*R1(1))*100.0*ROWS*CD*VSURF*(TS1-TSA)
c          XLH = ALV*100.0*ROWS*CD*(VSURF*(RSS-R1(1))-VSPRIME*(RP(1)-R2(1)))
	   ELSE
C
C  If radiative cooling is specified, calculate net surface radiative
C   flux that is consistent with vertically integrated cooling
C
          RADFLUX=0.0
          DO 2002 I=1,NP
           RADFLUX=RADFLUX+(CPD*(1.-R1(I))+CPV*R1(I))*RCOOL(I)*
     1      100.*(PH(I)-PH(I+1))/G
 2002     CONTINUE 
         END IF
         FTS=(RADFLUX-SFLUX)/(DEPTH*ROWL*CL0)
c
c  Add in feedback from hurricanes
c
c         FTS=(RADFLUX-SFLUX-2.0E-10*VMAX**6)/(DEPTH*ROWL*CL0)
c
        END IF
C
C  Calculate vertical advection and add to net tendencies
C
        TFAC=COS(VFREQ*TIME)
        VVP=MAX(OMEGA(2),0.0)*TFAC
        ALPHA=RD*T2(1)*(1.+R2(1)*EPSI)/(P(1)*(1.+R2(1)))
        FTADV=VVP*ALPHA/CPD-VVP*(T2(1)-T2(2))/(P(1)-P(2))
        FTNET(1)=FTNET(1)+FTADV
        FRADV=-VVP*(R2(1)-R2(2))/(P(1)-P(2))
        FRNET(1)=FRNET(1)+FRADV
        DO 140 I=2, NP-1
         VVP=MAX(OMEGA(I+1),0.0)*TFAC
         VVM=MIN(OMEGA(I),0.0)*TFAC
         ALPHA=RD*T2(I)*(1.+R2(I)*EPSI)/(P(I)*(1.+R2(I)))
         FTADV=(VVP+VVM)*ALPHA/CPD-
     1    VVP*(T2(I)-T2(I+1))/(P(I)-P(I+1))-
     2	  VVM*(T2(I-1)-T2(I))/(P(I-1)-P(I))
C
         FTNET(I)=FTNET(I)+FTADV
C
         FRADV=-VVP*(R2(I)-R2(I+1))/(P(I)-P(I+1))-    
     1	   VVM*(R2(I-1)-R2(I))/(P(I-1)-P(I))
C
         FRNET(I)=FRNET(I)+FRADV
C
C  Calculate certain average quantities for output graphs
C
         IF(TIME.GT.(ENDTIME-AVTIME-10.0))THEN
          FTBAR(I,4)=FTBAR(I,4)+FTADV
          FRBAR(I,4)=FRBAR(I,4)+FRADV
         END IF
  140	CONTINUE
C
C  Create smoothed quantities for time-height sections
C
	DO I=1,NP
	 TTEMPBAR(I)=TTEMPBAR(I)+T2(I)
	 QTEMPBAR(I)=QTEMPBAR(I)+1000.*R2(I)
	 RHTEMPBAR(I)=RHTEMPBAR(I)+R2(I)/(MAX(RS(I),1.0E-6))
	 MTEMPBAR(I)=MTEMPBAR(I)+MUP(I)
	 MDTEMPBAR(I)=MDTEMPBAR(I)-MDOWN(I)
	 MPTEMPBAR(I)=MPTEMPBAR(I)-MP(I)
	END DO
C
C  Put values into time arrays
C
        IF(PRINTIME.GT.(PRINFREQ-10.0))THEN
	   FINVNP=1./FLOAT(NPBAR)
         PRINTIME=0.0
         NPR=NPR+1
         GTIME(NPR)=TIME/(3600.*24.)
         GTIMEM(NPR)=(TIME-0.5*PRINFREQ)/(3600.*24.)
	   SWTOPTIME(NPR)=topsw(1)
	   LWTOPTIME(NPR)=toplw(1)
c	   SWTOPTIME(NPR)=solsw(1)
c	   LWTOPTIME(NPR)=sollw(1)
C
         RAIN_EVAP(NPR,1)=PREBAR*FINVNP
         PREBAR=0.0
         RAIN_EVAP(NPR,2)=100.*ROWS*CD*(VSURF*
     1    (RSS-R1(1))-VSPRIME*QPRIME)*1000.0*3600.0*24./ROWL
C
         SST_TS(NPR,1)=TSA-273.15
         SST_TS(NPR,2)=TS1-273.15
         TQ(NPR,1)=T2(21)-273.
         TQ(NPR,2)=1000.*R2(21)
         ENTS=0.0
	   TS1C=TS1-273.15
	   DO I=1,NP
	    CLDHOV(NPR,I)=CTEMPBAR(I)*FINVNP
	    CTEMPBAR(I)=0.0
	    BTIME(NPR,I)=BTEMPBAR(I)*FINVNP
	    BTEMPBAR(I)=0.0
	    THOV(NPR,I)=TTEMPBAR(I)*FINVNP
	    TTEMPBAR(I)=0.0
	    QHOV(NPR,I)=QTEMPBAR(I)*FINVNP
	    QTEMPBAR(I)=0.0
	    RHHOV(NPR,I)=100.*RHTEMPBAR(I)*FINVNP
	    rhhov(npr,i)=max(rhhov(npr,i),-99.9)
	    rhhov(npr,i)=min(rhhov(npr,i),100.0)
	    RHTEMPBAR(I)=0.0
	    MHOV(NPR,I)=100.*MTEMPBAR(I)*FINVNP
	    MTEMPBAR(I)=0.0
	    MDHOV(NPR,I)=100.*MDTEMPBAR(I)*FINVNP
	    MDTEMPBAR(I)=0.0
	    MPHOV(NPR,I)=100.*MPTEMPBAR(I)*FINVNP
	    MPTEMPBAR(I)=0.0
	   END DO
         NPBAR=0
C
         CALL PCMIN(TS1C,PG,P,T2,R2,NA,N,PMIN,VMAX,IFL,TOOUT)
	   VPOT(NPR)=VMAX
	   PPOT(NPR)=PMIN
	   TOUT(NPR)=TOOUT
C
         DO 138 I=1,NP
          ALV=LV0-CPVMCL*(T2(I)-273.15)
          ENTS=ENTS+((CPD*(1.-R2(I))+CPV*R2(I))*FTNET(I)+ALV
     1     *FRNET(I))*(PH(I)-PH(I+1))
C	   BTIME(NPR,I)=TVP(I)-T2(I)*(1.+R2(I)*EPSI-R2(I))
         WTIME(NPR,I)=R2(I)
         RHM=R2(I)/RS(I)
         RHM=MIN(RHM,1.0)
  138	 CONTINUE
         ENTS=ENTS+RADFLUX-SFLUX
         ENTHALPY(NPR)=ENTS*3600.*24./(CPD*(PH(1)-PH(N+1)))
        END IF
C
C  Calculate certain quantities averaged at end of integration
C
        IF(TIME.GT.(ENDTIME-AVTIME))THEN
         NAV=NAV+1
         FTBAR(1,3)=FTBAR(1,3)+FTSURF
         FRBAR(1,3)=FRBAR(1,3)+FRSURF
         COEF=0.5*(1.-2.*AFILT)/DT
         GZ(1)=0.0	
         ALV=LV0-CPVMCL*(T2(1)-273.15)
         DHH=(CPD*(1.-R2(1))+CL*R2(1))*(T2(1)-250.0)+ALV*R2(1)
         DO 150 I=1,NP
          FTBAR(I,1)=FTBAR(I,1)+FT(I)
          FTBAR(I,3)=FTBAR(I,3)+COEF*(T2(I)-TTEMP(I))
          FTBAR(I,2)=FTBAR(I,2)-RCOOL(I)
          FRBAR(I,1)=FRBAR(I,1)+FRDET(I)
          FRBAR(I,2)=FRBAR(I,2)+FRSUB(I)
          FRBAR(I,3)=FRBAR(I,3)+COEF*(R2(I)-RTEMP(I))
          BUOYBAR(I,1)=BUOYBAR(I,1)+TVP(I)-T2(I)*(1.+R2(I)*EPSI-R2(I))
c	  BUOYBAR(I)=BUOYBAR(I)+TP(I)*(1.+RS(I)*EPSI-R2(1))-
c     1     T2(I)*(1.+R2(I)*EPSI-R2(I))
          RHBAR(I)=RHBAR(I)+MAX(R2(I),0.0)/MAX(RS(I),1.0E-6)
          IF(I.GT.1)GZ(I)=GZ(I-1)+0.5*RD*(T2(I-1)*(1.+R2(I-1)*
     1      EPSI-R2(I-1))+T2(I)*(1.+R2(I)*EPSI-R2(I)))*
     2	    (P(I-1)-P(I))/PH(I)
          ALV=LV0-CPVMCL*(T2(I)-273.15)
          HBAR(I,1)=HBAR(I,1)+(CPD*(1.-R2(I))+CL*R2(I))*(T2(I)
     1	   -250.0)+ALV*R2(I)+GZ(I)
          DHS=(CPD*(1.-R2(I))+CL*R2(I))*(T2(I)-250.0)+ALV*RS(I)+GZ(I)
          HBAR(I,2)=HBAR(I,2)+DHS
c          IF(I.GT.1)THEN
c           BUOYBAR(I,2)=BUOYBAR(I,2)+P(I)*(TP(I-1)-TP(I+1))*(DHH-DHS)/
c     1      (RD*T2(I)*(P(I-1)-P(I+1)))
c          END IF
          MBAR(I,1)=MBAR(I,1)+MUP(I)
          MBAR(I,2)=MBAR(I,2)-MDOWN(I)
          MBAR(I,3)=MBAR(I,3)-MP(I)
          MDET(I,1)=MDET(I,1)+ENT(I)
          MDET(I,2)=MDET(I,2)+DET(I)
          WATERBAR(I,1)=WATERBAR(I,1)+WATER(I)*SIGD*WT(I)
          WATERBAR(I,2)=WATERBAR(I,2)+RP(I+1)-R2(I)
  150	 CONTINUE
        END IF
C
C  Advance quantities one time step, applying Asselin filter
C
        IF(TSINT.EQ.'y')THEN
         TS3=TS1+2.*DT*FTS
         TS1=TS2(1)+AFILT*(TS1+TS3-2.*TS2(1))
         TS2(1)=TS3
        END IF
C
        DO 160 I=1,NP
         IF(P(I).GE.PFIX)THEN
c         
c          if(p(i).lt.890.0.and.p(i).gt.85.0)then
c            ALPHA=RD*T2(I)*(1.+R2(I)*EPSI)/(P(I)*(1.+R2(I)))
c            alpha=alpha/cpd
c            if(ftnet(i).gt.0.0)then
c             frnet(i)=frnet(i)-ftnet(i)*(r2(i)-r2(i-1))/(alpha*(
c     1        p(i-1)-p(i))+t2(i)-t2(i-1)) 
c            else
c              frnet(i)=frnet(i)-ftnet(i)*(r2(i+1)-r2(i))/(alpha*(
c     1        p(i)-p(i+1))+t2(i+1)-t2(i))
c            end if  
c            ftnet(i)=0.0  
c          end if
 
          T3(I)=T1(I)+2.*DT*FTNET(I)
          R3(I)=R1(I)+2.*DT*FRNET(I)
	    R3(I)=MAX(R3(I),3.0E-6)
          U3(I)=U1(I)+2.*DT*FUNET(I)
          V3(I)=V1(I)+2.*DT*FVNET(I)
          TRA3(I,1)=TRA1(I,1)+2.*DT*FTRANET(I)
          T1(I)=T2(I)+AFILT*(T1(I)+T3(I)-2.*T2(I))
          R1(I)=R2(I)+AFILT*(R1(I)+R3(I)-2.*R2(I))
          U1(I)=U2(I)+AFILT*(U1(I)+U3(I)-2.*U2(I))
          V1(I)=V2(I)+AFILT*(V1(I)+V3(I)-2.*V2(I))
          TRA1(I,1)=TRA2(I,1)+AFILT*(TRA1(I,1)+TRA3(I,1)-2.*TRA2(I,1))
          T2(I)=T3(I)
          R2(I)=R3(I)
          U2(I)=U3(I)
          V2(I)=V3(I)
          TRA2(I,1)=TRA3(I,1)
         END IF
  160	CONTINUE
C
C  Return to beginning of time loop
C
        GOTO 100
C
C  Output processing starts here
C
  200	CONTINUE
C
C  Output to file sounding.out
C
      OPEN(UNIT=10,FILE='output/sounding.out',
     1   STATUS='UNKNOWN')
C
      TS2C=TS2(1)-273.15
      WRITE(10,210) NP,1000.*CBMF,TS2C
  210	FORMAT(2X, 'N levels = ',T15,I3,
     1    T28,'CBMF=',T35,F8.3,T48,'SST=',T55,F5.2,' C',/)
C
      WRITE(10,220)
  220	FORMAT(2X,'Pressure',T14,'Temp.',T25,'Spec. humid.',T40,
     1    'O3 mr',T48, 'Omega',T56,'Rad. cooling')
      WRITE(10,230)
  230	FORMAT(4X,'(mb)',T15,'(C)',T28,'(g/kg)',T40,'10**-6',T47,
     1    '(mb/hr)',T57,'(deg/day)')
  	WRITE(10,240)
  240	FORMAT(2X,'-------',T14,'-----',T25,'------------',T40,
     1   '-----',T48, '-----',T56,'------------')
C
        DO 260 I=1,NP
         OZ=1.0E6*O3(I)
         TC=T2(I)-273.15
         RM=1000.0*R2(I)
c	   RM=1000.*RHTEMP(I)*RS(I)
         AOMEG=OMEGA(I)*3600.0
         ARADC=-FTBAR(I,2)*3600.*24./FLOAT(NAV)
         WRITE(10,250)P(I),TC,RM,OZ,AOMEG,ARADC    
  250	 FORMAT(3X,F6.1,T12,F8.3,T25,F10.6,T39,F6.3,T47,F6.2,
     1   T58,F7.3)
  260	CONTINUE	
C
      OPEN(UNIT=12,FILE='output/time.out',
     1   STATUS='UNKNOWN')
      OPEN(UNIT=14,FILE='output/profile.out',
     1   STATUS='UNKNOWN')
c	OPEN(UNIT=15,FILE='output/bhov.out',
c     1   STATUS='UNKNOWN')
c	OPEN(UNIT=16,FILE='output/rhov.out',
c     1   STATUS='UNKNOWN')
	OPEN(UNIT=17,FILE='output/plabel.out',
     1   STATUS='UNKNOWN')
C
C  Normalize averaged quantities
C
        FAC=1./FLOAT(NAV)
        BUOYBAR(1,2)=0.0
        DO 300 I=1,NP
         RHBAR(I)=FAC*RHBAR(I)
	   CFRACBAR(I)=FAC*CFRACBAR(I)
	   CLQBAR(I)=FAC*1000.*CLQBAR(I)
         TRAGRAPH(I)=TRA2(I,1)
         DO 301 J=1,NPR    
C         BTIME(J,I)=BTIME(J,I)
	    WTIME(J,I)=1000.*(WTIME(J,I)-WTIME(NPR,I))
  301	 CONTINUE
         DO 270 J=1,2
          BUOYBAR(I,J)=FAC*BUOYBAR(I,J)
          BUOYBAR(I,J)=MAX(BUOYBAR(I,J),-4.0)
          HBAR(I,J)=FAC*HBAR(I,J)/CPD
          HBAR(I,J)=MIN(HBAR(I,J),HBAR(I,2))
          MDET(I,J)=FAC*MDET(I,J)*1000.0
          WATERBAR(I,J)=FAC*WATERBAR(I,J)*1000.0
  270	 CONTINUE
c	 BUOYBAR(I,2)=0.5*BUOYBAR(I,2)
         BUOYBAR(I,2)=MIN(BUOYBAR(I,2),5.0)
         DO 280 J=1,3
          MBAR(I,J)=FAC*MBAR(I,J)*1000.0
  280	 CONTINUE
         DO 290 J=1,4
          FTBAR(I,J)=FAC*FTBAR(I,J)*24.*3600.
          FRBAR(I,J)=FAC*FRBAR(I,J)*1000.*24.*3600.
  290	 CONTINUE
C
C  Clean up SIJ array for display purposes
C
         DO 295 J=1,NP
          IF(SIJ(I,J).GE.1.0)SIJ(I,J)=0.0
          SIJ(I,J)=MAX(SIJ(I,J),0.0)
          SIJ(I,J)=SIJ(I,J)-0.001
  295	 CONTINUE
  300	CONTINUE
C
C  Headers for output files
C
c     WRITE(12,3)
    3 FORMAT(1X,'Time (days),','Precipitation,','Evaporation,',
     * 'Surface T,','Sea surface temperature,',
     1 '500 mb T,','500 mb q,','Top of atmosphere shortwave,', 
     2 'Top of atmosphere longwave')
C
C      WRITE(14,4)
    4 FORMAT(1X,'Pressure (mb),','Temperature (C),',
     * 'Specific humidity (g/Kg),',
     * 'Buoyancy,','Relative humidity,',
     1 'Convective heating rate,','Radiative heating rate,',
     * 'Turbulent heating rate,','Adiabatic warming rate,',
     * 'Upward Mass Flux,',
     2 'Penetrative Mass Flux,','Unsaturated downdraft M,',
     3 'Entrainment rate,','Detrainment rate,','Net radiative flux,',
     4 'Moist static energy')
C
C  Write to output files
C
      DO 310 I=1,NPR
       WRITE(12,305)GTIME(I),RAIN_EVAP(I,1),RAIN_EVAP(I,2),
     1  SST_TS(I,1),SST_TS(I,2),TQ(I,1),TQ(I,2),SWTOPTIME(I),
     1  LWTOPTIME(I),VPOT(I), PPOT(I), GTIMEM(I), TOUT(I)
  305  FORMAT(1X,12(F9.3,' '),F9.3)
  310 CONTINUE
C
      DO 320 I=1,NP
        TC=T2(I)-273.15
        QG=1000.*R2(I)
        ALV=LV0-CPVMCL*(T2(I)-273.15)
       WRITE(14,315)P(I),TC,QG,BUOYBAR(I,1),RHBAR(I),FTBAR(I,1),
     1  FTBAR(I,2),FTBAR(I,3),FTBAR(I,4),MBAR(I,1),MBAR(I,2),
     2  MBAR(I,3),MDET(I,1),MDET(I,2),hbar(i,1),CFRACBAR(I),CLQBAR(I),
	3  WATERBAR(I,2),WATERBAR(I,1),TP(I)
  315  FORMAT(1X,16(F9.3,' '),E10.2,3(1X,F9.3))
  320 CONTINUE
C
C   Time-height sections
C
	OPEN(UNIT=15,FILE='output/cldhov.out',STATUS='UNKNOWN')
	DO I=1,NPR
	 WRITE(15,322)(CLDHOV(I,J),J=1,NP)
	END DO
  322 FORMAT(50(1X,F6.4))
	CLOSE(15)
C
	OPEN(UNIT=15,FILE='output/thov.out',STATUS='UNKNOWN')
	DO I=1,NPR
	 WRITE(15,323)(THOV(I,J),J=1,NP)
	END DO
	CLOSE(15)
  323 format(50(1x,f6.2))
C
	OPEN(UNIT=15,FILE='output/qhov.out',STATUS='UNKNOWN')
	DO I=1,NPR
	 WRITE(15,323)(QHOV(I,J),J=1,NP)
	END DO
	CLOSE(15)
C
	OPEN(UNIT=15,FILE='output/rhhov.out',STATUS='UNKNOWN')
	DO I=1,NPR
	 WRITE(15,323)(RHHOV(I,J),J=1,NP)
	END DO
	CLOSE(15)
C
	OPEN(UNIT=15,FILE='output/mhov.out',STATUS='UNKNOWN')
	DO I=1,NPR
	 WRITE(15,323)(MHOV(I,J),J=1,NP)
	END DO
	CLOSE(15)
C
	OPEN(UNIT=15,FILE='output/mdhov.out',STATUS='UNKNOWN')
	DO I=1,NPR
	 WRITE(15,323)(MDHOV(I,J),J=1,NP)
	END DO
	CLOSE(15)
C
	OPEN(UNIT=15,FILE='output/mphov.out',STATUS='UNKNOWN')
	DO I=1,NPR
	 WRITE(15,323)(MPHOV(I,J),J=1,NP)
	END DO
	CLOSE(15)
C
	OPEN(UNIT=15,FILE='output/lwbar.out',STATUS='UNKNOWN')
	 WRITE(15,*)(topbar/float(nradc))
	CLOSE(15)
C
	OPEN(UNIT=15,FILE='output/botbar.out',STATUS='UNKNOWN')
	 WRITE(15,*)(botbar/float(nradc))
	CLOSE(15)
C
      OPEN(UNIT=15,FILE='output/buoy.out',STATUS='UNKNOWN')
	DO 330 J=1,NPR
	 WRITE(15,325)(BTIME(J,I),I=1,NP)
c	 WRITE(16,325)(WTIME(J,I),J=1,NPR)
c	 WRITE(17,326)P(I)
  330	CONTINUE
  325	FORMAT(2X,50(F9.3,X))
  326 FORMAT(2X,F9.3)
      CLOSE(15)
C
        STOP 
        END
C
***************************************************************************
*****                       SUBROUTINE CONVECT                        *****
*****                          VERSION 4.3b                           *****
*****                         20 August, 2000                         *****
*****                          Kerry Emanuel                          *****
***************************************************************************
C
        SUBROUTINE CONVECT
     *    (T,   Q,    QS, FS,FL,    U,    V,      TRA,    P,    PH,
     *     ND,  NL,   NTRA,   DELT, IFLAG,  FT,     FQ,   FU,
     *     FV,  FTRA, PRECIP, WD,   TPRIME, QPRIME, QCONDC,CBMF,
     *     ENTPRO,ENTPRO2)
C
C-----------------------------------------------------------------------------
C    *** On input:      ***
C
C     T:   Array of absolute temperature (K) of dimension ND, with first
C           index corresponding to lowest model level. Note that this array
C           will be altered by the subroutine if dry convective adjustment
C           occurs and if IPBL is not equal to 0.
C
C     Q:   Array of specific humidity (gm/gm) of dimension ND, with first
C            index corresponding to lowest model level. Must be defined
C            at same grid levels as T. Note that this array will be altered
C            if dry convective adjustment occurs and if IPBL is not equal to 0.
C
C     QS:  Array of saturation specific humidity of dimension ND, with first
C            index corresponding to lowest model level. Must be defined
C            at same grid levels as T. Note that this array will be altered
C            if dry convective adjustment occurs and if IPBL is not equal to 0.
C
C     U:   Array of zonal wind velocity (m/s) of dimension ND, witth first
C            index corresponding with the lowest model level. Defined at
C            same levels as T. Note that this array will be altered if
C            dry convective adjustment occurs and if IPBL is not equal to 0.
C
C     V:   Same as U but for meridional velocity.
C
C     TRA: Array of passive tracer mixing ratio, of dimensions (ND,NTRA),
C            where NTRA is the number of different tracers. If no
C            convective tracer transport is needed, define a dummy
C            input array of dimension (ND,1). Tracers are defined at
C            same vertical levels as T. Note that this array will be altered
C            if dry convective adjustment occurs and if IPBL is not equal to 0.
C
C     P:   Array of pressure (mb) of dimension ND, with first
C            index corresponding to lowest model level. Must be defined
C            at same grid levels as T.
C
C     PH:  Array of pressure (mb) of dimension ND+1, with first index
C            corresponding to lowest level. These pressures are defined at
C            levels intermediate between those of P, T, Q and QS. The first
C            value of PH should be greater than (i.e. at a lower level than)
C            the first value of the array P.
C
C     ND:  The dimension of the arrays T,Q,QS,P,PH,FT and FQ
C
C     NL:  The maximum number of levels to which convection can
C            penetrate, plus 1.
C            NL MUST be less than or equal to ND-1.
C
C     NTRA:The number of different tracers. If no tracer transport
C            is needed, set this equal to 1. (On most compilers, setting
C            NTRA to 0 will bypass tracer calculation, saving some CPU.)  
C
C     DELT: The model time step (sec) between calls to CONVECT
C
C -- sb: interface with the cloud parameterization:
C
C     QCONDC: mixing ratio of condensed water within clouds (kg/kg) 
C               For use in the Bony-Emanuel cloud parameterization 
C sb --
C
C----------------------------------------------------------------------------
C    ***   On Output:         ***
C
C     IFLAG: An output integer whose value denotes the following:
C
C                VALUE                        INTERPRETATION
C                -----                        --------------
C                  0               No moist convection; atmosphere is not
C                                  unstable, or surface temperature is less
C                                  than 250 K or surface specific humidity
C                                  is non-positive.
C
C                  1               Moist convection occurs.
C
C                  2               No moist convection: lifted condensation
C                                  level is above the 200 mb level.
C
C                  3               No moist convection: cloud base is higher
C                                  then the level NL-1.
C
C                  4               Moist convection occurs, but a CFL condition
C                                  on the subsidence warming is violated. This
C                                  does not cause the scheme to terminate.
C
C     FT:   Array of temperature tendency (K/s) of dimension ND, defined at same
C             grid levels as T, Q, QS and P.
C
C     FQ:   Array of specific humidity tendencies ((gm/gm)/s) of dimension ND,
C             defined at same grid levels as T, Q, QS and P.
C
C     FU:   Array of forcing of zonal velocity (m/s^2) of dimension ND,
C             defined at same grid levels as T.
C
C     FV:   Same as FU, but for forcing of meridional velocity.
C
C     FTRA: Array of forcing of tracer content, in tracer mixing ratio per
C             second, defined at same levels as T. Dimensioned (ND,NTRA).
C
C     PRECIP: Scalar convective precipitation rate (mm/day).
C
C     WD:    A convective downdraft velocity scale. For use in surface
C             flux parameterizations. See convect.ps file for details.
C
C     TPRIME: A convective downdraft temperature perturbation scale (K).
C              For use in surface flux parameterizations. See convect.ps
C              file for details.
C
C     QPRIME: A convective downdraft specific humidity
C              perturbation scale (gm/gm).
C              For use in surface flux parameterizations. See convect.ps
C              file for details.
C
C     CBMF:   The cloud base mass flux ((kg/m**2)/s). THIS SCALAR VALUE MUST
C              BE STORED BY THE CALLING PROGRAM AND RETURNED TO CONVECT AT
C              ITS NEXT CALL. That is, the value of CBMF must be "remembered"
C              by the calling program between calls to CONVECT.
C
C------------------------------------------------------------------------------
C
C    ***  THE PARAMETER NA SHOULD IN GENERAL BE GREATER THAN   ***
C    ***                OR EQUAL TO  ND + 1                    ***
C
      PARAMETER (NA=120)           
C
      INTEGER NENT(NA)
      REAL T(ND),Q(ND),QS(ND),U(ND),V(ND),TRA(ND,NTRA),P(ND),PH(ND)
      REAL FT(ND),FQ(ND),FU(ND),FV(ND),FTRA(ND,NTRA)
      REAL UENT(NA,NA),VENT(NA,NA),TRAENT(NA,NA,NTRA),TRATM(NA)
      REAL UP(NA),VP(NA),TRAP(NA,NTRA)
      REAL M(NA),MP(NA),MENT(NA,NA),QENT(NA,NA),ELIJ(NA,NA)
      REAL SIJ(NA,NA),TVP(NA),TV(NA),WATER(NA)
      REAL QP(NA),EP(NA),TH(NA),WT(NA),EVAP(NA),CLW(NA)
      REAL SIGP(NA),TP(NA),TOLD(NA),CPN(NA)
      REAL LV(NA),LVCP(NA),LV0,H(NA),HP(NA),GZ(NA),HM(NA)
C -- sb:
      REAL QCONDC(ND)
      REAL QCOND(NA),NQCOND(NA),WA(NA),MA(NA),SIGA(NA),AX(NA)
C
      REAL MUP(NA),MDOWN(NA),ENT(NA),DET(NA),FQDET(NA),FQSUB(NA)
      COMMON / CVT / TVP,MP,QP,SIJ,MUP,MDOWN,ENT,DET,FQDET,FQSUB,
     1     WATER, INB, TP,WT
C
C -----------------------------------------------------------------------
C
C   ***                     Specify Switches                         ***
C
C   ***   IPBL: Set to zero to bypass dry adiabatic adjustment       ***
C   ***    Any other value results in dry adiabatic adjustment       ***
C   ***     (Zero value recommended for use in models with           ***
C   ***                   boundary layer schemes)                    ***
C
C   ***   MINORIG: Lowest level from which convection may originate  ***
C   ***     (Should be first model level at which T is defined       ***
C   ***      for models using bulk PBL schemes; otherwise, it should ***
C   ***      be the first model level at which T is defined above    ***
C   ***                      the surface layer)                      ***
C
        IPBL=1
        MINORIG=1
C
C------------------------------------------------------------------------------
C
C   ***                    SPECIFY PARAMETERS                        ***
C
C   *** ELCRIT IS THE AUTOCONVERSION THERSHOLD WATER CONTENT (gm/gm) ***
C   ***  TLCRIT IS CRITICAL TEMPERATURE BELOW WHICH THE AUTO-        ***
C   ***       CONVERSION THRESHOLD IS ASSUMED TO BE ZERO             ***
C   ***     (THE AUTOCONVERSION THRESHOLD VARIES LINEARLY            ***
C   ***               BETWEEN 0 C AND TLCRIT)                        ***
C   ***   ENTP IS THE COEFFICIENT OF MIXING IN THE ENTRAINMENT       ***
C   ***                       FORMULATION                            ***
C   ***  SIGD IS THE FRACTIONAL AREA COVERED BY UNSATURATED DNDRAFT  ***
C   ***  SIGS IS THE FRACTION OF PRECIPITATION FALLING OUTSIDE       ***
C   ***                        OF CLOUD                              ***
C   ***        OMTRAIN IS THE ASSUMED FALL SPEED (P/s) OF RAIN       ***
C   ***     OMTSNOW IS THE ASSUMED FALL SPEED (P/s) OF SNOW          ***
C   ***  COEFFR IS A COEFFICIENT GOVERNING THE RATE OF EVAPORATION   ***
C   ***                          OF RAIN                             ***
C   ***  COEFFS IS A COEFFICIENT GOVERNING THE RATE OF EVAPORATION   ***
C   ***                          OF SNOW                             ***
C   ***     CU IS THE COEFFICIENT GOVERNING CONVECTIVE MOMENTUM      ***
C   ***                         TRANSPORT                            ***
C   ***    DTMAX IS THE MAXIMUM NEGATIVE TEMPERATURE PERTURBATION    ***
C   ***        A LIFTED PARCEL IS ALLOWED TO HAVE BELOW ITS LFC      ***
C   ***    ALPHA AND DAMP ARE PARAMETERS THAT CONTROL THE RATE OF    ***
C   ***                 APPROACH TO QUASI-EQUILIBRIUM                ***
C   ***   (THEIR STANDARD VALUES ARE  0.20 AND 0.1, RESPECTIVELY)    ***
C   ***                   (DAMP MUST BE LESS THAN 1)                 ***
C
        ELCRIT=.0011
        TLCRIT=-55.0
        ENTP=1.5
        SIGD=0.05
        SIGS=0.12
	  SIGMIN=0.0
	  SIGMAX=0.9 
        OMTRAIN=50.0
        OMTSNOW=5.5 
        COEFFR=0.9
        COEFFS=0.6
        CU=0.7
        BETA=10.0
        DTMAX=0.9 
        ALPHA=0.015  
        DAMP=0.05
        DELTA=0.01       ! sb (for cloud parameterization)
C
c	DTMAX=0.1*DTMAX1*(T(1)-273.15)
c	DTMAX=MIN(DTMAX,DTMAX1)
c	DTMAX=MAX(DTMAX,0.0)
C
C   ***        ASSIGN VALUES OF THERMODYNAMIC CONSTANTS,        ***
C   ***            GRAVITY, AND LIQUID WATER DENSITY.           ***
C   ***             THESE SHOULD BE CONSISTENT WITH             ***
C   ***              THOSE USED IN CALLING PROGRAM              ***
C   ***     NOTE: THESE ARE ALSO SPECIFIED IN SUBROUTINE TLIFT  ***
C
      CPD=1005.7
      CPV=1870.0
c      CL=4190.0 
      CL=2500.0 
      RV=461.5
      RD=287.04
      LV0=2.501E6
      G=9.8  
      ROWL=1000.0
C
      CPVMCL=CL-CPV 
      EPS=RD/RV
      EPSI=1./EPS
      GINV=1.0/G
      DELTI=1.0/DELT
C
C           ***  INITIALIZE OUTPUT ARRAYS AND PARAMETERS  ***
C
        DO 5 I=1,ND
         FT(I)=0.0
         FQ(I)=0.0
         FU(I)=0.0
         FV(I)=0.0
C -- sb:
         QCONDC(I)=0.0
         QCOND(I)=0.0
         NQCOND(I)=0.0
         DO 4 J=1,NTRA
          FTRA(I,J)=0.0
    4    CONTINUE
         MDOWN(I)=0.0
         MUP(I)=0.0
         ENT(I)=0.0
         DET(I)=0.0
         FQDET(I)=0.0
         FQSUB(I)=0.0
         WATER(I)=0.0
         QP(I)=Q(I)
         MP(I)=0.0
         TVP(I)=T(I)*(1.+Q(I)*EPSI-Q(I))
    5   CONTINUE
        DO 7 I=1,NL+1
         RDCP=(RD*(1.-Q(I))+Q(I)*RV)/
     1    (CPD*(1.-Q(I))+Q(I)*CPV)
         TH(I)=T(I)*(1000.0/P(I))**RDCP
    7   CONTINUE
        PRECIP=0.0
        WD=0.0
        TPRIME=0.0
        QPRIME=0.0
        ENTPRO=0.0
        ENTPRO2=0.0
        IFLAG=0
        INB=ND
C
        IF(IPBL.NE.0)THEN
C
C     ***            PERFORM DRY ADIABATIC ADJUSTMENT            ***
C
        JC=0
        DO 30 I=NL-1,1,-1
         JN=0
          SUM=TH(I)*(1.+Q(I)*EPSI-Q(I))
         DO 10 J=I+1,NL
          SUM=SUM+TH(J)*(1.+Q(J)*EPSI-Q(J))
          THBAR=SUM/FLOAT(J+1-I)
          IF((TH(J)*(1.+Q(J)*EPSI-Q(J))).LT.THBAR)JN=J
          IF(I.EQ.1.AND.P(J).GT.948.)JN=MAX(JN,J)
   10    CONTINUE
c         IF(I.EQ.1)JN=MAX(JN,3)
         IF(JN.EQ.0)GOTO 30
   12    CONTINUE
         AHM=0.0
         RM=0.0
         UM=0.0
         VM=0.0
         DO K=1,NTRA
          TRATM(K)=0.0
         END DO
         DO 15 J=I,JN
          AHM=AHM+(CPD*(1.-Q(J))+Q(J)*CPV)*T(J)*(PH(J)-PH(J+1))
          RM=RM+Q(J)*(PH(J)-PH(J+1))
          UM=UM+U(J)*(PH(J)-PH(J+1))
          VM=VM+V(J)*(PH(J)-PH(J+1))
          DO K=1,NTRA
           TRATM(K)=TRATM(K)+TRA(J,K)*(PH(J)-PH(J+1))
          END DO
   15    CONTINUE
         DPHINV=1./(PH(I)-PH(JN+1))
         RM=RM*DPHINV
         UM=UM*DPHINV
         VM=VM*DPHINV
         DO K=1,NTRA
          TRATM(K)=TRATM(K)*DPHINV
         END DO
         A2=0.0
         DO 20 J=I,JN
          Q(J)=RM
          U(J)=UM
          V(J)=VM
          DO K=1,NTRA
           TRA(J,K)=TRATM(K)
          END DO
          RDCP=(RD*(1.-Q(J))+Q(J)*RV)/
     1     (CPD*(1.-Q(J))+Q(J)*CPV)  
          X=(0.001*P(J))**RDCP
          TOLD(J)=T(J)
          T(J)=X
          A2=A2+(CPD*(1.-Q(J))+Q(J)*CPV)*X*(PH(J)-PH(J+1))
   20    CONTINUE
         DO 25 J=I,JN
          TH(J)=AHM/A2
          T(J)=T(J)*TH(J)
          TC=TOLD(J)-273.15
          ALV=LV0-CPVMCL*TC
          QS(J)=QS(J)+QS(J)*(1.+QS(J)*0.608)*ALV*(T(J)-
     1     TOLD(J))/(RV*TOLD(J)*TOLD(J))
   25    CONTINUE
         IF((TH(JN+1)*(1.+Q(JN+1)*EPSI-Q(JN+1))).LT.
     1    (TH(JN)*(1.+Q(JN)*EPSI-Q(JN))))THEN
          JN=JN+1
          GOTO 12
         END IF
         IF(I.EQ.1)JC=JN 
   30   CONTINUE
C
C   ***   Remove any supersaturation that results from adjustment ***
C
      IF(JC.GT.1)THEN
       DO 38 J=1,JC
          IF(QS(J).LT.Q(J))THEN 
           ALV=LV0-CPVMCL*(T(J)-273.15)  
           TNEW=T(J)+ALV*(Q(J)-QS(J))/(CPD*(1.-Q(J))+
     1      CL*Q(J)+QS(J)*(CPV-CL+ALV*ALV/(RV*T(J)*T(J))))
           ALVNEW=LV0-CPVMCL*(TNEW-273.15)
           QNEW=(ALV*Q(J)-(TNEW-T(J))*(CPD*(1.-Q(J))+CL*Q(J)))/ALVNEW
           PRECIP=PRECIP+24.*3600.*1.0E5*(PH(J)-PH(J+1))*
     1      (Q(J)-QNEW)/(G*DELT*ROWL)
           T(J)=TNEW
           Q(J)=QNEW
           QS(J)=QNEW
          END IF     
   38  CONTINUE  
      END IF
C
      END IF
C
C  *** CALCULATE ARRAYS OF GEOPOTENTIAL, HEAT CAPACITY AND STATIC ENERGY
C  
        GZ(1)=0.0
        CPN(1)=CPD*(1.-Q(1))+Q(1)*CPV
        H(1)=T(1)*CPN(1)
        LV(1)=LV0-CPVMCL*(T(1)-273.15)
        HM(1)=LV(1)*Q(1)
        TV(1)=T(1)*(1.+Q(1)*EPSI-Q(1))
        AHMIN=1.0E12
        IHMIN=NL
        DO 40 I=2,NL+1
          TVX=T(I)*(1.+Q(I)*EPSI-Q(I))
          TVY=T(I-1)*(1.+Q(I-1)*EPSI-Q(I-1))
          GZ(I)=GZ(I-1)+0.5*RD*(TVX+TVY)*(P(I-1)-P(I))/PH(I)
          CPN(I)=CPD*(1.-Q(I))+CPV*Q(I)
          H(I)=T(I)*CPN(I)+GZ(I)
          LV(I)=LV0-CPVMCL*(T(I)-273.15)
          HM(I)=(CPD*(1.-Q(I))+CL*Q(I))*(T(I)-T(1))+
     1     LV(I)*Q(I)+GZ(I)
          TV(I)=T(I)*(1.+Q(I)*EPSI-Q(I))
C
C  ***  Find level of minimum moist static energy    ***
C
          IF(I.GE.MINORIG.AND.HM(I).LT.AHMIN.AND.HM(I).LT.HM(I-1))THEN
           AHMIN=HM(I)
           IHMIN=I
          END IF
   40   CONTINUE
        IHMIN=MIN(IHMIN, NL-1)
C
C  ***     Find that model level below the level of minimum moist       ***
C  ***  static energy that has the maximum value of moist static energy ***
C
        AHMAX=0.0
        DO 42 I=MINORIG,IHMIN
         IF(HM(I).GT.AHMAX)THEN
          NK=I
          AHMAX=HM(I)
         END IF
   42   CONTINUE
C
C  ***  CHECK WHETHER PARCEL LEVEL TEMPERATURE AND SPECIFIC HUMIDITY   ***
C  ***                          ARE REASONABLE                         ***
C  ***      Skip convection if HM increases monotonically upward       ***
C
        IF(T(NK).LT.250.0.OR.Q(NK).LE.0.0.OR.IHMIN.EQ.(NL-1))THEN
         IFLAG=0
         CBMF=0.0
         RETURN
        END IF
C
C   ***  CALCULATE LIFTED CONDENSATION LEVEL OF AIR AT PARCEL ORIGIN LEVEL ***
C   ***       (WITHIN 0.2% OF FORMULA OF BOLTON, MON. WEA. REV.,1980)      ***
C
        RH=Q(NK)/QS(NK)
        CHI=T(NK)/(1669.0-122.0*RH-T(NK))
        PLCL=P(NK)*(RH**CHI)
        IF(PLCL.LT.200.0.OR.PLCL.GE.2000.0)THEN
         IFLAG=2
         CBMF=0.0
         RETURN
        END IF
c	  PBLH=100.0
c	  IF(JC.GT.1)THEN
	   PBLH=(RD*T(NK)/G)*LOG(P(1)/PLCL)
c	  END IF
	  ROWS=100.*P(1)/(RD*T(1))
	  W3=PBLH*G*(FS/(CPD*T(1))+0.608*FL/2.5E6)/ROWS
	  W3=MAX(W3,0.0)
	  WSTAR=W3**(1./3.)
	  DTMAX=1.737*WSTAR
	  dtmax=max(dtmax,0.9)
c	  dtmax=2.4
C
C   ***  CALCULATE FIRST LEVEL ABOVE LCL (=ICB)  ***
C
        ICB=NL-1
        DO 50 I=NK+1,NL
         IF(P(I).LT.PLCL)THEN
          ICB=MIN(ICB,I)
         END IF
   50   CONTINUE
        IF(ICB.GE.(NL-1))THEN
         IFLAG=3
         CBMF=0.0
         RETURN
        END IF
C
C   *** FIND TEMPERATURE UP THROUGH ICB AND TEST FOR INSTABILITY           ***
C
C   *** SUBROUTINE TLIFT CALCULATES PART OF THE LIFTED PARCEL VIRTUAL      ***
C   ***  TEMPERATURE, THE ACTUAL TEMPERATURE AND THE ADIABATIC             ***
C   ***                   LIQUID WATER CONTENT                             ***
C
        CALL TLIFT(P,T,Q,QS,GZ,ICB,NK,TVP,TP,CLW,ND,NL,1)
        DO 54 I=NK,ICB
         TVP(I)=TVP(I)-TP(I)*Q(NK)
   54   CONTINUE
C
C   ***  If there was no convection at last time step and parcel    ***
C   ***       is stable at ICB then skip rest of calculation        ***
C
        IF(CBMF.EQ.0.0.AND.TVP(ICB).LE.(TV(ICB)-DTMAX))THEN
         IFLAG=0
         RETURN
        END IF
C
C   ***  IF THIS POINT IS REACHED, MOIST CONVECTIVE ADJUSTMENT IS NECESSARY ***
C
        IF(IFLAG.NE.4)IFLAG=1
C
C   ***  FIND THE REST OF THE LIFTED PARCEL TEMPERATURES          ***
C
        CALL TLIFT(P,T,Q,QS,GZ,ICB,NK,TVP,TP,CLW,ND,NL,2)
C
C   ***  SET THE PRECIPITATION EFFICIENCIES AND THE FRACTION OF   ***
C   ***          PRECIPITATION FALLING OUTSIDE OF CLOUD           ***
C   ***      THESE MAY BE FUNCTIONS OF TP(I), P(I) AND CLW(I)     ***
C                 
        DO 57 I=1,NK
         EP(I)=0.0
         SIGP(I)=SIGS
   57   CONTINUE
        DO 60 I=NK+1,NL
         TCA=TP(I)-273.15
         IF(TCA.GE.0.0)THEN
          ELACRIT=ELCRIT
         ELSE
          ELACRIT=ELCRIT*(1.0-TCA/TLCRIT)
         END IF
         ELACRIT=MAX(ELACRIT,0.0)
           EPMAX=0.999
         EP(I)=EPMAX*(1.0-ELACRIT/MAX(CLW(I),1.0E-8))
         EP(I)=MAX(EP(I),0.0)
         EP(I)=MIN(EP(I),EPMAX)
         SIGP(I)=SIGS
   60   CONTINUE
C
C   ***       CALCULATE VIRTUAL TEMPERATURE AND LIFTED PARCEL     ***
C   ***                    VIRTUAL TEMPERATURE                    ***
C
        DO 64 I=ICB+1,NL
c        TVP(I)=TVP(I)*(1.-Q(NK)+EP(I)*CLW(I))
        TVP(I)=TVP(I)-Q(NK)*TP(I)
   64   CONTINUE
        TVP(NL+1)=TVP(NL)-(GZ(NL+1)-GZ(NL))/CPD
C
C   ***        NOW INITIALIZE VARIOUS ARRAYS USED IN THE COMPUTATIONS       ***
C
        DO 70 I=1,NL+1
         HP(I)=H(I)
         NENT(I)=0
         WATER(I)=0.0
         EVAP(I)=0.0
         WT(I)=OMTSNOW
         MP(I)=0.0
         M(I)=0.0
         LVCP(I)=LV(I)/CPN(I)
         DO 70 J=1,NL+1
          QENT(I,J)=Q(J)
          ELIJ(I,J)=0.0
          MENT(I,J)=0.0
          SIJ(I,J)=0.0
          UENT(I,J)=U(J)
          VENT(I,J)=V(J)
          DO 70 K=1,NTRA
           TRAENT(I,J,K)=TRA(J,K)
   70   CONTINUE
        QP(1)=Q(1)
        UP(1)=U(1)
        VP(1)=V(1)
        DO 71 I=1,NTRA
         TRAP(1,I)=TRA(1,I)
   71	CONTINUE
        DO 72 I=2,NL+1
         QP(I)=Q(I-1)
         UP(I)=U(I-1)
         VP(I)=V(I-1)
         DO 72 J=1,NTRA
          TRAP(I,J)=TRA(I-1,J)
   72	CONTINUE
C
C  ***  FIND THE FIRST MODEL LEVEL (INB1) ABOVE THE PARCEL'S      ***
C  ***          HIGHEST LEVEL OF NEUTRAL BUOYANCY                 ***
C  ***     AND THE HIGHEST LEVEL OF POSITIVE CAPE (INB)           ***
C
        CAPE=0.0
        CAPEM=0.0
        INB=ICB+1
        INB1=INB
	  BYP=0.0
        DO 82 I=ICB+1,NL-1
         BY=(TVP(I)-TV(I))*(PH(I)-PH(I+1))/P(I)
         CAPE=CAPE+BY
         IF(BY.GE.0.0)INB1=I+1
         IF(CAPE.GT.0.0)THEN
          INB=I+1
          CAPEM=CAPE
          BYP=(TVP(I+1)-TV(I+1))*(PH(I+1)-PH(I+2))/P(I+1)
         END IF
   82   CONTINUE
	  INB=MAX(INB,INB1)
        CAPE=CAPEM+BYP
        DEFRAC=CAPEM-CAPE
        DEFRAC=MAX(DEFRAC,0.001)
        FRAC=-CAPE/DEFRAC
        FRAC=MIN(FRAC,1.0)
        FRAC=MAX(FRAC,0.0)
	  JMAX=INB
C
C   ***   CALCULATE LIQUID WATER STATIC ENERGY OF LIFTED PARCEL   ***
C
        DO 95 I=ICB,INB
         HP(I)=H(NK)+(LV(I)+(CPD-CPV)*T(I))*EP(I)*CLW(I)
   95   CONTINUE                  
C
C   ***  CALCULATE CLOUD BASE MASS FLUX AND RATES OF MIXING, M(I),  ***
c   ***                   AT EACH MODEL LEVEL                       ***
C
        DBOSUM=0.0
C   
C   ***     INTERPOLATE DIFFERENCE BETWEEN LIFTED PARCEL AND      ***
C   ***  ENVIRONMENTAL TEMPERATURES TO LIFTED CONDENSATION LEVEL  ***
C	
        TVPPLCL=TVP(ICB-1)-RD*TVP(ICB-1)*(P(ICB-1)-PLCL)/
     1    (CPN(ICB-1)*P(ICB-1))
        TVAPLCL=TV(ICB)+(TVP(ICB)-TVP(ICB+1))*(PLCL-P(ICB))/
     1    (P(ICB)-P(ICB+1))
        DTPBL=0.0
        DO 96 I=NK,ICB-1
         DTPBL=DTPBL+(TVP(I)-TV(I))*(PH(I)-PH(I+1))
   96   CONTINUE
        DTPBL=DTPBL/(PH(NK)-PH(ICB))
        DTMIN=TVPPLCL-TVAPLCL+DTMAX+DTPBL
        DTMA=DTMIN
C
C   ***  ADJUST CLOUD BASE MASS FLUX   ***
C
      CBMFOLD=CBMF
      DELT0=300.0
      DAMPS=DAMP*DELT/DELT0 
      CBMF=(1.-DAMPS)*CBMF+0.1*ALPHA*DTMA 
c      CBMF=0.02*DTMA
      CBMF=MAX(CBMF,0.0)
C
C   *** If cloud base mass flux is zero, skip rest of calculation  ***
C
      IF(CBMF.EQ.0.0.AND.CBMFOLD.EQ.0.0)THEN
       RETURN
      END IF
C
C   ***   CALCULATE RATES OF MIXING,  M(I)   ***
C
      M(ICB)=0.0
      DO 103 I=ICB+1,INB
       K=MIN(I,INB1)
c       DBO=ABS(TV(K+1)-TVP(K+1)-TV(K-1)+TVP(K-1))+
c     1  ENTP*0.04*(PH(K)-PH(K+1))
        dbo = abs(tv(k)-tvp(k))
     &        + entp*0.02*(ph(k)-ph(k+1))
       DBOSUM=DBOSUM+DBO
       M(I)=CBMF*DBO
  103 CONTINUE
      DO 110 I=ICB+1,INB
       M(I)=M(I)/DBOSUM  
  110 CONTINUE     
C
C   ***  CALCULATE ENTRAINED AIR MASS FLUX (MENT), TOTAL WATER MIXING  ***
C   ***     RATIO (QENT), TOTAL CONDENSED WATER (ELIJ), AND MIXING     ***
C   ***                        FRACTION (SIJ)                          ***
C
        DO 170 I=ICB+1,INB
         QTI=Q(NK)-EP(I)*CLW(I)
         DO 160 J=ICB,JMAX
          BF2=1.+LV(J)*LV(J)*QS(J)/(RV*T(J)*T(J)*CPD)
          ANUM=H(J)-HP(I)+(CPV-CPD)*T(J)*(QTI-Q(J))
c          ANUM=H(J)-HP(I)+(CPV-CPD)*T(J)*(QTI-Q(J))-CPN(J)*Q(NK)*TP(J)
          DENOM=H(I)-HP(I)+(CPD-CPV)*(Q(I)-QTI)*T(J)
          DEI=DENOM
          IF(ABS(DEI).LT.0.01)DEI=0.01
          SIJ(I,J)=ANUM/DEI
          SIJ(I,I)=1.0
          ALTEM=SIJ(I,J)*Q(I)+(1.-SIJ(I,J))*QTI-QS(J)
          ALTEM=ALTEM/BF2
          CWAT=CLW(J)*(1.-EP(J))
          STEMP=SIJ(I,J)
          IF((STEMP.LT.0.0.OR.STEMP.GT.1.0.OR.
     1      ALTEM.GT.CWAT).AND.J.GT.I)THEN
           ANUM=ANUM-LV(J)*(QTI-QS(J)-CWAT*BF2)
           DENOM=DENOM+LV(J)*(Q(I)-QTI)
           IF(ABS(DENOM).LT.0.01)DENOM=0.01
           SIJ(I,J)=ANUM/DENOM
           ALTEM=SIJ(I,J)*Q(I)+(1.-SIJ(I,J))*QTI-QS(J)
           ALTEM=ALTEM-(BF2-1.)*CWAT
          END IF
          IF(SIJ(I,J).GT.SIGMIN.AND.SIJ(I,J).LT.SIGMAX)THEN
           QENT(I,J)=SIJ(I,J)*Q(I)+(1.-SIJ(I,J))*QTI
           UENT(I,J)=SIJ(I,J)*U(I)+(1.-SIJ(I,J))*U(NK)
           VENT(I,J)=SIJ(I,J)*V(I)+(1.-SIJ(I,J))*V(NK)
           DO K=1,NTRA
            TRAENT(I,J,K)=SIJ(I,J)*TRA(I,K)+(1.-SIJ(I,J))*
     1       TRA(NK,K)
           END DO
           ELIJ(I,J)=ALTEM
           ELIJ(I,J)=MAX(0.0,ELIJ(I,J))
           MENT(I,J)=M(I)/(1.-SIJ(I,J))
           NENT(I)=NENT(I)+1
          END IF
          SIJ(I,J)=MAX(0.0,SIJ(I,J))
          SIJ(I,J)=MIN(1.0,SIJ(I,J))
  160    CONTINUE
C
C   ***   IF NO AIR CAN ENTRAIN AT LEVEL I ASSUME THAT UPDRAFT DETRAINS  ***
C   ***   AT THAT LEVEL AND CALCULATE DETRAINED AIR FLUX AND PROPERTIES  ***
C
         IF(NENT(I).EQ.0)THEN
          MENT(I,I)=M(I)
          QENT(I,I)=Q(NK)-EP(I)*CLW(I)
          UENT(I,I)=U(NK)
          VENT(I,I)=V(NK)
          DO J=1,NTRA
           TRAENT(I,I,J)=TRA(NK,J)
          END DO
          ELIJ(I,I)=CLW(I)
          SIJ(I,I)=1.0
         END IF 
  170   CONTINUE
        SIJ(INB,INB)=1.0
C
C   ***  NORMALIZE ENTRAINED AIR MASS FLUXES TO REPRESENT EQUAL  ***
C   ***              PROBABILITIES OF MIXING                     ***
C
        DO 200 I=ICB+1,INB
        IF(NENT(I).NE.0)THEN
         QP1=Q(NK)-EP(I)*CLW(I)
         ANUM=H(I)-HP(I)-LV(I)*(QP1-QS(I))
         DENOM=H(I)-HP(I)+LV(I)*(Q(I)-QP1)
         IF(ABS(DENOM).LT.0.01)DENOM=0.01
         SCRIT=ANUM/DENOM
         ALT=QP1-QS(I)+SCRIT*(Q(I)-QP1)
         IF(ALT.LT.0.0)SCRIT=1.0
	   SCRIT=MAX(SCRIT,0.0)
         ASIJ=0.0
         SMIN=1.0
         DO 175 J=ICB,JMAX
          IF(SIJ(I,J).GT.SIGMIN.AND.SIJ(I,J).LT.SIGMAX)THEN
           IF(J.GT.I)THEN
            SMID=MIN(SIJ(I,J),SCRIT)
            SJMAX=SMID
            SJMIN=SMID
            IF(SMID.LT.SMIN.AND.SIJ(I,J+1).LT.SMID)THEN
             SMIN=SMID
             SJMAX=MIN(SIJ(I,J+1),SIJ(I,J),SCRIT)
             SJMIN=MAX(SIJ(I,J-1),SIJ(I,J))
             SJMIN=MIN(SJMIN,SCRIT)
            END IF
           ELSE
            SJMAX=MAX(SIJ(I,J+1),SCRIT)
            SMID=MAX(SIJ(I,J),SCRIT)
            SJMIN=0.0
            IF(J.GT.1)SJMIN=SIJ(I,J-1)
            SJMIN=MAX(SJMIN,SCRIT)
           END IF
           DELP=ABS(SJMAX-SMID)
           DELM=ABS(SJMIN-SMID)
           ASIJ=ASIJ+(DELP+DELM)*(PH(J)-PH(J+1))
           MENT(I,J)=MENT(I,J)*(DELP+DELM)*(PH(J)-PH(J+1))
          END IF
  175    CONTINUE
         ASIJ=MAX(1.0E-21,ASIJ)
         ASIJ=1.0/ASIJ
         DO 180 J=ICB,INB
          MENT(I,J)=MENT(I,J)*ASIJ
  180    CONTINUE
         BSUM=0.0
         DO 190 J=ICB,INB
          BSUM=BSUM+MENT(I,J)
  190    CONTINUE
         IF(BSUM.LT.1.0E-18)THEN
          NENT(I)=0
          MENT(I,I)=M(I)
          QENT(I,I)=Q(NK)-EP(I)*CLW(I)
          UENT(I,I)=U(NK)
          VENT(I,I)=V(NK)
          DO J=1,NTRA
           TRAENT(I,I,J)=TRA(NK,J)
          END DO
          ELIJ(I,I)=CLW(I)
          SIJ(I,I)=1.0
         END IF
        END IF
  200   CONTINUE
C
C   ***  CHECK WHETHER EP(INB)=0, IF SO, SKIP PRECIPITATING    ***
C   ***             DOWNDRAFT CALCULATION                      ***
C
        IF(EP(INB).LT.0.0001)GOTO 405
C
C   ***  INTEGRATE LIQUID WATER EQUATION TO FIND CONDENSED WATER   ***
C   ***                AND CONDENSED WATER FLUX                    ***
C
        JTT=2
C
C    ***                    BEGIN DOWNDRAFT LOOP                    ***
C
        DO 400 I=INB,1,-1
C
C    ***              CALCULATE DETRAINED PRECIPITATION             ***
C
        WDTRAIN=G*EP(I)*M(I)*CLW(I)
        IF(I.GT.1)THEN
         DO 320 J=1,I-1
         AWAT=ELIJ(J,I)-(1.-EP(I))*CLW(I)
         AWAT=MAX(0.0,AWAT)
  320    WDTRAIN=WDTRAIN+G*AWAT*MENT(J,I)
        END IF
C
C    ***    FIND RAIN WATER AND EVAPORATION USING PROVISIONAL   ***
C    ***              ESTIMATES OF QP(I)AND QP(I-1)             ***
C     
c
c  ***  Value of terminal velocity and coeffecient of evaporation for snow   ***
c 
        COEFF=COEFFS
        WT(I)=OMTSNOW
c      
c  ***  Value of terminal velocity and coeffecient of evaporation for rain   ***
c
        IF(T(I).GT.273.0)THEN
         COEFF=COEFFR
         WT(I)=OMTRAIN
        END IF
        QSM=0.5*(Q(I)+QP(I+1))
        AFAC=COEFF*PH(I)*(QS(I)-QSM)/(1.0E4+2.0E3*PH(I)*QS(I))
        AFAC=MAX(AFAC,0.0)
        SIGT=SIGP(I)
        SIGT=MAX(0.0,SIGT)
        SIGT=MIN(1.0,SIGT)
        B6=100.*(PH(I)-PH(I+1))*SIGT*AFAC/WT(I)
        C6=(WATER(I+1)*WT(I+1)+WDTRAIN/SIGD)/WT(I)
        REVAP=0.5*(-B6+SQRT(B6*B6+4.*C6))
        EVAP(I)=SIGT*AFAC*REVAP
        WATER(I)=REVAP*REVAP
C
C    ***  CALCULATE PRECIPITATING DOWNDRAFT MASS FLUX UNDER     ***
C    ***              HYDROSTATIC APPROXIMATION                 ***
C   
        IF(I.EQ.1)GOTO 360
        DHDP=(H(I)-H(I-1))/(P(I-1)-P(I))
        DHDP=MAX(DHDP,10.0)
        MP(I)=100.*GINV*LV(I)*SIGD*EVAP(I)/DHDP
        MP(I)=MAX(MP(I),0.0)
C
C   ***   ADD SMALL AMOUNT OF INERTIA TO DOWNDRAFT              ***
C
        FAC=20.0/(PH(I-1)-PH(I))
        MP(I)=(FAC*MP(I+1)+MP(I))/(1.+FAC)
C   
C    ***      FORCE MP TO DECREASE LINEARLY TO ZERO                 ***
C    ***      BETWEEN ABOUT 950 MB AND THE SURFACE                  ***
C
          IF(P(I).GT.(0.949*P(1)))THEN
           JTT=MAX(JTT,I)
           MP(I)=MP(JTT)*(P(1)-P(I))/(P(1)-P(JTT))
          END IF              
  360   CONTINUE
C
C    ***       FIND MIXING RATIO OF PRECIPITATING DOWNDRAFT     ***
C
        IF(I.EQ.INB)GOTO 400
        IF(I.EQ.1)THEN
         QSTM=QS(1)
        ELSE
         QSTM=QS(I-1)
        END IF
        IF(MP(I).GT.MP(I+1))THEN
          RAT=MP(I+1)/MP(I)
          QP(I)=QP(I+1)*RAT+Q(I)*(1.0-RAT)+100.*GINV*
     1       SIGD*(PH(I)-PH(I+1))*(EVAP(I)/MP(I))
          UP(I)=UP(I+1)*RAT+U(I)*(1.-RAT)
          VP(I)=VP(I+1)*RAT+V(I)*(1.-RAT)
          DO J=1,NTRA
           TRAP(I,J)=TRAP(I+1,J)*RAT+TRAP(I,J)*(1.-RAT)
          END DO
         ELSE
          IF(MP(I+1).GT.0.0)THEN
            QP(I)=(GZ(I+1)-GZ(I)+QP(I+1)*(LV(I+1)+T(I+1)*(
     1        CL-CPD))+CPD*(T(I+1)-T(I)))/(LV(I)+T(I)*(CL-CPD))
            UP(I)=UP(I+1)
            VP(I)=VP(I+1)
            DO J=1,NTRA
             TRAP(I,J)=TRAP(I+1,J)
            END DO
          END IF
        END IF
        QP(I)=MIN(QP(I),QSTM)
        QP(I)=MAX(QP(I),0.0)
  400   CONTINUE
C
C   ***  CALCULATE SURFACE PRECIPITATION IN MM/DAY     ***
C
        PRECIP=PRECIP+WT(1)*SIGD*WATER(1)*3600.*24000./(ROWL*G)
C
  405   CONTINUE
C
C   ***  CALCULATE DOWNDRAFT VELOCITY SCALE AND SURFACE TEMPERATURE AND  ***
c   ***                    WATER VAPOR FLUCTUATIONS                      ***
C
      WD=BETA*ABS(MP(ICB))*0.01*RD*T(ICB)/(SIGD*P(ICB))
      QPRIME=0.5*(QP(1)-Q(1))
      TPRIME=LV0*QPRIME/CPD
C
C   ***  CALCULATE TENDENCIES OF LOWEST LEVEL POTENTIAL TEMPERATURE  ***
C   ***                      AND MIXING RATIO                        ***
C
        DPINV=0.01/(PH(1)-PH(2))
        AM=0.0
        IF(NK.EQ.1)THEN
         DO 410 K=2,INB
  410    AM=AM+M(K)
        END IF
        IF((2.*G*DPINV*AM).GE.DELTI)IFLAG=4
        FT(1)=FT(1)+G*DPINV*AM*(T(2)-T(1)+(GZ(2)-GZ(1))/CPN(1))
        FT(1)=FT(1)-LVCP(1)*SIGD*EVAP(1)
        FT(1)=FT(1)+SIGD*WT(2)*(CL-CPD)*WATER(2)*(T(2)-
     1   T(1))*DPINV/CPN(1)
c	FT(1)=FT(1)+0.01*SIGD*WT(1)*WATER(1)*RD*TV(1)/(PH(1)*CPN(1))
c	ENTPRO=SIGD*WT(1)*WATER(1)*RD*(P(1)-P(2))/(PH(1)*G)
c	ENTPRO=SIGD*WT(2)*WATER(2)*(GZ(2)-GZ(1))/(G*T(1))
        ENTPRO=0.0
        ENTPRO2=AM*(QS(2)-Q(2))*LV(1)*(T(1)-T(2))/T(1)**2
        ENTPRO2=ENTPRO2-AM*RV*(QS(2)-Q(2))*(PH(1)-PH(2))/P(1)
c	ENTPRO=SIGD*WT(2)*WATER(2)*(GZ(2)-GZ(1))/(G*T(1))-AM*
c     1   (QS(2)-Q(2))*RV*(LOG(Q(1)/QS(1))-LOG(Q(2)/QS(2)))
c	ENTPRO2=0.0
        FQ(1)=FQ(1)+G*MP(2)*(QP(2)-Q(1))*
     1    DPINV+SIGD*EVAP(1)
        FQDET(1)=FQ(1)
        FQ(1)=FQ(1)+G*AM*(Q(2)-Q(1))*DPINV
        FQSUB(1)=G*AM*(Q(2)-Q(1))*DPINV
        FU(1)=FU(1)+G*DPINV*(MP(2)*(UP(2)-U(1))+AM*(U(2)-U(1)))
        FV(1)=FV(1)+G*DPINV*(MP(2)*(VP(2)-V(1))+AM*(V(2)-V(1)))
        DO J=1,NTRA
         FTRA(1,J)=FTRA(1,J)+G*DPINV*(MP(2)*(TRAP(2,J)-TRA(1,J))+
     1    AM*(TRA(2,J)-TRA(1,J)))
        END DO
        AMDE=0.0
        DO 415 J=2,INB
         FQ(1)=FQ(1)+G*DPINV*MENT(J,1)*(QENT(J,1)-Q(1))
         FQSUB(1)=FQSUB(1)+G*DPINV*MENT(J,1)*(QENT(J,1)-Q(1))
         FU(1)=FU(1)+G*DPINV*MENT(J,1)*(UENT(J,1)-U(1))
         FV(1)=FV(1)+G*DPINV*MENT(J,1)*(VENT(J,1)-V(1))
         DO K=1,NTRA
          FTRA(1,K)=FTRA(1,K)+G*DPINV*MENT(J,1)*(TRAENT(J,1,K)-
     1     TRA(1,K))
         END DO
  415   CONTINUE
C
C   ***  CALCULATE TENDENCIES OF POTENTIAL TEMPERATURE AND MIXING RATIO  ***
C   ***               AT LEVELS ABOVE THE LOWEST LEVEL                   ***
C
C   ***  FIRST FIND THE NET SATURATED UPDRAFT AND DOWNDRAFT MASS FLUXES  ***
C   ***                      THROUGH EACH LEVEL                          ***
C
        DO 500 I=2,INB
        DPINV=0.01/(PH(I)-PH(I+1))
        CPINV=1.0/CPN(I)
        AMP1=0.0
        AD=0.0
        IF(I.GE.NK)THEN
         DO 440 K=I+1,INB+1
  440    AMP1=AMP1+M(K)
        END IF
        DO 450 K=1,I
        DO 450 J=I+1,INB+1
         AMP1=AMP1+MENT(K,J)
         MUP(I)=AMP1
  450   CONTINUE
        IF((2.*G*DPINV*AMP1).GE.DELTI)THEN
	 IFLAG=4
c	 CBMF=CBMFOLD-0.1*(2.*G*DPINV*AMP1*DELT-1.)
c         CBMF=MAX(CBMF,0.0)
	END IF
        DO 470 K=1,I-1
        DO 470 J=I,INB
        AD=AD+MENT(J,K)
        MDOWN(I)=AD
  470   CONTINUE
        FT(I)=FT(I)+G*DPINV*(AMP1*(T(I+1)-T(I)+(GZ(I+1)-GZ(I))*
     1   CPINV)-AD*(T(I)-T(I-1)+(GZ(I)-GZ(I-1))*CPINV))
     2   -SIGD*LVCP(I)*EVAP(I)
        FT(I)=FT(I)+G*DPINV*MENT(I,I)*(HP(I)-H(I)+
     1    T(I)*(CPV-CPD)*(Q(I)-QENT(I,I)))*CPINV
        FT(I)=FT(I)+SIGD*WT(I+1)*(CL-CPD)*WATER(I+1)*
     1    (T(I+1)-T(I))*DPINV*CPINV
c        ENTPRO=ENTPRO+SIGD*WT(I)*WATER(I)*RD*(P(I)-P(I+1))/
c     1    (PH(I)*G)
c	ENTPRO=ENTPRO+SIGD*WT(I+1)*WATER(I+1)*(GZ(I+1)-GZ(I))/(G*T(I))
        ENTPRO=ENTPRO+(AMP1-AD)*RD*(Q(NK)-EP(I)*CLW(I)-Q(I))*
     1   (PH(I)-PH(I+1))/P(I)-MP(I)*RD*(QP(I)-Q(I))*(PH(I)-PH(I+1))/
     2   P(I)
c	ENTPRO=ENTPRO-MP(I)*LOG(T(I-1)/T(I))*(QP(I+1)-Q(I+1))*(CL-CPD+
c     1   LV(I)/T(I))-MP(I)*(QP(I+1)-Q(I+1))*(RD-RV)*(P(I)-P(I+1))/
c     2    PH(I+1)
        ENTPRO2=ENTPRO2+(AMP1-AD)*(QS(I+1)-Q(I+1))*LV(I)*(T(I)-T(I+1))/
     1   T(I)**2-MP(I)*(QP(I+1)-Q(I+1))*LV(I)*(T(I)-T(I+1))/T(I)**2
        ENTPRO2=ENTPRO2-(AMP1-AD)*(QS(I+1)-Q(I+1))*RV*(PH(I)-PH(I+1))/
     2   P(I)+MP(I)*(QP(I+1)-Q(I+1))*RV*(PH(I)-PH(I+1))/P(I)
c	ENTPRO2=ENTPRO2+AMP2*RD*(TP(I)-T(I))*(PH(I)-PH(I+1))/(T(I)*P(I))
        FQ(I)=FQ(I)+G*DPINV*(AMP1*(Q(I+1)-Q(I))-
     1    AD*(Q(I)-Q(I-1)))
        FU(I)=FU(I)+G*DPINV*(AMP1*(U(I+1)-U(I))-
     1    AD*(U(I)-U(I-1)))
        FV(I)=FV(I)+G*DPINV*(AMP1*(V(I+1)-V(I))-
     1    AD*(V(I)-V(I-1)))
        DO K=1,NTRA
         FTRA(I,K)=FTRA(I,K)+G*DPINV*(AMP1*(TRA(I+1,K)-
     1    TRA(I,K))-AD*(TRA(I,K)-TRA(I-1,K)))
        END DO
        FQSUB(I)=FQ(I)
        DO 480 K=1,I-1
         AWAT=ELIJ(K,I)-(1.-EP(I))*CLW(I)
         AWAT=MAX(AWAT,0.0)
         FQ(I)=FQ(I)+G*DPINV*MENT(K,I)*(QENT(K,I)-AWAT-Q(I))
         FQDET(I)=G*DPINV*MENT(K,I)*(QENT(K,I)-AWAT-Q(I))
         FU(I)=FU(I)+G*DPINV*MENT(K,I)*(UENT(K,I)-U(I))
         FV(I)=FV(I)+G*DPINV*MENT(K,I)*(VENT(K,I)-V(I))
c -- sb:
C (saturated updrafts resulting from mixing)
         QCOND(I)=QCOND(I)+(ELIJ(K,I)-AWAT)
         NQCOND(I)=NQCOND(I)+1.
c sb --
         DO J=1,NTRA
          FTRA(I,J)=FTRA(I,J)+G*DPINV*MENT(K,I)*(TRAENT(K,I,J)-
     1     TRA(I,J))
         END DO
  480   CONTINUE
        DO 490 K=I,INB
         FQ(I)=FQ(I)+G*DPINV*MENT(K,I)*(QENT(K,I)-Q(I))
         FQDET(I)=FQDET(I)+G*DPINV*MENT(K,I)*(QENT(K,I)-Q(I))
         FU(I)=FU(I)+G*DPINV*MENT(K,I)*(UENT(K,I)-U(I))
         FV(I)=FV(I)+G*DPINV*MENT(K,I)*(VENT(K,I)-V(I))
         DO J=1,NTRA
          FTRA(I,J)=FTRA(I,J)+G*DPINV*MENT(K,I)*(TRAENT(K,I,J)-
     1     TRA(I,J))
         END DO
  490   CONTINUE
        FQ(I)=FQ(I)+SIGD*EVAP(I)+G*(MP(I+1)*
     1    (QP(I+1)-Q(I))-MP(I)*(QP(I)-Q(I-1)))*DPINV
        FQDET(I)=FQDET(I)+SIGD*EVAP(I)+G*(MP(I+1)*
     1    (QP(I+1)-Q(I))-MP(I)*(QP(I)-Q(I-1)))*DPINV
        FU(I)=FU(I)+G*(MP(I+1)*(UP(I+1)-U(I))-MP(I)*
     1    (UP(I)-U(I-1)))*DPINV
        FV(I)=FV(I)+G*(MP(I+1)*(VP(I+1)-V(I))-MP(I)*
     1    (VP(I)-V(I-1)))*DPINV
        DO J=1,NTRA
         FTRA(I,J)=FTRA(I,J)+G*DPINV*(MP(I+1)*(TRAP(I+1,J)-TRA(I,J))-
     1    MP(I)*(TRAP(I,J)-TRAP(I-1,J)))
        END DO
c -- sb:
C (saturated downdrafts resulting from mixing)
        DO K=I+1,INB
         QCOND(I)=QCOND(I)+ELIJ(K,I)
         NQCOND(I)=NQCOND(I)+1.
        ENDDO
C (particular case: no detraining level is found)
        IF (NENT(I).EQ.0) THEN
         QCOND(I)=QCOND(I)+(1-EP(I))*CLW(I)
         NQCOND(I)=NQCOND(I)+1.
        ENDIF
        IF (NQCOND(I).NE.0.) THEN
         QCOND(I)=QCOND(I)/NQCOND(I)
        ENDIF
c sb --
  500   CONTINUE
C
C   *** Adjust tendencies at top of convection layer to reflect  ***
C   ***       actual position of the level zero CAPE             ***
C
        FQOLD=FQ(INB)
        FQ(INB)=FQ(INB)*(1.-FRAC)
        FQ(INB-1)=FQ(INB-1)+FRAC*FQOLD*((PH(INB)-PH(INB+1))/
     1   (PH(INB-1)-PH(INB)))*LV(INB)/LV(INB-1)
        FTOLD=FT(INB)
        FT(INB)=FT(INB)*(1.-FRAC)
        FT(INB-1)=FT(INB-1)+FRAC*FTOLD*((PH(INB)-PH(INB+1))/
     1   (PH(INB-1)-PH(INB)))*CPN(INB)/CPN(INB-1)
        FUOLD=FU(INB)
        FU(INB)=FU(INB)*(1.-FRAC)
        FU(INB-1)=FU(INB-1)+FRAC*FUOLD*((PH(INB)-PH(INB+1))/
     1   (PH(INB-1)-PH(INB)))
        FVOLD=FV(INB)
        FV(INB)=FV(INB)*(1.-FRAC)
        FV(INB-1)=FV(INB-1)+FRAC*FVOLD*((PH(INB)-PH(INB+1))/
     1   (PH(INB-1)-PH(INB)))
        DO K=1,NTRA
         FTRAOLD=FTRA(INB,K)
         FTRA(INB,K)=FTRA(INB,K)*(1.-FRAC)
         FTRA(INB-1,K)=FTRA(INB-1,K)+FRAC*FTRAOLD*(PH(INB)-PH(INB+1))/
     1    (PH(INB-1)-PH(INB))
        END DO
C
C   ***   Very slightly adjust tendencies to force exact   ***
C   ***     enthalpy, momentum and tracer conservation     ***
C
        ENTS=0.0
        UAV=0.0
        VAV=0.0
        DO 680 I=1,INB
         ENTS=ENTS+(CPN(I)*FT(I)+LV(I)*FQ(I))*(PH(I)-PH(I+1))	
         UAV=UAV+FU(I)*(PH(I)-PH(I+1))
         VAV=VAV+FV(I)*(PH(I)-PH(I+1))
  680	CONTINUE
        ENTS=ENTS/(PH(1)-PH(INB+1))
        UAV=UAV/(PH(1)-PH(INB+1))
        VAV=VAV/(PH(1)-PH(INB+1))
        DO 640 I=1,INB
         FT(I)=FT(I)-ENTS/CPN(I)
         FU(I)=(1.-CU)*(FU(I)-UAV)
         FV(I)=(1.-CU)*(FV(I)-VAV)
  640	CONTINUE
        DO 700 K=1,NTRA
         TRAAV=0.0
         DO 690 I=1,INB
          TRAAV=TRAAV+FTRA(I,K)*(PH(I)-PH(I+1))
  690    CONTINUE
         TRAAV=TRAAV/(PH(1)-PH(INB+1))
         DO 695 I=1,INB
          FTRA(I,K)=FTRA(I,K)-TRAAV
  695    CONTINUE
  700	CONTINUE
        DO 750 I=1,INB
         DO 710 K=1,INB
          DET(I)=DET(I)+MENT(K,I)
  710    CONTINUE
         DO 720 K=1,INB
          ENT(I)=ENT(I)+MENT(I,K)
  720    CONTINUE
         ENT(I)=ENT(I)-MENT(I,I)
  750   CONTINUE
C    In-cloud mixing ratio of condensed water :

       DO I=1,ND
        MA(I)=0.0
        WA(I)=0.0
        SIGA(I)=0.0
       ENDDO

       DO I=NK,INB
       DO K=I+1,INB+1
        MA(I)=MA(I)+M(K)
       ENDDO
       ENDDO

       DO I=ICB,INB-1
        AX(I)=0.
        DO J=ICB,I
         AX(I)=AX(I)+RD*(TVP(J)-TV(J))*(PH(J)-PH(J+1))/P(J)
        ENDDO
        IF (AX(I).GT.0.) THEN
         WA(I)=SQRT(2.*AX(I))
        ENDIF
       ENDDO

       DO I=1,NL
       IF (WA(I).GT.0. )
     :    SIGA(I)=MA(I)/WA(I)*RD*TVP(I)/P(I)/100./DELTA
       SIGA(I) = MIN(SIGA(I),1.0) 
       QCONDC(I)=SIGA(I)*CLW(I)*(1.-EP(I))
     :          + (1.-SIGA(I))*QCOND(I)
       ENDDO
C
C   ***           RETURN           ***
C
        RETURN
C
        END
C
C ---------------------------------------------------------------------------
C
        SUBROUTINE TLIFT(P,T,Q,QS,GZ,ICB,NK,TVP,TPK,CLW,ND,NL,KK)
        REAL GZ(ND),TPK(ND),CLW(ND),P(ND)
        REAL T(ND),Q(ND),QS(ND),TVP(ND),LV0
C
C   ***   ASSIGN VALUES OF THERMODYNAMIC CONSTANTS     ***
C
        CPD=1005.7
        CPV=1870.0
c        CL=4190.0 
        CL=2500.0
        RV=461.5
        RD=287.04
        LV0=2.501E6
C
        CPVMCL=CL-CPV
        EPS=RD/RV
        EPSI=1./EPS
C
C   ***  CALCULATE CERTAIN PARCEL QUANTITIES, INCLUDING STATIC ENERGY   ***
C
	  ALV=LV0-CPVMCL*(T(NK)-273.15)
        AH0=(CPD*(1.-Q(NK))+CL*Q(NK))*T(NK)+Q(NK)*ALV+GZ(NK)
        CPP=CPD*(1.-Q(NK))+Q(NK)*CPV
        CPINV=1./CPP
C
        IF(KK.EQ.1)THEN
C
C   ***   CALCULATE LIFTED PARCEL QUANTITIES BELOW CLOUD BASE   ***
C
        DO 50 I=1,ICB-1
         CLW(I)=0.0
   50   CONTINUE
        DO 100 I=NK,ICB-1
         TPK(I)=T(NK)-(GZ(I)-GZ(NK))*CPINV
         TVP(I)=TPK(I)*(1.+Q(NK)*EPSI)
  100   CONTINUE
        END IF
C
C    ***  FIND LIFTED PARCEL QUANTITIES ABOVE CLOUD BASE    ***
C
        NST=ICB
        NSB=ICB
        IF(KK.EQ.2)THEN  
         NST=NL
         NSB=ICB+1
        END IF
        DO 300 I=NSB,NST
         TG=T(I)
         QG=QS(I)
         ALV=LV0-CPVMCL*(T(I)-273.15)
         DO 200 J=1,2
          S=CPD+ALV*ALV*QG/(RV*T(I)*T(I))
          S=1./S
          AHG=CPD*TG+(CL-CPD)*Q(NK)*T(I)+ALV*QG+GZ(I)
          TG=TG+S*(AH0-AHG)
          TG=MAX(TG,35.0)
          TC=TG-273.15
          DENOM=243.5+TC
          IF(TC.GE.0.0)THEN  
           ES=6.112*EXP(17.67*TC/DENOM)
          ELSE  
           ES=EXP(23.33086-6111.72784/TG+0.15215*LOG(TG))
          END IF  
          QG=EPS*ES/(P(I)-ES*(1.-EPS))
  200    CONTINUE
         ALV=LV0-CPVMCL*(T(I)-273.15)
         TPK(I)=(AH0-(CL-CPD)*Q(NK)*T(I)-GZ(I)-ALV*QG)/CPD
         CLW(I)=Q(NK)-QG
         CLW(I)=MAX(0.0,CLW(I))
         RG=QG/(1.-Q(NK))
         TVP(I)=TPK(I)*(1.+RG*EPSI)
  300   CONTINUE
        RETURN
        END
C
       SUBROUTINE PCMIN(SST,PSL,P,T,R,NA,N,PMIN,VMAX,IFL,TOOUT)
C
C   Revised on 9/24/2005 to fix convergence problems at high pressure
C
C   ***   This subroutine calculates the maximum wind speed        ***
C   ***             and mimimum central pressure                   ***
C   ***    achievable in tropical cyclones, given a sounding       ***
C   ***             and a sea surface temperature.                 ***
C
C  INPUT:   SST: Sea surface temperature in C
C
C           PSL: Sea level pressure (mb)
C
C           P,T,R: One-dimensional arrays of dimension NA
C             containing pressure (mb), temperature (C),
C             and mixing ratio (g/kg). The arrays MUST be
C             arranged so that the lowest index corresponds
C             to the lowest model level, with increasing index
C             corresponding to decreasing pressure. The temperature
C             sounding should extend to at least the tropopause and 
C             preferably to the lower stratosphere, however the
C             mixing ratios are not important above the boundary
C             layer. Missing mixing ratios can be replaced by zeros.
C
C           NA: The dimension of P,T and R
C
C           N:  The actual number of points in the sounding
C                (N is less than or equal to NA)
C
C  OUTPUT:  PMIN is the minimum central pressure, in mb
C
C           VMAX is the maximum surface wind speed, in m/s
C                  (reduced to reflect surface drag)
C
C           IFL is a flag: A value of 1 means OK; a value of 0
C              indicates no convergence (hypercane); a value of 2
C              means that the CAPE routine failed.
C
C-----------------------------------------------------------------------------
	REAL T(NA), P(NA), R(NA)
C
C   ***   Adjustable constant: Ratio of C_k to C_D    ***
C
	CKCD=0.9
C
C   ***   Adjustable constant for buoyancy of displaced parcels:  ***
C   ***    0=Reversible ascent;  1=Pseudo-adiabatic ascent        ***
C
      SIG=0.0
C
C   ***  Adjustable switch: if IDISS = 0, no dissipative heating is   ***
C   ***     allowed; otherwise, it is                                 ***
C
	IDISS=1
C
C   ***  Exponent, b, in assumed profile of azimuthal velocity in eye,   ***
C   ***   V=V_m(r/r_m)^b. Used only in calculation of central pressure   ***
C
	b=2.0
C
C   *** Set level from which parcels lifted   ***
C
	NK=1
C
C   *** Factor to reduce gradient wind to 10 m wind
C
	VREDUC=0.8
C
C   ***   Normalize certain quantities   ***
C
	SSTK=SST+273.15
	ES0=6.112*EXP(17.67*SST/(243.5+SST))
c	DO 40 I=1,N
c	 R(I)=R(I)*0.001
c	 T(I)=T(I)+273.15
c  40	CONTINUE
C
C   ***   Default values   ***
C
      VMAX=0.0
	PMIN=PSL 
	IFL=1
C
	NP=0
	PM=950.0
C
C   ***   Find environmental CAPE *** 
C
      TP=T(NK)
      RP=R(NK)
      PP=P(NK)
      CALL CAPE(TP,RP,PP,T,R,P,NA,N,SIG,CAPEA,TOA,IFLAG)
      IF(IFLAG.NE.1)IFL=2
C
C   ***   Begin iteration to find mimimum pressure   ***
C
  100 CONTINUE
C
C   ***  Find CAPE at radius of maximum winds   ***
C
      TP=T(NK)
      PP=MIN(PM,1000.0)
      RP=0.622*R(NK)*PSL/(PP*(0.622+R(NK))-R(NK)*PSL)
      CALL CAPE(TP,RP,PP,T,R,P,NA,N,SIG,CAPEM,TOM,IFLAG) 
      IF(IFLAG.NE.1)IFL=2
      RAT=SSTK/TOM
      IF(IDISS.EQ.0)RAT=1.0
C
C  ***  Find saturation CAPE at radius of maximum winds   ***
C
      TP=SSTK
      PP=MIN(PM,1000.0)
      RP=0.622*ES0/(PP-ES0)
      CALL CAPE(TP,RP,PP,T,R,P,NA,N,SIG,CAPEMS,TOMS,IFLAG)
	TOOUT=0.5*(TOMS+TOM)
      IF(IFLAG.NE.1)IFL=2
C
C  ***  Initial estimate of minimum pressure   ***
C
      RS0=RP
      TV1=T(1)*(1.+R(1)/0.622)/(1.+R(1))
	TVAV=0.5*(TV1+SSTK*(1.+RS0/0.622)/(1.+RS0))
C	CAT=0.5*CKCD*RAT*(CAPEMS-CAPEM)
	CAT=CAPEM-CAPEA+0.5*CKCD*RAT*(CAPEMS-CAPEM)
	CAT=MAX(CAT,0.0)
	PNEW=PSL*EXP(-CAT/(287.04*TVAV))
C
C   ***  Test for convergence   ***
C
	IF(ABS(PNEW-PM).GT.0.2)THEN
	 PM=PNEW
	 NP=NP+1
	 IF(NP.GT.1000.OR.PM.LT.400.0)THEN
	  PMIN=PSL
	  IFL=0
	  GOTO 900
	 END IF
	 GOTO 100
	ELSE
	 CATFAC=0.5*(1.+1./b)
C	 CAT=CKCD*RAT*CATFAC*(CAPEMS-CAPEM)
	 CAT=CAPEM-CAPEA+CKCD*RAT*CATFAC*(CAPEMS-CAPEM)
	 CAT=MAX(CAT,0.0)
	 PMIN=PSL*EXP(-CAT/(287.04*TVAV))
	END IF
  900	CONTINUE
	FAC=MAX(0.0,(CAPEMS-CAPEM))
	VMAX=VREDUC*SQRT(CKCD*RAT*FAC)
C
C   ***  Renormalize sounding arrays   ***
C	
c	DO 910 I=1,N
c	 R(I)=R(I)*1000.0
c	 T(I)=T(I)-273.15
c910	CONTINUE
C
	RETURN
	END
C        
      SUBROUTINE CAPE(TP,RP,PP,T,R,P,ND,N,SIG,CAPED,TOB,IFLAG)
C
C     This subroutine calculates the CAPE of a parcel with pressure PP (mb), 
C       temperature TP (K) and mixing ratio RP (gm/gm) given a sounding
C       of temperature (T in K) and mixing ratio (R in gm/gm) as a function
C       of pressure (P in mb). ND is the dimension of the arrays T,R and P,
C       while N is the actual number of points in the sounding. CAPED is
C       the calculated value of CAPE and TOB is the temperature at the
C       level of neutral buoyancy.  IFLAG is a flag
C       integer. If IFLAG = 1, routine is successful; if it is 0, routine did
C       not run owing to improper sounding (e.g.no water vapor at parcel level).
C       IFLAG=2 indicates that routine did not converge.                 
C
      REAL T(ND),R(ND),P(ND),TVRDIF(100)   
      REAL NA
C
C   ***   Default values   ***
C      
      CAPED=0.0
      TOB=T(1)
      IFLAG=1
C
C   ***   Check that sounding is suitable    ***
C
      IF(RP.LT.1.0E-6.OR.TP.LT.200.0)THEN
       IFLAG=0
       RETURN
      END IF            
C
C   ***   Assign values of thermodynamic constants     ***
C
      CPD=1005.7
      CPV=1870.0
C      CL=4190.0
      CL=2500.0
      CPVMCL=CPV-CL
      RV=461.5
      RD=287.04
      EPS=RD/RV
      ALV0=2.501E6
C
C   ***  Define various parcel quantities, including reversible   ***
C   ***                       entropy, S.                         ***
C                           
      TPC=TP-273.15
      ESP=6.112*EXP(17.67*TPC/(243.5+TPC))
      EVP=RP*PP/(EPS+RP)
      RH=EVP/ESP
	RH=MIN(RH,1.0)
      ALV=ALV0-CPVMCL*TPC
      S=(CPD+RP*CL)*LOG(TP)-RD*LOG(PP-EVP)+
     1   ALV*RP/TP-RP*RV*LOG(RH)            
C
C   ***  Find lifted condensation pressure, PLCL   ***
C     
	CHI=TP/(1669.0-122.0*RH-TP)
	PLCL=PP*(RH**CHI)
C
C   ***  Begin updraft loop   ***
C
	NCMAX=0
	DO J=1,N
	 TVRDIF(J)=0.0
	END DO
C
	JMIN=1E6
	DO 200 J=1,N
C
C    ***   Don't bother lifting parcel above 60 mb and skip sections of sounding below parcel level  ***
C
      IF(P(J).LT.59.0.OR.P(J).GE.PP)GOTO 200
C
	JMIN=MIN(JMIN,J)
C
C    ***  Parcel quantities below lifted condensation level   ***
C	 
	 IF(P(J).GE.PLCL)THEN
	  TG=TP*(P(J)/PP)**(RD/CPD)
	  RG=RP
C
C   ***   Calculate buoyancy   ***
C  
	  TLVR=TG*(1.+RG/EPS)/(1.+RG)
	  TVRDIF(J)=TLVR-T(J)*(1.+R(J)/EPS)/(1.+R(J))
	 ELSE
C
C   ***  Parcel quantities above lifted condensation level  ***
C	 
	  TG=T(J)          
	  TJC=T(J)-273.15 
	  ES=6.112*EXP(17.67*TJC/(243.5+TJC)) 
	  RG=EPS*ES/(P(J)-ES)
C
C   ***  Iteratively calculate lifted parcel temperature and mixing   ***
C   ***                ratio for reversible ascent                    ***
C
	  NC=0
  120	  CONTINUE
	  NC=NC+1
C
C   ***  Calculate estimates of the rates of change of the entropy    ***
C   ***           with temperature at constant pressure               ***
C  
	  ALV=ALV0-CPVMCL*(TG-273.15)
	  SL=(CPD+RP*CL+ALV*ALV*RG/(RV*TG*TG))/TG
	  EM=RG*P(J)/(EPS+RG)
	  SG=(CPD+RP*CL)*LOG(TG)-RD*LOG(P(J)-EM)+
     1      ALV*RG/TG
	  IF(NC.LT.3)THEN
	   AP=0.3
	  ELSE
	   AP=1.0
	  END IF
	  TGNEW=TG+AP*(S-SG)/SL  
C
C   ***   Test for convergence   ***
C
	  IF(ABS(TGNEW-TG).GT.0.001)THEN
	   TG=TGNEW
	   TC=TG-273.15
	   ENEW=6.112*EXP(17.67*TC/(243.5+TC))
C
C   ***   Bail out if things get out of hand   ***
C
	   IF(NC.GT.500.OR.ENEW.GT.(P(J)-1.0))THEN
            IFLAG=2
            RETURN
	   END IF
	   RG=EPS*ENEW/(P(J)-ENEW)           
	   GOTO 120
	  END IF
	  NCMAX=MAX(NC,NCMAX)
C
C   *** Calculate buoyancy   ***
C
        RMEAN=SIG*RG+(1.-SIG)*RP
	  TLVR=TG*(1.+RG/EPS)/(1.+RMEAN)
	  TVRDIF(J)=TLVR-T(J)*(1.+R(J)/EPS)/(1.+R(J))
	 END IF
  200	CONTINUE
C
C  ***  Begin loop to find NA, PA, and CAPE from reversible ascent ***
C
	NA=0.0
	PA=0.0
C
C   ***  Find maximum level of positive buoyancy, INB    ***
C
	INB=1
	DO 550 J=N,JMIN,-1
	 IF(TVRDIF(J).GT.0.0)INB=MAX(INB,J)
  550	CONTINUE
	IF(INB.EQ.1)RETURN
C
C   ***  Find positive and negative areas and CAPE  ***
C
	IF(INB.GT.1)THEN
	 DO 600 J=JMIN+1,INB
	  PFAC=RD*(TVRDIF(J)+TVRDIF(J-1))*(P(J-1)-P(J))/(P(J)+P(J-1))
	  PA=PA+MAX(PFAC,0.0)
	  NA=NA-MIN(PFAC,0.0)
  600	 CONTINUE
C
C   ***   Find area between parcel pressure and first level above it ***
C
	PMA=(PP+P(JMIN)) 
	PFAC=RD*(PP-P(JMIN))/PMA
	PA=PA+PFAC*MAX(TVRDIF(JMIN),0.0)
	NA=NA-PFAC*MIN(TVRDIF(JMIN),0.0)
C
C   ***   Find residual positive area above INB and TO  ***
C
       PAT=0.0
       TOB=T(INB)
       IF(INB.LT.N)THEN
        PINB=(P(INB+1)*TVRDIF(INB)-P(INB)*TVRDIF(INB+1))/
     1   (TVRDIF(INB)-TVRDIF(INB+1))
        PAT=RD*TVRDIF(INB)*(P(INB)-PINB)/(P(INB)+PINB)
	  TOB=(T(INB)*(PINB-P(INB+1))+T(INB+1)*(P(INB)-PINB))/
     1    (P(INB)-P(INB+1))
       END IF
C
C   ***   Find CAPE  ***
C            
	 CAPED=PA+PAT-NA
	 CAPED=MAX(CAPED,0.0)
	END IF
C
	RETURN
	END
C
      SUBROUTINE CLOUDS_SUB_LS_40(ND,R,RS,T,P,PH,DT,QSUBGRID
     :                        ,CLDF,CLDQ,PRADJ,FTADJ,FRADJ,SIGSUB)
      implicit none
C
C==========================================================================
C
C               CLOUDS_SUB_LS  version 4.0
C
C Purpose:
C --------
C
C   1) Call of the cloud parameterization
C
C   2) Large-scale super-saturation adjustment:
C     - condense water that exceeds saturation
C     - precipitate a fraction of that condensed water
C     - adjust the final in-cloud water content 
C
C Inputs:
C -------
C  ND----------: Number of vertical levels
C  R--------ND-: Grid-box average of the total water mixing ratio [kg/kg]
C  RS-------ND-: Mean saturation humidity mixing ratio within the gridbox [kg/kg]
C  T--------ND-: Grid-box average temperature [K]
C  P-------ND+1: Pressure at mid-levels [mb]
C  PH------ND+1: Pressure at interface levels [mb]
C  DT----------: Timestep [seconds]
C  QSUBGRID-ND-: in-cloud mixing ratio of cloud condensate [kg/kg] from CONVECT
C
C Outputs:
C --------
C  CLDF-----ND-: cloud fraction [0-1]
C  CLDQ-----ND-: in-cloud mixing ratio of condensed water [kg/kg]
C  PRADJ----ND-: precipitation associated with the LS super-saturation [mm/day]
C  FTADJ----ND-: temperature tendency associated with the LS adjustment [K/s]
C  FRADJ----ND-: total water tendency associated with the LS adjustment [kg/kg/s]
C  SIGSUB---ND-: SQRT(variance) of total water [kg/kg] (diagnostic only)
C
C Written by:
C -----------
C  Sandrine Bony (MIT & LMD/CNRS)   -    August 1999
C
C Modified by:
C ------------
C
C============================================================================
      integer NA,ND,I
      parameter (NA=120)
      real DT, PRADJ
      real TCA,ELACRIT,ALV,CPN,TNEW,RNEW,EP
      real R(ND),RS(ND),T(ND),P(ND),PH(ND+1)
     :    ,QSUBGRID(ND),CLDF(ND),CLDQ(ND)
     :    ,FTADJ(ND),FRADJ(ND),RNEWLS(NA),TNEWLS(NA),QLSP(NA)
     :    ,EPLS(NA)
c a retirer dans version pub:
      real RSSURF,CLDL(NA),CLDA(NA),CLDK(NA),CLDS(NA)
      real BIDON(ND)
      real SIGSUB(ND)
      COMMON / cldpara / CLDL, CLDA, CLDK, CLDS
C--------------------------------------------------------------------
C Thermodynamical constants:
C
        REAL CPD,CPV,CL,RD,RV,LV0,G,ROWL,EPS,EPSI,CPVMCL
        PARAMETER (CPD=1005.7, CPV=1870.0, CL=2500.0, RD=287.04)
        PARAMETER (RV=461.5, LV0=2.501E6, G=9.8, ROWL=1000.0 )
        PARAMETER (EPS=RD/RV, EPSI=1./EPS, CPVMCL=CL-CPV)
C--------------------------------------------------------------------
C Microphysical parameters:
C (here, we use the same values as in the convection scheme)
C
        REAL TLCRIT,ELCRIT,EPMAX
        PARAMETER (TLCRIT=-55.0, ELCRIT=0.0011, EPMAX=0.999)
C--------------------------------------------------------------------
C  Initialize output arrays
C
        PRADJ=0.0
C
        DO I = 1, ND
         FRADJ(I)  = 0.0
         FTADJ(I)  = 0.0
         TNEWLS(I) = T(I)
         RNEWLS(I) = R(I)
         QLSP(I)   = 0.0
         EPLS(I)   = 0.0
         BIDON(I)  = 0.0 ! used for all-or-nothing only
         SIGSUB(I) = 0.0 
        ENDDO
C-------------------------------------------------------------------
C If the cloud parameterization has a "convective" type of closure,
C then compute precipitation efficiencies associated with
C large-scale precipitation (independent on the presence of
C subgrid-scale variability), and then compute the cloud
C fraction and in-cloud water content associated with the
C combination of subgrid+large-scale condensation.
C
C
        DO 9999 I = ND, 1, -1
C
C Calculate large-scale condensation and precipitation:
C
         IF(R(I).GT.RS(I))THEN
C
C  Precipitation efficiencies:
C
         EP=0.0
         TCA=T(I)-273.15
         IF(TCA.GE.0.0)THEN
          ELACRIT=ELCRIT
         ELSE
          ELACRIT=ELCRIT*(1.0-TCA/TLCRIT)
         END IF
         ELACRIT=MAX(ELACRIT,0.0)
         EP= EPMAX * (1.0-ELACRIT/MAX(R(I)-RS(I),1.0E-8))
         EP=MAX(EP,0.0)
         EP=MIN(EP,EPMAX)
         EPLS(I) = EP ! large-scale precipitation efficiency
         QLSP(I) = EP*(R(I)-RS(I)) ! precipitated cloud water
C
C  Adjust temperature and humidity profiles:
C
          TCA=T(I)-273.15
          ALV=LV0-CPVMCL*TCA
          CPN=CPD*(1.-R(I))+CPV*R(I)
          TNEW=(ALV*(EP*R(I)+RS(I)*(ALV/(RV*T(I))-EP))
     :       +CPN*T(I))/
     1     (CPN+ALV*ALV*RS(I)/(RV*T(I)*T(I)))
          RNEW=RS(I)*(1.+(TNEW-T(I))*ALV/(RV*T(I)*T(I)))
     :       + (1.-EP)*(R(I)-RS(I))
          TNEWLS(I) = TNEW
          RNEWLS(I) = RNEW
          FRADJ(I)=FRADJ(I)+(RNEW-R(I))/DT
          FTADJ(I)=FTADJ(I)-ALV*(RNEW-R(I))/DT/CPN
          PRADJ=PRADJ-100.0*(PH(I)-PH(I+1))*(RNEW-R(I))/DT
     :            *1000.0*3600.0*24.0/(ROWL*G)
         END IF  ! R>RS
9999    CONTINUE
C-------------------------------------------------------------
C  Cloud parameterization:
C
c GNO PDF and convect-closure: 
         CALL CLOUDS_GNO_40(ND,R,RS,QSUBGRID,CLDF,CLDQ)
C--------------------------------------------------------------------
C  Remove large-scale precipitation from the cloud water content:
C
           DO I = 1, ND
            CLDQ(I) = CLDQ(I) - QLSP(I)
           ENDDO

         RETURN
         END
C
      SUBROUTINE CLOUDS_GNO_40(ND,R,RS,QSUB,CLDF,CLDQ)
      IMPLICIT NONE
C	   
C===================================================================================
C
C				CLOUDS_GNO 	version 4.0
C
C Purpose:
C --------
C
C   Parameterization of the cloudiness (cloud amount, cloud water content)
C   associated with cumulus convection.
C
C Principle:
C ----------
C
C   This cloud parameterization predicts the cloudiness that is associated with 
C   the presence of condensation within a large-scale domain: this condensation 
C   may be produced at the subgrid-scale by cumulus convection and at the 
C   large-scale by super-saturation.
C
C   IMPORTANT: in the present version of the scheme, the only source of subgrid-scale 
C   condensation that is considered is cumulus convection (condensation associated
C   with boundary layer turbulence, for instance, is not considered). 
C
C   The cloud fraction and the in-cloud water content are predicted by a
C   statistical approach. The subgrid-scale variability of total water 
C   (vapor + condensed) within the gridbox is described by a generalized 
C   log-normal Probability Distribution Function (PDF) whose mean, variance 
C   and skewness coefficient are predicted. The predictors are: 
C     1) the local concentration of condensed water that is produced at
C        the subgrid-scale by convection (output of the convection scheme) 
C     2) the saturation deficit or excess of the environment 
C     3) the domain-averaged mixing ratio of total water
C   Note that we impose the distribution of total water to be bounded by zero. 
C   On the other hand, no upper bound of the distribution is considered in this 
C   version of the scheme.
C
C   If no subgrid-scale condensation occurs within the domain, the scheme
C   becomes equivalent to an "all-or-nothing" large-scale saturation scheme.
C
C Inputs:
C -------
C
C  ND----------: Number of vertical levels
C  R--------ND-: Domain-averaged mixing ratio of total water 
C  RS-------ND-: Mean saturation humidity mixing ratio within the gridbox
C  QSUB-----ND-: Mixing ratio of condensed water within clouds associated
C                with SUBGRID-SCALE condensation processes (here, it is
C                predicted by the convection scheme)
C Outputs:
C --------
C
C  CLDF-----ND-: cloud fractional area (0-1)
C  CLDQ-----ND-: in-cloud mixing ratio of condensed water (kg/kg)
C
C CALL command:
C -------------
C
C     CALL CLOUDS_GNO(ND,R,RS,QSUBGRID,CLDF,CLDQ)
C
C Reference:
C ----------
C
C   Bony, S and K A Emanuel, 2001: A parameterization of the cloudiness
C       associated with cumulus convection; Evaluation using TOGA COARE data.
C	    J. Atmos. Sci., accepted.
C
C Written by:
C -----------
C
C  Sandrine Bony (MIT & LMD/CNRS; bony@wind.mit.edu) -  July 2000
C
C  Difference with version 1.0:
C	numerical method of resolution of equation 9
C	version 1.0: use a Gaussian PDF when erf(v)->1
C       version 2.0: use an asymptotic expression of erf(v) instead of a Gaussian PDF
C=====================================================================================
c
c -- input/output arguments of the subroutine:
      INTEGER ND
      REAL R(ND),  RS(ND), QSUB(ND), CLDF(ND), CLDQ(ND)
      REAL alpha, lambda, kew, skew, sigs
c -- lower bound of the PDF of total water:
c
      REAL PB
      PARAMETER ( PB = 0.0 )
c
c -- parameters controlling the iteration:
c --    nmax    : maximum nb of iterations (hopefully never reached!)
c --    epsilon : accuracy of the numerical resolution (here 2.0%)
c --    vmax    : v-value above which we use an asymptotic expression for ERF(v)
      INTEGER nmax, niter
      PARAMETER ( nmax = 10)
      REAL epsilon, vmax0, vmax
      PARAMETER ( epsilon = 0.02, vmax0 = 2.0 ) 
c -- gardes-fou:
      REAL min_mu, min_Q
      PARAMETER ( min_mu =  1.e-12, min_Q=1.e-12 )   
c -- misc:
      INTEGER K, n, m
      REAL*8 mu, qsat, delta, beta 
      REAL*8 xx, aux, coeff, block, dist, fprime, det
      REAL*8 pi, u, v, erfu, erfv, xx1, xx2, erf
      LOGICAL lconv
c ----------------------------------------------------------------------------------
c ----------------------------------------------------------------------------------
      pi = ACOS(-1.)
c
c -- loop over vertical levels :
c
      DO 500 K = 1, ND
c
      mu = R(K)
      mu = MAX(mu,min_mu)
      qsat = RS(K) 
      qsat = MAX(qsat,min_mu)
      delta = log(mu/qsat)
      IF ( QSUB(K) .lt. min_Q ) THEN
C===========================================================================
C  If no condensation is produced at the subgrid-scale:
C
C  -> the scheme becomes equivalent to a "large-scale condensation scheme"
C     ie: cldf = H(mu-qsat) and cldq = (mu-qsat)*H(mu-qsat)
C     where H is the Heaviside function.
C     (in the absence of subgrid-scale condensation, the generalized
C      log-normal PDF becomes equivalent to a gaussian PDF of variance
C      zero, i.e. it becomes equivalent to a Dirac function and the 
C      cumulative distribution function becomes an Heaviside function).
C===========================================================================
        CLDQ(K) = MAX( 0.0, mu-qsat )
        CLDF(K) = CLDQ(K) / MAX( CLDQ(K), min_mu )
        lambda = mu
        alpha  = 0.0
        kew    = 0.0   
        skew   = 0.0   
        sigs   = 0.0   
      ELSE 
C===========================================================================
C  Some condensation is produced at the subgrid-scale: 
C  (presence of subgrid-scale variability):
C
C Use the (iterative) numerical method of Newton to determine the parameters 
C that characterize the PDF of total water.
C
C Remark 1: the accuracy of the resolution is controlled by "epsilon"
C Remark 2: in GCMs, this numerical method may be too much CPU-time consuming.
C In that case, it may be more appropriate to substitute it by a tabulation 
C of equations 9 and 11 (see the Bony-Emanuel article cited in introduction).
C===========================================================================
         lconv  = .FALSE. ! flag for numerical convergence
         niter = 0
c +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c        PDF = generalized log-normal distribution (GNO):
c   (k<0 if a lower bound is considered for the PDF of total water)
c 
c   -> determine x (the parameter k of the GNO PDF) 
c      such that the contribution of subgrid-scale processes to the 
c      in-cloud water content is equal to QSUB(K)
c
c NB: the "error function" is called ERF or DERF (in double precision)
c +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        beta = QSUB(K)/mu + EXP( -MIN(0.0,delta) )
        vmax = vmax0
        IF ( .NOT. lconv ) then ! new     
c --  roots of equation v > vmax:
        det = delta + vmax**2.
        if (det.LE.0.0) vmax = vmax0 + 1.0
        det = delta + vmax**2.
        if (det.LE.0.) then
          xx = -0.0001
        else 
         xx1 = -SQRT(2.)*vmax*(1.0-SQRT(1.0+delta/(vmax**2.)))
         xx2 = -SQRT(2.)*vmax*(1.0+SQRT(1.0+delta/(vmax**2.)))
         xx = 1.01 * xx1
         if ( xx1 .GE. 0.0 ) xx = 0.5*xx2
        endif
        if (delta.LT.0.)
     :    xx = -0.5*SQRT(log(2.)) ! test comme avant
        DO n = 1, nmax ! iteration loop
          u = delta/(xx*sqrt(2.)) + xx/(2.*sqrt(2.))
          v = delta/(xx*sqrt(2.)) - xx/(2.*sqrt(2.))
          IF ( v .GT. vmax ) THEN 
            IF (     ABS(u)  .GT. vmax 
     :          .AND.  delta .LT. 0. ) THEN
c -- use asymptotic expression of erf for u and v large:
c ( -> analytic solution for xx )
             aux = 2.0*delta*(1.-beta*EXP(delta))
     :                       /(1.+beta*EXP(delta))
             xx = -SQRT(aux)
             block = EXP(-v*v) / v / sqrt(pi)
             dist = 0.0
             fprime = 1.0
            ELSE
c -- erfv -> 1.0, use an asymptotic expression of erfv for v large:
             erfu = ERF(u)
             aux = sqrt(pi) * (1.0-erfu) * EXP(v*v)
             coeff = 1.0 - 1./2./(v**2.) + 3./4./(v**4.)
             block = coeff * EXP(-v*v) / v / sqrt(pi)
             dist = v * aux / coeff - beta
             fprime = 2.0 / xx * (v**2.)
     :           * ( coeff*EXP(-delta) - u * aux )
     :           / coeff / coeff         
            ENDIF ! ABS(u)
          ELSE
c -- general case:
           erfu = ERF(u)
           erfv = ERF(v)
           block = 1.0-erfv
           dist = (1.0 - erfu) / (1.0 - erfv) - beta
           fprime = 2. /sqrt(pi) /xx /(1.0-erfv)**2.
     :           * (   (1.0-erfv)*v*EXP(-u*u)
     :               - (1.0-erfu)*u*EXP(-v*v) )
          ENDIF ! x
c -- numerical convergence reached?
        if ( ABS(dist/beta) .LT. epsilon ) then ! convergence 
          lconv = .TRUE. ! numerical convergence reached with GNO PDF
ccc          write(*,*) 'CV NEWTON GNO FOR K, N = ',K,n
c parameters that characterize the PDF:
           kew    = xx
           aux    = EXP(kew*kew)
           sigs   = mu * SQRT( aux-1.0 )
           sigs   = MAX(sigs,0.0)
           lambda = mu / SQRT (aux)
           alpha  = - lambda * kew
           skew   = ( 2.0 + aux**3. - 3*aux ) * (mu/sigs)**3.
c deduce the cloud fraction and the in-cloud water content:
c --------------------------------------------------------------
        CLDF(K) = 0.5 * block
        CLDQ(K) = QSUB(K) + MAX(mu-qsat,0.0)
        GOTO 100
        else
           xx = xx - dist/fprime
        endif
        ENDDO ! n
 100    continue
        if (.NOT. lconv) then 
          write(*,*) 'NO CV in CLOUDS_GNO: K,mu,qsub,qsat,n,niter: '
     :            ,K,mu*1000.,QSUB(K)*1000.,qsat*1000.,n,niter
          write(*,*) 'gamma,beta,error = ',delta,beta,ABS(dist/beta)
c all-or-nothing scheme in that (exceptional) case, may be improved later on:
          CLDQ(K) = MAX( 0.0, mu-qsat )
          CLDF(K) = CLDQ(K) / MAX( CLDQ(K), min_mu )
        endif
        ENDIF ! lconv
      ENDIF ! qsub
 500  CONTINUE  ! K
       RETURN
       END
C
      function erfc(x)
      parameter(
     &  pv= 9.20888710e+00,
     &  ph= 5.07254732e+00,
     &  p0= 3.86642217e-01,
     &  p1= 1.52430177e-01,
     &  p2= 2.38149125e-02,
     &  p3= 1.30227291e-03,
     &  q0= 1.16381965e-01,
     &  q1= 1.04753802e+00,
     &  q2= 2.92132156e+00,
     &  q3= 6.02608434e+00)
      y=x*x
      y=exp(-y)*x*(p3/(y+q3)+p2/(y+q2)
     &  +p1/(y+q1)+p0/(y+q0))
      if(x.lt.ph) y=y+2/(exp(pv*x)+1)
      erfc=y
      end
!
! error function in double precision
!
      function erf(x)
      implicit real*8 (a - h, o - z)
      dimension a(0 : 64), b(0 : 64)
      data (a(i), i = 0, 12) / 
     &    0.00000000005958930743d0, -0.00000000113739022964d0, 
     &    0.00000001466005199839d0, -0.00000016350354461960d0, 
     &    0.00000164610044809620d0, -0.00001492559551950604d0, 
     &    0.00012055331122299265d0, -0.00085483269811296660d0, 
     &    0.00522397762482322257d0, -0.02686617064507733420d0, 
     &    0.11283791670954881569d0, -0.37612638903183748117d0, 
     &    1.12837916709551257377d0 / 
      data (a(i), i = 13, 25) / 
     &    0.00000000002372510631d0, -0.00000000045493253732d0, 
     &    0.00000000590362766598d0, -0.00000006642090827576d0, 
     &    0.00000067595634268133d0, -0.00000621188515924000d0, 
     &    0.00005103883009709690d0, -0.00037015410692956173d0, 
     &    0.00233307631218880978d0, -0.01254988477182192210d0, 
     &    0.05657061146827041994d0, -0.21379664776456006580d0, 
     &    0.84270079294971486929d0 / 
      data (a(i), i = 26, 38) / 
     &    0.00000000000949905026d0, -0.00000000018310229805d0, 
     &    0.00000000239463074000d0, -0.00000002721444369609d0, 
     &    0.00000028045522331686d0, -0.00000261830022482897d0, 
     &    0.00002195455056768781d0, -0.00016358986921372656d0, 
     &    0.00107052153564110318d0, -0.00608284718113590151d0, 
     &    0.02986978465246258244d0, -0.13055593046562267625d0, 
     &    0.67493323603965504676d0 / 
      data (a(i), i = 39, 51) / 
     &    0.00000000000382722073d0, -0.00000000007421598602d0, 
     &    0.00000000097930574080d0, -0.00000001126008898854d0, 
     &    0.00000011775134830784d0, -0.00000111992758382650d0, 
     &    0.00000962023443095201d0, -0.00007404402135070773d0, 
     &    0.00050689993654144881d0, -0.00307553051439272889d0, 
     &    0.01668977892553165586d0, -0.08548534594781312114d0, 
     &    0.56909076642393639985d0 / 
      data (a(i), i = 52, 64) / 
     &    0.00000000000155296588d0, -0.00000000003032205868d0, 
     &    0.00000000040424830707d0, -0.00000000471135111493d0, 
     &    0.00000005011915876293d0, -0.00000048722516178974d0, 
     &    0.00000430683284629395d0, -0.00003445026145385764d0, 
     &    0.00024879276133931664d0, -0.00162940941748079288d0, 
     &    0.00988786373932350462d0, -0.05962426839442303805d0, 
     &    0.49766113250947636708d0 / 
      data (b(i), i = 0, 12) / 
     &    -0.00000000029734388465d0, 0.00000000269776334046d0, 
     &    -0.00000000640788827665d0, -0.00000001667820132100d0, 
     &    -0.00000021854388148686d0, 0.00000266246030457984d0, 
     &    0.00001612722157047886d0, -0.00025616361025506629d0, 
     &    0.00015380842432375365d0, 0.00815533022524927908d0, 
     &    -0.01402283663896319337d0, -0.19746892495383021487d0, 
     &    0.71511720328842845913d0 / 
      data (b(i), i = 13, 25) / 
     &    -0.00000000001951073787d0, -0.00000000032302692214d0, 
     &    0.00000000522461866919d0, 0.00000000342940918551d0, 
     &    -0.00000035772874310272d0, 0.00000019999935792654d0, 
     &    0.00002687044575042908d0, -0.00011843240273775776d0, 
     &    -0.00080991728956032271d0, 0.00661062970502241174d0, 
     &    0.00909530922354827295d0, -0.20160072778491013140d0, 
     &    0.51169696718727644908d0 / 
      data (b(i), i = 26, 38) / 
     &    0.00000000003147682272d0, -0.00000000048465972408d0, 
     &    0.00000000063675740242d0, 0.00000003377623323271d0, 
     &    -0.00000015451139637086d0, -0.00000203340624738438d0, 
     &    0.00001947204525295057d0, 0.00002854147231653228d0, 
     &    -0.00101565063152200272d0, 0.00271187003520095655d0, 
     &    0.02328095035422810727d0, -0.16725021123116877197d0, 
     &    0.32490054966649436974d0 / 
      data (b(i), i = 39, 51) / 
     &    0.00000000002319363370d0, -0.00000000006303206648d0, 
     &    -0.00000000264888267434d0, 0.00000002050708040581d0, 
     &    0.00000011371857327578d0, -0.00000211211337219663d0, 
     &    0.00000368797328322935d0, 0.00009823686253424796d0, 
     &    -0.00065860243990455368d0, -0.00075285814895230877d0, 
     &    0.02585434424202960464d0, -0.11637092784486193258d0, 
     &    0.18267336775296612024d0 / 
      data (b(i), i = 52, 64) / 
     &    -0.00000000000367789363d0, 0.00000000020876046746d0, 
     &    -0.00000000193319027226d0, -0.00000000435953392472d0, 
     &    0.00000018006992266137d0, -0.00000078441223763969d0, 
     &    -0.00000675407647949153d0, 0.00008428418334440096d0, 
     &    -0.00017604388937031815d0, -0.00239729611435071610d0, 
     &    0.02064129023876022970d0, -0.06905562880005864105d0, 
     &    0.09084526782065478489d0 / 
      w = abs(x)
      if (w .lt. 2.2d0) then
          t = w * w
          k = int(t)
          t = t - k
          k = k * 13
          y = ((((((((((((a(k) * t + a(k + 1)) * t + 
     &        a(k + 2)) * t + a(k + 3)) * t + a(k + 4)) * t + 
     &        a(k + 5)) * t + a(k + 6)) * t + a(k + 7)) * t + 
     &        a(k + 8)) * t + a(k + 9)) * t + a(k + 10)) * t + 
     &        a(k + 11)) * t + a(k + 12)) * w
      else if (w .lt. 6.9d0) then
          k = int(w)
          t = w - k
          k = 13 * (k - 2)
          y = (((((((((((b(k) * t + b(k + 1)) * t + 
     &        b(k + 2)) * t + b(k + 3)) * t + b(k + 4)) * t + 
     &        b(k + 5)) * t + b(k + 6)) * t + b(k + 7)) * t + 
     &        b(k + 8)) * t + b(k + 9)) * t + b(k + 10)) * t + 
     &        b(k + 11)) * t + b(k + 12)
          y = y * y
          y = y * y
          y = y * y
          y = 1 - y * y
      else
          y = 1
      end if
      if (x .lt. 0) y = -y
      erf = y
      end
!
C
       SUBROUTINE OPTICAL(ND,T,PH,P,CLDF,CLDQ
     :                   ,CLDEMI,CLDTAU,CLDFICE,CLDT,CLDWP)
       implicit none
C============================================================================
C Purpose:
C -----------
C   Compute cloud optical properties (emissivity, optical thickness)
C
C CALL command:
C -------------
C   CALL OPTICAL(NA,T2,PH,P,CLDF,CLDQ,CLDEMI,CLDTAU,CLDFICE,CLDT,CLDWP)
C
C Inputs:
C -----------
C  ND-----------: Number of vertical levels
C  T--------ND--: Grid-box average temperature (K)
C  PH-------ND+1: Pressure at interface levels (Pa)
C  P--------ND--: Pressure at mid-levels (Pa)
C  CLDF-----ND--: Cloud fraction at each vertical level (0-1)
C  CLDQ-----ND--: In-cloud condensate mixing ratio (kg/kg)
C
C Outputs:
C -----------
C  CLDEMI---ND--: Cloud longwave emissivity
C  CLDTAU---ND--: Cloud optical thickness
C  CLDFICE--ND--: Ice fraction within the grid-box
C  CLDT------1--: Total cloud cover 
C  CLDWP-----1--: Gridbox averaged cloud water path (kg/m2)
C
C Author & date:
C ---------------
C  Sandrine Bony (MIT & LMD/CNRS)   -    June 1999
C
C  Modification Feb 2000: if newmicro=TRUE, use cloud microphysical properties 
C                         suggested by Iacobellis and Somerville (1999) for
C                         rad_liq and rad_ice, and by Ebert and Curry (1992)
C                         for coef_liq and coef_ice.
C
C Caution: if newmicro=TRUE and if this routine is used in a GCM, information 
C          about the surface (land or ocean) should be used to specify the
C          droplet concentration (Noc or Nland); for the moment we
C          assume that this is an ocean surface.
C   
C============================================================================
C
      integer ND,k,i,n,no
      real T(ND),CLDF(ND),CLDQ(ND),CLDEMI(ND),CLDTAU(ND)
     :    ,PH(ND+1),P(ND),CLDFICE(ND),CLDT,CLDWP
c
c -- equivalent cloud droplet radius for liquid/ice clouds:
      real radius, rad_ice, rad_liq
c      parameter (rad_liq=10.0, rad_ice=30.0) !sb values
      parameter (rad_liq=2.0, rad_ice=30.0)  !KE after sb   1/7/04
c
c -- IR absorption coefficient for liquid/ice clouds (includes
c    the diffusivity factor):
      real coef, coef_ice, coef_liq
      parameter (coef_liq=0.13, coef_ice=0.09)
c
c -- Thresholds for liquid/ice clouds distinction:
      real seuil_neb, T_ice, T_o
      parameter (seuil_neb=0.001)
c from LMD ("CTRL" and "ICE-OPT" expts):
       parameter (T_o=273.15, T_ice=273.15-15.0)
c Zender and Kiehl (95): ("PHASE" expt)
c      parameter (T_o=273.15-10, T_ice=273.15-30.0)
c
c -- for the new microphysical parameterization:
      real Noc, Nland ! droplet nb concentration (cm^{-3})
      parameter (Noc=150., Nland=600.)! Bower et al. 1994
      real k_liq, k_ice0, k_ice, DF, k0,k1,k2,k3
      parameter (k_liq=0.0903, k_ice0=0.005) ! units = m2/g
      parameter (DF=1.66) ! diffusivity factor
      parameter (k0=60.75, k1=-2.47, k2=-0.11, k3=-0.001)
      real pi, aux, rel, rei, kabs, tc
c -- Misc:
      real RG, undef, zflwp, zfiwp
      parameter (RG=9.8,undef=999.999)
      logical lo, newmicro
C      parameter (newmicro=.FALSE.) !if T: new opt prop for ice clds, refnew
      parameter (newmicro=.TRUE.) !if T: new opt prop for ice clds
      pi = ACOS(-1.)
c
c Cloud optical thickness and emissivity:
c
      DO k = 1, ND
         CLDF(k) = MAX(CLDF(k), seuil_neb)
c liquid/ice fraction:
         CLDFICE(k) = 1.0 - (T(k)-T_ice) / (T_o-T_ice)
         CLDFICE(k) = MIN(MAX(CLDFICE(k),0.0),1.0)
c liquid/ice cloud water paths:
         zflwp = 1000.*(1.-CLDFICE(k))*CLDQ(k)
     :          *(PH(k)-PH(k+1))/RG
         zfiwp = 1000.*CLDFICE(k)*CLDQ(k)
     :          *(PH(k)-PH(k+1))/RG
c optical properties:
         IF (newmicro) THEN
c***********************************************************
c NB: "no" means that we modify only the optical properties
c      of ice clouds, not those of liquid clouds
c***********************************************************

c -- parameterization of the effective cloud droplet radius (microns):
c++++ for liquid water clouds: 
             rel = rad_liq ! no
c++++ for ice clouds: as a function of the ambiant temperature
         tc = T(k)-273.15
c ...... formula used by Iacobellis and Somerville (2000):
c ...... (with an asymptotical value of 3.5 microns at T<-81.4 C
c ......  added to be consistent with observations of Heymsfield et al. 1986):
         rei = 0.71*tc + 61.29
         if (tc.le.-81.4) rei = 3.5 ! only micronew, phasenew expts
c cloud optical thickness:
c for liquid clouds, LMD-like, for ice clouds, Ebert & Curry (1992):
c ------------------------------------------------------------------
         if (zflwp.eq.0.) rel = 1. ! no influence
         if (zfiwp.eq.0. .or. rei.le.0.) rei = 1. ! no influence
         CLDTAU(k) = 3.0/2.0 * ( zflwp/rel )
     .             + zfiwp * (3.448e-03  + 2.431/rei)
c cloud infrared emissivity:
c ---------------------------
c the broadband infrared absorption coefficient is parameterized as a
c function of the effective cld droplet radius:
c ... Ebert and Curry (1992) formula:
c as used by Kiehl & Zender (1995):  
         k_ice = k_ice0 + 1.0/rei ! iceopt
c cloud emissivity:
         CLDEMI(k) = 1.0
     .      - EXP( - coef_liq*zflwp - DF*k_ice*zfiwp )
         ELSE
c -- version LMD:
c ---------------
         CLDTAU(k) = 3.0/2.0
     .      * ( zflwp/rad_liq + zfiwp/rad_ice )
         CLDEMI(k) = 1.0
     .      - EXP( - coef_liq*zflwp - coef_ice*zfiwp )
         ENDIF ! newmicro

         lo = (CLDF(k) .LE. seuil_neb)
         IF (lo) CLDF(k) = 0.0
         IF (lo) CLDTAU(k) = 0.0
         IF (lo) CLDEMI(k) = 0.0
         lo = (CLDQ(k) .EQ. undef)
         IF (lo) write(*,*) 'PB: CLDQ EQ UNDEF '
         IF (lo) CLDTAU(k) = undef
         IF (lo) CLDEMI(k) = undef
      ENDDO
C
C Cloud liquid path and total cloudiness:
C
      CLDT  = 1.
      CLDWP = 0.
      no = 0
      DO k = ND, 1, -1
      if (CLDQ(k) .NE. undef .and. CLDF(k).ne.undef) then
       CLDWP = CLDWP
     .   + CLDQ(k)*CLDF(k)*(PH(k)-PH(k+1))/RG
        CLDT = CLDT * (1.0-CLDF(k)) ! random overlap
      else
       no = 1
       write(*,*) 'PB: k, CLDF, CLDQ: ',CLDF(k),CLDQ(k)
      endif
      ENDDO
      CLDT = 1.0 - CLDT ! random overlap only
C
      if (no.eq.1) then
       CLDT  = undef
       CLDWP = undef
      endif
      RETURN
      END
c======================================================================
      SUBROUTINE orbite(xjour,longi,dist)
      IMPLICIT none
c======================================================================
c Auteur(s): Z.X. Li (LMD/CNRS) (adapte du GCM du LMD) date: 19930818
c Objet: pour un jour donne, calculer la longitude vraie de la terre
c        (par rapport au point vernal-21 mars) dans son orbite solaire
c        calculer aussi la distance terre-soleil (unite astronomique)
c======================================================================
c Arguments:
c xjour--INPUT--R- jour de l'annee a compter du 1er janvier
c longi--OUTPUT-R- longitude vraie en degres par rapport au point
c                  vernal (21 mars) en degres
c dist---OUTPUT-R- distance terre-soleil (par rapport a la moyenne)
      REAL xjour, longi, dist
c======================================================================
      include "YOMCST.h"
C
C  -- Variables dynamiques locales
      REAL pir,xl,xllp,xee,xse,xlam,dlamm,anm,ranm,anv,ranv
C
C -- sb:
c call suphec to fill YOMSCT:
       call suphec
C sb --
      pir = 4.0*ATAN(1.0) / 180.0
      xl=R_peri+180.0
      xllp=xl*pir
      xee=R_ecc*R_ecc
      xse=SQRT(1.0-xee)
      xlam = (R_ecc/2.0+R_ecc*xee/8.0)*(1.0+xse)*SIN(xllp)
     .     - xee/4.0*(0.5+xse)*SIN(2.0*xllp)
     .     + R_ecc*xee/8.0*(1.0/3.0+xse)*SIN(3.0*xllp)
      xlam=2.0*xlam/pir
      dlamm=xlam+(xjour-81.0)
      anm=dlamm-xl
      ranm=anm*pir
      xee=xee*R_ecc
      ranv=ranm+(2.0*R_ecc-xee/4.0)*SIN(ranm)
     .         +5.0/4.0*R_ecc*R_ecc*SIN(2.0*ranm)
     .         +13.0/12.0*xee*SIN(3.0*ranm)
c
      anv=ranv/pir
      longi=anv+xl
C
      dist = (1-R_ecc*R_ecc)
     .      /(1+R_ecc*COS(pir*(longi-(R_peri+180.0))))
      RETURN
      END
c======================================================================
      SUBROUTINE angle(longi, lati, frac, muzero)
      IMPLICIT none
c======================================================================
c Auteur(s): Z.X. Li (LMD/CNRS) date: 19930818
c Objet: Calculer la duree d'ensoleillement pour un jour et la hauteur
c        du soleil (cosinus de l'angle zinithal) moyenne sur la journee
c======================================================================
c Arguments:
c longi----INPUT-R- la longitude vraie de la terre dans son plan 
c                   solaire a partir de l'equinoxe de printemps (degre)
c lati-----INPUT-R- la latitude d'un point sur la terre (degre)
c frac-----OUTPUT-R la duree d'ensoleillement dans la journee divisee
c                   par 24 heures (unite en fraction de 0 a 1)
c muzero---OUTPUT-R la moyenne du cosinus de l'angle zinithal sur
c                   la journee (0 a 1)
c======================================================================
      include "dimensions.h"
      include "dimphy.h"
      REAL longi
      REAL lati(klon), frac(klon), muzero(klon)
      include "YOMCST.h"
      REAL lat, omega, lon_sun, lat_sun
      REAL pi_local, incl
      INTEGER i
c
      pi_local = 4.0 * ATAN(1.0)
      incl=R_incl * pi_local / 180.
c
      lon_sun = longi * pi_local / 180.0
      lat_sun = ASIN (sin(lon_sun)*SIN(incl) )
c
      DO i = 1, klon
      lat = lati(i) * pi_local / 180.0
c
      IF ( lat .GE. (pi_local/2.+lat_sun)
     .    .OR. lat.LE.(-pi_local/2.+lat_sun)) THEN
         omega = 0.0   ! nuit polaire
      ELSE IF ( lat.GE.(pi_local/2.-lat_sun)
     .          .OR. lat.LE.(-pi_local/2.-lat_sun)) THEN
         omega = pi_local   ! journee polaire
      ELSE
         omega = -TAN(lat)*TAN(lat_sun)
         omega = ACOS (omega)
      ENDIF
c
      frac(i) = omega / pi_local
c
      IF (omega .GT. 0.0) THEN
         muzero(i) = SIN(lat)*SIN(lat_sun)
     .          + COS(lat)*COS(lat_sun)*SIN(omega) / omega
      ELSE
         muzero(i) = 0.0
      ENDIF
      ENDDO
c
      RETURN
      END
c====================================================================
      SUBROUTINE zenang(longi,gmtime,pdtrad,lat,long,
     s                  pmu0,frac)
      IMPLICIT none
c=============================================================
c Auteur : O. Boucher (LMD/CNRS)
c          d'apres les routines zenith et angle de Z.X. Li 
c Objet  : calculer les valeurs moyennes du cos de l'angle zenithal
c          et l'ensoleillement moyen entre gmtime1 et gmtime2 
c          connaissant la declinaison, la latitude et la longitude.
c Rque   : Different de la routine angle en ce sens que zenang 
c          fournit des moyennes de pmu0 et non des valeurs 
c          instantanees, du coup frac prend toutes les valeurs 
c          entre 0 et 1.
c Date   : premiere version le 13 decembre 1994
c          revu pour  GCM  le 30 septembre 1996
c===============================================================
c longi----INPUT : la longitude vraie de la terre dans son plan
c                  solaire a partir de l'equinoxe de printemps (degre)
c gmtime---INPUT : temps universel en fraction de jour
c pdtrad---INPUT : pas de temps du rayonnement (secondes)
c lat------INPUT : latitude en degres
c long-----INPUT : longitude en degres
c pmu0-----OUTPUT: angle zenithal moyen entre gmtime et gmtime+pdtrad
c frac-----OUTPUT: ensoleillement moyen entre gmtime et gmtime+pdtrad
c================================================================
      include "dimensions.h"
      include "dimphy.h"
      include "YOMCST.h"
c================================================================
      real longi, gmtime, pdtrad
      real lat(klon), long(klon), pmu0(klon), frac(klon)
c================================================================
      integer i
      real gmtime1, gmtime2
      real pi_local, deux_pi_local, incl
      real omega1, omega2, omega
c omega1, omega2 : temps 1 et 2 exprime en radian avec 0 a midi.
c omega : heure en radian du coucher de soleil 
c -omega est donc l'heure en radian de lever du soleil
      real omegadeb, omegafin
      real zfrac1, zfrac2, z1_mu, z2_mu
      real lat_sun          ! declinaison en radian
      real lon_sun          ! longitude solaire en radian
      real latr             ! latitude du pt de grille en radian
c================================================================
c
      pi_local = 4.0 * ATAN(1.0)
      deux_pi_local = 2.0 * pi_local
      incl=R_incl * pi_local / 180.
c
      lon_sun = longi * pi_local / 180.0
      lat_sun = ASIN (SIN(lon_sun)*SIN(incl) )
c
      gmtime1=gmtime*86400.
      gmtime2=gmtime*86400.+pdtrad
c
      DO i = 1, klon
c
      latr = lat(i) * pi_local / 180.
c
c--pose probleme quand lat=+/-90 degres
c
c      omega = -TAN(latr)*TAN(lat_sun)
c      omega = ACOS(omega)
c      IF (latr.GE.(pi_local/2.+lat_sun)
c     .    .OR. latr.LE.(-pi_local/2.+lat_sun)) THEN
c         omega = 0.0       ! nuit polaire
c      ENDIF
c      IF (latr.GE.(pi_local/2.-lat_sun)
c     .          .OR. latr.LE.(-pi_local/2.-lat_sun)) THEN
c         omega = pi_local  ! journee polaire
c      ENDIF
c
c--remplace par cela (le cas par defaut est different)
c
      omega=0.0  !--nuit polaire
      IF (latr.GE.(pi_local/2.-lat_sun)
     .          .OR. latr.LE.(-pi_local/2.-lat_sun)) THEN
         omega = pi_local  ! journee polaire
      ENDIF
      IF (latr.LT.(pi_local/2.+lat_sun).AND.
     .    latr.GT.(-pi_local/2.+lat_sun).AND.
     .    latr.LT.(pi_local/2.-lat_sun).AND.
     .    latr.GT.(-pi_local/2.-lat_sun)) THEN
      omega = -TAN(latr)*TAN(lat_sun)
      omega = ACOS(omega)
      ENDIF
c
         omega1 = gmtime1 + long(i)*86400.0/360.0
         omega1 = omega1 / 86400.0*deux_pi_local
         omega1 = MOD (omega1+deux_pi_local, deux_pi_local)
         omega1 = omega1 - pi_local
c
         omega2 = gmtime2 + long(i)*86400.0/360.0
         omega2 = omega2 / 86400.0*deux_pi_local
         omega2 = MOD (omega2+deux_pi_local, deux_pi_local)
         omega2 = omega2 - pi_local
c
      IF (omega1.LE.omega2) THEN  !--on est dans la meme journee locale
c
      IF (omega2.LE.-omega .OR. omega1.GE.omega
     .                     .OR. omega.LT.1e-5) THEN   !--nuit
         frac(i)=0.0
         pmu0(i)=0.0
      ELSE                                              !--jour+nuit/jour
        omegadeb=MAX(-omega,omega1)
        omegafin=MIN(omega,omega2)
        frac(i)=(omegafin-omegadeb)/(omega2-omega1)
        pmu0(i)=SIN(latr)*SIN(lat_sun) + 
     .          COS(latr)*COS(lat_sun)*
     .          (SIN(omegafin)-SIN(omegadeb))/
     .          (omegafin-omegadeb)        
      ENDIF
c
      ELSE  !---omega1 GT omega2 -- a cheval sur deux journees
c
c-------------------entre omega1 et pi
      IF (omega1.GE.omega) THEN  !--nuit
         zfrac1=0.0
         z1_mu =0.0
      ELSE                       !--jour+nuit
        omegadeb=MAX(-omega,omega1)
        omegafin=omega
        zfrac1=omegafin-omegadeb
        z1_mu =SIN(latr)*SIN(lat_sun) +
     .          COS(latr)*COS(lat_sun)*
     .          (SIN(omegafin)-SIN(omegadeb))/
     .          (omegafin-omegadeb)
      ENDIF 
c---------------------entre -pi et omega2
      IF (omega2.LE.-omega) THEN   !--nuit
         zfrac2=0.0
         z2_mu =0.0
      ELSE                         !--jour+nuit
         omegadeb=-omega
         omegafin=MIN(omega,omega2)
         zfrac2=omegafin-omegadeb
         z2_mu =SIN(latr)*SIN(lat_sun) +
     .           COS(latr)*COS(lat_sun)*
     .           (SIN(omegafin)-SIN(omegadeb))/
     .           (omegafin-omegadeb)
c
      ENDIF
c-----------------------moyenne 
      frac(i)=(zfrac1+zfrac2)/(omega2+deux_pi_local-omega1)
      pmu0(i)=(zfrac1*z1_mu+zfrac2*z2_mu)/MAX(zfrac1+zfrac2,1.E-10)
c
      ENDIF   !---comparaison omega1 et omega2
c
      ENDDO
c
      END
c===================================================================
      SUBROUTINE zenith (longi, gmtime, lat, long,
     s                   pmu0, fract)
      IMPLICIT none
c
c Auteur(s): Z.X. Li (LMD/ENS)
c
c Objet: calculer le cosinus de l'angle zenithal du soleil en
c        connaissant la declinaison du soleil, la latitude et la
c        longitude du point sur la terre, et le temps universel
c
c Arguments d'entree:
c     longi  : declinaison du soleil (en degres)
c     gmtime : temps universel en second qui varie entre 0 et 86400
c     lat    : latitude en degres
c     long   : longitude en degres
c Arguments de sortie:
c     pmu0   : cosinus de l'angle zenithal
c
c====================================================================
      include "dimensions.h"
      include "dimphy.h"
      include "YOMCST.h"
c====================================================================
      REAL longi, gmtime
      REAL lat(klon), long(klon), pmu0(klon), fract(klon)
c=====================================================================
      INTEGER n
      REAL zpi, zpir, omega, zgmtime
      REAL incl, lat_sun, lon_sun
c----------------------------------------------------------------------
      zpi = 4.0*ATAN(1.0)
      zpir = zpi / 180.0
      zgmtime=gmtime*86400.
c
      incl=R_incl * zpir
c
      lon_sun = longi * zpir
      lat_sun = ASIN (SIN(lon_sun)*SIN(incl) )
c
c--initialisation a la nuit
c
      DO n =1, klon
        pmu0(n)=0.
        fract(n)=0.0
      ENDDO
c
c 1 degre en longitude = 240 secondes en temps
c
      DO n = 1, klon
         omega = zgmtime + long(n)*86400.0/360.0
         omega = omega / 86400.0 * 2.0 * zpi
         omega = MOD(omega + 2.0 * zpi, 2.0 * zpi)
         omega = omega - zpi
         pmu0(n) = sin(lat(n)*zpir) * sin(lat_sun)
     .           + cos(lat(n)*zpir) * cos(lat_sun)
     .           * cos(omega)
         pmu0(n) = MAX (pmu0(n), 0.0)
         IF (pmu0(n).GT.1.E-6) fract(n)=1.0
      ENDDO
c
      RETURN
      END
      SUBROUTINE radlwsw(dist, rmu0, fract, co2_ppm, solaire,
     .                  paprs, pplay,tsol,albedo, t,q,wo,
     .                  cldfra, cldemi, cldtau,
     .                  heat,heat0,cool,cool0,radsol,albpla,
     .                  topsw,toplw,solsw,sollw,
     .                  topsw0,toplw0,solsw0,sollw0)
      IMPLICIT none
c======================================================================
c Auteur(s): Z.X. Li (LMD/CNRS) date: 19960719
c Objet: interface entre le modele et les rayonnements
c Arguments:
c dist-----input-R- distance astronomique terre-soleil
c rmu0-----input-R- cosinus de l'angle zenithal
c fract----input-R- duree d'ensoleillement normalisee
c co2_ppm--input-R- concentration du gaz carbonique (en ppm)
c solaire--input-R- constante solaire (W/m**2)
c paprs----input-R- pression a inter-couche (Pa)
c pplay----input-R- pression au milieu de couche (Pa)
c tsol-----input-R- temperature du sol (en K)
c albedo---input-R- albedo du sol (entre 0 et 1)
c t--------input-R- temperature (K)
c q--------input-R- vapeur d'eau (en kg/kg)
c wo-------input-R- contenu en ozone (en kg/kg)
c cldfra---input-R- fraction nuageuse (entre 0 et 1)
c cldtau---input-R- epaisseur optique des nuages dans le visible
c cldemi---input-R- emissivite des nuages dans l'IR (entre 0 et 1)
c
c heat-----output-R- echauffement atmospherique (visible) (K/jour)
c cool-----output-R- refroidissement dans l'IR (K/jour)
c radsol---output-R- bilan radiatif net au sol (W/m**2) (+ vers le bas)
c albpla---output-R- albedo planetaire (entre 0 et 1)
c topsw----output-R- ray. solaire absorbe au sommet de l'atm.
c toplw----output-R- ray. IR emis au sommet de l'atmosphere
c solsw----output-R- ray. solaire net absorbe a la surface
c sollw----output-R- ray. IR net emis par la surface
c======================================================================
      include "dimensions.h"
      include "dimphy.h"
      include "raddim.h"
c
      real rmu0(klon), fract(klon), dist
      real co2_ppm
      real solaire
c
      real paprs(klon,klev+1), pplay(klon,klev)
      real albedo(klon), tsol(klon)
      real t(klon,klev), q(klon,klev), wo(klon,klev)
      real cldfra(klon,klev), cldemi(klon,klev), cldtau(klon,klev)
      real heat(klon,klev), cool(klon,klev)
      real heat0(klon,klev), cool0(klon,klev)
      real radsol(klon), topsw(klon), toplw(klon)
      real solsw(klon), sollw(klon), albpla(klon)
      real topsw0(klon), toplw0(klon), solsw0(klon), sollw0(klon)
c
      REAL*8 zx_alpha1, zx_alpha2
c
      include "YOMCST.h"
c
      INTEGER k, kk, i, j, iof, nb_gr
      EXTERNAL lw, sw
c
      REAL*8 RCO2, RCH4, RN2O, RCFC11, RCFC12
      REAL*8 PSCT
c
      REAL*8 PALBD(kdlon,2), PALBP(kdlon,2)
      REAL*8 PEMIS(kdlon), PDT0(kdlon), PVIEW(kdlon)
      REAL*8 PPSOL(kdlon), PDP(kdlon,klev)
      REAL*8 PTL(kdlon,kflev+1), PPMB(kdlon,kflev+1)
      REAL*8 PTAVE(kdlon,kflev)
      REAL*8 PWV(kdlon,kflev), PQS(kdlon,kflev), POZON(kdlon,kflev)
      REAL*8 PAER(kdlon,kflev,5)
      REAL*8 PCLDLD(kdlon,kflev)
      REAL*8 PCLDLU(kdlon,kflev)
      REAL*8 PCLDSW(kdlon,kflev)
      REAL*8 PTAU(kdlon,2,kflev)
      REAL*8 POMEGA(kdlon,2,kflev)
      REAL*8 PCG(kdlon,2,kflev)
c
      REAL*8 zfract(kdlon), zrmu0(kdlon), zdist
c
      REAL*8 zheat(kdlon,kflev), zcool(kdlon,kflev)
      REAL*8 zheat0(kdlon,kflev), zcool0(kdlon,kflev)
      REAL*8 ztopsw(kdlon), ztoplw(kdlon)
      REAL*8 zsolsw(kdlon), zsollw(kdlon), zalbpla(kdlon)
      REAL*8 ztopsw0(kdlon), ztoplw0(kdlon)
      REAL*8 zsolsw0(kdlon), zsollw0(kdlon)
c
c-------------------------------------------
      nb_gr = klon / kdlon
      IF (nb_gr*kdlon .NE. klon) THEN
         PRINT*, "kdlon mauvais:", klon, kdlon, nb_gr
         CALL abort
      ENDIF
      IF (kflev .NE. klev) THEN
          PRINT*, "kflev differe de klev, kflev, klev"
          CALL abort
      ENDIF
c-------------------------------------------
      DO k = 1, klev
      DO i = 1, klon
         heat(i,k)=0.
         cool(i,k)=0.
         heat0(i,k)=0.
         cool0(i,k)=0.
      ENDDO
      ENDDO
c
      zdist = dist
c
      RCO2 = co2_ppm * 1.0e-06  * 44.011/28.97
      RCH4 = 1.72E-06* 16.043/28.97
      RN2O = 310.E-09* 44.013/28.97
      RCFC11 = 280.E-12* 137.3686/28.97
      RCFC12 = 484.E-12* 120.9140/28.97
      PSCT = solaire/zdist/zdist
c
      DO 99999 j = 1, nb_gr
      iof = kdlon*(j-1)
c
      DO i = 1, kdlon
         zfract(i) = fract(iof+i)
         zrmu0(i) = rmu0(iof+i)
         PALBD(i,1) = albedo(iof+i)
         PALBD(i,2) = albedo(iof+i)
         PALBP(i,1) = albedo(iof+i)
         PALBP(i,2) = albedo(iof+i)
ccc SB jul99         PEMIS(i) = 0.96
         PEMIS(i) = 0.999
         PVIEW(i) = 1.66
         PPSOL(i) = paprs(iof+i,1)
         zx_alpha1 = (paprs(iof+i,1)-pplay(iof+i,2)) 
     .             / (pplay(iof+i,1)-pplay(iof+i,2))
         zx_alpha2 = 1.0 - zx_alpha1
         PTL(i,1) = t(iof+i,1) * zx_alpha1 + t(iof+i,2) * zx_alpha2
         PTL(i,klev+1) = t(iof+i,klev)
         PDT0(i) = tsol(iof+i) - PTL(i,1)
      ENDDO
      DO k = 2, kflev
      DO i = 1, kdlon
         PTL(i,k) = (t(iof+i,k)+t(iof+i,k-1))*0.5
      ENDDO
      ENDDO
      DO k = 1, kflev
      DO i = 1, kdlon
         PDP(i,k) = paprs(iof+i,k)-paprs(iof+i,k+1)
         PTAVE(i,k) = t(iof+i,k)
         PWV(i,k) = MAX (q(iof+i,k), 1.0e-12)
         PQS(i,k) = PWV(i,k)
c wo:    cm.atm (epaisseur en cm dans la situation standard)
c POZON: kg/kg
c! sb        POZON(i,k) = MAX(wo(iof+i,k),1.0e-12)*RG/46.6968
c! sb    .               /(paprs(iof+i,k)-paprs(iof+i,k+1))
c! sb     .               *(paprs(iof+i,1)/101325.0)
         POZON(i,k) = wo(iof+i,k)
         PCLDLD(i,k) = cldfra(iof+i,k)*cldemi(iof+i,k)
         PCLDLU(i,k) = cldfra(iof+i,k)*cldemi(iof+i,k)
         PCLDSW(i,k) = cldfra(iof+i,k)
         PTAU(i,1,k) = MAX(cldtau(iof+i,k), 1.0e-05)! 1e-12 serait instable
         PTAU(i,2,k) = MAX(cldtau(iof+i,k), 1.0e-05)! pour 32-bit machines
         POMEGA(i,1,k) = 0.9999 - 5.0e-04 * EXP(-0.5 * PTAU(i,1,k))
         POMEGA(i,2,k) = 0.9988 - 2.5e-03 * EXP(-0.05 * PTAU(i,2,k))
         PCG(i,1,k) = 0.865
         PCG(i,2,k) = 0.910
      ENDDO
      ENDDO
c
      DO k = 1, kflev+1
      DO i = 1, kdlon
         PPMB(i,k) = paprs(iof+i,k)/100.0
      ENDDO
      ENDDO
c
      DO kk = 1, 5
      DO k = 1, kflev
      DO i = 1, kdlon
         PAER(i,k,kk) = 1.0E-15
      ENDDO
      ENDDO
      ENDDO
c
c======================================================================
      CALL LW(RCO2,RCH4,RN2O,RCFC11,RCFC12,
     .        PPMB, PDP,
     .        PPSOL,PDT0,PEMIS,
     .        PTL, PTAVE, PWV, POZON, PAER,
     .        PCLDLD,PCLDLU,
     .        PVIEW,
     .        zcool, zcool0,
     .        ztoplw,zsollw,ztoplw0,zsollw0)
      CALL SW(PSCT, RCO2, zrmu0, zfract,
     S        PPMB, PDP,
     S        PPSOL, PALBD, PALBP,
     S        PTAVE, PWV, PQS, POZON, PAER,
     S        PCLDSW, PTAU, POMEGA, PCG,
     S        zheat, zheat0,
     S        zalbpla,ztopsw,zsolsw,ztopsw0,zsolsw0)
c======================================================================
      DO i = 1, kdlon
         radsol(iof+i) = zsolsw(i) + zsollw(i)
         topsw(iof+i) = ztopsw(i)
         toplw(iof+i) = ztoplw(i)
         solsw(iof+i) = zsolsw(i)
         sollw(iof+i) = zsollw(i)
         topsw0(iof+i) = ztopsw0(i)
         toplw0(iof+i) = ztoplw0(i)
         solsw0(iof+i) = zsolsw0(i)
         sollw0(iof+i) = zsollw0(i)
         albpla(iof+i) = zalbpla(i)
      ENDDO
      DO k = 1, kflev
      DO i = 1, kdlon
         heat(iof+i,k) = zheat(i,k)
         cool(iof+i,k) = zcool(i,k)
         heat0(iof+i,k) = zheat0(i,k)
         cool0(iof+i,k) = zcool0(i,k)
      ENDDO
      ENDDO
c
99999 CONTINUE
      RETURN
      END
      SUBROUTINE SW(PSCT, RCO2, PRMU0, PFRAC, 
     S              PPMB, PDP, 
     S              PPSOL, PALBD, PALBP,
     S              PTAVE, PWV, PQS, POZON, PAER,
     S              PCLDSW, PTAU, POMEGA, PCG,
     S              PHEAT, PHEAT0,
     S              PALBPLA,PTOPSW,PSOLSW,PTOPSW0,PSOLSW0)
      IMPLICIT none
      include "dimensions.h"
      include "dimphy.h"
      include "raddim.h"
      include "YOMCST.h"
C
C     ------------------------------------------------------------------
C
C     PURPOSE.
C     --------
C
C          THIS ROUTINE COMPUTES THE SHORTWAVE RADIATION FLUXES IN TWO
C     SPECTRAL INTERVALS FOLLOWING FOUQUART AND BONNEL (1980).
C
C     METHOD.
C     -------
C
C          1. COMPUTES ABSORBER AMOUNTS                 (SWU)
C          2. COMPUTES FLUXES IN 1ST SPECTRAL INTERVAL  (SW1S)
C          3. COMPUTES FLUXES IN 2ND SPECTRAL INTERVAL  (SW2S)
C
C     REFERENCE.
C     ----------
C
C        SEE RADIATION'S PART OF THE ECMWF RESEARCH DEPARTMENT
C        DOCUMENTATION, AND FOUQUART AND BONNEL (1980)
C
C     AUTHOR.
C     -------
C        JEAN-JACQUES MORCRETTE  *ECMWF*
C
C     MODIFICATIONS.
C     --------------
C        ORIGINAL : 89-07-14
C        95-01-01   J.-J. MORCRETTE  Direct/Diffuse Albedo
C     ------------------------------------------------------------------
C
C* ARGUMENTS:
C
      REAL*8 PSCT  ! constante solaire (valeur conseillee: 1370)
      REAL*8 RCO2  ! concentration CO2 (IPCC: 353.E-06*44.011/28.97)
C
      REAL*8 PPSOL(KDLON)        ! SURFACE PRESSURE (PA)
      REAL*8 PDP(KDLON,KFLEV)    ! LAYER THICKNESS (PA)
      REAL*8 PPMB(KDLON,KFLEV+1) ! HALF-LEVEL PRESSURE (MB)
C
      REAL*8 PRMU0(KDLON)  ! COSINE OF ZENITHAL ANGLE
      REAL*8 PFRAC(KDLON)  ! fraction de la journee
C
      REAL*8 PTAVE(KDLON,KFLEV)  ! LAYER TEMPERATURE (K)
      REAL*8 PWV(KDLON,KFLEV)    ! SPECIFIC HUMIDITY (KG/KG)
      REAL*8 PQS(KDLON,KFLEV)    ! SATURATED WATER VAPOUR (KG/KG)
      REAL*8 POZON(KDLON,KFLEV)  ! OZONE CONCENTRATION (KG/KG)
      REAL*8 PAER(KDLON,KFLEV,5) ! AEROSOLS' OPTICAL THICKNESS
C
      REAL*8 PALBD(KDLON,2)  ! albedo du sol (lumiere diffuse)
      REAL*8 PALBP(KDLON,2)  ! albedo du sol (lumiere parallele)
C
      REAL*8 PCLDSW(KDLON,KFLEV)    ! CLOUD FRACTION
      REAL*8 PTAU(KDLON,2,KFLEV)    ! CLOUD OPTICAL THICKNESS
      REAL*8 PCG(KDLON,2,KFLEV)     ! ASYMETRY FACTOR
      REAL*8 POMEGA(KDLON,2,KFLEV)  ! SINGLE SCATTERING ALBEDO
C
      REAL*8 PHEAT(KDLON,KFLEV) ! SHORTWAVE HEATING (K/DAY)
      REAL*8 PHEAT0(KDLON,KFLEV)! SHORTWAVE HEATING (K/DAY) clear-sky
      REAL*8 PALBPLA(KDLON)     ! PLANETARY ALBEDO
      REAL*8 PTOPSW(KDLON)      ! SHORTWAVE FLUX AT T.O.A.
      REAL*8 PSOLSW(KDLON)      ! SHORTWAVE FLUX AT SURFACE
      REAL*8 PTOPSW0(KDLON)     ! SHORTWAVE FLUX AT T.O.A. (CLEAR-SKY)
      REAL*8 PSOLSW0(KDLON)     ! SHORTWAVE FLUX AT SURFACE (CLEAR-SKY)
C
C* LOCAL VARIABLES:
C
      REAL*8 ZOZ(KDLON,KFLEV)
      REAL*8 ZAKI(KDLON,2)     
      REAL*8 ZCLD(KDLON,KFLEV)
      REAL*8 ZCLEAR(KDLON) 
      REAL*8 ZDSIG(KDLON,KFLEV)
      REAL*8 ZFACT(KDLON)
      REAL*8 ZFD(KDLON,KFLEV+1)
      REAL*8 ZFDOWN(KDLON,KFLEV+1)
      REAL*8 ZFU(KDLON,KFLEV+1)
      REAL*8 ZFUP(KDLON,KFLEV+1)
      REAL*8 ZRMU(KDLON)
      REAL*8 ZSEC(KDLON)
      REAL*8 ZUD(KDLON,5,KFLEV+1)
      REAL*8 ZCLDSW0(KDLON,KFLEV)
c
      REAL*8 ZFSUP(KDLON,KFLEV+1)
      REAL*8 ZFSDN(KDLON,KFLEV+1)
      REAL*8 ZFSUP0(KDLON,KFLEV+1)
      REAL*8 ZFSDN0(KDLON,KFLEV+1)
      SAVE ZFSUP, ZFSDN, ZFSUP0, ZFSDN0
C
      INTEGER inu, jl, jk, i, k, kpl1
c
      INTEGER swpas  ! Every swpas steps, sw is calculated
      PARAMETER(swpas=1)
c
      INTEGER itapsw
      LOGICAL appel1er
      DATA itapsw /0/
      DATA appel1er /.TRUE./
c
      IF (appel1er) THEN
         PRINT*, 'SW calling frequency : ', swpas
         PRINT*, "   In general, it should be 1"
         appel1er = .FALSE.
      ENDIF
C     ------------------------------------------------------------------
      IF (MOD(itapsw,swpas).EQ.0) THEN
c
      DO JK = 1 , KFLEV
      DO JL = 1, KDLON
         ZCLDSW0(JL,JK) = 0.0
         ZOZ(JL,JK) = POZON(JL,JK)*46.6968/RG
     .               *PDP(JL,JK)*(101325.0/PPSOL(JL))
      ENDDO
      ENDDO
C
C
c clear-sky:
      CALL SWU(PSCT,RCO2,ZCLDSW0,PPMB,PPSOL,
     S         PRMU0,PFRAC,PTAVE,PWV,
     S         ZAKI,ZCLD,ZCLEAR,ZDSIG,ZFACT,ZRMU,ZSEC,ZUD)
      INU = 1
      CALL SW1S(INU,
     S     PAER, PALBD, PALBP, PCG, ZCLD, ZCLEAR, ZCLDSW0,
     S     ZDSIG, POMEGA, ZOZ, ZRMU, ZSEC, PTAU, ZUD,
     S     ZFD, ZFU)
      INU = 2
      CALL SW2S(INU,
     S     PAER, ZAKI, PALBD, PALBP, PCG, ZCLD, ZCLEAR, ZCLDSW0,
     S     ZDSIG, POMEGA, ZOZ, ZRMU, ZSEC, PTAU, ZUD,
     S     PWV, PQS,
     S     ZFDOWN, ZFUP)
      DO JK = 1 , KFLEV+1
      DO JL = 1, KDLON
         ZFSUP0(JL,JK) = (ZFUP(JL,JK)   + ZFU(JL,JK)) * ZFACT(JL)
         ZFSDN0(JL,JK) = (ZFDOWN(JL,JK) + ZFD(JL,JK)) * ZFACT(JL)
      ENDDO
      ENDDO
c cloudy-sky:
      CALL SWU(PSCT,RCO2,PCLDSW,PPMB,PPSOL,
     S         PRMU0,PFRAC,PTAVE,PWV,
     S         ZAKI,ZCLD,ZCLEAR,ZDSIG,ZFACT,ZRMU,ZSEC,ZUD)
      INU = 1
      CALL SW1S(INU,
     S     PAER, PALBD, PALBP, PCG, ZCLD, ZCLEAR, PCLDSW,
     S     ZDSIG, POMEGA, ZOZ, ZRMU, ZSEC, PTAU, ZUD,
     S     ZFD, ZFU)
      INU = 2
      CALL SW2S(INU,
     S     PAER, ZAKI, PALBD, PALBP, PCG, ZCLD, ZCLEAR, PCLDSW,
     S     ZDSIG, POMEGA, ZOZ, ZRMU, ZSEC, PTAU, ZUD,
     S     PWV, PQS,
     S    ZFDOWN, ZFUP)
      DO JK = 1 , KFLEV+1
      DO JL = 1, KDLON
         ZFSUP(JL,JK) = (ZFUP(JL,JK)   + ZFU(JL,JK)) * ZFACT(JL)
         ZFSDN(JL,JK) = (ZFDOWN(JL,JK) + ZFD(JL,JK)) * ZFACT(JL)
      ENDDO
      ENDDO
c
      itapsw = 0
      ENDIF
      itapsw = itapsw + 1
C
      DO k = 1, KFLEV
         kpl1 = k+1
         DO i = 1, KDLON
            PHEAT(i,k) = -(ZFSUP(i,kpl1)-ZFSUP(i,k))
     .                     -(ZFSDN(i,k)-ZFSDN(i,kpl1))
            PHEAT(i,k) = PHEAT(i,k) * RDAY*RG/RCPD / PDP(i,k)
            PHEAT0(i,k) = -(ZFSUP0(i,kpl1)-ZFSUP0(i,k))
     .                     -(ZFSDN0(i,k)-ZFSDN0(i,kpl1))
            PHEAT0(i,k) = PHEAT0(i,k) * RDAY*RG/RCPD / PDP(i,k)
         ENDDO
      ENDDO
      DO i = 1, KDLON
         PALBPLA(i) = ZFSUP(i,KFLEV+1)/(ZFSDN(i,KFLEV+1)+1.0e-20)
c
         PSOLSW(i) = ZFSDN(i,1) - ZFSUP(i,1)
         PTOPSW(i) = ZFSDN(i,KFLEV+1) - ZFSUP(i,KFLEV+1)
c
         PSOLSW0(i) = ZFSDN0(i,1) - ZFSUP0(i,1)
         PTOPSW0(i) = ZFSDN0(i,KFLEV+1) - ZFSUP0(i,KFLEV+1)
      ENDDO
C
      RETURN
      END
c
      SUBROUTINE SWU (PSCT,RCO2,PCLDSW,PPMB,PPSOL,PRMU0,PFRAC,
     S                PTAVE,PWV,PAKI,PCLD,PCLEAR,PDSIG,PFACT,
     S                PRMU,PSEC,PUD)
      IMPLICIT none
      include "dimensions.h"
      include "dimphy.h"
      include "raddim.h"
      include "radepsi.h"
      include "radopt.h"
      include "YOMCST.h"
C
C* ARGUMENTS:
C
      REAL*8 PSCT
      REAL*8 RCO2
      REAL*8 PCLDSW(KDLON,KFLEV)
      REAL*8 PPMB(KDLON,KFLEV+1)
      REAL*8 PPSOL(KDLON)
      REAL*8 PRMU0(KDLON)
      REAL*8 PFRAC(KDLON)
      REAL*8 PTAVE(KDLON,KFLEV)
      REAL*8 PWV(KDLON,KFLEV)
C
      REAL*8 PAKI(KDLON,2)
      REAL*8 PCLD(KDLON,KFLEV)
      REAL*8 PCLEAR(KDLON)
      REAL*8 PDSIG(KDLON,KFLEV)
      REAL*8 PFACT(KDLON)
      REAL*8 PRMU(KDLON)
      REAL*8 PSEC(KDLON)
      REAL*8 PUD(KDLON,5,KFLEV+1)
C
C* LOCAL VARIABLES:
C
      INTEGER IIND(2)
      REAL*8 ZC1J(KDLON,KFLEV+1)
      REAL*8 ZCLEAR(KDLON)
      REAL*8 ZCLOUD(KDLON)
      REAL*8 ZN175(KDLON)
      REAL*8 ZN190(KDLON)
      REAL*8 ZO175(KDLON)
      REAL*8 ZO190(KDLON)
      REAL*8 ZSIGN(KDLON)
      REAL*8 ZR(KDLON,2) 
      REAL*8 ZSIGO(KDLON)
      REAL*8 ZUD(KDLON,2)
      REAL*8 ZRTH, ZRTU, ZWH2O, ZDSCO2, ZDSH2O, ZFPPW
      INTEGER jl, jk, jkp1, jkl, jklp1, ja
C
C* Prescribed Data:
c
      REAL*8 ZPDH2O,ZPDUMG
      SAVE ZPDH2O,ZPDUMG
      REAL*8 ZPRH2O,ZPRUMG
      SAVE ZPRH2O,ZPRUMG
      REAL*8 RTDH2O,RTDUMG
      SAVE RTDH2O,RTDUMG
      REAL*8 RTH2O ,RTUMG
      SAVE RTH2O ,RTUMG
      DATA ZPDH2O,ZPDUMG / 0.8   , 0.75 /
      DATA ZPRH2O,ZPRUMG / 30000., 30000. /
      DATA RTDH2O,RTDUMG /  0.40  , 0.375 /
      DATA RTH2O ,RTUMG  /  240.  , 240.  /
C     ------------------------------------------------------------------
C
C*         1.     COMPUTES AMOUNTS OF ABSORBERS
C                 -----------------------------
C
 100  CONTINUE
C
      IIND(1)=1
      IIND(2)=2
C      
C
C*         1.1    INITIALIZES QUANTITIES
C                 ----------------------
C
 110  CONTINUE
C
      DO 111 JL = 1, KDLON
      PUD(JL,1,KFLEV+1)=0.
      PUD(JL,2,KFLEV+1)=0.
      PUD(JL,3,KFLEV+1)=0.
      PUD(JL,4,KFLEV+1)=0.
      PUD(JL,5,KFLEV+1)=0.
      PFACT(JL)= PRMU0(JL) * PFRAC(JL) * PSCT
      PRMU(JL)=SQRT(1224.* PRMU0(JL) * PRMU0(JL) + 1.) / 35.
      PSEC(JL)=1./PRMU(JL)
      ZC1J(JL,KFLEV+1)=0.
 111  CONTINUE
C
C*          1.3    AMOUNTS OF ABSORBERS
C                  --------------------
C
 130  CONTINUE
C
      DO 131 JL= 1, KDLON
      ZUD(JL,1) = 0.
      ZUD(JL,2) = 0.
      ZO175(JL) = PPSOL(JL)** (ZPDUMG+1.)
      ZO190(JL) = PPSOL(JL)** (ZPDH2O+1.)
      ZSIGO(JL) = PPSOL(JL)
      ZCLEAR(JL)=1.
      ZCLOUD(JL)=0.
 131  CONTINUE
C
      DO 133 JK = 1 , KFLEV
      JKP1 = JK + 1
      JKL = KFLEV+1 - JK
      JKLP1 = JKL+1
      DO 132 JL = 1, KDLON
      ZRTH=(RTH2O/PTAVE(JL,JK))**RTDH2O
      ZRTU=(RTUMG/PTAVE(JL,JK))**RTDUMG
      ZWH2O = MAX (PWV(JL,JK) , ZEPSCQ )
      ZSIGN(JL) = 100. * PPMB(JL,JKP1)
      PDSIG(JL,JK) = (ZSIGO(JL) - ZSIGN(JL))/PPSOL(JL)
      ZN175(JL) = ZSIGN(JL) ** (ZPDUMG+1.)
      ZN190(JL) = ZSIGN(JL) ** (ZPDH2O+1.)
      ZDSCO2 = ZO175(JL) - ZN175(JL)
      ZDSH2O = ZO190(JL) - ZN190(JL)
      PUD(JL,1,JK) = 1./( 10.* RG * (ZPDH2O+1.) )/(ZPRH2O**ZPDH2O)
     .             * ZDSH2O * ZWH2O  * ZRTH
      PUD(JL,2,JK) = 1./( 10.* RG * (ZPDUMG+1.) )/(ZPRUMG**ZPDUMG)
     .             * ZDSCO2 * RCO2 * ZRTU
      ZFPPW=1.6078*ZWH2O/(1.+0.608*ZWH2O)
      PUD(JL,4,JK)=PUD(JL,1,JK)*ZFPPW
      PUD(JL,5,JK)=PUD(JL,1,JK)*(1.-ZFPPW)
      ZUD(JL,1) = ZUD(JL,1) + PUD(JL,1,JK)
      ZUD(JL,2) = ZUD(JL,2) + PUD(JL,2,JK)
      ZSIGO(JL) = ZSIGN(JL)
      ZO175(JL) = ZN175(JL)
      ZO190(JL) = ZN190(JL)
C      
      IF (NOVLP.EQ.1) THEN
         ZCLEAR(JL)=ZCLEAR(JL)
     S               *(1.-MAX(PCLDSW(JL,JKL),ZCLOUD(JL)))
     S               /(1.-MIN(ZCLOUD(JL),1.-ZEPSEC))
         ZC1J(JL,JKL)= 1.0 - ZCLEAR(JL)
         ZCLOUD(JL) = PCLDSW(JL,JKL)
      ELSE IF (NOVLP.EQ.2) THEN
         ZCLOUD(JL) = MAX(PCLDSW(JL,JKL),ZCLOUD(JL))
         ZC1J(JL,JKL) = ZCLOUD(JL)
      ELSE IF (NOVLP.EQ.3) THEN
         ZCLEAR(JL) = ZCLEAR(JL)*(1.-PCLDSW(JL,JKL))
         ZCLOUD(JL) = 1.0 - ZCLEAR(JL)
         ZC1J(JL,JKL) = ZCLOUD(JL)
      END IF
 132  CONTINUE
 133  CONTINUE
      DO 134 JL=1, KDLON
      PCLEAR(JL)=1.-ZC1J(JL,1)
 134  CONTINUE
      DO 136 JK=1,KFLEV
      DO 135 JL=1, KDLON
      IF (PCLEAR(JL).LT.1.) THEN
         PCLD(JL,JK)=PCLDSW(JL,JK)/(1.-PCLEAR(JL))
      ELSE
         PCLD(JL,JK)=0.
      END IF
 135  CONTINUE
 136  CONTINUE           
C      
C
C*         1.4    COMPUTES CLEAR-SKY GREY ABSORPTION COEFFICIENTS
C                 -----------------------------------------------
C
 140  CONTINUE
C
      DO 142 JA = 1,2
      DO 141 JL = 1, KDLON
      ZUD(JL,JA) = ZUD(JL,JA) * PSEC(JL)
 141  CONTINUE
 142  CONTINUE
C
      CALL SWTT1(2, 2, IIND, ZUD, ZR)
C
      DO 144 JA = 1,2
      DO 143 JL = 1, KDLON
      PAKI(JL,JA) = -LOG( ZR(JL,JA) ) / ZUD(JL,JA)
 143  CONTINUE
 144  CONTINUE
C
C
C     ------------------------------------------------------------------
C
      RETURN
      END
      SUBROUTINE SW1S ( KNU
     S  ,  PAER  , PALBD , PALBP, PCG  , PCLD , PCLEAR, PCLDSW
     S  ,  PDSIG , POMEGA, POZ  , PRMU , PSEC , PTAU  , PUD  
     S  ,  PFD   , PFU)
      IMPLICIT none
      include "dimensions.h"
      include "dimphy.h"
      include "raddim.h"
C
C     ------------------------------------------------------------------
C     PURPOSE.
C     --------
C
C          THIS ROUTINE COMPUTES THE SHORTWAVE RADIATION FLUXES IN TWO
C     SPECTRAL INTERVALS FOLLOWING FOUQUART AND BONNEL (1980).
C
C     METHOD.
C     -------
C
C          1. COMPUTES UPWARD AND DOWNWARD FLUXES CORRESPONDING TO
C     CONTINUUM SCATTERING
C          2. MULTIPLY BY OZONE TRANSMISSION FUNCTION
C
C     REFERENCE.
C     ----------
C
C        SEE RADIATION'S PART OF THE ECMWF RESEARCH DEPARTMENT
C        DOCUMENTATION, AND FOUQUART AND BONNEL (1980)
C
C     AUTHOR.
C     -------
C        JEAN-JACQUES MORCRETTE  *ECMWF*
C
C     MODIFICATIONS.
C     --------------
C        ORIGINAL : 89-07-14
C        94-11-15   J.-J. MORCRETTE    DIRECT/DIFFUSE ALBEDO
C     ------------------------------------------------------------------
C
C* ARGUMENTS:
C
      INTEGER KNU
      REAL*8 PAER(KDLON,KFLEV,5)
      REAL*8 PALBD(KDLON,2)
      REAL*8 PALBP(KDLON,2)
      REAL*8 PCG(KDLON,2,KFLEV)  
      REAL*8 PCLD(KDLON,KFLEV)
      REAL*8 PCLDSW(KDLON,KFLEV)
      REAL*8 PCLEAR(KDLON)
      REAL*8 PDSIG(KDLON,KFLEV)
      REAL*8 POMEGA(KDLON,2,KFLEV)
      REAL*8 POZ(KDLON,KFLEV)
      REAL*8 PRMU(KDLON)
      REAL*8 PSEC(KDLON)
      REAL*8 PTAU(KDLON,2,KFLEV)
      REAL*8 PUD(KDLON,5,KFLEV+1)
C
      REAL*8 PFD(KDLON,KFLEV+1)
      REAL*8 PFU(KDLON,KFLEV+1)
C
C* LOCAL VARIABLES:
C
      INTEGER IIND(4)
C      
      REAL*8 ZCGAZ(KDLON,KFLEV) 
      REAL*8 ZDIFF(KDLON)
      REAL*8 ZDIRF(KDLON)        
      REAL*8 ZPIZAZ(KDLON,KFLEV)
      REAL*8 ZRAYL(KDLON)
      REAL*8 ZRAY1(KDLON,KFLEV+1)
      REAL*8 ZRAY2(KDLON,KFLEV+1)
      REAL*8 ZREFZ(KDLON,2,KFLEV+1)
      REAL*8 ZRJ(KDLON,6,KFLEV+1)
      REAL*8 ZRJ0(KDLON,6,KFLEV+1)
      REAL*8 ZRK(KDLON,6,KFLEV+1)
      REAL*8 ZRK0(KDLON,6,KFLEV+1)
      REAL*8 ZRMUE(KDLON,KFLEV+1)
      REAL*8 ZRMU0(KDLON,KFLEV+1)
      REAL*8 ZR(KDLON,4)
      REAL*8 ZTAUAZ(KDLON,KFLEV)
      REAL*8 ZTRA1(KDLON,KFLEV+1)
      REAL*8 ZTRA2(KDLON,KFLEV+1)
      REAL*8 ZW(KDLON,4)
C
      INTEGER jl, jk, k, jaj, ikm1, ikl
c
c Prescribed Data:
c
      REAL*8 RSUN(2)
      SAVE RSUN
      REAL*8 RRAY(2,6)
      SAVE RRAY
      DATA RSUN(1) / 0.441676 /
      DATA RSUN(2) / 0.558324 /
      DATA (RRAY(1,K),K=1,6) /
     S .428937E-01, .890743E+00,-.288555E+01,
     S .522744E+01,-.469173E+01, .161645E+01/
      DATA (RRAY(2,K),K=1,6) /
     S .697200E-02, .173297E-01,-.850903E-01,
     S .248261E+00,-.302031E+00, .129662E+00/
C     ------------------------------------------------------------------
C
C*         1.     FIRST SPECTRAL INTERVAL (0.25-0.68 MICRON)
C                 ----------------------- ------------------
C
 100  CONTINUE
C
C
C*         1.1    OPTICAL THICKNESS FOR RAYLEIGH SCATTERING
C                 -----------------------------------------
C
 110  CONTINUE
C
      DO 111 JL = 1, KDLON
      ZRAYL(JL) =  RRAY(KNU,1) + PRMU(JL) * (RRAY(KNU,2) + PRMU(JL)
     S          * (RRAY(KNU,3) + PRMU(JL) * (RRAY(KNU,4) + PRMU(JL)
     S          * (RRAY(KNU,5) + PRMU(JL) *  RRAY(KNU,6)       ))))
 111  CONTINUE
C
C
C     ------------------------------------------------------------------
C
C*         2.    CONTINUUM SCATTERING CALCULATIONS
C                ---------------------------------
C
 200  CONTINUE
C
C*         2.1   CLEAR-SKY FRACTION OF THE COLUMN
C                --------------------------------
C  
 210  CONTINUE
C
      CALL SWCLR ( KNU
     S  , PAER   , PALBP , PDSIG , ZRAYL, PSEC
     S  , ZCGAZ  , ZPIZAZ, ZRAY1 , ZRAY2, ZREFZ, ZRJ0
     S  , ZRK0   , ZRMU0 , ZTAUAZ, ZTRA1, ZTRA2)
C
C
C*         2.2   CLOUDY FRACTION OF THE COLUMN
C                -----------------------------
C
 220  CONTINUE
C
      CALL SWR ( KNU
     S  , PALBD ,PCG   ,PCLD  ,PDSIG ,POMEGA,ZRAYL
     S  , PSEC  ,PTAU
     S  , ZCGAZ ,ZPIZAZ,ZRAY1 ,ZRAY2 ,ZREFZ ,ZRJ  ,ZRK,ZRMUE
     S  , ZTAUAZ,ZTRA1 ,ZTRA2)
C
C
C     ------------------------------------------------------------------
C
C*         3.    OZONE ABSORPTION
C                ----------------
C
 300  CONTINUE
C
      IIND(1)=1
      IIND(2)=3
      IIND(3)=1
      IIND(4)=3
C      
C
C*         3.1   DOWNWARD FLUXES
C                ---------------
C
 310  CONTINUE
C
      JAJ = 2
C
      DO 311 JL = 1, KDLON
      ZW(JL,1)=0.
      ZW(JL,2)=0.
      ZW(JL,3)=0.
      ZW(JL,4)=0.
      PFD(JL,KFLEV+1)=((1.-PCLEAR(JL))*ZRJ(JL,JAJ,KFLEV+1)
     S     + PCLEAR(JL) *ZRJ0(JL,JAJ,KFLEV+1)) * RSUN(KNU)
 311  CONTINUE
      DO 314 JK = 1 , KFLEV
      IKL = KFLEV+1-JK
      DO 312 JL = 1, KDLON
      ZW(JL,1)=ZW(JL,1)+PUD(JL,1,IKL)/ZRMUE(JL,IKL)
      ZW(JL,2)=ZW(JL,2)+POZ(JL,  IKL)/ZRMUE(JL,IKL)
      ZW(JL,3)=ZW(JL,3)+PUD(JL,1,IKL)/ZRMU0(JL,IKL)
      ZW(JL,4)=ZW(JL,4)+POZ(JL,  IKL)/ZRMU0(JL,IKL)
 312  CONTINUE
C
      CALL SWTT1(KNU, 4, IIND, ZW, ZR)
C
      DO 313 JL = 1, KDLON
      ZDIFF(JL) = ZR(JL,1)*ZR(JL,2)*ZRJ(JL,JAJ,IKL)
      ZDIRF(JL) = ZR(JL,3)*ZR(JL,4)*ZRJ0(JL,JAJ,IKL)
      PFD(JL,IKL) = ((1.-PCLEAR(JL)) * ZDIFF(JL)
     S                  +PCLEAR(JL)  * ZDIRF(JL)) * RSUN(KNU)
 313  CONTINUE
 314  CONTINUE
C
C
C*         3.2   UPWARD FLUXES
C                -------------
C
 320  CONTINUE
C
      DO 325 JL = 1, KDLON
      PFU(JL,1) = ((1.-PCLEAR(JL))*ZDIFF(JL)*PALBD(JL,KNU)
     S               + PCLEAR(JL) *ZDIRF(JL)*PALBP(JL,KNU))
     S          * RSUN(KNU)
 325  CONTINUE
C
      DO 328 JK = 2 , KFLEV+1
      IKM1=JK-1
      DO 326 JL = 1, KDLON
      ZW(JL,1)=ZW(JL,1)+PUD(JL,1,IKM1)*1.66
      ZW(JL,2)=ZW(JL,2)+POZ(JL,  IKM1)*1.66
      ZW(JL,3)=ZW(JL,3)+PUD(JL,1,IKM1)*1.66
      ZW(JL,4)=ZW(JL,4)+POZ(JL,  IKM1)*1.66
 326  CONTINUE
C
      CALL SWTT1(KNU, 4, IIND, ZW, ZR)
C
      DO 327 JL = 1, KDLON
      ZDIFF(JL) = ZR(JL,1)*ZR(JL,2)*ZRK(JL,JAJ,JK)
      ZDIRF(JL) = ZR(JL,3)*ZR(JL,4)*ZRK0(JL,JAJ,JK)
      PFU(JL,JK) = ((1.-PCLEAR(JL)) * ZDIFF(JL)
     S                 +PCLEAR(JL)  * ZDIRF(JL)) * RSUN(KNU)
 327  CONTINUE
 328  CONTINUE
C
C     ------------------------------------------------------------------
C
      RETURN
      END
      SUBROUTINE SW2S ( KNU
     S  ,  PAER  ,PAKI, PALBD, PALBP, PCG   , PCLD, PCLEAR, PCLDSW
     S  ,  PDSIG ,POMEGA,POZ , PRMU , PSEC  , PTAU
     S  ,  PUD   ,PWV , PQS
     S  ,  PFDOWN,PFUP                                            )
      IMPLICIT none
      include "dimensions.h"
      include "dimphy.h"
      include "raddim.h"
      include "radepsi.h"
C
C     ------------------------------------------------------------------
C     PURPOSE.
C     --------
C
C          THIS ROUTINE COMPUTES THE SHORTWAVE RADIATION FLUXES IN THE
C     SECOND SPECTRAL INTERVAL FOLLOWING FOUQUART AND BONNEL (1980).
C
C     METHOD.
C     -------
C
C          1. COMPUTES REFLECTIVITY/TRANSMISSIVITY CORRESPONDING TO
C     CONTINUUM SCATTERING
C          2. COMPUTES REFLECTIVITY/TRANSMISSIVITY CORRESPONDING FOR
C     A GREY MOLECULAR ABSORPTION
C          3. LAPLACE TRANSFORM ON THE PREVIOUS TO GET EFFECTIVE AMOUNTS
C     OF ABSORBERS
C          4. APPLY H2O AND U.M.G. TRANSMISSION FUNCTIONS
C          5. MULTIPLY BY OZONE TRANSMISSION FUNCTION
C
C     REFERENCE.
C     ----------
C
C        SEE RADIATION'S PART OF THE ECMWF RESEARCH DEPARTMENT
C        DOCUMENTATION, AND FOUQUART AND BONNEL (1980)
C
C     AUTHOR.
C     -------
C        JEAN-JACQUES MORCRETTE  *ECMWF*
C
C     MODIFICATIONS.
C     --------------
C        ORIGINAL : 89-07-14
C        94-11-15   J.-J. MORCRETTE    DIRECT/DIFFUSE ALBEDO
C     ------------------------------------------------------------------
C* ARGUMENTS:
C
      INTEGER KNU
      REAL*8 PAER(KDLON,KFLEV,5)
      REAL*8 PAKI(KDLON,2)
      REAL*8 PALBD(KDLON,2)
      REAL*8 PALBP(KDLON,2)
      REAL*8 PCG(KDLON,2,KFLEV)
      REAL*8 PCLD(KDLON,KFLEV)
      REAL*8 PCLDSW(KDLON,KFLEV)
      REAL*8 PCLEAR(KDLON)
      REAL*8 PDSIG(KDLON,KFLEV)
      REAL*8 POMEGA(KDLON,2,KFLEV)
      REAL*8 POZ(KDLON,KFLEV)
      REAL*8 PQS(KDLON,KFLEV)
      REAL*8 PRMU(KDLON)
      REAL*8 PSEC(KDLON)
      REAL*8 PTAU(KDLON,2,KFLEV)
      REAL*8 PUD(KDLON,5,KFLEV+1)
      REAL*8 PWV(KDLON,KFLEV)
C
      REAL*8 PFDOWN(KDLON,KFLEV+1)
      REAL*8 PFUP(KDLON,KFLEV+1)
C
C* LOCAL VARIABLES:
C
      INTEGER IIND2(2), IIND3(3)
      REAL*8 ZCGAZ(KDLON,KFLEV)
      REAL*8 ZFD(KDLON,KFLEV+1)
      REAL*8 ZFU(KDLON,KFLEV+1) 
      REAL*8 ZG(KDLON)
      REAL*8 ZGG(KDLON)
      REAL*8 ZPIZAZ(KDLON,KFLEV)
      REAL*8 ZRAYL(KDLON)
      REAL*8 ZRAY1(KDLON,KFLEV+1)
      REAL*8 ZRAY2(KDLON,KFLEV+1)
      REAL*8 ZREF(KDLON)
      REAL*8 ZREFZ(KDLON,2,KFLEV+1)
      REAL*8 ZRE1(KDLON)
      REAL*8 ZRE2(KDLON)
      REAL*8 ZRJ(KDLON,6,KFLEV+1)
      REAL*8 ZRJ0(KDLON,6,KFLEV+1)
      REAL*8 ZRK(KDLON,6,KFLEV+1)
      REAL*8 ZRK0(KDLON,6,KFLEV+1)
      REAL*8 ZRL(KDLON,8)
      REAL*8 ZRMUE(KDLON,KFLEV+1)
      REAL*8 ZRMU0(KDLON,KFLEV+1)
      REAL*8 ZRMUZ(KDLON)
      REAL*8 ZRNEB(KDLON)
      REAL*8 ZRUEF(KDLON,8)
      REAL*8 ZR1(KDLON) 
      REAL*8 ZR2(KDLON,2)
      REAL*8 ZR3(KDLON,3)
      REAL*8 ZR4(KDLON)
      REAL*8 ZR21(KDLON)
      REAL*8 ZR22(KDLON)
      REAL*8 ZS(KDLON)
      REAL*8 ZTAUAZ(KDLON,KFLEV)
      REAL*8 ZTO1(KDLON)
      REAL*8 ZTR(KDLON,2,KFLEV+1)
      REAL*8 ZTRA1(KDLON,KFLEV+1)
      REAL*8 ZTRA2(KDLON,KFLEV+1)
      REAL*8 ZTR1(KDLON)
      REAL*8 ZTR2(KDLON)
      REAL*8 ZW(KDLON)   
      REAL*8 ZW1(KDLON)
      REAL*8 ZW2(KDLON,2)
      REAL*8 ZW3(KDLON,3)
      REAL*8 ZW4(KDLON)
      REAL*8 ZW5(KDLON)
C
      INTEGER jl, jk, k, jaj, ikm1, ikl, jn, jabs, jkm1
      INTEGER jref, jkl, jklp1, jajp, jkki, jkkp4, jn2j, iabs
      REAL*8 ZRMUM1, ZWH2O, ZCNEB, ZAA, ZBB, ZRKI, ZRE11
C
C* Prescribed Data:
C
      REAL*8 RSUN(2)
      SAVE RSUN
      REAL*8 RRAY(2,6)
      SAVE RRAY
      DATA RSUN(1) / 0.441676 /
      DATA RSUN(2) / 0.558324 /
      DATA (RRAY(1,K),K=1,6) /
     S .428937E-01, .890743E+00,-.288555E+01,
     S .522744E+01,-.469173E+01, .161645E+01/
      DATA (RRAY(2,K),K=1,6) /
     S .697200E-02, .173297E-01,-.850903E-01,
     S .248261E+00,-.302031E+00, .129662E+00/
C
C     ------------------------------------------------------------------
C
C*         1.     SECOND SPECTRAL INTERVAL (0.68-4.00 MICRON)
C                 -------------------------------------------
C
 100  CONTINUE
C
C
C*         1.1    OPTICAL THICKNESS FOR RAYLEIGH SCATTERING
C                 -----------------------------------------
C
 110  CONTINUE
C
      DO 111 JL = 1, KDLON
      ZRMUM1 = 1. - PRMU(JL)
      ZRAYL(JL) =  RRAY(KNU,1) + ZRMUM1   * (RRAY(KNU,2) + ZRMUM1
     S          * (RRAY(KNU,3) + ZRMUM1   * (RRAY(KNU,4) + ZRMUM1
     S          * (RRAY(KNU,5) + ZRMUM1   *  RRAY(KNU,6)     ))))
 111  CONTINUE
C
C
C     ------------------------------------------------------------------
C
C*         2.    CONTINUUM SCATTERING CALCULATIONS
C                ---------------------------------
C
 200  CONTINUE
C
C*         2.1   CLEAR-SKY FRACTION OF THE COLUMN
C                --------------------------------
C  
 210  CONTINUE
C
      CALL SWCLR ( KNU
     S  , PAER   , PALBP , PDSIG , ZRAYL, PSEC 
     S  , ZCGAZ  , ZPIZAZ, ZRAY1 , ZRAY2, ZREFZ, ZRJ0
     S  , ZRK0   , ZRMU0 , ZTAUAZ, ZTRA1, ZTRA2)
C
C
C*         2.2   CLOUDY FRACTION OF THE COLUMN
C                -----------------------------
C
 220  CONTINUE
C
      CALL SWR ( KNU
     S  , PALBD , PCG   , PCLD , PDSIG, POMEGA, ZRAYL
     S  , PSEC  , PTAU
     S  , ZCGAZ , ZPIZAZ, ZRAY1, ZRAY2, ZREFZ , ZRJ  , ZRK, ZRMUE
     S  , ZTAUAZ, ZTRA1 , ZTRA2)
C
C
C     ------------------------------------------------------------------
C
C*         3.    SCATTERING CALCULATIONS WITH GREY MOLECULAR ABSORPTION
C                ------------------------------------------------------
C
 300  CONTINUE
C
      JN = 2
C
      DO 361 JABS=1,2
C
C
C*         3.1  SURFACE CONDITIONS
C               ------------------
C
 310  CONTINUE
C
      DO 311 JL = 1, KDLON
      ZREFZ(JL,2,1) = PALBD(JL,KNU)
      ZREFZ(JL,1,1) = PALBD(JL,KNU)
 311  CONTINUE
C
C
C*         3.2  INTRODUCING CLOUD EFFECTS
C               -------------------------
C
 320  CONTINUE
C
      DO 324 JK = 2 , KFLEV+1
      JKM1 = JK - 1
      IKL=KFLEV+1-JKM1
      DO 322 JL = 1, KDLON
      ZRNEB(JL) = PCLD(JL,JKM1)
      IF (JABS.EQ.1 .AND. ZRNEB(JL).GT.2.*ZEELOG) THEN
         ZWH2O=MAX(PWV(JL,JKM1),ZEELOG)
         ZCNEB=MAX(ZEELOG,MIN(ZRNEB(JL),1.-ZEELOG))
         ZBB=PUD(JL,JABS,JKM1)*PQS(JL,JKM1)/ZWH2O
         ZAA=MAX((PUD(JL,JABS,JKM1)-ZCNEB*ZBB)/(1.-ZCNEB),ZEELOG)
      ELSE
         ZAA=PUD(JL,JABS,JKM1)
         ZBB=ZAA
      END IF
      ZRKI = PAKI(JL,JABS)
      ZS(JL) = EXP(-ZRKI * ZAA * 1.66)
      ZG(JL) = EXP(-ZRKI * ZAA / ZRMUE(JL,JK))
      ZTR1(JL) = 0.
      ZRE1(JL) = 0.
      ZTR2(JL) = 0.
      ZRE2(JL) = 0.
C
      ZW(JL)= POMEGA(JL,KNU,JKM1)
      ZTO1(JL) = PTAU(JL,KNU,JKM1) / ZW(JL)
     S               + ZTAUAZ(JL,JKM1) / ZPIZAZ(JL,JKM1)
     S               + ZBB * ZRKI
      ZR21(JL) = PTAU(JL,KNU,JKM1) + ZTAUAZ(JL,JKM1)
      ZR22(JL) = PTAU(JL,KNU,JKM1) / ZR21(JL)
      ZGG(JL) = ZR22(JL) * PCG(JL,KNU,JKM1)
     S              + (1. - ZR22(JL)) * ZCGAZ(JL,JKM1)
      ZW(JL) = ZR21(JL) / ZTO1(JL)
      ZREF(JL) = ZREFZ(JL,1,JKM1)
      ZRMUZ(JL) = ZRMUE(JL,JK)
 322  CONTINUE
C
      CALL SWDE(ZGG, ZREF, ZRMUZ, ZTO1, ZW,
     S          ZRE1, ZRE2, ZTR1, ZTR2)
C
      DO 323 JL = 1, KDLON
C
      ZREFZ(JL,2,JK) = (1.-ZRNEB(JL)) * (ZRAY1(JL,JKM1)
     S               + ZREFZ(JL,2,JKM1) * ZTRA1(JL,JKM1)
     S               * ZTRA2(JL,JKM1) ) * ZG(JL) * ZS(JL)
     S               + ZRNEB(JL) * ZRE1(JL)
C
      ZTR(JL,2,JKM1)=ZRNEB(JL)*ZTR1(JL)
     S              + (ZTRA1(JL,JKM1)) * ZG(JL) * (1.-ZRNEB(JL))
C
      ZREFZ(JL,1,JK)=(1.-ZRNEB(JL))*(ZRAY1(JL,JKM1)
     S                  +ZREFZ(JL,1,JKM1)*ZTRA1(JL,JKM1)*ZTRA2(JL,JKM1)
     S             /(1.-ZRAY2(JL,JKM1)*ZREFZ(JL,1,JKM1)))*ZG(JL)*ZS(JL)
     S             + ZRNEB(JL) * ZRE2(JL)
C
      ZTR(JL,1,JKM1)= ZRNEB(JL) * ZTR2(JL)
     S              + (ZTRA1(JL,JKM1)/(1.-ZRAY2(JL,JKM1)
     S              * ZREFZ(JL,1,JKM1)))
     S              * ZG(JL) * (1. -ZRNEB(JL))
C
 323  CONTINUE
 324  CONTINUE
C
C*         3.3  REFLECT./TRANSMISSIVITY BETWEEN SURFACE AND LEVEL
C               -------------------------------------------------
C
 330  CONTINUE
C
      DO 351 JREF=1,2
C
      JN = JN + 1
C
      DO 331 JL = 1, KDLON
      ZRJ(JL,JN,KFLEV+1) = 1.
      ZRK(JL,JN,KFLEV+1) = ZREFZ(JL,JREF,KFLEV+1)
 331  CONTINUE
C
      DO 333 JK = 1 , KFLEV
      JKL = KFLEV+1 - JK
      JKLP1 = JKL + 1
      DO 332 JL = 1, KDLON
      ZRE11 = ZRJ(JL,JN,JKLP1) * ZTR(JL,JREF,JKL)
      ZRJ(JL,JN,JKL) = ZRE11
      ZRK(JL,JN,JKL) = ZRE11 * ZREFZ(JL,JREF,JKL)
 332  CONTINUE
 333  CONTINUE
 351  CONTINUE
 361  CONTINUE
C
C
C     ------------------------------------------------------------------
C
C*         4.    INVERT GREY AND CONTINUUM FLUXES
C                --------------------------------
C
 400  CONTINUE
C
C
C*         4.1   UPWARD (ZRK) AND DOWNWARD (ZRJ) PSEUDO-FLUXES
C                ---------------------------------------------
C
 410  CONTINUE
C
      DO 414 JK = 1 , KFLEV+1
      DO 413 JAJ = 1 , 5 , 2
      JAJP = JAJ + 1
      DO 412 JL = 1, KDLON
      ZRJ(JL,JAJ,JK)=        ZRJ(JL,JAJ,JK) - ZRJ(JL,JAJP,JK)
      ZRK(JL,JAJ,JK)=        ZRK(JL,JAJ,JK) - ZRK(JL,JAJP,JK)
      ZRJ(JL,JAJ,JK)= MAX( ZRJ(JL,JAJ,JK) , ZEELOG )
      ZRK(JL,JAJ,JK)= MAX( ZRK(JL,JAJ,JK) , ZEELOG )
 412  CONTINUE
 413  CONTINUE
 414  CONTINUE
C
      DO 417 JK = 1 , KFLEV+1
      DO 416 JAJ = 2 , 6 , 2
      DO 415 JL = 1, KDLON
      ZRJ(JL,JAJ,JK)= MAX( ZRJ(JL,JAJ,JK) , ZEELOG )
      ZRK(JL,JAJ,JK)= MAX( ZRK(JL,JAJ,JK) , ZEELOG )
 415  CONTINUE
 416  CONTINUE
 417  CONTINUE
C
C*         4.2    EFFECTIVE ABSORBER AMOUNTS BY INVERSE LAPLACE
C                 ---------------------------------------------
C
 420  CONTINUE
C
      DO 437 JK = 1 , KFLEV+1
      JKKI = 1
      DO 425 JAJ = 1 , 2
      IIND2(1)=JAJ
      IIND2(2)=JAJ
      DO 424 JN = 1 , 2
      JN2J = JN + 2 * JAJ
      JKKP4 = JKKI + 4
C
C*         4.2.1  EFFECTIVE ABSORBER AMOUNTS
C                 --------------------------
C
 4210 CONTINUE
C
      DO 4211 JL = 1, KDLON
      ZW2(JL,1) = LOG( ZRJ(JL,JN,JK) / ZRJ(JL,JN2J,JK))
     S                               / PAKI(JL,JAJ)
      ZW2(JL,2) = LOG( ZRK(JL,JN,JK) / ZRK(JL,JN2J,JK))
     S                               / PAKI(JL,JAJ)
 4211 CONTINUE
C
C*         4.2.2  TRANSMISSION FUNCTION
C                 ---------------------
C
 4220 CONTINUE
C
      CALL SWTT1(KNU, 2, IIND2, ZW2, ZR2)
C
      DO 4221 JL = 1, KDLON
      ZRL(JL,JKKI) = ZR2(JL,1)
      ZRUEF(JL,JKKI) = ZW2(JL,1)
      ZRL(JL,JKKP4) = ZR2(JL,2)
      ZRUEF(JL,JKKP4) = ZW2(JL,2)
 4221 CONTINUE
C
      JKKI=JKKI+1
 424  CONTINUE
 425  CONTINUE
C
C*         4.3    UPWARD AND DOWNWARD FLUXES WITH H2O AND UMG ABSORPTION
C                 ------------------------------------------------------
C
 430  CONTINUE
C
      DO 431 JL = 1, KDLON
      PFDOWN(JL,JK) = ZRJ(JL,1,JK) * ZRL(JL,1) * ZRL(JL,3)
     S              + ZRJ(JL,2,JK) * ZRL(JL,2) * ZRL(JL,4)
      PFUP(JL,JK)   = ZRK(JL,1,JK) * ZRL(JL,5) * ZRL(JL,7)
     S              + ZRK(JL,2,JK) * ZRL(JL,6) * ZRL(JL,8)
 431  CONTINUE
 437  CONTINUE
C
C
C     ------------------------------------------------------------------
C
C*         5.    MOLECULAR ABSORPTION ON CLEAR-SKY FLUXES
C                ----------------------------------------
C
 500  CONTINUE
C
C
C*         5.1   DOWNWARD FLUXES
C                ---------------
C
 510  CONTINUE
C
      JAJ = 2
      IIND3(1)=1
      IIND3(2)=2
      IIND3(3)=3
C      
      DO 511 JL = 1, KDLON
      ZW3(JL,1)=0.
      ZW3(JL,2)=0.
      ZW3(JL,3)=0.
      ZW4(JL)  =0.
      ZW5(JL)  =0.
      ZR4(JL)  =1.
      ZFD(JL,KFLEV+1)= ZRJ0(JL,JAJ,KFLEV+1)
 511  CONTINUE
      DO 514 JK = 1 , KFLEV
      IKL = KFLEV+1-JK
      DO 512 JL = 1, KDLON
      ZW3(JL,1)=ZW3(JL,1)+PUD(JL,1,IKL)/ZRMU0(JL,IKL)
      ZW3(JL,2)=ZW3(JL,2)+PUD(JL,2,IKL)/ZRMU0(JL,IKL)
      ZW3(JL,3)=ZW3(JL,3)+POZ(JL,  IKL)/ZRMU0(JL,IKL)
      ZW4(JL)  =ZW4(JL)  +PUD(JL,4,IKL)/ZRMU0(JL,IKL)
      ZW5(JL)  =ZW5(JL)  +PUD(JL,5,IKL)/ZRMU0(JL,IKL)
 512  CONTINUE
C
      CALL SWTT1(KNU, 3, IIND3, ZW3, ZR3)
C
      DO 513 JL = 1, KDLON
C     ZR4(JL) = EXP(-RSWCE*ZW4(JL)-RSWCP*ZW5(JL))
      ZFD(JL,IKL) = ZR3(JL,1)*ZR3(JL,2)*ZR3(JL,3)*ZR4(JL)
     S            * ZRJ0(JL,JAJ,IKL)
 513  CONTINUE
 514  CONTINUE
C
C
C*         5.2   UPWARD FLUXES
C                -------------
C
 520  CONTINUE
C
      DO 525 JL = 1, KDLON
      ZFU(JL,1) = ZFD(JL,1)*PALBP(JL,KNU)
 525  CONTINUE
C
      DO 528 JK = 2 , KFLEV+1
      IKM1=JK-1
      DO 526 JL = 1, KDLON
      ZW3(JL,1)=ZW3(JL,1)+PUD(JL,1,IKM1)*1.66
      ZW3(JL,2)=ZW3(JL,2)+PUD(JL,2,IKM1)*1.66
      ZW3(JL,3)=ZW3(JL,3)+POZ(JL,  IKM1)*1.66
      ZW4(JL)  =ZW4(JL)  +PUD(JL,4,IKM1)*1.66
      ZW5(JL)  =ZW5(JL)  +PUD(JL,5,IKM1)*1.66
 526  CONTINUE
C
      CALL SWTT1(KNU, 3, IIND3, ZW3, ZR3)
C
      DO 527 JL = 1, KDLON
C     ZR4(JL) = EXP(-RSWCE*ZW4(JL)-RSWCP*ZW5(JL))
      ZFU(JL,JK) = ZR3(JL,1)*ZR3(JL,2)*ZR3(JL,3)*ZR4(JL)
     S           * ZRK0(JL,JAJ,JK)
 527  CONTINUE
 528  CONTINUE
C
C
C     ------------------------------------------------------------------
C
C*         6.     INTRODUCTION OF OZONE AND H2O CONTINUUM ABSORPTION
C                 --------------------------------------------------
C
 600  CONTINUE
      IABS=3
C
C*         6.1    DOWNWARD FLUXES
C                 ---------------
C
 610  CONTINUE
      DO 611 JL = 1, KDLON
      ZW1(JL)=0.
      ZW4(JL)=0.
      ZW5(JL)=0.
      ZR1(JL)=0.
      PFDOWN(JL,KFLEV+1) = ((1.-PCLEAR(JL))*PFDOWN(JL,KFLEV+1)
     S                   + PCLEAR(JL) * ZFD(JL,KFLEV+1)) * RSUN(KNU)
 611  CONTINUE
C
      DO 614 JK = 1 , KFLEV
      IKL=KFLEV+1-JK
      DO 612 JL = 1, KDLON
      ZW1(JL) = ZW1(JL)+POZ(JL,  IKL)/ZRMUE(JL,IKL)
      ZW4(JL) = ZW4(JL)+PUD(JL,4,IKL)/ZRMUE(JL,IKL)
      ZW5(JL) = ZW5(JL)+PUD(JL,5,IKL)/ZRMUE(JL,IKL)
C     ZR4(JL) = EXP(-RSWCE*ZW4(JL)-RSWCP*ZW5(JL))
 612  CONTINUE
C
      CALL SWTT(KNU, IABS, ZW1, ZR1)
C
      DO 613 JL = 1, KDLON
      PFDOWN(JL,IKL) = ((1.-PCLEAR(JL))*ZR1(JL)*ZR4(JL)*PFDOWN(JL,IKL)
     S                     +PCLEAR(JL)*ZFD(JL,IKL)) * RSUN(KNU)
 613  CONTINUE
 614  CONTINUE
C
C
C*         6.2    UPWARD FLUXES
C                 -------------
C
 620  CONTINUE
      DO 621 JL = 1, KDLON
      PFUP(JL,1) = ((1.-PCLEAR(JL))*ZR1(JL)*ZR4(JL) * PFUP(JL,1)
     S                 +PCLEAR(JL)*ZFU(JL,1)) * RSUN(KNU)
 621  CONTINUE
C
      DO 624 JK = 2 , KFLEV+1
      IKM1=JK-1
      DO 622 JL = 1, KDLON
      ZW1(JL) = ZW1(JL)+POZ(JL  ,IKM1)*1.66
      ZW4(JL) = ZW4(JL)+PUD(JL,4,IKM1)*1.66
      ZW5(JL) = ZW5(JL)+PUD(JL,5,IKM1)*1.66
C     ZR4(JL) = EXP(-RSWCE*ZW4(JL)-RSWCP*ZW5(JL))
 622  CONTINUE
C
      CALL SWTT(KNU, IABS, ZW1, ZR1)
C
      DO 623 JL = 1, KDLON
      PFUP(JL,JK) = ((1.-PCLEAR(JL))*ZR1(JL)*ZR4(JL) * PFUP(JL,JK)
     S                 +PCLEAR(JL)*ZFU(JL,JK)) * RSUN(KNU)
 623  CONTINUE
 624  CONTINUE
C
C     ------------------------------------------------------------------
C
      RETURN
      END
      SUBROUTINE SWCLR  ( KNU
     S  , PAER  , PALBP , PDSIG , PRAYL , PSEC
     S  , PCGAZ , PPIZAZ, PRAY1 , PRAY2 , PREFZ , PRJ  
     S  , PRK   , PRMU0 , PTAUAZ, PTRA1 , PTRA2                   )
      IMPLICIT none
      include "dimensions.h"
      include "dimphy.h"
      include "raddim.h"
      include "radepsi.h"
      include "radopt.h"
C
C     ------------------------------------------------------------------
C     PURPOSE.
C     --------
C           COMPUTES THE REFLECTIVITY AND TRANSMISSIVITY IN CASE OF
C     CLEAR-SKY COLUMN
C
C     REFERENCE.
C     ----------
C
C        SEE RADIATION'S PART OF THE ECMWF RESEARCH DEPARTMENT
C        DOCUMENTATION, AND FOUQUART AND BONNEL (1980)
C
C     AUTHOR.
C     -------
C        JEAN-JACQUES MORCRETTE  *ECMWF*
C
C     MODIFICATIONS.
C     --------------
C        ORIGINAL : 94-11-15
C     ------------------------------------------------------------------
C* ARGUMENTS:
C
      INTEGER KNU
      REAL*8 PAER(KDLON,KFLEV,5)
      REAL*8 PALBP(KDLON,2)
      REAL*8 PDSIG(KDLON,KFLEV)
      REAL*8 PRAYL(KDLON)
      REAL*8 PSEC(KDLON)
C
      REAL*8 PCGAZ(KDLON,KFLEV)     
      REAL*8 PPIZAZ(KDLON,KFLEV)
      REAL*8 PRAY1(KDLON,KFLEV+1)
      REAL*8 PRAY2(KDLON,KFLEV+1)
      REAL*8 PREFZ(KDLON,2,KFLEV+1)
      REAL*8 PRJ(KDLON,6,KFLEV+1)
      REAL*8 PRK(KDLON,6,KFLEV+1)
      REAL*8 PRMU0(KDLON,KFLEV+1)
      REAL*8 PTAUAZ(KDLON,KFLEV)
      REAL*8 PTRA1(KDLON,KFLEV+1)
      REAL*8 PTRA2(KDLON,KFLEV+1)
C
C* LOCAL VARIABLES:
C
      REAL*8 ZC0I(KDLON,KFLEV+1)       
      REAL*8 ZCLE0(KDLON,KFLEV)
      REAL*8 ZCLEAR(KDLON)
      REAL*8 ZR21(KDLON)
      REAL*8 ZR23(KDLON)
      REAL*8 ZSS0(KDLON)
      REAL*8 ZSCAT(KDLON)
      REAL*8 ZTR(KDLON,2,KFLEV+1)
C
      INTEGER jl, jk, ja, jae, jkl, jklp1, jaj, jkm1, in
      REAL*8 ZTRAY, ZGAR, ZRATIO, ZFF, ZFACOA, ZCORAE
      REAL*8 ZMUE, ZGAP, ZWW, ZTO, ZDEN, ZMU1, ZDEN1
      REAL*8 ZBMU0, ZBMU1, ZRE11
C
C* Prescribed Data for Aerosols:
C
      REAL*8 TAUA(2,5), RPIZA(2,5), RCGA(2,5)
      SAVE TAUA, RPIZA, RCGA
      DATA ((TAUA(IN,JA),JA=1,5),IN=1,2) /
     S .730719, .912819, .725059, .745405, .682188 ,
     S .730719, .912819, .725059, .745405, .682188 /
      DATA ((RPIZA(IN,JA),JA=1,5),IN=1,2) /
     S .872212, .982545, .623143, .944887, .997975 ,
     S .872212, .982545, .623143, .944887, .997975 /
      DATA ((RCGA (IN,JA),JA=1,5),IN=1,2) /
     S .647596, .739002, .580845, .662657, .624246 ,
     S .647596, .739002, .580845, .662657, .624246 /
C     ------------------------------------------------------------------
C
C*         1.    OPTICAL PARAMETERS FOR AEROSOLS AND RAYLEIGH
C                --------------------------------------------
C
 100  CONTINUE
C
      DO 103 JK = 1 , KFLEV+1
      DO 102 JA = 1 , 6
      DO 101 JL = 1, KDLON
      PRJ(JL,JA,JK) = 0.
      PRK(JL,JA,JK) = 0.
 101  CONTINUE
 102  CONTINUE
 103  CONTINUE
C
      DO 108 JK = 1 , KFLEV
      DO 104 JL = 1, KDLON
      PCGAZ(JL,JK) = 0.
      PPIZAZ(JL,JK) =  0.
      PTAUAZ(JL,JK) = 0.
 104  CONTINUE
      DO 106 JAE=1,5
      DO 105 JL = 1, KDLON
      PTAUAZ(JL,JK)=PTAUAZ(JL,JK)
     S        +PAER(JL,JK,JAE)*TAUA(KNU,JAE)
      PPIZAZ(JL,JK)=PPIZAZ(JL,JK)+PAER(JL,JK,JAE)
     S        * TAUA(KNU,JAE)*RPIZA(KNU,JAE)
      PCGAZ(JL,JK) =  PCGAZ(JL,JK) +PAER(JL,JK,JAE)
     S        * TAUA(KNU,JAE)*RPIZA(KNU,JAE)*RCGA(KNU,JAE)
 105  CONTINUE
 106  CONTINUE
C
      DO 107 JL = 1, KDLON
      IF (KAER.NE.0) THEN
         PCGAZ(JL,JK)=PCGAZ(JL,JK)/PPIZAZ(JL,JK)
         PPIZAZ(JL,JK)=PPIZAZ(JL,JK)/PTAUAZ(JL,JK)
         ZTRAY = PRAYL(JL) * PDSIG(JL,JK)
         ZRATIO = ZTRAY / (ZTRAY + PTAUAZ(JL,JK))
         ZGAR = PCGAZ(JL,JK)
         ZFF = ZGAR * ZGAR
         PTAUAZ(JL,JK)=ZTRAY+PTAUAZ(JL,JK)*(1.-PPIZAZ(JL,JK)*ZFF)
         PCGAZ(JL,JK) = ZGAR * (1. - ZRATIO) / (1. + ZGAR)
         PPIZAZ(JL,JK) =ZRATIO+(1.-ZRATIO)*PPIZAZ(JL,JK)*(1.-ZFF)
     S                       / (1. - PPIZAZ(JL,JK) * ZFF)
      ELSE
         ZTRAY = PRAYL(JL) * PDSIG(JL,JK)
         PTAUAZ(JL,JK) = ZTRAY
         PCGAZ(JL,JK) = 0.
         PPIZAZ(JL,JK) = 1.-REPSCT
      END IF
 107  CONTINUE
c      PRINT 9107,JK,((PAER(JL,JK,JAE),JAE=1,5)
c     $ ,PTAUAZ(JL,JK),PPIZAZ(JL,JK),PCGAZ(JL,JK),JL=1,KDLON)
c 9107 FORMAT(1X,'SWCLR_107',I3,8E12.5)
C
 108  CONTINUE
C
C     ------------------------------------------------------------------
C
C*         2.    TOTAL EFFECTIVE CLOUDINESS ABOVE A GIVEN LEVEL
C                ----------------------------------------------
C
 200  CONTINUE
C
      DO 201 JL = 1, KDLON
      ZR23(JL) = 0.
      ZC0I(JL,KFLEV+1) = 0.
      ZCLEAR(JL) = 1.
      ZSCAT(JL) = 0.
 201  CONTINUE
C
      JK = 1
      JKL = KFLEV+1 - JK
      JKLP1 = JKL + 1
      DO 202 JL = 1, KDLON
      ZFACOA = 1. - PPIZAZ(JL,JKL)*PCGAZ(JL,JKL)*PCGAZ(JL,JKL)
      ZCORAE = ZFACOA * PTAUAZ(JL,JKL) * PSEC(JL)
      ZR21(JL) = EXP(-ZCORAE   )
      ZSS0(JL) = 1.-ZR21(JL)
      ZCLE0(JL,JKL) = ZSS0(JL)
C
      IF (NOVLP.EQ.1) THEN
c* maximum-random
         ZCLEAR(JL) = ZCLEAR(JL)
     S                  *(1.0-MAX(ZSS0(JL),ZSCAT(JL)))
     S                  /(1.0-MIN(ZSCAT(JL),1.-ZEPSEC))
         ZC0I(JL,JKL) = 1.0 - ZCLEAR(JL)
         ZSCAT(JL) = ZSS0(JL)
      ELSE IF (NOVLP.EQ.2) THEN
C* maximum
         ZSCAT(JL) = MAX( ZSS0(JL) , ZSCAT(JL) )
         ZC0I(JL,JKL) = ZSCAT(JL)
      ELSE IF (NOVLP.EQ.3) THEN
c* random
         ZCLEAR(JL)=ZCLEAR(JL)*(1.0-ZSS0(JL))
         ZSCAT(JL) = 1.0 - ZCLEAR(JL)
         ZC0I(JL,JKL) = ZSCAT(JL)
      END IF
 202  CONTINUE
C
      DO 205 JK = 2 , KFLEV
      JKL = KFLEV+1 - JK
      JKLP1 = JKL + 1
      DO 204 JL = 1, KDLON
      ZFACOA = 1. - PPIZAZ(JL,JKL)*PCGAZ(JL,JKL)*PCGAZ(JL,JKL)
      ZCORAE = ZFACOA * PTAUAZ(JL,JKL) * PSEC(JL)
      ZR21(JL) = EXP(-ZCORAE   )
      ZSS0(JL) = 1.-ZR21(JL)
      ZCLE0(JL,JKL) = ZSS0(JL)
c     
      IF (NOVLP.EQ.1) THEN
c* maximum-random
         ZCLEAR(JL) = ZCLEAR(JL)
     S                  *(1.0-MAX(ZSS0(JL),ZSCAT(JL)))
     S                  /(1.0-MIN(ZSCAT(JL),1.-ZEPSEC))
         ZC0I(JL,JKL) = 1.0 - ZCLEAR(JL)
         ZSCAT(JL) = ZSS0(JL)
      ELSE IF (NOVLP.EQ.2) THEN
C* maximum
         ZSCAT(JL) = MAX( ZSS0(JL) , ZSCAT(JL) )
         ZC0I(JL,JKL) = ZSCAT(JL)
      ELSE IF (NOVLP.EQ.3) THEN
c* random
         ZCLEAR(JL)=ZCLEAR(JL)*(1.0-ZSS0(JL))
         ZSCAT(JL) = 1.0 - ZCLEAR(JL)
         ZC0I(JL,JKL) = ZSCAT(JL)
      END IF                  
 204  CONTINUE
 205  CONTINUE
C
C     ------------------------------------------------------------------
C
C*         3.    REFLECTIVITY/TRANSMISSIVITY FOR PURE SCATTERING
C                -----------------------------------------------
C
 300  CONTINUE
C
      DO 301 JL = 1, KDLON
      PRAY1(JL,KFLEV+1) = 0.
      PRAY2(JL,KFLEV+1) = 0.
      PREFZ(JL,2,1) = PALBP(JL,KNU)
      PREFZ(JL,1,1) = PALBP(JL,KNU)
      PTRA1(JL,KFLEV+1) = 1.
      PTRA2(JL,KFLEV+1) = 1.
 301  CONTINUE
C
      DO 346 JK = 2 , KFLEV+1
      JKM1 = JK-1
      DO 342 JL = 1, KDLON
C
C
C     ------------------------------------------------------------------
C
C*         3.1  EQUIVALENT ZENITH ANGLE
C               -----------------------
C
 310  CONTINUE
C
      ZMUE = (1.-ZC0I(JL,JK)) * PSEC(JL)
     S            + ZC0I(JL,JK) * 1.66
      PRMU0(JL,JK) = 1./ZMUE
C
C
C     ------------------------------------------------------------------
C
C*         3.2  REFLECT./TRANSMISSIVITY DUE TO RAYLEIGH AND AEROSOLS
C               ----------------------------------------------------
C
 320  CONTINUE
C
      ZGAP = PCGAZ(JL,JKM1)
      ZBMU0 = 0.5 - 0.75 * ZGAP / ZMUE
      ZWW = PPIZAZ(JL,JKM1)
      ZTO = PTAUAZ(JL,JKM1)
      ZDEN = 1. + (1. - ZWW + ZBMU0 * ZWW) * ZTO * ZMUE
     S       + (1-ZWW) * (1. - ZWW +2.*ZBMU0*ZWW)*ZTO*ZTO*ZMUE*ZMUE
      PRAY1(JL,JKM1) = ZBMU0 * ZWW * ZTO * ZMUE / ZDEN
      PTRA1(JL,JKM1) = 1. / ZDEN
C
      ZMU1 = 0.5
      ZBMU1 = 0.5 - 0.75 * ZGAP * ZMU1
      ZDEN1= 1. + (1. - ZWW + ZBMU1 * ZWW) * ZTO / ZMU1
     S       + (1-ZWW) * (1. - ZWW +2.*ZBMU1*ZWW)*ZTO*ZTO/ZMU1/ZMU1
      PRAY2(JL,JKM1) = ZBMU1 * ZWW * ZTO / ZMU1 / ZDEN1
      PTRA2(JL,JKM1) = 1. / ZDEN1
C
C
C
      PREFZ(JL,1,JK) = (PRAY1(JL,JKM1)
     S               + PREFZ(JL,1,JKM1) * PTRA1(JL,JKM1)
     S               * PTRA2(JL,JKM1)
     S               / (1.-PRAY2(JL,JKM1)*PREFZ(JL,1,JKM1)))
C
      ZTR(JL,1,JKM1) = (PTRA1(JL,JKM1)
     S               / (1.-PRAY2(JL,JKM1)*PREFZ(JL,1,JKM1)))
C
      PREFZ(JL,2,JK) = (PRAY1(JL,JKM1)
     S               + PREFZ(JL,2,JKM1) * PTRA1(JL,JKM1)
     S               * PTRA2(JL,JKM1) )
C
      ZTR(JL,2,JKM1) = PTRA1(JL,JKM1) 
C
 342  CONTINUE
 346  CONTINUE
      DO 347 JL = 1, KDLON
      ZMUE = (1.-ZC0I(JL,1))*PSEC(JL)+ZC0I(JL,1)*1.66
      PRMU0(JL,1)=1./ZMUE
 347  CONTINUE
C
C
C     ------------------------------------------------------------------
C
C*         3.5    REFLECT./TRANSMISSIVITY BETWEEN SURFACE AND LEVEL
C                 -------------------------------------------------
C
 350  CONTINUE
C
      IF (KNU.EQ.1) THEN
      JAJ = 2
      DO 351 JL = 1, KDLON
      PRJ(JL,JAJ,KFLEV+1) = 1.
      PRK(JL,JAJ,KFLEV+1) = PREFZ(JL, 1,KFLEV+1)
 351  CONTINUE
C
      DO 353 JK = 1 , KFLEV
      JKL = KFLEV+1 - JK
      JKLP1 = JKL + 1
      DO 352 JL = 1, KDLON
      ZRE11= PRJ(JL,JAJ,JKLP1) * ZTR(JL,  1,JKL)
      PRJ(JL,JAJ,JKL) = ZRE11
      PRK(JL,JAJ,JKL) = ZRE11 * PREFZ(JL,  1,JKL)
 352  CONTINUE
 353  CONTINUE
 354  CONTINUE
C
      ELSE
C
      DO 358 JAJ = 1 , 2
      DO 355 JL = 1, KDLON
      PRJ(JL,JAJ,KFLEV+1) = 1.
      PRK(JL,JAJ,KFLEV+1) = PREFZ(JL,JAJ,KFLEV+1)
 355  CONTINUE
C
      DO 357 JK = 1 , KFLEV
      JKL = KFLEV+1 - JK
      JKLP1 = JKL + 1
      DO 356 JL = 1, KDLON
      ZRE11= PRJ(JL,JAJ,JKLP1) * ZTR(JL,JAJ,JKL)
      PRJ(JL,JAJ,JKL) = ZRE11
      PRK(JL,JAJ,JKL) = ZRE11 * PREFZ(JL,JAJ,JKL)
 356  CONTINUE
 357  CONTINUE
 358  CONTINUE
C
      END IF
C
C     ------------------------------------------------------------------
C
      RETURN
      END
      SUBROUTINE SWR ( KNU
     S  , PALBD , PCG   , PCLD , PDSIG, POMEGA, PRAYL
     S  , PSEC  , PTAU
     S  , PCGAZ , PPIZAZ, PRAY1, PRAY2, PREFZ , PRJ  , PRK , PRMUE
     S  , PTAUAZ, PTRA1 , PTRA2 )
      IMPLICIT none
      include "dimensions.h"
      include "dimphy.h"
      include "raddim.h"
      include "radepsi.h"
      include "radopt.h"
C
C     ------------------------------------------------------------------
C     PURPOSE.
C     --------
C           COMPUTES THE REFLECTIVITY AND TRANSMISSIVITY IN CASE OF
C     CONTINUUM SCATTERING
C
C     METHOD.
C     -------
C
C          1. COMPUTES CONTINUUM FLUXES CORRESPONDING TO AEROSOL
C     OR/AND RAYLEIGH SCATTERING (NO MOLECULAR GAS ABSORPTION)
C
C     REFERENCE.
C     ----------
C
C        SEE RADIATION'S PART OF THE ECMWF RESEARCH DEPARTMENT
C        DOCUMENTATION, AND FOUQUART AND BONNEL (1980)
C
C     AUTHOR.
C     -------
C        JEAN-JACQUES MORCRETTE  *ECMWF*
C
C     MODIFICATIONS.
C     --------------
C        ORIGINAL : 89-07-14
C     ------------------------------------------------------------------
C* ARGUMENTS:
C
      INTEGER KNU
      REAL*8 PALBD(KDLON,2)
      REAL*8 PCG(KDLON,2,KFLEV)
      REAL*8 PCLD(KDLON,KFLEV)
      REAL*8 PDSIG(KDLON,KFLEV)
      REAL*8 POMEGA(KDLON,2,KFLEV)
      REAL*8 PRAYL(KDLON)
      REAL*8 PSEC(KDLON)
      REAL*8 PTAU(KDLON,2,KFLEV)
C
      REAL*8 PRAY1(KDLON,KFLEV+1)
      REAL*8 PRAY2(KDLON,KFLEV+1)
      REAL*8 PREFZ(KDLON,2,KFLEV+1)
      REAL*8 PRJ(KDLON,6,KFLEV+1)
      REAL*8 PRK(KDLON,6,KFLEV+1)
      REAL*8 PRMUE(KDLON,KFLEV+1)
      REAL*8 PCGAZ(KDLON,KFLEV)
      REAL*8 PPIZAZ(KDLON,KFLEV)
      REAL*8 PTAUAZ(KDLON,KFLEV)
      REAL*8 PTRA1(KDLON,KFLEV+1)
      REAL*8 PTRA2(KDLON,KFLEV+1)
C
C* LOCAL VARIABLES:
C
      REAL*8 ZC1I(KDLON,KFLEV+1)
      REAL*8 ZCLEQ(KDLON,KFLEV)
      REAL*8 ZCLEAR(KDLON)
      REAL*8 ZCLOUD(KDLON)
      REAL*8 ZGG(KDLON)
      REAL*8 ZREF(KDLON)
      REAL*8 ZRE1(KDLON)
      REAL*8 ZRE2(KDLON)
      REAL*8 ZRMUZ(KDLON)
      REAL*8 ZRNEB(KDLON)
      REAL*8 ZR21(KDLON)
      REAL*8 ZR22(KDLON)
      REAL*8 ZR23(KDLON)
      REAL*8 ZSS1(KDLON)
      REAL*8 ZTO1(KDLON)
      REAL*8 ZTR(KDLON,2,KFLEV+1)
      REAL*8 ZTR1(KDLON)
      REAL*8 ZTR2(KDLON)
      REAL*8 ZW(KDLON)
C
      INTEGER jk, jl, ja, jkl, jklp1, jkm1, jaj
      REAL*8 ZFACOA, ZFACOC, ZCORAE, ZCORCD
      REAL*8 ZMUE, ZGAP, ZWW, ZTO, ZDEN, ZDEN1
      REAL*8 ZMU1, ZRE11, ZBMU0, ZBMU1
C
C     ------------------------------------------------------------------
C
C*         1.    INITIALIZATION
C                --------------
C
 100  CONTINUE
C
      DO 103 JK = 1 , KFLEV+1
      DO 102 JA = 1 , 6
      DO 101 JL = 1, KDLON
      PRJ(JL,JA,JK) = 0.
      PRK(JL,JA,JK) = 0.
 101  CONTINUE
 102  CONTINUE
 103  CONTINUE
C
C
C     ------------------------------------------------------------------
C
C*         2.    TOTAL EFFECTIVE CLOUDINESS ABOVE A GIVEN LEVEL
C                ----------------------------------------------
C
 200  CONTINUE
C
      DO 201 JL = 1, KDLON
      ZR23(JL) = 0.
      ZC1I(JL,KFLEV+1) = 0.
      ZCLEAR(JL) = 1.
      ZCLOUD(JL) = 0.
 201  CONTINUE
C
      JK = 1
      JKL = KFLEV+1 - JK
      JKLP1 = JKL + 1
      DO 202 JL = 1, KDLON
      ZFACOA = 1. - PPIZAZ(JL,JKL)*PCGAZ(JL,JKL)*PCGAZ(JL,JKL)
      ZFACOC = 1. - POMEGA(JL,KNU,JKL) * PCG(JL,KNU,JKL)
     S                                 * PCG(JL,KNU,JKL)
      ZCORAE = ZFACOA * PTAUAZ(JL,JKL) * PSEC(JL)
      ZCORCD = ZFACOC * PTAU(JL,KNU,JKL) * PSEC(JL)
      ZR21(JL) = EXP(-ZCORAE   )
      ZR22(JL) = EXP(-ZCORCD   )
      ZSS1(JL) = PCLD(JL,JKL)*(1.0-ZR21(JL)*ZR22(JL))
     S               + (1.0-PCLD(JL,JKL))*(1.0-ZR21(JL))
      ZCLEQ(JL,JKL) = ZSS1(JL)
C
      IF (NOVLP.EQ.1) THEN
c* maximum-random
         ZCLEAR(JL) = ZCLEAR(JL)
     S                  *(1.0-MAX(ZSS1(JL),ZCLOUD(JL)))
     S                  /(1.0-MIN(ZCLOUD(JL),1.-ZEPSEC))
         ZC1I(JL,JKL) = 1.0 - ZCLEAR(JL)
         ZCLOUD(JL) = ZSS1(JL)
      ELSE IF (NOVLP.EQ.2) THEN
C* maximum
         ZCLOUD(JL) = MAX( ZSS1(JL) , ZCLOUD(JL) )
         ZC1I(JL,JKL) = ZCLOUD(JL)
      ELSE IF (NOVLP.EQ.3) THEN
c* random
         ZCLEAR(JL) = ZCLEAR(JL)*(1.0 - ZSS1(JL))
         ZCLOUD(JL) = 1.0 - ZCLEAR(JL)
         ZC1I(JL,JKL) = ZCLOUD(JL)
      END IF
 202  CONTINUE
C
      DO 205 JK = 2 , KFLEV
      JKL = KFLEV+1 - JK
      JKLP1 = JKL + 1
      DO 204 JL = 1, KDLON
      ZFACOA = 1. - PPIZAZ(JL,JKL)*PCGAZ(JL,JKL)*PCGAZ(JL,JKL)
      ZFACOC = 1. - POMEGA(JL,KNU,JKL) * PCG(JL,KNU,JKL)
     S                                 * PCG(JL,KNU,JKL)
      ZCORAE = ZFACOA * PTAUAZ(JL,JKL) * PSEC(JL)
      ZCORCD = ZFACOC * PTAU(JL,KNU,JKL) * PSEC(JL)
      ZR21(JL) = EXP(-ZCORAE   )
      ZR22(JL) = EXP(-ZCORCD   )
      ZSS1(JL) = PCLD(JL,JKL)*(1.0-ZR21(JL)*ZR22(JL))
     S               + (1.0-PCLD(JL,JKL))*(1.0-ZR21(JL))
      ZCLEQ(JL,JKL) = ZSS1(JL)
c     
      IF (NOVLP.EQ.1) THEN
c* maximum-random
         ZCLEAR(JL) = ZCLEAR(JL)
     S                  *(1.0-MAX(ZSS1(JL),ZCLOUD(JL)))
     S                  /(1.0-MIN(ZCLOUD(JL),1.-ZEPSEC))
         ZC1I(JL,JKL) = 1.0 - ZCLEAR(JL)
         ZCLOUD(JL) = ZSS1(JL)
      ELSE IF (NOVLP.EQ.2) THEN
C* maximum
         ZCLOUD(JL) = MAX( ZSS1(JL) , ZCLOUD(JL) )
         ZC1I(JL,JKL) = ZCLOUD(JL)
      ELSE IF (NOVLP.EQ.3) THEN
c* random
         ZCLEAR(JL) = ZCLEAR(JL)*(1.0 - ZSS1(JL))
         ZCLOUD(JL) = 1.0 - ZCLEAR(JL)
         ZC1I(JL,JKL) = ZCLOUD(JL)
      END IF
 204  CONTINUE
 205  CONTINUE
C
C     ------------------------------------------------------------------
C
C*         3.    REFLECTIVITY/TRANSMISSIVITY FOR PURE SCATTERING
C                -----------------------------------------------
C
 300  CONTINUE
C
      DO 301 JL = 1, KDLON
      PRAY1(JL,KFLEV+1) = 0.
      PRAY2(JL,KFLEV+1) = 0.
      PREFZ(JL,2,1) = PALBD(JL,KNU)
      PREFZ(JL,1,1) = PALBD(JL,KNU)
      PTRA1(JL,KFLEV+1) = 1.
      PTRA2(JL,KFLEV+1) = 1.
 301  CONTINUE
C
      DO 346 JK = 2 , KFLEV+1
      JKM1 = JK-1
      DO 342 JL = 1, KDLON
      ZRNEB(JL)= PCLD(JL,JKM1)
      ZRE1(JL)=0.
      ZTR1(JL)=0.
      ZRE2(JL)=0.
      ZTR2(JL)=0.
C
C
C     ------------------------------------------------------------------
C
C*         3.1  EQUIVALENT ZENITH ANGLE
C               -----------------------
C
 310  CONTINUE
C
      ZMUE = (1.-ZC1I(JL,JK)) * PSEC(JL)
     S            + ZC1I(JL,JK) * 1.66
      PRMUE(JL,JK) = 1./ZMUE
C
C
C     ------------------------------------------------------------------
C
C*         3.2  REFLECT./TRANSMISSIVITY DUE TO RAYLEIGH AND AEROSOLS
C               ----------------------------------------------------
C
 320  CONTINUE
C
      ZGAP = PCGAZ(JL,JKM1)
      ZBMU0 = 0.5 - 0.75 * ZGAP / ZMUE
      ZWW = PPIZAZ(JL,JKM1)
      ZTO = PTAUAZ(JL,JKM1)
      ZDEN = 1. + (1. - ZWW + ZBMU0 * ZWW) * ZTO * ZMUE
     S       + (1-ZWW) * (1. - ZWW +2.*ZBMU0*ZWW)*ZTO*ZTO*ZMUE*ZMUE
      PRAY1(JL,JKM1) = ZBMU0 * ZWW * ZTO * ZMUE / ZDEN
      PTRA1(JL,JKM1) = 1. / ZDEN
c      PRINT *,' LOOP 342 ** 3 ** JL=',JL,PRAY1(JL,JKM1),PTRA1(JL,JKM1)
C
      ZMU1 = 0.5
      ZBMU1 = 0.5 - 0.75 * ZGAP * ZMU1
      ZDEN1= 1. + (1. - ZWW + ZBMU1 * ZWW) * ZTO / ZMU1
     S       + (1-ZWW) * (1. - ZWW +2.*ZBMU1*ZWW)*ZTO*ZTO/ZMU1/ZMU1
      PRAY2(JL,JKM1) = ZBMU1 * ZWW * ZTO / ZMU1 / ZDEN1
      PTRA2(JL,JKM1) = 1. / ZDEN1
C
C
C     ------------------------------------------------------------------
C
C*         3.3  EFFECT OF CLOUD LAYER
C               ---------------------
C
 330  CONTINUE
C
      ZW(JL) = POMEGA(JL,KNU,JKM1)
      ZTO1(JL) = PTAU(JL,KNU,JKM1)/ZW(JL)
     S         + PTAUAZ(JL,JKM1)/PPIZAZ(JL,JKM1)
      ZR21(JL) = PTAU(JL,KNU,JKM1) + PTAUAZ(JL,JKM1)
      ZR22(JL) = PTAU(JL,KNU,JKM1) / ZR21(JL)
      ZGG(JL) = ZR22(JL) * PCG(JL,KNU,JKM1)
     S              + (1. - ZR22(JL)) * PCGAZ(JL,JKM1)
C Modif PhD - JJM 19/03/96 pour erreurs arrondis
C machine
C PHD PROTECTION ZW(JL) = ZR21(JL) / ZTO1(JL)
      IF (ZW(JL).EQ.1. .AND. PPIZAZ(JL,JKM1).EQ.1.) THEN
         ZW(JL)=1.
      ELSE
         ZW(JL) = ZR21(JL) / ZTO1(JL)
      END IF
      ZREF(JL) = PREFZ(JL,1,JKM1)
      ZRMUZ(JL) = PRMUE(JL,JK)
 342  CONTINUE
C
      CALL SWDE(ZGG  , ZREF  , ZRMUZ , ZTO1 , ZW,
     S          ZRE1 , ZRE2  , ZTR1  , ZTR2)
C
      DO 345 JL = 1, KDLON
C
      PREFZ(JL,1,JK) = (1.-ZRNEB(JL)) * (PRAY1(JL,JKM1)
     S               + PREFZ(JL,1,JKM1) * PTRA1(JL,JKM1)
     S               * PTRA2(JL,JKM1)
     S               / (1.-PRAY2(JL,JKM1)*PREFZ(JL,1,JKM1)))
     S               + ZRNEB(JL) * ZRE2(JL)
C
      ZTR(JL,1,JKM1) = ZRNEB(JL) * ZTR2(JL) + (PTRA1(JL,JKM1)
     S               / (1.-PRAY2(JL,JKM1)*PREFZ(JL,1,JKM1)))
     S               * (1.-ZRNEB(JL))
C
      PREFZ(JL,2,JK) = (1.-ZRNEB(JL)) * (PRAY1(JL,JKM1)
     S               + PREFZ(JL,2,JKM1) * PTRA1(JL,JKM1)
     S               * PTRA2(JL,JKM1) )
     S               + ZRNEB(JL) * ZRE1(JL)
C
      ZTR(JL,2,JKM1) = ZRNEB(JL) * ZTR1(JL)
     S               + PTRA1(JL,JKM1) * (1.-ZRNEB(JL))
C
 345  CONTINUE
 346  CONTINUE
      DO 347 JL = 1, KDLON
      ZMUE = (1.-ZC1I(JL,1))*PSEC(JL)+ZC1I(JL,1)*1.66
      PRMUE(JL,1)=1./ZMUE
 347  CONTINUE
C
C
C     ------------------------------------------------------------------
C
C*         3.5    REFLECT./TRANSMISSIVITY BETWEEN SURFACE AND LEVEL
C                 -------------------------------------------------
C
 350  CONTINUE
C
      IF (KNU.EQ.1) THEN
      JAJ = 2
      DO 351 JL = 1, KDLON
      PRJ(JL,JAJ,KFLEV+1) = 1.
      PRK(JL,JAJ,KFLEV+1) = PREFZ(JL, 1,KFLEV+1)
 351  CONTINUE
C
      DO 353 JK = 1 , KFLEV
      JKL = KFLEV+1 - JK
      JKLP1 = JKL + 1
      DO 352 JL = 1, KDLON
      ZRE11= PRJ(JL,JAJ,JKLP1) * ZTR(JL,  1,JKL)
      PRJ(JL,JAJ,JKL) = ZRE11
      PRK(JL,JAJ,JKL) = ZRE11 * PREFZ(JL,  1,JKL)
 352  CONTINUE
 353  CONTINUE
 354  CONTINUE
C
      ELSE
C
      DO 358 JAJ = 1 , 2
      DO 355 JL = 1, KDLON
      PRJ(JL,JAJ,KFLEV+1) = 1.
      PRK(JL,JAJ,KFLEV+1) = PREFZ(JL,JAJ,KFLEV+1)
 355  CONTINUE
C
      DO 357 JK = 1 , KFLEV
      JKL = KFLEV+1 - JK
      JKLP1 = JKL + 1
      DO 356 JL = 1, KDLON
      ZRE11= PRJ(JL,JAJ,JKLP1) * ZTR(JL,JAJ,JKL)
      PRJ(JL,JAJ,JKL) = ZRE11
      PRK(JL,JAJ,JKL) = ZRE11 * PREFZ(JL,JAJ,JKL)
 356  CONTINUE
 357  CONTINUE
 358  CONTINUE
C
      END IF
C
C     ------------------------------------------------------------------
C
      RETURN
      END
      SUBROUTINE SWDE (PGG,PREF,PRMUZ,PTO1,PW,
     S                 PRE1,PRE2,PTR1,PTR2)
      IMPLICIT none
      include "dimensions.h"
      include "dimphy.h"
      include "raddim.h"
C
C     ------------------------------------------------------------------
C     PURPOSE.
C     --------
C           COMPUTES THE REFLECTIVITY AND TRANSMISSIVITY OF A CLOUDY
C     LAYER USING THE DELTA-EDDINGTON'S APPROXIMATION.
C
C     METHOD.
C     -------
C
C          STANDARD DELTA-EDDINGTON LAYER CALCULATIONS.
C
C     REFERENCE.
C     ----------
C
C        SEE RADIATION'S PART OF THE MODEL'S DOCUMENTATION AND
C        ECMWF RESEARCH DEPARTMENT DOCUMENTATION OF THE IFS
C
C     AUTHOR.
C     -------
C        JEAN-JACQUES MORCRETTE  *ECMWF*
C
C     MODIFICATIONS.
C     --------------
C        ORIGINAL : 88-12-15
C     ------------------------------------------------------------------
C* ARGUMENTS:
C
      REAL*8 PGG(KDLON)   ! ASSYMETRY FACTOR
      REAL*8 PREF(KDLON)  ! REFLECTIVITY OF THE UNDERLYING LAYER
      REAL*8 PRMUZ(KDLON) ! COSINE OF SOLAR ZENITH ANGLE
      REAL*8 PTO1(KDLON)  ! OPTICAL THICKNESS
      REAL*8 PW(KDLON)    ! SINGLE SCATTERING ALBEDO
      REAL*8 PRE1(KDLON)  ! LAYER REFLECTIVITY (NO UNDERLYING-LAYER REFLECTION)
      REAL*8 PRE2(KDLON)  ! LAYER REFLECTIVITY
      REAL*8 PTR1(KDLON)  ! LAYER TRANSMISSIVITY (NO UNDERLYING-LAYER REFLECTION)
      REAL*8 PTR2(KDLON)  ! LAYER TRANSMISSIVITY
C
C* LOCAL VARIABLES:
C
      INTEGER jl
      REAL*8 ZFF, ZGP, ZTOP, ZWCP, ZDT, ZX1, ZWM
      REAL*8 ZRM2, ZRK, ZX2, ZRP, ZALPHA, ZBETA, ZARG
      REAL*8 ZEXMU0, ZARG2, ZEXKP, ZEXKM, ZXP2P, ZXM2P, ZAP2B, ZAM2B
      REAL*8 ZA11, ZA12, ZA13, ZA21, ZA22, ZA23
      REAL*8 ZDENA, ZC1A, ZC2A, ZRI0A, ZRI1A
      REAL*8 ZRI0B, ZRI1B
      REAL*8 ZB21, ZB22, ZB23, ZDENB, ZC1B, ZC2B
      REAL*8 ZRI0C, ZRI1C, ZRI0D, ZRI1D
C     ------------------------------------------------------------------
C
C*         1.      DELTA-EDDINGTON CALCULATIONS
C
 100  CONTINUE
C
      DO 131 JL   =   1, KDLON
C
C*         1.1     SET UP THE DELTA-MODIFIED PARAMETERS
C
 110  CONTINUE
C
      ZFF = PGG(JL)*PGG(JL)
      ZGP = PGG(JL)/(1.+PGG(JL))
      ZTOP = (1.- PW(JL) * ZFF) * PTO1(JL)
      ZWCP = (1-ZFF)* PW(JL) /(1.- PW(JL) * ZFF)
      ZDT = 2./3.
      ZX1 = 1.-ZWCP*ZGP
      ZWM = 1.-ZWCP
      ZRM2 =  PRMUZ(JL) * PRMUZ(JL)
      ZRK = SQRT(3.*ZWM*ZX1)
      ZX2 = 4.*(1.-ZRK*ZRK*ZRM2)
      ZRP=ZRK/ZX1
      ZALPHA = 3.*ZWCP*ZRM2*(1.+ZGP*ZWM)/ZX2
      ZBETA = 3.*ZWCP* PRMUZ(JL) *(1.+3.*ZGP*ZRM2*ZWM)/ZX2
      ZARG=MIN(ZTOP/PRMUZ(JL),200.)
      ZEXMU0=EXP(-ZARG)
      ZARG2=MIN(ZRK*ZTOP,200.)
      ZEXKP=EXP(ZARG2)
      ZEXKM = 1./ZEXKP
      ZXP2P = 1.+ZDT*ZRP
      ZXM2P = 1.-ZDT*ZRP
      ZAP2B = ZALPHA+ZDT*ZBETA
      ZAM2B = ZALPHA-ZDT*ZBETA
C
C*         1.2     WITHOUT REFLECTION FROM THE UNDERLYING LAYER
C
 120  CONTINUE
C
      ZA11 = ZXP2P
      ZA12 = ZXM2P
      ZA13 = ZAP2B
      ZA22 = ZXP2P*ZEXKP
      ZA21 = ZXM2P*ZEXKM
      ZA23 = ZAM2B*ZEXMU0
      ZDENA = ZA11 * ZA22 - ZA21 * ZA12
      ZC1A = (ZA22*ZA13-ZA12*ZA23)/ZDENA
      ZC2A = (ZA11*ZA23-ZA21*ZA13)/ZDENA
      ZRI0A = ZC1A+ZC2A-ZALPHA
      ZRI1A = ZRP*(ZC1A-ZC2A)-ZBETA
      PRE1(JL) = (ZRI0A-ZDT*ZRI1A)/ PRMUZ(JL)
      ZRI0B = ZC1A*ZEXKM+ZC2A*ZEXKP-ZALPHA*ZEXMU0
      ZRI1B = ZRP*(ZC1A*ZEXKM-ZC2A*ZEXKP)-ZBETA*ZEXMU0
      PTR1(JL) = ZEXMU0+(ZRI0B+ZDT*ZRI1B)/ PRMUZ(JL)
C
C*         1.3     WITH REFLECTION FROM THE UNDERLYING LAYER
C
 130  CONTINUE
C
      ZB21 = ZA21- PREF(JL) *ZXP2P*ZEXKM
      ZB22 = ZA22- PREF(JL) *ZXM2P*ZEXKP
      ZB23 = ZA23- PREF(JL) *ZEXMU0*(ZAP2B - PRMUZ(JL) )
      ZDENB = ZA11 * ZB22 - ZB21 * ZA12
      ZC1B = (ZB22*ZA13-ZA12*ZB23)/ZDENB
      ZC2B = (ZA11*ZB23-ZB21*ZA13)/ZDENB
      ZRI0C = ZC1B+ZC2B-ZALPHA
      ZRI1C = ZRP*(ZC1B-ZC2B)-ZBETA
      PRE2(JL) = (ZRI0C-ZDT*ZRI1C) / PRMUZ(JL)
      ZRI0D = ZC1B*ZEXKM + ZC2B*ZEXKP - ZALPHA*ZEXMU0
      ZRI1D = ZRP * (ZC1B*ZEXKM - ZC2B*ZEXKP) - ZBETA*ZEXMU0
      PTR2(JL) = ZEXMU0 + (ZRI0D + ZDT*ZRI1D) / PRMUZ(JL)
C
 131  CONTINUE
      RETURN
      END
      SUBROUTINE SWTT (KNU,KA,PU,PTR)
      IMPLICIT none
      include "dimensions.h"
      include "dimphy.h"
      include "raddim.h"
C
C-----------------------------------------------------------------------
C     PURPOSE.
C     --------
C           THIS ROUTINE COMPUTES THE TRANSMISSION FUNCTIONS FOR ALL THE
C     ABSORBERS (H2O, UNIFORMLY MIXED GASES, AND O3) IN THE TWO SPECTRAL
C     INTERVALS.
C
C     METHOD.
C     -------
C
C          TRANSMISSION FUNCTION ARE COMPUTED USING PADE APPROXIMANTS
C     AND HORNER'S ALGORITHM.
C
C     REFERENCE.
C     ----------
C
C        SEE RADIATION'S PART OF THE MODEL'S DOCUMENTATION AND
C        ECMWF RESEARCH DEPARTMENT DOCUMENTATION OF THE IFS
C
C     AUTHOR.
C     -------
C        JEAN-JACQUES MORCRETTE  *ECMWF*
C
C     MODIFICATIONS.
C     --------------
C        ORIGINAL : 88-12-15
C-----------------------------------------------------------------------
C
C* ARGUMENTS
C
      INTEGER KNU     ! INDEX OF THE SPECTRAL INTERVAL
      INTEGER KA      ! INDEX OF THE ABSORBER
      REAL*8 PU(KDLON)  ! ABSORBER AMOUNT
C
      REAL*8 PTR(KDLON) ! TRANSMISSION FUNCTION
C
C* LOCAL VARIABLES:
C
      REAL*8 ZR1(KDLON), ZR2(KDLON)
      INTEGER jl, i,j
C
C* Prescribed Data:
C
      REAL*8 APAD(2,3,7), BPAD(2,3,7), D(2,3)
      SAVE APAD, BPAD, D
      DATA ((APAD(1,I,J),I=1,3),J=1,7) /
     S 0.912418292E+05, 0.000000000E-00, 0.925887084E-04,
     S 0.723613782E+05, 0.000000000E-00, 0.129353723E-01,
     S 0.596037057E+04, 0.000000000E-00, 0.800821928E+00,
     S 0.000000000E-00, 0.000000000E-00, 0.242715973E+02,
     S 0.000000000E-00, 0.000000000E-00, 0.878331486E+02,
     S 0.000000000E-00, 0.000000000E-00, 0.191559725E+02,
     S 0.000000000E-00, 0.000000000E-00, 0.000000000E+00 /
      DATA ((APAD(2,I,J),I=1,3),J=1,7) /
     S 0.376655383E-08, 0.739646016E-08, 0.410177786E+03,
     S 0.978576773E-04, 0.131849595E-03, 0.672595424E+02,
     S 0.387714006E+00, 0.437772681E+00, 0.000000000E-00,
     S 0.118461660E+03, 0.151345118E+03, 0.000000000E-00,
     S 0.119079797E+04, 0.233628890E+04, 0.000000000E-00,
     S 0.293353397E+03, 0.797219934E+03, 0.000000000E-00,
     S 0.000000000E+00, 0.000000000E+00, 0.000000000E+00 /
C
      DATA ((BPAD(1,I,J),I=1,3),J=1,7) /
     S 0.912418292E+05, 0.000000000E-00, 0.925887084E-04,
     S 0.724555318E+05, 0.000000000E-00, 0.131812683E-01,
     S 0.602593328E+04, 0.000000000E-00, 0.812706117E+00,
     S 0.100000000E+01, 0.000000000E-00, 0.249863591E+02,
     S 0.000000000E-00, 0.000000000E-00, 0.931071925E+02,
     S 0.000000000E-00, 0.000000000E-00, 0.252233437E+02,
     S 0.000000000E-00, 0.000000000E-00, 0.100000000E+01 /
      DATA ((BPAD(2,I,J),I=1,3),J=1,7) /
     S 0.376655383E-08, 0.739646016E-08, 0.410177786E+03,
     S 0.979023421E-04, 0.131861712E-03, 0.731185438E+02,
     S 0.388611139E+00, 0.437949001E+00, 0.100000000E+01,
     S 0.120291383E+03, 0.151692730E+03, 0.000000000E+00,
     S 0.130531005E+04, 0.237071130E+04, 0.000000000E+00,
     S 0.415049409E+03, 0.867914360E+03, 0.000000000E+00,
     S 0.100000000E+01, 0.100000000E+01, 0.000000000E+00 /
c
      DATA (D(1,I),I=1,3) / 0.00, 0.00, 0.00 /
      DATA (D(2,I),I=1,3) / 0.000000000, 0.000000000, 0.800000000 /
C
C-----------------------------------------------------------------------
C
C*         1.      HORNER'S ALGORITHM TO COMPUTE TRANSMISSION FUNCTION
C
 100  CONTINUE
C
      DO 201 JL = 1, KDLON
      ZR1(JL) = APAD(KNU,KA,1) + PU(JL) * (APAD(KNU,KA,2) + PU(JL)
     S      * ( APAD(KNU,KA,3) + PU(JL) * (APAD(KNU,KA,4) + PU(JL)
     S      * ( APAD(KNU,KA,5) + PU(JL) * (APAD(KNU,KA,6) + PU(JL)
     S      * ( APAD(KNU,KA,7) ))))))
C
      ZR2(JL) = BPAD(KNU,KA,1) + PU(JL) * (BPAD(KNU,KA,2) + PU(JL)
     S      * ( BPAD(KNU,KA,3) + PU(JL) * (BPAD(KNU,KA,4) + PU(JL)
     S      * ( BPAD(KNU,KA,5) + PU(JL) * (BPAD(KNU,KA,6) + PU(JL)
     S      * ( BPAD(KNU,KA,7) ))))))
C     
C
C*         2.      ADD THE BACKGROUND TRANSMISSION
C
 200  CONTINUE
C
C
      PTR(JL) = (ZR1(JL) / ZR2(JL)) * (1. - D(KNU,KA)) + D(KNU,KA)
 201  CONTINUE
C
      RETURN
      END
      SUBROUTINE SWTT1(KNU,KABS,KIND, PU, PTR)
      IMPLICIT none
      include "dimensions.h"
      include "dimphy.h"
      include "raddim.h"
C
C-----------------------------------------------------------------------
C     PURPOSE.
C     --------
C           THIS ROUTINE COMPUTES THE TRANSMISSION FUNCTIONS FOR ALL THE
C     ABSORBERS (H2O, UNIFORMLY MIXED GASES, AND O3) IN THE TWO SPECTRAL
C     INTERVALS.
C
C     METHOD.
C     -------
C
C          TRANSMISSION FUNCTION ARE COMPUTED USING PADE APPROXIMANTS
C     AND HORNER'S ALGORITHM.
C
C     REFERENCE.
C     ----------
C
C        SEE RADIATION'S PART OF THE MODEL'S DOCUMENTATION AND
C        ECMWF RESEARCH DEPARTMENT DOCUMENTATION OF THE IFS
C
C     AUTHOR.
C     -------
C        JEAN-JACQUES MORCRETTE  *ECMWF*
C
C     MODIFICATIONS.
C     --------------
C        ORIGINAL : 95-01-20
C-----------------------------------------------------------------------
C* ARGUMENTS:
C
      INTEGER KNU          ! INDEX OF THE SPECTRAL INTERVAL
      INTEGER KABS         ! NUMBER OF ABSORBERS
      INTEGER KIND(KABS)   ! INDICES OF THE ABSORBERS
      REAL*8 PU(KDLON,KABS)  ! ABSORBER AMOUNT
C
      REAL*8 PTR(KDLON,KABS) ! TRANSMISSION FUNCTION
C
C* LOCAL VARIABLES:
C
      REAL*8 ZR1(KDLON)
      REAL*8 ZR2(KDLON)
      REAL*8 ZU(KDLON)
      INTEGER jl, ja, i, j, ia
C
C* Prescribed Data:
C
      REAL*8 APAD(2,3,7), BPAD(2,3,7), D(2,3)
      SAVE APAD, BPAD, D
      DATA ((APAD(1,I,J),I=1,3),J=1,7) /
     S 0.912418292E+05, 0.000000000E-00, 0.925887084E-04,
     S 0.723613782E+05, 0.000000000E-00, 0.129353723E-01,
     S 0.596037057E+04, 0.000000000E-00, 0.800821928E+00,
     S 0.000000000E-00, 0.000000000E-00, 0.242715973E+02,
     S 0.000000000E-00, 0.000000000E-00, 0.878331486E+02,
     S 0.000000000E-00, 0.000000000E-00, 0.191559725E+02,
     S 0.000000000E-00, 0.000000000E-00, 0.000000000E+00 /
      DATA ((APAD(2,I,J),I=1,3),J=1,7) /
     S 0.376655383E-08, 0.739646016E-08, 0.410177786E+03,
     S 0.978576773E-04, 0.131849595E-03, 0.672595424E+02,
     S 0.387714006E+00, 0.437772681E+00, 0.000000000E-00,
     S 0.118461660E+03, 0.151345118E+03, 0.000000000E-00,
     S 0.119079797E+04, 0.233628890E+04, 0.000000000E-00,
     S 0.293353397E+03, 0.797219934E+03, 0.000000000E-00,
     S 0.000000000E+00, 0.000000000E+00, 0.000000000E+00 /
C
      DATA ((BPAD(1,I,J),I=1,3),J=1,7) /
     S 0.912418292E+05, 0.000000000E-00, 0.925887084E-04,
     S 0.724555318E+05, 0.000000000E-00, 0.131812683E-01,
     S 0.602593328E+04, 0.000000000E-00, 0.812706117E+00,
     S 0.100000000E+01, 0.000000000E-00, 0.249863591E+02,
     S 0.000000000E-00, 0.000000000E-00, 0.931071925E+02,
     S 0.000000000E-00, 0.000000000E-00, 0.252233437E+02,
     S 0.000000000E-00, 0.000000000E-00, 0.100000000E+01 /
      DATA ((BPAD(2,I,J),I=1,3),J=1,7) /
     S 0.376655383E-08, 0.739646016E-08, 0.410177786E+03,
     S 0.979023421E-04, 0.131861712E-03, 0.731185438E+02,
     S 0.388611139E+00, 0.437949001E+00, 0.100000000E+01,
     S 0.120291383E+03, 0.151692730E+03, 0.000000000E+00,
     S 0.130531005E+04, 0.237071130E+04, 0.000000000E+00,
     S 0.415049409E+03, 0.867914360E+03, 0.000000000E+00,
     S 0.100000000E+01, 0.100000000E+01, 0.000000000E+00 /
c
      DATA (D(1,I),I=1,3) / 0.00, 0.00, 0.00 /
      DATA (D(2,I),I=1,3) / 0.000000000, 0.000000000, 0.800000000 /
C-----------------------------------------------------------------------
C
C*         1.      HORNER'S ALGORITHM TO COMPUTE TRANSMISSION FUNCTION
C
 100  CONTINUE
C
      DO 202 JA = 1,KABS
      IA=KIND(JA)
      DO 201 JL = 1, KDLON
      ZU(JL) = PU(JL,JA)
      ZR1(JL) = APAD(KNU,IA,1) + ZU(JL) * (APAD(KNU,IA,2) + ZU(JL)
     S      * ( APAD(KNU,IA,3) + ZU(JL) * (APAD(KNU,IA,4) + ZU(JL)
     S      * ( APAD(KNU,IA,5) + ZU(JL) * (APAD(KNU,IA,6) + ZU(JL)
     S      * ( APAD(KNU,IA,7) ))))))
C
      ZR2(JL) = BPAD(KNU,IA,1) + ZU(JL) * (BPAD(KNU,IA,2) + ZU(JL)
     S      * ( BPAD(KNU,IA,3) + ZU(JL) * (BPAD(KNU,IA,4) + ZU(JL)
     S      * ( BPAD(KNU,IA,5) + ZU(JL) * (BPAD(KNU,IA,6) + ZU(JL)
     S      * ( BPAD(KNU,IA,7) ))))))
C     
C
C*         2.      ADD THE BACKGROUND TRANSMISSION
C
 200  CONTINUE
C
      PTR(JL,JA) = (ZR1(JL)/ZR2(JL)) * (1.-D(KNU,IA)) + D(KNU,IA) 
 201  CONTINUE
 202  CONTINUE
C
      RETURN
      END
      SUBROUTINE LW(RCO2,RCH4,RN2O,RCFC11,RCFC12,
     .              PPMB, PDP,
     .              PPSOL,PDT0,PEMIS,
     .              PTL, PTAVE, PWV, POZON, PAER,
     .              PCLDLD,PCLDLU,
     .              PVIEW,
     .              PCOLR, PCOLR0,
     .              PTOPLW,PSOLLW,PTOPLW0,PSOLLW0)
      IMPLICIT none
      include "dimensions.h"
      include "dimphy.h"
      include "raddim.h"
      include "raddimlw.h"
      include "YOMCST.h"
C
C-----------------------------------------------------------------------
C     METHOD.
C     -------
C
C          1. COMPUTES THE PRESSURE AND TEMPERATURE WEIGHTED AMOUNTS OF
C     ABSORBERS.
C          2. COMPUTES THE PLANCK FUNCTIONS ON THE INTERFACES AND THE
C     GRADIENT OF PLANCK FUNCTIONS IN THE LAYERS.
C          3. PERFORMS THE VERTICAL INTEGRATION DISTINGUISHING THE CON-
C     TRIBUTIONS OF THE ADJACENT AND DISTANT LAYERS AND THOSE FROM THE
C     BOUNDARIES.
C          4. COMPUTES THE CLEAR-SKY DOWNWARD AND UPWARD EMISSIVITIES.
C          5. INTRODUCES THE EFFECTS OF THE CLOUDS ON THE FLUXES.
C
C
C     REFERENCE.
C     ----------
C
C        SEE RADIATION'S PART OF THE MODEL'S DOCUMENTATION AND
C        ECMWF RESEARCH DEPARTMENT DOCUMENTATION OF THE IFS
C
C     AUTHOR.
C     -------
C        JEAN-JACQUES MORCRETTE  *ECMWF*
C
C     MODIFICATIONS.
C     --------------
C        ORIGINAL : 89-07-14
C-----------------------------------------------------------------------
      REAL*8 RCO2   ! CO2 CONCENTRATION (IPCC:353.E-06* 44.011/28.97)
      REAL*8 RCH4   ! CH4 CONCENTRATION (IPCC: 1.72E-06* 16.043/28.97)
      REAL*8 RN2O   ! N2O CONCENTRATION (IPCC: 310.E-09* 44.013/28.97)
      REAL*8 RCFC11 ! CFC11 CONCENTRATION (IPCC: 280.E-12* 137.3686/28.97)
      REAL*8 RCFC12 ! CFC12 CONCENTRATION (IPCC: 484.E-12* 120.9140/28.97)
      REAL*8 PCLDLD(KDLON,KFLEV)  ! DOWNWARD EFFECTIVE CLOUD COVER
      REAL*8 PCLDLU(KDLON,KFLEV)  ! UPWARD EFFECTIVE CLOUD COVER
      REAL*8 PDP(KDLON,KFLEV)     ! LAYER PRESSURE THICKNESS (Pa)
      REAL*8 PDT0(KDLON)          ! SURFACE TEMPERATURE DISCONTINUITY (K)
      REAL*8 PEMIS(KDLON)         ! SURFACE EMISSIVITY
      REAL*8 PPMB(KDLON,KFLEV+1)  ! HALF LEVEL PRESSURE (mb)
      REAL*8 PPSOL(KDLON)         ! SURFACE PRESSURE (Pa)
      REAL*8 POZON(KDLON,KFLEV)   ! O3 CONCENTRATION (kg/kg)
      REAL*8 PTL(KDLON,KFLEV+1)   ! HALF LEVEL TEMPERATURE (K)
      REAL*8 PAER(KDLON,KFLEV,5)  ! OPTICAL THICKNESS OF THE AEROSOLS
      REAL*8 PTAVE(KDLON,KFLEV)   ! LAYER TEMPERATURE (K)
      REAL*8 PVIEW(KDLON)         ! COSECANT OF VIEWING ANGLE
      REAL*8 PWV(KDLON,KFLEV)     ! SPECIFIC HUMIDITY (kg/kg)
C
      REAL*8 PCOLR(KDLON,KFLEV)   ! LONG-WAVE TENDENCY (K/day)
      REAL*8 PCOLR0(KDLON,KFLEV)  ! LONG-WAVE TENDENCY (K/day) clear-sky
      REAL*8 PTOPLW(KDLON)        ! LONGWAVE FLUX AT T.O.A.
      REAL*8 PSOLLW(KDLON)        ! LONGWAVE FLUX AT SURFACE
      REAL*8 PTOPLW0(KDLON)       ! LONGWAVE FLUX AT T.O.A. (CLEAR-SKY)
      REAL*8 PSOLLW0(KDLON)       ! LONGWAVE FLUX AT SURFACE (CLEAR-SKY)
C
C-------------------------------------------------------------------------
      REAL*8 ZABCU(KDLON,NUA,3*KFLEV+1)
      REAL*8 ZOZ(KDLON,KFLEV)
c
      REAL*8 ZFLUX(KDLON,2,KFLEV+1) ! RADIATIVE FLUXES (1:up; 2:down)
      REAL*8 ZFLUC(KDLON,2,KFLEV+1) ! CLEAR-SKY RADIATIVE FLUXES
      REAL*8 ZBINT(KDLON,KFLEV+1)            ! Intermediate variable
      REAL*8 ZBSUI(KDLON)                    ! Intermediate variable
      REAL*8 ZCTS(KDLON,KFLEV)               ! Intermediate variable
      REAL*8 ZCNTRB(KDLON,KFLEV+1,KFLEV+1)   ! Intermediate variable
      SAVE ZFLUX, ZFLUC, ZBINT, ZBSUI, ZCTS, ZCNTRB
c
      INTEGER ilim, i, k, kpl1
C
      INTEGER lw0pas ! Every lw0pas steps, clear-sky is done
      PARAMETER (lw0pas=1)
      INTEGER lwpas  ! Every lwpas steps, cloudy-sky is done
      PARAMETER (lwpas=1)
c
      INTEGER itaplw0, itaplw
      LOGICAL appel1er
      SAVE appel1er, itaplw0, itaplw
      DATA appel1er /.TRUE./
      DATA itaplw0,itaplw /0,0/
C     ------------------------------------------------------------------
      IF (appel1er) THEN
         PRINT*, "LW clear-sky calling frequency: ", lw0pas
         PRINT*, "LW cloudy-sky calling frequency: ", lwpas
         PRINT*, "   In general, they should be 1"
         appel1er=.FALSE.
      ENDIF
C
      IF (MOD(itaplw0,lw0pas).EQ.0) THEN
      DO k = 1, KFLEV  ! convertir ozone de kg/kg en pa/pa
      DO i = 1, KDLON
         ZOZ(i,k) = POZON(i,k)*PDP(i,k) * 28.9644/47.9942
      ENDDO
      ENDDO
      CALL LWU(RCO2,RCH4, RN2O, RCFC11, RCFC12,
     S         PAER,PDP,PPMB,PPSOL,ZOZ,PTAVE,PVIEW,PWV,ZABCU)
      CALL LWBV(ILIM,PDP,PDT0,PEMIS,PPMB,PTL,PTAVE,ZABCU,
     S          ZFLUC,ZBINT,ZBSUI,ZCTS,ZCNTRB)
      itaplw0 = 0
      ENDIF
      itaplw0 = itaplw0 + 1
C
      IF (MOD(itaplw,lwpas).EQ.0) THEN
      CALL LWC(ILIM,PCLDLD,PCLDLU,PEMIS,
     S         ZFLUC,ZBINT,ZBSUI,ZCTS,ZCNTRB,
     S         ZFLUX)
      itaplw = 0
      ENDIF
      itaplw = itaplw + 1
C
      DO k = 1, KFLEV
         kpl1 = k+1
         DO i = 1, KDLON
            PCOLR(i,k) = ZFLUX(i,1,kpl1)+ZFLUX(i,2,kpl1)
     .                 - ZFLUX(i,1,k)-   ZFLUX(i,2,k)
            PCOLR(i,k) = PCOLR(i,k) * RDAY*RG/RCPD / PDP(i,k)
            PCOLR0(i,k) = ZFLUC(i,1,kpl1)+ZFLUC(i,2,kpl1)
     .                 - ZFLUC(i,1,k)-   ZFLUC(i,2,k)
            PCOLR0(i,k) = PCOLR0(i,k) * RDAY*RG/RCPD / PDP(i,k)
         ENDDO
      ENDDO
      DO i = 1, KDLON
         PSOLLW(i) = -ZFLUX(i,1,1)-ZFLUX(i,2,1)
         PTOPLW(i) = ZFLUX(i,1,KFLEV+1) + ZFLUX(i,2,KFLEV+1)
c
         PSOLLW0(i) = -ZFLUC(i,1,1)-ZFLUC(i,2,1)
         PTOPLW0(i) = ZFLUC(i,1,KFLEV+1) + ZFLUC(i,2,KFLEV+1)
      ENDDO
C     ------------------------------------------------------------------
      RETURN
      END
      SUBROUTINE LWU(RCO2, RCH4, RN2O, RCFC11, RCFC12,
     S               PAER,PDP,PPMB,PPSOL,POZ,PTAVE,PVIEW,PWV,
     S               PABCU)
      IMPLICIT none
      include "dimensions.h"
      include "dimphy.h"
      include "raddim.h"
      include "raddimlw.h"
      include "YOMCST.h"
      include "radepsi.h"
      include "radopt.h"
C
C     PURPOSE.
C     --------
C           COMPUTES ABSORBER AMOUNTS INCLUDING PRESSURE AND
C           TEMPERATURE EFFECTS
C
C     METHOD.
C     -------
C
C          1. COMPUTES THE PRESSURE AND TEMPERATURE WEIGHTED AMOUNTS OF
C     ABSORBERS.
C
C
C     REFERENCE.
C     ----------
C
C        SEE RADIATION'S PART OF THE MODEL'S DOCUMENTATION AND
C        ECMWF RESEARCH DEPARTMENT DOCUMENTATION OF THE IFS
C
C     AUTHOR.
C     -------
C        JEAN-JACQUES MORCRETTE  *ECMWF*
C
C     MODIFICATIONS.
C     --------------
C        ORIGINAL : 89-07-14
C        Voigt lines (loop 404 modified) - JJM & PhD - 01/96
C-----------------------------------------------------------------------
C* ARGUMENTS:
      REAL*8 RCO2
      REAL*8 RCH4, RN2O, RCFC11, RCFC12
      REAL*8 PAER(KDLON,KFLEV,5)
      REAL*8 PDP(KDLON,KFLEV)
      REAL*8 PPMB(KDLON,KFLEV+1)
      REAL*8 PPSOL(KDLON)
      REAL*8 POZ(KDLON,KFLEV)
      REAL*8 PTAVE(KDLON,KFLEV)
      REAL*8 PVIEW(KDLON)
      REAL*8 PWV(KDLON,KFLEV)
C
      REAL*8 PABCU(KDLON,NUA,3*KFLEV+1) ! EFFECTIVE ABSORBER AMOUNTS
C
C-----------------------------------------------------------------------
C* LOCAL VARIABLES:
      REAL*8 ZABLY(KDLON,NUA,3*KFLEV+1)
      REAL*8 ZDUC(KDLON,3*KFLEV+1)
      REAL*8 ZPHIO(KDLON)
      REAL*8 ZPSC2(KDLON)
      REAL*8 ZPSC3(KDLON)
      REAL*8 ZPSH1(KDLON)
      REAL*8 ZPSH2(KDLON)
      REAL*8 ZPSH3(KDLON)
      REAL*8 ZPSH4(KDLON)
      REAL*8 ZPSH5(KDLON)
      REAL*8 ZPSH6(KDLON)
      REAL*8 ZPSIO(KDLON)
      REAL*8 ZTCON(KDLON)
      REAL*8 ZPHM6(KDLON)
      REAL*8 ZPSM6(KDLON)
      REAL*8 ZPHN6(KDLON)
      REAL*8 ZPSN6(KDLON)
      REAL*8 ZSSIG(KDLON,3*KFLEV+1)
      REAL*8 ZTAVI(KDLON)
      REAL*8 ZUAER(KDLON,Ninter)
      REAL*8 ZXOZ(KDLON)
      REAL*8 ZXWV(KDLON)
C
      INTEGER jl, jk, jkj, jkjr, jkjp, ig1
      INTEGER jki, jkip1, ja, jj
      INTEGER jkl, jkp1, jkk, jkjpn
      INTEGER jae1, jae2, jae3, jae, jjpn
      INTEGER ir, jc, jcp1
      REAL*8 zdpm, zupm, zupmh2o, zupmco2, zupmo3, zu6, zup
      REAL*8 zfppw, ztx, ztx2, zzably
      REAL*8 zcah1, zcbh1, zcah2, zcbh2, zcah3, zcbh3
      REAL*8 zcah4, zcbh4, zcah5, zcbh5, zcah6, zcbh6
      REAL*8 zcac8, zcbc8
      REAL*8 zalup, zdiff
c
      REAL*8 PVGCO2, PVGH2O, PVGO3
C
      REAL*8 R10E  ! DECIMAL/NATURAL LOG.FACTOR
      PARAMETER (R10E=0.4342945)
c
c Used Data Block:
c
      REAL*8 TREF
      SAVE TREF
      REAL*8 RT1(2)
      SAVE RT1
      REAL*8 RAER(5,5)
      SAVE RAER
      REAL*8 AT(8,3), BT(8,3)
      SAVE AT, BT
      REAL*8 OCT(4)
      SAVE OCT
      DATA TREF /250.0/
      DATA (RT1(IG1),IG1=1,2) / -0.577350269, +0.577350269 /
      DATA RAER / .038520, .037196, .040532, .054934, .038520
     1          , .12613 , .18313 , .10357 , .064106, .126130
     2          , .012579, .013649, .018652, .025181, .012579
     3          , .011890, .016142, .021105, .028908, .011890
     4          , .013792, .026810, .052203, .066338, .013792 /
      DATA (AT(1,IR),IR=1,3) /
     S 0.298199E-02,-.394023E-03,0.319566E-04 /
      DATA (BT(1,IR),IR=1,3) /
     S-0.106432E-04,0.660324E-06,0.174356E-06 /
      DATA (AT(2,IR),IR=1,3) /
     S 0.143676E-01,0.366501E-02,-.160822E-02 /
      DATA (BT(2,IR),IR=1,3) /
     S-0.553979E-04,-.101701E-04,0.920868E-05 /
      DATA (AT(3,IR),IR=1,3) /
     S 0.197861E-01,0.315541E-02,-.174547E-02 /
      DATA (BT(3,IR),IR=1,3) /
     S-0.877012E-04,0.513302E-04,0.523138E-06 /
      DATA (AT(4,IR),IR=1,3) /
     S 0.289560E-01,-.208807E-02,-.121943E-02 /
      DATA (BT(4,IR),IR=1,3) /
     S-0.165960E-03,0.157704E-03,-.146427E-04 /
      DATA (AT(5,IR),IR=1,3) /
     S 0.103800E-01,0.436296E-02,-.161431E-02 /
      DATA (BT(5,IR),IR=1,3) /
     S -.276744E-04,-.327381E-04,0.127646E-04 /
      DATA (AT(6,IR),IR=1,3) /
     S 0.868859E-02,-.972752E-03,0.000000E-00 /
      DATA (BT(6,IR),IR=1,3) /
     S -.278412E-04,-.713940E-06,0.117469E-05 /
      DATA (AT(7,IR),IR=1,3) /
     S 0.250073E-03,0.455875E-03,0.109242E-03 /
      DATA (BT(7,IR),IR=1,3) /
     S 0.199846E-05,-.216313E-05,0.175991E-06 /
      DATA (AT(8,IR),IR=1,3) /
     S 0.307423E-01,0.110879E-02,-.322172E-03 /
      DATA (BT(8,IR),IR=1,3) /
     S-0.108482E-03,0.258096E-05,-.814575E-06 /
c
      DATA OCT /-.326E-03, -.102E-05, .137E-02, -.535E-05/
C-----------------------------------------------------------------------
c
      IF (LEVOIGT) THEN
         PVGCO2= 60.
         PVGH2O= 30.
         PVGO3 =400.
      ELSE
         PVGCO2= 0.
         PVGH2O= 0.
         PVGO3 = 0.
      ENDIF
C
C
C*         2.    PRESSURE OVER GAUSS SUB-LEVELS
C                ------------------------------
C
 200  CONTINUE
C
      DO 201 JL = 1, KDLON
      ZSSIG(JL, 1 ) = PPMB(JL,1) * 100.
 201  CONTINUE
C
      DO 206 JK = 1 , KFLEV
      JKJ=(JK-1)*NG1P1+1
      JKJR = JKJ
      JKJP = JKJ + NG1P1
      DO 203 JL = 1, KDLON
      ZSSIG(JL,JKJP)=PPMB(JL,JK+1)* 100.
 203  CONTINUE
      DO 205 IG1=1,NG1
      JKJ=JKJ+1
      DO 204 JL = 1, KDLON
      ZSSIG(JL,JKJ)= (ZSSIG(JL,JKJR)+ZSSIG(JL,JKJP))*0.5
     S  + RT1(IG1) * (ZSSIG(JL,JKJP) - ZSSIG(JL,JKJR)) * 0.5
 204  CONTINUE
 205  CONTINUE
 206  CONTINUE
C
C-----------------------------------------------------------------------
C
C
C*         4.    PRESSURE THICKNESS AND MEAN PRESSURE OF SUB-LAYERS
C                --------------------------------------------------
C
 400  CONTINUE
C
      DO 402 JKI=1,3*KFLEV
      JKIP1=JKI+1
      DO 401 JL = 1, KDLON
      ZABLY(JL,5,JKI)=(ZSSIG(JL,JKI)+ZSSIG(JL,JKIP1))*0.5
      ZABLY(JL,3,JKI)=(ZSSIG(JL,JKI)-ZSSIG(JL,JKIP1))
     S                                 /(10.*RG)
 401  CONTINUE
 402  CONTINUE
C
      DO 406 JK = 1 , KFLEV
      JKP1=JK+1
      JKL = KFLEV+1 - JK
      DO 403 JL = 1, KDLON
      ZXWV(JL) = MAX (PWV(JL,JK) , ZEPSCQ )
      ZXOZ(JL) = MAX (POZ(JL,JK) / PDP(JL,JK) , ZEPSCO )
 403  CONTINUE
      JKJ=(JK-1)*NG1P1+1
      JKJPN=JKJ+NG1
      DO 405 JKK=JKJ,JKJPN
      DO 404 JL = 1, KDLON
      ZDPM = ZABLY(JL,3,JKK)
      ZUPM = ZABLY(JL,5,JKK)             * ZDPM / 101325.
      ZUPMCO2 = ( ZABLY(JL,5,JKK) + PVGCO2 ) * ZDPM / 101325.
      ZUPMH2O = ( ZABLY(JL,5,JKK) + PVGH2O ) * ZDPM / 101325.
      ZUPMO3  = ( ZABLY(JL,5,JKK) + PVGO3  ) * ZDPM / 101325.
      ZDUC(JL,JKK) = ZDPM
      ZABLY(JL,12,JKK) = ZXOZ(JL) * ZDPM
      ZABLY(JL,13,JKK) = ZXOZ(JL) * ZUPMO3
      ZU6 = ZXWV(JL) * ZUPM
      ZFPPW = 1.6078 * ZXWV(JL) / (1.+0.608*ZXWV(JL))
      ZABLY(JL,6,JKK) = ZXWV(JL) * ZUPMH2O
      ZABLY(JL,11,JKK) = ZU6 * ZFPPW
      ZABLY(JL,10,JKK) = ZU6 * (1.-ZFPPW)
      ZABLY(JL,9,JKK) = RCO2 * ZUPMCO2
      ZABLY(JL,8,JKK) = RCO2 * ZDPM
 404  CONTINUE
 405  CONTINUE
 406  CONTINUE
C
C-----------------------------------------------------------------------
C
C
C*         5.    CUMULATIVE ABSORBER AMOUNTS FROM TOP OF ATMOSPHERE
C                --------------------------------------------------
C
 500  CONTINUE
C
      DO 502 JA = 1, NUA
      DO 501 JL = 1, KDLON
      PABCU(JL,JA,3*KFLEV+1) = 0.
  501 CONTINUE
  502 CONTINUE
C
      DO 529 JK = 1 , KFLEV
      JJ=(JK-1)*NG1P1+1
      JJPN=JJ+NG1
      JKL=KFLEV+1-JK
C
C
C*         5.1  CUMULATIVE AEROSOL AMOUNTS FROM TOP OF ATMOSPHERE
C               --------------------------------------------------
C
 510  CONTINUE
C
      JAE1=3*KFLEV+1-JJ
      JAE2=3*KFLEV+1-(JJ+1)
      JAE3=3*KFLEV+1-JJPN
      DO 512 JAE=1,5
      DO 511 JL = 1, KDLON
      ZUAER(JL,JAE) = (RAER(JAE,1)*PAER(JL,JKL,1)
     S      +RAER(JAE,2)*PAER(JL,JKL,2)+RAER(JAE,3)*PAER(JL,JKL,3)
     S      +RAER(JAE,4)*PAER(JL,JKL,4)+RAER(JAE,5)*PAER(JL,JKL,5))
     S      /(ZDUC(JL,JAE1)+ZDUC(JL,JAE2)+ZDUC(JL,JAE3))
 511  CONTINUE
 512  CONTINUE
C
C
C
C*         5.2  INTRODUCES TEMPERATURE EFFECTS ON ABSORBER AMOUNTS
C               --------------------------------------------------
C
 520  CONTINUE
C
      DO 521 JL = 1, KDLON
      ZTAVI(JL)=PTAVE(JL,JKL)
      ZTCON(JL)=EXP(6.08*(296./ZTAVI(JL)-1.))
      ZTX=ZTAVI(JL)-TREF
      ZTX2=ZTX*ZTX
      ZZABLY = ZABLY(JL,6,JAE1)+ZABLY(JL,6,JAE2)+ZABLY(JL,6,JAE3)
      ZUP=MIN( MAX( 0.5*R10E*LOG( ZZABLY ) + 5., 0.), 6.0)
      ZCAH1=AT(1,1)+ZUP*(AT(1,2)+ZUP*(AT(1,3)))
      ZCBH1=BT(1,1)+ZUP*(BT(1,2)+ZUP*(BT(1,3)))
      ZPSH1(JL)=EXP( ZCAH1 * ZTX + ZCBH1 * ZTX2 )
      ZCAH2=AT(2,1)+ZUP*(AT(2,2)+ZUP*(AT(2,3)))
      ZCBH2=BT(2,1)+ZUP*(BT(2,2)+ZUP*(BT(2,3)))
      ZPSH2(JL)=EXP( ZCAH2 * ZTX + ZCBH2 * ZTX2 )
      ZCAH3=AT(3,1)+ZUP*(AT(3,2)+ZUP*(AT(3,3)))
      ZCBH3=BT(3,1)+ZUP*(BT(3,2)+ZUP*(BT(3,3)))
      ZPSH3(JL)=EXP( ZCAH3 * ZTX + ZCBH3 * ZTX2 )
      ZCAH4=AT(4,1)+ZUP*(AT(4,2)+ZUP*(AT(4,3)))
      ZCBH4=BT(4,1)+ZUP*(BT(4,2)+ZUP*(BT(4,3)))
      ZPSH4(JL)=EXP( ZCAH4 * ZTX + ZCBH4 * ZTX2 )
      ZCAH5=AT(5,1)+ZUP*(AT(5,2)+ZUP*(AT(5,3)))
      ZCBH5=BT(5,1)+ZUP*(BT(5,2)+ZUP*(BT(5,3)))
      ZPSH5(JL)=EXP( ZCAH5 * ZTX + ZCBH5 * ZTX2 )
      ZCAH6=AT(6,1)+ZUP*(AT(6,2)+ZUP*(AT(6,3)))
      ZCBH6=BT(6,1)+ZUP*(BT(6,2)+ZUP*(BT(6,3)))
      ZPSH6(JL)=EXP( ZCAH6 * ZTX + ZCBH6 * ZTX2 )
      ZPHM6(JL)=EXP(-5.81E-4 * ZTX - 1.13E-6 * ZTX2 )
      ZPSM6(JL)=EXP(-5.57E-4 * ZTX - 3.30E-6 * ZTX2 )
      ZPHN6(JL)=EXP(-3.46E-5 * ZTX + 2.05E-7 * ZTX2 )
      ZPSN6(JL)=EXP( 3.70E-3 * ZTX - 2.30E-6 * ZTX2 )
 521  CONTINUE
C
      DO 522 JL = 1, KDLON
      ZTAVI(JL)=PTAVE(JL,JKL)
      ZTX=ZTAVI(JL)-TREF
      ZTX2=ZTX*ZTX
      ZZABLY = ZABLY(JL,9,JAE1)+ZABLY(JL,9,JAE2)+ZABLY(JL,9,JAE3)
      ZALUP = R10E * LOG ( ZZABLY )
      ZUP   = MAX( 0.0 , 5.0 + 0.5 * ZALUP )
      ZPSC2(JL) = (ZTAVI(JL)/TREF) ** ZUP
      ZCAC8=AT(8,1)+ZUP*(AT(8,2)+ZUP*(AT(8,3)))
      ZCBC8=BT(8,1)+ZUP*(BT(8,2)+ZUP*(BT(8,3)))
      ZPSC3(JL)=EXP( ZCAC8 * ZTX + ZCBC8 * ZTX2 )
      ZPHIO(JL) = EXP( OCT(1) * ZTX + OCT(2) * ZTX2)
      ZPSIO(JL) = EXP( 2.* (OCT(3)*ZTX+OCT(4)*ZTX2))
 522  CONTINUE
C
      DO 524 JKK=JJ,JJPN
      JC=3*KFLEV+1-JKK
      JCP1=JC+1
      DO 523 JL = 1, KDLON
      ZDIFF = PVIEW(JL)
      PABCU(JL,10,JC)=PABCU(JL,10,JCP1)
     S                +ZABLY(JL,10,JC)           *ZDIFF
      PABCU(JL,11,JC)=PABCU(JL,11,JCP1)
     S                +ZABLY(JL,11,JC)*ZTCON(JL)*ZDIFF
C
      PABCU(JL,12,JC)=PABCU(JL,12,JCP1)
     S                +ZABLY(JL,12,JC)*ZPHIO(JL)*ZDIFF
      PABCU(JL,13,JC)=PABCU(JL,13,JCP1)
     S                +ZABLY(JL,13,JC)*ZPSIO(JL)*ZDIFF
C
      PABCU(JL,7,JC)=PABCU(JL,7,JCP1)
     S               +ZABLY(JL,9,JC)*ZPSC2(JL)*ZDIFF
      PABCU(JL,8,JC)=PABCU(JL,8,JCP1)
     S               +ZABLY(JL,9,JC)*ZPSC3(JL)*ZDIFF
      PABCU(JL,9,JC)=PABCU(JL,9,JCP1)
     S               +ZABLY(JL,9,JC)*ZPSC3(JL)*ZDIFF
C
      PABCU(JL,1,JC)=PABCU(JL,1,JCP1)
     S               +ZABLY(JL,6,JC)*ZPSH1(JL)*ZDIFF
      PABCU(JL,2,JC)=PABCU(JL,2,JCP1)
     S               +ZABLY(JL,6,JC)*ZPSH2(JL)*ZDIFF
      PABCU(JL,3,JC)=PABCU(JL,3,JCP1)
     S               +ZABLY(JL,6,JC)*ZPSH5(JL)*ZDIFF
      PABCU(JL,4,JC)=PABCU(JL,4,JCP1)
     S               +ZABLY(JL,6,JC)*ZPSH3(JL)*ZDIFF
      PABCU(JL,5,JC)=PABCU(JL,5,JCP1)
     S               +ZABLY(JL,6,JC)*ZPSH4(JL)*ZDIFF
      PABCU(JL,6,JC)=PABCU(JL,6,JCP1)
     S               +ZABLY(JL,6,JC)*ZPSH6(JL)*ZDIFF
C
      PABCU(JL,14,JC)=PABCU(JL,14,JCP1)
     S                +ZUAER(JL,1)    *ZDUC(JL,JC)*ZDIFF
      PABCU(JL,15,JC)=PABCU(JL,15,JCP1)
     S                +ZUAER(JL,2)    *ZDUC(JL,JC)*ZDIFF
      PABCU(JL,16,JC)=PABCU(JL,16,JCP1)
     S                +ZUAER(JL,3)    *ZDUC(JL,JC)*ZDIFF
      PABCU(JL,17,JC)=PABCU(JL,17,JCP1)
     S                +ZUAER(JL,4)    *ZDUC(JL,JC)*ZDIFF
      PABCU(JL,18,JC)=PABCU(JL,18,JCP1)
     S                +ZUAER(JL,5)    *ZDUC(JL,JC)*ZDIFF
C
      PABCU(JL,19,JC)=PABCU(JL,19,JCP1)
     S               +ZABLY(JL,8,JC)*RCH4/RCO2*ZPHM6(JL)*ZDIFF
      PABCU(JL,20,JC)=PABCU(JL,20,JCP1)
     S               +ZABLY(JL,9,JC)*RCH4/RCO2*ZPSM6(JL)*ZDIFF
      PABCU(JL,21,JC)=PABCU(JL,21,JCP1)
     S               +ZABLY(JL,8,JC)*RN2O/RCO2*ZPHN6(JL)*ZDIFF
      PABCU(JL,22,JC)=PABCU(JL,22,JCP1)
     S               +ZABLY(JL,9,JC)*RN2O/RCO2*ZPSN6(JL)*ZDIFF
C
      PABCU(JL,23,JC)=PABCU(JL,23,JCP1)
     S               +ZABLY(JL,8,JC)*RCFC11/RCO2         *ZDIFF
      PABCU(JL,24,JC)=PABCU(JL,24,JCP1)
     S               +ZABLY(JL,8,JC)*RCFC12/RCO2         *ZDIFF
 523  CONTINUE
 524  CONTINUE
C
 529  CONTINUE
C
C
      RETURN
      END
      SUBROUTINE LWBV(KLIM,PDP,PDT0,PEMIS,PPMB,PTL,PTAVE,PABCU,
     S                PFLUC,PBINT,PBSUI,PCTS,PCNTRB)
      IMPLICIT none
      include "dimensions.h"
      include "dimphy.h"
      include "raddim.h"
      include "raddimlw.h"
      include "YOMCST.h"
C
C     PURPOSE.
C     --------
C           TO COMPUTE THE PLANCK FUNCTION AND PERFORM THE
C           VERTICAL INTEGRATION. SPLIT OUT FROM LW FOR MEMORY
C           SAVING
C
C     METHOD.
C     -------
C
C          1. COMPUTES THE PLANCK FUNCTIONS ON THE INTERFACES AND THE
C     GRADIENT OF PLANCK FUNCTIONS IN THE LAYERS.
C          2. PERFORMS THE VERTICAL INTEGRATION DISTINGUISHING THE CON-
C     TRIBUTIONS OF THE ADJACENT AND DISTANT LAYERS AND THOSE FROM THE
C     BOUNDARIES.
C          3. COMPUTES THE CLEAR-SKY COOLING RATES.
C
C     REFERENCE.
C     ----------
C
C        SEE RADIATION'S PART OF THE MODEL'S DOCUMENTATION AND
C        ECMWF RESEARCH DEPARTMENT DOCUMENTATION OF THE IFS
C
C     AUTHOR.
C     -------
C        JEAN-JACQUES MORCRETTE  *ECMWF*
C
C     MODIFICATIONS.
C     --------------
C        ORIGINAL : 89-07-14
C        MODIFICATION : 93-10-15 M.HAMRUD (SPLIT OUT FROM LW TO SAVE
C                                          MEMORY)
C-----------------------------------------------------------------------
C* ARGUMENTS:
      INTEGER KLIM
C
      REAL*8 PDP(KDLON,KFLEV)
      REAL*8 PDT0(KDLON)
      REAL*8 PEMIS(KDLON)
      REAL*8 PPMB(KDLON,KFLEV+1)
      REAL*8 PTL(KDLON,KFLEV+1)
      REAL*8 PTAVE(KDLON,KFLEV)
C
      REAL*8 PFLUC(KDLON,2,KFLEV+1)
C     
      REAL*8 PABCU(KDLON,NUA,3*KFLEV+1)
      REAL*8 PBINT(KDLON,KFLEV+1)
      REAL*8 PBSUI(KDLON)
      REAL*8 PCTS(KDLON,KFLEV)
      REAL*8 PCNTRB(KDLON,KFLEV+1,KFLEV+1)
C
C-------------------------------------------------------------------------
C
C* LOCAL VARIABLES:
      REAL*8 ZB(KDLON,Ninter,KFLEV+1)
      REAL*8 ZBSUR(KDLON,Ninter)
      REAL*8 ZBTOP(KDLON,Ninter)
      REAL*8 ZDBSL(KDLON,Ninter,KFLEV*2)
      REAL*8 ZGA(KDLON,8,2,KFLEV)
      REAL*8 ZGB(KDLON,8,2,KFLEV)
      REAL*8 ZGASUR(KDLON,8,2)
      REAL*8 ZGBSUR(KDLON,8,2)
      REAL*8 ZGATOP(KDLON,8,2)
      REAL*8 ZGBTOP(KDLON,8,2)
C
      INTEGER nuaer, ntraer
C     ------------------------------------------------------------------
C* COMPUTES PLANCK FUNCTIONS:
       CALL LWB(PDT0,PTAVE,PTL,
     S          ZB,PBINT,PBSUI,ZBSUR,ZBTOP,ZDBSL,
     S          ZGA,ZGB,ZGASUR,ZGBSUR,ZGATOP,ZGBTOP)
C     ------------------------------------------------------------------
C* PERFORMS THE VERTICAL INTEGRATION:
      NUAER = NUA
      NTRAER = NTRA
      CALL LWV(NUAER,NTRAER, KLIM
     R  , PABCU,ZB,PBINT,PBSUI,ZBSUR,ZBTOP,ZDBSL,PEMIS,PPMB,PTAVE
     R  , ZGA,ZGB,ZGASUR,ZGBSUR,ZGATOP,ZGBTOP
     S  , PCNTRB,PCTS,PFLUC)
C     ------------------------------------------------------------------
      RETURN
      END
      SUBROUTINE LWC(KLIM,PCLDLD,PCLDLU,PEMIS,PFLUC,
     R               PBINT,PBSUIN,PCTS,PCNTRB,
     S               PFLUX)
      IMPLICIT none
      include "dimensions.h"
      include "dimphy.h"
      include "raddim.h"
      include "radepsi.h"
      include "radopt.h"
C
C     PURPOSE.
C     --------
C           INTRODUCES CLOUD EFFECTS ON LONGWAVE FLUXES OR
C           RADIANCES
C
C        EXPLICIT ARGUMENTS :
C        --------------------
C     ==== INPUTS ===
C PBINT  : (KDLON,0:KFLEV)     ; HALF LEVEL PLANCK FUNCTION
C PBSUIN : (KDLON)             ; SURFACE PLANCK FUNCTION
C PCLDLD : (KDLON,KFLEV)       ; DOWNWARD EFFECTIVE CLOUD FRACTION
C PCLDLU : (KDLON,KFLEV)       ; UPWARD EFFECTIVE CLOUD FRACTION
C PCNTRB : (KDLON,KFLEV+1,KFLEV+1); CLEAR-SKY ENERGY EXCHANGE
C PCTS   : (KDLON,KFLEV)       ; CLEAR-SKY LAYER COOLING-TO-SPACE
C PEMIS  : (KDLON)             ; SURFACE EMISSIVITY
C PFLUC
C     ==== OUTPUTS ===
C PFLUX(KDLON,2,KFLEV)         ; RADIATIVE FLUXES :
C                     1  ==>  UPWARD   FLUX TOTAL
C                     2  ==>  DOWNWARD FLUX TOTAL
C
C     METHOD.
C     -------
C
C          1. INITIALIZES ALL FLUXES TO CLEAR-SKY VALUES
C          2. EFFECT OF ONE OVERCAST UNITY EMISSIVITY CLOUD LAYER
C          3. EFFECT OF SEMI-TRANSPARENT, PARTIAL OR MULTI-LAYERED
C     CLOUDS
C
C     REFERENCE.
C     ----------
C
C        SEE RADIATION'S PART OF THE MODEL'S DOCUMENTATION AND
C        ECMWF RESEARCH DEPARTMENT DOCUMENTATION OF THE IFS
C
C     AUTHOR.
C     -------
C        JEAN-JACQUES MORCRETTE  *ECMWF*
C
C     MODIFICATIONS.
C     --------------
C        ORIGINAL : 89-07-14
C        Voigt lines (loop 231 to 233)  - JJM & PhD - 01/96
C-----------------------------------------------------------------------
C* ARGUMENTS:
      INTEGER klim
      REAL*8 PFLUC(KDLON,2,KFLEV+1) ! CLEAR-SKY RADIATIVE FLUXES
      REAL*8 PBINT(KDLON,KFLEV+1)   ! HALF LEVEL PLANCK FUNCTION
      REAL*8 PBSUIN(KDLON)          ! SURFACE PLANCK FUNCTION
      REAL*8 PCNTRB(KDLON,KFLEV+1,KFLEV+1) !CLEAR-SKY ENERGY EXCHANGE
      REAL*8 PCTS(KDLON,KFLEV)      ! CLEAR-SKY LAYER COOLING-TO-SPACE
c
      REAL*8 PCLDLD(KDLON,KFLEV)
      REAL*8 PCLDLU(KDLON,KFLEV)
      REAL*8 PEMIS(KDLON)
C
      REAL*8 PFLUX(KDLON,2,KFLEV+1)
C-----------------------------------------------------------------------
C* LOCAL VARIABLES:
      INTEGER IMX(KDLON), IMXP(KDLON)
C
      REAL*8 ZCLEAR(KDLON),ZCLOUD(KDLON),ZDNF(KDLON,KFLEV+1,KFLEV+1)
     S  , ZFD(KDLON), ZFN10(KDLON), ZFU(KDLON)
     S  , ZUPF(KDLON,KFLEV+1,KFLEV+1)
      REAL*8 ZCLM(KDLON,KFLEV+1,KFLEV+1)
C
      INTEGER jk, jl, imaxc, imx1, imx2, jkj, jkp1, jkm1
      INTEGER jk1, jk2, jkc, jkcp1, jcloud
      INTEGER imxm1, imxp1
      REAL*8 zcfrac
C     ------------------------------------------------------------------
C
C*         1.     INITIALIZATION
C                 --------------
C
 100  CONTINUE
C
      IMAXC = 0
C
      DO 101 JL = 1, KDLON
      IMX(JL)=0
      IMXP(JL)=0
      ZCLOUD(JL) = 0.
 101  CONTINUE
C
C*         1.1    SEARCH THE LAYER INDEX OF THE HIGHEST CLOUD
C                 -------------------------------------------
C
 110  CONTINUE
C
      DO 112 JK = 1 , KFLEV
      DO 111 JL = 1, KDLON
      IMX1=IMX(JL)
      IMX2=JK
      IF (PCLDLU(JL,JK).GT.ZEPSC) THEN
         IMXP(JL)=IMX2
      ELSE
         IMXP(JL)=IMX1
      END IF
      IMAXC=MAX(IMXP(JL),IMAXC)
      IMX(JL)=IMXP(JL)
 111  CONTINUE
 112  CONTINUE
CGM*******
      IMAXC=KFLEV
CGM*******
C
      DO 114 JK = 1 , KFLEV+1
      DO 113 JL = 1, KDLON
      PFLUX(JL,1,JK) = PFLUC(JL,1,JK)
      PFLUX(JL,2,JK) = PFLUC(JL,2,JK)
 113  CONTINUE
 114  CONTINUE
C
C     ------------------------------------------------------------------
C
C*         2.      EFFECT OF CLOUDINESS ON LONGWAVE FLUXES
C                  ---------------------------------------
C
      IF (IMAXC.GT.0) THEN
C
         IMXP1 = IMAXC + 1
         IMXM1 = IMAXC - 1
C
C*         2.0     INITIALIZE TO CLEAR-SKY FLUXES
C                  ------------------------------
C
 200  CONTINUE
C
         DO 203 JK1=1,KFLEV+1
         DO 202 JK2=1,KFLEV+1
         DO 201 JL = 1, KDLON
         ZUPF(JL,JK2,JK1)=PFLUC(JL,1,JK1)
         ZDNF(JL,JK2,JK1)=PFLUC(JL,2,JK1)
 201     CONTINUE
 202     CONTINUE
 203     CONTINUE
C
C*         2.1     FLUXES FOR ONE OVERCAST UNITY EMISSIVITY CLOUD
C                  ----------------------------------------------
C
 210  CONTINUE
C
         DO 213 JKC = 1 , IMAXC
         JCLOUD=JKC
         JKCP1=JCLOUD+1
C
C*         2.1.1   ABOVE THE CLOUD
C                  ---------------
C
 2110 CONTINUE
C
         DO 2115 JK=JKCP1,KFLEV+1
         JKM1=JK-1
         DO 2111 JL = 1, KDLON
         ZFU(JL)=0.
 2111    CONTINUE
         IF (JK .GT. JKCP1) THEN
            DO 2113 JKJ=JKCP1,JKM1
            DO 2112 JL = 1, KDLON
            ZFU(JL) = ZFU(JL) + PCNTRB(JL,JK,JKJ)
 2112       CONTINUE
 2113       CONTINUE
         END IF
C
         DO 2114 JL = 1, KDLON
         ZUPF(JL,JKCP1,JK)=PBINT(JL,JK)-ZFU(JL)
 2114    CONTINUE
 2115    CONTINUE
C
C*         2.1.2   BELOW THE CLOUD
C                  ---------------
C
 2120 CONTINUE
C
         DO 2125 JK=1,JCLOUD
         JKP1=JK+1
         DO 2121 JL = 1, KDLON
         ZFD(JL)=0.
 2121    CONTINUE
C
         IF (JK .LT. JCLOUD) THEN
            DO 2123 JKJ=JKP1,JCLOUD
            DO 2122 JL = 1, KDLON
            ZFD(JL) = ZFD(JL) + PCNTRB(JL,JK,JKJ)
 2122       CONTINUE
 2123       CONTINUE
         END IF
         DO 2124 JL = 1, KDLON
         ZDNF(JL,JKCP1,JK)=-PBINT(JL,JK)-ZFD(JL)
 2124    CONTINUE
 2125    CONTINUE
C
 213     CONTINUE
C
C
C*         2.2     CLOUD COVER MATRIX
C                  ------------------
C
C*    ZCLM(JK1,JK2) IS THE OBSCURATION FACTOR BY CLOUD LAYERS BETWEEN
C     HALF-LEVELS JK1 AND JK2 AS SEEN FROM JK1
C
 220  CONTINUE
C
      DO 223 JK1 = 1 , KFLEV+1
      DO 222 JK2 = 1 , KFLEV+1
      DO 221 JL = 1, KDLON
      ZCLM(JL,JK1,JK2) = 0.
 221  CONTINUE
 222  CONTINUE
 223  CONTINUE
C
C
C
C*         2.4     CLOUD COVER BELOW THE LEVEL OF CALCULATION
C                  ------------------------------------------
C
 240  CONTINUE
C
      DO 244 JK1 = 2 , KFLEV+1
      DO 241 JL = 1, KDLON
      ZCLEAR(JL)=1.
      ZCLOUD(JL)=0.
 241  CONTINUE
      DO 243 JK = JK1 - 1 , 1 , -1
      DO 242 JL = 1, KDLON
      IF (NOVLP.EQ.1) THEN
c* maximum-random       
         ZCLEAR(JL)=ZCLEAR(JL)*(1.0-MAX(PCLDLU(JL,JK),ZCLOUD(JL)))
     *                        /(1.0-MIN(ZCLOUD(JL),1.-ZEPSEC))
         ZCLM(JL,JK1,JK) = 1.0 - ZCLEAR(JL)
         ZCLOUD(JL) = PCLDLU(JL,JK)
      ELSE IF (NOVLP.EQ.2) THEN 
c* maximum      
         ZCLOUD(JL) = MAX(ZCLOUD(JL) , PCLDLU(JL,JK))
         ZCLM(JL,JK1,JK) = ZCLOUD(JL)
      ELSE IF (NOVLP.EQ.3) THEN
c* random      
         ZCLEAR(JL) = ZCLEAR(JL)*(1.0 - PCLDLU(JL,JK))
         ZCLOUD(JL) = 1.0 - ZCLEAR(JL)
         ZCLM(JL,JK1,JK) = ZCLOUD(JL)
      END IF
 242  CONTINUE
 243  CONTINUE
 244  CONTINUE
C
C
C*         2.5     CLOUD COVER ABOVE THE LEVEL OF CALCULATION
C                  ------------------------------------------
C
 250  CONTINUE
C
      DO 254 JK1 = 1 , KFLEV
      DO 251 JL = 1, KDLON
      ZCLEAR(JL)=1.
      ZCLOUD(JL)=0.
 251  CONTINUE
      DO 253 JK = JK1 , KFLEV
      DO 252 JL = 1, KDLON
      IF (NOVLP.EQ.1) THEN
c* maximum-random       
         ZCLEAR(JL)=ZCLEAR(JL)*(1.0-MAX(PCLDLD(JL,JK),ZCLOUD(JL)))
     *                        /(1.0-MIN(ZCLOUD(JL),1.-ZEPSEC))
         ZCLM(JL,JK1,JK) = 1.0 - ZCLEAR(JL)
         ZCLOUD(JL) = PCLDLD(JL,JK)
      ELSE IF (NOVLP.EQ.2) THEN 
c* maximum      
         ZCLOUD(JL) = MAX(ZCLOUD(JL) , PCLDLD(JL,JK))
         ZCLM(JL,JK1,JK) = ZCLOUD(JL)
      ELSE IF (NOVLP.EQ.3) THEN
c* random      
         ZCLEAR(JL) = ZCLEAR(JL)*(1.0 - PCLDLD(JL,JK))
         ZCLOUD(JL) = 1.0 - ZCLEAR(JL)
         ZCLM(JL,JK1,JK) = ZCLOUD(JL)
      END IF
 252  CONTINUE
 253  CONTINUE
 254  CONTINUE
C
C
C
C*         3.      FLUXES FOR PARTIAL/MULTIPLE LAYERED CLOUDINESS
C                  ----------------------------------------------
C
 300  CONTINUE
C
C*         3.1     DOWNWARD FLUXES
C                  ---------------
C
 310  CONTINUE
C
      DO 311 JL = 1, KDLON
      PFLUX(JL,2,KFLEV+1) = 0.
 311  CONTINUE
C
      DO 317 JK1 = KFLEV , 1 , -1
C
C*                 CONTRIBUTION FROM CLEAR-SKY FRACTION
C
      DO 312 JL = 1, KDLON
      ZFD (JL) = (1. - ZCLM(JL,JK1,KFLEV)) * ZDNF(JL,1,JK1)
 312  CONTINUE
C
C*                 CONTRIBUTION FROM ADJACENT CLOUD
C
      DO 313 JL = 1, KDLON
      ZFD(JL) = ZFD(JL) + ZCLM(JL,JK1,JK1) * ZDNF(JL,JK1+1,JK1)
 313  CONTINUE
C
C*                 CONTRIBUTION FROM OTHER CLOUDY FRACTIONS
C
      DO 315 JK = KFLEV-1 , JK1 , -1
      DO 314 JL = 1, KDLON
      ZCFRAC = ZCLM(JL,JK1,JK+1) - ZCLM(JL,JK1,JK)
      ZFD(JL) =  ZFD(JL) + ZCFRAC * ZDNF(JL,JK+2,JK1)
 314  CONTINUE
 315  CONTINUE
C
      DO 316 JL = 1, KDLON
      PFLUX(JL,2,JK1) = ZFD (JL)
 316  CONTINUE
C
 317  CONTINUE
C
C
C
C
C*         3.2     UPWARD FLUX AT THE SURFACE
C                  --------------------------
C
 320  CONTINUE
C
      DO 321 JL = 1, KDLON
      PFLUX(JL,1,1) = PEMIS(JL)*PBSUIN(JL)-(1.-PEMIS(JL))*PFLUX(JL,2,1)
 321  CONTINUE
C
C
C
C*         3.3     UPWARD FLUXES
C                  -------------
C
 330  CONTINUE
C
      DO 337 JK1 = 2 , KFLEV+1
C
C*                 CONTRIBUTION FROM CLEAR-SKY FRACTION
C
      DO 332 JL = 1, KDLON
      ZFU (JL) = (1. - ZCLM(JL,JK1,1)) * ZUPF(JL,1,JK1)
 332  CONTINUE
C
C*                 CONTRIBUTION FROM ADJACENT CLOUD
C
      DO 333 JL = 1, KDLON
      ZFU(JL) =  ZFU(JL) + ZCLM(JL,JK1,JK1-1) * ZUPF(JL,JK1,JK1)
 333  CONTINUE
C
C*                 CONTRIBUTION FROM OTHER CLOUDY FRACTIONS
C
      DO 335 JK = 2 , JK1-1
      DO 334 JL = 1, KDLON
      ZCFRAC = ZCLM(JL,JK1,JK-1) - ZCLM(JL,JK1,JK)
      ZFU(JL) =  ZFU(JL) + ZCFRAC * ZUPF(JL,JK  ,JK1)
 334  CONTINUE
 335  CONTINUE
C
      DO 336 JL = 1, KDLON
      PFLUX(JL,1,JK1) = ZFU (JL)
 336  CONTINUE
C
 337  CONTINUE
C
C
      END IF
C
C
C*         2.3     END OF CLOUD EFFECT COMPUTATIONS
C
 230  CONTINUE
C
      IF (.NOT.LEVOIGT) THEN
        DO 231 JL = 1, KDLON
        ZFN10(JL) = PFLUX(JL,1,KLIM) + PFLUX(JL,2,KLIM)
 231    CONTINUE
        DO 233 JK = KLIM+1 , KFLEV+1
        DO 232 JL = 1, KDLON
        ZFN10(JL) = ZFN10(JL) + PCTS(JL,JK-1)
        PFLUX(JL,1,JK) = ZFN10(JL)
        PFLUX(JL,2,JK) = 0.0
 232    CONTINUE
 233    CONTINUE
      ENDIF
C
      RETURN
      END
      SUBROUTINE LWB(PDT0,PTAVE,PTL
     S  , PB,PBINT,PBSUIN,PBSUR,PBTOP,PDBSL
     S  , PGA,PGB,PGASUR,PGBSUR,PGATOP,PGBTOP)
      IMPLICIT none
      include "dimensions.h"
      include "dimphy.h"
      include "raddim.h"
      include "raddimlw.h"
C
C-----------------------------------------------------------------------
C     PURPOSE.
C     --------
C           COMPUTES PLANCK FUNCTIONS
C
C        EXPLICIT ARGUMENTS :
C        --------------------
C     ==== INPUTS ===
C PDT0   : (KDLON)             ; SURFACE TEMPERATURE DISCONTINUITY
C PTAVE  : (KDLON,KFLEV)       ; TEMPERATURE
C PTL    : (KDLON,0:KFLEV)     ; HALF LEVEL TEMPERATURE
C     ==== OUTPUTS ===
C PB     : (KDLON,Ninter,KFLEV+1); SPECTRAL HALF LEVEL PLANCK FUNCTION
C PBINT  : (KDLON,KFLEV+1)     ; HALF LEVEL PLANCK FUNCTION
C PBSUIN : (KDLON)             ; SURFACE PLANCK FUNCTION
C PBSUR  : (KDLON,Ninter)        ; SURFACE SPECTRAL PLANCK FUNCTION
C PBTOP  : (KDLON,Ninter)        ; TOP SPECTRAL PLANCK FUNCTION
C PDBSL  : (KDLON,Ninter,KFLEV*2); SUB-LAYER PLANCK FUNCTION GRADIENT
C PGA    : (KDLON,8,2,KFLEV); dB/dT-weighted LAYER PADE APPROXIMANTS
C PGB    : (KDLON,8,2,KFLEV); dB/dT-weighted LAYER PADE APPROXIMANTS
C PGASUR, PGBSUR (KDLON,8,2)   ; SURFACE PADE APPROXIMANTS
C PGATOP, PGBTOP (KDLON,8,2)   ; T.O.A. PADE APPROXIMANTS
C
C        IMPLICIT ARGUMENTS :   NONE
C        --------------------
C
C     METHOD.
C     -------
C
C          1. COMPUTES THE PLANCK FUNCTION ON ALL LEVELS AND HALF LEVELS
C     FROM A POLYNOMIAL DEVELOPMENT OF PLANCK FUNCTION
C
C     REFERENCE.
C     ----------
C
C        SEE RADIATION'S PART OF THE MODEL'S DOCUMENTATION AND
C        ECMWF RESEARCH DEPARTMENT DOCUMENTATION OF THE IFS           "
C
C     AUTHOR.
C     -------
C        JEAN-JACQUES MORCRETTE  *ECMWF*
C
C     MODIFICATIONS.
C     --------------
C        ORIGINAL : 89-07-14
C
C-----------------------------------------------------------------------
C
C ARGUMENTS:
C
      REAL*8 PDT0(KDLON)
      REAL*8 PTAVE(KDLON,KFLEV)
      REAL*8 PTL(KDLON,KFLEV+1)
C
      REAL*8 PB(KDLON,Ninter,KFLEV+1) ! SPECTRAL HALF LEVEL PLANCK FUNCTION
      REAL*8 PBINT(KDLON,KFLEV+1) ! HALF LEVEL PLANCK FUNCTION
      REAL*8 PBSUIN(KDLON) ! SURFACE PLANCK FUNCTION
      REAL*8 PBSUR(KDLON,Ninter) ! SURFACE SPECTRAL PLANCK FUNCTION
      REAL*8 PBTOP(KDLON,Ninter) ! TOP SPECTRAL PLANCK FUNCTION
      REAL*8 PDBSL(KDLON,Ninter,KFLEV*2) ! SUB-LAYER PLANCK FUNCTION GRADIENT
      REAL*8 PGA(KDLON,8,2,KFLEV) ! dB/dT-weighted LAYER PADE APPROXIMANTS
      REAL*8 PGB(KDLON,8,2,KFLEV) ! dB/dT-weighted LAYER PADE APPROXIMANTS
      REAL*8 PGASUR(KDLON,8,2) ! SURFACE PADE APPROXIMANTS
      REAL*8 PGBSUR(KDLON,8,2) ! SURFACE PADE APPROXIMANTS
      REAL*8 PGATOP(KDLON,8,2) ! T.O.A. PADE APPROXIMANTS
      REAL*8 PGBTOP(KDLON,8,2) ! T.O.A. PADE APPROXIMANTS
C
C-------------------------------------------------------------------------
C*  LOCAL VARIABLES:
      INTEGER INDB(KDLON),INDS(KDLON)
      REAL*8 ZBLAY(KDLON,KFLEV),ZBLEV(KDLON,KFLEV+1)
      REAL*8 ZRES(KDLON),ZRES2(KDLON),ZTI(KDLON),ZTI2(KDLON)
c
      INTEGER jk, jl, ic, jnu, jf, jg
      INTEGER jk1, jk2
      INTEGER k, j, ixtox, indto, ixtx, indt
      INTEGER indsu, indtp
      REAL*8 zdsto1, zdstox, zdst1, zdstx
c
C* Quelques parametres:
      REAL*8 TSTAND
      PARAMETER (TSTAND=250.0)
      REAL*8 TSTP
      PARAMETER (TSTP=12.5)
      INTEGER MXIXT
      PARAMETER (MXIXT=10)
C
C* Used Data Block:
      REAL*8 TINTP(11)
      SAVE TINTP
      REAL*8 GA(11,16,3), GB(11,16,3)
      SAVE GA, GB
      REAL*8 XP(6,6)
      SAVE XP
c
      DATA TINTP / 187.5, 200., 212.5, 225., 237.5, 250.,
     S             262.5, 275., 287.5, 300., 312.5 /
C-----------------------------------------------------------------------
C-- WATER VAPOR -- INT.1 -- 0- 500 CM-1 -- FROM ABS225 ----------------
C
C
C
C
C-- R.D. -- G = - 0.2 SLA
C
C
C----- INTERVAL = 1 ----- T =  187.5
C
C-- INDICES FOR PADE APPROXIMATION     1   15   29   45
      DATA (GA( 1, 1,IC),IC=1,3) /
     S 0.63499072E-02,-0.99506586E-03, 0.00000000E+00/
      DATA (GB( 1, 1,IC),IC=1,3) /
     S 0.63499072E-02, 0.97222852E-01, 0.10000000E+01/
      DATA (GA( 1, 2,IC),IC=1,3) /
     S 0.77266491E-02,-0.11661515E-02, 0.00000000E+00/
      DATA (GB( 1, 2,IC),IC=1,3) /
     S 0.77266491E-02, 0.10681591E+00, 0.10000000E+01/
C
C----- INTERVAL = 1 ----- T =  200.0
C
C-- INDICES FOR PADE APPROXIMATION     1   15   29   45
      DATA (GA( 2, 1,IC),IC=1,3) /
     S 0.65566348E-02,-0.10184169E-02, 0.00000000E+00/
      DATA (GB( 2, 1,IC),IC=1,3) /
     S 0.65566348E-02, 0.98862238E-01, 0.10000000E+01/
      DATA (GA( 2, 2,IC),IC=1,3) /
     S 0.81323287E-02,-0.11886130E-02, 0.00000000E+00/
      DATA (GB( 2, 2,IC),IC=1,3) /
     S 0.81323287E-02, 0.10921298E+00, 0.10000000E+01/
C
C----- INTERVAL = 1 ----- T =  212.5
C
C-- INDICES FOR PADE APPROXIMATION     1   15   29   45
      DATA (GA( 3, 1,IC),IC=1,3) /
     S 0.67849730E-02,-0.10404730E-02, 0.00000000E+00/
      DATA (GB( 3, 1,IC),IC=1,3) /
     S 0.67849730E-02, 0.10061504E+00, 0.10000000E+01/
      DATA (GA( 3, 2,IC),IC=1,3) /
     S 0.86507620E-02,-0.12139929E-02, 0.00000000E+00/
      DATA (GB( 3, 2,IC),IC=1,3) /
     S 0.86507620E-02, 0.11198225E+00, 0.10000000E+01/
C
C----- INTERVAL = 1 ----- T =  225.0
C
C-- INDICES FOR PADE APPROXIMATION     1   15   29   45
      DATA (GA( 4, 1,IC),IC=1,3) /
     S 0.70481947E-02,-0.10621792E-02, 0.00000000E+00/
      DATA (GB( 4, 1,IC),IC=1,3) /
     S 0.70481947E-02, 0.10256222E+00, 0.10000000E+01/
      DATA (GA( 4, 2,IC),IC=1,3) /
     S 0.92776391E-02,-0.12445811E-02, 0.00000000E+00/
      DATA (GB( 4, 2,IC),IC=1,3) /
     S 0.92776391E-02, 0.11487826E+00, 0.10000000E+01/
C
C----- INTERVAL = 1 ----- T =  237.5
C
C-- INDICES FOR PADE APPROXIMATION     1   15   29   45
      DATA (GA( 5, 1,IC),IC=1,3) /
     S 0.73585943E-02,-0.10847662E-02, 0.00000000E+00/
      DATA (GB( 5, 1,IC),IC=1,3) /
     S 0.73585943E-02, 0.10475952E+00, 0.10000000E+01/
      DATA (GA( 5, 2,IC),IC=1,3) /
     S 0.99806312E-02,-0.12807672E-02, 0.00000000E+00/
      DATA (GB( 5, 2,IC),IC=1,3) /
     S 0.99806312E-02, 0.11751113E+00, 0.10000000E+01/
C
C----- INTERVAL = 1 ----- T =  250.0
C
C-- INDICES FOR PADE APPROXIMATION     1   15   29   45
      DATA (GA( 6, 1,IC),IC=1,3) /
     S 0.77242818E-02,-0.11094726E-02, 0.00000000E+00/
      DATA (GB( 6, 1,IC),IC=1,3) /
     S 0.77242818E-02, 0.10720986E+00, 0.10000000E+01/
      DATA (GA( 6, 2,IC),IC=1,3) /
     S 0.10709803E-01,-0.13208251E-02, 0.00000000E+00/
      DATA (GB( 6, 2,IC),IC=1,3) /
     S 0.10709803E-01, 0.11951535E+00, 0.10000000E+01/
C
C----- INTERVAL = 1 ----- T =  262.5
C
C-- INDICES FOR PADE APPROXIMATION     1   15   29   45
      DATA (GA( 7, 1,IC),IC=1,3) /
     S 0.81472693E-02,-0.11372949E-02, 0.00000000E+00/
      DATA (GB( 7, 1,IC),IC=1,3) /
     S 0.81472693E-02, 0.10985370E+00, 0.10000000E+01/
      DATA (GA( 7, 2,IC),IC=1,3) /
     S 0.11414739E-01,-0.13619034E-02, 0.00000000E+00/
      DATA (GB( 7, 2,IC),IC=1,3) /
     S 0.11414739E-01, 0.12069945E+00, 0.10000000E+01/
C
C----- INTERVAL = 1 ----- T =  275.0
C
C-- INDICES FOR PADE APPROXIMATION     1   15   29   45
      DATA (GA( 8, 1,IC),IC=1,3) /
     S 0.86227527E-02,-0.11687683E-02, 0.00000000E+00/
      DATA (GB( 8, 1,IC),IC=1,3) /
     S 0.86227527E-02, 0.11257633E+00, 0.10000000E+01/
      DATA (GA( 8, 2,IC),IC=1,3) /
     S 0.12058772E-01,-0.14014165E-02, 0.00000000E+00/
      DATA (GB( 8, 2,IC),IC=1,3) /
     S 0.12058772E-01, 0.12108524E+00, 0.10000000E+01/
C
C----- INTERVAL = 1 ----- T =  287.5
C
C-- INDICES FOR PADE APPROXIMATION     1   15   29   45
      DATA (GA( 9, 1,IC),IC=1,3) /
     S 0.91396814E-02,-0.12038314E-02, 0.00000000E+00/
      DATA (GB( 9, 1,IC),IC=1,3) /
     S 0.91396814E-02, 0.11522980E+00, 0.10000000E+01/
      DATA (GA( 9, 2,IC),IC=1,3) /
     S 0.12623992E-01,-0.14378639E-02, 0.00000000E+00/
      DATA (GB( 9, 2,IC),IC=1,3) /
     S 0.12623992E-01, 0.12084229E+00, 0.10000000E+01/
C
C----- INTERVAL = 1 ----- T =  300.0
C
C-- INDICES FOR PADE APPROXIMATION     1   15   29   45
      DATA (GA(10, 1,IC),IC=1,3) /
     S 0.96825438E-02,-0.12418367E-02, 0.00000000E+00/
      DATA (GB(10, 1,IC),IC=1,3) /
     S 0.96825438E-02, 0.11766343E+00, 0.10000000E+01/
      DATA (GA(10, 2,IC),IC=1,3) /
     S 0.13108146E-01,-0.14708488E-02, 0.00000000E+00/
      DATA (GB(10, 2,IC),IC=1,3) /
     S 0.13108146E-01, 0.12019005E+00, 0.10000000E+01/
C
C----- INTERVAL = 1 ----- T =  312.5
C
C-- INDICES FOR PADE APPROXIMATION     1   15   29   45
      DATA (GA(11, 1,IC),IC=1,3) /
     S 0.10233955E-01,-0.12817135E-02, 0.00000000E+00/
      DATA (GB(11, 1,IC),IC=1,3) /
     S 0.10233955E-01, 0.11975320E+00, 0.10000000E+01/
      DATA (GA(11, 2,IC),IC=1,3) /
     S 0.13518390E-01,-0.15006791E-02, 0.00000000E+00/
      DATA (GB(11, 2,IC),IC=1,3) /
     S 0.13518390E-01, 0.11932684E+00, 0.10000000E+01/
C
C
C
C--- WATER VAPOR --- INTERVAL 2 -- 500-800 CM-1--- FROM ABS225 ---------
C
C
C
C
C--- R.D.  ---  G = 0.02 + 0.50 / ( 1 + 4.5 U )
C
C
C----- INTERVAL = 2 ----- T =  187.5
C
C-- INDICES FOR PADE APPROXIMATION     1   28   37   45
      DATA (GA( 1, 3,IC),IC=1,3) /
     S 0.11644593E+01, 0.41243390E+00, 0.00000000E+00/
      DATA (GB( 1, 3,IC),IC=1,3) /
     S 0.11644593E+01, 0.10346097E+01, 0.10000000E+01/
      DATA (GA( 1, 4,IC),IC=1,3) /
     S 0.12006968E+01, 0.48318936E+00, 0.00000000E+00/
      DATA (GB( 1, 4,IC),IC=1,3) /
     S 0.12006968E+01, 0.10626130E+01, 0.10000000E+01/
C
C----- INTERVAL = 2 ----- T =  200.0
C
C-- INDICES FOR PADE APPROXIMATION     1   28   37   45
      DATA (GA( 2, 3,IC),IC=1,3) /
     S 0.11747203E+01, 0.43407282E+00, 0.00000000E+00/
      DATA (GB( 2, 3,IC),IC=1,3) /
     S 0.11747203E+01, 0.10433655E+01, 0.10000000E+01/
      DATA (GA( 2, 4,IC),IC=1,3) /
     S 0.12108196E+01, 0.50501827E+00, 0.00000000E+00/
      DATA (GB( 2, 4,IC),IC=1,3) /
     S 0.12108196E+01, 0.10716026E+01, 0.10000000E+01/
C
C----- INTERVAL = 2 ----- T =  212.5
C
C-- INDICES FOR PADE APPROXIMATION     1   28   37   45
      DATA (GA( 3, 3,IC),IC=1,3) /
     S 0.11837872E+01, 0.45331413E+00, 0.00000000E+00/
      DATA (GB( 3, 3,IC),IC=1,3) /
     S 0.11837872E+01, 0.10511933E+01, 0.10000000E+01/
      DATA (GA( 3, 4,IC),IC=1,3) /
     S 0.12196717E+01, 0.52409502E+00, 0.00000000E+00/
      DATA (GB( 3, 4,IC),IC=1,3) /
     S 0.12196717E+01, 0.10795108E+01, 0.10000000E+01/
C
C----- INTERVAL = 2 ----- T =  225.0
C
C-- INDICES FOR PADE APPROXIMATION     1   28   37   45
      DATA (GA( 4, 3,IC),IC=1,3) /
     S 0.11918561E+01, 0.47048604E+00, 0.00000000E+00/
      DATA (GB( 4, 3,IC),IC=1,3) /
     S 0.11918561E+01, 0.10582150E+01, 0.10000000E+01/
      DATA (GA( 4, 4,IC),IC=1,3) /
     S 0.12274493E+01, 0.54085277E+00, 0.00000000E+00/
      DATA (GB( 4, 4,IC),IC=1,3) /
     S 0.12274493E+01, 0.10865006E+01, 0.10000000E+01/
C
C----- INTERVAL = 2 ----- T =  237.5
C
C-- INDICES FOR PADE APPROXIMATION     1   28   37   45
      DATA (GA( 5, 3,IC),IC=1,3) /
     S 0.11990757E+01, 0.48586286E+00, 0.00000000E+00/
      DATA (GB( 5, 3,IC),IC=1,3) /
     S 0.11990757E+01, 0.10645317E+01, 0.10000000E+01/
      DATA (GA( 5, 4,IC),IC=1,3) /
     S 0.12343189E+01, 0.55565422E+00, 0.00000000E+00/
      DATA (GB( 5, 4,IC),IC=1,3) /
     S 0.12343189E+01, 0.10927103E+01, 0.10000000E+01/
C
C----- INTERVAL = 2 ----- T =  250.0
C
C-- INDICES FOR PADE APPROXIMATION     1   28   37   45
      DATA (GA( 6, 3,IC),IC=1,3) /
     S 0.12055643E+01, 0.49968044E+00, 0.00000000E+00/
      DATA (GB( 6, 3,IC),IC=1,3) /
     S 0.12055643E+01, 0.10702313E+01, 0.10000000E+01/
      DATA (GA( 6, 4,IC),IC=1,3) /
     S 0.12404147E+01, 0.56878618E+00, 0.00000000E+00/
      DATA (GB( 6, 4,IC),IC=1,3) /
     S 0.12404147E+01, 0.10982489E+01, 0.10000000E+01/
C
C----- INTERVAL = 2 ----- T =  262.5
C
C-- INDICES FOR PADE APPROXIMATION     1   28   37   45
      DATA (GA( 7, 3,IC),IC=1,3) /
     S 0.12114186E+01, 0.51214132E+00, 0.00000000E+00/
      DATA (GB( 7, 3,IC),IC=1,3) /
     S 0.12114186E+01, 0.10753907E+01, 0.10000000E+01/
      DATA (GA( 7, 4,IC),IC=1,3) /
     S 0.12458431E+01, 0.58047395E+00, 0.00000000E+00/
      DATA (GB( 7, 4,IC),IC=1,3) /
     S 0.12458431E+01, 0.11032019E+01, 0.10000000E+01/
C
C----- INTERVAL = 2 ----- T =  275.0
C
C-- INDICES FOR PADE APPROXIMATION     1   28   37   45
      DATA (GA( 8, 3,IC),IC=1,3) /
     S 0.12167192E+01, 0.52341830E+00, 0.00000000E+00/
      DATA (GB( 8, 3,IC),IC=1,3) /
     S 0.12167192E+01, 0.10800762E+01, 0.10000000E+01/
      DATA (GA( 8, 4,IC),IC=1,3) /
     S 0.12506907E+01, 0.59089894E+00, 0.00000000E+00/
      DATA (GB( 8, 4,IC),IC=1,3) /
     S 0.12506907E+01, 0.11076379E+01, 0.10000000E+01/
C
C----- INTERVAL = 2 ----- T =  287.5
C
C-- INDICES FOR PADE APPROXIMATION     1   28   37   45
      DATA (GA( 9, 3,IC),IC=1,3) /
     S 0.12215344E+01, 0.53365803E+00, 0.00000000E+00/
      DATA (GB( 9, 3,IC),IC=1,3) /
     S 0.12215344E+01, 0.10843446E+01, 0.10000000E+01/
      DATA (GA( 9, 4,IC),IC=1,3) /
     S 0.12550299E+01, 0.60021475E+00, 0.00000000E+00/
      DATA (GB( 9, 4,IC),IC=1,3) /
     S 0.12550299E+01, 0.11116160E+01, 0.10000000E+01/
C
C----- INTERVAL = 2 ----- T =  300.0
C
C-- INDICES FOR PADE APPROXIMATION     1   28   37   45
      DATA (GA(10, 3,IC),IC=1,3) /
     S 0.12259226E+01, 0.54298448E+00, 0.00000000E+00/
      DATA (GB(10, 3,IC),IC=1,3) /
     S 0.12259226E+01, 0.10882439E+01, 0.10000000E+01/
      DATA (GA(10, 4,IC),IC=1,3) /
     S 0.12589256E+01, 0.60856112E+00, 0.00000000E+00/
      DATA (GB(10, 4,IC),IC=1,3) /
     S 0.12589256E+01, 0.11151910E+01, 0.10000000E+01/
C
C----- INTERVAL = 2 ----- T =  312.5
C
C-- INDICES FOR PADE APPROXIMATION     1   28   37   45
      DATA (GA(11, 3,IC),IC=1,3) /
     S 0.12299344E+01, 0.55150227E+00, 0.00000000E+00/
      DATA (GB(11, 3,IC),IC=1,3) /
     S 0.12299344E+01, 0.10918144E+01, 0.10000000E+01/
      DATA (GA(11, 4,IC),IC=1,3) /
     S 0.12624402E+01, 0.61607594E+00, 0.00000000E+00/
      DATA (GB(11, 4,IC),IC=1,3) /
     S 0.12624402E+01, 0.11184188E+01, 0.10000000E+01/
C
C
C
C
C
C
C- WATER VAPOR - INT. 3 -- 800-970 + 1110-1250 CM-1 -- FIT FROM 215 IS -
C
C
C-- WATER VAPOR LINES IN THE WINDOW REGION (800-1250 CM-1)
C
C
C
C--- G = 3.875E-03 ---------------
C
C----- INTERVAL = 3 ----- T =  187.5
C
C-- INDICES FOR PADE APPROXIMATION     1   28   37   45
      DATA (GA( 1, 7,IC),IC=1,3) /
     S 0.10192131E+02, 0.80737799E+01, 0.00000000E+00/
      DATA (GB( 1, 7,IC),IC=1,3) /
     S 0.10192131E+02, 0.82623280E+01, 0.10000000E+01/
      DATA (GA( 1, 8,IC),IC=1,3) /
     S 0.92439050E+01, 0.77425778E+01, 0.00000000E+00/
      DATA (GB( 1, 8,IC),IC=1,3) /
     S 0.92439050E+01, 0.79342219E+01, 0.10000000E+01/
C
C----- INTERVAL = 3 ----- T =  200.0
C
C-- INDICES FOR PADE APPROXIMATION     1   28   37   45
      DATA (GA( 2, 7,IC),IC=1,3) /
     S 0.97258602E+01, 0.79171158E+01, 0.00000000E+00/
      DATA (GB( 2, 7,IC),IC=1,3) /
     S 0.97258602E+01, 0.81072291E+01, 0.10000000E+01/
      DATA (GA( 2, 8,IC),IC=1,3) /
     S 0.87567422E+01, 0.75443460E+01, 0.00000000E+00/
      DATA (GB( 2, 8,IC),IC=1,3) /
     S 0.87567422E+01, 0.77373458E+01, 0.10000000E+01/
C
C----- INTERVAL = 3 ----- T =  212.5
C
C-- INDICES FOR PADE APPROXIMATION     1   28   37   45
      DATA (GA( 3, 7,IC),IC=1,3) /
     S 0.92992890E+01, 0.77609605E+01, 0.00000000E+00/
      DATA (GB( 3, 7,IC),IC=1,3) /
     S 0.92992890E+01, 0.79523834E+01, 0.10000000E+01/
      DATA (GA( 3, 8,IC),IC=1,3) /
     S 0.83270144E+01, 0.73526151E+01, 0.00000000E+00/
      DATA (GB( 3, 8,IC),IC=1,3) /
     S 0.83270144E+01, 0.75467334E+01, 0.10000000E+01/
C
C----- INTERVAL = 3 ----- T =  225.0
C
C-- INDICES FOR PADE APPROXIMATION     1   28   37   45
      DATA (GA( 4, 7,IC),IC=1,3) /
     S 0.89154021E+01, 0.76087371E+01, 0.00000000E+00/
      DATA (GB( 4, 7,IC),IC=1,3) /
     S 0.89154021E+01, 0.78012527E+01, 0.10000000E+01/
      DATA (GA( 4, 8,IC),IC=1,3) /
     S 0.79528337E+01, 0.71711188E+01, 0.00000000E+00/
      DATA (GB( 4, 8,IC),IC=1,3) /
     S 0.79528337E+01, 0.73661786E+01, 0.10000000E+01/
C
C----- INTERVAL = 3 ----- T =  237.5
C
C-- INDICES FOR PADE APPROXIMATION     1   28   37   45
      DATA (GA( 5, 7,IC),IC=1,3) /
     S 0.85730084E+01, 0.74627112E+01, 0.00000000E+00/
      DATA (GB( 5, 7,IC),IC=1,3) /
     S 0.85730084E+01, 0.76561458E+01, 0.10000000E+01/
      DATA (GA( 5, 8,IC),IC=1,3) /
     S 0.76286839E+01, 0.70015571E+01, 0.00000000E+00/
      DATA (GB( 5, 8,IC),IC=1,3) /
     S 0.76286839E+01, 0.71974319E+01, 0.10000000E+01/
C
C----- INTERVAL = 3 ----- T =  250.0
C
C-- INDICES FOR PADE APPROXIMATION     1   28   37   45
      DATA (GA( 6, 7,IC),IC=1,3) /
     S 0.82685838E+01, 0.73239981E+01, 0.00000000E+00/
      DATA (GB( 6, 7,IC),IC=1,3) /
     S 0.82685838E+01, 0.75182174E+01, 0.10000000E+01/
      DATA (GA( 6, 8,IC),IC=1,3) /
     S 0.73477879E+01, 0.68442532E+01, 0.00000000E+00/
      DATA (GB( 6, 8,IC),IC=1,3) /
     S 0.73477879E+01, 0.70408543E+01, 0.10000000E+01/
C
C----- INTERVAL = 3 ----- T =  262.5
C
C-- INDICES FOR PADE APPROXIMATION     1   28   37   45
      DATA (GA( 7, 7,IC),IC=1,3) /
     S 0.79978921E+01, 0.71929934E+01, 0.00000000E+00/
      DATA (GB( 7, 7,IC),IC=1,3) /
     S 0.79978921E+01, 0.73878952E+01, 0.10000000E+01/
      DATA (GA( 7, 8,IC),IC=1,3) /
     S 0.71035818E+01, 0.66987996E+01, 0.00000000E+00/
      DATA (GB( 7, 8,IC),IC=1,3) /
     S 0.71035818E+01, 0.68960649E+01, 0.10000000E+01/
C
C----- INTERVAL = 3 ----- T =  275.0
C
C-- INDICES FOR PADE APPROXIMATION     1   28   37   45
      DATA (GA( 8, 7,IC),IC=1,3) /
     S 0.77568055E+01, 0.70697065E+01, 0.00000000E+00/
      DATA (GB( 8, 7,IC),IC=1,3) /
     S 0.77568055E+01, 0.72652133E+01, 0.10000000E+01/
      DATA (GA( 8, 8,IC),IC=1,3) /
     S 0.68903312E+01, 0.65644820E+01, 0.00000000E+00/
      DATA (GB( 8, 8,IC),IC=1,3) /
     S 0.68903312E+01, 0.67623672E+01, 0.10000000E+01/
C
C----- INTERVAL = 3 ----- T =  287.5
C
C-- INDICES FOR PADE APPROXIMATION     1   28   37   45
      DATA (GA( 9, 7,IC),IC=1,3) /
     S 0.75416266E+01, 0.69539626E+01, 0.00000000E+00/
      DATA (GB( 9, 7,IC),IC=1,3) /
     S 0.75416266E+01, 0.71500151E+01, 0.10000000E+01/
      DATA (GA( 9, 8,IC),IC=1,3) /
     S 0.67032875E+01, 0.64405267E+01, 0.00000000E+00/
      DATA (GB( 9, 8,IC),IC=1,3) /
     S 0.67032875E+01, 0.66389989E+01, 0.10000000E+01/
C
C----- INTERVAL = 3 ----- T =  300.0
C
C-- INDICES FOR PADE APPROXIMATION     1   28   37   45
      DATA (GA(10, 7,IC),IC=1,3) /
     S 0.73491694E+01, 0.68455144E+01, 0.00000000E+00/
      DATA (GB(10, 7,IC),IC=1,3) /
     S 0.73491694E+01, 0.70420667E+01, 0.10000000E+01/
      DATA (GA(10, 8,IC),IC=1,3) /
     S 0.65386461E+01, 0.63262376E+01, 0.00000000E+00/
      DATA (GB(10, 8,IC),IC=1,3) /
     S 0.65386461E+01, 0.65252707E+01, 0.10000000E+01/
C
C----- INTERVAL = 3 ----- T =  312.5
C
C-- INDICES FOR PADE APPROXIMATION     1   28   37   45
      DATA (GA(11, 7,IC),IC=1,3) /
     S 0.71767400E+01, 0.67441020E+01, 0.00000000E+00/
      DATA (GB(11, 7,IC),IC=1,3) /
     S 0.71767400E+01, 0.69411177E+01, 0.10000000E+01/
      DATA (GA(11, 8,IC),IC=1,3) /
     S 0.63934377E+01, 0.62210701E+01, 0.00000000E+00/
      DATA (GB(11, 8,IC),IC=1,3) /
     S 0.63934377E+01, 0.64206412E+01, 0.10000000E+01/
C
C
C-- WATER VAPOR -- 970-1110 CM-1 ----------------------------------------
C
C-- G = 3.6E-03
C
C----- INTERVAL = 4 ----- T =  187.5
C
C-- INDICES FOR PADE APPROXIMATION     1   28   37   45
      DATA (GA( 1, 9,IC),IC=1,3) /
     S 0.24870635E+02, 0.10542131E+02, 0.00000000E+00/
      DATA (GB( 1, 9,IC),IC=1,3) /
     S 0.24870635E+02, 0.10656640E+02, 0.10000000E+01/
      DATA (GA( 1,10,IC),IC=1,3) /
     S 0.24586283E+02, 0.10490353E+02, 0.00000000E+00/
      DATA (GB( 1,10,IC),IC=1,3) /
     S 0.24586283E+02, 0.10605856E+02, 0.10000000E+01/
C
C----- INTERVAL = 4 ----- T =  200.0
C
C-- INDICES FOR PADE APPROXIMATION     1   28   37   45
      DATA (GA( 2, 9,IC),IC=1,3) /
     S 0.24725591E+02, 0.10515895E+02, 0.00000000E+00/
      DATA (GB( 2, 9,IC),IC=1,3) /
     S 0.24725591E+02, 0.10630910E+02, 0.10000000E+01/
      DATA (GA( 2,10,IC),IC=1,3) /
     S 0.24441465E+02, 0.10463512E+02, 0.00000000E+00/
      DATA (GB( 2,10,IC),IC=1,3) /
     S 0.24441465E+02, 0.10579514E+02, 0.10000000E+01/
C
C----- INTERVAL = 4 ----- T =  212.5
C
C-- INDICES FOR PADE APPROXIMATION     1   28   37   45
      DATA (GA( 3, 9,IC),IC=1,3) /
     S 0.24600320E+02, 0.10492949E+02, 0.00000000E+00/
      DATA (GB( 3, 9,IC),IC=1,3) /
     S 0.24600320E+02, 0.10608399E+02, 0.10000000E+01/
      DATA (GA( 3,10,IC),IC=1,3) /
     S 0.24311657E+02, 0.10439183E+02, 0.00000000E+00/
      DATA (GB( 3,10,IC),IC=1,3) /
     S 0.24311657E+02, 0.10555632E+02, 0.10000000E+01/
C
C----- INTERVAL = 4 ----- T =  225.0
C
C-- INDICES FOR PADE APPROXIMATION     1   28   37   45
      DATA (GA( 4, 9,IC),IC=1,3) /
     S 0.24487300E+02, 0.10472049E+02, 0.00000000E+00/
      DATA (GB( 4, 9,IC),IC=1,3) /
     S 0.24487300E+02, 0.10587891E+02, 0.10000000E+01/
      DATA (GA( 4,10,IC),IC=1,3) /
     S 0.24196167E+02, 0.10417324E+02, 0.00000000E+00/
      DATA (GB( 4,10,IC),IC=1,3) /
     S 0.24196167E+02, 0.10534169E+02, 0.10000000E+01/
C
C----- INTERVAL = 4 ----- T =  237.5
C
C-- INDICES FOR PADE APPROXIMATION     1   28   37   45
      DATA (GA( 5, 9,IC),IC=1,3) /
     S 0.24384935E+02, 0.10452961E+02, 0.00000000E+00/
      DATA (GB( 5, 9,IC),IC=1,3) /
     S 0.24384935E+02, 0.10569156E+02, 0.10000000E+01/
      DATA (GA( 5,10,IC),IC=1,3) /
     S 0.24093406E+02, 0.10397704E+02, 0.00000000E+00/
      DATA (GB( 5,10,IC),IC=1,3) /
     S 0.24093406E+02, 0.10514900E+02, 0.10000000E+01/
C
C----- INTERVAL = 4 ----- T =  250.0
C
C-- INDICES FOR PADE APPROXIMATION     1   28   37   45
      DATA (GA( 6, 9,IC),IC=1,3) /
     S 0.24292341E+02, 0.10435562E+02, 0.00000000E+00/
      DATA (GB( 6, 9,IC),IC=1,3) /
     S 0.24292341E+02, 0.10552075E+02, 0.10000000E+01/
      DATA (GA( 6,10,IC),IC=1,3) /
     S 0.24001597E+02, 0.10380038E+02, 0.00000000E+00/
      DATA (GB( 6,10,IC),IC=1,3) /
     S 0.24001597E+02, 0.10497547E+02, 0.10000000E+01/
C
C----- INTERVAL = 4 ----- T =  262.5
C
C-- INDICES FOR PADE APPROXIMATION     1   28   37   45
      DATA (GA( 7, 9,IC),IC=1,3) /
     S 0.24208572E+02, 0.10419710E+02, 0.00000000E+00/
      DATA (GB( 7, 9,IC),IC=1,3) /
     S 0.24208572E+02, 0.10536510E+02, 0.10000000E+01/
      DATA (GA( 7,10,IC),IC=1,3) /
     S 0.23919098E+02, 0.10364052E+02, 0.00000000E+00/
      DATA (GB( 7,10,IC),IC=1,3) /
     S 0.23919098E+02, 0.10481842E+02, 0.10000000E+01/
C
C----- INTERVAL = 4 ----- T =  275.0
C
C-- INDICES FOR PADE APPROXIMATION     1   28   37   45
      DATA (GA( 8, 9,IC),IC=1,3) /
     S 0.24132642E+02, 0.10405247E+02, 0.00000000E+00/
      DATA (GB( 8, 9,IC),IC=1,3) /
     S 0.24132642E+02, 0.10522307E+02, 0.10000000E+01/
      DATA (GA( 8,10,IC),IC=1,3) /
     S 0.23844511E+02, 0.10349509E+02, 0.00000000E+00/
      DATA (GB( 8,10,IC),IC=1,3) /
     S 0.23844511E+02, 0.10467553E+02, 0.10000000E+01/
C
C----- INTERVAL = 4 ----- T =  287.5
C
C-- INDICES FOR PADE APPROXIMATION     1   28   37   45
      DATA (GA( 9, 9,IC),IC=1,3) /
     S 0.24063614E+02, 0.10392022E+02, 0.00000000E+00/
      DATA (GB( 9, 9,IC),IC=1,3) /
     S 0.24063614E+02, 0.10509317E+02, 0.10000000E+01/
      DATA (GA( 9,10,IC),IC=1,3) /
     S 0.23776708E+02, 0.10336215E+02, 0.00000000E+00/
      DATA (GB( 9,10,IC),IC=1,3) /
     S 0.23776708E+02, 0.10454488E+02, 0.10000000E+01/
C
C----- INTERVAL = 4 ----- T =  300.0
C
C-- INDICES FOR PADE APPROXIMATION     1   28   37   45
      DATA (GA(10, 9,IC),IC=1,3) /
     S 0.24000649E+02, 0.10379892E+02, 0.00000000E+00/
      DATA (GB(10, 9,IC),IC=1,3) /
     S 0.24000649E+02, 0.10497402E+02, 0.10000000E+01/
      DATA (GA(10,10,IC),IC=1,3) /
     S 0.23714816E+02, 0.10324018E+02, 0.00000000E+00/
      DATA (GB(10,10,IC),IC=1,3) /
     S 0.23714816E+02, 0.10442501E+02, 0.10000000E+01/
C
C----- INTERVAL = 4 ----- T =  312.5
C
C-- INDICES FOR PADE APPROXIMATION     1   28   37   45
      DATA (GA(11, 9,IC),IC=1,3) /
     S 0.23943021E+02, 0.10368736E+02, 0.00000000E+00/
      DATA (GB(11, 9,IC),IC=1,3) /
     S 0.23943021E+02, 0.10486443E+02, 0.10000000E+01/
      DATA (GA(11,10,IC),IC=1,3) /
     S 0.23658197E+02, 0.10312808E+02, 0.00000000E+00/
      DATA (GB(11,10,IC),IC=1,3) /
     S 0.23658197E+02, 0.10431483E+02, 0.10000000E+01/
C
C
C
C-- H2O -- WEAKER PARTS OF THE STRONG BANDS  -- FROM ABS225 ----
C
C-- WATER VAPOR --- 350 - 500 CM-1
C
C-- G = - 0.2*SLA, 0.0 +0.5/(1+0.5U)
C
C----- INTERVAL = 5 ----- T =  187.5
C
C-- INDICES FOR PADE APPROXIMATION   1 35 40 45
      DATA (GA( 1, 5,IC),IC=1,3) /
     S 0.15750172E+00,-0.22159303E-01, 0.00000000E+00/
      DATA (GB( 1, 5,IC),IC=1,3) /
     S 0.15750172E+00, 0.38103212E+00, 0.10000000E+01/
      DATA (GA( 1, 6,IC),IC=1,3) /
     S 0.17770551E+00,-0.24972399E-01, 0.00000000E+00/
      DATA (GB( 1, 6,IC),IC=1,3) /
     S 0.17770551E+00, 0.41646579E+00, 0.10000000E+01/
C
C----- INTERVAL = 5 ----- T =  200.0
C
C-- INDICES FOR PADE APPROXIMATION   1 35 40 45
      DATA (GA( 2, 5,IC),IC=1,3) /
     S 0.16174076E+00,-0.22748917E-01, 0.00000000E+00/
      DATA (GB( 2, 5,IC),IC=1,3) /
     S 0.16174076E+00, 0.38913800E+00, 0.10000000E+01/
      DATA (GA( 2, 6,IC),IC=1,3) /
     S 0.18176757E+00,-0.25537247E-01, 0.00000000E+00/
      DATA (GB( 2, 6,IC),IC=1,3) /
     S 0.18176757E+00, 0.42345095E+00, 0.10000000E+01/
C
C----- INTERVAL = 5 ----- T =  212.5
C
C-- INDICES FOR PADE APPROXIMATION   1 35 40 45
      DATA (GA( 3, 5,IC),IC=1,3) /
     S 0.16548628E+00,-0.23269898E-01, 0.00000000E+00/
      DATA (GB( 3, 5,IC),IC=1,3) /
     S 0.16548628E+00, 0.39613651E+00, 0.10000000E+01/
      DATA (GA( 3, 6,IC),IC=1,3) /
     S 0.18527967E+00,-0.26025624E-01, 0.00000000E+00/
      DATA (GB( 3, 6,IC),IC=1,3) /
     S 0.18527967E+00, 0.42937476E+00, 0.10000000E+01/
C
C----- INTERVAL = 5 ----- T =  225.0
C
C-- INDICES FOR PADE APPROXIMATION   1 35 40 45
      DATA (GA( 4, 5,IC),IC=1,3) /
     S 0.16881124E+00,-0.23732392E-01, 0.00000000E+00/
      DATA (GB( 4, 5,IC),IC=1,3) /
     S 0.16881124E+00, 0.40222421E+00, 0.10000000E+01/
      DATA (GA( 4, 6,IC),IC=1,3) /
     S 0.18833348E+00,-0.26450280E-01, 0.00000000E+00/
      DATA (GB( 4, 6,IC),IC=1,3) /
     S 0.18833348E+00, 0.43444062E+00, 0.10000000E+01/
C
C----- INTERVAL = 5 ----- T =  237.5
C
C-- INDICES FOR PADE APPROXIMATION   1 35 40 45
      DATA (GA( 5, 5,IC),IC=1,3) /
     S 0.17177839E+00,-0.24145123E-01, 0.00000000E+00/
      DATA (GB( 5, 5,IC),IC=1,3) /
     S 0.17177839E+00, 0.40756010E+00, 0.10000000E+01/
      DATA (GA( 5, 6,IC),IC=1,3) /
     S 0.19100108E+00,-0.26821236E-01, 0.00000000E+00/
      DATA (GB( 5, 6,IC),IC=1,3) /
     S 0.19100108E+00, 0.43880316E+00, 0.10000000E+01/
C
C----- INTERVAL = 5 ----- T =  250.0
C
C-- INDICES FOR PADE APPROXIMATION   1 35 40 45
      DATA (GA( 6, 5,IC),IC=1,3) /
     S 0.17443933E+00,-0.24515269E-01, 0.00000000E+00/
      DATA (GB( 6, 5,IC),IC=1,3) /
     S 0.17443933E+00, 0.41226954E+00, 0.10000000E+01/
      DATA (GA( 6, 6,IC),IC=1,3) /
     S 0.19334122E+00,-0.27146657E-01, 0.00000000E+00/
      DATA (GB( 6, 6,IC),IC=1,3) /
     S 0.19334122E+00, 0.44258354E+00, 0.10000000E+01/
C
C----- INTERVAL = 5 ----- T =  262.5
C
C-- INDICES FOR PADE APPROXIMATION   1 35 40 45
      DATA (GA( 7, 5,IC),IC=1,3) /
     S 0.17683622E+00,-0.24848690E-01, 0.00000000E+00/
      DATA (GB( 7, 5,IC),IC=1,3) /
     S 0.17683622E+00, 0.41645142E+00, 0.10000000E+01/
      DATA (GA( 7, 6,IC),IC=1,3) /
     S 0.19540288E+00,-0.27433354E-01, 0.00000000E+00/
      DATA (GB( 7, 6,IC),IC=1,3) /
     S 0.19540288E+00, 0.44587882E+00, 0.10000000E+01/
C
C----- INTERVAL = 5 ----- T =  275.0
C
C-- INDICES FOR PADE APPROXIMATION   1 35 40 45
      DATA (GA( 8, 5,IC),IC=1,3) /
     S 0.17900375E+00,-0.25150210E-01, 0.00000000E+00/
      DATA (GB( 8, 5,IC),IC=1,3) /
     S 0.17900375E+00, 0.42018474E+00, 0.10000000E+01/
      DATA (GA( 8, 6,IC),IC=1,3) /
     S 0.19722732E+00,-0.27687065E-01, 0.00000000E+00/
      DATA (GB( 8, 6,IC),IC=1,3) /
     S 0.19722732E+00, 0.44876776E+00, 0.10000000E+01/
C
C----- INTERVAL = 5 ----- T =  287.5
C
C-- INDICES FOR PADE APPROXIMATION   1 35 40 45
      DATA (GA( 9, 5,IC),IC=1,3) /
     S 0.18097099E+00,-0.25423873E-01, 0.00000000E+00/
      DATA (GB( 9, 5,IC),IC=1,3) /
     S 0.18097099E+00, 0.42353379E+00, 0.10000000E+01/
      DATA (GA( 9, 6,IC),IC=1,3) /
     S 0.19884918E+00,-0.27912608E-01, 0.00000000E+00/
      DATA (GB( 9, 6,IC),IC=1,3) /
     S 0.19884918E+00, 0.45131451E+00, 0.10000000E+01/
C
C----- INTERVAL = 5 ----- T =  300.0
C
C-- INDICES FOR PADE APPROXIMATION   1 35 40 45
      DATA (GA(10, 5,IC),IC=1,3) /
     S 0.18276283E+00,-0.25673139E-01, 0.00000000E+00/
      DATA (GB(10, 5,IC),IC=1,3) /
     S 0.18276283E+00, 0.42655211E+00, 0.10000000E+01/
      DATA (GA(10, 6,IC),IC=1,3) /
     S 0.20029696E+00,-0.28113944E-01, 0.00000000E+00/
      DATA (GB(10, 6,IC),IC=1,3) /
     S 0.20029696E+00, 0.45357095E+00, 0.10000000E+01/
C
C----- INTERVAL = 5 ----- T =  312.5
C
C-- INDICES FOR PADE APPROXIMATION   1 35 40 45
      DATA (GA(11, 5,IC),IC=1,3) /
     S 0.18440117E+00,-0.25901055E-01, 0.00000000E+00/
      DATA (GB(11, 5,IC),IC=1,3) /
     S 0.18440117E+00, 0.42928533E+00, 0.10000000E+01/
      DATA (GA(11, 6,IC),IC=1,3) /
     S 0.20159300E+00,-0.28294180E-01, 0.00000000E+00/
      DATA (GB(11, 6,IC),IC=1,3) /
     S 0.20159300E+00, 0.45557797E+00, 0.10000000E+01/
C
C
C
C
C- WATER VAPOR - WINGS OF VIBRATION-ROTATION BAND - 1250-1450+1880-2820 -
C--- G = 0.0
C
C
C----- INTERVAL = 6 ----- T =  187.5
C
C-- INDICES FOR PADE APPROXIMATION   1 35 40 45
      DATA (GA( 1,11,IC),IC=1,3) /
     S 0.11990218E+02,-0.12823142E+01, 0.00000000E+00/
      DATA (GB( 1,11,IC),IC=1,3) /
     S 0.11990218E+02, 0.26681588E+02, 0.10000000E+01/
      DATA (GA( 1,12,IC),IC=1,3) /
     S 0.79709806E+01,-0.74805226E+00, 0.00000000E+00/
      DATA (GB( 1,12,IC),IC=1,3) /
     S 0.79709806E+01, 0.18377807E+02, 0.10000000E+01/
C
C----- INTERVAL = 6 ----- T =  200.0
C
C-- INDICES FOR PADE APPROXIMATION   1 35 40 45
      DATA (GA( 2,11,IC),IC=1,3) /
     S 0.10904073E+02,-0.10571588E+01, 0.00000000E+00/
      DATA (GB( 2,11,IC),IC=1,3) /
     S 0.10904073E+02, 0.24728346E+02, 0.10000000E+01/
      DATA (GA( 2,12,IC),IC=1,3) /
     S 0.75400737E+01,-0.56252739E+00, 0.00000000E+00/
      DATA (GB( 2,12,IC),IC=1,3) /
     S 0.75400737E+01, 0.17643148E+02, 0.10000000E+01/
C
C----- INTERVAL = 6 ----- T =  212.5
C
C-- INDICES FOR PADE APPROXIMATION   1 35 40 45
      DATA (GA( 3,11,IC),IC=1,3) /
     S 0.89126838E+01,-0.74864953E+00, 0.00000000E+00/
      DATA (GB( 3,11,IC),IC=1,3) /
     S 0.89126838E+01, 0.20551342E+02, 0.10000000E+01/
      DATA (GA( 3,12,IC),IC=1,3) /
     S 0.81804377E+01,-0.46188072E+00, 0.00000000E+00/
      DATA (GB( 3,12,IC),IC=1,3) /
     S 0.81804377E+01, 0.19296161E+02, 0.10000000E+01/
C
C----- INTERVAL = 6 ----- T =  225.0
C
C-- INDICES FOR PADE APPROXIMATION   1 35 40 45
      DATA (GA( 4,11,IC),IC=1,3) /
     S 0.85622405E+01,-0.58705980E+00, 0.00000000E+00/
      DATA (GB( 4,11,IC),IC=1,3) /
     S 0.85622405E+01, 0.19955244E+02, 0.10000000E+01/
      DATA (GA( 4,12,IC),IC=1,3) /
     S 0.10564339E+02,-0.40712065E+00, 0.00000000E+00/
      DATA (GB( 4,12,IC),IC=1,3) /
     S 0.10564339E+02, 0.24951120E+02, 0.10000000E+01/
C
C----- INTERVAL = 6 ----- T =  237.5
C
C-- INDICES FOR PADE APPROXIMATION   1 35 40 45
      DATA (GA( 5,11,IC),IC=1,3) /
     S 0.94892164E+01,-0.49305772E+00, 0.00000000E+00/
      DATA (GB( 5,11,IC),IC=1,3) /
     S 0.94892164E+01, 0.22227100E+02, 0.10000000E+01/
      DATA (GA( 5,12,IC),IC=1,3) /
     S 0.46896789E+02,-0.15295996E+01, 0.00000000E+00/
      DATA (GB( 5,12,IC),IC=1,3) /
     S 0.46896789E+02, 0.10957372E+03, 0.10000000E+01/
C
C----- INTERVAL = 6 ----- T =  250.0
C
C-- INDICES FOR PADE APPROXIMATION   1 35 40 45
      DATA (GA( 6,11,IC),IC=1,3) /
     S 0.13580937E+02,-0.51461431E+00, 0.00000000E+00/
      DATA (GB( 6,11,IC),IC=1,3) /
     S 0.13580937E+02, 0.31770288E+02, 0.10000000E+01/
      DATA (GA( 6,12,IC),IC=1,3) /
     S-0.30926524E+01, 0.43555255E+00, 0.00000000E+00/
      DATA (GB( 6,12,IC),IC=1,3) /
     S-0.30926524E+01,-0.67432659E+01, 0.10000000E+01/
C
C----- INTERVAL = 6 ----- T =  262.5
C
C-- INDICES FOR PADE APPROXIMATION   1 35 40 45
      DATA (GA( 7,11,IC),IC=1,3) /
     S-0.32050918E+03, 0.12373350E+02, 0.00000000E+00/
      DATA (GB( 7,11,IC),IC=1,3) /
     S-0.32050918E+03,-0.74061287E+03, 0.10000000E+01/
      DATA (GA( 7,12,IC),IC=1,3) /
     S 0.85742941E+00, 0.50380874E+00, 0.00000000E+00/
      DATA (GB( 7,12,IC),IC=1,3) /
     S 0.85742941E+00, 0.24550746E+01, 0.10000000E+01/
C
C----- INTERVAL = 6 ----- T =  275.0
C
C-- INDICES FOR PADE APPROXIMATION   1 35 40 45
      DATA (GA( 8,11,IC),IC=1,3) /
     S-0.37133165E+01, 0.44809588E+00, 0.00000000E+00/
      DATA (GB( 8,11,IC),IC=1,3) /
     S-0.37133165E+01,-0.81329826E+01, 0.10000000E+01/
      DATA (GA( 8,12,IC),IC=1,3) /
     S 0.19164038E+01, 0.68537352E+00, 0.00000000E+00/
      DATA (GB( 8,12,IC),IC=1,3) /
     S 0.19164038E+01, 0.49089917E+01, 0.10000000E+01/
C
C----- INTERVAL = 6 ----- T =  287.5
C
C-- INDICES FOR PADE APPROXIMATION   1 35 40 45
      DATA (GA( 9,11,IC),IC=1,3) /
     S 0.18890836E+00, 0.46548918E+00, 0.00000000E+00/
      DATA (GB( 9,11,IC),IC=1,3) /
     S 0.18890836E+00, 0.90279822E+00, 0.10000000E+01/
      DATA (GA( 9,12,IC),IC=1,3) /
     S 0.23513199E+01, 0.89437630E+00, 0.00000000E+00/
      DATA (GB( 9,12,IC),IC=1,3) /
     S 0.23513199E+01, 0.59008712E+01, 0.10000000E+01/
C
C----- INTERVAL = 6 ----- T =  300.0
C
C-- INDICES FOR PADE APPROXIMATION   1 35 40 45
      DATA (GA(10,11,IC),IC=1,3) /
     S 0.14209226E+01, 0.59121475E+00, 0.00000000E+00/
      DATA (GB(10,11,IC),IC=1,3) /
     S 0.14209226E+01, 0.37532746E+01, 0.10000000E+01/
      DATA (GA(10,12,IC),IC=1,3) /
     S 0.25566644E+01, 0.11127003E+01, 0.00000000E+00/
      DATA (GB(10,12,IC),IC=1,3) /
     S 0.25566644E+01, 0.63532616E+01, 0.10000000E+01/
C
C----- INTERVAL = 6 ----- T =  312.5
C
C-- INDICES FOR PADE APPROXIMATION   1 35 40 45
      DATA (GA(11,11,IC),IC=1,3) /
     S 0.19817679E+01, 0.74676119E+00, 0.00000000E+00/
      DATA (GB(11,11,IC),IC=1,3) /
     S 0.19817679E+01, 0.50437916E+01, 0.10000000E+01/
      DATA (GA(11,12,IC),IC=1,3) /
     S 0.26555181E+01, 0.13329782E+01, 0.00000000E+00/
      DATA (GB(11,12,IC),IC=1,3) /
     S 0.26555181E+01, 0.65558627E+01, 0.10000000E+01/
C
C
C
C
C
C-- END WATER VAPOR
C
C
C-- CO2 -- INT.2 -- 500-800 CM-1 --- FROM ABS225 ----------------------
C
C
C
C-- FIU = 0.8 + MAX(0.35,(7-IU)*0.9)  , X/T,  9
C
C----- INTERVAL = 2 ----- T =  187.5
C
C-- INDICES FOR PADE APPROXIMATION   1 30 38 45
      DATA (GA( 1,13,IC),IC=1,3) /
     S 0.87668459E-01, 0.13845511E+01, 0.00000000E+00/
      DATA (GB( 1,13,IC),IC=1,3) /
     S 0.87668459E-01, 0.23203798E+01, 0.10000000E+01/
      DATA (GA( 1,14,IC),IC=1,3) /
     S 0.74878820E-01, 0.11718758E+01, 0.00000000E+00/
      DATA (GB( 1,14,IC),IC=1,3) /
     S 0.74878820E-01, 0.20206726E+01, 0.10000000E+01/
C
C----- INTERVAL = 2 ----- T =  200.0
C
C-- INDICES FOR PADE APPROXIMATION   1 30 38 45
      DATA (GA( 2,13,IC),IC=1,3) /
     S 0.83754276E-01, 0.13187042E+01, 0.00000000E+00/
      DATA (GB( 2,13,IC),IC=1,3) /
     S 0.83754276E-01, 0.22288925E+01, 0.10000000E+01/
      DATA (GA( 2,14,IC),IC=1,3) /
     S 0.71650966E-01, 0.11216131E+01, 0.00000000E+00/
      DATA (GB( 2,14,IC),IC=1,3) /
     S 0.71650966E-01, 0.19441824E+01, 0.10000000E+01/
C
C----- INTERVAL = 2 ----- T =  212.5
C
C-- INDICES FOR PADE APPROXIMATION   1 30 38 45
      DATA (GA( 3,13,IC),IC=1,3) /
     S 0.80460283E-01, 0.12644396E+01, 0.00000000E+00/
      DATA (GB( 3,13,IC),IC=1,3) /
     S 0.80460283E-01, 0.21515593E+01, 0.10000000E+01/
      DATA (GA( 3,14,IC),IC=1,3) /
     S 0.68979615E-01, 0.10809473E+01, 0.00000000E+00/
      DATA (GB( 3,14,IC),IC=1,3) /
     S 0.68979615E-01, 0.18807257E+01, 0.10000000E+01/
C
C----- INTERVAL = 2 ----- T =  225.0
C
C-- INDICES FOR PADE APPROXIMATION   1 30 38 45
      DATA (GA( 4,13,IC),IC=1,3) /
     S 0.77659686E-01, 0.12191543E+01, 0.00000000E+00/
      DATA (GB( 4,13,IC),IC=1,3) /
     S 0.77659686E-01, 0.20855896E+01, 0.10000000E+01/
      DATA (GA( 4,14,IC),IC=1,3) /
     S 0.66745345E-01, 0.10476396E+01, 0.00000000E+00/
      DATA (GB( 4,14,IC),IC=1,3) /
     S 0.66745345E-01, 0.18275618E+01, 0.10000000E+01/
C
C----- INTERVAL = 2 ----- T =  237.5
C
C-- INDICES FOR PADE APPROXIMATION   1 30 38 45
      DATA (GA( 5,13,IC),IC=1,3) /
     S 0.75257056E-01, 0.11809511E+01, 0.00000000E+00/
      DATA (GB( 5,13,IC),IC=1,3) /
     S 0.75257056E-01, 0.20288489E+01, 0.10000000E+01/
      DATA (GA( 5,14,IC),IC=1,3) /
     S 0.64857571E-01, 0.10200373E+01, 0.00000000E+00/
      DATA (GB( 5,14,IC),IC=1,3) /
     S 0.64857571E-01, 0.17825910E+01, 0.10000000E+01/
C
C----- INTERVAL = 2 ----- T =  250.0
C
C-- INDICES FOR PADE APPROXIMATION   1 30 38 45
      DATA (GA( 6,13,IC),IC=1,3) /
     S 0.73179175E-01, 0.11484154E+01, 0.00000000E+00/
      DATA (GB( 6,13,IC),IC=1,3) /
     S 0.73179175E-01, 0.19796791E+01, 0.10000000E+01/
      DATA (GA( 6,14,IC),IC=1,3) /
     S 0.63248495E-01, 0.99692726E+00, 0.00000000E+00/
      DATA (GB( 6,14,IC),IC=1,3) /
     S 0.63248495E-01, 0.17442308E+01, 0.10000000E+01/
C
C----- INTERVAL = 2 ----- T =  262.5
C
C-- INDICES FOR PADE APPROXIMATION   1 30 38 45
      DATA (GA( 7,13,IC),IC=1,3) /
     S 0.71369063E-01, 0.11204723E+01, 0.00000000E+00/
      DATA (GB( 7,13,IC),IC=1,3) /
     S 0.71369063E-01, 0.19367778E+01, 0.10000000E+01/
      DATA (GA( 7,14,IC),IC=1,3) /
     S 0.61866970E-01, 0.97740923E+00, 0.00000000E+00/
      DATA (GB( 7,14,IC),IC=1,3) /
     S 0.61866970E-01, 0.17112809E+01, 0.10000000E+01/
C
C----- INTERVAL = 2 ----- T =  275.0
C
C-- INDICES FOR PADE APPROXIMATION   1 30 38 45
      DATA (GA( 8,13,IC),IC=1,3) /
     S 0.69781812E-01, 0.10962918E+01, 0.00000000E+00/
      DATA (GB( 8,13,IC),IC=1,3) /
     S 0.69781812E-01, 0.18991112E+01, 0.10000000E+01/
      DATA (GA( 8,14,IC),IC=1,3) /
     S 0.60673632E-01, 0.96080188E+00, 0.00000000E+00/
      DATA (GB( 8,14,IC),IC=1,3) /
     S 0.60673632E-01, 0.16828137E+01, 0.10000000E+01/
C
C----- INTERVAL = 2 ----- T =  287.5
C
C-- INDICES FOR PADE APPROXIMATION   1 30 38 45
      DATA (GA( 9,13,IC),IC=1,3) /
     S 0.68381606E-01, 0.10752229E+01, 0.00000000E+00/
      DATA (GB( 9,13,IC),IC=1,3) /
     S 0.68381606E-01, 0.18658501E+01, 0.10000000E+01/
      DATA (GA( 9,14,IC),IC=1,3) /
     S 0.59637277E-01, 0.94657562E+00, 0.00000000E+00/
      DATA (GB( 9,14,IC),IC=1,3) /
     S 0.59637277E-01, 0.16580908E+01, 0.10000000E+01/
C
C----- INTERVAL = 2 ----- T =  300.0
C
C-- INDICES FOR PADE APPROXIMATION   1 30 38 45
      DATA (GA(10,13,IC),IC=1,3) /
     S 0.67139539E-01, 0.10567474E+01, 0.00000000E+00/
      DATA (GB(10,13,IC),IC=1,3) /
     S 0.67139539E-01, 0.18363226E+01, 0.10000000E+01/
      DATA (GA(10,14,IC),IC=1,3) /
     S 0.58732178E-01, 0.93430511E+00, 0.00000000E+00/
      DATA (GB(10,14,IC),IC=1,3) /
     S 0.58732178E-01, 0.16365014E+01, 0.10000000E+01/
C
C----- INTERVAL = 2 ----- T =  312.5
C
C-- INDICES FOR PADE APPROXIMATION   1 30 38 45
      DATA (GA(11,13,IC),IC=1,3) /
     S 0.66032012E-01, 0.10404465E+01, 0.00000000E+00/
      DATA (GB(11,13,IC),IC=1,3) /
     S 0.66032012E-01, 0.18099779E+01, 0.10000000E+01/
      DATA (GA(11,14,IC),IC=1,3) /
     S 0.57936092E-01, 0.92363528E+00, 0.00000000E+00/
      DATA (GB(11,14,IC),IC=1,3) /
     S 0.57936092E-01, 0.16175164E+01, 0.10000000E+01/
C
C
C
C
C
C
C
C
C
C
C-- CARBON DIOXIDE LINES IN THE WINDOW REGION (800-1250 CM-1)
C
C
C-- G = 0.0
C
C
C----- INTERVAL = 4 ----- T =  187.5
C
C-- INDICES FOR PADE APPROXIMATION     1   15   29   45
      DATA (GA( 1,15,IC),IC=1,3) /
     S 0.13230067E+02, 0.22042132E+02, 0.00000000E+00/
      DATA (GB( 1,15,IC),IC=1,3) /
     S 0.13230067E+02, 0.22051750E+02, 0.10000000E+01/
      DATA (GA( 1,16,IC),IC=1,3) /
     S 0.13183816E+02, 0.22169501E+02, 0.00000000E+00/
      DATA (GB( 1,16,IC),IC=1,3) /
     S 0.13183816E+02, 0.22178972E+02, 0.10000000E+01/
C
C----- INTERVAL = 4 ----- T =  200.0
C
C-- INDICES FOR PADE APPROXIMATION     1   15   29   45
      DATA (GA( 2,15,IC),IC=1,3) /
     S 0.13213564E+02, 0.22107298E+02, 0.00000000E+00/
      DATA (GB( 2,15,IC),IC=1,3) /
     S 0.13213564E+02, 0.22116850E+02, 0.10000000E+01/
      DATA (GA( 2,16,IC),IC=1,3) /
     S 0.13189991E+02, 0.22270075E+02, 0.00000000E+00/
      DATA (GB( 2,16,IC),IC=1,3) /
     S 0.13189991E+02, 0.22279484E+02, 0.10000000E+01/
C
C----- INTERVAL = 4 ----- T =  212.5
C
C-- INDICES FOR PADE APPROXIMATION     1   15   29   45
      DATA (GA( 3,15,IC),IC=1,3) /
     S 0.13209140E+02, 0.22180915E+02, 0.00000000E+00/
      DATA (GB( 3,15,IC),IC=1,3) /
     S 0.13209140E+02, 0.22190410E+02, 0.10000000E+01/
      DATA (GA( 3,16,IC),IC=1,3) /
     S 0.13209485E+02, 0.22379193E+02, 0.00000000E+00/
      DATA (GB( 3,16,IC),IC=1,3) /
     S 0.13209485E+02, 0.22388551E+02, 0.10000000E+01/
C
C----- INTERVAL = 4 ----- T =  225.0
C
C-- INDICES FOR PADE APPROXIMATION     1   15   29   45
      DATA (GA( 4,15,IC),IC=1,3) /
     S 0.13213894E+02, 0.22259478E+02, 0.00000000E+00/
      DATA (GB( 4,15,IC),IC=1,3) /
     S 0.13213894E+02, 0.22268925E+02, 0.10000000E+01/
      DATA (GA( 4,16,IC),IC=1,3) /
     S 0.13238789E+02, 0.22492992E+02, 0.00000000E+00/
      DATA (GB( 4,16,IC),IC=1,3) /
     S 0.13238789E+02, 0.22502309E+02, 0.10000000E+01/
C
C----- INTERVAL = 4 ----- T =  237.5
C
C-- INDICES FOR PADE APPROXIMATION     1   15   29   45
      DATA (GA( 5,15,IC),IC=1,3) /
     S 0.13225963E+02, 0.22341039E+02, 0.00000000E+00/
      DATA (GB( 5,15,IC),IC=1,3) /
     S 0.13225963E+02, 0.22350445E+02, 0.10000000E+01/
      DATA (GA( 5,16,IC),IC=1,3) /
     S 0.13275017E+02, 0.22608508E+02, 0.00000000E+00/
      DATA (GB( 5,16,IC),IC=1,3) /
     S 0.13275017E+02, 0.22617792E+02, 0.10000000E+01/
C
C----- INTERVAL = 4 ----- T =  250.0
C
C-- INDICES FOR PADE APPROXIMATION     1   15   29   45
      DATA (GA( 6,15,IC),IC=1,3) /
     S 0.13243806E+02, 0.22424247E+02, 0.00000000E+00/
      DATA (GB( 6,15,IC),IC=1,3) /
     S 0.13243806E+02, 0.22433617E+02, 0.10000000E+01/
      DATA (GA( 6,16,IC),IC=1,3) /
     S 0.13316096E+02, 0.22723843E+02, 0.00000000E+00/
      DATA (GB( 6,16,IC),IC=1,3) /
     S 0.13316096E+02, 0.22733099E+02, 0.10000000E+01/
C
C----- INTERVAL = 4 ----- T =  262.5
C
C-- INDICES FOR PADE APPROXIMATION     1   15   29   45
      DATA (GA( 7,15,IC),IC=1,3) /
     S 0.13266104E+02, 0.22508089E+02, 0.00000000E+00/
      DATA (GB( 7,15,IC),IC=1,3) /
     S 0.13266104E+02, 0.22517429E+02, 0.10000000E+01/
      DATA (GA( 7,16,IC),IC=1,3) /
     S 0.13360555E+02, 0.22837837E+02, 0.00000000E+00/
      DATA (GB( 7,16,IC),IC=1,3) /
     S 0.13360555E+02, 0.22847071E+02, 0.10000000E+01/
C
C----- INTERVAL = 4 ----- T =  275.0
C
C-- INDICES FOR PADE APPROXIMATION     1   15   29   45
      DATA (GA( 8,15,IC),IC=1,3) /
     S 0.13291782E+02, 0.22591771E+02, 0.00000000E+00/
      DATA (GB( 8,15,IC),IC=1,3) /
     S 0.13291782E+02, 0.22601086E+02, 0.10000000E+01/
      DATA (GA( 8,16,IC),IC=1,3) /
     S 0.13407324E+02, 0.22949751E+02, 0.00000000E+00/
      DATA (GB( 8,16,IC),IC=1,3) /
     S 0.13407324E+02, 0.22958967E+02, 0.10000000E+01/
C
C----- INTERVAL = 4 ----- T =  287.5
C
C-- INDICES FOR PADE APPROXIMATION     1   15   29   45
      DATA (GA( 9,15,IC),IC=1,3) /
     S 0.13319961E+02, 0.22674661E+02, 0.00000000E+00/
      DATA (GB( 9,15,IC),IC=1,3) /
     S 0.13319961E+02, 0.22683956E+02, 0.10000000E+01/
      DATA (GA( 9,16,IC),IC=1,3) /
     S 0.13455544E+02, 0.23059032E+02, 0.00000000E+00/
      DATA (GB( 9,16,IC),IC=1,3) /
     S 0.13455544E+02, 0.23068234E+02, 0.10000000E+01/
C
C----- INTERVAL = 4 ----- T =  300.0
C
C-- INDICES FOR PADE APPROXIMATION     1   15   29   45
      DATA (GA(10,15,IC),IC=1,3) /
     S 0.13349927E+02, 0.22756246E+02, 0.00000000E+00/
      DATA (GB(10,15,IC),IC=1,3) /
     S 0.13349927E+02, 0.22765522E+02, 0.10000000E+01/
      DATA (GA(10,16,IC),IC=1,3) /
     S 0.13504450E+02, 0.23165146E+02, 0.00000000E+00/
      DATA (GB(10,16,IC),IC=1,3) /
     S 0.13504450E+02, 0.23174336E+02, 0.10000000E+01/
C
C----- INTERVAL = 4 ----- T =  312.5
C
C-- INDICES FOR PADE APPROXIMATION     1   15   29   45
      DATA (GA(11,15,IC),IC=1,3) /
     S 0.13381108E+02, 0.22836093E+02, 0.00000000E+00/
      DATA (GB(11,15,IC),IC=1,3) /
     S 0.13381108E+02, 0.22845354E+02, 0.10000000E+01/
      DATA (GA(11,16,IC),IC=1,3) /
     S 0.13553282E+02, 0.23267456E+02, 0.00000000E+00/
      DATA (GB(11,16,IC),IC=1,3) /
     S 0.13553282E+02, 0.23276638E+02, 0.10000000E+01/
C     ------------------------------------------------------------------
      DATA (( XP(  J,K),J=1,6),       K=1,6) /
     S 0.46430621E+02, 0.12928299E+03, 0.20732648E+03,
     S 0.31398411E+03, 0.18373177E+03,-0.11412303E+03,
     S 0.73604774E+02, 0.27887914E+03, 0.27076947E+03,
     S-0.57322111E+02,-0.64742459E+02, 0.87238280E+02,
     S 0.37050866E+02, 0.20498759E+03, 0.37558029E+03,
     S 0.17401171E+03,-0.13350302E+03,-0.37651795E+02,
     S 0.14930141E+02, 0.89161160E+02, 0.17793062E+03,
     S 0.93433860E+02,-0.70646020E+02,-0.26373150E+02,
     S 0.40386780E+02, 0.10855270E+03, 0.50755010E+02,
     S-0.31496190E+02, 0.12791300E+00, 0.18017770E+01,
     S 0.90811926E+01, 0.75073923E+02, 0.24654438E+03,
     S 0.39332612E+03, 0.29385281E+03, 0.89107921E+02 /
C
C
C*         1.0     PLANCK FUNCTIONS AND GRADIENTS
C                  ------------------------------
C
 100  CONTINUE
C
      DO 102 JK = 1 , KFLEV+1
      DO 101 JL = 1, KDLON
      PBINT(JL,JK) = 0.
 101  CONTINUE
 102  CONTINUE
      DO 103 JL = 1, KDLON
      PBSUIN(JL) = 0.
 103  CONTINUE
C
      DO 141 JNU=1,Ninter
C
C
C*         1.1   LEVELS FROM SURFACE TO KFLEV
C                ----------------------------
C
 110  CONTINUE
C
      DO 112 JK = 1 , KFLEV
      DO 111 JL = 1, KDLON
      ZTI(JL)=(PTL(JL,JK)-TSTAND)/TSTAND
      ZRES(JL) = XP(1,JNU)+ZTI(JL)*(XP(2,JNU)+ZTI(JL)*(XP(3,JNU)
     S       +ZTI(JL)*(XP(4,JNU)+ZTI(JL)*(XP(5,JNU)+ZTI(JL)*(XP(6,JNU)
     S       )))))
      PBINT(JL,JK)=PBINT(JL,JK)+ZRES(JL)
      PB(JL,JNU,JK)= ZRES(JL)
      ZBLEV(JL,JK) = ZRES(JL)
      ZTI2(JL)=(PTAVE(JL,JK)-TSTAND)/TSTAND
      ZRES2(JL)=XP(1,JNU)+ZTI2(JL)*(XP(2,JNU)+ZTI2(JL)*(XP(3,JNU)
     S     +ZTI2(JL)*(XP(4,JNU)+ZTI2(JL)*(XP(5,JNU)+ZTI2(JL)*(XP(6,JNU)
     S       )))))
      ZBLAY(JL,JK) = ZRES2(JL)
 111  CONTINUE
 112  CONTINUE
C
C
C*         1.2   TOP OF THE ATMOSPHERE AND SURFACE
C                ---------------------------------
C
 120  CONTINUE
C
      DO 121 JL = 1, KDLON
      ZTI(JL)=(PTL(JL,KFLEV+1)-TSTAND)/TSTAND
      ZTI2(JL) = (PTL(JL,1) + PDT0(JL) - TSTAND) / TSTAND
      ZRES(JL) = XP(1,JNU)+ZTI(JL)*(XP(2,JNU)+ZTI(JL)*(XP(3,JNU)
     S    +ZTI(JL)*(XP(4,JNU)+ZTI(JL)*(XP(5,JNU)+ZTI(JL)*(XP(6,JNU)
     S       )))))
      ZRES2(JL) = XP(1,JNU)+ZTI2(JL)*(XP(2,JNU)+ZTI2(JL)*(XP(3,JNU)
     S    +ZTI2(JL)*(XP(4,JNU)+ZTI2(JL)*(XP(5,JNU)+ZTI2(JL)*(XP(6,JNU)
     S       )))))
      PBINT(JL,KFLEV+1) = PBINT(JL,KFLEV+1)+ZRES(JL)
      PB(JL,JNU,KFLEV+1)= ZRES(JL)
      ZBLEV(JL,KFLEV+1) = ZRES(JL)
      PBTOP(JL,JNU) = ZRES(JL)
      PBSUR(JL,JNU) = ZRES2(JL)
      PBSUIN(JL) = PBSUIN(JL) + ZRES2(JL)
 121  CONTINUE
C
C
C*         1.3   GRADIENTS IN SUB-LAYERS
C                -----------------------
C
 130  CONTINUE
C
      DO 132 JK = 1 , KFLEV
      JK2 = 2 * JK
      JK1 = JK2 - 1
      DO 131 JL = 1, KDLON
      PDBSL(JL,JNU,JK1) = ZBLAY(JL,JK  ) - ZBLEV(JL,JK)
      PDBSL(JL,JNU,JK2) = ZBLEV(JL,JK+1) - ZBLAY(JL,JK)
 131  CONTINUE
 132  CONTINUE
C
 141  CONTINUE
C
C*         2.0   CHOOSE THE RELEVANT SETS OF PADE APPROXIMANTS
C                ---------------------------------------------
C
 200  CONTINUE
C
C
 210  CONTINUE
C
      DO 211 JL=1, KDLON
      ZDSTO1 = (PTL(JL,KFLEV+1)-TINTP(1)) / TSTP
      IXTOX = MAX( 1, MIN( MXIXT, INT( ZDSTO1 + 1. ) ) )
      ZDSTOX = (PTL(JL,KFLEV+1)-TINTP(IXTOX))/TSTP
      IF (ZDSTOX.LT.0.5) THEN
         INDTO=IXTOX
      ELSE
         INDTO=IXTOX+1
      END IF
      INDB(JL)=INDTO
      ZDST1 = (PTL(JL,1)-TINTP(1)) / TSTP
      IXTX = MAX( 1, MIN( MXIXT, INT( ZDST1 + 1. ) ) )
      ZDSTX = (PTL(JL,1)-TINTP(IXTX))/TSTP
      IF (ZDSTX.LT.0.5) THEN
         INDT=IXTX
      ELSE
         INDT=IXTX+1
      END IF
      INDS(JL)=INDT
 211  CONTINUE
C
      DO 214 JF=1,2
      DO 213 JG=1, 8
      DO 212 JL=1, KDLON
      INDSU=INDS(JL)
      PGASUR(JL,JG,JF)=GA(INDSU,2*JG-1,JF)
      PGBSUR(JL,JG,JF)=GB(INDSU,2*JG-1,JF)
      INDTP=INDB(JL)
      PGATOP(JL,JG,JF)=GA(INDTP,2*JG-1,JF)
      PGBTOP(JL,JG,JF)=GB(INDTP,2*JG-1,JF)
 212  CONTINUE
 213  CONTINUE
 214  CONTINUE
C
 220  CONTINUE
C
      DO 225 JK=1,KFLEV
      DO 221 JL=1, KDLON
      ZDST1 = (PTAVE(JL,JK)-TINTP(1)) / TSTP
      IXTX = MAX( 1, MIN( MXIXT, INT( ZDST1 + 1. ) ) )
      ZDSTX = (PTAVE(JL,JK)-TINTP(IXTX))/TSTP
      IF (ZDSTX.LT.0.5) THEN
         INDT=IXTX
      ELSE
         INDT=IXTX+1
      END IF
      INDB(JL)=INDT
 221  CONTINUE
C
      DO 224 JF=1,2
      DO 223 JG=1, 8
      DO 222 JL=1, KDLON
      INDT=INDB(JL)
      PGA(JL,JG,JF,JK)=GA(INDT,2*JG,JF)
      PGB(JL,JG,JF,JK)=GB(INDT,2*JG,JF)
 222  CONTINUE
 223  CONTINUE
 224  CONTINUE
 225  CONTINUE
C
C     ------------------------------------------------------------------
C
      RETURN
      END
      SUBROUTINE LWV(KUAER,KTRAER, KLIM
     R  , PABCU,PB,PBINT,PBSUIN,PBSUR,PBTOP,PDBSL,PEMIS,PPMB,PTAVE
     R  , PGA,PGB,PGASUR,PGBSUR,PGATOP,PGBTOP
     S  , PCNTRB,PCTS,PFLUC)
      IMPLICIT none
      include "dimensions.h"
      include "dimphy.h"
      include "raddim.h"
      include "raddimlw.h"
      include "YOMCST.h"
C
C-----------------------------------------------------------------------
C     PURPOSE.
C     --------
C           CARRIES OUT THE VERTICAL INTEGRATION TO GIVE LONGWAVE
C           FLUXES OR RADIANCES
C
C     METHOD.
C     -------
C
C          1. PERFORMS THE VERTICAL INTEGRATION DISTINGUISHING BETWEEN
C     CONTRIBUTIONS BY -  THE NEARBY LAYERS
C                      -  THE DISTANT LAYERS
C                      -  THE BOUNDARY TERMS
C          2. COMPUTES THE CLEAR-SKY DOWNWARD AND UPWARD EMISSIVITIES.
C
C     REFERENCE.
C     ----------
C
C        SEE RADIATION'S PART OF THE MODEL'S DOCUMENTATION AND
C        ECMWF RESEARCH DEPARTMENT DOCUMENTATION OF THE IFS
C
C     AUTHOR.
C     -------
C        JEAN-JACQUES MORCRETTE  *ECMWF*
C
C     MODIFICATIONS.
C     --------------
C        ORIGINAL : 89-07-14
C-----------------------------------------------------------------------
C
C* ARGUMENTS:
      INTEGER KUAER,KTRAER, KLIM
C
      REAL*8 PABCU(KDLON,NUA,3*KFLEV+1) ! EFFECTIVE ABSORBER AMOUNTS
      REAL*8 PB(KDLON,Ninter,KFLEV+1) ! SPECTRAL HALF-LEVEL PLANCK FUNCTIONS
      REAL*8 PBINT(KDLON,KFLEV+1) ! HALF-LEVEL PLANCK FUNCTIONS
      REAL*8 PBSUR(KDLON,Ninter) ! SURFACE SPECTRAL PLANCK FUNCTION
      REAL*8 PBSUIN(KDLON) ! SURFACE PLANCK FUNCTION
      REAL*8 PBTOP(KDLON,Ninter) ! T.O.A. SPECTRAL PLANCK FUNCTION
      REAL*8 PDBSL(KDLON,Ninter,KFLEV*2) ! SUB-LAYER PLANCK FUNCTION GRADIENT
      REAL*8 PEMIS(KDLON) ! SURFACE EMISSIVITY
      REAL*8 PPMB(KDLON,KFLEV+1) ! HALF-LEVEL PRESSURE (MB)
      REAL*8 PTAVE(KDLON,KFLEV) ! TEMPERATURE
      REAL*8 PGA(KDLON,8,2,KFLEV) ! PADE APPROXIMANTS
      REAL*8 PGB(KDLON,8,2,KFLEV) ! PADE APPROXIMANTS
      REAL*8 PGASUR(KDLON,8,2) ! PADE APPROXIMANTS
      REAL*8 PGBSUR(KDLON,8,2) ! PADE APPROXIMANTS
      REAL*8 PGATOP(KDLON,8,2) ! PADE APPROXIMANTS
      REAL*8 PGBTOP(KDLON,8,2) ! PADE APPROXIMANTS
C
      REAL*8 PCNTRB(KDLON,KFLEV+1,KFLEV+1) ! CLEAR-SKY ENERGY EXCHANGE MATRIX
      REAL*8 PCTS(KDLON,KFLEV) ! COOLING-TO-SPACE TERM
      REAL*8 PFLUC(KDLON,2,KFLEV+1) ! CLEAR-SKY RADIATIVE FLUXES
C-----------------------------------------------------------------------
C LOCAL VARIABLES:
      REAL*8 ZADJD(KDLON,KFLEV+1)
      REAL*8 ZADJU(KDLON,KFLEV+1)
      REAL*8 ZDBDT(KDLON,Ninter,KFLEV)
      REAL*8 ZDISD(KDLON,KFLEV+1)
      REAL*8 ZDISU(KDLON,KFLEV+1)
C
      INTEGER jk, jl
C-----------------------------------------------------------------------
C
      DO 112 JK=1,KFLEV+1
      DO 111 JL=1, KDLON
      ZADJD(JL,JK)=0.
      ZADJU(JL,JK)=0.
      ZDISD(JL,JK)=0.
      ZDISU(JL,JK)=0.
 111  CONTINUE
 112  CONTINUE
C
      DO 114 JK=1,KFLEV
      DO 113 JL=1, KDLON
      PCTS(JL,JK)=0.
 113  CONTINUE
 114  CONTINUE
C
C* CONTRIBUTION FROM ADJACENT LAYERS
C
      CALL LWVN(KUAER,KTRAER
     R  , PABCU,PDBSL,PGA,PGB
     S  , ZADJD,ZADJU,PCNTRB,ZDBDT)
C* CONTRIBUTION FROM DISTANT LAYERS
C
      CALL LWVD(KUAER,KTRAER
     R  , PABCU,ZDBDT,PGA,PGB
     S  , PCNTRB,ZDISD,ZDISU)
C
C* EXCHANGE WITH THE BOUNDARIES
C
      CALL LWVB(KUAER,KTRAER, KLIM
     R  , PABCU,ZADJD,ZADJU,PB,PBINT,PBSUIN,PBSUR,PBTOP
     R  , ZDISD,ZDISU,PEMIS,PPMB
     R  , PGA,PGB,PGASUR,PGBSUR,PGATOP,PGBTOP
     S  , PCTS,PFLUC)
C
C
      RETURN
      END
      SUBROUTINE LWVB(KUAER,KTRAER, KLIM
     R  , PABCU,PADJD,PADJU,PB,PBINT,PBSUI,PBSUR,PBTOP
     R  , PDISD,PDISU,PEMIS,PPMB
     R  , PGA,PGB,PGASUR,PGBSUR,PGATOP,PGBTOP
     S  , PCTS,PFLUC)
      IMPLICIT none
      include "dimensions.h"
      include "dimphy.h"
      include "raddim.h"
      include "raddimlw.h"
      include "radopt.h"
C
C-----------------------------------------------------------------------
C     PURPOSE.
C     --------
C           INTRODUCES THE EFFECTS OF THE BOUNDARIES IN THE VERTICAL
C           INTEGRATION
C
C     METHOD.
C     -------
C
C          1. COMPUTES THE ENERGY EXCHANGE WITH TOP AND SURFACE OF THE
C     ATMOSPHERE
C          2. COMPUTES THE COOLING-TO-SPACE AND HEATING-FROM-GROUND
C     TERMS FOR THE APPROXIMATE COOLING RATE ABOVE 10 HPA
C          3. ADDS UP ALL CONTRIBUTIONS TO GET THE CLEAR-SKY FLUXES
C
C     REFERENCE.
C     ----------
C
C        SEE RADIATION'S PART OF THE MODEL'S DOCUMENTATION AND
C        ECMWF RESEARCH DEPARTMENT DOCUMENTATION OF THE IFS
C
C     AUTHOR.
C     -------
C        JEAN-JACQUES MORCRETTE  *ECMWF*
C
C     MODIFICATIONS.
C     --------------
C        ORIGINAL : 89-07-14
C        Voigt lines (loop 2413 to 2427)  - JJM & PhD - 01/96
C-----------------------------------------------------------------------
C
C*       0.1   ARGUMENTS
C              ---------
C
      INTEGER KUAER,KTRAER, KLIM
C
      REAL*8 PABCU(KDLON,NUA,3*KFLEV+1) ! ABSORBER AMOUNTS
      REAL*8 PADJD(KDLON,KFLEV+1) ! CONTRIBUTION BY ADJACENT LAYERS
      REAL*8 PADJU(KDLON,KFLEV+1) ! CONTRIBUTION BY ADJACENT LAYERS
      REAL*8 PB(KDLON,Ninter,KFLEV+1) ! SPECTRAL HALF-LEVEL PLANCK FUNCTIONS
      REAL*8 PBINT(KDLON,KFLEV+1) ! HALF-LEVEL PLANCK FUNCTIONS
      REAL*8 PBSUR(KDLON,Ninter) ! SPECTRAL SURFACE PLANCK FUNCTION
      REAL*8 PBSUI(KDLON) ! SURFACE PLANCK FUNCTION
      REAL*8 PBTOP(KDLON,Ninter) ! SPECTRAL T.O.A. PLANCK FUNCTION
      REAL*8 PDISD(KDLON,KFLEV+1) ! CONTRIBUTION BY DISTANT LAYERS
      REAL*8 PDISU(KDLON,KFLEV+1) ! CONTRIBUTION BY DISTANT LAYERS
      REAL*8 PEMIS(KDLON) ! SURFACE EMISSIVITY
      REAL*8 PPMB(KDLON,KFLEV+1) ! PRESSURE MB
      REAL*8 PGA(KDLON,8,2,KFLEV) ! PADE APPROXIMANTS
      REAL*8 PGB(KDLON,8,2,KFLEV) ! PADE APPROXIMANTS
      REAL*8 PGASUR(KDLON,8,2) ! SURFACE PADE APPROXIMANTS
      REAL*8 PGBSUR(KDLON,8,2) ! SURFACE PADE APPROXIMANTS
      REAL*8 PGATOP(KDLON,8,2) ! T.O.A. PADE APPROXIMANTS
      REAL*8 PGBTOP(KDLON,8,2) ! T.O.A. PADE APPROXIMANTS
C
      REAL*8 PFLUC(KDLON,2,KFLEV+1) ! CLEAR-SKY RADIATIVE FLUXES
      REAL*8 PCTS(KDLON,KFLEV) ! COOLING-TO-SPACE TERM
C
C* LOCAL VARIABLES:
C
      REAL*8 ZBGND(KDLON)
      REAL*8 ZFD(KDLON)
      REAL*8  ZFN10(KDLON)
      REAL*8 ZFU(KDLON)
      REAL*8  ZTT(KDLON,NTRA)
      REAL*8 ZTT1(KDLON,NTRA)
      REAL*8 ZTT2(KDLON,NTRA)
      REAL*8  ZUU(KDLON,NUA) 
      REAL*8 ZCNSOL(KDLON)
      REAL*8 ZCNTOP(KDLON)
C
      INTEGER jk, jl, ja
      INTEGER jstra, jstru
      INTEGER ind1, ind2, ind3, ind4, in, jlim
      REAL*8 zctstr
C-----------------------------------------------------------------------
C
C*         1.    INITIALIZATION
C                --------------
C
 100  CONTINUE
C
C
C*         1.2     INITIALIZE TRANSMISSION FUNCTIONS
C                  ---------------------------------
C
 120  CONTINUE
C
      DO 122 JA=1,NTRA
      DO 121 JL=1, KDLON
      ZTT (JL,JA)=1.0
      ZTT1(JL,JA)=1.0
      ZTT2(JL,JA)=1.0
 121  CONTINUE
 122  CONTINUE
C
      DO 124 JA=1,NUA
      DO 123 JL=1, KDLON
      ZUU(JL,JA)=1.0
 123  CONTINUE
 124  CONTINUE
C
C     ------------------------------------------------------------------
C
C*         2.      VERTICAL INTEGRATION
C                  --------------------
C
 200  CONTINUE
C
      IND1=0
      IND3=0
      IND4=1
      IND2=1
C
C
C*         2.3     EXCHANGE WITH TOP OF THE ATMOSPHERE
C                  -----------------------------------
C
 230  CONTINUE
C
      DO 235 JK = 1 , KFLEV
      IN=(JK-1)*NG1P1+1
C
      DO 232 JA=1,KUAER
      DO 231 JL=1, KDLON
      ZUU(JL,JA)=PABCU(JL,JA,IN)
 231  CONTINUE
 232  CONTINUE
C
C
      CALL LWTT(PGATOP(1,1,1), PGBTOP(1,1,1), ZUU, ZTT)
C
      DO 234 JL = 1, KDLON
      ZCNTOP(JL)=PBTOP(JL,1)*ZTT(JL,1)          *ZTT(JL,10)
     2      +PBTOP(JL,2)*ZTT(JL,2)*ZTT(JL,7)*ZTT(JL,11)
     3      +PBTOP(JL,3)*ZTT(JL,4)*ZTT(JL,8)*ZTT(JL,12)
     4      +PBTOP(JL,4)*ZTT(JL,5)*ZTT(JL,9)*ZTT(JL,13)
     5      +PBTOP(JL,5)*ZTT(JL,3)          *ZTT(JL,14)
     6      +PBTOP(JL,6)*ZTT(JL,6)          *ZTT(JL,15)
      ZFD(JL)=ZCNTOP(JL)-PBINT(JL,JK)-PDISD(JL,JK)-PADJD(JL,JK)
      PFLUC(JL,2,JK)=ZFD(JL)
 234  CONTINUE
C
 235  CONTINUE
C
      JK = KFLEV+1
      IN=(JK-1)*NG1P1+1
C
      DO 236 JL = 1, KDLON
      ZCNTOP(JL)= PBTOP(JL,1)
     1   + PBTOP(JL,2)
     2   + PBTOP(JL,3)
     3   + PBTOP(JL,4)
     4   + PBTOP(JL,5)
     5   + PBTOP(JL,6)
      ZFD(JL)=ZCNTOP(JL)-PBINT(JL,JK)-PDISD(JL,JK)-PADJD(JL,JK)
      PFLUC(JL,2,JK)=ZFD(JL)
 236  CONTINUE
C
C*         2.4     COOLING-TO-SPACE OF LAYERS ABOVE 10 HPA
C                  ---------------------------------------
C
 240  CONTINUE
C
C
C*         2.4.1   INITIALIZATION
C                  --------------
C
 2410 CONTINUE
C
      JLIM = KFLEV
C
      IF (.NOT.LEVOIGT) THEN
      DO 2412 JK = KFLEV,1,-1
      IF(PPMB(1,JK).LT.10.0) THEN
         JLIM=JK
      ENDIF   
 2412 CONTINUE
      ENDIF
      KLIM=JLIM
C
      IF (.NOT.LEVOIGT) THEN
        DO 2414 JA=1,KTRAER
        DO 2413 JL=1, KDLON
        ZTT1(JL,JA)=1.0
 2413   CONTINUE
 2414   CONTINUE
C
C*         2.4.2   LOOP OVER LAYERS ABOVE 10 HPA
C                  -----------------------------
C
 2420   CONTINUE
C
        DO 2427 JSTRA = KFLEV,JLIM,-1
        JSTRU=(JSTRA-1)*NG1P1+1
C
        DO 2423 JA=1,KUAER
        DO 2422 JL=1, KDLON
        ZUU(JL,JA)=PABCU(JL,JA,JSTRU)
 2422   CONTINUE
 2423   CONTINUE
C
C
        CALL LWTT(PGA(1,1,1,JSTRA), PGB(1,1,1,JSTRA), ZUU, ZTT)
C
        DO 2424 JL = 1, KDLON
        ZCTSTR =
     1   (PB(JL,1,JSTRA)+PB(JL,1,JSTRA+1))
     1       *(ZTT1(JL,1)           *ZTT1(JL,10)
     1       - ZTT (JL,1)           *ZTT (JL,10))
     2  +(PB(JL,2,JSTRA)+PB(JL,2,JSTRA+1))
     2       *(ZTT1(JL,2)*ZTT1(JL,7)*ZTT1(JL,11)
     2       - ZTT (JL,2)*ZTT (JL,7)*ZTT (JL,11))
     3  +(PB(JL,3,JSTRA)+PB(JL,3,JSTRA+1))
     3       *(ZTT1(JL,4)*ZTT1(JL,8)*ZTT1(JL,12)
     3       - ZTT (JL,4)*ZTT (JL,8)*ZTT (JL,12))
     4  +(PB(JL,4,JSTRA)+PB(JL,4,JSTRA+1))
     4       *(ZTT1(JL,5)*ZTT1(JL,9)*ZTT1(JL,13)
     4       - ZTT (JL,5)*ZTT (JL,9)*ZTT (JL,13))
     5  +(PB(JL,5,JSTRA)+PB(JL,5,JSTRA+1))
     5       *(ZTT1(JL,3)           *ZTT1(JL,14)
     5       - ZTT (JL,3)           *ZTT (JL,14))
     6  +(PB(JL,6,JSTRA)+PB(JL,6,JSTRA+1))
     6       *(ZTT1(JL,6)           *ZTT1(JL,15)
     6       - ZTT (JL,6)           *ZTT (JL,15))
        PCTS(JL,JSTRA)=ZCTSTR*0.5
 2424   CONTINUE
        DO 2426 JA=1,KTRAER
        DO 2425 JL=1, KDLON
        ZTT1(JL,JA)=ZTT(JL,JA)
 2425   CONTINUE
 2426   CONTINUE
 2427   CONTINUE
      ENDIF
C Mise a zero de securite pour PCTS en cas de LEVOIGT
      IF(LEVOIGT)THEN
        DO 2429 JSTRA = 1,KFLEV
        DO 2428 JL = 1, KDLON
          PCTS(JL,JSTRA)=0.
 2428   CONTINUE
 2429   CONTINUE
      ENDIF
C
C
C*         2.5     EXCHANGE WITH LOWER LIMIT
C                  -------------------------
C
 250  CONTINUE
C
      DO 251 JL = 1, KDLON
      ZBGND(JL)=PBSUI(JL)*PEMIS(JL)-(1.-PEMIS(JL))
     S               *PFLUC(JL,2,1)-PBINT(JL,1)
 251  CONTINUE
C
      JK = 1
      IN=(JK-1)*NG1P1+1
C
      DO 252 JL = 1, KDLON
      ZCNSOL(JL)=PBSUR(JL,1)
     1 +PBSUR(JL,2)
     2 +PBSUR(JL,3)
     3 +PBSUR(JL,4)
     4 +PBSUR(JL,5)
     5 +PBSUR(JL,6)
      ZCNSOL(JL)=ZCNSOL(JL)*ZBGND(JL)/PBSUI(JL)
      ZFU(JL)=ZCNSOL(JL)+PBINT(JL,JK)-PDISU(JL,JK)-PADJU(JL,JK)
      PFLUC(JL,1,JK)=ZFU(JL)
 252  CONTINUE
C
      DO 257 JK = 2 , KFLEV+1
      IN=(JK-1)*NG1P1+1
C
C
      DO 255 JA=1,KUAER
      DO 254 JL=1, KDLON
      ZUU(JL,JA)=PABCU(JL,JA,1)-PABCU(JL,JA,IN)
 254  CONTINUE
 255  CONTINUE
C
C
      CALL LWTT(PGASUR(1,1,1), PGBSUR(1,1,1), ZUU, ZTT)
C
      DO 256 JL = 1, KDLON
      ZCNSOL(JL)=PBSUR(JL,1)*ZTT(JL,1)          *ZTT(JL,10)
     2      +PBSUR(JL,2)*ZTT(JL,2)*ZTT(JL,7)*ZTT(JL,11)
     3      +PBSUR(JL,3)*ZTT(JL,4)*ZTT(JL,8)*ZTT(JL,12)
     4      +PBSUR(JL,4)*ZTT(JL,5)*ZTT(JL,9)*ZTT(JL,13)
     5      +PBSUR(JL,5)*ZTT(JL,3)          *ZTT(JL,14)
     6      +PBSUR(JL,6)*ZTT(JL,6)          *ZTT(JL,15)
      ZCNSOL(JL)=ZCNSOL(JL)*ZBGND(JL)/PBSUI(JL)
      ZFU(JL)=ZCNSOL(JL)+PBINT(JL,JK)-PDISU(JL,JK)-PADJU(JL,JK)
      PFLUC(JL,1,JK)=ZFU(JL)
 256  CONTINUE
C
C
 257  CONTINUE
C
C
C
C*         2.7     CLEAR-SKY FLUXES
C                  ----------------
C
 270  CONTINUE
C
      IF (.NOT.LEVOIGT) THEN
      DO 271 JL = 1, KDLON
      ZFN10(JL) = PFLUC(JL,1,JLIM) + PFLUC(JL,2,JLIM)
 271  CONTINUE
      DO 273 JK = JLIM+1,KFLEV+1
      DO 272 JL = 1, KDLON
      ZFN10(JL) = ZFN10(JL) + PCTS(JL,JK-1)
      PFLUC(JL,1,JK) = ZFN10(JL)
      PFLUC(JL,2,JK) = 0.
 272  CONTINUE
 273  CONTINUE
      ENDIF
C
C     ------------------------------------------------------------------
C
      RETURN
      END
      SUBROUTINE LWVD(KUAER,KTRAER
     S  , PABCU,PDBDT
     R  , PGA,PGB
     S  , PCNTRB,PDISD,PDISU)
      IMPLICIT none
      include "dimensions.h"
      include "dimphy.h"
      include "raddim.h"
      include "raddimlw.h"
C
C-----------------------------------------------------------------------
C     PURPOSE.
C     --------
C           CARRIES OUT THE VERTICAL INTEGRATION ON THE DISTANT LAYERS
C
C     METHOD.
C     -------
C
C          1. PERFORMS THE VERTICAL INTEGRATION CORRESPONDING TO THE
C     CONTRIBUTIONS OF THE DISTANT LAYERS USING TRAPEZOIDAL RULE
C
C     REFERENCE.
C     ----------
C
C        SEE RADIATION'S PART OF THE MODEL'S DOCUMENTATION AND
C        ECMWF RESEARCH DEPARTMENT DOCUMENTATION OF THE IFS
C
C     AUTHOR.
C     -------
C        JEAN-JACQUES MORCRETTE  *ECMWF*
C
C     MODIFICATIONS.
C     --------------
C        ORIGINAL : 89-07-14
C-----------------------------------------------------------------------
C* ARGUMENTS:
C
      INTEGER KUAER,KTRAER
C
      REAL*8 PABCU(KDLON,NUA,3*KFLEV+1) ! ABSORBER AMOUNTS
      REAL*8 PDBDT(KDLON,Ninter,KFLEV) ! LAYER PLANCK FUNCTION GRADIENT
      REAL*8 PGA(KDLON,8,2,KFLEV) ! PADE APPROXIMANTS
      REAL*8 PGB(KDLON,8,2,KFLEV) ! PADE APPROXIMANTS
C
      REAL*8 PCNTRB(KDLON,KFLEV+1,KFLEV+1) ! ENERGY EXCHANGE MATRIX
      REAL*8 PDISD(KDLON,KFLEV+1) !  CONTRIBUTION BY DISTANT LAYERS
      REAL*8 PDISU(KDLON,KFLEV+1) !  CONTRIBUTION BY DISTANT LAYERS
C
C* LOCAL VARIABLES:
C
      REAL*8 ZGLAYD(KDLON)
      REAL*8 ZGLAYU(KDLON)
      REAL*8 ZTT(KDLON,NTRA)
      REAL*8 ZTT1(KDLON,NTRA)
      REAL*8 ZTT2(KDLON,NTRA)
C
      INTEGER jl, jk, ja, ikp1, ikn, ikd1, jkj, ikd2
      INTEGER ikjp1, ikm1, ikj, jlk, iku1, ijkl, iku2
      INTEGER ind1, ind2, ind3, ind4, itt
      REAL*8 zww, zdzxdg, zdzxmg
C
C*         1.    INITIALIZATION
C                --------------
C
 100  CONTINUE
C
C*         1.1     INITIALIZE LAYER CONTRIBUTIONS
C                  ------------------------------
C
 110  CONTINUE
C
      DO 112 JK = 1, KFLEV+1
      DO 111 JL = 1, KDLON
      PDISD(JL,JK) = 0.
      PDISU(JL,JK) = 0.
  111 CONTINUE
  112 CONTINUE
C
C*         1.2     INITIALIZE TRANSMISSION FUNCTIONS
C                  ---------------------------------
C
 120  CONTINUE
C
C
      DO 122 JA = 1, NTRA
      DO 121 JL = 1, KDLON
      ZTT (JL,JA) = 1.0
      ZTT1(JL,JA) = 1.0
      ZTT2(JL,JA) = 1.0
  121 CONTINUE
  122 CONTINUE
C
C     ------------------------------------------------------------------
C
C*         2.      VERTICAL INTEGRATION
C                  --------------------
C
 200  CONTINUE
C
      IND1=0
      IND3=0
      IND4=1
      IND2=1
C
C
C*         2.2     CONTRIBUTION FROM DISTANT LAYERS
C                  ---------------------------------
C
 220  CONTINUE
C
C
C*         2.2.1   DISTANT AND ABOVE LAYERS
C                  ------------------------
C
 2210 CONTINUE
C
C
C
C*         2.2.2   FIRST UPPER LEVEL
C                  -----------------
C
 2220 CONTINUE
C
      DO 225 JK = 1 , KFLEV-1
      IKP1=JK+1
      IKN=(JK-1)*NG1P1+1
      IKD1= JK  *NG1P1+1
C
      CALL LWTTM(PGA(1,1,1,JK), PGB(1,1,1,JK)
     2          , PABCU(1,1,IKN),PABCU(1,1,IKD1),ZTT1)
C
C
C
C*         2.2.3   HIGHER UP
C                  ---------
C
 2230 CONTINUE
C
      ITT=1
      DO 224 JKJ=IKP1,KFLEV
      IF(ITT.EQ.1) THEN
         ITT=2
      ELSE
         ITT=1
      ENDIF
      IKJP1=JKJ+1
      IKD2= JKJ  *NG1P1+1
C
      IF(ITT.EQ.1) THEN
         CALL LWTTM(PGA(1,1,1,JKJ),PGB(1,1,1,JKJ)
     2             , PABCU(1,1,IKN),PABCU(1,1,IKD2),ZTT1)
      ELSE
         CALL LWTTM(PGA(1,1,1,JKJ),PGB(1,1,1,JKJ)
     2             , PABCU(1,1,IKN),PABCU(1,1,IKD2),ZTT2)
      ENDIF
C
      DO 2235 JA = 1, KTRAER
      DO 2234 JL = 1, KDLON
      ZTT(JL,JA) = (ZTT1(JL,JA)+ZTT2(JL,JA))*0.5
 2234 CONTINUE
 2235 CONTINUE
C
      DO 2236 JL = 1, KDLON
      ZWW=PDBDT(JL,1,JKJ)*ZTT(JL,1)          *ZTT(JL,10)
     S   +PDBDT(JL,2,JKJ)*ZTT(JL,2)*ZTT(JL,7)*ZTT(JL,11)
     S   +PDBDT(JL,3,JKJ)*ZTT(JL,4)*ZTT(JL,8)*ZTT(JL,12)
     S   +PDBDT(JL,4,JKJ)*ZTT(JL,5)*ZTT(JL,9)*ZTT(JL,13)
     S   +PDBDT(JL,5,JKJ)*ZTT(JL,3)          *ZTT(JL,14)
     S   +PDBDT(JL,6,JKJ)*ZTT(JL,6)          *ZTT(JL,15)
      ZGLAYD(JL)=ZWW
      ZDZXDG=ZGLAYD(JL)
      PDISD(JL,JK)=PDISD(JL,JK)+ZDZXDG
      PCNTRB(JL,JK,IKJP1)=ZDZXDG
 2236 CONTINUE
C
C
 224  CONTINUE
 225  CONTINUE
C
C
C*         2.2.4   DISTANT AND BELOW LAYERS
C                  ------------------------
C
 2240 CONTINUE
C
C
C
C*         2.2.5   FIRST LOWER LEVEL
C                  -----------------
C
 2250 CONTINUE
C
      DO 228 JK=3,KFLEV+1
      IKN=(JK-1)*NG1P1+1
      IKM1=JK-1
      IKJ=JK-2
      IKU1= IKJ  *NG1P1+1
C
C
      CALL LWTTM(PGA(1,1,1,IKJ),PGB(1,1,1,IKJ)
     2          , PABCU(1,1,IKU1),PABCU(1,1,IKN),ZTT1)
C
C
C
C*         2.2.6   DOWN BELOW
C                  ----------
C
 2260 CONTINUE
C
      ITT=1
      DO 227 JLK=1,IKJ
      IF(ITT.EQ.1) THEN
         ITT=2
      ELSE
         ITT=1
      ENDIF
      IJKL=IKM1-JLK
      IKU2=(IJKL-1)*NG1P1+1
C
C
      IF(ITT.EQ.1) THEN
         CALL LWTTM(PGA(1,1,1,IJKL),PGB(1,1,1,IJKL)
     2             , PABCU(1,1,IKU2),PABCU(1,1,IKN),ZTT1)
      ELSE
         CALL LWTTM(PGA(1,1,1,IJKL),PGB(1,1,1,IJKL)
     2             , PABCU(1,1,IKU2),PABCU(1,1,IKN),ZTT2)
      ENDIF
C
      DO 2265 JA = 1, KTRAER
      DO 2264 JL = 1, KDLON
      ZTT(JL,JA) = (ZTT1(JL,JA)+ZTT2(JL,JA))*0.5
 2264 CONTINUE
 2265 CONTINUE
C
      DO 2266 JL = 1, KDLON
      ZWW=PDBDT(JL,1,IJKL)*ZTT(JL,1)          *ZTT(JL,10)
     S   +PDBDT(JL,2,IJKL)*ZTT(JL,2)*ZTT(JL,7)*ZTT(JL,11)
     S   +PDBDT(JL,3,IJKL)*ZTT(JL,4)*ZTT(JL,8)*ZTT(JL,12)
     S   +PDBDT(JL,4,IJKL)*ZTT(JL,5)*ZTT(JL,9)*ZTT(JL,13)
     S   +PDBDT(JL,5,IJKL)*ZTT(JL,3)          *ZTT(JL,14)
     S   +PDBDT(JL,6,IJKL)*ZTT(JL,6)          *ZTT(JL,15)
      ZGLAYU(JL)=ZWW
      ZDZXMG=ZGLAYU(JL)
      PDISU(JL,JK)=PDISU(JL,JK)+ZDZXMG
      PCNTRB(JL,JK,IJKL)=ZDZXMG
 2266 CONTINUE
C
C
 227  CONTINUE
 228  CONTINUE
C
      RETURN
      END
      SUBROUTINE LWVN(KUAER,KTRAER
     R  , PABCU,PDBSL,PGA,PGB
     S  , PADJD,PADJU,PCNTRB,PDBDT)
      IMPLICIT none
      include "dimensions.h"
      include "dimphy.h"
      include "raddim.h"
      include "raddimlw.h"
C
C-----------------------------------------------------------------------
C     PURPOSE.
C     --------
C           CARRIES OUT THE VERTICAL INTEGRATION ON NEARBY LAYERS
C           TO GIVE LONGWAVE FLUXES OR RADIANCES
C
C     METHOD.
C     -------
C
C          1. PERFORMS THE VERTICAL INTEGRATION CORRESPONDING TO THE
C     CONTRIBUTIONS OF THE ADJACENT LAYERS USING A GAUSSIAN QUADRATURE
C
C     REFERENCE.
C     ----------
C
C        SEE RADIATION'S PART OF THE MODEL'S DOCUMENTATION AND
C        ECMWF RESEARCH DEPARTMENT DOCUMENTATION OF THE IFS
C
C     AUTHOR.
C     -------
C        JEAN-JACQUES MORCRETTE  *ECMWF*
C
C     MODIFICATIONS.
C     --------------
C        ORIGINAL : 89-07-14
C-----------------------------------------------------------------------
C
C* ARGUMENTS:
C
      INTEGER KUAER,KTRAER
C
      REAL*8 PABCU(KDLON,NUA,3*KFLEV+1) ! ABSORBER AMOUNTS
      REAL*8 PDBSL(KDLON,Ninter,KFLEV*2) ! SUB-LAYER PLANCK FUNCTION GRADIENT
      REAL*8 PGA(KDLON,8,2,KFLEV) ! PADE APPROXIMANTS
      REAL*8 PGB(KDLON,8,2,KFLEV) ! PADE APPROXIMANTS
C
      REAL*8 PADJD(KDLON,KFLEV+1) ! CONTRIBUTION OF ADJACENT LAYERS
      REAL*8 PADJU(KDLON,KFLEV+1) ! CONTRIBUTION OF ADJACENT LAYERS
      REAL*8 PCNTRB(KDLON,KFLEV+1,KFLEV+1) ! CLEAR-SKY ENERGY EXCHANGE MATRIX
      REAL*8 PDBDT(KDLON,Ninter,KFLEV) !  LAYER PLANCK FUNCTION GRADIENT
C
C* LOCAL ARRAYS:
C
      REAL*8 ZGLAYD(KDLON)
      REAL*8 ZGLAYU(KDLON)
      REAL*8 ZTT(KDLON,NTRA)
      REAL*8 ZTT1(KDLON,NTRA)
      REAL*8 ZTT2(KDLON,NTRA)
      REAL*8 ZUU(KDLON,NUA)
C
      INTEGER jk, jl, ja, im12, ind, inu, ixu, jg
      INTEGER ixd, ibs, idd, imu, jk1, jk2, jnu
      REAL*8 zwtr
c
C* Data Block:
c
      REAL*8 WG1(2)
      SAVE WG1
      DATA (WG1(jk),jk=1,2) /1.0, 1.0/
C-----------------------------------------------------------------------
C
C*         1.    INITIALIZATION
C                --------------
C
 100  CONTINUE
C
C*         1.1     INITIALIZE LAYER CONTRIBUTIONS
C                  ------------------------------
C
 110  CONTINUE
C
      DO 112 JK = 1 , KFLEV+1
      DO 111 JL = 1, KDLON
      PADJD(JL,JK) = 0.
      PADJU(JL,JK) = 0.
 111  CONTINUE
 112  CONTINUE
C
C*         1.2     INITIALIZE TRANSMISSION FUNCTIONS
C                  ---------------------------------
C
 120  CONTINUE
C
      DO 122 JA = 1 , NTRA
      DO 121 JL = 1, KDLON
      ZTT (JL,JA) = 1.0
      ZTT1(JL,JA) = 1.0
      ZTT2(JL,JA) = 1.0
 121  CONTINUE
 122  CONTINUE
C
      DO 124 JA = 1 , NUA
      DO 123 JL = 1, KDLON
      ZUU(JL,JA) = 0.
 123  CONTINUE
 124  CONTINUE
C
C     ------------------------------------------------------------------
C
C*         2.      VERTICAL INTEGRATION
C                  --------------------
C
 200  CONTINUE
C
C
C*         2.1     CONTRIBUTION FROM ADJACENT LAYERS
C                  ---------------------------------
C
 210  CONTINUE
C
      DO 215 JK = 1 , KFLEV
C
C*         2.1.1   DOWNWARD LAYERS
C                  ---------------
C
 2110 CONTINUE
C
      IM12 = 2 * (JK - 1)
      IND = (JK - 1) * NG1P1 + 1
      IXD = IND
      INU = JK * NG1P1 + 1
      IXU = IND
C
      DO 2111 JL = 1, KDLON
      ZGLAYD(JL) = 0.
      ZGLAYU(JL) = 0.
 2111 CONTINUE
C
      DO 213 JG = 1 , NG1
      IBS = IM12 + JG
      IDD = IXD + JG
      DO 2113 JA = 1 , KUAER
      DO 2112 JL = 1, KDLON
      ZUU(JL,JA) = PABCU(JL,JA,IND) - PABCU(JL,JA,IDD)
 2112 CONTINUE
 2113 CONTINUE
C
C
      CALL LWTT(PGA(1,1,1,JK), PGB(1,1,1,JK), ZUU, ZTT)
C
      DO 2114 JL = 1, KDLON
      ZWTR=PDBSL(JL,1,IBS)*ZTT(JL,1)          *ZTT(JL,10)
     S    +PDBSL(JL,2,IBS)*ZTT(JL,2)*ZTT(JL,7)*ZTT(JL,11)
     S    +PDBSL(JL,3,IBS)*ZTT(JL,4)*ZTT(JL,8)*ZTT(JL,12)
     S    +PDBSL(JL,4,IBS)*ZTT(JL,5)*ZTT(JL,9)*ZTT(JL,13)
     S    +PDBSL(JL,5,IBS)*ZTT(JL,3)          *ZTT(JL,14)
     S    +PDBSL(JL,6,IBS)*ZTT(JL,6)          *ZTT(JL,15)
      ZGLAYD(JL)=ZGLAYD(JL)+ZWTR*WG1(JG)
 2114 CONTINUE
C
C*         2.1.2   DOWNWARD LAYERS
C                  ---------------
C
 2120 CONTINUE
C
      IMU = IXU + JG
      DO 2122 JA = 1 , KUAER
      DO 2121 JL = 1, KDLON
      ZUU(JL,JA) = PABCU(JL,JA,IMU) - PABCU(JL,JA,INU)
 2121 CONTINUE
 2122 CONTINUE
C
C
      CALL LWTT(PGA(1,1,1,JK), PGB(1,1,1,JK), ZUU, ZTT)
C
      DO 2123 JL = 1, KDLON
      ZWTR=PDBSL(JL,1,IBS)*ZTT(JL,1)          *ZTT(JL,10)
     S    +PDBSL(JL,2,IBS)*ZTT(JL,2)*ZTT(JL,7)*ZTT(JL,11)
     S    +PDBSL(JL,3,IBS)*ZTT(JL,4)*ZTT(JL,8)*ZTT(JL,12)
     S    +PDBSL(JL,4,IBS)*ZTT(JL,5)*ZTT(JL,9)*ZTT(JL,13)
     S    +PDBSL(JL,5,IBS)*ZTT(JL,3)          *ZTT(JL,14)
     S    +PDBSL(JL,6,IBS)*ZTT(JL,6)          *ZTT(JL,15)
      ZGLAYU(JL)=ZGLAYU(JL)+ZWTR*WG1(JG)
 2123 CONTINUE
C
 213  CONTINUE
C
      DO 214 JL = 1, KDLON
      PADJD(JL,JK) = ZGLAYD(JL)
      PCNTRB(JL,JK,JK+1) = ZGLAYD(JL)
      PADJU(JL,JK+1) = ZGLAYU(JL)
      PCNTRB(JL,JK+1,JK) = ZGLAYU(JL)
      PCNTRB(JL,JK  ,JK) = 0.0
 214  CONTINUE
C
 215  CONTINUE
C
      DO 218 JK = 1 , KFLEV
      JK2 = 2 * JK
      JK1 = JK2 - 1
      DO 217 JNU = 1 , Ninter
      DO 216 JL = 1, KDLON
      PDBDT(JL,JNU,JK) = PDBSL(JL,JNU,JK1) + PDBSL(JL,JNU,JK2)
 216  CONTINUE
 217  CONTINUE
 218  CONTINUE
C
      RETURN
C
      END
      SUBROUTINE LWTT(PGA,PGB,PUU, PTT)
      IMPLICIT none
      include "dimensions.h"
      include "dimphy.h"
      include "raddim.h"
      include "raddimlw.h"
C
C-----------------------------------------------------------------------
C     PURPOSE.
C     --------
C           THIS ROUTINE COMPUTES THE TRANSMISSION FUNCTIONS FOR ALL THE
C     ABSORBERS (H2O, UNIFORMLY MIXED GASES, AND O3) IN ALL SIX SPECTRAL
C     INTERVALS.
C
C     METHOD.
C     -------
C
C          1. TRANSMISSION FUNCTION BY H2O AND UNIFORMLY MIXED GASES ARE
C     COMPUTED USING PADE APPROXIMANTS AND HORNER'S ALGORITHM.
C          2. TRANSMISSION BY O3 IS EVALUATED WITH MALKMUS'S BAND MODEL.
C          3. TRANSMISSION BY H2O CONTINUUM AND AEROSOLS FOLLOW AN
C     A SIMPLE EXPONENTIAL DECREASE WITH ABSORBER AMOUNT.
C
C     REFERENCE.
C     ----------
C
C        SEE RADIATION'S PART OF THE MODEL'S DOCUMENTATION AND
C        ECMWF RESEARCH DEPARTMENT DOCUMENTATION OF THE IFS
C
C     AUTHOR.
C     -------
C        JEAN-JACQUES MORCRETTE  *ECMWF*
C
C     MODIFICATIONS.
C     --------------
C        ORIGINAL : 88-12-15
C
C-----------------------------------------------------------------------
      REAL*8 O1H, O2H
      PARAMETER (O1H=2230.)
      PARAMETER (O2H=100.)
      REAL*8 RPIALF0
      PARAMETER (RPIALF0=2.0)
C
C* ARGUMENTS:
C
      REAL*8 PUU(KDLON,NUA)
      REAL*8 PTT(KDLON,NTRA)
      REAL*8 PGA(KDLON,8,2)
      REAL*8 PGB(KDLON,8,2)
C
C* LOCAL VARIABLES:
C
      REAL*8 zz, zxd, zxn
      REAL*8 zpu, zpu10, zpu11, zpu12, zpu13
      REAL*8 zeu, zeu10, zeu11, zeu12, zeu13
      REAL*8 zx, zy, zsq1, zsq2, zvxy, zuxy
      REAL*8 zaercn, zto1, zto2, zxch4, zych4, zxn2o, zyn2o
      REAL*8 zsqn21, zodn21, zsqh42, zodh42
      REAL*8 zsqh41, zodh41, zsqn22, zodn22, zttf11, zttf12
      REAL*8 zuu11, zuu12, za11, za12
      INTEGER jl, ja
C     ------------------------------------------------------------------
C
C*         1.     HORNER'S ALGORITHM FOR H2O AND CO2 TRANSMISSION
C                 -----------------------------------------------
C
 100  CONTINUE
C
C
      DO 130 JA = 1 , 8
      DO 120 JL = 1, KDLON
      ZZ      =SQRT(PUU(JL,JA))
c     ZXD(JL,1)=PGB( JL, 1,1) + ZZ(JL, 1)*(PGB( JL, 1,2) + ZZ(JL, 1))
c     ZXN(JL,1)=PGA( JL, 1,1) + ZZ(JL, 1)*(PGA( JL, 1,2) )
c     PTT(JL,1)=ZXN(JL,1)/ZXD(JL,1)
      ZXD      =PGB( JL,JA,1) + ZZ       *(PGB( JL,JA,2) + ZZ       )
      ZXN      =PGA( JL,JA,1) + ZZ       *(PGA( JL,JA,2) )
      PTT(JL,JA)=ZXN      /ZXD
  120 CONTINUE
  130 CONTINUE
C
C     ------------------------------------------------------------------
C
C*         2.     CONTINUUM, OZONE AND AEROSOL TRANSMISSION FUNCTIONS
C                 ---------------------------------------------------
C
 200  CONTINUE
C
      DO 201 JL = 1, KDLON
      PTT(JL, 9) = PTT(JL, 8)
C
C-  CONTINUUM ABSORPTION: E- AND P-TYPE
C
      ZPU   = 0.002 * PUU(JL,10)
      ZPU10 = 112. * ZPU
      ZPU11 = 6.25 * ZPU
      ZPU12 = 5.00 * ZPU
      ZPU13 = 80.0 * ZPU
      ZEU   =  PUU(JL,11)
      ZEU10 =  12. * ZEU
      ZEU11 = 6.25 * ZEU
      ZEU12 = 5.00 * ZEU
      ZEU13 = 80.0 * ZEU
C
C-  OZONE ABSORPTION
C
      ZX = PUU(JL,12)
      ZY = PUU(JL,13)
      ZUXY = 4. * ZX * ZX / (RPIALF0 * ZY)
      ZSQ1 = SQRT(1. + O1H * ZUXY ) - 1.
      ZSQ2 = SQRT(1. + O2H * ZUXY ) - 1.
      ZVXY = RPIALF0 * ZY / (2. * ZX)
      ZAERCN = PUU(JL,17) + ZEU12 + ZPU12
      ZTO1 = EXP( - ZVXY * ZSQ1 - ZAERCN )
      ZTO2 = EXP( - ZVXY * ZSQ2 - ZAERCN )
C
C-- TRACE GASES (CH4, N2O, CFC-11, CFC-12)
C
C* CH4 IN INTERVAL 800-970 + 1110-1250 CM-1
C
c     NEXOTIC=1
c     IF (NEXOTIC.EQ.1) THEN
      ZXCH4 = PUU(JL,19)
      ZYCH4 = PUU(JL,20)
      ZUXY = 4. * ZXCH4*ZXCH4/(0.103*ZYCH4)
      ZSQH41 = SQRT(1. + 33.7 * ZUXY) - 1.
      ZVXY = 0.103 * ZYCH4 / (2. * ZXCH4)
      ZODH41 = ZVXY * ZSQH41
C
C* N2O IN INTERVAL 800-970 + 1110-1250 CM-1
C
      ZXN2O = PUU(JL,21)
      ZYN2O = PUU(JL,22)
      ZUXY = 4. * ZXN2O*ZXN2O/(0.416*ZYN2O)
      ZSQN21 = SQRT(1. + 21.3 * ZUXY) - 1.
      ZVXY = 0.416 * ZYN2O / (2. * ZXN2O)
      ZODN21 = ZVXY * ZSQN21
C
C* CH4 IN INTERVAL 1250-1450 + 1880-2820 CM-1
C
      ZUXY = 4. * ZXCH4*ZXCH4/(0.113*ZYCH4)
      ZSQH42 = SQRT(1. + 400. * ZUXY) - 1.
      ZVXY = 0.113 * ZYCH4 / (2. * ZXCH4)
      ZODH42 = ZVXY * ZSQH42
C
C* N2O IN INTERVAL 1250-1450 + 1880-2820 CM-1
C
      ZUXY = 4. * ZXN2O*ZXN2O/(0.197*ZYN2O)
      ZSQN22 = SQRT(1. + 2000. * ZUXY) - 1.
      ZVXY = 0.197 * ZYN2O / (2. * ZXN2O)
      ZODN22 = ZVXY * ZSQN22
C
C* CFC-11 IN INTERVAL 800-970 + 1110-1250 CM-1
C
      ZA11 = 2. * PUU(JL,23) * 4.404E+05
      ZTTF11 = 1. - ZA11 * 0.003225
C
C* CFC-12 IN INTERVAL 800-970 + 1110-1250 CM-1
C
      ZA12 = 2. * PUU(JL,24) * 6.7435E+05
      ZTTF12 = 1. - ZA12 * 0.003225
C
      ZUU11 = - PUU(JL,15) - ZEU10 - ZPU10
      ZUU12 = - PUU(JL,16) - ZEU11 - ZPU11 - ZODH41 - ZODN21
      PTT(JL,10) = EXP( - PUU(JL,14) )
      PTT(JL,11) = EXP( ZUU11 )
      PTT(JL,12) = EXP( ZUU12 ) * ZTTF11 * ZTTF12
      PTT(JL,13) = 0.7554 * ZTO1 + 0.2446 * ZTO2
      PTT(JL,14) = PTT(JL,10) * EXP( - ZEU13 - ZPU13 )
      PTT(JL,15) = EXP ( - PUU(JL,14) - ZODH42 - ZODN22 )
 201  CONTINUE
C
      RETURN
      END
      SUBROUTINE LWTTM(PGA,PGB,PUU1,PUU2, PTT)
      IMPLICIT none
      include "dimensions.h"
      include "dimphy.h"
      include "raddim.h"
      include "raddimlw.h"
C
C     ------------------------------------------------------------------
C     PURPOSE.
C     --------
C           THIS ROUTINE COMPUTES THE TRANSMISSION FUNCTIONS FOR ALL THE
C     ABSORBERS (H2O, UNIFORMLY MIXED GASES, AND O3) IN ALL SIX SPECTRAL
C     INTERVALS.
C
C     METHOD.
C     -------
C
C          1. TRANSMISSION FUNCTION BY H2O AND UNIFORMLY MIXED GASES ARE
C     COMPUTED USING PADE APPROXIMANTS AND HORNER'S ALGORITHM.
C          2. TRANSMISSION BY O3 IS EVALUATED WITH MALKMUS'S BAND MODEL.
C          3. TRANSMISSION BY H2O CONTINUUM AND AEROSOLS FOLLOW AN
C     A SIMPLE EXPONENTIAL DECREASE WITH ABSORBER AMOUNT.
C
C     REFERENCE.
C     ----------
C
C        SEE RADIATION'S PART OF THE MODEL'S DOCUMENTATION AND
C        ECMWF RESEARCH DEPARTMENT DOCUMENTATION OF THE IFS
C
C     AUTHOR.
C     -------
C        JEAN-JACQUES MORCRETTE  *ECMWF*
C
C     MODIFICATIONS.
C     --------------
C        ORIGINAL : 88-12-15
C
C-----------------------------------------------------------------------
      REAL*8 O1H, O2H
      PARAMETER (O1H=2230.)
      PARAMETER (O2H=100.)
      REAL*8 RPIALF0
      PARAMETER (RPIALF0=2.0)
C
C* ARGUMENTS:
C
      REAL*8 PGA(KDLON,8,2) ! PADE APPROXIMANTS
      REAL*8 PGB(KDLON,8,2) ! PADE APPROXIMANTS
      REAL*8 PUU1(KDLON,NUA) ! ABSORBER AMOUNTS FROM TOP TO LEVEL 1
      REAL*8 PUU2(KDLON,NUA) ! ABSORBER AMOUNTS FROM TOP TO LEVEL 2
      REAL*8 PTT(KDLON,NTRA) ! TRANSMISSION FUNCTIONS
C
C* LOCAL VARIABLES:
C
      INTEGER ja, jl
      REAL*8 zz, zxd, zxn
      REAL*8 zpu, zpu10, zpu11, zpu12, zpu13
      REAL*8 zeu, zeu10, zeu11, zeu12, zeu13
      REAL*8 zx, zy, zuxy, zsq1, zsq2, zvxy, zaercn, zto1, zto2
      REAL*8 zxch4, zych4, zsqh41, zodh41
      REAL*8 zxn2o, zyn2o, zsqn21, zodn21, zsqh42, zodh42
      REAL*8 zsqn22, zodn22, za11, zttf11, za12, zttf12
      REAL*8 zuu11, zuu12
C     ------------------------------------------------------------------
C
C*         1.     HORNER'S ALGORITHM FOR H2O AND CO2 TRANSMISSION
C                 -----------------------------------------------
C
 100  CONTINUE
C
C
      DO 130 JA = 1 , 8
      DO 120 JL = 1, KDLON
      ZZ      =SQRT(PUU1(JL,JA) - PUU2(JL,JA))
      ZXD      =PGB( JL,JA,1) + ZZ       *(PGB( JL,JA,2) + ZZ       )
      ZXN      =PGA( JL,JA,1) + ZZ       *(PGA( JL,JA,2) )
      PTT(JL,JA)=ZXN      /ZXD
  120 CONTINUE
  130 CONTINUE
C
C     ------------------------------------------------------------------
C
C*         2.     CONTINUUM, OZONE AND AEROSOL TRANSMISSION FUNCTIONS
C                 ---------------------------------------------------
C
 200  CONTINUE
C
      DO 201 JL = 1, KDLON
      PTT(JL, 9) = PTT(JL, 8)
C
C-  CONTINUUM ABSORPTION: E- AND P-TYPE
C
      ZPU   = 0.002 * (PUU1(JL,10) - PUU2(JL,10))
      ZPU10 = 112. * ZPU
      ZPU11 = 6.25 * ZPU
      ZPU12 = 5.00 * ZPU
      ZPU13 = 80.0 * ZPU
      ZEU   = (PUU1(JL,11) - PUU2(JL,11))
      ZEU10 =  12. * ZEU
      ZEU11 = 6.25 * ZEU
      ZEU12 = 5.00 * ZEU
      ZEU13 = 80.0 * ZEU
C
C-  OZONE ABSORPTION
C
      ZX = (PUU1(JL,12) - PUU2(JL,12))
      ZY = (PUU1(JL,13) - PUU2(JL,13))
      ZUXY = 4. * ZX * ZX / (RPIALF0 * ZY)
      ZSQ1 = SQRT(1. + O1H * ZUXY ) - 1.
      ZSQ2 = SQRT(1. + O2H * ZUXY ) - 1.
      ZVXY = RPIALF0 * ZY / (2. * ZX)
      ZAERCN = (PUU1(JL,17) -PUU2(JL,17)) + ZEU12 + ZPU12
      ZTO1 = EXP( - ZVXY * ZSQ1 - ZAERCN )
      ZTO2 = EXP( - ZVXY * ZSQ2 - ZAERCN )
C
C-- TRACE GASES (CH4, N2O, CFC-11, CFC-12)
C
C* CH4 IN INTERVAL 800-970 + 1110-1250 CM-1
C
      ZXCH4 = (PUU1(JL,19) - PUU2(JL,19))
      ZYCH4 = (PUU1(JL,20) - PUU2(JL,20))
      ZUXY = 4. * ZXCH4*ZXCH4/(0.103*ZYCH4)
      ZSQH41 = SQRT(1. + 33.7 * ZUXY) - 1.
      ZVXY = 0.103 * ZYCH4 / (2. * ZXCH4)
      ZODH41 = ZVXY * ZSQH41
C
C* N2O IN INTERVAL 800-970 + 1110-1250 CM-1
C
      ZXN2O = (PUU1(JL,21) - PUU2(JL,21))
      ZYN2O = (PUU1(JL,22) - PUU2(JL,22))
      ZUXY = 4. * ZXN2O*ZXN2O/(0.416*ZYN2O)
      ZSQN21 = SQRT(1. + 21.3 * ZUXY) - 1.
      ZVXY = 0.416 * ZYN2O / (2. * ZXN2O)
      ZODN21 = ZVXY * ZSQN21
C
C* CH4 IN INTERVAL 1250-1450 + 1880-2820 CM-1
C
      ZUXY = 4. * ZXCH4*ZXCH4/(0.113*ZYCH4)
      ZSQH42 = SQRT(1. + 400. * ZUXY) - 1.
      ZVXY = 0.113 * ZYCH4 / (2. * ZXCH4)
      ZODH42 = ZVXY * ZSQH42
C
C* N2O IN INTERVAL 1250-1450 + 1880-2820 CM-1
C
      ZUXY = 4. * ZXN2O*ZXN2O/(0.197*ZYN2O)
      ZSQN22 = SQRT(1. + 2000. * ZUXY) - 1.
      ZVXY = 0.197 * ZYN2O / (2. * ZXN2O)
      ZODN22 = ZVXY * ZSQN22
C
C* CFC-11 IN INTERVAL 800-970 + 1110-1250 CM-1
C
      ZA11 = (PUU1(JL,23) - PUU2(JL,23)) * 4.404E+05
      ZTTF11 = 1. - ZA11 * 0.003225
C
C* CFC-12 IN INTERVAL 800-970 + 1110-1250 CM-1
C
      ZA12 = (PUU1(JL,24) - PUU2(JL,24)) * 6.7435E+05
      ZTTF12 = 1. - ZA12 * 0.003225
C
      ZUU11 = - (PUU1(JL,15) - PUU2(JL,15)) - ZEU10 - ZPU10
      ZUU12 = - (PUU1(JL,16) - PUU2(JL,16)) - ZEU11 - ZPU11 -
     S         ZODH41 - ZODN21
      PTT(JL,10) = EXP( - (PUU1(JL,14)- PUU2(JL,14)) )
      PTT(JL,11) = EXP( ZUU11 )
      PTT(JL,12) = EXP( ZUU12 ) * ZTTF11 * ZTTF12
      PTT(JL,13) = 0.7554 * ZTO1 + 0.2446 * ZTO2
      PTT(JL,14) = PTT(JL,10) * EXP( - ZEU13 - ZPU13 )
      PTT(JL,15) = EXP ( - (PUU1(JL,14) - PUU2(JL,14)) - ZODH42-ZODN22 )
 201  CONTINUE
C
      RETURN
      END
c==============================================================
      SUBROUTINE suphec
C
      include "YOMCST.h"
      include "YOETHF.h"
C      -----------------------------------------------------------------
C
C*       1.    DEFINE FUNDAMENTAL CONSTANTS.
C              -----------------------------
C
c sb      WRITE(UNIT=6,FMT='(''0*** Constants of the ICM   ***'')')
      RPI=2.*ASIN(1.)
      RCLUM=299792458.
      RHPLA=6.6260755E-34
      RKBOL=1.380658E-23
      RNAVO=6.0221367E+23
c sb      WRITE(UNIT=6,FMT='('' *** Fundamental constants ***'')')
c sb      WRITE(UNIT=6,FMT='(''           PI = '',E13.7,'' -'')')RPI
c sb      WRITE(UNIT=6,FMT='(''            c = '',E13.7,''m s-1'')')
c sb     S RCLUM
c sb      WRITE(UNIT=6,FMT='(''            h = '',E13.7,''J s'')')
c sb     S RHPLA
c sb      WRITE(UNIT=6,FMT='(''            K = '',E13.7,''J K-1'')')
c sb     S RKBOL
c sb      WRITE(UNIT=6,FMT='(''            N = '',E13.7,''mol-1'')')
c sb     S RNAVO
C
C     ----------------------------------------------------------------
C
C*       2.    DEFINE ASTRONOMICAL CONSTANTS.
C              ------------------------------
C
      RDAY=86400.
      REA=149597870000.
      REPSM=0.409093
C
      RSIYEA=365.25*RDAY*2.*RPI/6.283076
      RSIDAY=RDAY/(1.+RDAY/RSIYEA)
      ROMEGA=2.*RPI/RSIDAY
c
c exp1      R_ecc = 0.05
c exp1      R_peri = 102.04
c exp1      R_incl = 22.5
c exp1      print*, 'Parametres orbitaux modifies'
c ref      R_ecc = 0.016724
c ref      R_peri = 102.04
c ref      R_incl = 23.5
      R_ecc = 0.016724
      R_peri = 102.04
      R_incl = 23.5
c
c sb      WRITE(UNIT=6,FMT='('' *** Astronomical constants ***'')')
c sb      WRITE(UNIT=6,FMT='(''          day = '',E13.7,'' s'')')RDAY
c sb      WRITE(UNIT=6,FMT='('' half g. axis = '',E13.7,'' m'')')REA
c sb      WRITE(UNIT=6,FMT='('' mean anomaly = '',E13.7,'' -'')')REPSM
c sb      WRITE(UNIT=6,FMT='('' sideral year = '',E13.7,'' s'')')RSIYEA
c sb      WRITE(UNIT=6,FMT='(''  sideral day = '',E13.7,'' s'')')RSIDAY
c sb      WRITE(UNIT=6,FMT='(''        omega = '',E13.7,'' s-1'')')
c sb     S                  ROMEGA
c sb      write(unit=6,fmt='('' excentricite = '',e13.7,''-'')')R_ecc
c sb      write(unit=6,fmt='(''     equinoxe = '',e13.7,''-'')')R_peri
c sb      write(unit=6,fmt='(''  inclinaison = '',e13.7,''-'')')R_incl
C
C     ------------------------------------------------------------------
C
C*       3.    DEFINE GEOIDE.
C              --------------
C
      RG=9.80665
      RA=6371229.
      R1SA=SNGL(1.D0/DBLE(RA))
c sb      WRITE(UNIT=6,FMT='('' ***         Geoide         ***'')')
c sb      WRITE(UNIT=6,FMT='(''      Gravity = '',E13.7,'' m s-2'')')
c sb     S      RG
c sb      WRITE(UNIT=6,FMT='('' Earth radius = '',E13.7,'' m'')')RA
c sb      WRITE(UNIT=6,FMT='('' Inverse E.R. = '',E13.7,'' m'')')R1SA
C
C     -----------------------------------------------------------------
C
C*       4.    DEFINE RADIATION CONSTANTS.
C              ---------------------------
C
c z.x.li      RSIGMA=2. * RPI**5 * RKBOL**4 /(15.* RCLUM**2 * RHPLA**3)
      rsigma = 2.*rpi**5 * (rkbol/rhpla)**3 * rkbol/rclum/rclum/15.
      RI0=1370.
c sb      WRITE(UNIT=6,FMT='('' ***        Radiation       ***'')')
c sb      WRITE(UNIT=6,FMT='('' Stefan-Bol.  = '',E13.7,'' W m-2 K-4''
c sb     S )')  RSIGMA
c sb      WRITE(UNIT=6,FMT='('' Solar const. = '',E13.7,'' W m-2'')')
c sb     S      RI0
C
C     -----------------------------------------------------------------
C
C*       5.    DEFINE THERMODYNAMIC CONSTANTS, GAS PHASE.
C              ------------------------------------------
C
      R=RNAVO*RKBOL
      RMD=28.9644
      RMV=18.0153
      RD=1000.*R/RMD
      RV=1000.*R/RMV
      RCPD=3.5*RD
      RCVD=RCPD-RD
      RCPV=4. *RV
      RCVV=RCPV-RV
      RKAPPA=RD/RCPD
      RETV=RV/RD-1.
c sb      WRITE(UNIT=6,FMT='('' *** Thermodynamic, gas     ***'')')
c sb      WRITE(UNIT=6,FMT='('' Perfect gas  = '',e13.7)') R
c sb      WRITE(UNIT=6,FMT='('' Dry air mass = '',e13.7)') RMD
c sb      WRITE(UNIT=6,FMT='('' Vapour  mass = '',e13.7)') RMV
c sb      WRITE(UNIT=6,FMT='('' Dry air cst. = '',e13.7)') RD
c sb      WRITE(UNIT=6,FMT='('' Vapour  cst. = '',e13.7)') RV
c sb      WRITE(UNIT=6,FMT='(''         Cpd  = '',e13.7)') RCPD
c sb      WRITE(UNIT=6,FMT='(''         Cvd  = '',e13.7)') RCVD
c sb      WRITE(UNIT=6,FMT='(''         Cpv  = '',e13.7)') RCPV
c sb      WRITE(UNIT=6,FMT='(''         Cvv  = '',e13.7)') RCVV
c sb      WRITE(UNIT=6,FMT='(''      Rd/Cpd  = '',e13.7)') RKAPPA
c sb      WRITE(UNIT=6,FMT='(''     Rv/Rd-1  = '',e13.7)') RETV
C
C     ----------------------------------------------------------------
C
C*       6.    DEFINE THERMODYNAMIC CONSTANTS, LIQUID PHASE.
C              ---------------------------------------------
C
      RCW=4218.
c sb      WRITE(UNIT=6,FMT='('' *** Thermodynamic, liquid  ***'')')
c sb      WRITE(UNIT=6,FMT='(''         Cw   = '',E13.7)') RCW
C
C     ----------------------------------------------------------------
C
C*       7.    DEFINE THERMODYNAMIC CONSTANTS, SOLID PHASE.
C              --------------------------------------------
C
      RCS=2106.
c sb      WRITE(UNIT=6,FMT='('' *** thermodynamic, solid   ***'')')
c sb      WRITE(UNIT=6,FMT='(''         Cs   = '',E13.7)') RCS
C
C     ----------------------------------------------------------------
C
C*       8.    DEFINE THERMODYNAMIC CONSTANTS, TRANSITION OF PHASE.
C              ----------------------------------------------------
C
      RTT=273.16
      RLVTT=2.5008E+6
      RLSTT=2.8345E+6
      RLMLT=RLSTT-RLVTT
      RATM=100000.
c sb      WRITE(UNIT=6,FMT='('' *** Thermodynamic, trans.  ***'')')
c sb      WRITE(UNIT=6,FMT='('' Fusion point  = '',E13.7)') RTT
c sb      WRITE(UNIT=6,FMT='(''        RLvTt  = '',E13.7)') RLVTT
c sb      WRITE(UNIT=6,FMT='(''        RLsTt  = '',E13.7)') RLSTT
c sb      WRITE(UNIT=6,FMT='(''        RLMlt  = '',E13.7)') RLMLT
c sb      WRITE(UNIT=6,FMT='('' Normal press. = '',E13.7)') RATM
c sb      WRITE(UNIT=6,FMT='('' Latent heat :  '')')
C
C     ----------------------------------------------------------------
C
C*       9.    SATURATED VAPOUR PRESSURE.
C              --------------------------
C
      RESTT=611.14
      RGAMW=(RCW-RCPV)/RV
      RBETW=RLVTT/RV+RGAMW*RTT
      RALPW=LOG(RESTT)+RBETW/RTT+RGAMW*LOG(RTT)
      RGAMS=(RCS-RCPV)/RV
      RBETS=RLSTT/RV+RGAMS*RTT
      RALPS=LOG(RESTT)+RBETS/RTT+RGAMS*LOG(RTT)
      RGAMD=RGAMS-RGAMW
      RBETD=RBETS-RBETW
      RALPD=RALPS-RALPW
C
C     ------------------------------------------------------------------
c
c calculer les constantes pour les fonctions thermodynamiques
c
      RVTMP2=RCPV/RCPD-1.
      RHOH2O=RATM/100.
c sb      WRITE(UNIT=6,FMT='(''        RHOL  = '',E13.7)') RHOH2O ! ajout sb
      R2ES=RESTT*RD/RV
      R3LES=17.269
      R3IES=21.875
      R4LES=35.86
      R4IES=7.66
      R5LES=R3LES*(RTT-R4LES)
      R5IES=R3IES*(RTT-R4IES)
C
      RETURN
      END
      SUBROUTINE alboc(rjour,rlat,albedo)
      IMPLICIT none
c======================================================================
c Auteur(s): Z.X. Li (LMD/CNRS) (adaptation du GCM du LMD)
c Date: le 16 mars 1995
c Objet: Calculer l'albedo sur l'ocean
c Methode: Integrer numeriquement l'albedo pendant une journee
c
c Arguments;
c rjour (in,R)  : jour dans l'annee (a compter du 1 janvier)
c rlat (in,R)   : latitude en degre
c albedo (out,R): albedo obtenu (de 0 a 1)
c
c Nov 1999: modification S. Bony:
c  use another formula for the dependance of the sfc ocean albedo
c  on the solar zenith angle (derived from ECHAM-3; the original
c  doesn't seem reasonable to me...).
c======================================================================
      include "dimensions.h"
      include "dimphy.h"
      include "YOMCST.h"
c
      REAL fmagic ! un facteur magique pour regler l'albedo
ccc      PARAMETER (fmagic=0.7)
      PARAMETER (fmagic=1.0)
      INTEGER npts ! il controle la precision de l'integration
      PARAMETER (npts=120) ! 120 correspond a l'interval 6 minutes
c
      REAL rlat(klon), rjour, albedo(klon)
      REAL zdist, zlonsun, zpi, zdeclin
      REAL rmu,alb, srmu, salb, fauxo, aa, bb
      INTEGER i, k
c
      zpi = 4. * ATAN(1.)
c
c Calculer la longitude vraie de l'orbite terrestre:
      CALL orbite(rjour,zlonsun,zdist)
c
c Calculer la declinaison du soleil (qui varie entre + et - R_incl):
      zdeclin = ASIN(SIN(zlonsun*zpi/180.0)*SIN(R_incl*zpi/180.0))
c
      DO 999 i=1,klon
      aa = SIN(rlat(i)*zpi/180.0) * SIN(zdeclin)
      bb = COS(rlat(i)*zpi/180.0) * COS(zdeclin)
c
c Midi local (angle du temps = 0.0):
      rmu = aa + bb * COS(0.0)
      rmu = MAX(0.0, rmu)
c -- sb:
ccc      fauxo = (1.47-ACOS(rmu))/.15
ccc      alb = 0.03+0.630/(1.+fauxo*fauxo)
         alb = 0.05/(rmu + 0.15) ! echam-3 formula
         alb = MIN(alb,0.15)     ! echam-3
c sb --
      srmu = rmu
      salb = alb * rmu
c
c Faire l'integration numerique de midi a minuit (le facteur 2
c prend en compte l'autre moitie de la journee):
      DO k = 1, npts
         rmu = aa + bb * COS(FLOAT(k)/FLOAT(npts)*zpi)
         rmu = MAX(0.0, rmu)
c -- sb: 
ccc         fauxo = (1.47-ACOS(rmu))/.15
ccc         alb = 0.03+0.630/(1.+fauxo*fauxo)
         alb = 0.05 / (rmu+0.15) ! echam-3
         alb = MIN(alb,0.15)     ! echam-3
c sb --
         srmu = srmu + rmu * 2.0
         salb = salb + alb*rmu * 2.0
      ENDDO
      IF (srmu .NE. 0.0) THEN
         albedo(i) = salb / srmu * fmagic
      ELSE ! nuit polaire (on peut prendre une valeur quelconque)
         albedo(i) = fmagic
      ENDIF
  999 CONTINUE
      RETURN
      END
c=====================================================================
      SUBROUTINE alboc_cd(rmu0,albedo)
      IMPLICIT none
c======================================================================
c Auteur(s): Z.X. Li (LMD/CNRS)
c date: 19940624
c Calculer l'albedo sur l'ocean en fonction de l'angle zenithal moyen
c Formule due a Larson and Barkstrom (1977) Proc. of the symposium
C on radiation in the atmosphere, 19-28 August 1976, science Press,
C 1977 pp 451-453, ou These de 3eme cycle de Sylvie Joussaume.
c
c Arguments
c rmu0    (in): cosinus de l'angle solaire zenithal
c albedo (out): albedo de surface de l'ocean
c======================================================================
      include "dimensions.h"
      include "dimphy.h"
      REAL rmu0(klon), albedo(klon)
c
      REAL fauxo
      INTEGER i
c
      DO i = 1, klon
c -- sb:
ccc         fauxo = ( 1.47 - ACOS( rmu0(i) ) )/0.15
ccc         albedo(i) = 1.1*( .03 + .630/( 1. + fauxo*fauxo))
ccc         albedo(i) = MAX(MIN(albedo(i),0.60),0.04)
c -- sb: use echam-3 formula:
         albedo(i) = 0.05/(rmu0(i) + 0.15) 
         albedo(i) = MIN(albedo(i),0.15)
c sb --
      ENDDO
c
      RETURN
      END
c========================================================================
      SUBROUTINE albsno(veget, agesno, alb_neig)
      IMPLICIT none
c
      include "dimensions.h"
      include "dimphy.h"
      INTEGER nvm
      PARAMETER (nvm=8)
      REAL veget(klon,nvm)
      REAL alb_neig(klon)
      REAL agesno(klon)
c
      INTEGER i, nv
c
      REAL init(nvm), decay(nvm), as
      SAVE init, decay
      DATA init /0.55, 0.14, 0.18, 0.29, 0.15, 0.15, 0.14, 0./
      DATA decay/0.30, 0.67, 0.63, 0.45, 0.40, 0.14, 0.06, 1./
c
      DO i = 1, klon
         alb_neig(i) = 0.0
      ENDDO
      DO nv = 1, nvm
         DO i = 1, klon
            as = init(nv)+decay(nv)*EXP(-agesno(i)/5.)
            alb_neig(i) = alb_neig(i) + veget(i,nv)*as
         ENDDO
      ENDDO
c
      RETURN
      END
