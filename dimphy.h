c-----------------------------------------------------------------------
      INTEGER KIDIA, KFDIA, KLON, KLEV
c sb      PARAMETER (KIDIA=1,KFDIA=iim*(jjm-1)+2,
      PARAMETER (KIDIA=1,KFDIA=iim*(jjm-1)+2-1/jjm,
     .           KLON=KFDIA-KIDIA+1,KLEV=llm)
c-----------------------------------------------------------------------
      INTEGER nbtr ! nombre de vrais traceurs
      PARAMETER (nbtr=nqmx-2+1/(nqmx-1))
c-----------------------------------------------------------------------
