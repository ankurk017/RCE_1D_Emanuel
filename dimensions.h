c-----------------------------------------------------------------------
c   INCLUDE 'dimensions.h'
c
c   dimensions.h contient les dimensions du modele
c   ndm est tel que iim=2**ndm
c   nqmx est la dimension de la variable traceur q
c-----------------------------------------------------------------------

      INTEGER iim,jjm,llm,ndm
      PARAMETER (iim= 1,jjm=1,llm=46,ndm=1) ! 25mb
c      PARAMETER (iim= 1,jjm=1,llm=21,ndm=1) ! 50mb
      integer nqmx
      parameter (nqmx=3)

c-----------------------------------------------------------------------