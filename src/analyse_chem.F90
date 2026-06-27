!=======================================================================
      SUBROUTINE ANALYSE_CHEMISTRY(POS,LOCRAD,LOCCR,GRIDPOINT,TIME,DENSITY,TEMPERATURE, &
     &                             NRAYS,AV,NSPEC,SPECIES,ABUNDANCE, &
     &                             NREAC,REACTANT,PRODUCT,RATE)

!T.Bell
!
!     Outputs the dominant formation and destruction pathways of EVERY
!     species in the active chemical network, for every grid point.
!
!     The original routine only examined a hardcoded list of 11 species
!     and wrote a verbose, human-readable report. Because we now emit the
!     whole network the file would be enormous in that format, so the
!     output below is deliberately compact and machine-parseable. It is
!     read back by PDR-studio (pdrstudio/chemanalysis.py) for the
!     Chemical Analysis feature.
!
!     File layout (whitespace-separated, '#' comment lines are headers):
!
!       # 3D-PDR chemanalysis v1
!       # nspec=<N> nreac=<M> time=<t> yr
!       R <id> <reactant tokens ...> > <product tokens ...>   (once, M lines)
!       P <point> <x_pc> <Av> <nH> <Tgas> <FUV> <CR_s-1>      (per grid point)
!       S <id> <species> <abundance> <Ptot_s-1> <Dtot_s-1>    (per species)
!       <reaction_id> <signed_percent>   (+ = formation, - = destruction)
!
!     The reaction legend (R lines) is written once, on the first grid
!     point, so the bulky reaction strings are never repeated. Each
!     reaction is referenced thereafter by its 1-based index <id>, which
!     is the reaction's record order in chemfiles/rates_*.d.

      USE DEFINITIONS
      USE HEALPIX_TYPES
      USE GLOBAL_MODULE

      IMPLICIT NONE

      INTEGER(KIND=I4B), INTENT(IN) :: GRIDPOINT,NRAYS,NSPEC,NREAC
      REAL(KIND=DP),     INTENT(IN) :: POS,LOCRAD,LOCCR,TIME,DENSITY,TEMPERATURE
      REAL(KIND=DP),     INTENT(IN) :: AV,ABUNDANCE(1:NSPEC),RATE(1:NREAC)
      CHARACTER(LEN=10), INTENT(IN) :: SPECIES(1:NSPEC),REACTANT(1:NREAC,1:3),PRODUCT(1:NREAC,1:4)

      INTEGER(KIND=I4B) :: I,J,K,M,N,IP,ID
      INTEGER(KIND=I4B) :: NPR(1:NREAC),NDR(1:NREAC)
      INTEGER(KIND=I4B) :: NMAX(1),TOTAL
      REAL(KIND=DP)     :: X1,X2,RMULT,PERCENT
      REAL(KIND=DP)     :: P(1:NREAC),PTOT
      REAL(KIND=DP)     :: D(1:NREAC),DTOT
      CHARACTER(LEN=256):: LINEBUF

!-----------------------------------------------------------------------
!     On the first grid point, write the file header and the reaction
!     legend. The legend maps each reaction index to its equation, so the
!     per-point data below only ever has to quote the (integer) index.
!-----------------------------------------------------------------------
      IF(GRIDPOINT.EQ.1) THEN
         WRITE(98,'(A)') '# 3D-PDR chemanalysis v1'
         WRITE(98,'(A,I0,A,I0,A,ES12.5,A)') &
     &        '# nspec=',NSPEC,' nreac=',NREAC,' time=',TIME,' yr'
         WRITE(98,'(A)') '# R <id> <reactants ...> > <products ...>'
         WRITE(98,'(A)') '# P <point> <x_pc> <Av> <nH> <Tgas> <FUV> <CR_s-1>'
         WRITE(98,'(A)') '# S <id> <species> <abundance> <Ptot_s-1> <Dtot_s-1>'
         WRITE(98,'(A)') '# <reaction_id> <signed_percent>  (+ formation, - destruction)'

         DO J=1,NREAC
            WRITE(LINEBUF,'(A,I0)') 'R ',J
            DO K=1,3
               IF(TRIM(ADJUSTL(REACTANT(J,K))).NE.'') &
     &            LINEBUF=TRIM(LINEBUF)//' '//TRIM(ADJUSTL(REACTANT(J,K)))
            ENDDO
            LINEBUF=TRIM(LINEBUF)//' >'
            DO K=1,4
               IF(TRIM(ADJUSTL(PRODUCT(J,K))).NE.'') &
     &            LINEBUF=TRIM(LINEBUF)//' '//TRIM(ADJUSTL(PRODUCT(J,K)))
            ENDDO
            WRITE(98,'(A)') TRIM(LINEBUF)
         ENDDO
      ENDIF

!-----------------------------------------------------------------------
!     Per grid point header: position, environment and physical state.
!-----------------------------------------------------------------------
      WRITE(98,'(A,I0,6(1X,ES12.5))') 'P ',GRIDPOINT, &
     &     POS,AV,DENSITY,TEMPERATURE,LOCRAD,LOCCR*1.3D-17

!-----------------------------------------------------------------------
!     Loop over EVERY species in the network.
!-----------------------------------------------------------------------
      DO I=1,NSPEC

!        Reset the formation/destruction rate counters
         IP=0
         ID=0

!        Reset the total formation/destruction rates
         PTOT=0.0D0
         DTOT=0.0D0

!        Check all reactions to find relevant ones
         DO J=1,NREAC
!           Reset the reactant abundances
            X1=0.0D0
            X2=0.0D0

!-----------------------------------------------------------------------
!           Formation routes for species I
!-----------------------------------------------------------------------
            IF(PRODUCT(J,1).EQ.SPECIES(I) .OR. PRODUCT(J,2).EQ.SPECIES(I) .OR. &
     &         PRODUCT(J,3).EQ.SPECIES(I) .OR. PRODUCT(J,4).EQ.SPECIES(I)) THEN

               IP=IP+1

!              Store the reactant abundances
               DO N=1,NSPEC
                  IF(SPECIES(N).EQ.REACTANT(J,1)) X1=ABUNDANCE(N)
                  IF(SPECIES(N).EQ.REACTANT(J,2)) X2=ABUNDANCE(N)
               ENDDO

               RMULT=0.0D0
               IF(PRODUCT(J,1).EQ.SPECIES(I)) RMULT=RMULT+1.0D0
               IF(PRODUCT(J,2).EQ.SPECIES(I)) RMULT=RMULT+1.0D0
               IF(PRODUCT(J,3).EQ.SPECIES(I)) RMULT=RMULT+1.0D0
               IF(PRODUCT(J,4).EQ.SPECIES(I)) RMULT=RMULT+1.0D0

!              Calculate the reaction rate
               IF(REACTANT(J,2).EQ."PHOTON" .OR. &
     &            REACTANT(J,2).EQ."CRP   " .OR. &
     &            REACTANT(J,2).EQ."CRPHOT" .OR. &
     &            REACTANT(J,2).EQ."FREEZE" .OR. &
     &            REACTANT(J,2).EQ."ELFRZE" .OR. &
     &            REACTANT(J,2).EQ."CRH   " .OR. &
     &            REACTANT(J,2).EQ."PHOTD " .OR. &
     &            REACTANT(J,2).EQ."THERM ") THEN
                  P(IP)=RATE(J)*X1*RMULT
               ELSE
                  P(IP)=RATE(J)*X1*X2*DENSITY*RMULT
               ENDIF

               PTOT=PTOT+P(IP)
               NPR(IP)=J
            ENDIF

!-----------------------------------------------------------------------
!           Destruction routes for species I
!-----------------------------------------------------------------------

            IF(REACTANT(J,1).EQ.SPECIES(I) .OR. REACTANT(J,2).EQ.SPECIES(I)) THEN

               ID=ID+1

!              Store the reactant abundances
               DO N=1,NSPEC
                  IF(SPECIES(N).EQ.REACTANT(J,1)) X1=ABUNDANCE(N)
                  IF(SPECIES(N).EQ.REACTANT(J,2)) X2=ABUNDANCE(N)
               ENDDO

               RMULT=0.0D0
               IF(REACTANT(J,1).EQ.SPECIES(I)) RMULT=RMULT+1.0D0
               IF(REACTANT(J,2).EQ.SPECIES(I)) RMULT=RMULT+1.0D0

!              Calculate the reaction rate
               IF(REACTANT(J,2).EQ."PHOTON" .OR. &
     &            REACTANT(J,2).EQ."CRP   " .OR. &
     &            REACTANT(J,2).EQ."CRPHOT" .OR. &
     &            REACTANT(J,2).EQ."FREEZE" .OR. &
     &            REACTANT(J,2).EQ."ELFRZE" .OR. &
     &            REACTANT(J,2).EQ."CRH   " .OR. &
     &            REACTANT(J,2).EQ."PHOTD " .OR. &
     &            REACTANT(J,2).EQ."THERM ") THEN
                  D(ID)=RATE(J)*X1*RMULT
               ELSE
                  D(ID)=RATE(J)*X1*X2*DENSITY*RMULT
               ENDIF

               DTOT=DTOT+D(ID)
               NDR(ID)=J
            ENDIF

!        End of loop over all reactions
         ENDDO

!        Prevent divide-by-zero errors by setting minimum finite
!        values for the total formation and destruction rates
         IF(PTOT.LT.1.0D-99) PTOT=1.0D-99
         IF(DTOT.LT.1.0D-99) DTOT=1.0D-99

!-----------------------------------------------------------------------
!        Species record: index, name, abundance and total rates.
!-----------------------------------------------------------------------
         WRITE(98,'(A,I0,1X,A,3(1X,ES12.5))') 'S ',I, &
     &        TRIM(ADJUSTL(SPECIES(I))),ABUNDANCE(I),PTOT,DTOT

!-----------------------------------------------------------------------
!        Output the formation reactions (positive percentages),
!        in order of decreasing importance, until they sum to 100%.
!-----------------------------------------------------------------------
         NMAX=0
         TOTAL=0
         DO M=1,IP
            NMAX=MAXLOC(P(1:IP))
            N=NMAX(1)
            IF(P(N).EQ.0.0D0) EXIT

            PERCENT=1.0D2*(P(N)/PTOT)
            IF(PERCENT.LT.0.5D0) PERCENT=1.0D0
            TOTAL=TOTAL+NINT(PERCENT)

            WRITE(98,'(I0,1X,I0)') NPR(N),NINT(PERCENT)

            IF(TOTAL.GE.100) EXIT
            P(N)=0.0D0
         ENDDO

!-----------------------------------------------------------------------
!        Output the destruction reactions (negative percentages),
!        in order of decreasing importance, until they sum to 100%.
!-----------------------------------------------------------------------
         NMAX=0
         TOTAL=0
         DO M=1,ID
            NMAX=MAXLOC(D(1:ID))
            N=NMAX(1)
            IF(D(N).EQ.0.0D0) EXIT

            PERCENT=1.0D2*(D(N)/DTOT)
            IF(PERCENT.LT.0.5D0) PERCENT=1.0D0
            TOTAL=TOTAL+NINT(PERCENT)

            WRITE(98,'(I0,1X,I0)') NDR(N),-NINT(PERCENT)

            IF(TOTAL.GE.100) EXIT
            D(N)=0.0D0
         ENDDO

!     End of loop over species
      ENDDO

      RETURN
      END SUBROUTINE ANALYSE_CHEMISTRY
!=======================================================================
