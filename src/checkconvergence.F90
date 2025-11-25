subroutine checkconvergence
use healpix_types
use maincode_module
use global_module

#ifdef OPENMP
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(p,ilevel,k,RELCH)
#endif 
do p=1,pdr_ptot
   pdr(p)%coolant(:)%isconverged=.true.
   do k=1,coo
      DO ilevel=1,coolant(k)%cnlev
         IF(pdr(p)%coolant(k)%solution(ilevel).GE.pdr(p)%abundance(coolant(k)%cspec)*1.0D-10) THEN
            IF(pdr(p)%coolant(k)%solution(ilevel).EQ.0.0D0 .AND. pdr(p)%coolant(k)%pop(ilevel).EQ.0.0D0) THEN
               RELCH=0.0D0
            ELSE
               RELCH=2.0D0*ABS((pdr(p)%coolant(k)%solution(ilevel)-pdr(p)%coolant(k)%pop(ilevel))&
                     &/(pdr(p)%coolant(k)%solution(ilevel)+pdr(p)%coolant(k)%pop(ilevel)))
            ENDIF
            IF(RELCH.GT.1.0D-2) then 
                RELCH_conv=.false.
                pdr(p)%coolant(k)%isconverged = .false.
            ENDIF        
         ENDIF
      ENDDO
   enddo
   pdr(p)%levelconverged = all([(pdr(p)%coolant(k)%isconverged, k = 1, coo)])
enddo !p=1,pdr_ptot
#ifdef OPENMP
!$OMP END PARALLEL DO
#endif

return
end subroutine
