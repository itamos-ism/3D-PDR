subroutine readcoolants
use maincode_module
#ifdef LEVMAX
real(kind=dp), pointer :: lm_p1(:), lm_p2(:,:), lm_p3(:,:,:)
#endif

allocate(coolant(1:coo))
do i = 1,coo
  open(unit=44,file=coolfile(i),status='old')
  read(44,'()') 
  read(44,*) coolant(i)%cname
  read(44,'()')
  read(44,*) coolant(i)%molweight
  read(44,'()')
  read(44,*) cur_nlev, cur_ntemp
  close(44)
  coolant(i)%cnlev = cur_nlev
  coolant(i)%cntemp = cur_ntemp
  allocate(coolant(i)%energies(1:cur_nlev))
  allocate(coolant(i)%weights(1:cur_nlev))
  allocate(coolant(i)%A_COEFFS(1:cur_nlev,1:cur_nlev))
  allocate(coolant(i)%B_COEFFS(1:cur_nlev,1:cur_nlev))
  allocate(coolant(i)%frequencies(1:cur_nlev,1:cur_nlev))
  allocate(coolant(i)%temperatures(1:7,1:cur_ntemp))
  allocate(coolant(i)%H_COL(1:cur_NLEV,1:cur_NLEV,1:cur_NTEMP))
  allocate(coolant(i)%HP_COL(1:cur_nlev,1:cur_nlev,1:cur_ntemp))
  allocate(coolant(i)%EL_COL(1:cur_NLEV,1:cur_NLEV,1:cur_NTEMP))
  allocate(coolant(i)%HE_COL(1:cur_NLEV,1:cur_NLEV,1:cur_NTEMP))
  allocate(coolant(i)%H2_COL(1:cur_NLEV,1:cur_NLEV,1:cur_NTEMP))
  allocate(coolant(i)%PH2_COL(1:cur_NLEV,1:cur_NLEV,1:cur_NTEMP))
  allocate(coolant(i)%OH2_COL(1:cur_NLEV,1:cur_NLEV,1:cur_NTEMP))
  call readinput(coolfile(i),cur_nlev,cur_ntemp,coolant(i)%energies,coolant(i)%weights,&
          coolant(i)%A_COEFFS,coolant(i)%B_COEFFS,coolant(i)%frequencies,coolant(i)%temperatures,&
          coolant(i)%H_COL,coolant(i)%HP_COL,coolant(i)%EL_COL,coolant(i)%HE_COL,coolant(i)%H2_COL,&
          coolant(i)%PH2_COL,coolant(i)%OH2_COL)
#ifdef LEVMAX
  if (lmax_levels.gt.0 .and. cur_nlev.gt.lmax_levels) then
    ! 1D: energies, weights
    allocate(lm_p1(1:lmax_levels))
    lm_p1 = coolant(i)%energies(1:lmax_levels)
    deallocate(coolant(i)%energies); coolant(i)%energies => lm_p1; nullify(lm_p1)
    allocate(lm_p1(1:lmax_levels))
    lm_p1 = coolant(i)%weights(1:lmax_levels)
    deallocate(coolant(i)%weights); coolant(i)%weights => lm_p1; nullify(lm_p1)
    ! 2D: A_COEFFS, B_COEFFS, frequencies
    allocate(lm_p2(1:lmax_levels,1:lmax_levels))
    lm_p2 = coolant(i)%A_COEFFS(1:lmax_levels,1:lmax_levels)
    deallocate(coolant(i)%A_COEFFS); coolant(i)%A_COEFFS => lm_p2; nullify(lm_p2)
    allocate(lm_p2(1:lmax_levels,1:lmax_levels))
    lm_p2 = coolant(i)%B_COEFFS(1:lmax_levels,1:lmax_levels)
    deallocate(coolant(i)%B_COEFFS); coolant(i)%B_COEFFS => lm_p2; nullify(lm_p2)
    allocate(lm_p2(1:lmax_levels,1:lmax_levels))
    lm_p2 = coolant(i)%frequencies(1:lmax_levels,1:lmax_levels)
    deallocate(coolant(i)%frequencies); coolant(i)%frequencies => lm_p2; nullify(lm_p2)
    ! 3D: all 7 collision partner tables
    allocate(lm_p3(1:lmax_levels,1:lmax_levels,1:cur_ntemp))
    lm_p3 = coolant(i)%H_COL(1:lmax_levels,1:lmax_levels,1:cur_ntemp)
    deallocate(coolant(i)%H_COL); coolant(i)%H_COL => lm_p3; nullify(lm_p3)
    allocate(lm_p3(1:lmax_levels,1:lmax_levels,1:cur_ntemp))
    lm_p3 = coolant(i)%HP_COL(1:lmax_levels,1:lmax_levels,1:cur_ntemp)
    deallocate(coolant(i)%HP_COL); coolant(i)%HP_COL => lm_p3; nullify(lm_p3)
    allocate(lm_p3(1:lmax_levels,1:lmax_levels,1:cur_ntemp))
    lm_p3 = coolant(i)%EL_COL(1:lmax_levels,1:lmax_levels,1:cur_ntemp)
    deallocate(coolant(i)%EL_COL); coolant(i)%EL_COL => lm_p3; nullify(lm_p3)
    allocate(lm_p3(1:lmax_levels,1:lmax_levels,1:cur_ntemp))
    lm_p3 = coolant(i)%HE_COL(1:lmax_levels,1:lmax_levels,1:cur_ntemp)
    deallocate(coolant(i)%HE_COL); coolant(i)%HE_COL => lm_p3; nullify(lm_p3)
    allocate(lm_p3(1:lmax_levels,1:lmax_levels,1:cur_ntemp))
    lm_p3 = coolant(i)%H2_COL(1:lmax_levels,1:lmax_levels,1:cur_ntemp)
    deallocate(coolant(i)%H2_COL); coolant(i)%H2_COL => lm_p3; nullify(lm_p3)
    allocate(lm_p3(1:lmax_levels,1:lmax_levels,1:cur_ntemp))
    lm_p3 = coolant(i)%PH2_COL(1:lmax_levels,1:lmax_levels,1:cur_ntemp)
    deallocate(coolant(i)%PH2_COL); coolant(i)%PH2_COL => lm_p3; nullify(lm_p3)
    allocate(lm_p3(1:lmax_levels,1:lmax_levels,1:cur_ntemp))
    lm_p3 = coolant(i)%OH2_COL(1:lmax_levels,1:lmax_levels,1:cur_ntemp)
    deallocate(coolant(i)%OH2_COL); coolant(i)%OH2_COL => lm_p3; nullify(lm_p3)
    coolant(i)%cnlev = lmax_levels
    cur_nlev = lmax_levels
    write(6,'(a,a,a,i4)') ' [LEVMAX] Coolant ', trim(coolant(i)%cname), &
        ': levels truncated to ', lmax_levels
  endif
#endif
enddo
return
end subroutine
