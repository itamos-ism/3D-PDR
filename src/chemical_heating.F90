!=======================================================================
! chemical_heating.F90
!
! Automatic identification of the exothermic chemical-heating reactions
! used by CHEMICAL_HEATING in heatingfunctions.F90.
!
! Previously those reactions were referenced by hard-wired RATE(index)
! numbers, with a separate #ifdef block per network (REDUCED/MEDIUM/
! MEDIUMS/SCUTUM/FULL/MYNETWORK). Those indices change whenever a reaction
! moves position in rates_NETWORK.d, and a channel that is absent in a new
! network had to be deleted by hand - so every new network meant editing
! this routine.
!
! Here the heating reactions are defined ONCE by their chemistry
! (reactants -> products, with the exothermicity in eV), independent of the
! network. init_chemical_heating scans the loaded network and tags every
! reaction whose reactants AND products match a table entry. Matching is
! done on UNORDERED multisets of (trimmed) species names, so the order in
! which products are written does not matter. Then:
!   - reactions absent from the network are simply skipped (handles the
!     "2 instead of 3 channels" case the user mentioned);
!   - reactions split over several temperature ranges (DUPLICATE entries,
!     identical reactants+products) are ALL tagged - at any gas temperature
!     only the active one has a non-zero RATE, so the sum stays correct;
!   - the heating rate uses the consistent physical formula
!         sum_r  [ prod_j x(reactant_j) ] * density**n_reactants * RATE(r) * E_r
!     (number density of each reactant = abundance * density).
!
! To add a heating channel for a future network you edit the readable table
! in init_chemical_heating (reactants, products, energy) - never an index.
!
! Gated by the AUTOCHEMHEAT flag; with AUTOCHEMHEAT=0 the original
! per-network hard-wired blocks in heatingfunctions.F90 are used instead.
!=======================================================================
module chemical_heating_module
  use healpix_types
  implicit none

#ifdef AUTOCHEMHEAT
  private
  public :: init_chemical_heating, chemical_heating_rate, chemheat_ready

  ! Packed list of the tagged heating reactions (built once at startup).
  integer(kind=i4b)              :: nheat = 0      ! number of tagged reactions
  integer(kind=i4b), allocatable :: hrxn(:)        ! reaction (RATE) index
  real(kind=dp),     allocatable :: henergy(:)     ! exothermicity [eV]
  integer(kind=i4b), allocatable :: hs(:,:)        ! reactant species indices (:,1:3)
  integer(kind=i4b), allocatable :: hnr(:)         ! number of reactants
  logical :: chemheat_ready = .false.

contains

  !---------------------------------------------------------------------
  ! Build the tagged-reaction list by matching the network against the
  ! canonical table. Call once after READ_RATES / read_species.
  !---------------------------------------------------------------------
  subroutine init_chemical_heating(nreac,reactant,product,nspec,species)
    integer(kind=i4b), intent(in) :: nreac,nspec
    character(len=10), intent(in) :: reactant(nreac,3), product(nreac,4)
    character(len=10), intent(in) :: species(nspec)

    integer(kind=i4b), parameter :: nch = 12
    character(len=10) :: tr(nch,3), tp(nch,4)
    real(kind=dp)     :: te(nch)
    logical           :: found(nch)
    integer(kind=i4b) :: i,j,ic,nr,sidx,hs_tmp(3)
    logical           :: ok

    ! ---- canonical exothermic heating reactions (network-independent) ----
    !   See Clavel et al. (1978); Hollenbach & McKee. Energies in eV.
    tr = "          "; tp = "          "; te = 0.0d0
    !                        reactants                products                       E[eV]
    call setrx(tr,tp,te,nch, 1, "H2+","e-"," ",  "H","H"," "," ",      10.90d0) ! H2+  + e- -> H + H
    call setrx(tr,tp,te,nch, 2, "H2+","H"," ",   "H2","H+"," "," ",     0.94d0) ! H2+  + H  -> H2 + H+
    call setrx(tr,tp,te,nch, 3, "HCO+","e-"," ", "CO","H"," "," ",      7.51d0) ! HCO+ + e- -> CO + H
    call setrx(tr,tp,te,nch, 4, "H3+","e-"," ",  "H","H","H"," ",       4.76d0) ! H3+  + e- -> 3H
    call setrx(tr,tp,te,nch, 5, "H3+","e-"," ",  "H2","H"," "," ",      9.23d0) ! H3+  + e- -> H2 + H
    call setrx(tr,tp,te,nch, 6, "H3O+","e-"," ", "OH","H","H"," ",      1.16d0) ! H3O+ + e- -> OH + H + H
    call setrx(tr,tp,te,nch, 7, "H3O+","e-"," ", "OH","H2"," "," ",     5.63d0) ! H3O+ + e- -> OH + H2
    call setrx(tr,tp,te,nch, 8, "H3O+","e-"," ", "H2O","H"," "," ",     6.27d0) ! H3O+ + e- -> H2O + H
    call setrx(tr,tp,te,nch, 9, "H3O+","e-"," ", "O","H2","H"," ",      1.23d0) ! H3O+ + e- -> O + H2 + H
    call setrx(tr,tp,te,nch,10, "He+","H2"," ",  "He","H+","H"," ",     6.51d0) ! He+  + H2 -> He + H+ + H
    call setrx(tr,tp,te,nch,11, "He+","H2"," ",  "He","H2+"," "," ",    9.16d0) ! He+  + H2 -> He + H2+ (charge transfer)
    call setrx(tr,tp,te,nch,12, "He+","CO"," ",  "O","C+","He"," ",     2.22d0) ! He+  + CO -> O + C+ + He

    if (allocated(hrxn)) deallocate(hrxn,henergy,hs,hnr)
    allocate(hrxn(nreac),henergy(nreac),hs(nreac,3),hnr(nreac))
    hrxn=0; henergy=0.0d0; hs=0; hnr=0
    nheat=0
    found=.false.

    do i=1,nreac
       do ic=1,nch
          if (multiset_eq(reactant(i,:),tr(ic,:),3) .and. &
            & multiset_eq(product(i,:), tp(ic,:),4)) then
             ! count reactants and map each to its species index
             nr=0; ok=.true.; hs_tmp=0
             do j=1,3
                if (len_trim(adjustl(reactant(i,j))).gt.0) then
                   nr=nr+1
                   sidx=spec_index(reactant(i,j),species,nspec)
                   if (sidx.le.0) ok=.false.
                   hs_tmp(nr)=sidx
                endif
             enddo
             if (.not.ok) then
                write(6,'(" [AUTOCHEMHEAT] WARNING: reactant not in species list; skipping reaction ",I0)') i
                cycle
             endif
             nheat=nheat+1
             hrxn(nheat)=i
             henergy(nheat)=te(ic)
             hnr(nheat)=nr
             hs(nheat,1:3)=hs_tmp(1:3)
             found(ic)=.true.
             write(6,'("   rxn ",I5,":  ",A,1X,A,1X,A," -> ",A,1X,A,1X,A,1X,A,3X,F6.2," eV")') &
                  & i, trim(reactant(i,1)),trim(reactant(i,2)),trim(reactant(i,3)), &
                  & trim(product(i,1)),trim(product(i,2)),trim(product(i,3)),trim(product(i,4)),te(ic)
             exit
          endif
       enddo
    enddo

    write(6,'(" [AUTOCHEMHEAT] tagged ",I0," chemical-heating reaction(s) in this network")') nheat
    do ic=1,nch
       if (.not.found(ic)) &
          write(6,'("   - channel not present in this network: ",A," + ",A," (",F6.2," eV)")') &
               & trim(tr(ic,1)),trim(tr(ic,2)),te(ic)
    enddo
    chemheat_ready=.true.

  end subroutine init_chemical_heating


  !---------------------------------------------------------------------
  ! Total exothermic chemical heating [erg cm^-3 s^-1] at a grid point.
  !---------------------------------------------------------------------
  function chemical_heating_rate(abundance,density,rate,nspec,nreac) result(ch)
    integer(kind=i4b), intent(in) :: nspec,nreac
    real(kind=dp),     intent(in) :: abundance(nspec), rate(nreac), density
    real(kind=dp) :: ch, term
    integer(kind=i4b) :: m,j
    ch=0.0d0
    if (.not.chemheat_ready) return
    do m=1,nheat
       term = rate(hrxn(m))*henergy(m)*EV
       do j=1,hnr(m)
          term = term*abundance(hs(m,j))
       enddo
       term = term*density**hnr(m)
       ch = ch + term
    enddo
  end function chemical_heating_rate


  !---- helpers -------------------------------------------------------
  subroutine setrx(tr,tp,te,nch,k,r1,r2,r3,p1,p2,p3,p4,e)
    integer(kind=i4b), intent(in)    :: nch,k
    character(len=10), intent(inout) :: tr(nch,3), tp(nch,4)
    real(kind=dp),     intent(inout) :: te(nch)
    character(len=*),  intent(in)    :: r1,r2,r3,p1,p2,p3,p4
    real(kind=dp),     intent(in)    :: e
    tr(k,1)=r1; tr(k,2)=r2; tr(k,3)=r3
    tp(k,1)=p1; tp(k,2)=p2; tp(k,3)=p3; tp(k,4)=p4
    te(k)=e
  end subroutine setrx

  integer(kind=i4b) function spec_index(name,species,nspec)
    character(len=*),  intent(in) :: name
    integer(kind=i4b), intent(in) :: nspec
    character(len=10), intent(in) :: species(nspec)
    integer(kind=i4b) :: i
    spec_index=0
    do i=1,nspec
       if (trim(adjustl(species(i))).eq.trim(adjustl(name))) then
          spec_index=i; return
       endif
    enddo
  end function spec_index

  ! .true. if the non-blank entries of a and b are equal as unordered multisets
  logical function multiset_eq(a,b,n)
    integer(kind=i4b), intent(in) :: n
    character(len=10), intent(in) :: a(n), b(n)
    logical           :: usedb(n)
    integer(kind=i4b) :: i,j,na,nb
    na=0; nb=0; usedb=.false.
    do i=1,n
       if (len_trim(adjustl(a(i))).gt.0) na=na+1
       if (len_trim(adjustl(b(i))).gt.0) nb=nb+1
    enddo
    multiset_eq=.false.
    if (na.ne.nb) return
    do i=1,n
       if (len_trim(adjustl(a(i))).eq.0) cycle
       do j=1,n
          if (.not.usedb(j)) then
             if (trim(adjustl(b(j))).eq.trim(adjustl(a(i)))) then
                usedb(j)=.true.
                go to 10
             endif
          endif
       enddo
       return        ! a(i) had no partner in b
10     continue
    enddo
    multiset_eq=.true.
  end function multiset_eq
#endif

end module chemical_heating_module
