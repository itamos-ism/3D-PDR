!=======================================================================
! restart_io.F90
!
! Checkpoint/restore for the RESTART flag.
!
! The original restart only stored {fullyconverged, Tgas, abundance} per
! grid point. That is NOT enough to resume a converged solution: on restart
! the level populations were rebuilt from LTE (chemicaliterations(1,...)) and
! had to re-converge through the whole level-population + thermal-balance
! cycle, so even a fully converged model needed many (>20) iterations.
!
! These routines checkpoint the full per-point state needed to continue from
! the last "level populations converged" step:
!   - gas temperature (Tgas, nTgas) and the thermal-balance bracket
!     (Tlow, Thigh, previouschange, dobinarychop, doleveltmin, and the
!      Illinois/REBRACKET search state),
!   - the chemical abundances,
!   - the per-coolant LEVEL POPULATIONS (the piece that was missing),
!   - the cooling rates and convergence flags.
!
! save_restart writes the file and CLOSES it every checkpoint, so a run that
! is interrupted always leaves a complete, flushed file (the original kept the
! unit open and only truncated/re-wrote at the next checkpoint, so the most
! recent checkpoint could be left unflushed).
!
! A small header (pdr_ptot, coo) is validated on load; a mismatched or
! old-format file is rejected and the run falls back to a normal cold start.
!=======================================================================

subroutine save_restart
  use healpix_types
  use maincode_module
  use global_module
  use maincode_local
  implicit none
#if defined RESTART && defined THERMALBALANCE
  ! p and k are the module-level loop counters (maincode_module)
  close(77)   ! no-op if not connected
  open(unit=77, file=restart_file, status='replace', form='unformatted')

  ! Header for validation on load
  write(77) pdr_ptot, coo

  do p = 1, pdr_ptot
     write(77) pdr(p)%fullyconverged, pdr(p)%levelconverged, pdr(p)%doleveltmin, &
             & pdr(p)%dobinarychop, pdr(p)%previouschange, &
             & pdr(p)%Tgas, pdr(p)%nTgas, pdr(p)%Tlow, pdr(p)%Thigh
#ifdef ILLINOIS
     write(77) pdr(p)%Flow, pdr(p)%Fhigh, pdr(p)%ill_last
#endif
#ifdef REBRACKET
     write(77) pdr(p)%nrebracket
#endif
     write(77) pdr(p)%abundance
     write(77) pdr(p)%cooling
     write(77) pdr(p)%heating          ! output-only; converged points skip changetemperature
     do k = 1, coo
        write(77) pdr(p)%coolant(k)%isconverged, pdr(p)%coolant(k)%pop
        write(77) pdr(p)%coolant(k)%line   ! output-only; converged points skip coolingfunctions
#ifdef FREEZECOOL
        write(77) pdr(p)%coolant(k)%coolprev
#endif
     enddo
  enddo

  close(77)   ! flush and close so an interrupted run keeps a complete file
#endif
  return
end subroutine save_restart


subroutine load_restart
  use healpix_types
  use maincode_module
  use global_module
  use maincode_local
  implicit none
#if defined RESTART && defined THERMALBALANCE
  ! p and k are the module-level loop counters (maincode_module)
  integer(kind=i4b) :: rp, rc

  open(unit=77, file=restart_file, status='old', form='unformatted')

  ! Validate the header. A mismatched (or old-format) file is rejected and the
  ! run continues as a normal cold start; pdr(:) is left untouched.
  read(77) rp, rc
  if (rp.ne.pdr_ptot .or. rc.ne.coo) then
     write(6,*) ' *** [RESTART] file does not match this model'// &
              & ' (points/coolants differ); ignoring it ***'
     close(77)
     restart = .false.
     return
  endif

  do p = 1, pdr_ptot
     read(77) pdr(p)%fullyconverged, pdr(p)%levelconverged, pdr(p)%doleveltmin, &
            & pdr(p)%dobinarychop, pdr(p)%previouschange, &
            & pdr(p)%Tgas, pdr(p)%nTgas, pdr(p)%Tlow, pdr(p)%Thigh
#ifdef ILLINOIS
     read(77) pdr(p)%Flow, pdr(p)%Fhigh, pdr(p)%ill_last
#endif
#ifdef REBRACKET
     read(77) pdr(p)%nrebracket
#endif
     read(77) pdr(p)%abundance
     read(77) pdr(p)%cooling
     read(77) pdr(p)%heating
     do k = 1, coo
        read(77) pdr(p)%coolant(k)%isconverged, pdr(p)%coolant(k)%pop
        read(77) pdr(p)%coolant(k)%line
#ifdef FREEZECOOL
        read(77) pdr(p)%coolant(k)%coolprev
#endif
        ! Make the solver state consistent with the restored populations so the
        ! first convergence check after restart compares like with like.
        pdr(p)%coolant(k)%solution = pdr(p)%coolant(k)%pop
#ifdef NGRELAX
        pdr(p)%coolant(k)%noscil = 0
#endif
     enddo
     pdr(p)%totalcooling = sum(pdr(p)%cooling)
  enddo

  close(77)

  ! The Ng-acceleration history is not checkpointed; start it fresh.
  levpop_iteration = 0
#ifdef NGACCEL
  ng_nhist = 0
#endif
#endif
  return
end subroutine load_restart
