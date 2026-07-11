!=======================================================================
! restart_io.F90
!
! Checkpoint/restore for the RESTART flag.
!
! ----------------------------------------------------------------------
! Two-file split (fully-converged cells vs. still-active cells)
! ----------------------------------------------------------------------
! The single-file checkpoint used to rewrite EVERY grid point on EVERY
! "level populations converged" checkpoint, including points that had
! already reached thermal balance (fullyconverged=.true.) many checkpoints
! earlier and will never change again. In a heavy 3D model this made the
! write cost scale with the FULL grid size on every checkpoint, right when
! the run is close to convergence and checkpoints are (relatively) frequent
! - the opposite of what should happen as fewer points remain active.
!
! This version splits the checkpoint into two files:
!
!   PREFIX_restart_A.bin  Cells with fullyconverged=.true. (thermal balance
!                          reached). Once a point enters this state,
!                          changetemperature/coolingfunctions/chemicaliterations
!                          all "cycle" past it unconditionally (see the
!                          fullyconverged guards in those routines), so its
!                          Tgas/abundance/cooling/heating/level-populations
!                          are frozen for the rest of the run. This file is
!                          therefore APPEND-ONLY: each point is written to it
!                          exactly once, the first time it converges, and is
!                          never rewritten. pdr(p)%restconverged (previously
!                          an unused field) marks "already appended to A" so
!                          repeat checkpoints don't duplicate it.
!
!   PREFIX_restart_B.bin  The full solver state (temperature bracket search,
!                          abundances, level populations, ...) for every
!                          point that is NOT YET fullyconverged. This file IS
!                          fully rewritten on every checkpoint - but it only
!                          ever contains the shrinking set of still-active
!                          points, so the rewrite gets cheaper as the model
!                          approaches full convergence, instead of staying
!                          pinned at the full grid size.
!
! Both files tag every record with the point index p (unformatted sequential
! records are read back with "do ... read(unit,iostat=ios) p ... if
! (ios/=0) exit"), so points do not need to be visited in order and neither
! file needs a separate index/manifest.
!
! Fields deliberately NOT stored in A (all are solver-internal bookkeeping
! for the not-yet-converged search; none are ever read again once
! fullyconverged, because every routine that would read them "cycle"s past
! fullyconverged points first):
!   doleveltmin, dobinarychop, previouschange, Tlow, Thigh, Illinois
!   Flow/Fhigh/ill_last, REBRACKET nrebracket, coolprev.
! coolprev in particular does not need to be reconstructed either: by the
! time a point stops being visited by coolingfunctions (the only place that
! writes coolprev), coolprev and cooling(k) already hold the same value (the
! last update wrote both from it), so it is simply re-derived on load as
! coolprev = cooling(k).
!
! levelconverged (per-point) is not stored in EITHER file: it is
! unconditionally reset to .false. and then recomputed from scratch by
! checkconvergence() for every point on every iteration (3DPDR.F90), so
! whatever value it has at load time is overwritten before anything reads
! it - keeping it in the checkpoint was dead weight.
!
! save_restart CLOSES both units after writing, so an interrupted run always
! leaves complete, flushed files (a half-written trailing record in A would
! only ever be the newest, not-yet-acknowledged batch: everything written in
! earlier checkpoints was already flushed and closed).
!
! A small header (pdr_ptot, coo) is validated - for BOTH files, before any
! part of pdr(:) is touched - on load; a mismatched or old-format pair is
! rejected and the run falls back to a normal cold start.
!=======================================================================

subroutine save_restart
  use healpix_types
  use maincode_module
  use global_module
  use maincode_local
  implicit none
#if defined RESTART && defined THERMALBALANCE
  ! p and k are the module-level loop counters (maincode_module)
  logical :: haveA
  integer(kind=i4b) :: nAnew, nBcur

  ! ---- A: append any cells that newly reached full (thermal-balance)
  !         convergence since the last checkpoint. Never rewritten. ----
  inquire(file=restart_file_A, exist=haveA)
  if (haveA) then
     open(unit=78, file=restart_file_A, status='old', form='unformatted', position='append')
  else
     open(unit=78, file=restart_file_A, status='replace', form='unformatted')
     write(78) pdr_ptot, coo
  endif

  nAnew = 0
  do p = 1, pdr_ptot
     if (pdr(p)%fullyconverged .and. .not.pdr(p)%restconverged) then
        write(78) p
        write(78) pdr(p)%Tgas, pdr(p)%nTgas
        write(78) pdr(p)%abundance
        write(78) pdr(p)%cooling
        write(78) pdr(p)%heating          ! output-only
        do k = 1, coo
           write(78) pdr(p)%coolant(k)%isconverged, pdr(p)%coolant(k)%pop
           write(78) pdr(p)%coolant(k)%line   ! output-only
        enddo
        pdr(p)%restconverged = .true.     ! don't write this point to A again
        nAnew = nAnew + 1
     endif
  enddo
  close(78)   ! flush and close so an interrupted run keeps a complete file

  ! ---- B: full rewrite, but only of the still-active (not yet
  !         fullyconverged) points - this set shrinks as the run proceeds. ----
  open(unit=77, file=restart_file_B, status='replace', form='unformatted')
  write(77) pdr_ptot, coo

  nBcur = 0
  do p = 1, pdr_ptot
     if (.not.pdr(p)%fullyconverged) then
        write(77) p
        write(77) pdr(p)%doleveltmin, pdr(p)%dobinarychop, pdr(p)%previouschange, &
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
           write(77) pdr(p)%coolant(k)%coolprev
        enddo
        nBcur = nBcur + 1
     endif
  enddo
  close(77)   ! flush and close so an interrupted run keeps a complete file

  write(6,'("  [RESTART] +",I0," cell(s) appended to _restart_A.bin (thermal-balance converged)")') nAnew
  write(6,'("  [RESTART] ",I0," active cell(s) written to _restart_B.bin")') nBcur
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
  ! p and k are the module-level loop counters (maincode_module), reused
  ! below only for the coolant (1:coo) loop; the point index read back from
  ! each record is kept in the local p_rec so this routine does not depend
  ! on visiting points in order.
  integer(kind=i4b) :: rpA, rcA, rpB, rcB, p_rec, ios, countA, countB
  logical :: haveA

  inquire(file=restart_file_A, exist=haveA)

  ! ---- Validate both headers BEFORE touching any pdr(:) state, so a
  !      mismatched (or old-format) pair falls back to a clean cold start
  !      rather than leaving pdr(:) partially overwritten. ----
  if (haveA) then
     open(unit=78, file=restart_file_A, status='old', form='unformatted')
     read(78) rpA, rcA
     if (rpA.ne.pdr_ptot .or. rcA.ne.coo) then
        write(6,*) ' *** [RESTART] _restart_A.bin does not match this model'// &
                 & ' (points/coolants differ); ignoring checkpoint ***'
        close(78)
        restart = .false.
        return
     endif
  endif

  open(unit=77, file=restart_file_B, status='old', form='unformatted')
  read(77) rpB, rcB
  if (rpB.ne.pdr_ptot .or. rcB.ne.coo) then
     write(6,*) ' *** [RESTART] _restart_B.bin does not match this model'// &
              & ' (points/coolants differ); ignoring checkpoint ***'
     close(77)
     if (haveA) close(78)
     restart = .false.
     return
  endif

  ! ---- A: fully-converged cells - restore the frozen final state. ----
  countA = 0
  if (haveA) then
     do
        read(78, iostat=ios) p_rec
        if (ios.ne.0) exit
        read(78) pdr(p_rec)%Tgas, pdr(p_rec)%nTgas
        read(78) pdr(p_rec)%abundance
        read(78) pdr(p_rec)%cooling
        read(78) pdr(p_rec)%heating
        do k = 1, coo
           read(78) pdr(p_rec)%coolant(k)%isconverged, pdr(p_rec)%coolant(k)%pop
           read(78) pdr(p_rec)%coolant(k)%line
           ! Make the solver state consistent with the restored populations
           ! (same as below for B) so nothing downstream sees a mismatch.
           pdr(p_rec)%coolant(k)%solution = pdr(p_rec)%coolant(k)%pop
           ! Re-derive rather than store: by the time a point stops being
           ! visited by coolingfunctions, coolprev already equals cooling(k)
           ! (the last call wrote both from the same value).
           pdr(p_rec)%coolant(k)%coolprev = pdr(p_rec)%cooling(k)
#ifdef NGRELAX
           pdr(p_rec)%coolant(k)%noscil = 0
#endif
        enddo
        pdr(p_rec)%totalcooling   = sum(pdr(p_rec)%cooling)
        pdr(p_rec)%fullyconverged = .true.
        pdr(p_rec)%restconverged  = .true.   ! already in A; don't re-append on next checkpoint
        countA = countA + 1
     enddo
     close(78)
  endif

  ! ---- B: still-active cells - restore the full solver state. ----
  countB = 0
  do
     read(77, iostat=ios) p_rec
     if (ios.ne.0) exit
     read(77) pdr(p_rec)%doleveltmin, pdr(p_rec)%dobinarychop, pdr(p_rec)%previouschange, &
            & pdr(p_rec)%Tgas, pdr(p_rec)%nTgas, pdr(p_rec)%Tlow, pdr(p_rec)%Thigh
#ifdef ILLINOIS
     read(77) pdr(p_rec)%Flow, pdr(p_rec)%Fhigh, pdr(p_rec)%ill_last
#endif
#ifdef REBRACKET
     read(77) pdr(p_rec)%nrebracket
#endif
     read(77) pdr(p_rec)%abundance
     read(77) pdr(p_rec)%cooling
     read(77) pdr(p_rec)%heating
     do k = 1, coo
        read(77) pdr(p_rec)%coolant(k)%isconverged, pdr(p_rec)%coolant(k)%pop
        read(77) pdr(p_rec)%coolant(k)%line
        read(77) pdr(p_rec)%coolant(k)%coolprev
        ! Make the solver state consistent with the restored populations so the
        ! first convergence check after restart compares like with like.
        pdr(p_rec)%coolant(k)%solution = pdr(p_rec)%coolant(k)%pop
#ifdef NGRELAX
        pdr(p_rec)%coolant(k)%noscil = 0
#endif
     enddo
     pdr(p_rec)%totalcooling = sum(pdr(p_rec)%cooling)
     countB = countB + 1
  enddo
  close(77)

  write(6,'("  [RESTART] loaded ",I0," fully-converged cell(s) from _restart_A.bin")') countA
  write(6,'("  [RESTART] loaded ",I0," active cell(s) from _restart_B.bin")') countB

  ! The Ng-acceleration history is not checkpointed; start it fresh.
  levpop_iteration = 0
#ifdef NGACCEL
  ng_nhist = 0
#endif
#endif
  return
end subroutine load_restart
