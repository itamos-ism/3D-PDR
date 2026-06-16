subroutine ng_accelerate(nlev, xnew, h1, h2, h3)
! Ng (1974) acceleration of the level-population fixed-point iteration.
! xnew is the latest plain iterate; h1,h2,h3 are the three previous iterates
! (h1 the most recent). On success xnew is replaced by the third-order
! accelerated estimate; on any degeneracy (singular least-squares system,
! negative or non-finite populations) xnew is left untouched so the plain
! iteration proceeds as before.
use definitions
use healpix_types
implicit none

integer(kind=i4b), intent(in) :: nlev
real(kind=dp), intent(inout)  :: xnew(1:nlev)
real(kind=dp), intent(in)     :: h1(1:nlev), h2(1:nlev), h3(1:nlev)

real(kind=dp) :: r0(1:nlev), r1(1:nlev), u(1:nlev), v(1:nlev), w(1:nlev), xacc(1:nlev)
real(kind=dp) :: a11, a12, a22, b1, b2, det, ca, cb, floorval, maxr
#ifdef NGCYCLE
real(kind=dp) :: res2
#endif
integer(kind=i4b) :: i

! Relative weighting, with levels far below the convergence-check threshold
! (1e-10 of the total population) capped so numerical noise cannot dominate
floorval = 1.0D-10*sum(xnew)
if (floorval.le.0.0D0) return

r0 = xnew - h1
r1 = h1 - h2
u  = xnew - 2.0D0*h1 + h2
v  = xnew - h1 - h2 + h3
do i=1,nlev
   w(i) = 1.0D0/max(xnew(i),floorval)**2
enddo

! Nothing to gain if the latest update is already below the 1% convergence
! threshold used by checkconvergence: the plain iteration is about to
! converge on its own, and extrapolating noise-level residuals can only
! un-converge it
maxr = 0.0D0
do i=1,nlev
   maxr = max(maxr, abs(r0(i))/max(xnew(i),floorval))
enddo
if (maxr.le.1.0D-2) return

! Ng assumes a (near-)linear contractive iteration; right after a
! temperature or chemistry change the sequence is transient and the
! extrapolation is unreliable, so require the residual to be shrinking
if (sum(w*r0*r0).gt.sum(w*r1*r1)) return

! Least-squares solve for (ca,cb) minimizing ||r0 - ca*u - cb*v||_w
a11 = sum(w*u*u)
a12 = sum(w*u*v)
a22 = sum(w*v*v)
b1  = sum(w*u*r0)
b2  = sum(w*v*r0)
det = a11*a22 - a12*a12
if (det.le.1.0D-14*a11*a22) then
#ifdef NGCYCLE
   ! Degenerate least-squares system. This happens in particular for a pure
   ! period-2 oscillation (x0=h2, h1=h3 makes v vanish identically), which
   ! the two-coefficient extrapolation cannot treat but the one-coefficient
   ! (Aitken-like) form solves exactly: for a linear flip-flop mode it
   ! returns the cycle average, i.e. the fixed point. Such 2-cycles otherwise
   ! persist just above the 1% convergence threshold and stall the level
   ! populations until FORCECONVERGENCE truncates the stretch.
   if (a11.le.0.0D0) return
   ca = b1/a11
   ! Trust a large coefficient only when the one-mode model actually
   ! explains the residual: a slow mode (contraction factor lambda close
   ! to 1) NEEDS ca = lambda/(lambda-1), which is large; capping it leaves
   ! points crawling at ~1%/iteration just above the convergence threshold
   ! for a hundred iterations. A poor fit with a large coefficient still
   ! signals noise and is rejected as before.
   res2 = sum(w*(r0-ca*u)**2)
   if (res2.gt.1.0D-2*sum(w*r0*r0) .and. abs(ca).gt.10.0D0) return
   if (abs(ca).gt.1.0D3) return
   xacc = (1.0D0-ca)*xnew + ca*h1
   do i=1,nlev
      if (xacc(i).ne.xacc(i)) return            ! NaN: reject the whole step
      ! A large (legitimate) extrapolation amplifies the noise of already-
      ! converged levels and can push them slightly negative; keep those
      ! levels un-extrapolated instead of rejecting the whole step, which
      ! would block exactly the slow-mode acceleration this branch is for.
      if (xacc(i).lt.0.0D0) xacc(i) = xnew(i)
   enddo
   xnew = xacc
#endif
   return
endif

ca = (b1*a22 - b2*a12)/det
cb = (b2*a11 - b1*a12)/det

! Reject wild extrapolations (large coefficients signal a noisy or
! non-linear regime where the least-squares fit is meaningless)
#ifdef NGCYCLE
! ...unless the two-mode model explains the residual well: slow modes
! (contraction factor close to 1) need large coefficients, and capping
! them leaves points crawling at ~1%/iteration just above the convergence
! threshold for a hundred iterations (see the one-coefficient branch above).
res2 = sum(w*(r0-ca*u-cb*v)**2)
if (res2.gt.1.0D-2*sum(w*r0*r0) .and. abs(ca)+abs(cb).gt.10.0D0) return
if (abs(ca)+abs(cb).gt.1.0D3) return
#else
if (abs(ca)+abs(cb).gt.10.0D0) return
#endif

xacc = (1.0D0-ca-cb)*xnew + ca*h1 + cb*h2
#ifdef NGCYCLE
do i=1,nlev
   if (xacc(i).ne.xacc(i)) return            ! NaN: reject the whole step
   ! see the one-coefficient branch: clamp noise-driven negatives instead of
   ! rejecting the whole step, which would block slow-mode acceleration
   if (xacc(i).lt.0.0D0) xacc(i) = xnew(i)
enddo
#else
do i=1,nlev
   if (.not.(xacc(i).ge.0.0D0)) return  ! rejects negatives and NaNs
enddo
#endif
xnew = xacc

return
end subroutine

#ifdef NGCYCLE
subroutine ng_damposc(nlev, xnew, h1, h2)
! Per-iteration damping of oscillating (flip-flop) level populations.
! xnew is the latest plain iterate; h1 and h2 the two previous ones.
! When the latest update reverses the previous one (weighted residuals
! anti-correlated), the iteration is cycling around its fixed point, which
! lies between the last two states; some points enter genuine period-2
! limit cycles whose amplitude sits just above the 1% convergence test, so
! they would otherwise stall every stretch until the FORCECONVERGENCE
! averaging starts at level-population iteration 75. Averaging the last two
! iterates (under-relaxation with weight 1/2) turns any oscillatory mode
! with multiplier in (-3,1) into a contraction. Each point/coolant is
! treated individually, every iteration, as soon as a reversal appears;
! monotonically converging sequences (no reversal) are left untouched for
! the Ng extrapolation to accelerate.
use definitions
use healpix_types
implicit none

integer(kind=i4b), intent(in) :: nlev
real(kind=dp), intent(inout)  :: xnew(1:nlev)
real(kind=dp), intent(in)     :: h1(1:nlev), h2(1:nlev)

real(kind=dp) :: r0(1:nlev), r1(1:nlev), w(1:nlev), xacc(1:nlev)
real(kind=dp) :: floorval, maxr
integer(kind=i4b) :: i

floorval = 1.0D-10*sum(xnew)
if (floorval.le.0.0D0) return

r0 = xnew - h1
r1 = h1 - h2
do i=1,nlev
   w(i) = 1.0D0/max(xnew(i),floorval)**2
enddo

! Leave near-converged points alone (same threshold as checkconvergence)
maxr = 0.0D0
do i=1,nlev
   maxr = max(maxr, abs(r0(i))/max(xnew(i),floorval))
enddo
if (maxr.le.1.0D-2) return

! No reversal: normally contracting/transient sequence, nothing to damp
if (sum(w*r0*r1).ge.0.0D0) return

xacc = 0.5D0*(xnew + h1)
do i=1,nlev
   if (.not.(xacc(i).ge.0.0D0)) return  ! rejects negatives and NaNs
enddo
xnew = xacc

return
end subroutine
#endif

#ifdef NGRELAX
subroutine ng_relax(nlev, xnew, h1, h2, noscil)
! Adaptive per-point under-relaxation for persistently oscillating level
! populations. Fixed-weight averaging (ng_damposc, weight 1/2) only
! stabilises an oscillatory mode whose multiplier lambda satisfies
! lambda > -3; near the CO(1-0) inversion boundary (excitation temperature
! -> infinity, source function S ~ 1/(R-1) with R=g1 n0/(g0 n1) -> 1) the
! population map is violently non-linear with lambda < -3, so weight 1/2 is
! marginal and the cell flip-flops between R<1 and R>1 for ~80 iterations
! until FORCECONVERGENCE truncates the stretch. A per-point counter (noscil)
! records consecutive reversals and shrinks the relaxation weight
! omega = 1/(2+noscil); the effective multiplier 1-omega*(1-lambda) is driven
! below 1 in magnitude within a few iterations for any lambda, so the cell
! settles on the running mean of its cycle (= the true marginal fixed point
! R~=1). Non-reversing (converging/transient) sequences reset the counter and
! are left untouched for the Ng extrapolation.
! xnew = latest plain iterate, h1 = previous, h2 = the one before that.
use definitions
use healpix_types
implicit none

integer(kind=i4b), intent(in)    :: nlev
real(kind=dp), intent(inout)     :: xnew(1:nlev)
real(kind=dp), intent(in)        :: h1(1:nlev), h2(1:nlev)
integer(kind=i4b), intent(inout) :: noscil

real(kind=dp) :: r0(1:nlev), r1(1:nlev), w(1:nlev), xacc(1:nlev)
real(kind=dp) :: floorval, maxr, omega
integer(kind=i4b) :: i

floorval = 1.0D-10*sum(xnew)
if (floorval.le.0.0D0) return

r0 = xnew - h1                 ! latest update
r1 = h1 - h2                   ! previous update
do i=1,nlev
   w(i) = 1.0D0/max(xnew(i),floorval)**2
enddo

! Near-converged: nothing to do; let the oscillation counter decay
maxr = 0.0D0
do i=1,nlev
   maxr = max(maxr, abs(r0(i))/max(xnew(i),floorval))
enddo
if (maxr.le.1.0D-2) then
   noscil = 0
   return
endif

! Reversal (the latest update undoes the previous one). If the sequence is
! not reversing it is converging or in a transient: reset the counter and
! leave it for the Ng extrapolation / plain iteration.
if (sum(w*r0*r1).ge.0.0D0) then
   noscil = 0
   return
endif

! Persistent reversal: tighten the relaxation each consecutive time.
noscil = noscil + 1
omega = 1.0D0/real(2+noscil,kind=dp)
xacc = omega*xnew + (1.0D0-omega)*h1
do i=1,nlev
   if (.not.(xacc(i).ge.0.0D0)) return  ! rejects negatives and NaNs
enddo
xnew = xacc

return
end subroutine
#endif
