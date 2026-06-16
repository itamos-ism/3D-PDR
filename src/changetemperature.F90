subroutine changetemperature
use maincode_module
use maincode_local
use global_module

#ifdef THERMALBALANCE
! The heating evaluated below is only consumed by the thermal-balance update
! (and the final output), so skip the whole rates+heating sweep on iterations
! where no temperature update can occur.
if (.not.level_conv) return
#endif

#ifdef OPENMP
#ifdef THERMALBALANCE
#ifdef ILLINOIS
!$OMP PARALLEL DO DEFAULT(SHARED) SCHEDULE(GUIDED) PRIVATE(p,allheating,temp_Tgas,xlo,xhi,xnew,denom,rate) &
#else
!$OMP PARALLEL DO DEFAULT(SHARED) SCHEDULE(GUIDED) PRIVATE(p,allheating,temp_Tgas,rate) &
#endif
#else
!$OMP PARALLEL DO DEFAULT(SHARED) SCHEDULE(GUIDED) PRIVATE(p,allheating,rate)&
#endif
!$OMP PRIVATE(NRGR,NRH2,NRHD,NRCO,NRCI,NRSI)
#endif
        do p=1,pdr_ptot
#ifdef THERMALBALANCE
           ! Skip this pdrpoint if the temperature has already converged
           if (pdr(p)%fullyconverged) cycle
#endif
           if (allocated(allheating)) deallocate(allheating); allocate(allheating(1:12))
           CALL CALCULATE_REACTION_RATES(pdr(p)%nTgas,pdr(p)%Tdust,nrays,pdr(p)%rad_surface(0:nrays-1),pdr(p)%AV(0:nrays-1),&
                  &pdr(p)%column_NH2(0:nrays-1),pdr(p)%column_NHD(0:nrays-1),pdr(p)%column_NCO(0:nrays-1),&
                  &pdr(p)%column_NC(0:nrays-1),pdr(p)%column_NS(0:nrays-1),&
                  &nreac, reactant, product, alpha, beta, gamma, rate, rtmin, rtmax, duplicate, nspec,&
                  &NRGR,NRH2,NRHD,NRCO,NRCI,NRSI,pdr(p)%abundance(NELECT)*pdr(p)%rho,pdr(p)%rho,pdr(p)%zetalocal)
           call calc_heating(pdr(p)%rho,pdr(p)%nTgas,pdr(p)%Tdust,pdr(p)%UVfield, &
                  &v_turb,nspec,pdr(p)%abundance(:),nreac,rate,allheating,&
                  &NRGR,NRH2,NRHD,NRCO,NRCI,NRSI,pdr(p)%zetalocal)

           pdr(p)%heating=allheating

#ifdef THERMALBALANCE
           ! Calculate the difference between the total heating and total cooling rates (Fmean)
           ! and the absolute value of the relative difference between the two rates (Fratio)
           pdr(p)%Fmean = pdr(p)%heating(12) - sum(pdr(p)%cooling(:))
           pdr(p)%Fratio = 2.0D0*abs(pdr(p)%Fmean)/abs(pdr(p)%heating(12) + sum(pdr(p)%cooling(:)))

           ! Store the current temperature in a dummy variable
           ! Do not start testing the thermal balance until enough iterations have passed for the level populations to begin to converge...
           if (level_conv.and.first_time) then
                   temp_Tgas = pdr(p)%nTgas
                   ! Determine the temperature bracket to begin searching within...
                   if (pdr(p)%Fmean.eq.0) then ! Handle the (very rare) case when the initial guess temperature is the correct value
                         pdr(p)%Tlow = pdr(p)%nTgas  ! Update the value of Tlow
                         pdr(p)%Thigh = pdr(p)%nTgas ! Update the value of Thigh
#ifdef ILLINOIS
                         pdr(p)%Flow = pdr(p)%Fmean
                         pdr(p)%Fhigh = pdr(p)%Fmean
#endif
                   else if (pdr(p)%Fmean.gt.0) then !---> HEATING
                         pdr(p)%Tlow = pdr(p)%nTgas  ! Update the value of Tlow
#ifdef ILLINOIS
                         pdr(p)%Flow = pdr(p)%Fmean  ! Record F at this just-evaluated point
#endif
                         pdr(p)%nTgas = 1.3D0*pdr(p)%nTgas !increase 30%
                         pdr(p)%previouschange = "H" !we increased
                   else if (pdr(p)%Fmean.lt.0) then !---> COOLING
                         pdr(p)%Thigh = pdr(p)%nTgas ! Update the value of Thigh
#ifdef ILLINOIS
                         pdr(p)%Fhigh = pdr(p)%Fmean ! Record F at this just-evaluated point
#endif
                         pdr(p)%nTgas = 0.7D0*pdr(p)%nTgas !decrease 30%
                         if (pdr(p)%nTgas.lt.Tmin) pdr(p)%nTgas=Tmin
                         pdr(p)%previouschange = "C" !we decreased
                   endif !(Fmean.eq.0)
                   pdr(p)%Tgas = temp_Tgas
                   if (pdr(p)%Tgas.lt.Tmin) pdr(p)%Tgas=Tmin
           else if (level_conv.and..not.first_time) then
                   temp_Tgas = pdr(p)%nTgas
                   ! Check for convergence in both the heating-cooling imbalance and the temperature difference between iterations
                   if (pdr(p)%Fratio.le.Fcrit) pdr(p)%fullyconverged = .true.
                   if (.not.pdr(p)%dobinarychop) then
                         !if we *still* need to heat, increase by 30% 
                         if (pdr(p)%Fmean.gt.0.and.pdr(p)%previouschange.eq."H") then
                               pdr(p)%Tlow = pdr(p)%nTgas
#ifdef ILLINOIS
                               pdr(p)%Flow = pdr(p)%Fmean
#endif
                               pdr(p)%nTgas = 1.3D0*pdr(p)%nTgas
                               pdr(p)%Thigh = pdr(p)%nTgas
                               pdr(p)%previouschange = "H"
                         endif
                         !if we *still* need to cool, decrease by 30%
                         if (pdr(p)%Fmean.lt.0.and.pdr(p)%previouschange.eq."C") then
                               pdr(p)%Thigh = pdr(p)%nTgas
#ifdef ILLINOIS
                               pdr(p)%Fhigh = pdr(p)%Fmean
#endif
                               pdr(p)%nTgas = 0.7D0*pdr(p)%nTgas
                               pdr(p)%Tlow = pdr(p)%nTgas
                               pdr(p)%previouschange = "C"
                               if (pdr(p)%nTgas.lt.Tmin) then
                                   pdr(p)%nTgas=Tmin
                                   pdr(p)%Tlow=Tmin
                               endif
                         endif !(Fmean.lt.0.and.pdr(p)%previouschange.eq."C")
                        !For all other cases do binary chop and flag the process as .true.
                        !Needs heating but previously it was decreased by 30%. 
                        if (pdr(p)%Fmean.gt.0.and.pdr(p)%previouschange.eq."C") then
#ifdef ILLINOIS
                               pdr(p)%Tlow = pdr(p)%nTgas
                               pdr(p)%Flow = pdr(p)%Fmean
#endif
                               pdr(p)%nTgas = (pdr(p)%Thigh + pdr(p)%Tlow)/2.0D0
                               pdr(p)%dobinarychop=.true.  !from now on
#ifdef ILLINOIS
                               pdr(p)%ill_last = "N"
#endif
                        endif
                        !Needs cooling but previously it was increased by 30%. Now do binary chop
                        if (pdr(p)%Fmean.lt.0.and.pdr(p)%previouschange.eq."H") then
#ifdef ILLINOIS
                               pdr(p)%Thigh = pdr(p)%nTgas
                               pdr(p)%Fhigh = pdr(p)%Fmean
#endif
                               pdr(p)%nTgas = (pdr(p)%Thigh + pdr(p)%Tlow)/2.0D0
                               pdr(p)%dobinarychop=.true.  !from now on
#ifdef ILLINOIS
                               pdr(p)%ill_last = "N"
#endif
                        endif
#ifdef ILLINOIS
                   else !from now on Illinois-modified regula-falsi in ln(T) (bracket already established)
                        if (pdr(p)%Fmean.gt.0) then
                               pdr(p)%Tlow = pdr(p)%nTgas
                               pdr(p)%Flow  = pdr(p)%Fmean
                               if (pdr(p)%ill_last.eq."L") pdr(p)%Fhigh = 0.5D0*pdr(p)%Fhigh
                               pdr(p)%ill_last = "L"
                        endif
                        if (pdr(p)%Fmean.lt.0) then
                               pdr(p)%Thigh = pdr(p)%nTgas
                               pdr(p)%Fhigh = pdr(p)%Fmean
                               if (pdr(p)%ill_last.eq."H") pdr(p)%Flow = 0.5D0*pdr(p)%Flow
                               pdr(p)%ill_last = "H"
                        endif
                        ! Fmean==0 already forced fullyconverged=.true. via Fratio<=Fcrit above; nothing to do.

                        if (pdr(p)%Fmean.ne.0) then
                              xlo   = log(pdr(p)%Tlow)
                              xhi   = log(pdr(p)%Thigh)
                              denom = pdr(p)%Fhigh - pdr(p)%Flow

                              if (abs(denom).lt.1.0D-300 .or. abs(xhi-xlo).lt.1.0D-12) then
                                    xnew = 0.5D0*(xlo + xhi)
                              else
                                    xnew = xhi - pdr(p)%Fhigh*(xhi-xlo)/denom
                                    ! keep strictly inside the bracket (avoid edge-sticking)
                                    if (xnew.le.min(xlo,xhi).or.xnew.ge.max(xlo,xhi)) then
                                          xnew = 0.5D0*(xlo + xhi)
                                    endif
                                    ! safeguard against landing too close to an endpoint
                                    if ((xnew-min(xlo,xhi)).lt.1.0D-3*abs(xhi-xlo) .or. &
                                       &(max(xlo,xhi)-xnew).lt.1.0D-3*abs(xhi-xlo)) then
                                          xnew = 0.5D0*(xlo + xhi)
                                    endif
                              endif

                              pdr(p)%nTgas = exp(xnew)
                        endif
#else
                   else !from now on only binary chop (we found the min-max by the 30% increase/decrease)
                        if (pdr(p)%Fmean.gt.0) then
                               pdr(p)%Tlow = pdr(p)%nTgas
                               pdr(p)%nTgas = (pdr(p)%nTgas + pdr(p)%Thigh) / 2.0D0
                        endif
                        if (pdr(p)%Fmean.lt.0) then
                               pdr(p)%Thigh = pdr(p)%nTgas
                               pdr(p)%nTgas = (pdr(p)%nTgas + pdr(p)%Tlow) / 2.0D0
                        endif
#endif
                   endif !dobinarychop

                   ! If the temperatures are converging in a value that thermal balance can't be reached
                   ! double the high temperature and half the low one. If the temperatures are still not converging, force convergence.
!                   if ((abs(pdr(p)%nTgas-pdr(p)%Tgas).le.Tdiff).and.(pdr(p)%Fratio.gt.Fcrit)) pdr(p)%fullyconverged=.true.
#ifdef REBRACKET
                   ! The bracket can become confined away from the true root when the
                   ! early F evaluations used level populations that were still relaxing
                   ! (their cooling was over/under-estimated). Instead of force-converging
                   ! a stalled point with a large residual imbalance, re-open the search
                   ! on the side the imbalance points to and resume the 30% expansion
                   ! phase; the next evaluations use the now-converged populations. After
                   ! 3 attempts force convergence as before (genuinely rootless points).
                   if ((abs(pdr(p)%nTgas-temp_Tgas).le.Tdiff).and.(pdr(p)%Fratio.gt.Fcrit)) then
                         if (pdr(p)%nrebracket.lt.3 .and. temp_Tgas.ge.Tmin .and. temp_Tgas.le.Tmax) then
                               pdr(p)%nrebracket = pdr(p)%nrebracket + 1
                               pdr(p)%dobinarychop = .false.
#ifdef ILLINOIS
                               pdr(p)%ill_last = "N"
#endif
                               if (pdr(p)%Fmean.gt.0) then !too cold: search upwards again
                                     pdr(p)%Tlow = temp_Tgas
#ifdef ILLINOIS
                                     pdr(p)%Flow = pdr(p)%Fmean
#endif
                                     pdr(p)%nTgas = 1.3D0*temp_Tgas
                                     pdr(p)%Thigh = pdr(p)%nTgas
                                     pdr(p)%previouschange = "H"
                               else                        !too hot: search downwards again
                                     pdr(p)%Thigh = temp_Tgas
#ifdef ILLINOIS
                                     pdr(p)%Fhigh = pdr(p)%Fmean
#endif
                                     pdr(p)%nTgas = 0.7D0*temp_Tgas
                                     pdr(p)%Tlow = pdr(p)%nTgas
                                     pdr(p)%previouschange = "C"
                                     if (pdr(p)%nTgas.lt.Tmin) then
                                         pdr(p)%nTgas=Tmin
                                         pdr(p)%Tlow=Tmin
                                     endif
                               endif
                         else
                               pdr(p)%fullyconverged=.true.
                         endif
                   endif
#else
                   if ((abs(pdr(p)%nTgas-temp_Tgas).le.Tdiff).and.(pdr(p)%Fratio.gt.Fcrit)) pdr(p)%fullyconverged=.true.
#endif

                   ! Replace the previous temperature with the current value
                   pdr(p)%Tgas = temp_Tgas

                   if ((temp_Tgas.lt.Tmin).and.(pdr(p)%Fmean.lt.0)) pdr(p)%fullyconverged=.true.
                   if ((temp_Tgas.gt.Tmax).and.(pdr(p)%Fmean.gt.0)) pdr(p)%fullyconverged=.true.
 
                   if (pdr(p)%fullyconverged) then
                         if (temp_Tgas.lt.Tmin) then 
                               pdr(p)%Tgas = Tmin
                               pdr(p)%nTgas = Tmin
                               if (pdr(p)%doleveltmin) then 
                                    pdr(p)%fullyconverged=.true.
                               else
                                    pdr(p)%fullyconverged=.false.
                                    pdr(p)%levelconverged=.false.
                                    temp_Tgas=Tmin
                               endif
                               pdr(p)%doleveltmin=.true.
                         endif !(temp_Tgas.lt.Tmin)
                         if (temp_Tgas.gt.Tmax) then
                              pdr(p)%Tgas = Tmax
                              pdr(p)%nTgas = Tmax
                         endif !(temp_Tgas.gt.Tmax)
                   endif !(pdr(p)%fullyconverged)
           endif !(level_conv.and.first_time)
#endif

        enddo !1,pdr_ptot
#ifdef OPENMP
!$OMP END PARALLEL DO
#endif

return
end subroutine
