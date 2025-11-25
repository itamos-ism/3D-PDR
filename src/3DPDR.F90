program threedpdr

use definitions
use healpix_types
use healpix_module
use maincode_module
use uclpdr_module, only : start
use global_module
#ifdef OPENMP
use omp_lib
#endif
use chemistry_module
use maincode_local
!use m_IOAndVisu

!call logo

real(kind=dp)::thfpix, phfpix

call read_command_line
call readparams

#if defined RESTART && defined THERMALBALANCE
restart_file = trim(adjustl(output))//"_restart.bin"
inquire(file=restart_file, exist=restart)
if (restart) then
  write(6,*) ' '
  write(6,*) ' *****************'
  write(6,*) ' *** [RESTART] ***'
  write(6,*) ' *****************'
  write(6,*) ' '
endif
#endif

v_turb=v_turb_inp*1.0D5
start=.true.
!call writeparams
call readdensity
call allocations
call readcoolants
call read_species(nspec, species, init_abundance, mass)
call READ_RATES(NREAC,REACTANT,PRODUCT,ALPHA,BETA,GAMMA,rate,DUPLICATE,RTMIN,RTMAX)

!HEALPix setup
nside=2**level
nrays=12*nside**2
ns_max=8192
allocate(vector(1:3))
allocate(vertex(1:3,1:4))
allocate(vectors(1:3,0:nrays-1))
call mk_xy2pix 

do i=1,nrays                                                                            
  ipix=i-1 !ipix is the ID of a HEALPix ray. Runs with values 0:nrays-1
  call pix2vec_nest(nside,ipix,pix2x,pix2y,vector,vertex)                               
  vectors(1:3,ipix)=vector(1:3) !Stores in memory
enddo

call initialization
!call InitVisu

#ifdef OPENMP
!$OMP PARALLEL
!$OMP MASTER
CPUs = OMP_GET_NUM_THREADS()
write(6,*) "Number of CPUs: ", OMP_GET_NUM_THREADS()
!$OMP END MASTER
!$OMP END PARALLEL
#else
write(6,*) 'Processing in SERIAL'
#endif

write(6,*) 'Performing ray-tracing'
call evaluation_points

write(6,*) 'Calculating UV field'
call calc_UVfield

out_file = trim(adjustl(output))//".UV.fin"
out_file2 = trim(adjustl(out_file))//"]"
write(6,'(" Writing file [",A)') trim(adjustl(out_file2))
open(unit=21,file=out_file,status='replace')
do p=1,pdr_ptot
  write(21,*) pdr(p)%x, pdr(p)%y, pdr(p)%z, pdr(p)%rho, pdr(p)%UVfield
enddo
close(21)

out_file = trim(adjustl(output))//".rayDir.fin"
out_file2 = trim(adjustl(out_file))//"]"
write(6,'(" Writing file [",A)') trim(adjustl(out_file2))
open(unit=21,file=out_file,status='replace')
do p=0,nrays-1
  call pix2ang_nest(nside, p, thfpix, phfpix)
  write(21, *) thfpix, phfpix
enddo
close(21)

out_file = trim(adjustl(output))//".rayAV.fin"
out_file2 = trim(adjustl(out_file))//"]"
write(6,'(" Writing file [",A)') trim(adjustl(out_file2))
open(unit=21,file=out_file,status='replace')

do p=1,pdr_ptot
  write(21, *) pdr(p)%AV
enddo
close(21)
#ifdef UVDEBUG
stop
#endif

#ifdef CRATTENUATION
write(6,*) 'Calculating CRIR attenuation'
#else
write(6,*) 'Setting CR ionization rate'
#endif
call calc_CRattenuation

write(6,*) 'Setting dust temperatures'
#ifdef DUST
call calculate_dust_temperatures
#else
do p=1,pdr_ptot
  pdr(p)%Tdust=dust_temperature
enddo
#endif

#ifdef THERMALBALANCE
first_time=.true.
#endif
write(6,*) 'Calculating column densities'
referee=0
call calc_columndens

referee=1

start_time = 0.0D0

#if defined RESTART && defined THERMALBALANCE
if (restart) then
   open(unit=77,file=restart_file,status='old',form='unformatted')
   do p=1,pdr_ptot
     read(77) pdr(p)%restconverged,pdr(p)%Tgas,pdr(p)%abundance
     pdr(p)%nTgas = pdr(p)%Tgas
     pdr(p)%Tlow = pdr(p)%Tgas/2.0
     pdr(p)%Thigh = pdr(p)%Tgas*1.5
   enddo
endif
#endif

ITERATION = 0
!======== LTE LEVEL POPULATIONS ============
write(6,*) ''; write(6,*) 'Calculating LTE level populations' 


call chemicaliterations(1,CHEMITERATIONS)

do p=1,pdr_ptot
  do i=1,coo
    pdr(p)%coolant(i)%isconverged = .false.
  enddo
  pdr(p)%levelconverged = .false.
#ifdef THERMALBALANCE
  pdr(p)%fullyconverged = .false.
#ifndef ONEDIMENSIONAL
  if (pdr(p)%rho<1.0) pdr(p)%fullyconverged = .true.
#endif
  pdr(p)%doleveltmin    = .false.
#endif
  do i=1,coo
    allocate(pdr(p)%cooling(i))
  enddo
  do i=1,12
    allocate(pdr(p)%heating(i))
  enddo
enddo

write(6,*) ''
write(6,*) '----------------'
write(6,*) 'Iterations begin'
write(6,*) '----------------'

levpop_iteration=0
do k=1,coo
  coolant(k)%percentage=0
enddo
  
DO ITERATION=1,ITERTOT
        write_output=.false.
        write(6,*) ''
        write(6,'(" Iteration = ",I4)') iteration
        
        levpop_iteration=levpop_iteration+1
        write(6,'(" Level population iteration = ",I3)') levpop_iteration

#if defined RESTART && defined THERMALBALANCE
        if (level_conv.and.first_time) then
            do p=1,pdr_ptot
               if (pdr(p)%restconverged) pdr(p)%fullyconverged = .true.
            enddo
        endif
#endif

        !LTE/LVG computations
        IF (iteration.gt.1.and.levpop_iteration.eq.1) THEN
               call chemicaliterations(2,3)
        ELSE  
               call coolingfunctions
        ENDIF 
   
        call changetemperature

#ifdef THERMALBALANCE
        if (level_conv.and.first_time) first_time=.false.
        i=0
        do p=1,pdr_ptot
           if (pdr(p)%fullyconverged) i=i+1
        enddo
        if (level_conv) then
           write(6,*) 'Resetting [level_conv=.false.]'
           level_conv=.false.
           write_output=.true.
        endif
        thermal_percentage = 100.D0*real(i,kind=dp)/real(pdr_ptot,kind=dp)
        write(*,'(" Thermal balance is ",F5.1,"% converged.")') thermal_percentage
        !write(*,'(" [",I6,"/",I6,"]")') i,pdr_ptot
        if (i.eq.pdr_ptot) then
             write(6,*) '#### Converged through thermal balance ####'
             goto 2
        endif
#endif

        RELCH_conv=.true.
        do p=1,pdr_ptot
          pdr(p)%levelconverged = .false.
        enddo
        
        call checkconvergence

        if (.not.relch_conv) then 
             goto 1
        else
             write(6,*) '#### Converged through level populations ####'
             levpop_iteration=0
#ifdef THERMALBALANCE
#ifdef RESTART
             close(77);open(unit=77,file=restart_file,status='replace',form='unformatted')
             out_file2 = trim(adjustl(restart_file))//"]"
             write(6,'(" Writing file [",A)') trim(adjustl(out_file2))
             do p=1,pdr_ptot
                 write(77) pdr(p)%fullyconverged,pdr(p)%Tgas,pdr(p)%abundance
             enddo
#endif
             write(6,*) 'Enabling thermal balance routine in next iteration'
             level_conv=.true.
             goto 1
#else
             goto 2
#endif
        endif

1 continue 

        i=0
        do k=1,coo
           coolant(k)%incr = 0
        enddo
        do p=1,pdr_ptot
           if (pdr(p)%levelconverged) i=i+1
           do k=1,coo
             if (pdr(p)%coolant(k)%isconverged) coolant(k)%incr = coolant(k)%incr + 1
           enddo
        enddo
        levpop_percentage  = 100.D0*real(i,kind=dp)/real(pdr_ptot,kind=dp)
        do k=1,coo
           coolant(k)%percentage = 100.D0*real(coolant(k)%incr)/real(pdr_ptot)
           write(*,'(A, F6.2, A)') ' Coolant '//trim(adjustl(coolant(k)%cname))//' is ', coolant(k)%percentage, '% converged'
         !  write(*,*) trim(adjustl(coolant(k)%cname)),' %convergence: ',coolant(k)%percentage
        enddo
        write(*,'(" Level populations are ",F5.1,"% converged.")') levpop_percentage
        !write(*,'(" [",I6,"/",I6,"]")') i,pdr_ptot
          
#ifdef THERMALBALANCE
        if (int(levpop_percentage,kind=i4b).ge.100) then
           do p=1,pdr_ptot
             pdr(p)%levelconverged = .false.
           enddo
           write(6,*) 'Resetting [level_converged=.false.] array'
        endif
#endif
        !write(6,*) 'Updating population densities...'
#ifdef OPENMP
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(p,ilevel)
#endif
       do p=1,pdr_ptot
         do k=1,coo
           DO ilevel=1,coolant(k)%cnlev
             pdr(p)%coolant(k)%pop(ilevel) = pdr(p)%coolant(k)%solution(ilevel)
           ENDDO
         enddo
       enddo
#ifdef OPENMP
!$OMP END PARALLEL DO
#endif

END DO !ITERATIONS

2  continue
if (iteration.ge.1) then

   if (iteration.lt.itertot) then
     WRITE(6,*) '3DPDR converged after ',ITERATION-1,' iterations'
   else
     write(6,*) 'Reached maximum number of iterations without convergence.'
     write(6,*) 'To reach convergence, increase the relative number in [params.dat]'
   endif
   write(6,*) 'Writing final outputs'
   
   call writeoutputs

endif

!#ifdef RAYTHEIA
!call dump_visu
!#endif

write(6,*) ''
#ifdef OPENMP
write(6,*) "Number of CPUs: ",CPUs
#endif
write(6,*) 'Finished !'

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~~~~~~~~END OF MAIN PROGRAM~~~~~~~~~~~~~~~~
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

end Program
