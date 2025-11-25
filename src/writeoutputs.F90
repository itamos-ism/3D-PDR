subroutine writeoutputs
!for ASCII outputs
use healpix_types
use maincode_module
use global_module

#if (defined(OUTRAYINFO) &&  !defined(ONEDIMENSIONAL))
character(len=100) :: outFormat
#endif

#ifdef PYWRAP
out_file = trim(adjustl(output))//".pywrap"
out_file2 = trim(adjustl(out_file))//"]"
write(6,'(" Writing file [",A)') trim(adjustl(out_file2))
open(unit=21,file=out_file,status='replace')

write(21,'("x Av Tgas Tdust nH fuv")',advance='no')
do i=1,nspec
  write(21,'(A6)',advance='no') trim(species(i))
enddo
write(21,'(" dummy1 photoelectric dummy2 Cionization H2formation H2photodis FUVpumping")',advance='no')
write(21,'(" cosmicray turbulent chemical gasgrain totheat ")',advance='no')
do i=1,coo
   write(21,'(A10)',advance='no') trim(adjustl(coolant(i)%cname))//"cool"
enddo
write(21,'(" totcool")',advance='no')
write(21,*) ''

do p=1,pdr_ptot
     write(21,'(300ES15.7)') pdr(p)%x, pdr(p)%AV(6), pdr(p)%Tgas,pdr(p)%Tdust,pdr(p)%rho,pdr(p)%UVfield,&
        &pdr(p)%abundance,pdr(p)%heating,pdr(p)%cooling,pdr(p)%totalcooling
enddo

close(21)
#endif


#ifdef CHEMANALYSIS
!-------------------------------------
!OUTPUT FOR CHEMICAL ANALYSIS
!-------------------------------------
out_file = trim(adjustl(output))//".rates.fin"
out_file2 = trim(adjustl(out_file))//"]"
write(6,'(" Writing file [",A)') trim(adjustl(out_file2))
open(unit=98,file=out_file,status='replace')
do p=1,pdr_ptot
   call analyse_chemistry(p, end_time, pdr(p)%rho, pdr(p)%Tgas, &
     &12, pdr(p)%AV(6), nspec, species,pdr(p)%abundance(1:nspec),nreac, reactant, &
     & product, temp_rate(:,p))
enddo
close(98)
#endif


!-------------------------------------
!OUTPUT FOR PDR PARAMETERS
!-------------------------------------
   out_file = trim(adjustl(output))//".params"
   out_file2 = trim(adjustl(out_file))//"]"
   write(6,'(" Writing file [",A)') trim(adjustl(out_file2))
   open(unit=21,file=out_file,status='replace')
   write(21,'(3ES11.3)') Gext(1), zeta*1.3d-17, metallicity
   close(21)


!-------------------------------------
!OUTPUT FOR ABUNDANCES AND TEMPERATURE
!-------------------------------------
   out_file = trim(adjustl(output))//".pdr.fin"
   out_file2 = trim(adjustl(out_file))//"]"
   write(6,'(" Writing file [",A)') trim(adjustl(out_file2))
   open(unit=21,file=out_file,status='replace')

  do p=1,pdr_ptot
#ifdef ONEDIMENSIONAL
     write(21,'(I7,4ES15.7,I5,300ES15.7)') p,pdr(p)%x, pdr(p)%AV(6), pdr(p)%Tgas,pdr(p)%Tdust,pdr(p)%etype,&
     &pdr(p)%rho,pdr(p)%UVfield,pdr(p)%abundance
#else
     write(21,'(I7,5ES15.7,I5,300ES15.7)') p,pdr(p)%x, pdr(p)%y, pdr(p)%z, pdr(p)%Tgas,pdr(p)%Tdust,&
     &pdr(p)%etype,pdr(p)%rho,pdr(p)%UVfield,pdr(p)%abundance!,pdr(p)%AV
#endif
  enddo

close(21)
!-----------------------------------------
!END OUTPUT FOR ABUNDANCES AND TEMPERATURE
!-----------------------------------------

#ifdef CRATTENUATION
!-------------------------
!OUTPUT FOR LOCAL CR VALUE
!-------------------------
   out_file = trim(adjustl(output))//".cr.fin"
   out_file2 = trim(adjustl(out_file))//"]"
   write(6,'(" Writing file [",A)') trim(adjustl(out_file2))
   open(unit=21,file=out_file,status='replace')

  do p=1,pdr_ptot
     write(21,'(ES15.7)') pdr(p)%zetalocal*1.3d-17
  enddo

close(21)
!-----------------------------
!END OUTPUT FOR LOCAL CR VALUE
!-----------------------------
#endif

#if (defined(OUTRAYINFO) &&  !defined(ONEDIMENSIONAL))
  write(outFormat, '( "(", I6, "ES15.7)" )') nrays

  out_file = trim(adjustl(output))//".rayAV.fin"
  out_file2 = trim(adjustl(out_file))//"]"
  write(6,'(" Writing file [",A)') trim(adjustl(out_file2))
  open(unit=21,file=out_file,status='replace')

  do p=1,pdr_ptot
    write(21, outFormat) pdr(p)%AV
  enddo

  out_file = trim(adjustl(output))//".rayH2.fin"
  out_file2 = trim(adjustl(out_file))//"]"
  write(6,'(" Writing file [",A)') trim(adjustl(out_file2))
  open(unit=21,file=out_file,status='replace')

  do p=1,pdr_ptot
    write(21, outFormat) pdr(p)%column_NH2
  enddo

  out_file = trim(adjustl(output))//".rayHD.fin"
  out_file2 = trim(adjustl(out_file))//"]"
  write(6,'(" Writing file [",A)') trim(adjustl(out_file2))
  open(unit=21,file=out_file,status='replace')

  do p=1,pdr_ptot
    write(21, outFormat) pdr(p)%column_NHD
  enddo
#endif


!---------------------------
!OUTPUT FOR COOLING FUNCTION
!---------------------------
   out_file = trim(adjustl(output))//trim(adjustl(".cool"))//".fin"
   out_file2 = trim(adjustl(out_file))//"]"
   write(6,'(" Writing file [",A)') trim(adjustl(out_file2))
   open(unit=13,file=out_file,status='replace')

   do p=1,pdr_ptot
#ifdef ONEDIMENSIONAL
      write(13,'(I7,200ES15.7)') p, pdr(p)%x, pdr(p)%AV(6), pdr(p)%cooling(:),pdr(p)%totalcooling
     
#else
      write(13,'(I7,200ES15.7)') p, pdr(p)%x, pdr(p)%y, pdr(p)%z, pdr(p)%cooling(:),pdr(p)%totalcooling!, pdr(p)%AV(:)
#endif
   enddo
!-------------------------------
!END OUTPUT FOR COOLING FUNCTION
!-------------------------------


!---------------------------
!OUTPUT FOR HEATING FUNCTION
!---------------------------
   out_file = trim(adjustl(output))//trim(adjustl(".heat"))//".fin"
   out_file2 = trim(adjustl(out_file))//"]"
   write(6,'(" Writing file [",A)') trim(adjustl(out_file2))
   open(unit=14,file=out_file,status='replace')


   do p=1,pdr_ptot
#ifdef ONEDIMENSIONAL
      write(14,'(I7,200ES15.7)') p, pdr(p)%x, pdr(p)%AV(6), pdr(p)%heating
#else
      write(14,'(I7,200ES15.7)') p, pdr(p)%x, pdr(p)%y, pdr(p)%z, pdr(p)%heating!, pdr(p)%AV
#endif
   enddo
   close(14)
!-------------------------------
!END OUTPUT FOR HEATING FUNCTION
!-------------------------------


#ifdef ONEDIMENSIONAL
!-----------------------
!OUTPUT FOR EMISSIVITIES
!-----------------------
 do k=1,coo
   out_file = trim(adjustl(output))//"."//trim(adjustl(coolant(k)%cname))//trim(adjustl(".line"))//".fin"
   out_file2 = trim(adjustl(out_file))//"]"
   write(6,'(" Writing file [",A)') trim(adjustl(out_file2))
   open(unit=16,file=out_file,status='replace')

   !experimental [1D models only]
   do p = 1, pdr_ptot
     write(16, '(I5, 1X, ES15.7, 1X, ES15.7, 1X)', advance='no') p,pdr(p)%x, pdr(p)%AV(6)
     do ilevel=1,coolant(k)%cnlev-1
       if (pdr(p)%coolant(k)%line(ilevel+1,ilevel).lt.1d-99) pdr(p)%coolant(k)%line(ilevel+1,ilevel)=0.0D0
       write(16, '(100ES15.7)', advance='no') pdr(p)%coolant(k)%line(ilevel+1,ilevel)
     end do
     write(16, *)
   end do
   close(16)
 enddo
!---------------------------
!END OUTPUT FOR EMISSIVITIES
!---------------------------

!----------------------------
!OUTPUT FOR LEVEL POPULATIONS
!----------------------------
 do k=1,coo
   out_file = trim(adjustl(output))//"."//trim(adjustl(coolant(k)%cname))//trim(adjustl(".spop"))//".fin"
   out_file2 = trim(adjustl(out_file))//"]"
   write(6,'(" Writing file [",A)') trim(adjustl(out_file2))
   open(unit=16,file=out_file,status='replace')

   !experimental [1D models only]
   do p = 1, pdr_ptot
     write(16, '(I5, 1X, ES15.7, 1X, ES15.7, 1X)', advance='no') p,pdr(p)%x, pdr(p)%AV(6)
     do ilevel=1,coolant(k)%cnlev
       if (pdr(p)%coolant(k)%pop(ilevel).lt.1d-99) pdr(p)%coolant(k)%pop(ilevel) = 0.0D0
     end do
     write(16, '(100ES15.7)', advance='no') pdr(p)%coolant(k)%pop
     write(16, *)
   end do
   close(16)
 enddo
!--------------------------------
!END OUTPUT FOR LEVEL POPULATIONS
!--------------------------------
#endif


!----------------------------
!OUTPUT FOR RTtool
!----------------------------
   out_file = trim(adjustl(output))//trim(adjustl(".RTspop"))//".fin"
   out_file2 = trim(adjustl(out_file))//"]"
   write(6,'(" Writing file [",A)') trim(adjustl(out_file2))
   open(unit=16,file=out_file,status='replace')

#ifdef REDUCED
   write(16,*) 'REDUCED'
#elif MEDIUM
   write(16,*) 'MEDIUM'
#elif FULL
   write(16,*) 'FULL'
#endif
   do k=1,coo
     write(16,*) coolfile(k)
   enddo
   write(16,*) 'ENDCOOLFILES'
   do p = 1, pdr_ptot
     do k = 1, coo
       do ilevel=1,coolant(k)%cnlev
         if (pdr(p)%coolant(k)%pop(ilevel).lt.1d-99) pdr(p)%coolant(k)%pop(ilevel) = 0.0D0
       end do
     end do
     do k = 1, coo
       write(16, '(100ES15.7)', advance='no') pdr(p)%coolant(k)%pop(1:coolant(k)%cnlev)
     end do
     write(16, *)
   end do

   close(16)
!--------------------------------
!END OUTPUT FOR LEVEL POPULATIONS
!--------------------------------


return
end subroutine

