!Written by Thomas G. Bisbas
program RTtool
use RTtool_module

call readparams
call readcoolants
call initializations
call readpdrfiles

call readspopfile
call excitation
call opticaldepth
call radiativetransfer
call writeoutputs

end program
