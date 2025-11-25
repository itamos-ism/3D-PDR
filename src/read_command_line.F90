SUBROUTINE read_command_line
!T.Bisbas, T.Bell
    use definitions
    use maincode_module, only : paramFile

    integer :: commandLineCount, iarg, argIndex
    character(len=50)::arg,param,val
    character :: delim

    delim = '='    
    commandLineCount = command_argument_count()
    
    if (commandLineCount.gt.0) then
        do iarg = 1, commandLineCount
            param = ''
            val = ''
            call get_command_argument(iarg, arg)
            arg = trim(arg)
            argIndex = index(arg,delim)
            if (argIndex.eq.0) then
                param = trim(arg)
            else
                param = trim(arg(:argIndex-1))
                val = trim(arg(argIndex+1:))
            endif
            select case (param)
            case ('-p', '--param')
                if (val.eq.'') then
                    write(6,*) "User must provide param file as -p=paramName.dat"
                    stop
                endif
                paramFile = val
                write(6,*) "Using user provided param file: ", paramFile
            case ('-h', '--help')
                write(6,*) "Allowed parameters to be passed into 3DPDR: "
                write(6,*) "============================================"
                write(6,*) "-h, --help        Write out the help menu"
                write(6,*) "-p, --params      Provide a parameter file (50 char max)"
                write(6,*) "                  e.g. -p=paramInput.dat"
                write(6,*) "============================================"
                stop
            case default
                write(6,*) "Unrecognized command line option: ", param
                write(6,*) "Use -h or --help to see allowed options"
                STOP
            end select
        end do
    else
        write(6,*) "Using default parameter file: params.dat"
        paramFile = 'params.dat'
    endif
end subroutine read_command_line