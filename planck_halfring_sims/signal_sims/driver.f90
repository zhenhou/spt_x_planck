subroutine mpi_assign(njobs_total, npc, rank, istart, iend)

    use healpix_types
    
    integer(I4B),   intent(in)  :: njobs_total, npc, rank
    integer(I4B),   intent(out) :: istart, iend

    integer(I4B)                :: nsims_tmp

    if (mod(njobs_total, npc) .eq. 0) then 
        nsims_tmp = njobs_total/npc
    else
        nsims_tmp = njobs_total/npc + 1
    endif
    
    istart = rank * nsims_tmp
    iend = min((rank+1)*nsims_tmp-1, njobs_total-1)

    return
end subroutine

program sxp_sims
    
    use healpix_types
    use extension,      only: nArguments, getArgument
    use sxp_signal_sims_utils

    use mpi
    
    implicit none

    character(256)          :: ini_file
    character(6)            :: s_str

    integer(I4B)            :: nsims, isim
    integer(I4B)            :: istart, iend
    integer(I4B)            :: ierr, npc, rank
    integer(I4B)            :: run_status

    if (nArguments() .eq. 1) then
        call getArgument(1, ini_file)
    else
        write(*,*) "./sxp_sims param.ini"
        stop
    endif

    call params_init(ini_file)
    call input_data_init()

    nsims = params%nsims

    call MPI_INIT(IERR)
    if (ierr .ne. MPI_SUCCESS) then
        print*, "MPI initialization error"
        stop
    endif
    call MPI_COMM_SIZE(MPI_COMM_WORLD, NPC, IERR)
    call MPI_COMM_RANK(MPI_COMM_WORLD, RANK, IERR)
    call mpi_assign(nsims, npc, rank, istart, iend)

    print*, rank, istart, iend

    do isim=istart, iend 
        call check_status(isim, run_status)
        
        if (run_status .eq. 0) then
            write(s_str,'(I6)') isim
            s_str = adjustl(s_str)

            write(*,*) "start isim = "//trim(s_str)
            call do_sxp_sims(isim)
            call write_prj(isim)
        endif
    enddo
end program
