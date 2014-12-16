module sxp_signal_sims_utils
    
    use healpix_types
    use paramfile_io
    use pix_tools,      only: nside2npix
    use alm_tools,      only: alm2map, alter_alm, pixel_window
    use misc_utils,     only: string
    use rngmod

    use healpix_taylor
    
    implicit none

    integer(I4B),   parameter   :: num_fields = 20
    character(50)               :: fields_name(0:num_fields-1)

    type ini_params
        logical         :: init = .false.

        integer(I4B)    :: nsims
        integer(I4B)    :: nside, lmax
        integer(I4B)    :: lmax_input_alms = 8000
        integer(I4B)    :: max_order = 4

        character(512)  :: input_alms_path, input_alms_root

        character(512)  :: beam_file
        logical         :: apply_pixel_window
        logical         :: vary_beam

        character(512)  :: input_gal_theta_path, input_gal_theta_root
        character(512)  :: input_gal_phi_path, input_gal_phi_root

        character(512)  :: output_path, output_root

    end type ini_params

    type data_params
        integer(I4B),   allocatable :: num_pixels(:)
        real(DP),       allocatable :: gal_theta(:), gal_phi(:)
        real(DP),       allocatable :: beam_data(:,:)
        real(DP),       allocatable :: wpix(:,:)
    end type data_params

    type(ini_params)            :: params
    type(data_params),  target  :: input_data

    real(DP),       allocatable,    target  :: beam_total(:,:)
    real(DP),       allocatable,    target  :: map_sim(:,:)
    complex(SPC),   allocatable             :: alms_read(:,:,:)
    complex(DPC),   allocatable             :: alms_sim(:,:,:)


    contains

    subroutine fields_info_init()

        integer(I4B)    :: i_field

        fields_name(0)   = 'ra5h30dec-55_2008'
        fields_name(1)   = 'ra23h30dec-55_2008'
        fields_name(2)   = 'ra23h30dec-55_2010'
        fields_name(3)   = 'ra21hdec-60'
        fields_name(4)   = 'ra3h30dec-60'
        fields_name(5)   = 'ra21hdec-50'
        fields_name(6)   = 'ra4h10dec-50'
        fields_name(7)   = 'ra0h50dec-50'
        fields_name(8)   = 'ra2h30dec-50'
        fields_name(9)   = 'ra1hdec-60'
        fields_name(10)  = 'ra5h30dec-45'
        fields_name(11)  = 'ra6h30dec-55'
        fields_name(12)  = 'ra23hdec-62.5'
        fields_name(13)  = 'ra21hdec-42.5'
        fields_name(14)  = 'ra22h30dec-55'
        fields_name(15)  = 'ra23hdec-45'
        fields_name(16)  = 'ra6hdec-62.5'
        fields_name(17)  = 'ra3h30dec-42.5'
        fields_name(18)  = 'ra1hdec-42.5'
        fields_name(19)  = 'ra6h30dec-45'
        !fields_name(20)  = 'ra5h30dec-55_2011'

        do i_field=0, num_fields-1
            fields_name(i_field) = adjustl(fields_name(i_field))
        enddo

        return
    end subroutine

    subroutine params_init(ini_file)

        implicit none

        character(*),     intent(in)        :: ini_file
        type(paramfile_handle)              :: ini_handle

        call fields_info_init()

        ini_handle = parse_init(ini_file, silent=.true.)

        params%input_gal_theta_path = parse_string(ini_handle, 'gal_theta_path')
        params%input_gal_theta_root = parse_string(ini_handle, 'gal_theta_root')

        params%input_gal_phi_path   = parse_string(ini_handle, 'gal_phi_path')
        params%input_gal_phi_root   = parse_string(ini_handle, 'gal_phi_root')
    
        params%input_alms_path      = parse_string(ini_handle, 'input_alms_path')
        params%input_alms_root      = parse_string(ini_handle, 'input_alms_root')

        params%output_path          = parse_string(ini_handle, 'output_path')
        params%output_root          = parse_string(ini_handle, 'output_root')

        params%lmax_input_alms      = parse_long(ini_handle, 'lmax_input_alms')

        params%nsims                = parse_long(ini_handle, 'num_sims', 1)
        params%nside                = parse_long(ini_handle, 'nside', 2048)
        params%lmax                 = parse_long(ini_handle, 'lmax', 4096)
        params%max_order            = parse_long(ini_handle, 'max_expand_order', 4)

        params%apply_pixel_window   = parse_lgt(ini_handle, 'apply_pixel_window', .true.)
        params%vary_beam            = parse_lgt(ini_handle, 'vary_beam', .false.)
        params%beam_file            = parse_string(ini_handle, 'beam_file')

        params%init = .true.

        write(*,*) "params_init finished"

        return
    end subroutine


    subroutine input_data_init()
        
        implicit none

        character(512)      :: filename

        integer(I4B)        :: nside, lmax
        integer(I4B)        :: num_modes

        logical             :: theta_alive, phi_alive
        integer(I4B)        :: theta_size, phi_size
        integer(I4B)        :: file_size_tot
        integer(I4B)        :: num_pixels, num_pixels_tot
        integer(I4B)        :: i_field, i_start, i_end, i, il

        integer(I4B)        :: npixtot, nmaps
        integer(I4B)        :: file_unit
        integer(I4B)        :: lmax_beam

        real(DP),   allocatable :: bell_tmp(:)
        real(DP),   pointer     :: theta(:), phi(:)

        nside = params%nside
        lmax  = params%lmax

        if (.not. params%init) then
            write(*,*) "INPUT_DATA_INIT: "
            write(*,*) "params should be initialized first"
            stop
        endif

        allocate(input_data%num_pixels(0:num_fields-1))

        file_size_tot = 0
        
        do i_field=0, num_fields-1
            theta_alive = .false.
            phi_alive   = .false.

            !! get the file size of theta and phi
            filename = trim(params%input_gal_theta_path)//'/'//trim(fields_name(i_field))//'/'//trim(params%input_gal_theta_root)
            inquire(file=trim(filename), size=theta_size, exist=theta_alive)

            filename = trim(params%input_gal_phi_path)//'/'//trim(fields_name(i_field))//'/'//trim(params%input_gal_phi_root)
            inquire(file=trim(filename), size=phi_size, exist=phi_alive)
            
            if (theta_alive .and. phi_alive .and. theta_size .eq. phi_size) then
                num_pixels = theta_size / 8
                input_data%num_pixels(i_field) = num_pixels
                file_size_tot = file_size_tot + theta_size
            else
                print*, i_field
                print*, filename
                print*, theta_alive, phi_alive
                write(*,*) "ERROR while reading theta or phi binary files"
                stop
            endif
        enddo

        num_pixels_tot = file_size_tot / 8

        allocate(input_data%gal_theta(0:num_pixels_tot-1))
        allocate(input_data%gal_phi(0:num_pixels_tot-1))

        do i_field=0, num_fields-1
            i_start = sum(input_data%num_pixels(0:i_field)) - input_data%num_pixels(i_field)
            i_end   = i_start + input_data%num_pixels(i_field) - 1

            theta => input_data%gal_theta(i_start:i_end)
            phi   => input_data%gal_phi(i_start:i_end)

            filename = trim(params%input_gal_theta_path)//'/'//trim(fields_name(i_field))//'/'//trim(params%input_gal_theta_root)
            open(newunit=file_unit, file=trim(filename), form='binary', status='old', action='read')
            read(file_unit) theta
            close(file_unit)
            
            do i=1, input_data%num_pixels(i_field)
                if (theta(i) .lt. 0.00_dp .or. theta(i) .gt. PI) then
                    write(*,*) "theta = "//trim(adjustl(string(theta(i),'(E16.7)')))
                    write(*,*) "out of the range [0, PI]"
                    stop
                endif
            enddo

            filename = trim(params%input_gal_phi_path)//'/'//trim(fields_name(i_field))//'/'//trim(params%input_gal_phi_root)
            open(newunit=file_unit, file=trim(filename), form='binary', status='old', action='read')
            read(file_unit) phi
            close(file_unit)

            do i=1, input_data%num_pixels(i_field)
                if (phi(i) .lt. 0.00_dp .or. phi(i) .gt. 2.00_dp*PI) then
                    write(*,*) "phi = "//trim(adjustl(string(theta(i),'(E16.7)')))
                    write(*,*) "out of the range [0, 2*PI]"
                    stop
                endif
            enddo
        enddo
        
        allocate(input_data%wpix(0:lmax,1:1))
        if (params%apply_pixel_window) then
            do i=0, lmax
                input_data%wpix(i,1) = 1.00_dp
            enddo
        else
            call pixel_window(input_data%wpix, nside)
        endif

        if (params%vary_beam) then
            num_modes = 6
        else
            num_modes = 1
        endif

        call check_lmax_beam(lmax_beam)

        if (lmax .gt. lmax_beam) then
            write(*,*) "WARNING: the lmax set in .ini file is higher than what is allowed in the beam file."
            write(*,*) "The lmax is set to be compatible with beam file."
            write(*,*) "lmax_beam = "//string(lmax_beam,'(I4)')
            params%lmax = lmax_beam
        endif

        lmax = params%lmax

        allocate(input_data%beam_data(0:lmax,1:num_modes))
        allocate(bell_tmp(1:num_modes))

        input_data%beam_data(0:lmax,1:num_modes) = 0.00_dp
        
        open(newunit=file_unit, file=trim(params%beam_file), status='old', action='read')
        do while (.not. eof(file_unit))
            read(file_unit,'(I6,6E16.7)') il, bell_tmp

            input_data%beam_data(il,:) = bell_tmp(:)
            if (il .eq. lmax) exit
        enddo
        close(file_unit)

        write(*,*) "input_data_init finished"

        deallocate(bell_tmp)

        return
    end subroutine


    subroutine check_lmax_beam(lmax_beam)

        integer(I4B),   intent(out) :: lmax_beam
        
        integer(I4B)    :: file_unit
        integer(I4B)    :: il
        real(DP)        :: b0, b1

        open(newunit=file_unit, file=trim(params%beam_file), status='old', action='read')
        do while (.not. eof(file_unit))
            if (params%vary_beam) then
                read(file_unit,'(I6,2E16.7)') il, b0, b1
                if (b1 .ne. b1) then
                    lmax_beam = il-1
                    return
                endif
            else
                read(file_unit,'(I6,E16.7)') il, b0
            endif
        enddo
        close(file_unit) 

        lmax_beam = il
        return
    end subroutine


    subroutine check_status(isim, run_status)

        integer(I4B),   intent(in)  :: isim
        integer(I4B),   intent(out) :: run_status

        integer(I4B)        :: i_field

        character(6)        :: s_str
        character(512)      :: filename, theta_filename
        logical             :: alive
        integer(I4B)        :: file_size, theta_file_size

        write(s_str,'(I6)') isim
        s_str = adjustl(s_str)

        do i_field=0, num_fields-1
            filename = trim(params%output_path)//'/'//trim(fields_name(i_field))//'/'//trim(params%output_root)//'_sim'//trim(s_str)
            inquire(file=trim(filename), size=file_size, exist=alive)

            if (.not. alive) then
                run_status = 0
                return
            else
                theta_filename = trim(params%input_gal_theta_path)//'/'//trim(fields_name(i_field))//'/'//trim(params%input_gal_theta_root)
                inquire(file=trim(theta_filename), size=theta_file_size)

                if (file_size .ne. theta_file_size) then
                    write(*,*) "WARNING: "
                    write(*,*) trim(filename)
                    write(*,*) "does exist. But the filesize is different from what is supposed to be."
                    write(*,*) "redo it."
                    run_status = 0
                    return
                endif
            endif
        enddo
        
        write(*,*) 'all fields done for isim = '//trim(s_str)
        run_status = 1
        return    
    end subroutine


    subroutine read_input_alms(isim)

        integer(I4B),   intent(in)  :: isim

        character(512)              :: filename
        character(4)                :: isim_str

        integer(I4B)                :: file_unit
        integer(I4B)                :: il1, il2, lmax_read, lmax

        write(isim_str,'(I4)') isim
        isim_str = adjustl(isim_str)

        lmax_read = params%lmax_input_alms
        lmax = params%lmax

        do il1=0, lmax_read
            do il2=0, lmax_read
                alms_read(1,il1,il2) = (0.0_spc,0.0_spc)
            enddo
        enddo

        filename = trim(params%input_alms_path)//'/'//trim(params%input_alms_root)//'_'//trim(isim_str)//'.bin'

        open(newunit=file_unit, file=trim(filename), status='old', form='binary')
        read(file_unit) alms_read
        close(file_unit)

        do il1=0, lmax
            do il2=0, lmax
                alms_sim(1,il1,il2) = dcmplx(alms_read(1,il1,il2))
            enddo
        enddo

        return
    end subroutine


    subroutine do_sxp_sims(isim)

        integer(I4B),   intent(in)  :: isim

        integer(I4B)            :: nside, npix, ipix
        integer(I4B)            :: lmax, lmax_input, il, i_mode

        real(DP)                :: tmp

        type(planck_rng)        :: rng_handle
        
        real(DP),       pointer :: map(:)

        nside = params%nside
        npix  = nside2npix(nside)

        lmax  = params%lmax
        lmax_input = params%lmax_input_alms

        if (.not. proj%init) call proj_init(nside, input_data%gal_theta, input_data%gal_phi)
        if (.not. allocated(beam_total)) allocate(beam_total(0:lmax,1:1))
        if (.not. allocated(map_sim))    allocate(map_sim(0:npix-1,1:1))
        if (.not. allocated(alms_read))  allocate(alms_read(1:1,0:lmax_input,0:lmax_input))
        if (.not. allocated(alms_sim))   allocate(alms_sim(1:1,0:lmax,0:lmax))

        !! this step remove the previously stored values in map_prj
        do ipix=0, proj%num_prj-1
            proj%map_prj(ipix) = 0.00_dp
        enddo

        do il=0, lmax
            beam_total(il,1) = input_data%beam_data(il,1) * input_data%wpix(il,1)
        enddo

        if (params%vary_beam) then 
            call rand_init(rng_handle, isim+1, -(isim+1)*10, (isim+1)*100, -(isim+1)*1000)
            do il=0, lmax
                tmp = 0.00_dp
                do i_mode=1, 5
                    tmp = tmp + rand_gauss(rng_handle) * input_data%beam_data(il,i_mode+1)
                enddo
                beam_total(il,1) = beam_total(il,1) * dexp(tmp)
            enddo
        endif
        
        call read_input_alms(isim)
        call alter_alm(nside, lmax, lmax, 0.0_dp, alms_sim, ' ', beam_total)
        
        map => map_sim(0:npix-1,1)
        call alm2map(nside, lmax, lmax, alms_sim, map)

        call deriv_init(map_sim, params%max_order, params%lmax)
        call taylor_exp()

        return
    end subroutine


    subroutine write_prj(isim)

        integer(I4B),   intent(in)  :: isim

        !type(ProjInfo),     intent(in), target  :: proj
        !type(data_params),  intent(in)          :: params

        character(256)      :: filename
        character(6)        :: s_str

        integer(I4B)        :: file_unit
        integer(I4B)        :: i_field, i_start, i_end

        real(DP),   pointer :: prj(:)

        write(s_str,'(I6)') isim
        s_str = adjustl(s_str)
        
        if (proj%init .and. params%init) then
            do i_field=0, num_fields-1
                i_start = sum(input_data%num_pixels(0:i_field)) - input_data%num_pixels(i_field)
                i_end   = i_start + input_data%num_pixels(i_field) - 1

                prj => proj%map_prj(i_start:i_end)
                
                !write(*,*) "writing projected map to field ", i_field
                
                filename = trim(params%output_path)//'/'//trim(fields_name(i_field))//'/'//trim(params%output_root)//'_sim'//trim(s_str)

                open(newunit=file_unit, file=trim(filename), status='unknown', action='write', form='binary')
                write(file_unit) prj
                close(file_unit)
            enddo
        else
            write(*,*) "params and proj should have been initialized"
            stop
        endif

    end subroutine
end module
