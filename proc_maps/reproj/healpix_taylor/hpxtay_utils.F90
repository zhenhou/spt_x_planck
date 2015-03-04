module hpxtay_utils

    use healpix_types
    use pix_tools,      only: nside2npix
    use fitstools,      only: getsize_fits, input_map
    use paramfile_io
    use misc_utils,     only: string

    use healpix_taylor

    implicit none

    type data_params
        logical         :: init = .false.
        integer(I4B)    :: hpx_map_nside
        character(512)  :: hpx_map_file
        character(512)  :: hpx_mask_file
        
        integer(I4B)    :: num_fields

        character(512), allocatable :: gal_theta_file(:)
        character(512), allocatable :: gal_phi_file(:)
        character(512), allocatable :: out_prj_file(:)
        
        integer(I4B)                :: lmax
        integer(I4B)                :: max_order = 2
        logical                     :: has_pol = .false.
        logical                     :: use_mask = .false.

        integer(I4B),   allocatable :: num_pixels(:)
        real(DP),       allocatable :: hpx_map(:,:)
        real(DP),       allocatable :: hpx_mask(:,:)

        real(DP),   allocatable     :: gal_theta(:)
        real(DP),   allocatable     :: gal_phi(:)

    end type

    contains

    subroutine params_init(ini_file, params)
        
        character(512),     intent(in)      :: ini_file
        type(data_params),  intent(out)     :: params

        type(paramfile_handle)              :: ini_handle
    
        character(3)                        :: field_str
        integer(I4B)                        :: num_fields, i_field

        ini_handle = parse_init(ini_file, silent=.true.)

        num_fields             = parse_int(ini_handle, 'num_fields')
        params%num_fields      = num_fields

        allocate(params%gal_theta_file(0:num_fields-1))
        allocate(params%gal_phi_file(0:num_fields-1))
        allocate(params%out_prj_file(0:num_fields-1))
        allocate(params%num_pixels(0:num_fields-1))

        params%hpx_map_nside   = parse_long(ini_handle, 'hpx_map_nside')
        params%hpx_map_file    = parse_string(ini_handle, 'hpx_map_file')

        do i_field=0, num_fields-1
            write(field_str, '(I3)') i_field
            field_str = adjustl(field_str)

            params%gal_theta_file(i_field)  = parse_string(ini_handle, 'gal_theta_file'//trim(field_str))
            params%gal_phi_file(i_field)    = parse_string(ini_handle, 'gal_phi_file'//trim(field_str))
            params%out_prj_file(i_field)    = parse_string(ini_handle, 'out_prj_file'//trim(field_str))
        enddo

        params%has_pol         = parse_lgt(ini_handle, 'has_pol', .false.)
        params%use_mask        = parse_lgt(ini_handle, 'use_mask', .false.)
        
        if (params%use_mask) then
            params%hpx_mask_file    = parse_string(ini_handle, 'hpx_mask_file')
        endif

        params%lmax            = parse_long(ini_handle, 'expand_lmax', 2*params%hpx_map_nside)
        params%max_order       = parse_int(ini_handle, 'max_expand_order', 2)

        if (params%has_pol) then
            write(*,*) "Polarization Q/U is not supported yet."
            stop
        endif

        params%init = .true.

        call parse_finish(ini_handle)

        return
    end subroutine


    subroutine read_data_params(params)
        
        type(data_params),  intent(inout),  target  :: params
        
        logical             :: theta_alive, phi_alive
        integer(I4B)        :: theta_size, phi_size
        integer(I4B)        :: file_size_tot

        integer(I4B)        :: file_unit
        integer(I4B)        :: nside, npix, p, nmaps
        integer(I4B)        :: ip, ipix, i
        integer(I4B)        :: npixtot, num_pixels, num_pixels_tot
        integer(I4B)        :: num_fields, i_field
        integer(I4B)        :: i_start, i_end

        real(DP),   pointer :: theta(:), phi(:)

        if (.not. params%init) then
            write(*,*) "do params_init first before reading data"
            stop
        endif

        num_fields = params%num_fields

        file_size_tot = 0
        
        do i_field=0, num_fields-1
            theta_alive = .false.
            phi_alive   = .false.

            !! get the file size of theta and phi
            inquire(file=trim(params%gal_theta_file(i_field)), size=theta_size, exist=theta_alive)
            inquire(file=trim(params%gal_phi_file(i_field)), size=phi_size, exist=phi_alive)
            
            if (theta_alive .and. phi_alive .and. theta_size .eq. phi_size) then
                num_pixels = theta_size / 8
                params%num_pixels(i_field) = num_pixels
                file_size_tot = file_size_tot + theta_size
            else
                write(*,*) "ERROR while reading theta or phi binary files"
                stop
            endif
        enddo

        num_pixels_tot = file_size_tot / 8
                
        allocate(params%gal_theta(0:num_pixels_tot-1))
        allocate(params%gal_phi(0:num_pixels_tot-1))

        do i_field=0, num_fields-1
            i_start = sum(params%num_pixels(0:i_field)) - params%num_pixels(i_field)
            i_end   = i_start + params%num_pixels(i_field) - 1

            theta => params%gal_theta(i_start:i_end)
            phi   => params%gal_phi(i_start:i_end)

            open(newunit=file_unit, file=trim(params%gal_theta_file(i_field)), form='binary', status='old', action='read')
            read(file_unit) theta
            close(file_unit)

            do i=1, params%num_pixels(i_field)
                if (theta(i) .lt. 0.00_dp .or. theta(i) .gt. PI) then
                    write(*,*) "theta = "//trim(adjustl(string(theta(i),'(E16.7)')))
                    write(*,*) "out of the range [0, PI]"
                    stop
                endif
            enddo

            open(newunit=file_unit, file=trim(params%gal_phi_file(i_field)), form='binary', status='old', action='read')
            read(file_unit) phi
            close(file_unit)

            do i=1, params%num_pixels(i_field)
                if (phi(i) .lt. 0.00_dp .or. phi(i) .gt. 2.00_dp*PI) then
                    write(*,*) "phi = "//trim(adjustl(string(theta(i),'(E16.7)')))
                    write(*,*) "out of the range [0, 2*PI]"
                    stop
                endif
            enddo
        enddo

        nside = params%hpx_map_nside
        npix = nside2npix(nside)
        if (params%has_pol) then
            p=3
        else
            p=1
        endif

        npixtot = getsize_fits(trim(params%hpx_map_file), nmaps=nmaps)
        if (npixtot .eq. npix) then
            if (nmaps .eq. p) then
                allocate(params%hpx_map(0:npix-1,1:p))
                call input_map(trim(params%hpx_map_file), params%hpx_map, npix, nmaps)
            else
                write(*,*) "# of maps in hpx_map_file not consistent with has_pol"
                stop
            endif
        else
            write(*,*) "# of pixels in fits file does not match nside setting"
            stop
        endif

        if (params%use_mask) then
            npixtot = getsize_fits(trim(params%hpx_mask_file), nmaps=nmaps)
            if (npixtot .eq. npix) then
                if (nmaps .eq. p) then
                    allocate(params%hpx_mask(0:npix-1,1:p))
                    call input_map(trim(params%hpx_mask_file), params%hpx_mask, npix, nmaps)
                else
                    write(*,*) "# of maps in hpx_map_file not consistent with has_pol"
                    stop
                endif
            else
                write(*,*) "# of pixels in fits file does not match nside setting"
                stop
            endif
            
            do ip=1, p
                do ipix=0, npix-1
                    params%hpx_map(ipix,ip) = params%hpx_map(ipix,ip) * params%hpx_mask(ipix,ip)
                enddo
            enddo
        endif
        
        return
    end subroutine


    subroutine hpx_tayexp(params, proj, deriv)
        
        type(data_params),  intent(in)      :: params
        type(ProjInfo),     intent(out)     :: proj
        type(DerivInfo),    intent(out)     :: deriv

        integer(I4B)            :: i_field

        real(DP),   pointer     :: theta
        real(DP),   pointer     :: phi

        call proj_init(params%hpx_map_nside, params%gal_theta, params%gal_phi, proj)
        call deriv_init(params%hpx_map, proj, deriv, params%max_order, params%lmax)
        call taylor_exp(proj, deriv)

        return
    end subroutine


    subroutine write_prj(params, proj)

        type(ProjInfo),     intent(in), target  :: proj
        type(data_params),  intent(in)          :: params

        integer(I4B)        :: file_unit
        integer(I4B)        :: i_field, i_start, i_end

        real(DP),   pointer :: prj(:)
        
        if (proj%init .and. params%init) then
            do i_field=0, params%num_fields-1
                i_start = sum(params%num_pixels(0:i_field)) - params%num_pixels(i_field)
                i_end   = i_start + params%num_pixels(i_field) - 1

                prj => proj%map_prj(i_start:i_end)
                
                write(*,*) "writing projected map to field ", i_field
                open(newunit=file_unit, file=trim(params%out_prj_file(i_field)), status='unknown', action='write', form='binary')
                write(file_unit) prj
                close(file_unit)
            enddo
        else
            write(*,*) "params and proj should have been initialized"
            stop
        endif

    end subroutine


    subroutine params_finalize(params)

        type(data_params),  intent(inout)   :: params

        deallocate(params%gal_theta_file)
        deallocate(params%gal_phi_file)
        deallocate(params%out_prj_file)

        deallocate(params%num_pixels)
        deallocate(params%hpx_map, params%hpx_mask)
        deallocate(params%gal_theta, params%gal_phi)
    
        params%hpx_map_nside = 0
        params%hpx_map_file  = ' '

        params%num_fields = 0
        params%lmax       = 0

        params%init = .false.

        return

    end subroutine

end module
