program test

    use healpix_types
    use pix_tools,  only: nside2npix, pix2ang_ring, ang2pix_ring
    use alm_tools,  only: create_alm, map2alm, alm2map

    use rngmod,     only: rand_init, planck_rng

    use healpix_taylor
    
    implicit none

    integer(I4B),   parameter   :: nside = 2048
    integer(I4B),   parameter   :: lmax = 4000
    character(*),   parameter   :: field_file = 'data_test/pixels_ring_150_ra5h30dec-55_2008.bin'
    
    character(80)               :: header(1:60)

    integer(I4B)                :: npix, ipix
    integer(I4B)                :: fsize, num_proj, i_proj

    real(DP)                    :: theta_tmp, phi_tmp
    real(DP)                    :: zbounds(1:2)
    real(DP)                    :: w8ring(1:2*8192,1:1)

    real(DP)                    :: vec(1:3)

    type(ProjInfo)              :: proj
    type(DerivInfo)             :: deriv
    type(planck_rng)            :: rng_handle
    
    
    integer(I4B),   allocatable :: pix_proj(:)
    real(DP),       allocatable :: theta(:), phi(:) 

    real(DP),       allocatable :: map(:)
    real(DP),       allocatable :: map_highres(:)
    
    real(DP),       allocatable :: map_dp(:,:)
    real(DP),       allocatable :: der1(:,:), der2(:,:)
    real(DP),       allocatable :: der_save(:,:)

    complex(DPC),   allocatable :: alms(:,:,:)

    header = " "
    zbounds = (/-1.0_dp, 1.0_dp/)
    call input_w8ring(8192, w8ring)

    npix = nside2npix(nside)

    inquire(file=field_file, size=fsize)
    num_proj = fsize/4
    
    allocate(pix_proj(0:num_proj-1))
    allocate(theta(0:num_proj-1), phi(0:num_proj-1))
    !allocate(der_save(0:num_proj-1,0:5)) !! test up to the second order
    
    open(unit=10, file=field_file, form='binary', status='old', action='read')
    read(10) pix_proj
    close(10)

    do i_proj=0, num_proj-1
        call pix2ang_ring(8192, pix_proj(i_proj), theta_tmp, phi_tmp)

        theta(i_proj) = theta_tmp
        phi(i_proj) = phi_tmp
    enddo

    allocate(map(0:npix-1), map_dp(0:npix-1,1:1))
    allocate(map_highres(0:nside2npix(8192)-1))
    allocate(alms(1:1,0:lmax,0:lmax))

    call rand_init(rng_handle, -1)
    call create_alm(8192, lmax, lmax, 0, 'data_test/test_lensedCls.fits', rng_handle, 0.00_dp, alms, header)
    write(*,*) "alm2map for nside=8192"
    call alm2map(8192, lmax, lmax, alms, map_highres)
    write(*,*) "alm2map for nside=",nside
    call alm2map(nside, lmax, lmax, alms, map)

    do ipix=0, npix-1
        map_dp(ipix,1) = map(ipix) 
    enddo
    call proj_init(nside, theta, phi, proj)
    call deriv_init(map_dp, proj, deriv, 4, lmax)
    call taylor_exp(proj, deriv)


    do i_proj=0, num_proj-1
        call pix2ang_ring(8192, pix_proj(i_proj), theta_tmp, phi_tmp)
        call ang2pix_ring(nside, theta_tmp, phi_tmp, ipix)
        
        print*, "check resulting map"
        print*, map(ipix), map_highres(pix_proj(i_proj))
        print*, proj%map_prj(i_proj)
        pause
        
    enddo

    write(*,*) "Hello world!"

end program

