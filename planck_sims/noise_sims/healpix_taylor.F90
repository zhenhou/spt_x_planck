module healpix_taylor

    use healpix_types
    use extension,  only: getEnvironment

    use alm_tools,  only: map2alm, map2alm_iterative, alm2map_der
    use fitstools,  only: input_map
    use pix_tools,  only: npix2nside, ang2vec, pix2ang_ring, vec2pix_ring

    use misc_utils, only: brag_openmp, string

    integer(I4B),   allocatable :: deriv_index(:)

    type ProjInfo
        logical                     :: init=.false.
        integer(I4B)                :: num_prj
        real(DP),       allocatable :: theta_prj(:)
        real(DP),       allocatable :: phi_prj(:)
        real(DP),       allocatable :: map_prj(:)

        real(DP),       allocatable :: vec_prj(:,:)
        real(DP),       allocatable :: theta_hpx(:)
        real(DP),       allocatable :: phi_hpx(:)
        integer(I4B),   allocatable :: pix_hpx(:)
    end type ProjInfo

    type DerivInfo
        logical                     :: init=.false.
        integer(I4B)                :: max_order
        real(DP),       allocatable :: deriv_maps(:,:)
        !! deriv_map stores the projection patch of the derivative maps
    end type DerivInfo

    type(ProjInfo),     target  :: proj
    type(DerivInfo),    target  :: deriv

    real(DP),       allocatable :: md_map_dp(:)
    real(DP),       allocatable :: md_maps_der1(:,:)
    real(DP),       allocatable :: md_maps_der2(:,:)
    real(DP),       allocatable :: md_der2_save(:,:)
    complex(DPC),   allocatable :: md_alms(:,:,:)

    contains
    
    !! for sims, proj_init can only be done once !!
    subroutine proj_init(nside, theta, phi)
        
        implicit none
        
        integer(I4B),   intent(in)  :: nside
        !real(DP),       intent(in)  :: map(0:)
        real(DP),       intent(in)  :: theta(0:), phi(0:)
        !type(ProjInfo), intent(out) :: proj

        integer(I4B)        :: num_proj, i_proj
        integer(I4B)        :: npix, ipix
        real(DP)            :: vec_tmp(1:3)
        real(DP)        :: theta_tmp, phi_tmp

        num_proj = size(theta)

        if (num_proj .ne. size(phi)) then
            write(*,*) "theta and phi should be of same dimension. stop"
            stop
        endif

        proj%num_prj = num_proj

        allocate(proj%vec_prj(1:3,0:num_proj-1))
        allocate(proj%theta_hpx(0:num_proj-1), proj%phi_hpx(0:num_proj-1))
        allocate(proj%theta_prj(0:num_proj-1), proj%phi_prj(0:num_proj-1))
        allocate(proj%pix_hpx(0:num_proj-1))
        allocate(proj%map_prj(0:num_proj-1))

        !npix = size(map)
        !nside = npix2nside(npix)
        
        do i_proj=0, num_proj-1
            proj%theta_prj(i_proj) = theta(i_proj)
            proj%phi_prj(i_proj)   = phi(i_proj)
            
            call ang2vec(theta(i_proj), phi(i_proj), vec_tmp)
            proj%vec_prj(1:3,i_proj) = vec_tmp(1:3)
            
            ipix = 0
            call vec2pix_ring(nside, vec_tmp, ipix)
            proj%pix_hpx(i_proj) = ipix
            proj%map_prj(i_proj) = 0.00_dp
            !proj%map_prj(i_proj) = map(ipix) !! zeroth order

            call pix2ang_ring(nside, ipix, theta_tmp, phi_tmp)
            proj%theta_hpx(i_proj) = theta_tmp
            proj%phi_hpx(i_proj) = phi_tmp
        enddo

        proj%init = .true.

        write(*,*) "proj_init finished"

        return
    end subroutine proj_init


    subroutine deriv_init(map, max_order, lmax)

        real(DP),           intent(in)  :: map(0:,1:)
        !type(ProjInfo),     intent(in)  :: proj
        !type(DerivInfo),    intent(out) :: deriv
        integer(I4B),       intent(in)  :: max_order
        integer(I4B),       intent(in)  :: lmax

        
        logical         :: alive

        integer(I4B)    :: num_proj, i_proj
        integer(I4B)    :: num_terms, i_term
        integer(I4B)    :: ct_der
        integer(I4B)    :: ct_der2_1, ct_der2_3, ct_der2

        integer(I4B)    :: num_der2sav
        
        integer(I4B)    :: ipix
        integer(I4B)    :: npix, nside

        real(DP)        :: zbounds(1:2)

        logical,        allocatable :: deriv_null(:)
        logical,        allocatable :: der2sav_null(:)
        real(DP),       allocatable :: w8ring(:,:)

        !real(DP),       allocatable :: map_dp(:)
        !real(DP),       allocatable :: maps_der1(:,:)
        !real(DP),       allocatable :: maps_der2(:,:)
        !real(DP),       allocatable :: der2_save(:,:)
        !complex(DPC),   allocatable :: alms(:,:,:)

        if (.not. proj%init) then
            write(*,*) "PROJ_INIT should be done prior DERIV_INIT"
            stop
        endif

        deriv%max_order = max_order

        num_proj  = proj%num_prj
        num_terms = (max_order+2)*(max_order+1)/2
        num_der2sav = ((max_order-2)/2+2)*((max_order-2)/2+1)/2
        

        if (.not. allocated(deriv%deriv_maps)) allocate(deriv%deriv_maps(0:num_proj-1,0:num_terms-1))

        allocate(deriv_null(0:num_terms-1))
        deriv_null(0:num_terms-1) = .true.

        npix = size(map(:,1))
        nside = npix2nside(npix)
        !lmax  = 2*nside
        
        !! for cut sky, re-consider the setting of zbounds
        !! to speed up
        zbounds = (/-1.0_dp, 1.0_dp/)
        allocate(w8ring(1:2*nside,1:1)) 
        call input_w8ring(nside, w8ring)

        if (.not. allocated(md_map_dp)) allocate(md_map_dp(0:npix-1))
        if (.not. allocated(md_maps_der1)) allocate(md_maps_der1(0:npix-1,1:2))
        if (.not. allocated(md_maps_der2)) allocate(md_maps_der2(0:npix-1,1:3))
        if (.not. allocated(md_der2_save)) allocate(md_der2_save(0:npix-1,0:num_der2sav-1))
        if (.not. allocated(md_alms)) allocate(md_alms(1:1,0:lmax,0:lmax))

        allocate(der2sav_null(0:num_der2sav-1))
        der2sav_null(0:num_der2sav-1) = .true.
        
        !! zero all the relevant variables !!
        do ipix=0, npix-1
            md_map_dp(ipix) = 0.0_dp
            md_maps_der1(ipix,1) = 0.0_dp
            md_maps_der1(ipix,2) = 0.0_dp
            md_maps_der2(ipix,1) = 0.0_dp
            md_maps_der2(ipix,2) = 0.0_dp
            md_maps_der2(ipix,3) = 0.0_dp
            md_der2_save(ipix,:) = 0.0_dp
        enddo

        do i=0, lmax
            do j=0, lmax
                md_alms(1,i,j) = 0.00_dpc
            enddo
        enddo


        
        if (max_order .gt. 0) then
            do ipix=0, npix-1
                md_der2_save(ipix,0) = map(ipix,1)
            enddo
            der2sav_null(0) = .false.
        endif

        do i_proj=0, num_proj-1
            deriv%deriv_maps(i_proj, 0) = map(proj%pix_hpx(i_proj),1)
        enddo
        
        i_term = 1
        do i=0, max_order, 2
            do j=0, max_order, 2
                if (i+j .ge. max_order) cycle

                ct_der2 = ((i+j)/2+1)*(i+j)/2/2 + j/2
                if (der2sav_null(ct_der2)) then
                    write(*,*) "FATAL ERROR: previous map missing for this step of derivatives"
                    stop
                endif

                do ipix=0, npix-1
                    md_map_dp(ipix) = md_der2_save(ipix, ct_der2)
                enddo

                call brag_openmp()

                call map2alm(nside, lmax, lmax, md_map_dp, md_alms, zbounds, w8ring)
                call alm2map_der(nside, lmax, lmax, md_alms, md_map_dp, md_maps_der1, md_maps_der2)
                
                !! save the der2 maps to der2_save for next steps !!
                ct_der2_1 = ((i+2+j)/2+1)*(i+2+j)/2/2 + j/2
                ct_der2_3 = ((i+j+2)/2+1)*(i+j+2)/2/2 + (j+2)/2

                if (ct_der2_1 .lt. num_der2sav .and. ct_der2_3 .lt. num_der2sav) then
                    if (der2sav_null(ct_der2_1)) then
                        do ipix=0, npix-1
                            md_der2_save(ipix, ct_der2_1) = md_maps_der2(ipix,1)
                        enddo
                        der2sav_null(ct_der2_1) = .false.
                    endif

                    if (der2sav_null(ct_der2_3)) then
                        do ipix=0, npix-1
                            md_der2_save(ipix, ct_der2_3) = md_maps_der2(ipix,3)
                        enddo
                        der2sav_null(ct_der2_3) = .false.
                    endif
                endif

                do i_der=i, i+2
                    do j_der=j, j+2
                        ct_der = i_der + j_der - i - j
                        if (ct_der .eq. 0 .or. ct_der .gt. 2 .or. i_der+j_der .gt. max_order) cycle

                        i_term = order2index(i_der,j_der)

                        if (deriv_null(i_term)) then 
                            do i_proj=0, num_proj-1
                                ipix = proj%pix_hpx(i_proj)
                
                                if (ct_der .eq. 1) deriv%deriv_maps(i_proj, i_term) = md_maps_der1(ipix, j_der-j+1)
                                if (ct_der .eq. 2) deriv%deriv_maps(i_proj, i_term) = md_maps_der2(ipix, j_der-j+1)
                            enddo
                            deriv_null(i_term) = .false.
                        endif
                    enddo
                enddo

            enddo
        enddo

        !deallocate(map_dp, maps_der1, maps_der2, der2_save, alms)
        deallocate(w8ring, deriv_null, der2sav_null)

        deriv%init = .true.

        write(*,*) "deriv_init finished"

        return

    end subroutine deriv_init


    subroutine proj_finalize()
        
        !type(ProjInfo),     intent(out)     :: proj

        deallocate(proj%vec_prj)
        deallocate(proj%theta_hpx, proj%phi_hpx)
        deallocate(proj%theta_prj, proj%phi_prj)
        deallocate(proj%pix_hpx)
        deallocate(proj%map_prj)

        proj%num_prj = 0
        proj%init = .false.

        return
    end subroutine


    subroutine deriv_finalize()

        !type(DerivInfo),    intent(out)     :: deriv

        deallocate(deriv%deriv_maps)
        
        deriv%max_order = 0
        deriv%init = .false.

        return
    end subroutine deriv_finalize


    subroutine taylor_exp()
        
        !type(ProjInfo),     intent(inout)   :: proj
        !type(DerivInfo),    intent(in)      :: deriv

        !! deriv%deriv_maps stores the derivatives starting from order 0
        !! proj%map_prj stores the map to be projected

        integer(I4B)            :: num_proj, i_proj
        integer(I4B)            :: max_order, num_terms
        
        integer(I4B)            :: i, j, ct
        real(DP)                :: d_theta, d_phi, theta0
        real(DP)                :: term
        
        integer(I4B),   allocatable :: i_terms(:)
        real(DP),       allocatable :: taylor_factor(:)

        if (.not. proj%init) then
            write(*,*) "Please do PROJ_INIT first before Taylor expansion"
            stop
        endif

        if (.not. deriv%init) then
            write(*,*) "Please do DERIV_INIT before Taylor expansion"
            stop
        endif

        num_proj = proj%num_prj
        max_order = deriv%max_order
        num_terms = (max_order+2)*(max_order+1)/2

        allocate(taylor_factor(0:num_terms-1))
        allocate(i_terms(0:num_terms-1))
        
        ct = 0
        do i=0, max_order
            do j=0, max_order
                if (i+j .gt. max_order) cycle
                i_terms(ct) = order2index(i,j)
                taylor_factor(ct) = 1.00_dp / dble(factor(i)*factor(j))

                ct = ct + 1
            enddo
        enddo


        do i_proj=0, num_proj-1
            theta0  = proj%theta_hpx(i_proj)
            d_theta = proj%theta_prj(i_proj) - proj%theta_hpx(i_proj)
            d_phi   = proj%phi_prj(i_proj) - proj%phi_hpx(i_proj)

            if (d_phi .gt. PI) d_phi = d_phi - 2.00_dp * PI
            if (d_phi .lt. -PI) d_phi = d_phi + 2.00_dp * PI
            
            ct = 0
            do i=0, max_order
                do j=0, max_order
                    if (i+j .gt. max_order) cycle
                    term = taylor_factor(ct) * d_theta**dble(i) * (dsin(theta0)*d_phi)**dble(j)
                    !! note the zero-th order
                    proj%map_prj(i_proj) = proj%map_prj(i_proj) + term*deriv%deriv_maps(i_proj, i_terms(ct))
                    
                    ct = ct+1
                enddo
            enddo
        enddo

        deallocate(taylor_factor, i_terms)
    end subroutine


    subroutine input_w8ring(nside, w8ring)

        integer(I4B),   intent(in)  :: nside
        real(DP),       intent(out) :: w8ring(0:,1:)
        
        logical         :: alive
        character(256)  :: hpxdir, w8file

        call getEnvironment('HEALPIX', hpxdir)

        hpxdir = adjustl(hpxdir)
        w8file = trim(hpxdir)//"/data/weight_ring_n"//trim(string(nside,"(i5.5)"))//".fits"
        inquire(file=trim(w8file), exist=alive)

        w8ring = 0.0_dp
        if (alive) then
            call input_map(trim(w8file), w8ring, 2*nside, 1, fmissval=0.0_dp)
        else
            write(*,*) "WARNING: "
            write(*,*) "RING weight file not found. Use uniform RING weight."
        endif
        w8ring(1:,1) = w8ring(1:,1) + 1.0_dp

        return
    end subroutine input_w8ring


    function order2index(i, j)

        integer(I4B),   intent(in)  :: i, j
        integer(I4B)                :: order2index

        order2index = (i+j+1)*(i+j)/2 + j

        return

    end function order2index


!      Compute the factorial of N
    function factor(n)
        
        integer(I4B),   intent(in)  :: n
        integer(I4B)                :: factor

        integer(I4B)                :: i

        factor = 1
        DO i=n, 1, -1
            factor = factor * i
        END DO

        return
    end function

    
end module healpix_taylor
