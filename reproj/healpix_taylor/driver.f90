program driver

    use healpix_types
    use healpix_taylor
    use hpxtay_utils
    
    use extension

    implicit none

    character(512)              :: ini_file
    
    type(data_params)           :: params
    type(ProjInfo)              :: proj
    type(DerivInfo)             :: deriv

    if (nArguments() == 1) call getArgument(1, ini_file)
    
    call params_init(ini_file, params)
    call read_data_params(params)
    call hpx_tayexp(params, proj, deriv)
    call write_prj(params, proj)

end program driver
