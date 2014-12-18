pro rewrite_planck_maps_2014

    freqs = ['100','143','217']
    types = ['missionhalf-1', 'missionhalf-2','year-1','year-2']

    spt_mask = '/home/hou/Projects/projects/spt_x_planck/healpix_pixels/mask_spt2500deg2_nside2048_nest.fits'
    read_fits_map, spt_mask, mask

    planck_fits_path = '/data23/hou/planck_data/2014/dpc_maps_hfi/'
    rewrite_path     = '/data23/hou/planck_data/2014/single_field_maps/'

    ;;HFI_SkyMap_100_2048_DX11d_missionhalf-1.fits
    
    num_freqs = n_elements(freqs)
    num_types = n_elements(types)

    T_field = 0
    TT_field = 4

    for i_freq=0, num_freqs-1 do begin
        for i_type=0, num_types-1 do begin
            input_file = planck_fits_path+'HFI_SkyMap_'+freqs[i_freq]+'_2048_DX11d_'+types[i_type]+'.fits'

            read_fits_map, input_file, maps

            map1 = maps[*,T_field]
            map2 = maps[*,TT_field]
            
            ;; I do median filter within 1 deg radius for the bad pixels
            ip_null = where(map1 eq -1.6375E30 and mask eq 1)
            print, freqs[i_freq], types[i_type], n_elements(ip_null)
            if (ip_null[0] ne -1) then begin
                num_null = n_elements(ip_null)

                for i_null=0L, num_null-1L do begin
                    pix2vec_nest, 2048L, ip_null[i_null], vec
                    query_disc, 2048L, vec, 1.00, list, /deg, /nest

                    if (list[0] ne -1) then begin
                        map1[ip_null[i_null]] = median(map1[list], /even)
                        map2[ip_null[i_null]] = median(map2[list], /even)
                    endif else begin
                        print, freqs[i_freq]+' '+types[i_type]
                        print, 'radius too small for median filter'
                        stop
                    endelse
                endfor
            endif

            map_ring = reorder(map1, in='nest', out='ring')
            output_file1 = rewrite_path+'HFI_SkyMap_RING_'+freqs[i_freq]+'_2048_DX11d_'+types[i_type]+'.fits'
            write_fits_map, output_file1, map_ring, /ring

            map_ring = reorder(map2, in='nest', out='ring')
            output_file2 = rewrite_path+'HFI_CovMap_RING_'+freqs[i_freq]+'_2048_DX11d_'+types[i_type]+'.fits'
            write_fits_map, output_file2, map_ring, /ring
        endfor
    endfor
end