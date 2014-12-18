pro do_proj
    
    spt_freq = 150
    planck_freq = 143
    planck_type = 'nominal_ringhalf_2'
    map_type    = 'CovMap'

    max_order = 4
    jackhalf  = 0


    spt_output_root = 'spt_'+strcompress(string(spt_freq),/remove)

    planck_map_file = '/data23/hou/planck_data/2013/all_sky_maps/single_field_maps/HFI_'+map_type+'_RING_'+ $
                      strcompress(string(planck_freq),/remove)+'_2048_R1.10_'+planck_type+'.fits'
    
    if planck_type eq 'nominal' then begin
        planck_output_root = 'hfi_'+map_type+'_'+strcompress(string(planck_freq),/remove)+'_'+planck_type+'_ringfull_maxOrder'+strcompress(string(max_order),/remove)
        jackhalf = 1
        half_roots = 'hfi_'+map_type+'_'+strcompress(string(planck_freq),/remove)+'_'+ $
                     ['nominal_ringhalf_1','nominal_ringhalf_2']+ $
                     '_maxOrder'+strcompress(string(max_order),/remove)
    endif else begin
        planck_output_root = 'hfi_'+map_type+'_'+strcompress(string(planck_freq),/remove)+'_'+planck_type+'_maxOrder'+strcompress(string(max_order),/remove)
    endelse

    proj_planck_sptsz, max_order, planck_map_file=planck_map_file, rewrite_thetaphi=rewrite_thetaphi, $
    spt_output_path=spt_output_path, spt_output_root=spt_output_root, spt_freq=spt_freq, $
    planck_output_path=planck_output_path, planck_output_root=planck_output_root, $
    jackhalf=jackhalf, half_roots=half_roots


end
