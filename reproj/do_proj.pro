pro do_proj, planck_freq, planck_type, version
    
    spt_freq = 150
    ;planck_freq = 143
    ;planck_type = 'nominal_ringhalf_2'
    ;version = 'DX11d'
    map_type    = 'SkyMap'

    max_order = 4
    jackhalf  = 0

    hostname = getenv('HOSTNAME')

    if (hostname eq 'spt') then data_path = '/data23/hou/'
    if (hostname eq 'midway') then data_path = '/home/zhenhou/scratch-midway2/'


    spt_output_root = 'spt_'+strcompress(string(spt_freq),/remove)

    planck_map_file = data_path+'planck_data/2014/all_sky_maps/single_field_maps/HFI_'+map_type+'_RING_'+ $
                      strcompress(string(planck_freq),/remove)+'_2048_'+version+'_'+planck_type+'.fits'
    
    planck_output_root = 'hfi_'+map_type+'_'+strcompress(string(planck_freq),/remove)+'_'+version+'_'+planck_type+'_O'+strcompress(string(max_order),/remove)
    if planck_type eq 'nominal' then begin 
        jackhalf = 1
        half_roots = 'hfi_'+map_type+'_'+strcompress(string(planck_freq),/remove)+'_'+version+'_'+ $
                     ['missionhalf-1','missionhalf-2']+ $
                     '_O'+strcompress(string(max_order),/remove)
    endif

    proj_planck_sptsz, max_order, planck_map_file=planck_map_file, rewrite_thetaphi=rewrite_thetaphi, $
    spt_output_path=spt_output_path, spt_output_root=spt_output_root, spt_freq=spt_freq, $
    planck_output_path=planck_output_path, planck_output_root=planck_output_root, $
    jackhalf=jackhalf, half_roots=half_roots


end
