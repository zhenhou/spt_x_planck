pro do_proj, version=version
    
    spt_freq = 150
    ;planck_freq = 143
    map_type    = 'SkyMap'

    freqs = ['100','143','217']
    types = ['full', 'halfmission-1', 'halfmission-2','year-1','year-2']

    nfreqs = n_elements(freqs)
    ntypes = n_elements(types)

    max_order = 4
    jackhalf  = 0

    hostname = getenv('HOSTNAME')

    if (hostname eq 'spt') then home = '/home/hou/'
    if (hostname eq 'midway') then home = '/home/zhenhou/'
    data_path = home+'data/'

    if (not keyword_set(version)) then version='R2.00'


    spt_output_root = 'spt_'+strcompress(string(spt_freq),/remove)

    for ifreq=0, nfreqs-1 do begin
        for itype=0, ntypes-1 do begin

            planck_map_file = data_path+'planck_data/2015/single_field_maps/HFI_'+map_type+'_RING_'+freqs[ifreq]+'_2048_'+version+'_'+types[itype]+'.fits'
            
            planck_output_path = data_path+'projects/spt_x_planck/planck_2015/reproj/'
            planck_output_root = 'hfi_'+map_type+'_'+freqs[ifreq]+'_'+version+'_'+types[itype]+'_O'+strcompress(string(max_order),/remove)
            if types[itype] eq 'nominal' then begin 
                jackhalf = 1
                half_roots = 'hfi_'+map_type+'_'+freqs[ifreq]+'_'+version+'_'+ $
                             ['halfmission-1','halfmission-2']+ $
                             '_O'+strcompress(string(max_order),/remove)
            endif

            proj_planck_sptsz, max_order, planck_map_file=planck_map_file, rewrite_thetaphi=rewrite_thetaphi, $
            spt_output_path=spt_output_path, spt_output_root=spt_output_root, spt_freq=spt_freq, $
            planck_output_path=planck_output_path, planck_output_root=planck_output_root, $
            jackhalf=jackhalf, half_roots=half_roots

        endfor
    endfor

end
