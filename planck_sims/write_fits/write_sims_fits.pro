pro write_sims_fits, freq, type=type, tgz=tgz, copy=copy

    resolve_routine, 'proj_planck_sptsz'
    ;resolve_routine, 'array_1dto2d'

    num_sims = 400
    num_fields = 20

    host = getenv('HOSTNAME')
    if host eq 'spt' then home='/home/hou/'
    if host eq  'midway' then home='/home/zhenhou/'
    if (not keyword_set(type)) then type='halfmission'

    freq_str = strcompress(string(freq), /remove)

    f = lps12_fieldstruct()
    signal_bin_path = home+'data/projects/spt_x_planck/sims/signal/hfi_'+freq_str+'_R2.00/'
    noise_bin_path  = home+'data/projects/spt_x_planck/sims/noise/hfi_'+freq_str+'_R2.00_'+type+'/'
    splusn_bin_path = home+'data/projects/spt_x_planck/sims/signal+noise/hfi_'+freq_str+'_R2.00_'+type+'/'

    for i_field=0, num_fields-1 do begin
        field_name = f[i_field].name

        print, field_name

        model_fits_file = home+'data/projects/spt_x_planck/planck_2015/reproj/'+field_name+'/hfi_SkyMap_'+freq_str+'_R2.00_'+type+'-1_O4_prj.fits'

        t = read_spt_fits(model_fits_file)
        nx = long(t.mapinfo.nsidex)
        ny = long(t.mapinfo.nsidey)
        
        map = fltarr(nx,ny)
        noise = create_struct('MAP',map)

        m = rem_tag(t,['WEIGHT','PROCESSING'])

        num_pixels = nx * ny
        noise1_tmp = dblarr(num_pixels)
        noise2_tmp = dblarr(num_pixels)
        signal_tmp = dblarr(num_pixels)

        for isim=0, num_sims-1 do begin
            sidx = strcompress(string(isim),/remove)

            rh1_bin_file = noise_bin_path+field_name+'/hfi_'+freq_str+'_R2.00_noise_'+type+'-1_sim'+sidx
            rh2_bin_file = noise_bin_path+field_name+'/hfi_'+freq_str+'_R2.00_noise_'+type+'-2_sim'+sidx
            sig_bin_file = signal_bin_path+field_name+'/hfi_'+freq_str+'_R2.00_signal_sim'+sidx
    
            res0 = file_info(sig_bin_file)
            res1 = file_info(rh1_bin_file)
            res2 = file_info(rh2_bin_file)

            if (res1.exists) then begin
                new_dir = noise_bin_path+field_name+'/bin_maps'
                file_mkdir, new_dir
                spawn, ['mv', rh1_bin_file, new_dir+'/'], /noshell
                rh1_bin_file = new_dir+'/hfi_'+freq_str+'_R2.00_noise_'+type+'-1_sim'+sidx
            endif

            if (res2.exists) then begin
                spawn, ['mv', rh2_bin_file, new_dir+'/'], /noshell
                rh2_bin_file = new_dir+'/hfi_'+freq_str+'_R2.00_noise_'+type+'-2_sim'+sidx
            endif

            ;if ((not res0.exists) or (not res1.exists) or (not res2.exists)) then continue

            fits_dir = splusn_bin_path+field_name+'/fits_maps'
            inf = file_info(fits_dir)
            if (not inf.exists) then file_mkdir, fits_dir

            rh1_fits_file = fits_dir+'/hfi_'+freq_str+'_R2.00_'+type+'-1_sim'+sidx+'.fits'
            rh2_fits_file = fits_dir+'/hfi_'+freq_str+'_R2.00_'+type+'-2_sim'+sidx+'.fits'

            res1 = file_info(rh1_fits_file)
            res2 = file_info(rh2_fits_file)

            if (res1.exists and res2.exists) then continue

            openr, 5, sig_bin_file
            readu, 5, signal_tmp
            close, 5
    
            openr, 5, rh1_bin_file
            readu, 5, noise1_tmp
            close, 5
            array_1dto2d, signal_tmp+noise1_tmp, m2d, nx=nx, ny=ny
            m.map.map = float(m2d)
            array_1dto2d, noise1_tmp, m2d, nx=nx, ny=ny
            m.dmap.map = float(m2d)
            write_processed_fits, m, rh1_fits_file

            openr, 5, rh2_bin_file
            readu, 5, noise2_tmp
            close, 5
            array_1dto2d, signal_tmp+noise2_tmp, m2d, nx=nx, ny=ny
            m.map.map = float(m2d)
            array_1dto2d, noise2_tmp, m2d, nx=nx, ny=ny
            m.dmap.map = float(m2d)
            write_processed_fits, m, rh2_fits_file
            
        endfor
        
        if (keyword_set(tgz)) then begin
            cd, splusn_bin_path
            spawn, 'tar -czf '+field_name+'.tgz '+field_name+'/fits_maps/*.fits'
        endif

        if (keyword_set(copy)) then begin
            spawn, 'scp '+field_name+'.tgz hou@spt.uchicago.edu:/home/hou/data/projects/spt_x_planck/sims/signal+noise/hfi_'+freq_str+'_R2.00_'+type+'/'
        endif
    endfor
end
