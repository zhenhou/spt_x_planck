pro write_sims_fits, tgz=tgz, copy=copy

    resolve_routine, 'proj_planck_sptsz'
    ;resolve_routine, 'array_1dto2d'

    num_sims = 400
    num_fields = 20

    f = lps12_fieldstruct()

    noise_bin_path = '/home/zhenhou/scratch-data/projects/spt_x_planck/halfring_noise_sims/'
    signal_bin_path = '/home/zhenhou/scratch-midway/projects/spt_x_planck/sims/signal/hfi_143_signal_sims/'

    for i_field=0, num_fields-1 do begin
        field_name = f[i_field].name

        print, field_name

        model_fits_file = '/home/zhenhou/scratch-data/projects/spt_x_planck/reproj/'+field_name+'/hfi_CovMap_143_nominal_ringhalf_1_maxOrder4_prj.fits'

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

            rh1_bin_file = noise_bin_path+field_name+'/hfi_143_noise_ringhalf_1_sim'+sidx
            rh2_bin_file = noise_bin_path+field_name+'/hfi_143_noise_ringhalf_2_sim'+sidx
            sig_bin_file = signal_bin_path+field_name+'/hfi_143_signal_sim'+sidx
    
            res0 = file_info(sig_bin_file)
            res1 = file_info(rh1_bin_file)
            res2 = file_info(rh2_bin_file)

            if ((not res0.exists) or (not res1.exists) or (not res2.exists)) then continue

            rh1_fits_file = noise_bin_path+field_name+'/sim_maps/hfi_143_ringhalf_1_sim'+sidx+'.fits'
            rh2_fits_file = noise_bin_path+field_name+'/sim_maps/hfi_143_ringhalf_2_sim'+sidx+'.fits'

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
            cd, bin_path
            spawn, 'tar -czf '+field_name+'.tgz '+field_name+'/sim_maps/*.fits'
        endif

        if (keyword_set(copy)) then begin
            spawn, 'scp '+field_name+'.tgz hou@spt.uchicago.edu:/data23/hou/projects/spt_x_planck/powspec/planck_2013_halfring_noise_sims_cross/'
        endif
    endfor
end
