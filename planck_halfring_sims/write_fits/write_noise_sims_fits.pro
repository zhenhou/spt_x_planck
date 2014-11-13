pro write_noise_sims_fits, tgz=tgz, copy=copy

    resolve_routine, 'proj_planck_sptsz'
    ;resolve_routine, 'array_1dto2d'

    num_sims = 400
    num_fields = 20

    f = lps12_fieldstruct()

    bin_path = '/home/zhenhou/scratch-data/projects/spt_x_planck/halfring_noise_sims/'

    for i_field=0, num_fields-1 do begin
        field_name = f[i_field].name

        print, field_name

        model_fits_file = '/home/zhenhou/scratch-data/projects/spt_x_planck/reproj/'+field_name+'/hfi_CovMap_143_nominal_ringhalf_1_maxOrder4_prj.fits'

        t = read_spt_fits(model_fits_file)
        nx = long(t.mapinfo.nsidex)
        ny = long( t.mapinfo.nsidey)

        m = rem_tag(t,['DMAP','WEIGHT','PROCESSING'])

        num_pixels = nx * ny

        tmp = dblarr(num_pixels)

        for isim=0, num_sims-1 do begin
            sidx = strcompress(string(isim),/remove)

            rh1_bin_file = bin_path+field_name+'/hfi_143_noise_ringhalf_1_sim'+sidx
            rh2_bin_file = bin_path+field_name+'/hfi_143_noise_ringhalf_2_sim'+sidx

            res1 = file_info(rh1_bin_file)
            res2 = file_info(rh2_bin_file)

            if ((not res1.exists) or (not res2.exists)) then continue

            rh1_fits_file = bin_path+field_name+'/sim_maps/hfi_143_noise_ringhalf_1_sim'+sidx+'.fits'
            rh2_fits_file = bin_path+field_name+'/sim_maps/hfi_143_noise_ringhalf_2_sim'+sidx+'.fits'

            res1 = file_info(rh1_fits_file)
            res2 = file_info(rh2_fits_file)

            if (res1.exists and res2.exists) then continue
    
            openr, 5, rh1_bin_file
            readu, 5, tmp
            close, 5
            array_1dto2d, tmp, m2d, nx=nx, ny=ny
            m.map.map = float(m2d)
            write_processed_fits, m, rh1_fits_file

            openr, 5, rh2_bin_file
            readu, 5, tmp
            close, 5
            array_1dto2d, tmp, m2d, nx=nx, ny=ny
            m.map.map = float(m2d)
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
