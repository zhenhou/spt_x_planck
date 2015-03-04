pro make_half_sims_xspec, ifield, freq, type=type

    num_sims = 400L

    host = getenv('HOSTNAME')
    if (host eq 'spt') then home='/home/hou/'
    if (host eq 'midway') then home = '/home/zhenhou/'

    if (not keyword_set(type)) then type = 'halfmission'

    freq_str = strcompress(string(freq),/remove)

    workpath = home+'data/projects/spt_x_planck/sims/signal+noise/hfi_'+freq_str+'_R2.00_'+type
    beamfile = home+'data/planck_data/2015/beams/hfi_beam_'+freq_str+'x'+freq_str+'_nominal_wpix_R2.00.txt'

    f = lps12_fieldstruct()
    field_name = f[ifield].name

    print, field_name
    
    sidx = strcompress(string(num_sims-1),/remove) 

    ;; run 1, for xspec of halfmission/year noise sims (beam uncert is not included in sims)
    mapname = 'DMAP.MAP'
    sim_map_root = 'hfi_'+freq_str+'_R2.00_'+type+'-'+['1','2']
    dls_sav_root = 'dls_hfi_'+freq_str+'_R2.00_noise_'+type+'_xspec'
    dls_sims_path    = workpath+'/'+field_name+'/sim_dls'
    inf = file_info(dls_sims_path)
    if (not inf.exists) then file_mkdir, dls_sims_path
    dls_sims_savfile = dls_sims_path+'/'+dls_sav_root+'_0sims'+sidx+'.sav'

    res = file_info(dls_sims_savfile)
    if (res.exists) then begin
        restore, dls_sims_savfile
    endif else begin
        run_sims_xspec, ifield, sim_map_root, num_sims, $
        mapname=mapname, $
        workpath=workpath, beamfile=beamfile, $
        bandcenters=bandcenters, dls_sims=dls_sims_noise, $
        dls_sav_root=dls_sav_root, sims_type=type, $
        /delete_intfile, /resume
        save, bandcenters, dls_sims_noise, filename=dls_sims_savfile
    endelse

    ;; run 2, for xspec of halfmission/year signal+noise sims
    mapname = 'MAP.MAP'
    dls_sav_root = 'dls_hfi_'+freq_str+'_R2.00_'+type+'_xspec'
    dls_sims_savfile = dls_sims_path+'/'+dls_sav_root+'_0sims'+sidx+'.sav'

    res = file_info(dls_sims_savfile)
    if (res.exists) then begin
        restore, dls_sims_savfile
    endif else begin
        run_sims_xspec, ifield, sim_map_root, num_sims, $
        mapname=mapname, $
        workpath=workpath, beamfile=beamfile, $
        bandcenters=bandcenters, dls_sims=dls_sims, $
        dls_sav_root=dls_sav_root, sims_type=type, $
        /delete_intfile, /resume
        save, bandcenters, dls_sims, filename=dls_sims_savfile
    endelse
end
