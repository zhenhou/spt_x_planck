pro make_halfring_sims_xspec, freq, planck_dr, ifield

    num_sims = 400L

    data_path = getenv('DATAPATH')
    freq_str = strcompress(string(freq),/remove)
    dr_str   = strcompress(string(planck_dr),/remove)

    workpath = data_path + 'projects/spt_x_planck/powspec/planck_'+dr_str+'_sims'
    beamfile = data_path + 'planck_data/'+dr_str+'/beams/hfi_beam_'+freq_str+'x'+freq_str+'_nominal_wpix_R1.10.txt'

    f = lps12_fieldstruct()
    field_name = f[ifield].name

    print, field_name
    
    sidx = strcompress(string(num_sims-1),/remove) 

    res = file_info(dls_sims_savfile)

    if res.exists then begin
        restore, dls_sims_savfile
    endif else begin
        ;; run 1, for xspec of halfring badkground+noise sims (beam uncert is probably included in sims)
        mapname = 'MAP.MAP'
        sim_map_root = 'hfi_'+freq_str+'_ringhalf_'+['1','2']
        dls_sav_root = 'dls_hfi_'+freq_str+'_ringhalf_xspec'
        dls_sims_savfile = workpath+field_name+'/sim_dls/'+dls_sav_root+'_0sims'+sidx+'.sav'

        run_halfring_sims_xspec, ifield, sim_map_root, num_sims, $
        mapname=mapname, $
        workpath=workpath, beamfile=beamfile, $
        bandcenters=bandcenters, dls_sims=dls_sims, $
        dls_sav_root=dls_sav_root, $
        /delete_intfile, /resume
        save, bandcenters, dls_sims, filename=dls_sims_savfile

        ;; run 2, for xspec of halfring noise sims
        mapname = 'DMAP.MAP'
        sim_map_root = 'hfi_'+freq_str+'_ringhalf_'+['1','2']
        dls_sav_root = 'dls_hfi_'+freq_str+'_ringhalf_noise_xspec'
        dls_sims_savfile = workpath+field_name+'/sim_dls/'+dls_sav_root+'_0sims'+sidx+'.sav'

        run_halfring_sims_xspec, ifield, sim_map_root, num_sims, $
        mapname=mapname, $
        workpath=workpath, beamfile=beamfile, $
        bandcenters=bandcenters, dls_sims=dls_sims_noise, $
        dls_sav_root=dls_sav_root, $
        /delete_intfile, /resume
        save, bandcenters, dls_sims_noise, filename=dls_sims_savfile
    endelse
end
