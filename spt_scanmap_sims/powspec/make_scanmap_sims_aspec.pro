pro make_scanmap_sims_aspec, freq, ifield
    
    num_sims = 400L

    data_path = getenv('DATAPATH')
    freq_str = strcompress(string(freq),/remove)

    f = lps12_fieldstruct()
    field_name = f[ifield].name
    field_year = strcompress(string(f[ifield].year),/remove)

    print, field_name

    workpath = data_path + 'projects/sptsz_lowl/scanmap_sim/coadd_maps_fits'
    beamfile = data_path + 'projects/spt_x_planck/lps12/beams/blgrid_'+field_year+'_'+freq_str+'.txt'

    sidx = strcompress(string(num_sims-1),/remove) 

    ;; run 1, for xspec of halfring badkground+noise sims (beam uncert is probably included in sims)
    mapname = 'MAP.MAP'
    sim_map_root = 'scanmap_sims_coadd_'+freq_str
    dls_sav_root = 'dls_scanmap_sims_coadd_'+freq_str+'_spec'
    dls_sims_savfile = workpath+'/'+field_name+'/sim_dls/'+dls_sav_root+'_0sims'+sidx+'.sav'

    res = file_info(dls_sims_savfile)
    if (res.exists) then begin
        restore, dls_sims_savfile
    endif else begin
        run_scanmap_sims_coadd_aspec, ifield, sim_map_root, num_sims, $
        mapname=mapname, $
        workpath=workpath, beamfile=beamfile, $
        bandcenters=bandcenters, dls_sims=dls_sims, $
        dls_sav_root=dls_sav_root, $
        /delete_intfile, /resume
        save, bandcenters, dls_sims, filename=dls_sims_savfile
    endelse

end
