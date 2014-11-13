pro make_halfring_noise_xpowspec, ifield, bandcenters, dls_noise_sims
    
    num_sims = 50L

    mapfile_path = '/data23/hou/projects/spt_x_planck/reproj/'
    mapfile_root = ['hfi_CovMap_143_nominal_ringhalf_1_maxOrder4_prj','hfi_CovMap_143_nominal_ringhalf_2_maxOrder4_prj']
    workpath     = '/data23/hou/projects/spt_x_planck/powspec/planck_2013_halfring_noise_sims_cross/'
    beamfile     = '/data23/hou/planck_data/2013/beams/hfi_beam_143x143_nominal_wpix_R1.10.txt'

    f = lps12_fieldstruct()
    field_name = f[ifield].name

    print, field_name
    
    sidx = strcompress(string(num_sims-1),/remove)
    dls_sims_savfile = workpath+field_name+'/sim_dls/dls_halfring_noise_cross_0sims'+sidx+'.sav'
    
    res = file_info(dls_sims_savfile)

    if res.exists then begin
        restore, dls_sims_savfile
    endif else begin
        covmap_fits_files = mapfile_path+field_name+'/'+mapfile_root+'.fits'

        run_noise_sims_cross, ifield, covmap_fits_files, num_sims, $
        workpath=workpath, beamfile=beamfile, $
        bandcenters=bandcenters, dls_noise_sims=dls_noise_sims, $
        intfile_ident=intfile_ident, $
        /delete_intfile, $
        /resume

        save, bandcenters, dls_noise_sims, filename=dls_sims_savfile
    endelse

    return
end
