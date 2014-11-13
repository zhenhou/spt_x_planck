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
    dls_sims_savfile = workpath+field_name+'/sim_dls/dls_halfring_noise_cross_0sims'+sidx+'.sav'

    res = file_info(dls_sims_savfile)

    if res.exists then begin
        restore, dls_sims_savfile
    endif else begin

    endelse
end
