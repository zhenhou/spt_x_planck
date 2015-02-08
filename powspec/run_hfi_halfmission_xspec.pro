pro run_hfi_halfmission_xspec
    
    freqs = [100, 143, 217]
    nfreqs = n_elements(freqs)
    
    f = lps12_fieldstruct()
    num_fields = 20

    host = getenv('HOSTNAME')
    if host eq 'spt' then home='/home/hou/'
    if host eq  'midway' then home='/home/zhenhou/'

    mapfile_path = home+'data/projects/spt_x_planck/planck_2015/reproj/'
    workpath     = home+'data/projects/spt_x_planck/powspec/halfmission_xspec/'
    
    for ifreq=0, nfreqs-1 do begin
        fstr = strcompress(string(freqs[ifreq]), /remove)
        mapfile_root = ['hfi_SkyMap_'+fstr+'_R2.00_halfmission-1_O4_prj','hfi_SkyMap_'+fstr+'_R2.00_halfmission-2_O4_prj']
        beamfile     = home+'data/planck_data/2015/beams/hfi_beam_'+fstr+'x'+fstr+'_nominal_wpix_R2.00.txt'

        for ifield=0, num_fields-1 do begin
            field_name = f[ifield].name

            dpath = workpath+'/'+field_name+'/data'
            inf   = file_info(dpath)
            if (not inf.exists) then file_mkdir, dpath

            savfile = dpath+'/dls_halfmission_xspec.sav'
            res = file_info(savfile)
            
            if not res.exists then begin
                run_end2end_planck_field, ifield, $
                mapfile_path=mapfile_path, mapfile_root=mapfile_root, mapname='MAP.MAP', $
                beamfile=beamfile, $
                workpath=workpath, $
                bandcenters=bandcenters, spectrum=spectrum, /resume

                dls_halfmission_xspec = spectrum[1:*,*]

                save, bandcenters, dls_halfmission_xspec, filename=savfile
            endif

        endfor

    endfor

end
