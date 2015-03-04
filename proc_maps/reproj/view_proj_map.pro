pro view_proj_map, i_field, max_order, prj_map=prj_map

    planck_freq = 143
    planck_type = 'nominal_ringfull'
    ;max_order = 4
    spt_freq = 150

    fields_info = lps12_fieldstruct()
    
    planck_output_path = '/data23/hou/projects/spt_x_planck/reproj/'
    planck_output_root = 'hfi_SkyMap_'+strcompress(string(planck_freq),/remove)+'_'+planck_type+'_maxOrder'+strcompress(string(max_order),/remove)
    ;planck_output_root = 'hfi_'+strcompress(string(planck_freq),/remove)+'_'+planck_type

    prj_file = planck_output_path+'/'+fields_info[i_field].name+'/'+planck_output_root+'_prj.fits'

    print, prj_file

    res = read_spt_fits(prj_file)
    prj = expand_fits_struct(res)

    tv_spt_map, prj.map.map*1e6, resolution=prj.mapinfo.reso_arcmin, minval=-400.0, maxval=400.0;, xra=[-200,200], yra=[-200,200]

    spt_map_path  = '/data23/kstory/lps12/maps/20120620/coadds/'
    spt_fits_file = spt_map_path+'coadd_'+fields_info[i_field].name+'_50mJy.fits'
    res = read_spt_fits(spt_fits_file)
    spt_map = expand_fits_struct(res)

    tv_spt_map, spt_map.map.map*1e6, resolution=spt_map.mapinfo.reso_arcmin, minval=-400.0, maxval=400.0;, xra=[-200,200], yra=[-200,200]

    prj_map = prj.map.map*1e6

    return
end
