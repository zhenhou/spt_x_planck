pro plot_diff_field, ifield
    nsims = 50
    sidx = strcompress(string(nsims-1),/remove)

    f = lps12_fieldstruct()

    field_name = f[ifield].name
    
    lps12_end2end_spectrum_sav = '/data/hou/projects/spt_x_planck/lps12/end2end/lps12_end2end_spectrum_allfields.sav'
    read_spt_end2end_spectrum, il_spt, dls_spt, savfile=lps12_end2end_spectrum_sav

    dls_spt_field = dls_spt[*,ifield] * 0.72
    
    planck_halfring_xspec_sav = '/data/hou/projects/spt_x_planck/powspec/planck_2013_halfring_cross/end2end/hfi_143_halfring_xspec_allfields.sav'
    read_planck_halfring_xspec, il_planck, dls_planck, savfile=planck_halfring_xspec_sav

    dls_planck_field = dls_planck[*,ifield]

    ip_spt = where(il_spt ge 500 and il_spt le 2500)
    ip_planck = where(il_planck ge 500 and il_planck le 2500)

    diff = dls_spt_field[ip_spt] - dls_planck_field[ip_planck]

    sim_sav = '/data23/hou/projects/spt_x_planck/powspec/planck_2013_halfring_noise_sims_cross/'+field_name+'/sim_dls/dls_halfring_noise_cross_0sims'+sidx+'.sav'
    restore, sim_sav
    il_sims = bandcenters
    ip_sims = where(il_sims ge 500 and il_sims le 2500)
    dls_sims = dls_noise_sims[ip_sims,*]
    
    ct_sims = n_elements(ip_sims)
    weight = dblarr(nsims)
    weight[*] = 1.00d0

    med_sims = dblarr(ct_sims)
    s1_sims = dblarr(ct_sims,2)
    s2_sims = dblarr(ct_sims,2)

    for i=0, ct_sims-1 do begin
        array = reform(dls_sims[i,*])
        get_conflev, array, weight, report, weight_accu=weight_accu

        med_sims[i] = report[0]
        s1_sims[i,0] = report[1]
        s1_sims[i,1] = report[2]
        s2_sims[i,0] = report[3]
        s2_sims[i,1] = report[4]
    endfor

    dl_str = '!13D!D!12l!N!6'
    
    window, xsize=800, ysize=600
    xyouts, 0, 0, '!6 '

    plot, il_spt[ip_spt], diff, xstyle=1, ystyle=1, position=[0.15,0.1,0.95,0.9], xtitle='!12l!6', ytitle=dl_str+' [!7l!6K!U2!N]', $
    charsize=2, thick=2, xrange=[600, 2500], yrange=[-500,500]
    
    errplot, il_sims[ip_sims], s1_sims[*,0], s1_sims[*,1], width=0.01
end
