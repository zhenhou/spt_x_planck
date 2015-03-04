pro rewrite_planck_beams_2015, convol_wpix=convol_wpix, nominal_only=nominal_only
        
    home = getenv('HOME')
    fits_file = home+'/data/planck_data/2015/RIMO/HFI_RIMO_Beams-075pc_R2.00.fits'
    output_path = home+'/data/planck_data/2015/beams/'
    
    if (keyword_set(wpix)) then begin
        wpix = fltarr(8001,1)
        wpix[*] = 1.00e0
    endif else begin
        wpix = healpixwindow(2048)
    endelse

    extno_start = 3
    extno_end = 49

    for iext=extno_start, extno_end do begin
        res=mrdfits(fits_file, iext, hdr)

        n = n_elements(hdr)
        pos = -1
        for i=0, n-1 do begin
            pos = strpos(hdr[i], 'qb_bl_DX11d_r3_75pcY5_freqs/Bl_')
            if pos ne -1 then break
        endfor
        
        if pos eq -1 then begin
            print, "something wrong while reading"
            print, str
            print, hdr
            stop
        endif

        pos = strpos(hdr[i],'Bl_')
        str = strcompress(strmid(hdr[i],pos+3),/remove)

        is_nominal = 0
        if (keyword_set(nominal_only) or n_tags(res) eq 1) then begin
            str+='_nominal'
            is_nominal = 1
        endif

        if (keyword_set(convol_wpix)) then $
        output_file = output_path+'hfi_beam_'+str+'_wpix_R2.00.txt' $
        else $
        output_file = output_path+'hfi_beam_'+str+'_R2.00.txt'
        openw, 5, output_file
        for il=0L, 4000L do begin
            if (is_nominal) then begin
                printf, 5, format='(I6,E16.7)', il, (res.nominal)[il]*wpix[il,0]
            endif else begin
                if (res.EIGEN_1)[il] ne (res.EIGEN_1)[il] then break
                printf, 5, format='(I6,6E16.7)', il, (res.nominal)[il]*wpix[il,0], (res.EIGEN_1)[il], (res.EIGEN_2)[il], (res.EIGEN_3)[il], (res.EIGEN_4)[il], (res.EIGEN_5)[il]
            endelse
        endfor
        close, 5

        
    endfor
    
end
