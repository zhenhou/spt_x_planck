pro rewrite_planck_beams_2015

    fits_file = '/data23/hou/planck_data/2013/RIMO/HFI_RIMO_R1.10.fits'
    output_path = '/data23/hou/planck_data/2013/beams/'
    
    freqs = [100, 143, 217, 353, 545, 857]
    
    if (keyword_set(wpix)) then begin
        wpix = fltarr(8001,1)
        wpix[*] = 1.00e0
    endif else begin
        wpix = healpixwindow(2048)
    endelse
end
