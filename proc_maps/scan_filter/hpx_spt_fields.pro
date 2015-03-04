pro convert_pixels, nside_in, pixels_in, nside_out, pixels_out
    pix2vec_nest, nside_in, pixels_in, vec
    vec2pix_nest, nside_out, vec, ipnest

    pixels_out = ipnest[UNIQ(ipnest, SORT(ipnest))]
end

pro hpx_spt_fields, map_2048

    npix = nside2npix(8192)
    map = lonarr(npix)

    f = lps12_fieldstruct()
    
    num_fields = 20
    output_path = '/home/zhenhou/data/projects/spt_x_planck/planck_2015/scan_filter/pixel_lists/'

    for i=0, num_fields-1 do begin
        read_pixels, i, nside, pixels_in

        if (nside ne -1) then begin
            convert_pixels, 8192, pixels_in, 2048, pixels

            pixels = long64(pixels)
            map[pixels] =1

            filename = output_path + 'pixels_nside2048_150_'+f[i].name
            npix_field = n_elements(pixels)
            
            a = '2048 '
            a += strcompress(string(npix_field),/remove)+' '
            a += strcompress(string(pixels[0]),/remove)+' '
            a += strcompress(string(pixels[npix_field-1]),/remove)+' G NESTED'

            openw, 10, filename
            printf, 10, a
            writeu, 10, pixels
            close, 10
        endif
    endfor

    pixels = where(map eq 1)
    
    map_2048 = fltarr(nside2npix(2048))
    map_2048[pixels] = 1.0

end
