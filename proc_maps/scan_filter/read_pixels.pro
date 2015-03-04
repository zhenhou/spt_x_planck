pro read_pixels, ifield, nside, pixels
    home = getenv('HOME')

    f = lps12_fieldstruct()

    pixel_path = '/home/zhenhou/scratch-data/projects_bak/sptsz_lowl/setup_healpix/'

    a = ' '
    filename = pixel_path+'pixels_nside8192_150_'+f[ifield].name

    res = file_info(filename)
    if (not res.exists) then begin
        nside = -1
        pixels = -1
        return
    endif

    openr, 5, filename
    readf, 5, a
    res = strsplit(a, /extract)

    nside = long(res[0])
    npix = long(res[1])

    pix_start = long(res[2])
    pix_end = long(res[3])

    pix = lon64arr(npix)

    readu, 5, pix
    close, 5

    if (pix[0] ne pix_start or pix[npix-1] ne pix_end) then begin
        print, "something wrong in the pixel file"
        stop
    endif

    pixels = long(pix)

    return

end
