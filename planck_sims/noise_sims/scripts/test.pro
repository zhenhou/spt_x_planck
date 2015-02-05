pro test,nside
    !PATH = !PATH +':'+EXPAND_PATH('+/home/zhenhou/Projects/CMBtools/healpix/Healpix_3.11_intel/src/idl')
    print, nside2npix(nside)
end
