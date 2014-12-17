from numpy import *
import healpy as hp
import matplotlib.pyplot as plt

path = '/home/zhenhou/scratch-midway2/planck_data/2014/single_field_maps/'
freqs = [100, 143, 217]
dmaps = ['CovMap']
types = ['missionhalf-1','missionhalf-2','year-1','year-2']

for fq in freqs:
    for ds in dmaps:
        for tp in types:
            # HFI_CovMap_RING_143_2048_DX11d_year-1.fits
            filename = path+'HFI_'+ds+'_RING_'+str(fq)+'_2048_DX11d_'+tp+'.fits'
            dat = hp.read_map(filename)

            hp.mollview(dat*1e12, xsize=1000, title='Cov II '+str(fq)+' '+tp, unit=r'$\mu K^2_{\rm{CMB}}$', min=0, max=5e3)
            outname = 'hfi_CovII_'+str(fq)+'_'+tp+'.png'
            plt.savefig(outname)
