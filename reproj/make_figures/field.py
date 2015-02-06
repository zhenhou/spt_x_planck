from make_figures import *
import numpy as np
import matplotlib.pyplot as plt

import os

field = 'ra5h30dec-55_2008'
si = ShowImage(field)

home = os.getenv('HOME')
hfi_path = home+'/data/projects/spt_x_planck/planck_2015/reproj/'
hfi_map = hfi_path + field + '/hfi_SkyMap_143_R2.00_full_O4_prj'

spt_path = home+'/data/projects/spt_x_planck/lps12/coadd_data/'
spt_map = spt_path + 'coadd_'+field+'_50mJy.dat'

si.read_binary_map(hfi_map,dtype=np.float64)
xra = (si.ra0 - 0.25*si.nsidex/60.00, si.ra0 + 0.25*si.nsidex/60.00)
yra = (si.dec0 - 0.25*si.nsidey/60.00, si.dec0 + 0.25*si.nsidey/60.00)
si.cut_map(xra, yra)
si.make_image(scale=1.0e6, xra=xra, yra=yra, xticks=(80.0,82.5,85.0),yticks=(-52.5,-55.0,-57.5),xtitle='RA [Deg]',ytitle='DEC [Deg]', xfontsize=14, yfontsize=14, cbticks=[-300,-200,-100,0,100,200,300], cbticks_fontsize=13, \
vmin=-350,vmax=350, cbformat='%4d', cbtitle=r'$T_{\rm CMB}\,[\,\mu\rm{K}\,]$', cbtitle_fontsize=17)
plt.savefig('hfi_'+field+'.png')

si.read_binary_map(spt_map,dtype=np.float32)
si.cut_map(xra, yra)
si.make_image(scale=1.0e6, xra=xra, yra=yra, xticks=(80.0,82.5,85.0),yticks=(-52.5,-55.0,-57.5),xtitle='RA [Deg]',ytitle='DEC [Deg]', xfontsize=14, yfontsize=14, cbticks=[-150,-100,-50,0,50,100,150], cbticks_fontsize=13, \
vmin=-150,vmax=150, cbformat='%4d', cbtitle=r'$T_{\rm CMB}\,[\,\mu\rm{K}\,]$', cbtitle_fontsize=17)
plt.savefig('spt_'+field+'.png')

