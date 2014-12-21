import numpy as np
import matplotlib.pyplot as plt

__all__ = ['ShowImage']

class ShowImage(object):

    def __init__(self, field_name, spt_info=None):
        if (spt_info is None): spt_info = \
        '/home/zhenhou/Projects/CMBtools/idl_routines/spt_tools/fields/sptsz_fields_info.txt'
        tmp_str1, tmp_str2 = np.loadtxt(spt_info,usecols=(0,0), unpack=True, dtype=np.str_, comments='#')
        field_names = tmp_str1
        del tmp_str2

        ra0_arr, dec0_arr = np.loadtxt(spt_info, usecols=(1,2), unpack=True, dtype=np.float32, comments='#')
        nsidex_arr, nsidey_arr = np.loadtxt(spt_info, usecols=(3,4), unpack=True, dtype=np.int32, comments='#')

        ip = field_names.tolist().index(field_name)
        
        self.ra0 = ra0_arr[ip]
        self.dec0 = dec0_arr[ip]
        self.nsidex = nsidex_arr[ip]
        self.nsidey = nsidey_arr[ip]

        self.xra = (self.ra0 - 0.5*self.nsidex/60.00, self.ra0 + 0.5*self.nsidex/60.00)
        self.yra = (self.dec0 - 0.5*self.nsidey/60.00, self.dec0 + 0.5*self.nsidey/60.00)

    def read_binary_map(self, binfile, dtype=None):
        if (dtype is None): dtype=np.float32

        tmp = np.fromfile(binfile, dtype=dtype)
        self.map2d = tmp.reshape((self.nsidex, self.nsidey))

        del tmp
        


    def make_image(self, scale=1.0, vmin=None, vmax=None, savefile=None, xticks=None, yticks=None, xtitle=None, ytitle=None, cbticks=None, cbformat=None, cbtitle=None):
        fig, ax = plt.subplots()

        im = plt.imshow(self.map2d*scale,cmap=plt.get_cmap('bone'), extent=(self.xra[0], self.xra[1], self.yra[0], self.yra[1]), vmin=vmin, vmax=vmax)
        cbar = fig.colorbar(im, ticks=cbticks, format=cbformat)

        ax = cbar.ax
        ax.text(1.3,0.5,cbtitle,rotation=270, fontsize=20)

        if xticks != None: plt.xticks(xticks)
        if yticks != None: plt.yticks(yticks)
        if xtitle != None: plt.xlabel(xtitle)
        if ytitle != None: plt.ylabel(ytitle)
