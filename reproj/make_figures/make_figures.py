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
        
    def cut_map(self, xra, yra):
        xarr = np.arange(0,self.nsidex)/60.0 + self.ra0 - 0.5*self.nsidex/60.00
        yarr = np.arange(0,self.nsidey)/60.0 + self.dec0 - 0.5*self.nsidey/60.00

        ipx = np.where((xarr >= min(xra)) & (xarr <= max(xra)))
        ipy = np.where((yarr >= min(yra)) & (yarr <= max(yra)))

        self.map2d_cut = self.map2d[ipx[0].min():(ipx[0].max()+1),ipy[0].min():(ipy[0].max()+1)]


    def make_image(self, scale=1.0, vmin=None, vmax=None, xra=None, yra=None, savefile=None, xticks=None, yticks=None, xtitle=None, ytitle=None, xfontsize=None, yfontsize=None, \
        cbticks=None, cbticks_fontsize=14, cbformat=None, cbtitle=None, cbtitle_xpos=4.0, cbtitle_fontsize=16):

        fig, ax = plt.subplots()
        
        if xra is None: xra = self.xra
        if yra is None: yra = self.yra

        im = plt.imshow(self.map2d_cut*scale,cmap=plt.get_cmap('bone'), extent=(xra[0],xra[1],yra[0],yra[1]), vmin=vmin, vmax=vmax, interpolation='bicubic')
        cbar = fig.colorbar(im, ticks=cbticks, format=cbformat)
        cbar.solids.set_edgecolor("face")
        for t in cbar.ax.get_yticklabels():
            t.set_fontsize(cbticks_fontsize)

        ax = cbar.ax
        ax.text(cbtitle_xpos,0.5,cbtitle,rotation=270, fontsize=cbtitle_fontsize, ha='center', va='center')

        if xticks != None: plt.xticks(xticks, fontsize=xfontsize)
        if yticks != None: plt.yticks(yticks, fontsize=xfontsize)
        if xtitle != None: plt.xlabel(xtitle, fontsize=xfontsize)
        if ytitle != None: plt.ylabel(ytitle, fontsize=yfontsize)
