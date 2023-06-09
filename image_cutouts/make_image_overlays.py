###make overlays for GRGs

import numpy as np, matplotlib.pyplot as plt, aplpy, argparse
from astropy.table import Table
from astropy.io import fits
from astropy.coordinates import SkyCoord
from astropy import units as u
from astropy.wcs import WCS

plt.interactive(True)

###to make a radio overlay on an optical image call:
###>python make_image_overlays.py optical.fits optical.jpg radio.fits

#############################################################
#############################################################
###functions

def remove_extradim_info_from_header(header):
    'clean up a fits header that has 4D info for a 2D object -- e.g. FIRST!'
    
    head = header.copy()
    
    hkeys = list(head.keys())
    for key in hkeys:
        if len(key)>0:
            lastchar = key[-1]
            if lastchar == '3' or lastchar == '4':
                del head[key]
    
    return head
    

def radio_image_as_2d(file):
    'takes 4D image file and outputs 2D HDU for ease of working with astropy on flat image'
    
    hdu = fits.open(file)
    
    imdata = hdu[0].data
    ###extract 2D image array
    if hdu[0].data.ndim == 4:
        imdata = imdata[0][0]
    
        ###extract header info for 2D
        head2d = hdu[0].header.copy()
    
        hkeys = list(head2d.keys())
    
        crkeys = ['CTYPE', 'CRVAL', 'CDELT', 'CRPIX', 'CUNIT']
        cr3 = [c + '3' for c in crkeys]
        cr4 = [c + '4' for c in crkeys]
        badkeys = cr3 + cr4 + ['NAXIS3', 'NAXIS4']
    
        for key in hkeys:
            if 'PC3' in key or 'PC4' in key or '_3' in key or '_4' in key:
                badkeys.append(key)
            if 'PC03' in key or 'PC04' in key or '_03' in key or '_04' in key:
                badkeys.append(key)
            if key in badkeys:
                del(head2d[key])

        head2d['NAXIS'] = 2

    else:
        head2d = remove_extradim_info_from_header(hdu[0].header)

    ###create new 2D hdu
    hdu2d = fits.PrimaryHDU(imdata)
    hdu2d.header = head2d
    hdulist2d = fits.HDUList(hdu2d)
    
    return(hdulist2d)



def cont_levs(n_levs=6, rms=0.0035, base=2.0):
    'create list of log scaled contout levels'
    clevs = rms * base ** np.arange(1, n_levs+1)
    
    return clevs


def update_lsdr9_hdu_for_aplpy(filename, image_idx=1,
                               return_hdulist=True,
                               banddict={0: 'g', 1: 'r', 2: 'z'},
                               remkeys=['NAXIS3', 'BANDS',
                                        'BAND0', 'BAND1',
                                        'BAND2']):
    'update the hdu from the legacy survey fits file for use with aplpy'
    
    ###load up fits file and extract image data
    old_hdu = fits.open(filename)
    header = old_hdu[0].header
    data = old_hdu[0].data[image_idx]
    
    ###remove superfluous header info
    for key in remkeys:
        del(header[key])
    
    ##create new hdu
    hdu = fits.PrimaryHDU(data)
    
    ###update necessary header keys to reflect single band data
    hdu.header = header
    hdu.header['BAND'] = banddict[image_idx]
    
    if return_hdulist == True:
        hdu = fits.HDUList(hdu)
    
    return hdu


def make_rgb_overlay(ofile, jfile, rfile,
                     outname='optical-radio.png',
                     host_ra=None, host_dec=None,
                     hostcolor='lime', contour_color='cyan',
                     figsize=(8,7), fontsize=14,
                     contour_levs=cont_levs(),
                     sbsize=1*u.arcmin, sbcolor='w', sbwidth=2,
                     savefig=False, hostmarker_scale=1/20,
                     contalpha=1):
    'make overlay image with radio contours and host pos on rgb optical'
    
    ###set up hdu and wcs
    ohdu = update_lsdr9_hdu_for_aplpy(ofile)
    wcs = WCS(ohdu[0].header)
    radio = radio_image_as_2d(rfile)
    
    fig = plt.figure(figsize=figsize)
    fig = aplpy.FITSFigure(ohdu, figure=fig)
    fig.tick_labels.set_font(size=fontsize-1)
    fig.axis_labels.set_font(size=fontsize)
    fig.show_rgb(jfile)
    fig.show_contour(radio, levels=contour_levs,
                     colors=contour_color, linewidths=1.2,
                     alpha=contalpha)

    if host_ra is not None:
        hpos = SkyCoord(ra=host_ra, dec=host_dec, unit='deg')
        xh, yh = wcs.world_to_pixel(hpos)
        mrad = int(ohdu[0].header['NAXIS1']*hostmarker_scale)
        fig.show_circles(xw=xh, yw=yh, coords_frame='pixel',
                         radius=mrad, edgecolor=hostcolor, lw=1.5)
                         
    ####add in scalebar
    sblab = ' '.join([str(int(sbsize.value)), str(sbsize.unit)])
    fig.add_scalebar(sbsize, lw=sbwidth, color=sbcolor)
    fig.scalebar.set_label(sblab)
    fig.scalebar.set_font_size(fontsize)
    
    if savefig==True:
        plt.savefig(outname, dpi=180)
        plt.close()
    
    return


def parse_args():
    'parse command line arguments'
    parser = argparse.ArgumentParser(description="make optical radio overlay")
    parser.add_argument("ofits",
                        help="optical fits file")
    parser.add_argument("ojpeg",
                        help="color optical jpeg file")
    parser.add_argument("radio",
                        help="radio fits file")
    parser.add_argument("--ra", action='store',
                        type=float, default=None,
                        help="RA to overlay [deg]")
    parser.add_argument("--dec", action='store',
                        type=str, default=None,
                        help="Decl. to overlay [deg]")
    parser.add_argument("--rms", action='store',
                        type=float, default=3.5,
                        help="noise level [mJy / beam]")
    parser.add_argument("--outname", action='store',
                        type=str, default='optical_radio_overlay.png',
                        help="output filename")
    
    args = parser.parse_args()
    
    ###convert rms to Jy/bm
    args.rms = args.rms/1000
    
    return args


#############################################################
#############################################################
###main


if __name__ == '__main__':
    args = parse_args()
    contour_levels = cont_levs(rms=args.rms)
    make_rgb_overlay(ofile=args.ofits,
                     jfile=args.ojpeg,
                     rfile=args.radio,
                     outname=args.outname,
                     host_ra=args.ra, host_dec=args.dec,
                     hostcolor='lime', contour_color='cyan',
                     figsize=(6,5), fontsize=13,
                     contour_levs=contour_levels,
                     sbsize=1*u.arcmin, sbcolor='w', sbwidth=2,
                     savefig=True, hostmarker_scale=1/20,
                     contalpha=1)



