##script to get cutout images from the legacy survey skyviewer

import numpy as np, argparse
from astropy.table import Table
from astroquery.skyview import SkyView
from astropy.coordinates import SkyCoord
from astropy import units as u


### for a list of targets in a file named 'target_list.fits' call as:
### >python get_TGSS_cutouts.py target_list.fits
### target file must include an ID/name column, RA and Dec columns

#############################################################
#############################################################
###functions


def fetch_tgss_cutout(skypos=None, size=4*u.arcmin):
    'get HDUlist of fits cutout for TGSS in native pixel scale'
    pxsize = 6.2
    npix = int(np.ceil(size.to('arcsec').value/pxsize))
    
    hdul = SkyView.get_images(position=skypos,
                              survey='TGSS ADR1',
                              pixels=npix)[0]
                              
    return hdul


def write_tgss_image_to_file(hdulist, name, outdir='.'):
    'write HDU list to fits file'
    outdir = outdir.removesuffix('/')
    fname = f'{outdir}/{name}_TGSS.fits'
    hdulist.writeto(fname)
    
    return


def batch_fetch_tgss(data, idcol='Name',
                     racol='RA', deccol='DEC',
                     size=4*u.arcmin, outdir='.',
                     posunits=('deg', 'deg')):
    'loop through input data and get TGSS cutout for all rows'
    for row in data:
        name = row[idcol]
        ra = row[racol]
        dec = row[deccol]
        pos = SkyCoord(ra=ra, dec=dec, unit=posunits)
        hdu = fetch_tgss_cutout(skypos=pos, size=size)
        write_tgss_image_to_file(hdulist=hdu, name=name,
                                 outdir=outdir)

    return


def parse_args():
    'parse command line arguments'
    parser = argparse.ArgumentParser(description="get cutouts from ls skyviewer for a list of targets")
    parser.add_argument("targets",
                        help="table file containing sky coords of targets")
    parser.add_argument("--RAcol", action='store',
                        type=str, default='RA',
                        help="RA column name in targets")
    parser.add_argument("--DEcol", action='store',
                        type=str, default='DEC',
                        help="Decl. column name in targets")
    parser.add_argument("--IDcol", action='store',
                        type=str, default='Name',
                        help="ID column name in targets")
    parser.add_argument("--size", action='store',
                        type=str, default='4arcmin',
                        help="cutout size (per side)")
    parser.add_argument("--outdir", action='store',
                        type=str, default='.',
                        help="directory to write files to")
    
    args = parser.parse_args()
    
    ###make size an astropy quantity
    args.size = u.Quantity(args.size)
    
    return args


#############################################################
#############################################################
###main

if __name__ == '__main__':
    args = parse_args()
    data = Table.read(args.targets)
    if data[args.RAcol].unit is None:
        print(f'WARNING: {args.RAcol} contains no unit information, assuming R.A. in decimal degrees')
        aunit = 'deg'
    else:
        aunit = data[args.RAcol].unit
    if data[args.DEcol].unit is None:
        print(f'WARNING: {args.DEcol} contains no unit information, assuming Decl. in decimal degrees')
        dunit = 'deg'
    else:
        dunit = data[args.DEcol].unit
    posunits = (aunit, dunit)
    batch_fetch_tgss(data=data, idcol=args.IDcol,
                     racol=args.RAcol, deccol=args.DEcol,
                     size=args.size, outdir=args.outdir,
                     posunits=posunits)
