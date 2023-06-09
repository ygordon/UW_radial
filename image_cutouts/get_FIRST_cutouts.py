###script to obtain cutout fits files from the FIRST survey

import numpy as np, argparse
from astropy.table import Table
from astropy import units as u
from astropy.coordinates import SkyCoord
from astroquery.image_cutouts.first import First

#####################################################################
#####################################################################



def first_cutout(position, size=3*u.arcmin, outname=None,
                 outdir='.', verbose=True):
    ###download FIRST cutout and write to file
    if outname is None:
        astring = position.ra.to_string(unit='hour', sep='', precision=3, pad=True)
        astring = astring[:len(astring)-2]
        dstring = position.dec.to_string(unit='deg', sep='', precision=2, pad=True,
                                             alwayssign=True)
        dstring = dstring[:len(dstring)-3] ##removes decimal
        
        outdir = outdir.removesuffix('/')
        outname = (f'{outdir}/J{astring}{dstring}_FIRST.fits')

    hdu = First.get_images(position, image_size=size)
    hdu.writeto(outname)
    
    if verbose==True:
        print(f'cutout file saved to {outname}')

    return


def parse_args():
    'parse command line arguments'
    parser = argparse.ArgumentParser(description="get FIRST cutouts from list of targets")
    parser.add_argument("targets",
                        help="table file containing sky coords of targets")
    parser.add_argument("--RAcol", action='store',
                        type=str, default='RA',
                        help="RA column name in targets")
    parser.add_argument("--DEcol", action='store',
                        type=str, default='DEC',
                        help="Decl. column name in targets")
    parser.add_argument("--RAunit", action='store',
                        type=str, default='deg',
                        help="units to assume for RA if not in table meta")
    parser.add_argument("--DEunit", action='store',
                        type=str, default='deg',
                        help="units to assume for Decl. if not in table meta")
    parser.add_argument("--size", action='store',
                        type=str, default='3arcmin',
                        help="cutout size (per side)")
    parser.add_argument("--outdir", action='store',
                        type=str, default='.',
                        help="directory to write files to")
    
    args = parser.parse_args()
    
    ###make size an astropy quantity
    args.size = u.Quantity(args.size)
    
    return args

#####################################################################
#####################################################################

if __name__ == '__main__':
    args = parse_args()
    data = Table.read(args.targets)
    ###ensure unit info in data
    if data[args.RAcol].unit is None:
        data[args.RAcol].unit = args.RAunit
    if data[args.DEcol].unit is None:
        data[args.DEcol].unit = args.DEunit
    poslist = SkyCoord(ra=data[args.RAcol], dec=data[args.DEcol])
    for pos in poslist:
        first_cutout(position=pos, size=args.size,
                     outdir=args.outdir)


