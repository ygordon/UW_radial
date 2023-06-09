##script to get cutout images from the legacy survey skyviewer
###ASSUMES COORDINATES IN DECIMAL DEGREES

import wget, argparse
from astropy.table import Table
from astropy import units as u

###e.g. for a list of targets in a file called test.fits where you want jpeg cutouts from sdss call:
### >python get_ls_cutouts.py targets.fits --survey=sdss --image_format=jpg

#############################################################
#############################################################
###functions


def get_ls_cutout(ra, dec, name, size=180, outdir='.',
                  survey='ls-dr10', fmt='fits'):
    'download legacy survey cutouts'
    ##use ls-dr10-early (fits and jpg) and vlass1.2
    ##
    outdir = outdir.removesuffix('/')
    
    pxscale = {'vlass1.2': 1.0,
               'ls-dr9': 0.262,
               'ls-dr10': 0.262,
               'sdss': 0.396,
               'galex': 1.5,
               'hsc-dr3': 0.168,
               'unwise-neo7': 2.75}
    pixelscale = pxscale[survey]
    
    pxsize = int(size/pixelscale)
    
    url = f'https://www.legacysurvey.org/viewer/cutout.{fmt}?ra={ra}&dec={dec}&layer={survey}&pixscale={pixelscale}&size={pxsize}'
    if survey == 'vlass1.2':
        survey = 'VLASS1'
    outfile = f'{outdir}/{name}_{survey}.{fmt}'
    print(outfile)
    
    try:
        wget.download(url=url, out=outfile)
    except:
        print(f'WARNING: {outfile} not downloaded')
        print(url)

    return


def iterate_through_input(data, idcol='Name',
                          racol='RA', deccol='DEC',
                          survey='ls-dr10', fmt='fits',
                          size_arcsec=180, outdir='.'):
    'iterate through the rows in a data table and get cutout for each row'
    for row in data:
        name = row[idcol]
        ra = row[racol]
        dec = row[deccol]
        get_ls_cutout(ra=ra, dec=dec, name=name,
                      size=size_arcsec, outdir=outdir,
                      survey=survey, fmt=fmt)
                      
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
    parser.add_argument("--survey", action='store',
                        type=str, default='ls-dr10',
                        help="survey to get cutout from")
    parser.add_argument("--image_format", action='store',
                        type=str, default='fits',
                        help="image format to download")
    parser.add_argument("--size", action='store',
                        type=str, default='3arcmin',
                        help="cutout size (per side)")
    parser.add_argument("--outdir", action='store',
                        type=str, default='.',
                        help="directory to write files to")
    
    args = parser.parse_args()
    
    ###make size an astropy quantity
    args.size = int(u.Quantity(args.size).to('arcsec').value)
    
    return args


#############################################################
#############################################################
###main

if __name__ == '__main__':
    args = parse_args()
    data = Table.read(args.targets)
    iterate_through_input(data=data, idcol=args.IDcol,
                          racol=args.RAcol, deccol=args.DEcol,
                          survey=args.survey, fmt=args.image_format,
                          size_arcsec=args.size, outdir=args.outdir)
    
