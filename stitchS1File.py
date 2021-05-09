#!/usr/bin/env python

''' stitchS1File.py by Rohan Maddamsetti.

This script uses the img2pdf python library
(https://gitlab.mister-muffin.de/josch/img2pdf)
to combine the pngs of I-modulon evolution in 
"../results/gene-modules/figures/I-modulon-plots"
into a single PDF: "../results/gene-modules/figures/S1File.pdf".

Usage: python stitchS1File.py

'''

import img2pdf
import os
import subprocess

def stitchFile(outfile, dirname):
    """ This snippet comes straight from  documentation
 at: https://gitlab.mister-muffin.de/josch/img2pdf
    """

    ''' hack to get img2pdf to work:
    these images contain transparency which cannot be retained in PDF.
    since img2pdf will not perform a lossy operator,
    remove the alpha channel using imagemagick:
    convert input.png -background white -alpha remove -alpha off output.png
    '''
    for fname in os.listdir(dirname):
        if not fname.endswith(".png"):
            continue
        if fname.startswith("no-alpha-"):
            continue
        path = os.path.join(dirname, fname)
        if os.path.isdir(path):
            continue
        new_fname = 'no-alpha-' + fname
        outpath = os.path.join(dirname, new_fname)
        subprocess.run(["convert", path, '-background', 'white','-alpha','remove','-alpha','off',outpath])

    with open(outfile,"wb") as f:
        imgs = []
        for fname in os.listdir(dirname):
            if not (fname.startswith("no-alpha-") and fname.endswith(".png")):
                continue
            path = os.path.join(dirname, fname)
            if os.path.isdir(path):
                continue
            imgs.append(path)
        ## sort the pages.
        imgs.sort()
        f.write(img2pdf.convert(imgs))
    return

S1outfile = "../results/gene-modules/figures/S1File.pdf"
# convert all files ending in .png inside the directory
S1dirname = "../results/gene-modules/figures/I-modulon-plots"
        
stitchFile(S1outfile, S1dirname)

S2outfile = "../results/gene-modules/figures/S2File.pdf"
# convert all files ending in .png inside the directory
S2dirname = "../results/gene-modules/figures/all-pops-I-modulon-plots"

stitchFile(S2outfile, S2dirname)
