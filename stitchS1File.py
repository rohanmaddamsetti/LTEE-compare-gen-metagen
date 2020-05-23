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

S1outfile = "../results/gene-modules/figures/S1File.pdf"

## This snippet comes straight from  documentation
## at: https://gitlab.mister-muffin.de/josch/img2pdf

# convert all files ending in .png inside the directory
dirname = "../results/gene-modules/figures/I-modulon-plots"
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

with open(S1outfile,"wb") as f:
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
