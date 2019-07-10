import os
import sys
import glob
import calc_sectors
import matplotlib.pyplot as plot
import numpy as np
import re

nruns = 80

if len(sys.argv) > 1 :
    folder_postfix = sys.argv[1]
else:
    print("usage: plot_volume_scaled.py folder_postfix")

scaleV = False
scaleL = False
xshift = False
yshift = False
if len(sys.argv) > 2 :
    for arg in argv[2:]:
        if arg == 'scaleV':
            scaleV = True
        elif arg == 'scaleL':
            scaleL = True
        elif arg == 'xshift':
            xshift = True
        elif arg == 'yshift':
            yshift = True
        else:
            print(f"Option {arg} not recognized")

colors = iter(['black','blue','red','green','magenta','yellow']*20)
symbols = iter(['o','x']*30)

folders = glob.glob("L*"+folder_postfix)

for folder in folders:
    os.chdir(folder)
    mean, sigma = calc_sectors.calc_sectors()
    os.chdir('..')

    L = int(re.search('L(.+?)'+folder_postfix, folder).group(1))
    x = np.linspace(0, sigma.shape[0]-1, sigma.shape[0])
    print(L)
    
    if xshift:
        maxindex = np.argmax(mean)
        first = max(maxindex-2,0)
        xw = x[first:first+5]
        mw = mean[first:first+5]
        poly = np.polyfit( xw, mw, 2)
        maxindex = -poly[1]/poly[0]/2
        x -= maxindex

    if scaleV:
        x = x/L**2

    if scaleL:
        x *= L
    
    if yshift:
        mean += np.log(L)

    plot.errorbar( x, mean, sigma, fmt='o', color=next(colors), label = f"L{L}" )


plot.legend(loc='upper right')
if rescale:
    plot.xlabel('Sector / L displaced')
    plot.ylabel('F + log(L)')
elif noscale:
    plot.xlabel('Sector')
    plot.ylabel('F')
else:
    plot.xlabel('Sector / Volume')
    plot.ylabel('F + log(L)')
plot.show()

