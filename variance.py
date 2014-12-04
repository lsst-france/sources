import pyfits
import re
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import glob
import os
import math
from scipy import spatial

#where = '/workspace/LSST/data/08AL01/D3'
where = '/afs/in2p3.fr/home/l/lsstprod/data/DC2014/CFHTLS/output/src/08AL01/D3/'

os.chdir(where)

# list of sources: 
#   key = exposure, 
#   for each exposure: coords (ra,dec) of reconstructed sources indexed by the source id
#   
# list of fluxes: 
#   key = exposure, 
#   for each exposure: flux of reconstructed sources
#   
# for each exposure, we add a KDTree from the array of coords
#
all_sources = {}
all_fluxes = {}
trees = {}

dates = 0

def R (coord):
    x = coord[0]
    y = coord[1]
    return (math.sqrt (x*x + y*y))

for date in glob.glob('*'):
    print date
    all_sources[date] = {}
    all_fluxes[date] = {}
    sources = all_sources[date]
    fluxes = all_fluxes[date]
    ns = 0
    # this is the simple array of coords to build the KDTree
    coords = []
    for file in glob.glob('%s/r/*00.fits' % date):
        print file
        h = pyfits.open (file)
        data = h[1].data
        for row in data:
            #id = row['id']
            coord = row['coord']
            flux = row['flux_psf']
            x = coord[0]
            y = coord[1]
            coords.append ((x,y))
            sources[ns] = [R((x,y))]
            fluxes[ns] = [flux]
            ns += 1
        
    all_sources[date] = sources
    all_fluxes[date] = fluxes
    
    tree = spatial.KDTree(coords)
    trees[date] = tree

    dates += 1
    #if dates > 1:
    #   break

d0 = trees.keys()[0]

combined = {}

def init_combined ():
    global combined
    
    combined = {}
    sources = all_sources[d0]
    for ns in sources:
        combined[ns] = sources[ns]
    
fluxes = all_fluxes[d0]

# arc-seconds from radians
sec = math.pi/(180*3600)

def distance (step):
    return step*sec/5.0

def accumulate_bin (histo, bin):
    if bin in histo:
        count = histo[bin]
    else:
        count = 0
    count += 1

    histo[bin] = count

def accumulate_list (histo, bin, add_values):
    if bin in histo:
        values = histo[bin]
    else:
        values = []
    for value in add_values:
        values.append (value)

    histo[bin] = values

    
#
# we associate sources coming from different exposures (each esposure was stored as a KDTree)
# one exposure is selected as the reference (date=d0)
#
# the set of associations is stored into combined
#
# return: the distribution of accumulated match counts for all exposures 
#
def associate (dist):
    print dist/sec, 'arcsec'
    
    t0 = trees[d0]
    ms = {}
    mc = {}
    for date in trees:
        if date == d0:
            continue        
        t = trees[date]
        matches = t0.query_ball_tree (t, dist)
        other = all_sources[date]
        otherf = all_fluxes[date]
        # there is one match information per source in the reference exposure
        # each match info is a list (possibly empty) of sources from the other tree
        for i_source in range(len(matches)):
            found = matches[i_source]

            # we count associations found for this source of the reference exposure.
            # one match may associate 0 to N sources from the other exposure
            # here we accumulate the distribution of number of associations

            n_associations = len(found)
            accumulate_bin (ms, n_associations)

            coords = []
            for j_source in found:
                for coord in other[j_source]:
                    coords.append (coord)

            accumulate_list (combined, i_source, coords)
            accumulate_list (mc, n_associations, coords)

    return ms, mc


#------------------------------------

#evaluate_kdtree ()
#variance ()

init_combined ()
dist = distance (6)
associate (dist)

i = 0

fluxes = all_fluxes[d0]

fig, ax1 = plt.subplots()
ax1.set_xscale("log")
ax1.set_yscale("log")

for s in combined:
    #print s, len (combined[s]), combined[s], fluxes[s]
    print s, len (combined[s]), fluxes[s]
    i += 1
    if i > 20:
        break
        pass
    continue
    if not np.isnan (fluxes[s]):
        for r in combined[s]:
            ax1.plot (np.var(combined[s]), fluxes[s], '.')

#ax1.legend()
ax1.grid()
plt.show()
