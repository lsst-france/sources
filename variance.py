import pyfits
import re
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import glob
import os
import math
from scipy import spatial

where = '/workspace/LSST/data/06AL01/D3'
where = '/afs/in2p3.fr/home/l/lsstprod/data/DC2014/CFHTLS/output/src/06AL01/D3/'

os.chdir(where)

all_sources = {}
all_fluxes = {}
trees = None

dates = None 
d0 = None

# arc-seconds from radians
sec = math.pi/(180*3600)

#------------------------------------------------------------------------
def read_sources ():
    global dates
    global all_sources
    global all_fluxes
    global trees

    dates = 0
    trees = {}
    for date in glob.glob('*'):
	print date
	all_sources[date] = {}
	all_fluxes[date] = {}
	sources = all_sources[date]
	fluxes = all_fluxes[date]
	ns = 0
	coords = []
        images = '%s/r/*00.fits' % date
        #images = '%s/r/*.fits' % date

	for file in glob.glob(images):
	    #print file
	    h = pyfits.open (file)
	    data = h[1].data
	    for row in data:
	        #id = row['id']
		coord = row['coord']
		flux = row['flux_psf']
		x = coord[0]
		y = coord[1]
		coords.append ((x,y))
		sources[ns] = [(x,y)]
		fluxes[ns] = [flux]
		ns += 1
        
	all_sources[date] = sources
	all_fluxes[date] = fluxes
    
	tree = spatial.KDTree(coords)
	trees[date] = tree

	dates += 1

#------------------------------------------------------------------------
def distance (step):
    global sec
    return step*sec/5.0

#------------------------------------------------------------------------
#
# we associate sources coming from different exposures (each esposure was stored as a KDTree)
# one exposure is selected as the reference (date=d0)
#
# the set of associations is stored into combined
#
def associate (dist):
    global dates
    global all_sources
    global all_fluxes
    global trees
    global d0
    global t0
    global combined

    print dist
    ms = {}
    for date in trees:
        if date == d0:
            continue        
        t = trees[date]
        matches = t0.query_ball_tree (t, dist)
        other = all_sources[date]
        otherf = all_fluxes[date]
        for i_source in range(len(matches)):
            found = matches[i_source]

            # we count associations found for this source of the ref exposure.
            # one match may associate 0 to N sources from the other exposure

            n_associations = len(found)
            if n_associations in ms:
                count = ms[n_associations]
            else:
                count = 0
            count += 1

            ms[n_associations] = count

            existing_coords = combined[i_source]
            for j_source in found:
                for coord in other[j_source]:
                    existing_coords.append (coord)
            combined[i_source] = existing_coords

    return ms


#------------------------------------------------------------------------
def evaluate_kdtree ():
    global dates
    global all_sources
    global all_fluxes
    global trees
    global d0
    global t0
    global combined

    # we select one of the dates to be the reference data -> the reference tree
    # each exposure is stored 
    #  as a KDTRee of the detected sources
    #  and a list of sources
    
    # let's copy the sources of the referenced exposure in the combined list 
    # where we will accumulate all source associations 

    combined = {}
    sources = all_sources[d0]
    for ns in sources:
	combined[ns] = sources[ns]
    
    fluxes = all_fluxes[d0]

    # steps to evaluate the evolution of the association quality vs the distance 
    steps = 24

    # we setup a list of fixed-length arrays of to store association counts 
    all_associations = {}

    fig, ax1 = plt.subplots()
    #ax1.set_yscale("log")
    colors = 'rgbyc' * 10

    max_found = 0
    for step in range (1, steps):
	dist = step*sec/5.0
	associations = associate (dist)
	if len(associations) > max_found:
	    max_found = len(associations)
        
	print 'step=', step, 'associations=', associations
	for k in associations:
	    if k in all_associations:
		sums = all_associations[k]
	    else:
		sums = {}
	    sums[step] = associations[k]
	    all_associations[k] = sums

    print 'max_found=', max_found
    for k in all_associations:
	sums = all_associations[k]
	for step in range (1, steps):
	    if not step in sums:
		sums[step] = 0
	all_associations[k] = sums
            
    # now count association counts at each step.
    # using negative numbers except for count = 1
    sum = []

    first = True

    for k in all_associations:
	associations = all_associations[k]
	if first:
	    first = False
	    for a in associations:
		sum.append (0)
            
	print associations
	ax1.plot (associations.values(), '-' + colors[k], label='%d' % (k))
	i = 0
	for a in associations.values():
	    sign = -1
	    if k == 1:
		sign = 1
        
	    sum[i] += sign * a
	    i += 1
    
    ax1.plot (sum, '-o')

    ax1.legend()
    ax1.grid()
    plt.show()


#------------------------------------------------------------------------
def variance ():
    global dates
    global all_sources
    global all_fluxes
    global trees
    global d0
    global t0
    global combined

    fluxes = all_fluxes[d0]
    #print fluxes

    # l'id de chaque source est l'index dans l'array de t0
    # chaque source contient un array des coordonnes associees pour chaque exposition 

    combined = {}
    sources = all_sources[d0]
    for ns in sources:
	combined[ns] = sources[ns]
    
    fig, ax1 = plt.subplots()
    ax1.set_xscale("log")
    ax1.set_yscale("log")

    for ns in combined:
	coords = combined[ns]
	flux = fluxes[ns]
	if len(coords) == 1:
	    continue
	rs = []
	for c in coords:
	    x = c[0]
	    y = c[1]
	    r = math.sqrt (x*x + y*y)
	    rs.append (r)
	rv = np.var (rs)
	ax1.plot(flux, rv, 'b.')
    plot.show()

def distrib ():
    global dates
    global all_sources
    global all_fluxes
    global trees
    global d0
    global t0
    global combined

    combined = {}
    sources = all_sources[d0]
    for ns in sources:
        combined[ns] = sources[ns]
    
    sec = math.pi/(180*3600)
    dist = 6*sec/5.0
    associate (dist)

    colors = 'rgbyc' * 10
    
    fig, ax1 = plt.subplots()

    k = {}
    for s in combined:
        associations = combined[s]
        #print associations
    
        n = len(associations)
        
        if n == 0:
            continue
        if n > len(colors):
            continue
    
        coords = associations[0]
        ra = coords[0]
        dec = coords[1]
    
        c = n-1
        #if c != 2:
        #    continue
        if c in k:
            sum = k[c]
        else:
            sum = 0
        sum += 1
        k[c] = sum
        
        #ax1.plot(ra, dec, colors[c] + '.', label='%d' % c)
        ax1.plot(ra, dec, colors[c] + '.')

    print k
    #ax1.legend()
    ax1.grid()
    plt.show()


#------------------------------------

read_sources ()

d0 = trees.keys()[0]
t0 = trees[d0]

evaluate_kdtree ()
#variance ()
distrib()

