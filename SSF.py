#!/usr/bin/env python3
import numpy as np
import sys,argparse
import glob,os.path
import itertools
from math import pi,sin,cos,ceil
from cmath import exp
np.set_printoptions(threshold=np.inf)

#***********************************************************************
# Parsing arguments
#***********************************************************************
parser = argparse.ArgumentParser(description='Calculate spin structure factor for a chain using <Si.Sj> output from dmrgpp. In the output: 2nd column is sin*sin, and 3rd column is cos.')
parser.add_argument('L', type=int, help='total number of sites')
parser.add_argument('files', type=str, nargs='+', help='corr_szsz.dat, corr_spsm.dat, corr_smsp.dat')
parser.add_argument('-q','--alternative-qx', dest='altqx', action='store_const', default=False, const=True, help='use discretization qx=k*pi/Lrung (inexact, but provides points at qx=pi,pi/2,...) in place of the correct qx=k*pi/(Lrung+1)')
parser.add_argument('-SDW','--subtract-SDW', type=str, nargs='+', metavar='SDW', dest='SDW', default='', help='Provide {ev_sz.dat, ev_Sz.dat} to subtract SDW from the <Si.Sj> matrix, i.e., we now define correlation function as C(i,j) = <Si.Sj> - <Si><Sj> = <Si.Sj> - 3<Sz_i><Sz_j>. (Note that the last line follows from SU(2) symmetry which can be falsely broken due to numerical errors...)')
args = parser.parse_args()

#***********************************************************************
# Loading files and initializing variables
#***********************************************************************
L = args.L # L = TotalNumberOfSites

with open(args.files[0]) as f:
	n_lines = sum(1 for line in f)
	f.seek(0)
	szsz = np.loadtxt(f, dtype='float', skiprows=n_lines-L)

if (len(args.files) > 1):
	with open(args.files[1]) as f:
		n_lines = sum(1 for line in f)
		f.seek(0)
		spsm = np.loadtxt(f, dtype='float', skiprows=n_lines-L)
		# It is not necessary to calculate <gs|sminus_i splus_j|gs>,
		# since we have the following relation:
		smsp = spsm 

# However, if one calculates <gs|S- S+|gs> explicitly, it can be used here:
if (len(args.files) > 2):
	with open(args.files[2]) as f:
		n_lines = sum(1 for line in f)
		f.seek(0)
		smsp = np.loadtxt(f, dtype='float', skiprows=n_lines-L)

#**********************************************************************
# <SiSj> = <SzSz + 0.5*(S+S- + S-S+)>
#**********************************************************************
ss = szsz + 0.5*(spsm+smsp)
ss = ss + np.tril(np.transpose(ss), k=-1) # reflecting the matrix

# if provided a file with <Sz>, we subtract SDW
if args.SDW: 
	# first, read sz, Sz into arrays
	with open(args.SDW[0]) as f:
		n_lines = sum(1 for line in f)
		f.seek(0)
		sz = np.loadtxt(f, dtype='float', skiprows=n_lines-L)[:,1]
	if (len(args.SDW) > 1):
		with open(args.SDW[1]) as f:
			n_lines = sum(1 for line in f)
			f.seek(0)
			Sz = np.loadtxt(f, dtype='float', skiprows=n_lines-L)[:,1]
	SzTot = sz+Sz
	# then, we subtract SDW, i.e., <SS> --> <SS> - <Sz><Sz>
	#-------------
	# Note that formally we should have:
	# <SS> --> <SS> -<Sx><Sx> -<Sy><Sy> -<Sz><Sz> = <SS> -3*<Sz><Sz> [due to SU2]
	# but for conserved Sz, we have <Sx>=<Sy>=0, since <S+>=<S->=0,
	# thus only <Sz><Sz> is nonzero.
	#-------------
	for i, j in itertools.product(range(L), range(L)):
		ss[i,j] -= SzTot[i]*SzTot[j]

q = np.arange(1, L+1) * pi/(L+1.0) # [alternatively --> q = np.arange(-L/2, L/2) * 2*pi/L]
if args.altqx == True:
	q = np.arange(1, L+1) * pi/(L) # diff. discretization to have points at pi, pi/2, etc. (correct in the limit L->infty)
ssf_s = np.zeros(L) # def. with sine
ssf_c = np.zeros(L) # def. with cosine [i.e., exp(iq(l-m)) for C(r)=C(-r)]
for i in range(L):
	for l, m in itertools.product(range(L), range(L)): # l,m = 0,..,L-1
		ssf_s[i] += 2.0/(L+1.0) * sin(q[i]*l+1) * sin(q[i]*m+1) * ss[l, m]
		ssf_c[i] += 1.0/L * cos(q[i]*(l-m)) * ss[l, m]

print("# %s" % L)
print("#              <S^2> = %g" % (np.sum(ss.diagonal())/L))
print("# (sin*sin) sum rule = %g" % (np.sum(ssf_s)/L))
print("#     (cos) sum rule = %g" % (np.sum(ssf_c)/L))
np.savetxt(sys.stdout, np.c_[q, ssf_s, ssf_c], newline='\n')				
