import numpy as np
import os
from scipy import optimize

filename_in  = "Record__CM-unwrapped"
filename_out = "fit-result"

# load data
t, x = np.loadtxt( filename_in, usecols=(0,3), unpack=True )

# fixed parameters
x0 = x[0]
v0 = 0.612372457981110

# fit
def fit_func_fit_a( t, a ):
    return x0 + v0*t + 0.5*a*t**2

def fit_func_fit_av( t, a, v ):
    return x0 + v*t + 0.5*a*t**2

def fit_func_fit_avx( t, a, v, x ):
    return x + v*t + 0.5*a*t**2

para_fit_a,   para_cov_fit_a   = optimize.curve_fit( fit_func_fit_a,   t, x, p0=[0.0] )
para_fit_av,  para_cov_fit_av  = optimize.curve_fit( fit_func_fit_av,  t, x, p0=[0.0,v0] )
para_fit_avx, para_cov_fit_avx = optimize.curve_fit( fit_func_fit_avx, t, x, p0=[0.0,v0,x0] )

print para_fit_a
print para_fit_av
print para_fit_avx

# output
file_out = open( filename_out, "w" )
file_out.write( "# fitting function: x(t) = x0 + v0*t + 0.5*a0*t^2\n" )
file_out.write( "# =========================================================\n" )
file_out.write( "# case1: fit a0\n" )
file_out.write( "x0 = %14.7e\n" % x0 )
file_out.write( "v0 = %14.7e\n" % v0 )
file_out.write( "a0 = %14.7e\n" % para_fit_a[0] )
file_out.write( "# =========================================================\n" )
file_out.write( "# case2: fit (a0, v0)\n" )
file_out.write( "x0 = %14.7e\n" % x0 )
file_out.write( "v0 = %14.7e\n" % para_fit_av[1] )
file_out.write( "a0 = %14.7e\n" % para_fit_av[0] )
file_out.write( "# =========================================================\n" )
file_out.write( "# case3: fit (a0, v0, x0)\n" )
file_out.write( "x0 = %14.7e\n" % para_fit_avx[2] )
file_out.write( "v0 = %14.7e\n" % para_fit_avx[1] )
file_out.write( "a0 = %14.7e\n" % para_fit_avx[0] )
file_out.close()
