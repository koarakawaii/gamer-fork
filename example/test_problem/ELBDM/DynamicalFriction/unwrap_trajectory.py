import numpy as np
import os

file_in  = "Record__CM"
file_out = "Record__CM-unwrapped"
nx       = 512
ny       = 512
nz       = 512

data_in = np.loadtxt( file_in, usecols=(0,1,2,3,4,5) )

File = open( file_out, "w" )
File.write( "#%12s   %9s   %5s   %13s   %13s   %13s\n" %
            ("Time", "Step", "NIter", "CM-x", "CM-y", "CM-z") )

x  = data_in[0][3]
y  = data_in[0][4]
z  = data_in[0][5]
cx = 0
cy = 0
cz = 0

for t in range( data_in.shape[0] ):
   x_prev = x
   y_prev = y
   z_prev = z
   x      = data_in[t][3] + cx*nx
   y      = data_in[t][4] + cy*ny
   z      = data_in[t][5] + cz*nz

   if np.abs( x - x_prev ) > 0.5*nx:
      cx = cx - np.sign( x - x_prev )
      x  = data_in[t][3] + cx*nx

   if np.abs( y - y_prev ) > 0.5*ny:
      cy = cy - np.sign( y - y_prev )
      y  = data_in[t][4] + cy*ny

   if np.abs( z - z_prev ) > 0.5*nz:
      cz = cz - np.sign( z - z_prev )
      z  = data_in[t][5] + cz*nz

   File.write( "%13.7e   %9ld   %5d   %13.7e   %13.7e   %13.7e\n" %
               (data_in[t][0], data_in[t][1], data_in[t][2], x, y, z) )

File.close()
