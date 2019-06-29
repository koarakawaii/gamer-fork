import argparse
import sys
import yt

# load the command-line parameters
parser = argparse.ArgumentParser( description='Plot the gas slices and projections' )

parser.add_argument( '-i', action='store', required=False, type=str, dest='prefix',
                     help='path prefix [%(default)s]', default='./' )
parser.add_argument( '-s', action='store', required=True,  type=int, dest='idx_start',
                     help='first data index' )
parser.add_argument( '-e', action='store', required=True,  type=int, dest='idx_end',
                     help='last data index' )
parser.add_argument( '-d', action='store', required=False, type=int, dest='didx',
                     help='delta data index [%(default)d]', default=1 )

args=parser.parse_args()

# take note
print( '\nCommand-line arguments:' )
print( '-------------------------------------------------------------------' )
for t in range( len(sys.argv) ):
   print str(sys.argv[t]),
print( '' )
print( '-------------------------------------------------------------------\n' )


idx_start = args.idx_start
idx_end   = args.idx_end
didx      = args.didx
prefix    = args.prefix

field     = ("deposit", "all_cic")
#center    = (356,256,256)
radius    = 10.0
dpi       = 150

yt.enable_parallelism()

ts = yt.load( [ prefix+'/Data_%06d'%idx for idx in range(idx_start, idx_end+1, didx) ] )

for ds in ts.piter():

   val, center = ds.find_max( field )
   my_sphere   = ds.sphere( center, radius )

   prof = yt.ProfilePlot( my_sphere, 'radius', field, weight_field='cell_volume', n_bins=32 )

   prof.set_unit( 'radius', 'code_length' )
   prof.set_unit( field, 'code_mass/code_length**3' )
#  prof.set_xlim( 1.0e0, 1.0e2 )
#  prof.set_ylim( field, 1.0e-6, 2.0e-2 )
   prof.save( mpl_kwargs={"dpi":dpi} )

