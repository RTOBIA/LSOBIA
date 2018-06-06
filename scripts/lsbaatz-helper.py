import sys
import math

if len(sys.argv)!=7:
   print("Usage: {} nb_bands w h max_iter nb_proc mem_per_proc(Mb)".format(sys.argv[0]))
   sys.exit(1)

nb_bands = int(sys.argv[1])
w = int(sys.argv[2])
h = int(sys.argv[3])
max_iter = int(sys.argv[4])
nb_procs = int(sys.argv[5])
mem_per_proc = int(sys.argv[6])


# Compute memory used by a single pixel
margin = 1.1
graph_node_size = margin*(4+305+(8*nb_bands)+0.3*4*16 +3*8) # (multithread cost)
img_node_size = (4+4+4)*nb_bands # gdal buffer + reader buffer + extract ROI buffer
node_size = graph_node_size + img_node_size

# Image aspect ratio
ratio = float(w)/float(h)

# Number of pixels per tile
nb_pixels_per_tile = float(w*h)/nb_procs

# Estimate a target number of divisions in width from number of pixels per tile and aspect ratio
target_nb_tiles_w = int(w/math.floor(math.sqrt(nb_pixels_per_tile/ratio)))

# Find the factor of nb procs closest to the target number of divisions in width
dist = abs(target_nb_tiles_w - 1)
best_factor=1

for factor in range(2,nb_procs+1):
   if nb_procs % factor == 0:
      new_dist = abs(target_nb_tiles_w-factor)
      if new_dist<dist:
         best_factor=factor
         dist = new_dist

# Compute number of tiles in w and h
nb_tiles_w = best_factor
nb_tiles_h = nb_procs/nb_tiles_w

# Compute tile width and height
tile_w = int(math.ceil(float(w)/nb_tiles_w))
tile_h = int(math.ceil(float(h)/nb_tiles_h))

print("Partionning: {}x{}={} tiles of {}x{} pixels".format(nb_tiles_w,nb_tiles_h,nb_tiles_h*nb_tiles_w,tile_w,tile_h))

# Loop on initial number of iterations
for it in range(1,max_iter):

   # Compute tile size with margins
   margin = 2**(it+1)-2
   tile_w_margin = int(math.floor(tile_w + 2*margin))
   tile_h_margin = int(math.floor(tile_h + 2*margin))

   # Compute memory usage for this margin
   used_memory_graph = (tile_w_margin+2*margin)*(tile_h_margin+2*margin)*graph_node_size/(10**6)
   used_memory_img = (tile_w_margin+2*margin)*(tile_h_margin+2*margin)*img_node_size/(10**6)

   # Check if required memory is higher than available memory
   if used_memory_graph+used_memory_img>mem_per_proc:
      print("For {} initial iterations: margin={}, \trequired memory={} Mb, not enough memory available".format(it,margin,int(used_memory_graph+used_memory_img)))
   else:
      print("For {} initial iterations: margin={}, \trequired memory={} Mb".format(it,margin,int(used_memory_graph+used_memory_img)))

