#! /usr/bin/env python

if __name__ == '__main__':
	
	print "Importing..."
	from enthought.mayavi import mlab
	from enthought.mayavi.modules.surface import Surface
	import numpy
	import sys
	print "  done."
	
	if(len(sys.argv) != 5): 
		raise RuntimeError("Wrong number of arguments!\nUsage: grid_to_m.py sol.vtu var_name width step")
	
	mlab.options.offscreen = True
	engine = mlab.get_engine()
	vtk_file_reader = engine.open(sys.argv[1])
	vtk_file_reader.point_scalars_name = sys.argv[2]

	#ug  = vtk_file_reader.outputs[0]

	surface = Surface()
	engine.add_filter(surface, vtk_file_reader)

	# getting the data
	# you need to set the vars
	w = int(sys.argv[3])
	s = int(sys.argv[4])
	
	x_min = -w
	x_max = w
	x_step = s*1j
	y_min = -w
	y_max = w
	y_step = s*1j
	x_g, y_g, z_g = numpy.mgrid[x_min:x_max:x_step, y_min:y_max:y_step, 0:1:1j]
	res = mlab.pipeline.probe_data(surface, x_g, y_g, z_g, type="scalars")
	
	print "Writing samples to %s_%i.txt..."%(sys.argv[2],w)
	f = open("%s_%i.txt"%(sys.argv[2],w), "w")
	# now the data are accessible for you, so you can easly construct your matrix (x, y, val):
	for i in range(int(x_step.imag)):
		for j in range(int(y_step.imag)):
		    f.write("%g\t" % res[i,j,0])
		f.write("\n")
	f.close()
	print "   done."
	#mlab.show()
