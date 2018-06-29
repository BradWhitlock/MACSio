import os, string, subprocess

description = """
This is a weak scaling study of MACSio that shows how different I/O interfaces perform
when writing data out for different mesh types. MACSio was run on LLNL's Surface cluster
to demonstrate the performance of the Conduit/Relay output plugin for MACsio. The performance
was measured for 4 Conduit protocols (including ADIOS) and compared to MACSio's HDF5 serial interface. MIFMAX
trials run MACSio so it will generate a separate file for each MPI rank at each time step.
SIF trials generate a single file for each time step. Note that in SIF mode, data are lacking
for certain protocols that are only able to add data to the file by reading back and appending
the full data in serial.
"""

# LLNL
bank="wbronze"
ppn = 16
Qdeb="pbatch" ## surface only has pbatch "pdebug"
Qreg="pbatch"

# Cori
# Add something for haswell or knl
#bank="m636"
#ppn = 68
#Qdeb="debug"
#Qreg="regular"

# dimensions to vary.
job_sizes = (
    {"cores":16,   "ppn":ppn, "queue":Qdeb, "time":60, "bank":bank},
    {"cores":64,   "ppn":ppn, "queue":Qdeb, "time":60, "bank":bank},
    {"cores":256,  "ppn":ppn, "queue":Qdeb, "time":60, "bank":bank},
    {"cores":800, "ppn":ppn, "queue":Qreg, "time":120, "bank":bank}#,
    #{"cores":4096, "ppn":ppn, "queue":Qreg, "time":30, "bank":bank}
)

interface_protocol = (
    ("hdf5", "",               "--plugin_args --no_collective --compression gzip level=5"),
    ("conduit", "json",        "--plugin_args --protocol json"),
    ("conduit", "conduit_bin", "--plugin_args --protocol conduit_bin"),
    ("conduit", "hdf5",        "--plugin_args --protocol hdf5 --option hdf5/chunking/compression/method gzip --option hdf5/chunking/compression/level 5"),
    ("conduit", "adios",       "--plugin_args --protocol adios --option adios/write/transform zfp --option adios/write/transform_options rate=0.25")
)

parallel_mode = (
   ("MIFMAX", 1),  # 1 file per core/part
#   ("MIF", 16), # let this mean cores/16 number of files
   ("SIF", 1)
)

part_types = ("uniform", "rectilinear", "curvilinear", "unstructured", "arbitrary")

part_dims = (2, 3)

def iterate(func):
    n = 0
    for jsize in job_sizes:
        for iface in interface_protocol:
            for mode in parallel_mode:
                if func != None:
                    func(jsize, iface, mode)
                n = n + 1
    return n

def mkdir(dname):
    try:
        os.mkdir(dname)
    except:
        pass

class JobMaker(object):
    def __init__(self, p, mio):
        super(JobMaker, self).__init__()
        self.path = p
        self.macsio = mio

    def count_jobs(self):
        n = iterate(None)  
        print n, " jobs"

    def make_job(self, jsize, iface, mode):
        # Make the job.
        casename = "%s_%d_%04d_%s_%s" % (mode[0], mode[1], jsize["cores"], iface[0], iface[1])
        if casename[-1:] == "_":
            casename = casename[:-1]
        jobdir = os.path.join(self.path, casename)

        mkdir(jobdir)
        jobname = os.path.join(jobdir, casename+".sh")
        print jobname
        print jsize
        f = open(jobname, "wt")
        f.write("#!/bin/sh\n")
        f.write("#SBATCH -N %d\n" % (jsize["cores"] / jsize["ppn"]))
        f.write("#SBATCH -p %s\n" % jsize["queue"])
        f.write("#SBATCH -J %s\n" % casename)
	if iface[1] == "json":
            f.write("#SBATCH -t 120\n")
	else:
	    t = jsize["time"]
	    if mode[0] == "SIF":
	        t = t*2
            f.write("#SBATCH -t %d\n" % t)
        f.write("#SBATCH -A %s\n" % jsize["bank"])
        f.write("#SBATCH -o %s\n" % (jobname[:-3]+".out"))
        f.write("#SBATCH -e %s\n" % (jobname[:-3]+".err"))
        f.write("#SBATCH --license=lscratchv\n")
        f.write("\n")
        f.write("eval `/usr/bin/modulecmd sh load intel/16.0`\n")
	f.write("eval `/usr/bin/modulecmd sh load mvapich2-intel-ofa/1.7`\n")
	f.write("export LD_LIBRARY_PATH=/opt/mvapich2-intel-ofa-1.7/lib:$LD_LIBRARY_PATH\n")
        f.write("\n")
        f.write("cd %s\n" % jobdir)
        f.write("\n")
	
        if mode[0] == "MIFMAX":
            nfiles = jsize["cores"]
        elif mode[0] == "MIF":
            nfiles = jsize["cores"] / mode[1]
        else:
            nfiles = 1

        # Iterate over parttype, partdim here so we have fewer jobs with more steps.
        f.write("# Make directories \n")
        for parttype in part_types:
            for partdim in part_dims:
                f.write("mkdir %s_%d\n" % (parttype, partdim))

        nsteps = 10
	if jsize["cores"] > 256:
	    nsteps = 1

        f.write("# Run job steps\n")
        for parttype in part_types:
            for partdim in part_dims:
                if partdim == 2:
                    dims = (1000,1000,0)
                else:
                    dims = (100,100,100)

                f.write("echo \"************** Running %s for %s %d **************\"\n" % (casename, parttype, partdim))
                f.write("cd %s/%s_%d\n" % (jobdir, parttype, partdim))
                f.write("date\n")
		f.write("if test -e macsio-timings.log ; then\n")
		f.write("    echo \"Skipping %s %s %d because it already ran.\"\n" % (casename,parttype,partdim))
		f.write("else\n")
                f.write("    srun -n %d %s --num_dumps %d --interface %s --part_type %s --part_dim %d --part_mesh_dims %d %d %d --parallel_file_mode %s %d %s\n" % (jsize["cores"], self.macsio, nsteps, iface[0], parttype, partdim, dims[0], dims[1], dims[2], mode[0], nfiles, iface[2]))
                f.write("fi\n\n")
        f.close()

    def make_jobs(self):
        n = iterate(self.make_job)  

    def get_file_size(self, filename):
        s = 0
	try:
	    sr = os.stat(filename)
	    s = sr.st_size
	except:
	    s = 0
        return s

    def gather_file_sizes(self, jobdir, exts):
        fs = 0
	return 0
        files = os.listdir(jobdir)
	for f in files:
	    #print f
	    for e in exts:
	        elen = len(e)
		if f[-elen:] == e:
		    #print "match ", f, e
		    fs = fs + self.get_file_size(os.path.join(jobdir, f))
        return fs

    def gather_data_one_case(self, jsize, iface, mode, parttype, partdim, keys):
        ok = 1
	d = {}
	# Make the job.
        casename = "%s_%d_%04d_%s_%s" % (mode[0], mode[1], jsize["cores"], iface[0], iface[1])
        if casename[-1:] == "_":
            casename = casename[:-1]
	pname = "%s_%d" % (parttype, partdim)
        jobdir = os.path.join(self.path, casename, pname)
        if mode[0] == "MIFMAX":
            nfiles = jsize["cores"]
        elif mode[0] == "MIF":
            nfiles = jsize["cores"] / mode[1]
        else:
            nfiles = 1

        timings_log = os.path.join(jobdir, "macsio-timings.log")
	try:
	    lines = open(timings_log, "rt").readlines()
            print timings_log
	    idx = string.find(lines[4], "TOT=")
	    if idx != -1:
	        s = lines[4][idx+4:]
		s = string.replace(s, "CNT=", "")
		s = string.replace(s, "MIN=", "")
		s = string.replace(s, "):", ",")
		s = string.replace(s, "(", ",")
		s = string.replace(s, "AVG=", "")
		s = string.replace(s, "MAX=", "")
		s = string.replace(s, "DEV=", "")
                s = string.replace(s, "\n", "")
		
		if partdim == 2:
                    dims = (1000,1000,0)
                else:
                    dims = (100,100,100)

                command = "srun -n %d %s --interface %s --part_type %s --part_dim %d --part_mesh_dims %d %d %d --parallel_file_mode %s %d %s" % (jsize["cores"], self.macsio, iface[0], parttype, partdim, dims[0], dims[1], dims[2], mode[0], nfiles, iface[2])
                fs = self.gather_file_sizes(jobdir, [".h5", ".root"])
		if iface[1] == "":
		    iname = iface[0]
		else:
		    iname = iface[0] + "/" + iface[1]

                # Try and save the data as a dictionary
                d["cores"] = jsize["cores"]
		d["interface"] = iname
		d["mode"] = mode[0]
		d["parttype"] = parttype
		d["partdim"] = partdim
		d["command"] = command
                fvalues = [float(string.replace(x,":","")) for x in string.split(s, ",")]
		i = 0
		for k in keys:
                    d[k] = fvalues[i]
		    i = i + 1

                self.dataout.write("%s,%d,%s,%s,%d,%s,%d,%s" % (mode[0],partdim, parttype,iname, jsize["cores"], s, fs, command))
	except:
            print "CANNOT read ", timings_log
            d = {}
	    ok = 0

        m_log = os.path.join(jobdir, "macsio-log.log")
        if ok:
	    try:
	        ok = 0
	        key = "Summed  BW:"
		lkey = len(key)
	        lines = open(m_log, "rt").readlines()
		for line in lines:
		    idx = string.find(line, key)
		    if idx != -1:
		        s = line[idx+lkey:]
			idx = string.find(s, ":")
			s = s[:idx]
			#print s
			bw = 0.
			idx = string.find(s, "Mi/sec")
			if idx != -1:
			    s = s[:idx]
			    bw = float(s)
			else:
			    idx = string.find(s, "Gi/sec")
			    if idx != -1:
			        s = s[:idx]
			        bw = float(s) * 1000.
			#print s, bw
		        self.dataout.write(",%g" % bw)
			d["SummedBW"] = bw
		        ok = 1
		        break
		self.dataout.write("\n")
            except:
	        ok = 0
                print "CANNOT read ", m_log

	return (ok,d)

    def gather_data_it(self):
        keys = string.split("TOT,CNT,MIN,MINp,MINrank,AVG,MAX,MAXp,MAXrank,DEV", ",")
        # iterate in a mode that is better for the spreadsheet.
        for mode in parallel_mode:
	    for partdim in part_dims:
	        for parttype in part_types:
		    for iface in interface_protocol:
                        for jsize in job_sizes:
			    ret = self.gather_data_one_case(jsize, iface, mode, parttype, partdim, keys)
			    if ret[0] == 0:
                                print "Error reading data for case: ", jsize["cores"], iface, mode, parttype, partdim
	
    def gather_data(self, filename):
        self.dataout = open(filename, "wt")
        self.dataout.write("ParallelMode,PartDims,PartType,Interface/Protocol,Cores,TOT,CNT,MIN,MINp,MINrank,AVG,MAX,MAXp,MAXrank,DEV,FileSizes,Command,SummedBW Mi/sec\n")
        self.gather_data_it()
        self.dataout.close()

    def generate_gnuplot_image(self, datafilename, commfilename, imgfilename, d, datakey, datalabel, title, resolution):
		    # Make a gnuplot data file for the bandwidth data.
		    gpdata = open(datafilename, "wt")
		    gpdata.write("# %s\n" % datafilename)
                    b = 0
		    pcomm = []
                    for iface in interface_protocol:
		        iname = "%s/%s" % (iface[0],iface[1])
			if iname[-1] == '/':
			    iname = iname[:-1]
			if len(d[iname].keys()) > 0:
			    gpdata.write("# data block (index %d)\n" % b)
			    gpdata.write("# X Y\n")
			    cores = sorted(d[iname].keys())
                            for c in cores:
			        gpdata.write("%d %g\n" % (c, d[iname][c][datakey]))
			    gpdata.write("\n\n")
			    
			    p = "'%s' index %d with linespoints ls %d title '%s'" % (datafilename, b, b+1, iname)
			    pcomm = pcomm + [p]
			    
			    b = b + 1
	            gpdata.close()

                    g = open(commfilename, "wt")
                    g.write("set style line 1 lc rgb '#0060ad' lt 1 lw 2 pt 7 ps 1.5 # -- blue\n")
                    g.write("set style line 2 lc rgb '#ad6060' lt 1 lw 2 pt 7 ps 1.5 # -- red\n")
                    g.write("set style line 3 lc rgb '#60ad60' lt 1 lw 2 pt 7 ps 1.5 # -- green\n")
                    g.write("set style line 4 lc rgb '#eeee60' lt 1 lw 2 pt 7 ps 1.5 # -- gold\n")
                    g.write("set style line 5 lc rgb '#ff8c00' lt 1 lw 2 pt 7 ps 1.5 # -- darkorange\n")
                    g.write("\n")
                    g.write("set terminal png size %d,%d\n" % resolution)
                    g.write("set output '%s'\n" % imgfilename)
                    g.write("\n")
                    g.write("set title '%s'\n" % title)
                    g.write("show title\n")
                    g.write("\n")
                    g.write("set xlabel 'Cores'\n")
                    g.write("show xlabel\n")
                    g.write("\n")
                    g.write("set ylabel '%s'\n" % datalabel)
                    g.write("show ylabel\n")
                    g.write("\n")
                    g.write("set mxtics 5\n")
                    g.write("set mytics 5\n")
                    g.write("set grid mxtics mytics\n")
                    g.write("set logscale x\n")
                    g.write("\n")
                    g.write("set key left top\n")
                    g.write("\n")
		    pc = string.join(pcomm, ", \\\n")
		    #print pcomm
		    #print pc
                    g.write("plot %s\n" % pc)
                    g.write("quit\n")
                    g.close()

                    # Generate the image
                    #cmd = ["/usr/local/bin/gnuplot", commfilename]
                    #ret = subprocess.call(cmd, shell=True)
                    ret = os.system("gnuplot %s" % commfilename)
                    print ret

    def gather_data_gnuplot(self, excelfilename, filename):
        self.dataout = open(excelfilename, "wt")
        self.dataout.write("ParallelMode,PartDims,PartType,Interface/Protocol,Cores,TOT,CNT,MIN,MINp,MINrank,AVG,MAX,MAXp,MAXrank,DEV,FileSizes,Command,SummedBW Mi/sec\n")

        html = open(filename, "wt")
        html.write("<html>\n")
        html.write("<head><title>MACSio Scaling</title></head>\n")
        html.write("<body bgcolor=\"#ffffff\">\n")
        html.write("<p>%s</p>\n" % description)

        html.write("<ul>\n")
        for mode in parallel_mode:
	    for partdim in part_dims:
	        for parttype in part_types:
                    aname = "%s %dD %s" % (mode[0], partdim, parttype)
                    alink = "%s_%dD_%s" % (mode[0], partdim, parttype)
                    html.write("<li><a href=\"#%s\">%s</a></li>\n" % (alink,aname))
        html.write("</ul>\n")

        keys = string.split("TOT,CNT,MIN,MINp,MINrank,AVG,MAX,MAXp,MAXrank,DEV", ",")
	i = 0
        # iterate in a mode that is better for gnu plot data creation. Populate a 
	# d dictionary that contains the macsio results for the selected
	# mode,partdim,parttype cases so we can plot them relative to one another.
        for mode in parallel_mode:
	    for partdim in part_dims:
	        for parttype in part_types:
		    d = {}
		    for iface in interface_protocol:
		        iname = "%s/%s" % (iface[0],iface[1])
			if iname[-1] == '/':
			    iname = iname[:-1]
		        d[iname] = {}
                        for jsize in job_sizes:
			    ret = self.gather_data_one_case(jsize, iface, mode, parttype, partdim, keys)
			    if ret[0] != 0:
#                                print "Error reading data for case: ", jsize["cores"], iface, mode, parttype, partdim
#		            else:
			        d[iname][jsize["cores"]] = ret[1]

                    # Generate some images.
                    res = (500,500)
                    bw_datafilename = "gp_bw%04d.dat" % i
                    bw_commfilename = "gp_bw%04d.bp" % i
                    bw_imgfilename = "graph_bw%04d.png" % i
                    bw_title = "%s %dD %s" % (mode[0], partdim, parttype)
                    self.generate_gnuplot_image(bw_datafilename, bw_commfilename, bw_imgfilename, d, "SummedBW", "Summed Bandwidth Mi/sec", bw_title, res)
                    i = i + 1
                    tot_datafilename = "gp_tot%04d.dat" % i
                    tot_commfilename = "gp_tot%04d.bp" % i
                    tot_imgfilename = "graph_tot%04d.png" % i
                    tot_title = "%s %dD %s" % (mode[0], partdim, parttype)
                    self.generate_gnuplot_image(tot_datafilename, tot_commfilename, tot_imgfilename, d, "TOT", "Total Reduced Time sec", tot_title, res)
                    i = i + 1

                    # Add to the table.
                    aname = "%s_%dD_%s" % (mode[0], partdim, parttype)
                    html.write("<a name=\"%s\">\n" % aname)
                    html.write("<h1>%s</h1>\n" % bw_title)

                    html.write("<img src=\"%s\"></img>\n" % bw_imgfilename)
                    html.write("<img src=\"%s\"></img>\n" % tot_imgfilename)
                    html.write("<br>\n")

                    html.write("Runs:<br><table border=\"1\">\n")
                    html.write("<tr><td><b>Interface</b></td><td><b>Cores</b></td><td><b>Command</b></td></tr>")
                    for iface in interface_protocol:
		        iname = "%s/%s" % (iface[0],iface[1])
			if iname[-1] == '/':
			    iname = iname[:-1]
			if len(d[iname].keys()) > 0:
			    cores = sorted(d[iname].keys())
                            for c in cores:
			        html.write("<tr><td>%s</td>" % d[iname][c]["interface"])
			        html.write("<td>%d</td>" % d[iname][c]["cores"])
			        html.write("<td>%s</td></tr>" % string.replace(d[iname][c]["command"], self.macsio, "macsio"))
                    html.write("</table><br>\n")

        html.write("</body>\n</html>\n")
        self.dataout.close()

def main():
    jobdir = os.path.join(os.path.abspath(os.curdir), "macsio_jobs")
#    jobdir = os.path.join(os.getenv("SCRATCH"), "macsio_jobs")
    mkdir(jobdir)

    macsio = "/g/g17/whitlocb/Development/AICR/MACSio_install/macsio"

    jobs = JobMaker(jobdir, macsio)
    jobs.count_jobs()
    jobs.make_jobs()
#    jobs.gather_data("results.dat")
#    jobs.gather_data_gnuplot("results.dat", "results.html")

main()
