#
# Class CEP2Info defines different parameters, locus names, necessary ENV variables, etc.
# related to pipeline processing
#
import os, sys
from pulp_logging import PulpLogger
import re

class CEP2Info:
	# setting different attributes
        def __init__(self):
		self.cluster_headnode=os.popen('hostname').readlines()[0].strip().split(".")[0]
                try:
		    self.user = os.environ['USER']
                except:
                    print "Env. variable 'USER' is not defined!"                    
                    raise
		try:
			self.uid  = os.environ['UID']
		except:
			self.uid = os.popen('id -u').readlines()[0]
                try:
		    self.home = os.environ['HOME']
                except:
                    print "Env. variable 'HOME' is not defined!"                    
                    raise
                self.lofarsoft = ""
		# Additional possible directory with parfiles
		self.parfile_dir_extra = "/home/kondratiev/parfiles"
		# maximum radial distance (in deg) to find pulsars in FOV
		self.fov_lba = 12. # IS fwhm for LBA
		self.fov_hba = 6. # IS fwhm for HBA
		self.fwhm_lba = 2. # CS fwhm for LBA
		self.fwhm_hba = 1. # CS fwhm for HBA
                try:
		    self.pythonpath = os.environ['PYTHONPATH']
                except: pass
		self.current_node = "" # current node
		self.current_dir = "" # current dir
		# set the logdir, the log-file from it will be copied to corresponding _red dir
		self.logdir = ""
		# ** ONLY for auto mode! ** the logfile directory for logs for specific beams, or summary
		self.local_logdir = "" # should be composed from "local_logdir_prefix" and PipeID
		self.local_logdir_prefix = "/data/scratch/lofarsys/Observation"
		# the logfile will be initialized later when we get info about the ObsID
		self.logfile = ""
		# pipeline feedback file to be used by MoM for archive ingest
		# will be initialized when info about ObsID (pipeline Id in the future) is known
		self.feedbackfile = ""
		# directory where all parset files live
		self.parset_dir="/globalhome/lofarsystem/log"
		# keyword prefix for feedback info
		self.feedback_prefix="LOFAR.ObsSW.Observation.DataProducts.Output_Pulsar_"
		# keyword prefix for feedback info about general LOG file
		self.feedback_log_prefix="LOFAR.ObsSW.Observation.Output_Pulsar_.log"
		# keyword for showing number of DataProducts in the Feebback file for MoM
		self.feedback_nrofoutputs_keyword="LOFAR.ObsSW.Observation.DataProducts.nrOfOutput_Pulsar_"
		# directory with raw data (can be changed by user from command line)
		self.rawdir = "/data"
		# introduced specifically for CEP4 as the raw data are stored not under rawdir/ObsID
		self.rawdir_suffix_specificator=""
		# prefix for default directory with processed data
		self.processed_dir_prefix="LOFAR_PULSAR_ARCHIVE"
		self.processed_dir_root = "/data" # in reality it can be /data1/LOFAR_PUL..., /data2/LOFAR_PUL..
		# Processing directory (used in auto mode), not for summaries
		self.processed_dir = ""
		# data directory prefix on hoover nodes
		self.hoover_data_dir="/cep2"
		# dictionary for all nodes
		self.cexec_nodes={}
		# list of alive processing nodes
		self.alive_nodes = []
		# list of down nodes
		self.down_nodes = []	
		# Pipeline ID (when run in auto mode), the same as ObsID in manual mode
		self.pipeid = ""
		#
		# SLURM related
		#
		# extra shell calls for Slurm (not needed for Dragnet though)
		self.start_shell="/bin/sh -c \""
		self.end_shell="\""
		#
		self.srun_general_opts="--exclusive -n 1"
		# extra options
		self.slurm_extra_opts=""
		# extra options for summary nodes
		self.slurm_summaries_extra_opts=""
		# job id
		self.slurm_jobid=""
		# job name
		self.slurm_jobname="Pulp"
		#
		# DOCKER related
		#
		self.docker_common_opts=""
		self.docker_cmd_prefix=""

		# settings for CEP2
		if self.cluster_headnode[:5] == "lhn00":
			self.ncores = 24 # number of cores in one locus nodes. Can be used to limit a number of simultaneous processes
                        try:
			    self.lofarsoft = os.environ['LOFARSOFT']
                        except:
                            print "Env. variable 'LOFARSOFT' is not defined!"
                            raise
			# Directory with existing par-files
			self.parfile_dir = self.lofarsoft + "/release/share/pulsar/data/parfile"
			# db file from Psrcat
			self.psrcatdb = self.lofarsoft + "/release/share/pulsar/data/psrcat.db"
			# Puma2 header file template for bf2puma2
			self.puma2header = self.lofarsoft + "/release/share/pulsar/data/header.puma2"
			# full list of nodes and its cexec corresponding table
			self.locus_nodes=["locus%03d" % (num+1) for num in range(100)]
			self.hoover_nodes=["locus101", "locus102"]   # first is used to process CS data (if files per beam distributed over many nodes), second - to process IS data
#			self.hoover_nodes=["locus101"]   # is used to process data (both CS and IS), if files per beam distributed over many nodes
			for num in range(100): # we have 100 locus nodes
				key="locus%03d" % (num+1)
				self.cexec_nodes[key] = "locus:%d" % (num)
			# adding hoover nodes as well
			self.cexec_nodes["locus101"] = "hoover:0"
			self.cexec_nodes["locus102"] = "hoover:1"
			# cexec command to run. Using this mapfile makes keep mapping of the locus to be always the same
			self.cexeccmd="cexec -f /etc/c3.conf.full"
			# summary nodes
			self.summary_nodes={"CS": "locus092", "CV": "locus093", "IS": "locus094"}

		# settings for Dragnet
		elif self.cluster_headnode[:4] == "drag" or self.cluster_headnode[:3] == "drg":
			self.ncores = 16 # number of cores in one dragnet node. Can be used to limit a number of simultaneous processes
			self.lofarsoft = "/usr/local/"
			# Directory with existing par-files
			self.parfile_dir = self.lofarsoft + "etc/parfiles"
			# db file from Psrcat
			self.psrcatdb = self.lofarsoft + "bin/psrcat.db"
			# Puma2 header file template for bf2puma2
			self.puma2header = self.lofarsoft + "etc/header.puma2"
			# full list of nodes and its cexec corresponding table
			self.locus_nodes=["drg%02d" % (num+1) for num in range(23)]
			self.hoover_nodes=["dragproc"]   # first is used to process CS data (if files per beam distributed over many nodes), second - to process IS data
			for num in range(23): # we have 23 dragnet nodes
				key="drg%02d" % (num+1)
				self.cexec_nodes[key] = "drg:%d" % (num)
			# adding hoover nodes as well
			self.cexec_nodes["dragproc"] = "dragnet:23"
			# cexec command to run. Using this mapfile makes keep mapping of the locus to be always the same
			self.cexeccmd="cexec -f /usr/local/etc/c3.conf"
			# summary nodes
			self.summary_nodes={"CS": "dragproc", "CV": "dragproc", "IS": "dragproc"}
			#
			# SLURM related
			#
			# extra shell calls for Slurm (not needed for Dragnet though)
			self.start_shell="/bin/sh -c "
			self.end_shell=""
			# extra options
			self.srun_general_opts="-n 1"
			#self.slurm_extra_opts="-p proc,workers"
			self.slurm_extra_opts="-N 1 --mem-per-cpu=8192"
			# extra options for summary nodes
			#self.slurm_summaries_extra_opts="-p proc"
			self.slurm_summaries_extra_opts="-N 1 -w dragproc --mem-per-cpu=8192"

		# settings for CEP3
		elif self.cluster_headnode[:5] == "lhd00":
			self.ncores = 40 # number of cores in one CEP3 node. Can be used to limit a number of simultaneous processes
                        try:
			    self.lofarsoft = os.environ['LOFARSOFT']
                        except:
                            print "Env. variable 'LOFARSOFT' is not defined!"
                            raise
			# Directory with existing par-files
			self.parfile_dir = self.lofarsoft + "/release/share/pulsar/data/parfile"
			# db file from Psrcat
			self.psrcatdb = self.lofarsoft + "/release/share/pulsar/data/psrcat.db"
			# Puma2 header file template for bf2puma2
			self.puma2header = self.lofarsoft + "/release/share/pulsar/data/header.puma2"
	                # prefix for default directory with processed data
			self.processed_dir_prefix="scratch/vlad/LOFAR_PULSAR_ARCHIVE"
			# full list of nodes and its cexec corresponding table
			self.locus_nodes=["lof%03d" % (num+1) for num in range(20)]
			self.hoover_nodes=["lof020"]   # first is used to process CS data (if files per beam distributed over many nodes), second - to process IS data
			for num in range(8): # we have 100 locus nodes
				key="lof%03d" % (num+1)
				self.cexec_nodes[key] = "lof:%d" % (num)
			# adding hoover nodes as well
			self.cexec_nodes["lof020"] = "lof:19"
			# cexec command to run. Using this mapfile makes keep mapping of the locus to be always the same
			self.cexeccmd="cexec"
			# summary nodes
			self.summary_nodes={"CS": "lhd002", "CV": "lhd002", "IS": "lhd002"}

		# settings for CEP4
		elif self.cluster_headnode[:5] == "head0" or self.cluster_headnode[:3] == "cpu":
			self.ncores = 20 # number of cores in one CEP4 node. Can be used to limit a number of simultaneous processes
                        try:
			    self.lofarsoft = os.environ['LOFARSOFT']
                        except:
			    self.lofarsoft = "/usr/local/"
                        # check if everything is in /usr/local or LOFARSOFT/...
                        if self.lofarsoft == "/usr/local/" or (self.lofarsoft != "/usr/local/" and \
                                not os.path.exists("%s/release/share/pulsar/data/parfile" % (self.lofarsoft,))):
			    # Directory with existing par-files
        		    self.parfile_dir = "%s/%s/parfiles" % (self.processed_dir_root, self.processed_dir_prefix)
        		    # db file from Psrcat
	        	    self.psrcatdb = self.lofarsoft + "bin/psrcat.db"
		            # Puma2 header file template for bf2puma2
			    self.puma2header = self.lofarsoft + "etc/header.puma2"
                        else: # if you have LOFARSOFT with full installation
			    # Directory with existing par-files
        		    self.parfile_dir = self.lofarsoft + "/release/share/pulsar/data/parfile"
			    # db file from Psrcat
        		    self.psrcatdb = self.lofarsoft + "/release/share/pulsar/data/psrcat.db"
			    # Puma2 header file template for bf2puma2
        		    self.puma2header = self.lofarsoft + "/release/share/pulsar/data/header.puma2"

			self.local_logdir_prefix = "/data/scratch/pipeline/Observation"
			# full list of nodes and its cexec corresponding table
			self.locus_nodes=["cpu%02d" % (num+1) for num in range(50)]
			self.hoover_nodes=[]   # first is used to process CS data (if files per beam distributed over many nodes), second - to process IS data
			for num in range(50): # we have 50 CEP4 compute nodes
				key="cpu%02d" % (num+1)
				self.cexec_nodes[key] = "cpu:%d" % (num)
			# adding hoover nodes as well
			# ...
			# cexec command to run. Using this mapfile makes keep mapping of the locus to be always the same
			self.cexeccmd="cexec"
			# summary nodes
			# on CEP4 these could be any node... I put these keywords in the dictionary for consistency, but with Slurm they will never be used
			# I've put "CEP4" as it was needed for CEP4 feedback files
			# summary nodes
			self.summary_nodes={"CS": "CEP4", "CV": "CEP4", "IS": "CEP4"}
			#	
			# introduced specifically for CEP4 as the raw data are stored not under rawdir/ObsID
			self.rawdir_suffix_specificator="*/"
			#
			# SLURM related
			#
			# extra options
			#self.slurm_extra_opts="--mem-per-cpu=8192"
			# extra options for summary nodes
			#self.slurm_summaries_extra_opts="--mem-per-cpu=8192"
			#
			# DOCKER related
			#
			#self.docker_common_opts="docker run --rm -u %s -e USER=%s -e HOME=%s -v %s/.ssh:%s/.ssh:ro -v /data:/data --net=host" % \
			self.docker_common_opts="/data/bin/docker-run-slurm.sh --rm -u %s -e USER=%s -e HOME=%s -v %s/.ssh:%s/.ssh:ro -v /data:/data --net=host" % \
				(self.uid, self.user, self.home, self.home, self.home)
			self.docker_cmd_prefix="ssh -n -tt -x localhost "
			#self.docker_cmd_prefix="ssh -t -x localhost "

		# undefined cluster
		else:
			print "Unknown cluster: %s. Add proper settings to pulp_sysinfo.py" % (self.cluster_headnode)
			raise Exception

	# checking connection to all processing nodes, to determine which ones are "alive"
	def check_connection(self, log=None):
		msg="\nChecking connection to processing nodes..."
		if log != None: log.info(msg)
		else: print msg
		# forming string with all processing nodes to check in one cexec command
		cexeclocus=self.cexec_nodes[self.locus_nodes[0]] # there is always at least one processing node
		if len(self.locus_nodes) > 1:
			for s in self.locus_nodes[1:]:
				cexeclocus += ",%s" % (self.cexec_nodes[s].split(":")[1])
		# adding hoover nodes to cexeclocus
		cexeclocus += " %s" % (self.cexec_nodes[self.hoover_nodes[0]])
		if len(self.hoover_nodes) > 1:
			for s in self.hoover_nodes[1:]:
				cexeclocus += ",%s" % (self.cexec_nodes[s].split(":")[1])
		cmd="%s %s 'date' | grep -v denied | grep -v xauth | grep -v connect | grep -v closed | grep -v Permission | grep -v The | grep -v maintenance | grep -v @ | egrep -v \'\\*\\*\\*\\*\\*\'" % (self.cexeccmd, cexeclocus)
		cexec_output=[line[:-1] for line in os.popen(cmd).readlines()]
		# finding all processing nodes that have the dir with raw data
		try:
			for l in range(len(cexec_output)):
				if re.match("^-----", cexec_output[l]) is None:
					self.alive_nodes.append(cexec_output[l-1].split(" ")[1].split("-")[0])
		except Exception:
			msg="Problem with connection to processing nodes...\nTry removing processing nodes' entries from your ~/.ssh/known_hosts file or try again later"
			if log != None: log.error(msg)
			else: print msg
			raise
		if len(self.alive_nodes) == 0:
			msg="The connection to all processing nodes is down. Try again later"
			if log != None: log.error(msg)
			else: print msg
			raise Exception

	# print nodes that are down
	def print_down_nodes(self, log=None):
		all_nodes = self.cexec_nodes.keys()
		self.down_nodes=list(set(all_nodes)-set(all_nodes).intersection(set(self.alive_nodes)))
		msg="Nodes are down [%d]: %s" % (len(self.down_nodes), ", ".join(self.down_nodes))
		if log != None: log.info(msg)
		else: print msg

	# return list of alive nodes
	def get_alive_nodes(self):
		return self.alive_nodes

	# return $LOFARSOFT
	def get_lofarsoft(self):
		return self.lofarsoft

	# return current node
	def get_current_node(self):
		return self.current_node

	# set current node
	def set_current_node(self):
		self.current_node = os.popen('hostname').readlines()[0].strip().split(".")[0]

	# return current dir
	def get_current_dir(self):
		return self.current_dir

	# set current dir
	def set_current_dir(self):
		self.current_dir = os.popen('pwd').readlines()[0][:-1]

	# return logfile
	def get_logfile(self):
		return self.logfile

	# return logdir
	def get_logdir(self):
		return self.logdir

	# return local logdir
	def get_local_logdir(self, cmdline):
		if self.local_logdir == "": return self.logdir 
		if cmdline.opts.is_auto and cmdline.opts.is_local: return self.local_logdir
		return self.logdir

	# return feedback file
	def get_feedbackfile(self):
		return self.feedbackfile

	# set the logdir
	def set_logdir(self, directory):
		self.logdir = directory
		cmd="mkdir -p %s" % (self.logdir)	
		os.system(cmd)

	# set local logdir
	def set_local_logdir(self, pipeid):
		self.local_logdir = "%s%s" % (self.local_logdir_prefix, pipeid)
		cmd="mkdir -p %s" % (self.local_logdir)	
		os.system(cmd)

	# set the logfile
	def set_logfile(self, f):
		self.logfile = f

	# set the feedback file
	def set_feedbackfile(self, f):
		self.feedbackfile = f

	# set the pipeline ID
	def set_pipeid(self, id):
		self.pipeid = id

	# get the pipeline ID
	def get_pipeid(self):
		return self.pipeid

	# set SLURM jobid
	def set_slurm_jobid(self):
		try:
			self.slurm_jobid = os.environ['SLURM_JOB_ID']
		except: raise

	# set SLURM job name
	def set_slurm_jobname(self):
		try:
			jobname = os.environ['SLURM_JOB_NAME']
			self.slurm_jobname = jobname
		#except: raise
		except: pass

	# set processing directory
	def set_processed_dir(self, procdir):
		self.processed_dir = procdir

	# print info of all set attributes
	def print_info(self, cmdline, log=None):
		if log != None:
			log.info("")
			log.info("USER = %s" % (self.user))
			log.info("Current node = %s" % (self.current_node))
			log.info("Current directory = %s" % (self.current_dir))
			log.info("LOFARSOFT = %s" % (self.lofarsoft))
			if cmdline.opts.is_debug:
				log.info("PYTHONPATH = %s" % (self.pythonpath))
			log.info("")
		
