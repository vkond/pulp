# Class CMDLine collects info about user input from the command line
#
import optparse as opt
import os, sys
import numpy as np
from pulp_sysinfo import CEP2Info
from pulp_logging import PulpLogger
import time
import re
import random, string
import logging

# checks if given pulsar is good ones, i.e. either in ATNF catalog
# or par-file do exist
# return value is True (psr is good), or False (psr is bad)
def check_pulsars(psr, cmdline, cep2, log=None):
	if psr not in cmdline.catpsrs:
		msg="Warning: Pulsar %s is not in the ATNF pulsar catalog! Checking for par-file..." % (psr)
		if log != None: log.warning(msg)
		# checking if par-file exist
		if cmdline.opts.parfile != "":
			if os.path.exists(cmdline.opts.parfile): 
				msg="Found parfile '%s'. Continue..." % (cmdline.opts.parfile)
				if log != None: log.info(msg)
				return True
			else: 
				msg="Can't find user parfile '%s'. Exiting..." % (cmdline.opts.parfile)
				if log != None: log.error(msg)
				else: print msg
				raise Exception
		# check parfile when parfile name does not have 'B|J'
		parfile="%s/%s.par" % (cep2.parfile_dir, re.sub(r'^[BJ]', '', psr))
		if not os.path.exists(parfile):
			# checking another parfile name
			parfile="%s/%s.par" % (cep2.parfile_dir, psr)
			if not os.path.exists(parfile): 
				# checking extra directory as well
				parfile="%s/%s.par" % (cep2.parfile_dir_extra, re.sub(r'^[BJ]', '', psr))
				if not os.path.exists(parfile):
					# checking another parfile name in the extra directory
					parfile="%s/%s.par" % (cep2.parfile_dir_extra, psr)
					if not os.path.exists(parfile):
						return False
		msg="Found parfile '%s'. Continue..." % (parfile)
		if log != None: log.info(msg)
	return True

class CMDLine:
	# parsing a command line
        def __init__(self, version="", wrapper=None):
		self.prg = sys.argv[0]
		self.options = sys.argv[1:]  # storing original cmd line
		self.version = version
		self.psrs = []  # list of pulsars to fold
		self.beams = [] # list of beams to process
		self.user_beams = [] # list of User beams to process
		self.user_excluded_beams = [] # list of User excluded beams
		self.catpsrs = [] # list of pulsars from ATNF catalog
		self.ras = self.decs = self.s400 = self.DMs = self.P0s = [] # lists of RA, DEC, S400, DM, and P0 of catalog pulsars
        	self.usage = "Usage: %prog <--id ObsID> [-h|--help] [OPTIONS]"
        	self.cmd = opt.OptionParser(self.usage, version="%prog " + self.version)
        	self.cmd.add_option('--id', '--obsid', dest='obsid', metavar='ObsID',
                           help="specify the Observation ID (i.e. L30251). This option is required", default="", type='str')
        	self.cmd.add_option('-p', '-P', '--pulsar', dest='psr', metavar='PSRS|word',
                           help="specify the Pulsar Name or comma-separated list of Pulsars for folding (w/o spaces) or \
give one of the 5 special words: \"parset\" or \"meta\" - to take pulsar name from the source field for each SAP \
separately from the parset file, or \"sapfind\", \"sapfind3\" to find the best (3 best) pulsars in FOV \
of the particular SAP, or \"tabfind\" to find the brightest pulsar for each TAB individually, or \
\"tabfind+\" to first get pulsar from the parset if pulsar name their is legitimate, then get the brightest \
one in the SAP (same as \"sapfind\"), and get another pulsar from the TAB (same as \"tabfind\"). \
If no pulsars are given and no special words used, then pipeline will try to take source names from \
parset file first, and then look for the best pulsars in SAP's FOV (same as \"sapfind\"). \
Word 'NONE' as a pulsar name is ignored", default="", type='str')
        	self.cmd.add_option('-o', '-O', '--output', dest='outdir', metavar='DIR',
                           help="specify the Output Processing Location relative to /data*/LOFAR_PULSAR_ARCHIVE_*. \
Default is corresponding *_red or *_redIS directory", default="", type='str')
        	self.cmd.add_option('--par', '--parfile', '--eph', dest='parfile', metavar='PARFILE',
                           help="specify the parfile for one pulsar to fold. Pulsar name should be given explicitely using --pulsar option \
and only one pulsar name should be given for --par option to work", default="", type='str')
        	self.cmd.add_option('--parfile-dir', dest='parfiledir', metavar='DIR',
                           help="specify the directory where to search for parfiles. By default, parfiles are looked for in the build directory \
unless parfile is given with the --par option", default="", type='str')
                self.cmd.add_option('--auto', action="store_true", dest='is_auto',
                           help="optional parameter to indicate this is an automated pipeline run from the LOFAR pipeline framework. DO NOT USE IT \
when running the pipeline manually from the terminal!", default=False)
        	self.cmd.add_option('--nodecode', action="store_true", dest='is_nodecode',
                           help="optional parameter to skip decoding the data (2bf2fits/bf2puma2)", default=False)
        	self.cmd.add_option('--norfi', action="store_true", dest='is_norfi',
                           help="optional parameter to skip rfifind and subdyn.py RFI checker and/or Coast_Guard's clean.py (used to be paz)", default=False)
        	self.cmd.add_option('--nofold', action="store_true", dest='is_nofold',
                           help="optional parameter to turn off folding of data (prepfold is not run)", default=False)
        	self.cmd.add_option('--nopdmp', action="store_true", dest='is_nopdmp',
                           help="turn off running pdmp in the pipeline", default=False)
        	self.cmd.add_option('--summary', action="store_true", dest='is_summary', 
                           help="making only summaries on already processed data", default=False)
        	self.cmd.add_option('--nosummary', action="store_true", dest='is_nosummary', 
                           help="do not make summaries of already processed data", default=False)
        	self.cmd.add_option('--feedback', action="store_true", dest='is_feedback', 
                           help="making only feedback files on already processed data", default=False)
        	self.cmd.add_option('--plots-only', action="store_true", dest='is_plots_only', 
                           help="creating diagnostic plots only on processing nodes assuming the data required for plots is there already", default=False)
        	self.cmd.add_option('--single-pulse', action="store_true", dest='is_single_pulse', 
                           help="running single-pulse analysis in addition to folding a profile (implemented only for CS/IS \
for PRESTO part only, and for CV data using digifil approach", default=False)
        	self.cmd.add_option('--rrats', action="store_true", dest='is_rrats', 
                           help="running prepsubband for a range of DMs (default - 1000 DM trials +/-5 around nominal DM of the pulsar) and prepdata for DM=0 \
followed by single_pulse_search.py. Prepsubband is run using -nsub option with the highest possible value of subbands, the highest common denominator of the number \
of channels that is <=1024. Use --prepsubband-extra-opts to set different DM range/step and/or nsub", default=False)
        	self.cmd.add_option('--beams', dest='beam_str', metavar='[^]SAP#:TAB#[,SAP#:TAB#,...]',
                           help="user-specified beams to process separated by commas and written as station beam number, colon, \
TA beam number, with no spaces. The argument can have leading hat character '^' to indicate that \
specified beams are to be excluded from processing", default="", type='str')
        	self.cmd.add_option('--stokes', dest='stokes', metavar='STOKES#',
                           help="only process specific Stokes file. Only relevant for IQUV data. Possible values are 0,1,2,3. Default - process all", default=-1, type='int')
        	self.cmd.add_option('--noIS', action="store_true", dest='is_noIS',
                           help="optional parameter to turn off processing of Incoherent sum (IS) data", default=False)
        	self.cmd.add_option('--noCS', action="store_true", dest='is_noCS',
                           help="optional parameter to turn off processing of Coherent sum (CS) data", default=False)
        	self.cmd.add_option('--noCV', action="store_true", dest='is_noCV',
                           help="optional parameter to turn off processing of Complex voltages (CV) data", default=False)
        	self.cmd.add_option('--noFE', action="store_true", dest='is_noFE',
                           help="optional parameter to turn off processing of Fly's Eye (FE) data", default=False)
        	self.cmd.add_option('--del', '--delete', action="store_true", dest='is_delete',
                           help="optional parameter to delete the previous entire output processing location if it exists. \
Otherwise, the new results will be overwritten/added to existing directory", default=False)
        	self.cmd.add_option('--parset', dest='parset', metavar='FILE',
                           help="specify explicitely the input parameter file (parset file). By default, it will be looked for in standard system directory", default="", type='str')
        	self.cmd.add_option('--raw', dest='rawdir', metavar='RAWDIR',
                           help="specify the location of input raw data. Directory structure is assumed as RAWDIR/<ObsID>.", default="/data", type='str')
        	self.cmd.add_option('--locate-rawdata', action="store_true", dest='is_locate_rawdata',
                           help="search for input raw data in all alive nodes instead of using the list of nodes from the parset file", default=False)
        	self.cmd.add_option('--skip-check-rawdata', action="store_true", dest='is_skip_check_rawdata',
                           help="skip checking at all if rawdata are present on the processing nodes", default=False)
        	self.cmd.add_option('--log-append', action="store_true", dest='is_log_append',
                           help="optional parameter to append log output to already existent log files. Default is overwrite", default=False)
        	self.cmd.add_option('--nthreads', dest='nthreads', metavar='#THREADS',
                           help="number of threads for all dspsr calls. Default: %default", default=2, type='int')
        	self.cmd.add_option('--tsubint', dest='tsubint', metavar='SECS',
                           help="set the length of each subintegration to SECS. Default is 60 secs for CS/IS and 5 secs for CV mode", default=-1, type='int')
        	self.cmd.add_option('--dspsr-extra-opts', dest='dspsr_extra_opts', metavar='STRING',
                           help="specify extra additional options for Dspsr command'", default="", type='str')
        	self.cmd.add_option('--digifil-extra-opts', dest='digifil_extra_opts', metavar='STRING',
                           help="specify extra additional options for digifil command", default="", type='str')
        	self.cmd.add_option('--prepdata-extra-opts', dest='prepdata_extra_opts', metavar='STRING',
                           help="specify extra additional options for prepdata command", default="", type='str')
        	self.cmd.add_option('--raw-to-8bit', action="store_true", dest='is_raw_to_8bit',
                           help="convert raw 32-bit data to 8 bits in addition to other processing", default=False)
        	self.cmd.add_option('--digitize-sigma', dest='digitize_sigma', metavar='FLOAT',
                           help="clip raw data above this threshold (in sigmas) for conversion raw data from 32 to 8 bits. Default: %default", default=5.0, type='float')
        	self.cmd.add_option('--hoover', action="store_false", dest='is_nohoover',
                           help="use hoover node locus101 for processing, rather than rsync data to target processing nodes", default=True)
		self.cmd.add_option('--all-parts-at-once', action="store_true", dest='is_all_parts_at_once',
			   help="process all parts in parallel when they all are on one node", default=False)
        	self.cmd.add_option('--first-frequency-split-coh', dest='first_freq_splitCS', metavar='SPLIT#',
                           help="start processing from this frequency split for CS beams. It works only for processing with DAL support. Default: %default", default=0, type='int')
        	self.cmd.add_option('--first-frequency-split-incoh', dest='first_freq_splitIS', metavar='SPLIT#',
                           help="start processing from this frequency split for IS beams. It works only for processing with DAL support. Default: %default", default=0, type='int')
        	self.cmd.add_option('--nsplits-coh', dest='nsplitsCS', metavar='#SPLITS',
                           help="only process the #SPLITS splits for CS beams starting from SPLIT# determined by --first-frequency-split-coh. \
It works only for processing with DAL support. Default: all splits", default=-1, type='int')
        	self.cmd.add_option('--nsplits-incoh', dest='nsplitsIS', metavar='#SPLITS',
                           help="only process the #SPLITS splits for IS beams starting from SPLIT# determined by --first-frequency-split-incoh. \
It works only for processing with DAL support. Default: all splits", default=-1, type='int')
        	self.cmd.add_option('--fwhm-CS', dest='fwhm_CS', metavar='FWHM (deg)',
                           help="set the full-width at half maximum (in degrees) for CS beams. Default is 1 deg for HBA and 2 deg for LBA that \
corresponds roughly to the beam sizes of Superterp at 120 and 60 MHz", default=-1., type='float')
        	self.cmd.add_option('--fwhm-IS', dest='fwhm_IS', metavar='FWHM (deg)',
                           help="set the full-width at half maximum (in degrees) for IS beams. Default is 6 deg for HBA and 12 deg for LBA that \
corresponds roughly to the beam sizes of Superterp at 120 and 60 MHz", default=-1., type='float')
        	self.cmd.add_option('--bgp', action="store_false", dest='is_cobalt',
                           help="use BG/P oriented pipeline, i.e. using parset file. Default is Cobalt oriented pipeline, i.e. read metadata from .h5 files", default=True)
        	self.cmd.add_option('--slurm', action="store_true", dest='is_slurm',
                           help="use Slurm to connect with other nodes, otherwise SSH is used", default=False)
        	self.cmd.add_option('--globalfs', action="store_true", dest='is_globalfs',
                           help="all nodes share the same global disk space without need to copy the data between nodes", default=False)
        	self.cmd.add_option('--docker', action="store_true", dest='is_docker',
                           help="pipeline is executed from within Docker container", default=False)
        	self.cmd.add_option('--docker-container', dest='docker_container', metavar='STRING',
                           help="specify the name of the docker image and its tag (default: %default)", default="pulp:latest", type='str')
        	self.cmd.add_option('--debug', action="store_true", dest='is_debug',
                           help="optional for testing: turns on debug level logging in Python and intermediate data files are not deleted", default=False)
        	self.cmd.add_option('-q', '--quiet', action="store_true", dest='is_quiet',
                           help="optional parameter to turn off user's warnings and waiting time of 10 seconds in the beginning", default=False)
        	self.cmd.add_option('--noinit', action="store_true", dest='is_noinit',
                           help="do not check for down nodes and available input raw data. Observation config is read from saved file \
rather then is initialized using parset file (mostly for _internal_ use only)", default=False)
        	self.cmd.add_option('--local', action="store_true", dest='is_local', 
                           help="to process the data locally on current processing node for one beam only. Should only be used together \
with --beams option and only the first beam will be used if there are several specified in --beams \
(mostly for _internal_ use only)", default=False)
        	self.cmd.add_option('--confdir', dest='confdir', metavar='CONFIGDIR',
                           help="specify the location of configuration directory. By default, ~/.pulp/<ObsID>", default="", type='str')
		# adding CS/IS/FE extra options
	        self.groupCS = opt.OptionGroup(self.cmd, "CS/IS/FE extra options")
        	self.groupCS.add_option('--skip-subdyn', action="store_true", dest='is_skip_subdyn',
                           help="optional parameter to skip subdyn.py only", default=False)
        	self.groupCS.add_option('--subdyn-time-average', dest='subdyn_time_average', metavar='TIME (s)',
                           help="for subdyn.py set the averaging time in seconds. Default = %default s", default=0.5, type='float')
        	self.groupCS.add_option('--skip-dspsr', action="store_true", dest='is_skip_dspsr',
                           help="optional parameter to turn off running dspsr part of the pipeline when running without DAL support (including pdmp and creation of corresponding plots)", default=False)
        	self.groupCS.add_option('--skip-prepfold', action="store_true", dest='is_skip_prepfold',
                           help="optional parameter to turn off running prepfold part of the pipeline", default=False)
        	self.groupCS.add_option('--with-dal', action="store_true", dest='is_with_dal',
                           help="use dspsr directly to read raw data instead of 2bf2fits. No PRESTO routines though...", default=False)
        	self.groupCS.add_option('--2bf2fits-extra-opts', dest='bf2fits_extra_opts', metavar='STRING',
                           help="specify extra additional options for 2bf2fits command", default="", type='str')
        	self.groupCS.add_option('--decode-nblocks', dest='decode_nblocks', metavar='#BLOCKS',
                           help="read #BLOCKS at once in 2bf2fits. Same as -A option in 2bf2fits. Default: %default", default=100, type='int')
        	self.groupCS.add_option('--decode-sigma', dest='decode_sigma', metavar='INTEGER',
                           help="sigma limit value used for packing in 2bf2fits. Same as -sigma option in 2bf2fits. Default: %default", default=3, type='int')
        	self.groupCS.add_option('--prepfold-extra-opts', dest='prepfold_extra_opts', metavar='STRING',
                           help="specify extra additional options for Prepfold command", default="", type='str')
        	self.groupCS.add_option('--prepsubband-extra-opts', dest='prepsubband_extra_opts', metavar='STRING',
                           help="specify extra additional options for Prepsubband command when --rrats is used", default="", type='str')
        	self.groupCS.add_option('--rfifind-extra-opts', dest='rfifind_extra_opts', metavar='STRING',
                           help="specify extra additional options for rfifind command", default="", type='str')
		self.cmd.add_option_group(self.groupCS)
		# adding CV extra options
	        self.groupCV = opt.OptionGroup(self.cmd, "Complex voltage (CV) extra options")
        	self.groupCV.add_option('--nodal', action="store_true", dest='is_nodal',
                           help="use bf2puma2 to read raw data instead of using dspsr to read *.h5 files directly", default=False)
        	self.groupCV.add_option('--skip-rmfit', action="store_true", dest='is_skip_rmfit',
                           help="skip running rmfit program", default=False)
        	self.groupCV.add_option('--hist-cutoff', dest='hist_cutoff', metavar='FRACTION',
                           help="clip FRACTION off the edges of the samples histogram. Be noted, it eliminates spiky RFI, but may also \
clip bright pulsar pulses. Default: %default (no clipping)", default=0.02, type='float')
        	self.groupCV.add_option('--nblocks', dest='nblocks', metavar='#BLOCKS',
                           help="only read the first #BLOCKS blocks. Default: all blocks", default=-1, type='int')
        	self.groupCV.add_option('--all-for-scaling', action="store_true", dest='is_all_times',
                           help="normalize the data based on entire data set. Otherwise, the scaling is updated after every data block", default=False)
        	self.groupCV.add_option('--write-ascii', action="store_true", dest='is_write_ascii',
                           help="write out also ascii files (.rv) containing complex values", default=False)
		self.cmd.add_option_group(self.groupCV)
        
		# reading cmd options
		(self.opts, self.args) = self.cmd.parse_args()

		# check if any input parameters are given
		if len(sys.argv[1:]) == 0:
			self.cmd.print_usage()
			quit(0)

		# Here is the place to call function that parses wrapper.pulsar_parms object to get pipeline tunable parameters
		# and sets proper self.opts and self.options
		if self.opts.is_auto and not self.opts.is_local:
			self.parse_tunable_parameters(wrapper)	

		# we need to run it here separately (in addition to sysinfo), because CEP2 object gets created after CMDLine
		# and we need to know where we running the code, on CEP4 or not
		cluster_headnode=os.popen('hostname').readlines()[0].strip().split(".")[0]
		if cluster_headnode[:5] == "head0" or cluster_headnode[:3] == "cpu": # CEP4
			# if option has spaces then we protect it back slashes
			for ii in range(len(self.options)):
				self.options[ii].replace(" ", "\\\ ")
		else: # NOT CEP4
			# if option has spaces then we protect it with quotes
			for ii in range(len(self.options)):
				if ' ' in self.options[ii]: 
					if "=" in self.options[ii]:
						# for long options with usage of "=" we should put quotes after this "="
						self.options[ii] = self.options[ii].split("=")[0] + "=" + '"' + "=".join(self.options[ii].split("=")[1:]) + '"'
					else:
						self.options[ii] = '"' + self.options[ii] + '"'


	# function that is used only in automatic mode to parse tunable parameters in the pipeline parset provided by the pipeline parset file
	def parse_tunable_parameters(self, wrapper=None):
		# check if wrapper.pulsar_params object is None or not. If it is None then raise an exception
		if wrapper.pulsar_parms == None:
			raise Exception("The object wrapper.pulsar_params is None!")
		else:
			# general parameters
			self.opts.psr=wrapper.pulsar_parms.getString("pulsar")
			if self.opts.psr != "": self.options.extend(["-p", "%s" % self.opts.psr])
			self.opts.is_single_pulse = wrapper.pulsar_parms.getBool("single_pulse")
			if self.opts.is_single_pulse: self.options.append("--single-pulse")
			self.opts.is_raw_to_8bit = wrapper.pulsar_parms.getBool("raw_to_8bit")
			if self.opts.is_raw_to_8bit: self.options.append("--raw-to-8bit")
			self.opts.dspsr_extra_opts = wrapper.pulsar_parms.getString("dspsr_extra_opts")
			if self.opts.dspsr_extra_opts != "": self.options.extend(["--dspsr-extra-opts", "%s" % self.opts.dspsr_extra_opts])
			self.opts.prepdata_extra_opts = wrapper.pulsar_parms.getString("prepdata_extra_opts")
			if self.opts.prepdata_extra_opts != "": self.options.extend(["--prepdata-extra-opts", "%s" % self.opts.prepdata_extra_opts])
			self.opts.digitize_sigma = wrapper.pulsar_parms.getFloat("8bit_conversion_sigma")
			if self.opts.digitize_sigma != 5.0: self.options.extend(["--digitize-sigma", "%f" % self.opts.digitize_sigma])
			self.opts.tsubint = wrapper.pulsar_parms.getInt("tsubint")
			if self.opts.tsubint != -1: self.options.extend(["--tsubint", "%d" % self.opts.tsubint])
			self.opts.is_norfi = wrapper.pulsar_parms.getBool("norfi")
			if self.opts.is_norfi: self.options.append("--norfi")
			self.opts.is_nofold = wrapper.pulsar_parms.getBool("nofold")
			if self.opts.is_nofold: self.options.append("--nofold")
			self.opts.is_nopdmp = wrapper.pulsar_parms.getBool("nopdmp")
			if self.opts.is_nopdmp: self.options.append("--nopdmp")
			self.opts.is_skip_dspsr = wrapper.pulsar_parms.getBool("skip_dspsr")
			if self.opts.is_skip_dspsr: self.options.append("--skip-dspsr")
			# CS/IS parameters
			self.opts.is_rrats = wrapper.pulsar_parms.getBool("rrats")
			if self.opts.is_rrats: self.options.append("--rrats")
			self.opts.bf2fits_extra_opts = wrapper.pulsar_parms.getString("2bf2fits_extra_opts")
			if self.opts.bf2fits_extra_opts != "": self.options.extend(["--2bf2fits-extra-opts", "%s" % self.opts.bf2fits_extra_opts])
			self.opts.decode_sigma = wrapper.pulsar_parms.getInt("decode_sigma")
			if self.opts.decode_sigma != 3: self.options.extend(["--decode-sigma", "%d" % self.opts.decode_sigma])
			self.opts.decode_nblocks = wrapper.pulsar_parms.getInt("decode_nblocks")
			if self.opts.decode_nblocks != 100: self.options.extend(["--decode-nblocks", "%d" % self.opts.decode_nblocks])
			self.opts.rfifind_extra_opts = wrapper.pulsar_parms.getString("rfifind_extra_opts")
			if self.opts.rfifind_extra_opts != "": self.options.extend(["--rfifind-extra-opts", "%s" % self.opts.rfifind_extra_opts])
			self.opts.prepfold_extra_opts = wrapper.pulsar_parms.getString("prepfold_extra_opts")
			if self.opts.prepfold_extra_opts != "": self.options.extend(["--prepfold-extra-opts", "%s" % self.opts.prepfold_extra_opts])
			self.opts.prepsubband_extra_opts = wrapper.pulsar_parms.getString("prepsubband_extra_opts")
			if self.opts.prepsubband_extra_opts != "": self.options.extend(["--prepsubband-extra-opts", "%s" % self.opts.prepsubband_extra_opts])
			self.opts.subdyn_time_average = wrapper.pulsar_parms.getFloat("dynamic_spectrum_time_average")
			if self.opts.subdyn_time_average != 0.5: self.options.extend(["--subdyn-time-average", "%f" % self.opts.subdyn_time_average])
			self.opts.is_skip_subdyn = wrapper.pulsar_parms.getBool("skip_dynamic_spectrum")
			if self.opts.is_skip_subdyn: self.options.append("--skip-subdyn")
			self.opts.is_skip_prepfold = wrapper.pulsar_parms.getBool("skip_prepfold")
			if self.opts.is_skip_prepfold: self.options.append("--skip-prepfold")
			# CV parameters
			self.opts.digifil_extra_opts = wrapper.pulsar_parms.getString("digifil_extra_opts")
			if self.opts.digifil_extra_opts != "": self.options.extend(["--digifil-extra-opts", "%s" % self.opts.digifil_extra_opts])


	# set version of the program
	def set_version(self, version):
		self.version=version

	# prints countdown (only to terminal)
	def press_controlc(self, fr, cur):
		msg="\bPress Control-C to stop in: %s" % (" ".join(str(i) for i in range(fr, cur-1, -1)))
		print "\b" * int("%d" % (29 + 3*(fr-cur))), msg,
		sys.stdout.flush()

	# waiting for user to make a decision about "Stop or continue"
	def waiting_for_user(self, secs_to_wait, log=None):
		msg="Waiting %d seconds before starting..." % (secs_to_wait)
		if log != None: log.info(msg)
		else: print msg
		for sec in range(secs_to_wait, 0, -1):
			self.press_controlc(secs_to_wait, sec)
			time.sleep(1)
		if log != None: log.info("")
		else: print ""

	# return value of the confdir
	def get_confdir(self):
		return self.opts.confdir

	# checks if given arguments are OK and not mutually exclusive, etc.
	def check_options(self, cep2, log=None):

		# with Cobalt there is no parset file, so for --pulsar option we the also allow to use "meta" special word
		# to read the metadata from .h5, similar to BG/P parsets. In the code though we are checking for the "parset"
		# word, thus we replace 'meta' (if given) to 'parset'
		self.opts.psr = self.opts.psr.replace("meta", "parset")

		# changing the logging level for debug if option is given
		if log != None and self.opts.is_debug:
			log.set_loglevel(logging.DEBUG)

		# checking if --local is used without specified beam
		if self.opts.is_local and self.opts.beam_str == "":
			msg="You have to use --beams option with --local!"
			if log != None: log.error(msg)
			else: print msg
			raise Exception

		# checking if both --summary and --nosummary are used
		if self.opts.is_summary and self.opts.is_nosummary:
			msg="Mutually exclusive options: --summary and --nosummary. Choose one!"
			if log != None: log.error(msg)
			else: print msg
			raise Exception

		# check if all required options are given
		if self.opts.obsid == "":
			if not self.opts.is_auto:
				msg="ObsID is not given. What do you want to process?"
				if log != None: log.error(msg)
				else: print msg
				raise Exception

		# check if confdir is specified or not
		if self.opts.confdir == "":
			self.opts.confdir = "%s" % (cep2.get_logdir())

		# check if parfile directory is specified. If it is, then we override the system directory with parfiles
		if self.opts.parfiledir != "":
			cep2.parfile_dir = self.opts.parfiledir

		# when do only summaries (or plots-only or --nodecode) then ignore --del option if given, otherwise 
		# everything will be deleted and if raw data are already erased then we are screwed
		if (self.opts.is_summary or self.opts.is_plots_only or self.opts.is_nodecode or self.opts.is_feedback) and self.opts.is_delete:
			self.opts.is_delete = False
			msg="***\n*** Warning: You give --del with one of other options (--summary or --plots-only or --nodecode).\n\
*** Deleting of previous results will be ignored and new results will be overwritten.\n***"
			if log != None: log.warning(msg)
			else: print msg

		# set to ignore checking for rawdata if one of the flags below is true:
		if self.opts.is_summary or self.opts.is_feedback or self.opts.is_plots_only or self.opts.is_nodecode:
			self.opts.is_skip_check_rawdata = True

		# checking that if --beams used then beams are specified correctly
		# we have this complicated "if" because we used --beams to pass summary locus node when --summary and --local
		if self.opts.beam_str != "" and (not self.opts.is_summary or (self.opts.is_summary and not self.opts.is_local)):
			# checking first if our list of beams is actually the list of excluded beams
			if self.opts.beam_str[0] == '^':
				is_excluded = True
				self.opts.beam_str = self.opts.beam_str[1:]
				if self.opts.beam_str == "":
					msg="Option --beams should have at least one excluded beam!"
					if log != None: log.error(msg)
					else: print msg
					raise Exception
			else: is_excluded = False
			if re.search(r'[^\,\:\d\/]+', self.opts.beam_str) is not None:
				msg="Option --beams can only has digits, colons, commas and in some cases / for parts!"
				if log != None: log.error(msg)
				else: print msg
				raise Exception
			elif re.search(r'[\:]+', self.opts.beam_str) is None:
				msg="Option --beams should have at least one colon!"
				if log != None: log.error(msg)
				else: print msg
				raise Exception
			else:   # forming array of beams
				beams=self.opts.beam_str.split(",")
				# also, we have to remove --beams option with argument from self.options
				# deleting here all instances of --beams (if several) and its arguments
				for jj in reversed([ii for ii in range(len(self.options)) if self.options[ii] == '--beams']):
					del(self.options[jj:jj+2])
				# checking if neither SAP or TAB are empty
				for bb in beams:
					(sap, tab) = bb.split(":")
					if sap == "" or tab == "":
						msg="Option --beams has at least one empty SAP or TAB value!"
						if log != None: log.error(msg)
						else: print msg
						raise Exception
				# defining proper lists of beams
				if is_excluded: self.user_excluded_beams = beams
				else: self.user_beams = beams

		# warning user that some of the results can still be overwritten, if --del is not used
		if not self.opts.is_delete:
			if not self.opts.is_quiet:
				msg="***\n*** Warning: Some of the previous results still can be overwritten.\n\
*** You may want to use --del to have clean run, or specify new output directory.\n***"
				if log != None: log.warning(msg)
				else: print msg
				self.waiting_for_user(10, log)

		# checking if rawdir is specified
		if self.opts.rawdir != "/data":
			cep2.rawdir = self.opts.rawdir

		# checking if stokes is specified and has a valid value
		if self.opts.stokes != -1:
			if self.opts.stokes != 0 and self.opts.stokes != 1 and self.opts.stokes != 2 and self.opts.stokes != 3:
				self.opts.stokes = -1

		# NONE is ignored as pulsar name
		if self.opts.psr != "":
       	        	self.psrs=[psr for psr in self.opts.psr.split(",") if psr != "NONE"]
			self.psrs=list(np.unique(self.psrs))

		# checking --nofold and pulsar list
		if not self.opts.is_nofold and len(self.psrs) == 0:
			if not self.opts.is_quiet:
				msg="***\n*** Warning: Pulsar is not specified and PULP will use source names\n\
*** from the parset file first if given, and then will look for the best\n\
*** pulsar in the SAP's FOV for each SAP separately! You also can use\n\
*** predefined words: \"parset\", \"sapfind\", \"sapfind3\", \"tabfind\", or \"tabfind+\".\n\
*** See help for more details.\n***"
				if log != None: log.warning(msg)
				else: print msg
				self.waiting_for_user(10, log)

		# checking if number of threads is good
		if self.opts.nthreads <= 0:
			msg="Number of threads should be positive! Given is %d. Exiting..." % (self.opts.nthreads)
			if log != None: log.error(msg)
			else: print msg
			raise Exception
		if self.opts.nthreads > 6:
			msg="\n*** Warning *** Number of threads %d is too high! Safe is 2-6. Best is 2-3\n" % (self.opts.nthreads)
			if log != None: log.warning(msg)
			else: print msg

		# do some basic checking that user numbers are ok
		if self.opts.first_freq_splitCS < 0: self.opts.first_freq_splitCS = 0
		if self.opts.first_freq_splitIS < 0: self.opts.first_freq_splitIS = 0
		if self.opts.nsplitsCS < -1 or self.opts.nsplitsCS == 0: self.opts.nsplitsCS = -1
		if self.opts.nsplitsIS < -1 or self.opts.nsplitsIS == 0: self.opts.nsplitsIS = -1

		if not self.opts.is_nofold:
			# creating temporary ATNF pulsar catalog listing to collect info about all psrs
			catname="%s/.tmpcat-%s" % (cep2.get_local_logdir(self), ''.join(random.choice(string.lowercase) for iii in range(7)))
			cmd="psrcat -db_file %s -nohead -o short -c 'name raj decj rajd decjd p0 dm s400' > %s" % (cep2.psrcatdb, catname)
			os.system(cmd)
			# reading B1950 and J2000 pulsar names from the catalog and checking if our pulsars are listed there
			# also reading coordinates and flux
			self.catpsrs, self.ras, self.decs, self.P0s, self.DMs, self.s400 = \
				np.loadtxt(catname, comments='#', usecols=(1,4,5,6,7,8), dtype=str, unpack=True)
			self.ras = np.array([float(i) for i in self.ras])
			self.decs = np.array([float(i) for i in self.decs])
			# removing this temporary catalog
			cmd="rm -f %s" % (catname)
			os.system(cmd)
			# filtering out those pulsars that do not have known values of either P0 or DM
			crit=(self.P0s != "*")&(self.DMs != "*")
			self.catpsrs = self.catpsrs[crit]
			self.P0s = self.P0s[crit]
			self.DMs = self.DMs[crit]
			self.s400 = self.s400[crit]
			self.ras = self.ras[crit]
			self.decs = self.decs[crit]

			# checking if given psr(s) names are valid, and these pulsars are in the catalog
			if len(self.psrs) != 0 and self.psrs[0] != "parset" and self.psrs[0] != "sapfind" and \
				self.psrs[0] != "sapfind3" and self.psrs[0] != "tabfind" and self.psrs[0] != "tabfind+":
				if self.opts.parfile != "" and len(self.psrs) > 1:
					msg="Parfile '%s' is given, but more than 1 pulsar are given to fold. Exiting..." % (self.opts.parfile)
					if log != None: log.error(msg)
					else: print msg
					raise Exception
				# copying parfile (if given) to temporary ~/.pulp/<obsid> dir
				# and re-defining parfile as located either in .pulp/<obsid> dir or dir with all par-files
				if self.opts.parfile != "":
					# checking first if parfile is in the current directory. If not then checking
					# if this parfile exists in the directory with all parfiles
					if os.path.exists(self.opts.parfile):
						# now checking if the path in parfile has '.pulp'. If not, then copy the file to ~/.pulp/<obsid>
						if re.search(r'\.pulp', self.opts.parfile) is None:
							msg="Copying parfile '%s' to %s..." % (self.opts.parfile, cep2.get_logdir())
							if log != None: log.info(msg)
							else: print msg
							cmd="cp -f %s %s" % (self.opts.parfile, cep2.get_logdir())
							os.system(cmd)
							# we should also change the parfile name in self.options in order to pass correct file to processing nodes
							newpar="%s/%s" % (cep2.get_logdir(), self.opts.parfile.split("/")[-1])
							self.options = [opt != self.opts.parfile and opt or newpar for opt in self.options]
							self.opts.parfile=newpar
					else:
						msg="Checking if given parfile '%s' exists in %s directory..." % (self.opts.parfile.split("/")[-1], cep2.parfile_dir)
						if log != None: log.info(msg)
						else: print msg
						newpar="%s/%s" % (cep2.parfile_dir, self.opts.parfile.split("/")[-1])
						self.options = [re.search(self.opts.parfile, opt) is None and opt or re.sub(self.opts.parfile, newpar, opt) for opt in self.options]
						self.opts.parfile=newpar
						if not os.path.exists(self.opts.parfile):
							msg="Can't find parfile '%s'. Exiting..." % (self.opts.parfile.split("/")[-1])
							if log != None: log.error(msg)
							else: print msg
							raise Exception

				if len(self.psrs) > 3:
					msg="%d pulsars are given, but only first 3 will be used for folding" % (len(self.psrs))
					if log != None: log.warning(msg)
					else: print msg
					self.psrs=self.psrs[:3]
				# checking if given pulsar names are good
				for psr in self.psrs:
					if not check_pulsars(psr, self, cep2, log):
						msg="No parfile found for pulsar %s. Exiting..." % (psr)
						if log != None: log.error(msg)
						else: print msg
						raise Exception
			else:
				msg="No pulsar names are given. PULP will find the proper pulsar(s) to fold..."
				if log != None: log.info(msg)
				else: print msg

        # checking if raw data for specified beams are located on one of the down nodes
        def is_rawdata_available(self, cep2, obs, log=None):
                # first forming the actual list of beams to process taking also into account
                # cmdline flags, like --noIS, --noCS, --noCV, --noFE
                if len(self.user_beams) > 0: self.beams = self.user_beams
                else:
                        self.beams = []
                        for sap in obs.saps:
                                for tab in sap.tabs:
                                        beam="%d:%d" % (sap.sapid, tab.tabid)
                                        # checking if this beam is already excluded
                                        if beam in self.user_excluded_beams: continue
                                        # ignoring IS beams
                                        if self.opts.is_noIS and not tab.is_coherent: continue
                                        # ignoring FE beams (both CS and CV)
                                        if self.opts.is_noFE and tab.is_coherent and tab.specificationType == "flyseye": continue
                                        # ignoring CS beams
                                        if self.opts.is_noCS and obs.CS and tab.is_coherent and tab.specificationType != "flyseye": continue
                                        # ignoring CV beams
                                        if self.opts.is_noCV and obs.CV and tab.is_coherent and tab.specificationType != "flyseye": continue
                                        self.beams.append(beam)

                # now we are checking if raw data are available for beams we want to process
		if  not self.opts.is_skip_check_rawdata:
	               	msg="Checking if all raw data/nodes are available for user-specified beams..."
        	        if log != None: log.info(msg)
       	        	else: print msg

	               	# if some TABs have raw data in several locations
        	       	# we also need to check if hoover nodes are up
               		avail_hoover_nodes=list(set(cep2.hoover_nodes).intersection(set(cep2.alive_nodes)))

	                excluded_beams_id=[]
       		        for ii in range(len(self.beams)):
               		        sapid=int(self.beams[ii].split(":")[0])
                       		tabid=int(self.beams[ii].split(":")[1])
	                        for ss in obs.saps:
        	                        if ss.sapid == sapid:
                	                        for tt in ss.tabs:
                        	                        if tt.tabid == tabid:
                                	                        tab = tt
                                        	                break

				# checking if for this beam we already know that either node is down or files are missing
				if not tab.is_data_available: 
					excluded_beams_id.append(ii)
					continue
				if len(tab.location) > 0:
					# if here, it means node is available for this beam
                	                if len(tab.location) > 1 and not self.opts.is_nohoover and len(avail_hoover_nodes) != len(cep2.hoover_nodes):
                        	       	        loc=""
                                	       	if tab.is_coherent and "locus101" not in avail_hoover_nodes: loc="locus101"
#						if not tab.is_coherent and "locus102" not in avail_hoover_nodes: loc="locus102"
						if not tab.is_coherent and "locus101" not in avail_hoover_nodes: loc="locus101"
						if loc != "":
        	                                        excluded_beams_id.append(ii)
                	                       	        msg="Hoover node %s is not available for the beam %d:%d [#locations = %d] - excluded" % (loc, sapid, tabid, len(tab.location))
                        	                       	if log != None: log.warning(msg)
							else: print msg
	               	        else: # no data available
        	               	        excluded_beams_id.append(ii)
               		                msg="No data available for the beam %d:%d - excluded" % (sapid, tabid)
                       		        if log != None: log.warning(msg)
                               		else: print msg

	                # now giving summary of excluded beams and deleted them from the list
       		        if len(excluded_beams_id) > 0:
               		        msg="Excluded beams [%d]: %s" % (len(excluded_beams_id), ", ".join([self.beams[id] for id in excluded_beams_id]))
                       		if log != None: log.info(msg)
				else: print msg
       		                # deleting these excluded beams from the cmdline.beams list
               		        for id in reversed(excluded_beams_id):
                       		        del(self.beams[id])
	                else:
       		                if len(self.beams) > 0:
               		                msg="All data/nodes are available"
                       		        if log != None: log.info(msg)
                               		else: print msg

	# updating cmdline default parameters based on obtained info about Observation
	# such as, number of frequency splits
	# and about FWHMs of CS and IS beams (depends on what observation is, HBA or LBA)
	def update_default_values(self, obs, cep2, log=None):
		# updating number of splits...
		if self.opts.first_freq_splitCS >= obs.nsplitsCS: self.opts.first_freq_splitCS = 0
		if self.opts.nsplitsCS == -1: self.opts.nsplitsCS = obs.nsplitsCS
		if self.opts.first_freq_splitCS + self.opts.nsplitsCS > obs.nsplitsCS:
			self.opts.nsplitsCS -= (self.opts.first_freq_splitCS + self.opts.nsplitsCS - obs.nsplitsCS)
		# for IS beams now
		if self.opts.first_freq_splitIS >= obs.nsplitsIS: self.opts.first_freq_splitIS = 0
		if self.opts.nsplitsIS == -1: self.opts.nsplitsIS = obs.nsplitsIS
		if self.opts.first_freq_splitIS + self.opts.nsplitsIS > obs.nsplitsIS:
			self.opts.nsplitsIS -= (self.opts.first_freq_splitIS + self.opts.nsplitsIS - obs.nsplitsIS)
		# updating the real values of FWHM to use
		if self.opts.fwhm_CS < 0.0:
			if obs.antenna == "HBA": self.opts.fwhm_CS = cep2.fwhm_hba
			if obs.antenna == "LBA": self.opts.fwhm_CS = cep2.fwhm_lba
		if self.opts.fwhm_IS < 0.0:
			if obs.antenna == "HBA": self.opts.fwhm_IS = cep2.fov_hba
			if obs.antenna == "LBA": self.opts.fwhm_IS = cep2.fov_lba
		# updating subintegration time
		if self.opts.tsubint == 0 or self.opts.tsubint <= -1:
			if obs.CV: self.opts.tsubint = 5
			else: self.opts.tsubint = 60  # CS/IS
		# updating stokes if it was set up but not relevant
		if self.opts.stokes != -1:
			if (obs.CS and obs.stokesCS == "IQUV") or (obs.IS and obs.stokesIS == "IQUV"): pass
			else: self.opts.stokes == -1
			

	# print summary of all set input parameters
	def print_summary(self, cep2, obs, log=None):
		if log != None:
			log.info("")
			if self.opts.is_auto:
				log.info("Pulsar Pipeline (automatic run)")
			else:
				log.info("Pulsar Pipeline, %s" % (self.version))
			log.info("Prg: %s" % (self.prg))
			log.info(" -> SLURM = %s, Global FS = %s, Docker = %s" % \
				(self.opts.is_slurm and "yes" or "no (SSH)", self.opts.is_globalfs and "yes" or "no", self.opts.is_docker and "yes" or "no"))
			if len(self.user_beams) == 0 and len(self.user_excluded_beams) == 0:
				log.info("Cmdline: %s %s" % (self.prg.split("/")[-1], " ".join(self.options)))
			elif len(self.user_beams) != 0:
				log.info("Cmdline: %s %s" % (self.prg.split("/")[-1], " ".join(self.options + ['--beams'] + [",".join(self.user_beams)])))
			else:
				log.info("Cmdline: %s %s" % (self.prg.split("/")[-1], " ".join(self.options + ['--beams'] + ["^"+",".join(self.user_excluded_beams)])))
			log.info("")
			log.info("ObsID = %s" % (obs.id))
			if self.opts.is_nofold: pulsar_status = "No folding"
			else:
				if len(self.psrs) == 0: 
					if self.opts.is_cobalt:
						pulsar_status = "default:meta -> sapfind"
					else:
						pulsar_status = "default:parset -> sapfind"
				else: 
					if self.psrs[0] == "parset" or self.psrs[0] == "sapfind" or self.psrs[0] == "sapfind3" \
						or self.psrs[0] == "tabfind" or self.psrs[0] == "tabfind+":
						pulsar_status = self.psrs[0]
						if self.opts.is_cobalt: 
							pulsar_status = pulsar_status.replace("parset", "meta")
					else: pulsar_status = ", ".join(self.psrs)
			log.info("PSR(s) = %s" % (pulsar_status))
			if not self.opts.is_nofold:
				if len(self.psrs) == 0:
					for sap in obs.saps:
						if sap.source != "" and check_pulsars(sap.source, self, cep2, None): log.info("SAP=%d   PSR: %s" % (sap.sapid, sap.source))
						else: log.info("SAP=%d   PSR(s): %s" % (sap.sapid, len(sap.psrs)>0 and sap.psrs[0] or ""))
				if len(self.psrs) != 0 and self.psrs[0] == "sapfind":
					for sap in obs.saps: log.info("SAP=%d   PSR(s): %s" % (sap.sapid, len(sap.psrs)>0 and sap.psrs[0] or ""))
				if len(self.psrs) != 0 and self.psrs[0] == "sapfind3":
					for sap in obs.saps: log.info("SAP=%d   PSR(s): %s" % (sap.sapid, ", ".join(sap.psrs)))
				if len(self.psrs) != 0 and self.psrs[0] == "parset":
					for sap in obs.saps: log.info("SAP=%d   PSR: %s" % (sap.sapid, sap.source))
				if len(self.psrs) != 0 and self.psrs[0] == "tabfind+":
					for sap in obs.saps:
						if sap.source != "" and check_pulsars(sap.source, self, cep2, None): 
							# check if there are pulsars in SAP and if the first one is not the same as in the parset
							if len(sap.psrs) > 0 and sap.source != sap.psrs[0]:
								log.info("SAP=%d   PSR(s): %s, %s" % (sap.sapid, sap.source, len(sap.psrs)>0 and sap.psrs[0] or ""))
							else: log.info("SAP=%d   PSR(s): %s" % (sap.sapid, sap.source))
						else: log.info("SAP=%d   PSR(s): %s" % (sap.sapid, len(sap.psrs)>0 and sap.psrs[0] or ""))
			if self.opts.parfile != "":
				log.info("User-specified Parfile = %s" % (self.opts.parfile))
			if self.opts.rawdir[:5] != "/data":
				log.info("User-specified Raw data directory = %s" % (self.opts.rawdir))
			if self.opts.parset != "":
				log.info("User-specified Parset file = %s" % (self.opts.parset))
			if not self.opts.is_auto:
				log.info("Output Dir = %s*/%s_*/%s" % (cep2.processed_dir_root, cep2.processed_dir_prefix, self.opts.outdir != "" and self.opts.outdir or obs.id + "_red*"))
			else:
				log.info("Output Dir = %s" % (cep2.processed_dir))
			log.info("Delete previous results = %s" % (self.opts.is_delete and "yes" or "no"))
			log.info("Log files mode = %s" % (self.opts.is_log_append and "append" or "overwrite"))
			if self.opts.is_nosummary: log.info("NO summaries")
			elif self.opts.is_summary: log.info("Summaries ONLY")
			elif self.opts.is_feedback: log.info("Feedbacks ONLY")
			else:
				skipped=""
				if obs.CS and self.opts.is_noCS: skipped += " CS"
				if obs.IS and self.opts.is_noIS: skipped += " IS"
				if obs.CV and self.opts.is_noCV: skipped += " CV"
				if obs.FE and self.opts.is_noFE: skipped += " FE"
				if skipped != "":
					log.info("Skipped processing of: %s" % (skipped))
				if len(self.user_beams) != 0:
					log.info("User-specified BEAMS to process: %s" % (", ".join(self.user_beams)))
				if len(self.user_excluded_beams) != 0:
					log.info("User-specified BEAMS to be excluded: %s" % (", ".join(self.user_excluded_beams)))
				if (obs.CS and obs.stokesCS == "IQUV") or (obs.IS and obs.stokesIS == "IQUV"):
					log.info("STOKES = %s" % (self.opts.stokes == -1 and "IQUV" or "IQUV"[self.opts.stokes]))
				if self.opts.is_plots_only: log.info("Diagnostic plots ONLY")
				else:
					log.info("USING HOOVER NODES = %s" % (self.opts.is_nohoover and "no" or "yes"))
					log.info("RAW DATA 32 -> 8 bits = %s" % (self.opts.is_raw_to_8bit and "yes (-s %g)" % (self.opts.digitize_sigma) or "no"))
					if obs.CV: 
						log.info("DSPSR with LOFAR DAL = %s" % (self.opts.is_nodal and "no" or "yes"))
					else:
						log.info("DSPSR with LOFAR DAL = %s" % (self.opts.is_with_dal and "yes" or "no"))
					if self.opts.first_freq_splitCS != 0:
						log.info("FIRST FREQUENCY SPLIT (CS) = %d" % (self.opts.first_freq_splitCS))
					if self.opts.first_freq_splitIS != 0:
						log.info("FIRST FREQUENCY SPLIT (IS) = %d" % (self.opts.first_freq_splitIS))
					if obs.CV or obs.CS:
						if self.opts.nsplitsCS != -1:
							log.info("NUMBER OF SPLITS (CS) = %d" % (self.opts.nsplitsCS))
					if obs.IS:
						if self.opts.nsplitsIS != -1:
							log.info("NUMBER OF SPLITS (IS) = %d" % (self.opts.nsplitsIS))
					log.info("Data decoding = %s" % (self.opts.is_nodecode and "no" or "yes (-A %d -sigma %d)" % (self.opts.decode_nblocks, self.opts.decode_sigma)))
					if not obs.CV and not self.opts.is_nodecode and not self.opts.is_with_dal and self.opts.bf2fits_extra_opts != "":
						log.info("2bf2fits user extra options: %s" % (self.opts.bf2fits_extra_opts))
					log.info("RFI Zapping = %s" % (self.opts.is_norfi and "no" or "yes"))
					if self.opts.rfifind_extra_opts != "" and not self.opts.is_norfi:
						log.info("Rfifind user extra options: %s" % (self.opts.rfifind_extra_opts))
					log.info("Subdyn.py = %s" % ((self.opts.is_skip_subdyn == False and self.opts.is_norfi == False) and \
						"yes (averaging time = %g s)" % (self.opts.subdyn_time_average) or "no"))
					log.info("Prepfold = %s" % (self.opts.is_skip_prepfold and "no" or "yes"))
					if self.opts.prepfold_extra_opts != "" and not self.opts.is_skip_prepfold:
						log.info("Prepfold user extra options: %s" % (self.opts.prepfold_extra_opts))
					log.info("Single-pulse analysis = %s" % (self.opts.is_single_pulse and "yes" or "no"))
					log.info("RRATs analysis = %s" % (self.opts.is_rrats and "yes" or "no"))
					if self.opts.prepsubband_extra_opts != "" and self.opts.is_rrats:
						log.info("Prepsubband user extra options: %s" % (self.opts.prepsubband_extra_opts))
					if self.opts.prepdata_extra_opts != "" and (self.opts.is_single_pulse or self.opts.is_rrats):
						log.info("Prepdata user extra options: %s" % (self.opts.prepdata_extra_opts))
					log.info("DSPSR = %s" % (self.opts.is_skip_dspsr and "no" or \
						(self.opts.nthreads == 2 and "yes, #threads = %d (default)" % (self.opts.nthreads) or "yes, #threads = %d" % (self.opts.nthreads))))
					if self.opts.dspsr_extra_opts != "" and not self.opts.is_skip_dspsr:
						log.info("DSPSR user extra options: %s" % (self.opts.dspsr_extra_opts))
					if self.opts.digifil_extra_opts != "":
						log.info("Digifil user extra options: %s" % (self.opts.digifil_extra_opts))
					log.info("pdmp = %s" % ((self.opts.is_nopdmp or self.opts.is_nofold) and "no" or "yes"))
					if obs.CV: log.info("rmfit = %s" % (self.opts.is_skip_rmfit and "no" or "yes"))
			log.info("")
