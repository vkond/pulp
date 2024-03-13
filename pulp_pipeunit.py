###################################################################
#
# Class PipeUnit - aka Processing unit per Beam 
# and other derivative classes
#

import os, sys, glob, time, re, os.path
import math
import numpy as np
import cPickle
import subprocess, shlex
from subprocess import PIPE, STDOUT, Popen
import psr_utils as pu
from pulp_parset import Observation, radial_distance, find_pulsars
from pulp_usercmd import CMDLine, check_pulsars
from pulp_sysinfo import CEP2Info
from pulp_logging import PulpLogger
from pulp_feedback import FeedbackUnit


# return list of dirs in the current base directory
def getDirs(base):
	return [x for x in glob.iglob(os.path.join(base, '*')) if os.path.isdir(x)]

# recursive glob
def rglob(base, pattern, maxdepth=0):
	flist=[]
	flist.extend(glob.glob(os.path.join(base, pattern)))
	dirs = getDirs(base)
	if maxdepth <= 0: return flist
	if len(dirs):
		for d in dirs:
			flist.extend(rglob(os.path.join(base, d), pattern, maxdepth-1))
	return flist

# common postprocessing routines after DSPSR calls for different Units and both DAL and non_DAL parts
# root - class instance to execute
# ref - class instance to get values of chans, etc.
def dspsr_postproc(root, ref, cmdline, obs, psr, total_chan, nsubs_eff, curdir, output_prefix):
	# removing corrupted freq channels
	if ref.nrChanPerSub > 1:
		root.log.info("Zapping every %d channel..." % (ref.nrChanPerSub))
		cmd="paz -z \"%s\" -m %s_%s.ar" % \
			(" ".join([str(jj) for jj in range(0, total_chan, ref.nrChanPerSub)]), psr, output_prefix)
		root.execute(cmd, workdir=curdir)

	# rfi zapping
	if not cmdline.opts.is_norfi:
		# first we check how large the dataset is (product of number of channels and subints)
		# this is necessary as clean.py uses a lot of memory. If som then we use old-fashioned paz -r
		if (obs.duration/cmdline.opts.tsubint) * total_chan >= 512000:
			root.log.info("Zapping channels using median smoothed difference algorithm...")
			cmd="paz -r -e paz.ar %s_%s.ar" % (psr, output_prefix)
			root.execute(cmd, workdir=curdir)
		else:
			root.log.info("Running Patrick Lazarus's COAST_GUARD RFI cleaner using surgical approach...")
			cmd="clean.py -F surgical -o %s_%s.paz.ar %s_%s.ar" % (psr, output_prefix, psr, output_prefix)
			try:
				root.execute(cmd, workdir=curdir)
			except:
				root.log.info("COAST_GUARD RFI cleaner failed! Will try to use paz -r...")
				root.log.info("Zapping channels using median smoothed difference algorithm...")
				cmd="paz -r -e paz.ar %s_%s.ar" % (psr, output_prefix)
				root.execute(cmd, workdir=curdir)

	# dedispersing
	# checking if there was already an option -K. That means we do not need to run dedispersion as all sub-integrations
	# have been aligned already
	if re.match("^\-K$", cmdline.opts.dspsr_extra_opts) or re.match("^.*\s+\-K$", cmdline.opts.dspsr_extra_opts) or re.match("^.*\s+\-K\s+.*$", cmdline.opts.dspsr_extra_opts) or re.match("^\-K\s+.*$", cmdline.opts.dspsr_extra_opts):
		cmd="mv %s_%s.ar %s_%s.dd" % (psr, output_prefix, psr, output_prefix)
		root.execute(cmd, workdir=curdir)
	else:
		root.log.info("Dedispersing...")
		if not cmdline.opts.is_norfi or os.path.exists("%s/%s_%s.paz.ar" % (curdir, psr, output_prefix)):
			cmd="pam -D -m %s_%s.paz.ar" % (psr, output_prefix)
			root.execute(cmd, workdir=curdir)
		cmd="pam -D -e dd %s_%s.ar" % (psr, output_prefix)
		root.execute(cmd, workdir=curdir)

	# scrunching in frequency
	root.log.info("Scrunching in frequency to have %d channels in the output ar-file..." % (nsubs_eff))
	if ref.nrChanPerSub > 1:
		# first, running fscrunch on zapped archive
		if not cmdline.opts.is_norfi or os.path.exists("%s/%s_%s.paz.ar" % (curdir, psr, output_prefix)):
			cmd="pam --setnchn %d -e fscr.AR %s_%s.paz.ar" % (nsubs_eff, psr, output_prefix)
			root.execute(cmd, workdir=curdir)
			# remove non-scrunched zapped archive (we will always have unzapped non-scrunched version)
			cmd="rm -f %s_%s.paz.ar" % (psr, output_prefix)
			root.execute(cmd, workdir=curdir)
		# running fscrunching on non-zapped archive
		cmd="pam --setnchn %d -e fscr.AR %s_%s.dd" % (nsubs_eff, psr, output_prefix)
		root.execute(cmd, workdir=curdir)
		# remove non-scrunched dedispersed archive (we will always have unzapped non-dedispersed non-scrunched version)
		cmd="rm -f %s_%s.dd" % (psr, output_prefix)
		root.execute(cmd, workdir=curdir)
	else: # if number of chans == number of subs, we will just rename .paz.ar to .paz.fscr.AR and .dd to .fscr.AR
		if not cmdline.opts.is_norfi or os.path.exists("%s/%s_%s.paz.ar" % (curdir, psr, output_prefix)):
			cmd="mv -f %s_%s.paz.ar %s_%s.paz.fscr.AR" % (psr, output_prefix, psr, output_prefix)
			root.execute(cmd, workdir=curdir)
		cmd="mv -f %s_%s.dd %s_%s.fscr.AR" % (psr, output_prefix, psr, output_prefix)
		root.execute(cmd, workdir=curdir)

# running pav's to make DSPSR diagnostic plots
# root - class instance to execute
# ref - class instance to get values of chans, etc.
def make_dspsr_plots(root, ref, cmdline, obs, nsubs_eff, curdir, output_prefix, is_pdmp=False):

	# first, calculating the proper max divisor for the number of subbands
	root.log.info("Getting proper value of nchans in pav -f between %d and %d..." % (1, min(nsubs_eff, 63)))
	# calculating the greatest common denominator of self.tab.nrSubbands starting from 63 down
	pav_nchans = ref.hcd(1, min(nsubs_eff, 63), nsubs_eff)
	fscrunch_factor = nsubs_eff / pav_nchans
	if is_pdmp: suffix=".pdmp.AR"
	else: suffix=".AR"
	root.log.info("Creating diagnostic plots...")
	for psr in ref.psrs:  # pulsar list is empty if --nofold is used
		if not cmdline.opts.is_norfi or os.path.exists("%s/%s_%s.paz.fscr%s" % (curdir, psr, output_prefix, suffix)):
			output_stem=".paz.fscr%s" % (suffix)
		else: output_stem=".fscr%s" % (suffix)

		# creating DSPSR diagnostic plots
		if obs.CV: plot_type = "SFTd"
		else: plot_type = "DFTpd"

		cmd="pav -%s -C -g %s_%s_%s.ps/cps %s_%s%s" % (plot_type, psr, output_prefix, plot_type, psr, output_prefix, output_stem)
		root.execute(cmd, workdir=curdir)
		cmd="pav -GTpdf%d -C -g %s_%s_GTpdf%d.ps/cps %s_%s%s" % (fscrunch_factor, psr, output_prefix, fscrunch_factor, psr, output_prefix, output_stem)
		root.execute(cmd, workdir=curdir)
		cmd="pav -YFpd -C -g %s_%s_YFpd.ps/cps %s_%s%s" % (psr, output_prefix, psr, output_prefix, output_stem)
		root.execute(cmd, workdir=curdir)
		cmd="pav -j -g %s_%s_j.ps/cps %s_%s%s" % (psr, output_prefix, psr, output_prefix, output_stem)
		root.execute(cmd, workdir=curdir)
		if not cmdline.opts.is_skip_rmfit and obs.CV:
			try:
				# running rmfit for negative and positive RMs
				cmd="rmfit -m -100,0,100 -D -K %s_%s.negRM.ps/cps %s_%s%s" % (psr, output_prefix, psr, output_prefix, output_stem)
				root.execute(cmd, workdir=curdir)
				cmd="rmfit -m 0,100,100 -D -K %s_%s.posRM.ps/cps %s_%s%s" % (psr, output_prefix, psr, output_prefix, output_stem)
				root.execute(cmd, workdir=curdir)
				cmd="convert \( %s_%s_GTpdf%d.ps %s_%s_j.ps %s_%s.posRM.ps +append \) \( %s_%s_%s.ps %s_%s_YFpd.ps %s_%s.negRM.ps +append \) -append -rotate 90 -background white -flatten %s_%s_%s.png" % \
					(psr, output_prefix, fscrunch_factor, psr, output_prefix, psr, output_prefix, psr, output_prefix, plot_type, \
					psr, output_prefix, psr, output_prefix, psr, output_prefix, is_pdmp and "diag_pdmp" or "diag")
				root.execute(cmd, workdir=curdir)
			except Exception:
				root.log.warning("***** Warning! Rmfit has failed. Diagnostic plots will be made without rmfit plots. *****")
				cmd="convert \( %s_%s_GTpdf%d.ps %s_%s_j.ps +append \) \( %s_%s_%s.ps %s_%s_YFpd.ps +append \) -append -rotate 90 -background white -flatten %s_%s_%s.png" % \
					(psr, output_prefix, fscrunch_factor, psr, output_prefix, psr, output_prefix, plot_type, \
					psr, output_prefix, psr, output_prefix, is_pdmp and "diag_pdmp" or "diag")
				root.execute(cmd, workdir=curdir)
		else:
			cmd="convert \( %s_%s_GTpdf%d.ps %s_%s_j.ps +append \) \( %s_%s_%s.ps %s_%s_YFpd.ps +append \) -append -rotate 90 -background white -flatten %s_%s_%s.png" % \
				(psr, output_prefix, fscrunch_factor, psr, output_prefix, psr, output_prefix, plot_type, \
				psr, output_prefix, psr, output_prefix, is_pdmp and "diag_pdmp" or "diag")
			root.execute(cmd, workdir=curdir)

# start pdmp calls
# root - class instance to execute
# ref - class instance to get values of chans, etc.
def start_pdmp(root, ref, cmdline, obs, nsubs_eff, curdir, output_prefix):

	# now running pdmp without waiting...
	root.log.info("Running pdmp...")
	pdmp_popens=[]  # list of pdmp Popen objects	
	for psr in ref.psrs:
		if not cmdline.opts.is_norfi or os.path.exists("%s/%s_%s.paz.fscr.AR" % (curdir, psr, output_prefix)):
			output_stem=".paz.fscr.AR"
		else: output_stem=".fscr.AR"
		# getting the number of bins in the ar-file (it can be different from self.get_best_nbins, because
		# we still provide our own number of bins in --dspsr-extra-opts
		try:
			cmd="psredit -q -Q -c nbin %s/%s_%s%s" % (curdir, psr, output_prefix, output_stem)
			binsline=os.popen(cmd).readlines()
			if np.size(binsline) > 0:
				best_nbins=int(binsline[0][:-1])
				cmd="pdmp -mc %d -mb %d -g %s_%s_pdmp.ps/cps %s_%s%s" % (nsubs_eff, min(128, best_nbins), psr, output_prefix, psr, output_prefix, output_stem)
				pdmp_popen = root.start_and_go(cmd, workdir=curdir)
				pdmp_popens.append(pdmp_popen)
		except Exception: pass
	return pdmp_popens

# finish pdmps
# root - class instance to execute
# ref - class instance to get values of chans, etc.
def finish_pdmp(root, ref, pdmp_popens, cmdline, obs, curdir, output_prefix):
	# waiting for pdmp to finish
	root.waiting_list("pdmp", pdmp_popens)
	# when pdmp is finished do extra actions with files...
	for psr in ref.psrs:
		if not cmdline.opts.is_norfi or os.path.exists("%s/%s_%s.paz.fscr.AR" % (curdir, psr, output_prefix)):
			output_stem=".paz.fscr.AR"
		else: output_stem=".fscr.AR"
		cmd="grep %s %s/pdmp.per > %s/%s_%s_pdmp.per" % (psr, curdir, curdir, psr, output_prefix)
		root.execute(cmd, is_os=True)
		cmd="grep %s %s/pdmp.posn > %s/%s_%s_pdmp.posn" % (psr, curdir, curdir, psr, output_prefix)
		root.execute(cmd, is_os=True)
		# getting new topo period
		logf = open(root.log.get_logfile(), 'r')
		newp0s = [str(float(ff.split()[5])/1000.) for ff in logf.read().splitlines() if re.search("^Best TC", ff) is not None]
		logf.close()
		if np.size(newp0s) > 1: newp0 = newp0s[-1]
		else: newp0 = newp0s[0]
		# reading new DM from the *.per file
		newdm = np.loadtxt("%s/%s_%s_pdmp.per" % (curdir, psr, output_prefix), comments='#', usecols=(3,3), dtype=float, unpack=True)[0]
		if np.size(newdm) > 1: cmd="pam -e pdmp.AR -d %f --period %s -D %s_%s%s" % (newdm[-1], newp0, psr, output_prefix, output_stem)
		else: cmd="pam -e pdmp.AR -d %f --period %s -D %s_%s%s" % (newdm, newp0, psr, output_prefix, output_stem)
		root.execute(cmd, workdir=curdir)

# Running spectar.py to calculate pulsar spectra
def calc_psr_spectra(root, ref, cmdline, obs, sapid, tabid, curdir, output_prefix):
	if not cmdline.opts.is_cobalt:
		root.log.info("Calculating pulsar(s) spectra...")
		try:
			for psr in ref.psrs:
				psrar=glob.glob("%s/%s_%s.paz.fscr.pdmp.AR" % (curdir, psr, output_prefix))
				if len(psrar) == 0:
					psrar=glob.glob("%s/%s_%s.paz.fscr.AR" % (curdir, psr, output_prefix))
					if len(psrar) == 0:
						psrar=glob.glob("%s/%s_%s.fscr.AR" % (curdir, psr, output_prefix))
				if len(psrar) == 0: continue
				cmd="spectar.py -f %s -p %s -o %s/%s_%s_%s --sap %d --tab %d" % (psrar[0], obs.parset, curdir, psr, output_prefix, "spectra.txt", sapid, tabid)
				root.execute(cmd, workdir=curdir)
		except: pass
	else:
		root.log.info("Calculating pulsar(s) spectra...skipped. No observation parset anymore.")

# patch to fill in source coordinates to the header of output ar-file
def fix_header_coords(root, ref, psr, curdir, output_prefix):
	# getting coordinates string of the pulsar
	psr_ra=pu.coord_to_string(*pu.rad_to_hms(ref.tab.rarad))
	psr_dec=pu.coord_to_string(*pu.rad_to_dms(ref.tab.decrad))
	if ref.tab.decrad < 0: psr_coords=psr_ra+psr_dec
	else: psr_coords=psr_ra+"+"+psr_dec
	cmd="psredit -m -c coord=%s %s_%s.ar" % (psr_coords, psr, output_prefix)
	root.execute(cmd, workdir=curdir)

# running single-pulse analysis for CV data (prepdata + single_pulse_search.py) for all specified pulsars
# for CV data means, that we already ran digifil with further rfifind on a .fil file
def run_prepdata_CV(root, ref, cmdline, outdir, curdir, output_prefix):
	root.log.info("Running single-pulse analysis (prepdata + single_pulse_search.py for a pulsar(s) DM")
	try:
		for psr in ref.psrs:
			psr2=re.sub(r'^[BJ]', '', psr)
			if cmdline.opts.is_nofold:
				psrdm = 0.0
			else:
				if os.path.exists("%s/%s.par" % (outdir, psr2)):
					psrdm=ref.get_psr_dm("%s/%s.par" % (outdir, psr2))
				else: # this check is necessary because when combining parts the parfile is in the curdir, not in the outdir
					psrdm=ref.get_psr_dm("%s/%s.par" % (curdir, psr2))
			# running prepdata with mask (if --norfi was not set)
			if not cmdline.opts.is_norfi or os.path.exists("%s/%s_%s_rfifind.mask" % (curdir, psr, output_prefix)):
				cmd="prepdata -noclip -nobary -dm %f -mask %s_%s_rfifind.mask -o %s_%s %s %s_%s.fil" % \
					(psrdm, psr, output_prefix, psr, output_prefix, cmdline.opts.prepdata_extra_opts, psr, output_prefix)
			else: # running prepdata without mask
				cmd="prepdata -noclip -nobary -dm %f -o %s_%s %s %s_%s.fil" % \
					(psrdm, psr, output_prefix, cmdline.opts.prepdata_extra_opts, psr, output_prefix)
			root.execute(cmd, workdir=curdir)
			# running single_pulse_search.py on .dat file (created either with mask or not)
			cmd="single_pulse_search_lotaas.py --no-width-label %s_%s.dat" % (psr, output_prefix)
			root.execute(cmd, workdir=curdir)
	except: pass

# base class for the single processing (a-ka beam)
class PipeUnit:
	def __init__(self, obs, cep2, cmdline, tab, log, curstokes, part=0):
		self.code = ""  # 2-letter code, CS, IS, CV
		self.stokes = "" # Stokes I, IQUV, or XXYY
		self.stokes_index = curstokes
		self.sapid = tab.parent_sapid
		self.tabid = tab.tabid
		self.tab = tab
		self.parent = None   # parent Popen project
		self.procs = []      # list of open processes
		# process.pid - pid, process.returncode - status
		# if returncode == None, it means that process is still running
		self.rawdir = ""     # actual rawdir for this beam (where the processing is happening, e.g. either /data1 or /data2)
		self.processed_dir_root = ""  # root directory relative to cep2.processed_dir_prefix, i.e. where LOFAR_PULSAR_ARCHIVE resides
					      # because it can be different from cep2.rawdir  
		self.outdir = ""     # root directory with processed data
		self.curdir = ""     # current processing directory
		self.beams_root_dir = ""  # 'stokes' for CS, 'incoherentstokes' for IS
		self.summary_node = "" # summary node, for CS - locus092, for IS - locus094
		self.summary_node_dir_suffix = ""  # for CS this is "_CSplots" and for IS - "_redIS"
		self.archive_suffix = "" # "_plotsCS.tar.gz" for CS and "_plotsIS.tar.gz" for IS
		self.outdir_suffix = "" #  for CS this is "_red" (if user did not specify) and "_redIS" for IS
		self.raw2fits_extra_options = ""  # extra options for 2bf2fits:  -CS -H for CS, -CS -H -IS for IS
		self.raw_8bit_dir = "raw-8bit" # directory with raw 8-bit data (when requested)
		self.nrChanPerSub = 0
		self.sampling = 0
		self.downsample_factor = 0
		self.log = None
		self.logfile = ""    # basename of the logfile
		self.start_time = 0  # start time of the processing (in s)
		self.end_time = 0    # end time (in s)
		self.total_time = 0  # total time in s 
		# extensions of the files to copy to summary node
		self.extensions=["*.log", "*.txt", "*.pdf", "*.ps", "*.bestprof", "*.rfirep", "*png", "*_DM*.inf", "*.singlepulse"]
		self.procdir = "BEAM%d" % (self.tabid)
		# extensions for the full archive, e.g. LTA
		self.full_archive_exts=["*.log", "*.txt", "*.pdf", "*.ps", "*.pfd", "*.bestprof", "*.polycos", "*.inf", "*.rfirep", "*png", "*.ar", "*.AR", "*pdmp*", "*_rfifind*", "*.dat", "*.singlepulse", "*.h5", "*.fil", "*.rv", "*.out", "*.raw", "*.fits"]
		# feedback unit
		self.fbunit = None
		# bandwidth part index
		self.part = part

		# pulsars to fold for this unit
		self.psrs = []
		if not cmdline.opts.is_nofold:
			self.psrs = self.get_pulsars_to_fold(obs, cep2, cmdline, log)
			# to be sure that we have unique list of pulsars (especially relevant for tabfind+ option)
			self.psrs = np.unique(self.psrs)

			# dspsr folding options
			# making choice between -L %d   and "-s"
			# by default -L is used, but if -s is given in the dspsr_extra_opts, then we should get rid of -L
			self.dspsr_folding_options="-L %d" % (cmdline.opts.tsubint)
			if re.match("^\-s$", cmdline.opts.dspsr_extra_opts) or re.match("^.*\s+\-s$", cmdline.opts.dspsr_extra_opts) or re.match("^.*\s+\-s\s+.*$", cmdline.opts.dspsr_extra_opts) or re.match("^\-s\s+.*$", cmdline.opts.dspsr_extra_opts):
				self.dspsr_folding_options=""

	# function to set outdir and curdir directories
	def set_outdir(self, obs, cep2, cmdline):
		if len(self.tab.location) == 0:    # when raw data are erased but still want to run only summaries
			self.location=cep2.hoover_nodes[0] # basically there is no need in this, as only summaries will be done
		elif len(self.tab.location) > 1:
			if (self.tab.is_coherent and obs.stokesCS == "IQUV" and obs.nsplitsCS == 1) or (self.tab.is_coherent == False and obs.stokesIS == "IQUV" and obs.nsplitsIS == 1):
				self.location=[key for key in self.tab.rawfiles.keys() if any("_S%d_" % (self.stokes_index) in ff for ff in self.tab.rawfiles[key])][0]
			else:
				if not cmdline.opts.is_nohoover: # are not using hoover nodes
					self.location=cep2.hoover_nodes[0]
				# otherwise self.location is already defined in the PartUnit
		else:
			self.location=self.tab.location[0]
		# setting output directory
		if cmdline.opts.is_auto: # auto pipeline mode
			if cep2.cluster_headnode[:5] == "head0" or cep2.cluster_headnode[:3] == "cpu": # CEP4
				self.outdir = "%s%s" % (cep2.processed_dir, self.outdir_suffix)
			else: # not CEP4
				self.outdir = "%s" % (cep2.processed_dir)

			if cep2.rawdir[:5] == "/data":
				if cep2.cluster_headnode[:5] == "head0" or cep2.cluster_headnode[:3] == "cpu": # CEP4
					self.rawdir = "%s%s" % (cep2.rawdir, "/".join(cep2.processed_dir.split("/data")[-1].split("/")[:-1]))
					self.processed_dir_root = cep2.processed_dir_root
				else: # CEP2
					self.rawdir = "%s%s" % (cep2.rawdir, cep2.processed_dir.split("/data")[-1].split("/")[0])
					self.processed_dir_root = "%s%s" % (cep2.processed_dir_root, cep2.processed_dir.split("/data")[-1].split("/")[0])
			else: 
				if cep2.cluster_headnode[:5] == "head0" or cep2.cluster_headnode[:3] == "cpu": # CEP4
					self.rawdir = "/".join(cep2.processed_dir.split("/")[:-1])
					self.processed_dir_root = ""
				else: # CEP2
					self.rawdir = cep2.rawdir
					self.processed_dir_root = cep2.processed_dir_root
		else:
			if cep2.cluster_headnode[:5] == "head0" or cep2.cluster_headnode[:3] == "cpu": # CEP4
				if cep2.rawdir[:5] == "/data":
					myfiles=[ff for ff in self.tab.rawfiles[self.location] if "_P%03d_" % (self.part) in ff]
					self.rawdir = "%s%s" % (cep2.rawdir, "/".join(myfiles[0].split("/data")[-1].split("/")[:-3]))
					#self.processed_dir_root = "%s%s" % (cep2.processed_dir_root, "/".join(myfiles[0].split("/data")[-1].split("/")[:-3]))
				else: 
					self.rawdir = cep2.rawdir
				self.processed_dir_root = cep2.processed_dir_root
			else: # not CEP4
				if cep2.rawdir[:5] == "/data":
					myfiles=[ff for ff in self.tab.rawfiles[self.location] if "_P%03d_" % (self.part) in ff]
					self.rawdir = "%s%s" % (cep2.rawdir, myfiles[0].split("/data")[-1].split("/")[0])
					self.processed_dir_root = "%s%s" % (cep2.processed_dir_root, myfiles[0].split("/data")[-1].split("/")[0])
				else:
					self.rawdir = cep2.rawdir
					self.processed_dir_root = cep2.processed_dir_root

			# if user specified output dir (relative to /data/LOFAR_PULSAR_....)
			if cmdline.opts.outdir != "":
				if not cmdline.opts.is_globalfs:
					self.outdir = "%s/%s_%s/%s" % (self.processed_dir_root, cep2.processed_dir_prefix, self.location, cmdline.opts.outdir)
				else:
					self.outdir = "%s/%s/%s" % (self.processed_dir_root, cep2.processed_dir_prefix, cmdline.opts.outdir)
			else: # if output dir was not set
				if not cmdline.opts.is_globalfs:
					self.outdir = "%s/%s_%s/%s%s" % (self.processed_dir_root, cep2.processed_dir_prefix, self.location, obs.id, self.outdir_suffix)
				else:
					self.outdir = "%s/%s/%s%s" % (self.processed_dir_root, cep2.processed_dir_prefix, obs.id, self.outdir_suffix)
		self.curdir = "%s/%s/SAP%d/%s" % (self.outdir, self.beams_root_dir, self.sapid, self.procdir)

	# function to get feedback index for the feedback unit
	def get_feedback_index(self, obs, cmdline):
		fbindex=0
		for sap in sorted(obs.saps, key=lambda x: x.sapid):
			for tab in sorted(sap.tabs, key=lambda x: x.tabid):
				if self.sapid == sap.sapid and self.tabid == tab.tabid:
					if tab.is_coherent: 
						if obs.stokesCS == "IQUV": 
							fbindex += 4 * (self.part - cmdline.opts.first_freq_splitCS)
							fbindex += self.stokes_index
						else: fbindex += (self.part - cmdline.opts.first_freq_splitCS)
					else:
						if obs.stokesIS == "IQUV": 
							fbindex += 4 * (self.part - cmdline.opts.first_freq_splitIS)
							fbindex += self.stokes_index
						else: fbindex += (self.part - cmdline.opts.first_freq_splitIS)
					return fbindex
				else:
					if tab.is_coherent:
						if obs.stokesCS == "IQUV": 
							fbindex += 4 * cmdline.opts.nsplitsCS
						else: fbindex += cmdline.opts.nsplitsCS
					else:
						if obs.stokesIS == "IQUV": 
							fbindex += 4 * cmdline.opts.nsplitsIS
						else: fbindex += cmdline.opts.nsplitsIS
				return fbindex

	# function to initialize feedback unit (should be called after self.outdir is set)
	def feedback_init(self, obs, cep2, cmdline):
		fbindex=self.get_feedback_index(obs, cmdline)
		self.fbunit = FeedbackUnit(fbindex, self.location, self.outdir, self.sapid, self.tabid, self.stokes_index, self.part)
		if (self.tab.is_coherent and obs.stokesCS == "IQUV") or (self.tab.is_coherent == False and obs.stokesIS == "IQUV"):
			self.fbunit.update("%s/" % (self.outdir),
				"%s/.%s_SAP%03d_B%03d_S%d_P%03d.fb" % (cep2.get_logdir(), cep2.pipeid, self.sapid, self.tabid, self.stokes_index, self.part), self.code, obs)
		else:
			self.fbunit.update("%s/" % (self.outdir),
				"%s/.%s_SAP%03d_B%03d_P%03d.fb" % (cep2.get_logdir(), cep2.pipeid, self.sapid, self.tabid, self.part), self.code, obs)
		# writing initial feedback file to disk, so we have something in case pipeline will fail
		self.fbunit.flush(0, cep2, False)

	# function to get the list of pulsars to fold for this TAB (unit)
	def get_pulsars_to_fold(self, obs, cep2, cmdline, log):
		try:
			# get pulsar name from the parset
			# if pulsar is not given in the command line, I also have to find pulsar if parset entry is empty
			if len(cmdline.psrs) == 0 or cmdline.psrs[0] == "parset":
				for sap in obs.saps:
					if self.sapid == sap.sapid and sap.source != "" and check_pulsars(sap.source, cmdline, cep2, None):
						self.psrs.append(sap.source)

			# if --pulsar is not given and source field in parset is empty
			if len(cmdline.psrs) == 0 and len(self.psrs) == 0:
				for sap in obs.saps:
					if self.sapid == sap.sapid:
						self.psrs[:] = sap.psrs
						break
				if len(self.psrs)>0: self.psrs = self.psrs[:1]  # leave only one pulsar
	
			# if special word "tabfind" or "tabfind+" is given
			# if it is "tabfind+", then we take pulsar from the parset (if exist), then one pulsar from the SAP (if different from the parset)
			# and then also search for another pulsar in the TAB
			# In case of "tabfind" we only search for pulsars in the TAB
			if len(cmdline.psrs) != 0 and (cmdline.psrs[0] == "tabfind" or cmdline.psrs[0] == "tabfind+"):
				# in the special case of "tabfind+"...
				if cmdline.psrs[0] == "tabfind+":
					for sap in obs.saps:
						if self.sapid == sap.sapid: 
							if sap.source != "" and check_pulsars(sap.source, cmdline, cep2, None):
								self.psrs.append(sap.source)
							if len(sap.psrs) > 0: self.psrs.append(sap.psrs[0])
				log.info("Searching for best pulsar for folding in SAP=%d TAB=%d..." % (self.sapid, self.tabid))
				tabpsrs = find_pulsars(self.tab.rarad, self.tab.decrad, cmdline, cmdline.opts.fwhm_CS/2.)
				if len(tabpsrs) > 0: 
					self.psrs.append(tabpsrs[0]) # use only one pulsar from those found in a TAB
					log.info("%s" % (tabpsrs[0]))

			# using pulsars from SAP
			if len(cmdline.psrs) != 0 and (cmdline.psrs[0] == "sapfind" or cmdline.psrs[0] == "sapfind3"):
				for sap in obs.saps:
					if self.sapid == sap.sapid:
						self.psrs[:] = sap.psrs
						break
				if cmdline.psrs[0] == "sapfind" and len(self.psrs)>0: self.psrs = self.psrs[:1]  # leave only one pulsar

			# if --pulsar is used but no special word
			if len(cmdline.psrs) != 0 and cmdline.psrs[0] != "parset" and cmdline.psrs[0] != "tabfind" and \
					cmdline.psrs[0] != "sapfind" and cmdline.psrs[0] != "sapfind3" and cmdline.psrs[0] != "tabfind+":
				self.psrs[:] = cmdline.psrs # copying all items

			# checking if pulsars are in ATNF catalog, or if not par-files do exist fo them, if not - exit
			for psr in self.psrs:
				if not check_pulsars(psr, cmdline, cep2, log):
					log.info("*** No parfile found for pulsar %s for SAP=%d TAB=%d. Exiting..." % (psr, self.sapid, self.tabid))
					raise Exception

			# if pulsar list is still empty, and we did not set --nofold then set is_nofold flag to True
			if len(self.psrs) == 0:
				log.warning("*** No pulsar found to fold for SAP=%d TAB=%d. Setting --nofold flag for this beam..." % (self.sapid, self.tabid))
			
			return self.psrs

		except Exception:
			self.log.exception("Oops... get_pulsars_to_fold has crashed!")
			self.kill()  # killing all open processes
			raise

	# getting proper parfile in the processing directory
	def get_parfile(self, cmdline, cep2):
		# local function to check the existence of the parfile
		def check_parfile (parfile, psrstem):
			if os.path.exists(parfile):
				try:
					cmd="rsync -a %s %s/%s.par" % (parfile, self.outdir, psrstem)
					self.execute(cmd)
					# also converting them to UNIX format (if were prepared on Macs for example)
					# if they are in UNIX format already, nothing happens
					cmd="dos2unix %s/%s.par" % (self.outdir, psrstem)
					self.execute(cmd)
				except: pass
				return True
			return False

		for psr in self.psrs:
			psr2=re.sub(r'^[BJ]', '', psr)
			if cmdline.opts.parfile != "":
				if os.path.exists(cmdline.opts.parfile): 
					self.log.info("Copying user-specified parfile '%s' to %s/%s.par" % \
						(cmdline.opts.parfile, self.outdir, psr2))
					try:
						cmd="rsync -a %s %s/%s.par" % (cmdline.opts.parfile, self.outdir, psr2)
						self.execute(cmd)
					except: pass
					continue
				else: 
					self.log.error("Can't find user-specified parfile '%s'. Exiting..." % (cmdline.opts.parfile))
					self.kill()
					raise Exception

			# checking repository dir withiout ^[BJ] in the name of the parfile
			parfile="%s/%s.par" % (cep2.parfile_dir, psr2)
			if check_parfile(parfile, psr2): continue
			# with the B or J in the name
			parfile="%s/%s.par" % (cep2.parfile_dir, psr)
			if check_parfile(parfile, psr2): continue
			# checking extra directory withiout ^[BJ] in the name of the parfile
			parfile="%s/%s.par" % (cep2.parfile_dir_extra, psr2)
			if check_parfile(parfile, psr2): continue
			# checking extra directory with the B or J in the name
			parfile="%s/%s.par" % (cep2.parfile_dir_extra, psr)
			if check_parfile(parfile, psr2): continue

			self.log.info("Parfile does not exist. Creating parfile base on pulsar ephemeris from ATNF catalog...")
			# for -e option, we need to use pulsar name with leading B or J, otherwise it is possible (I cama across this!)
			# that there are two pulsars with the same name, one with leading J, another with leading B,
			# in this case psrcat returns records for both pulsars, and output parfile gets messed up
			if not cmdline.opts.is_globalfs:
				cmd="psrcat -db_file %s -e %s > %s/%s.par" % (cep2.psrcatdb, psr, self.outdir, psr2)
				self.execute(cmd, is_os=True)
			else:
				import tempfile
				tmpdir = tempfile.mkdtemp()
				parf = os.path.join(tmpdir, "%s.par" % (psr2))
				cmd="psrcat -db_file %s -e %s > %s" % (cep2.psrcatdb, psr, parf)
				self.execute(cmd, is_os=True)
				cmd="rsync -a %s %s" % (parf, self.outdir)
				self.execute(cmd)
				cmd="rm -rf %s" % (tmpdir)
				os.system(cmd)

		# To avoid potential simultaneous access (as probably happens on CEP4 with global file system)
		# each PipeUnit will have also copies of the parfiles inside the temporary directories
		# which will be deleted at the end
		import tempfile
		tmpdir = tempfile.mkdtemp()

		# Now we check the par-files if they have non-appropriate flags that can cause prepfold to crash
		toremove_psrs=[] # list of pulsars to remove from folding in case there is a problem with the parfile (e.g. PEPOCH is absent)
		for psr in self.psrs:
			psr2=re.sub(r'^[BJ]', '', psr)
			parfile="%s/%s.par" % (self.outdir, psr2)
			cmd="rsync -a %s %s" % (parfile, tmpdir)
			self.execute(cmd)
			parf="%s/%s" % (tmpdir, parfile.split("/")[-1])
			# check first that PEPOCH is in the parfile
			cmd="grep 'PEPOCH' %s" % (parf,)
			status=os.popen(cmd).readlines()
			if np.size(status)==0:
				self.log.warning("WARNING: Par-file %s has no PEPOCH keyword, so this pulsar will be excluded from folding." % (parfile,))
				toremove_psrs.append(psr)
				continue
			# check first for PSRB name. It should be changed to PSRJ
			cmd="grep 'PSRB' %s" % (parf,)
			status=os.popen(cmd).readlines()
			if np.size(status)>0:
				self.log.warning("WARNING: Par-file %s has PSRB keyword that can cause prepfold to crash!\n\
If your pipeline run calls prepfold you might need to change PSRB to PSRJ." % (parfile,))
			# checking CLK flag
			cmd="grep 'CLK' %s" % (parf,)
			status=os.popen(cmd).readlines()
			if np.size(status)>0:
				self.log.warning("WARNING: Par-file %s has CLK keyword that can cause prepfold to crash!\n\
CLK line will be removed from the parfile!" % (parfile,))
				cmd="sed -i '/^CLK/d' %s" % (parf,)
				self.execute(cmd, self.log, is_os=True)
				cmd="rsync -a %s %s" % (parf, self.outdir)
				self.execute(cmd)
			# checking for UNITS flag
			cmd="grep 'UNITS' %s" % (parf,)
			status=os.popen(cmd).readlines()
			if np.size(status)>0:
				self.log.warning("WARNING: Par-file %s has UNITS keyword that can cause prepfold to crash!\n\
UNITS line will be removed from the parfile!" % (parfile,))
				cmd="sed -i '/^UNITS/d' %s" % (parf,)
				self.execute(cmd, self.log, is_os=True)
				cmd="rsync -a %s %s" % (parf, self.outdir)
				self.execute(cmd)
		# removing pulsars from folding if their parfile is bad (e.g. no PEPOCH)
		if len(toremove_psrs) > 0:
			indices_to_remove=set()
			for psr in np.unique(toremove_psrs):
				for ii in xrange(len(self.psrs)):
					if psr == self.psrs[ii]: indices_to_remove.add(ii) 
			indices_to_keep=sorted(set(np.arange(len(self.psrs)))-indices_to_remove)
			self.psrs = self.psrs[indices_to_keep]
		# return tmpdir
		return tmpdir

	def execute(self, cmd, workdir=None, shell=False, is_os=False):
		"""
		Execute the command 'cmd' after logging the command
		and the wall-clock amount of time the command took to execute.
		This function waits for process to finish
		"""
		try:
			self.log.info(re.sub("\n", "\\\\n", cmd))  # also escaping \n to show it as it is
			job_start = time.time()
			self.log.info("Start at UTC %s" % (time.asctime(time.gmtime())))
			status = 1
			if is_os: status = os.system(cmd)
			else:
				process = Popen(shlex.split(cmd), stdout=PIPE, stderr=STDOUT, cwd=workdir, shell=shell, bufsize=1048576)
				self.procs.append(process)
				self.log.process2log(process)
				process.communicate()
				status=process.poll()
				self.procs.remove(process)
			job_end = time.time()
			job_total_time = job_end - job_start
			self.log.info("Finished at UTC %s, status=%s, Total running time: %.1f s (%.2f hrs)" % (time.asctime(time.gmtime()), status, job_total_time, job_total_time/3600.))
			self.log.info("")
			# if job is not successful
			if status != 0:
				raise Exception
		except Exception:
			self.log.exception("Oops... job has crashed!\n%s\nStatus=%s" % (re.sub("\n", "\\\\n", cmd), status))
			raise Exception

	def start_and_go(self, cmd, workdir=None, shell=False, immediate_status_check=False):
		"""
		Execute the command 'cmd' after logging the command
		This function start the cmd and leaves the function
		returning the Popen object, it does not wait for process to finish
		"""
		status=1
		try:
			self.log.info(re.sub("\n", "\\\\n", cmd))
			self.log.info("Start at UTC %s" % (time.asctime(time.gmtime())))
			if immediate_status_check:
				process = Popen(shlex.split(cmd), cwd=workdir, shell=shell, bufsize=1048576)
				time.sleep(10)  # waiting 10 secs to see if process crashes right away
				if process.poll() is not None and process.poll() != 0:
					raise Exception
				else: process.kill()  # if process is still running, it means that cmd is good, so we kill it in order to
			# restart it with proper stdout/stderr and add it to the list
			process = Popen(shlex.split(cmd), stdout=PIPE, stderr=STDOUT, cwd=workdir, shell=shell, bufsize=1048576)
			status=process.returncode
			self.procs.append(process)
			self.log.info("Job pid=%d, not waiting for it to finish." % (process.pid))
			return process
		except Exception:
			self.log.exception("Oops... job has crashed!\n%s\nStatus=%s" % (re.sub("\n", "\\\\n", cmd), status))
			raise Exception

	def waiting(self, prg, popen):
		"""
		Waiting for process to finish
		"""
		try:
			job_start = time.time()
			self.log.info("Waiting for %s to finish, pid=%d" % (prg, popen.pid))
			(sout, serr) = popen.communicate()
			# we pipe serr to sout, so no need to log serr
#			self.log.stdout2log(sout)
			self.log.info(sout)
			job_end = time.time()
			job_total_time = job_end - job_start
			self.log.info("Process pid=%d (%s) has finished at UTC %s, status=%d, Waiting time: %.1f s (%.2f hrs)" % \
				(popen.pid, prg, time.asctime(time.gmtime()), popen.returncode, job_total_time, job_total_time/3600.))
			self.log.info("")
			self.procs.remove(popen)
			# if job is not successful
			if popen.poll() != 0:
				raise Exception
		except Exception:
			self.log.exception("Oops... %s has crashed!\npid=%d, Status=%s" % (prg, popen.pid, popen.returncode))
			raise Exception

	def waiting_list(self, prg, popen_list):
		"""
		Waiting for list of processes to finish
		"""
		try:
			job_start = time.time()
			self.log.info("Waiting for %s processes to finish..." % (prg))
			run_units = [u.pid for u in popen_list if u.poll() is None]
			self.log.info("Still running [%d]: %s" % (len(run_units), run_units))
			for unit in popen_list:
				self.waiting(prg, unit)
				run_units = [u.pid for u in popen_list if u.poll() is None]
				finished_units = [u for u in popen_list if u.poll() is not None]
				for fu in finished_units:
					if fu.returncode != 0:
						self.log.exception("Oops... %s has crashed!\npid=%d, Status=%s" % (prg, fu.pid, fu.returncode))
				if len(run_units) > 0: self.log.info("Still running [%d]: %s" % (len(run_units), run_units))
			job_end = time.time()
			job_total_time = job_end - job_start
			self.log.info("Processes of %s have finished at UTC %s, Waiting time: %.1f s (%.2f hrs)" % (prg, time.asctime(time.gmtime()), job_total_time, job_total_time/3600.))
			self.log.info("")
		except Exception:
			self.log.exception("Oops... %s has crashed!\npids = %s" % (prg, ",".join(["%d" % (fu.pid) for fu in popen_list if fu.poll() is not None])))
			raise Exception

	"""
	def waiting_list(self, prg, popen_list):
	"""
	"""
	Waiting for list of processes to finish
	"""
	"""
		try:
			job_start = time.time()
			self.log.info("Waiting for %s processes to finish..." % (prg))
			run_popens = popen_list
			while True:
				if len(run_popens) == 0: break
				finished_units = [u for u in run_popens if u.poll() is not None]
				run_units = [u.pid for u in run_popens if u.poll() is None]
				run_popens = [u for u in run_popens if u.poll() is None]
				for fu in finished_units:
					if fu.returncode != 0:
						self.log.exception("Oops... %s has crashed!\npid=%d, Status=%s" % (prg, fu.pid, fu.returncode))
					else: self.waiting(prg, fu)
				if len(run_units) > 0 and len(finished_units) > 0: 
					self.log.info("Still running [%d]: %s" % (len(run_units), run_units))
			job_end = time.time()
			job_total_time = job_end - job_start
			self.log.info("Processes of %s have finished at UTC %s, Waiting time: %.1f s (%.2f hrs)" % \
				(prg, time.asctime(time.gmtime()), job_total_time, job_total_time/3600.))
			self.log.info("")
		except Exception:
			self.log.exception("Oops... %s has crashed!\npids = %s" % (prg, ",".join(["%d" % (fu.pid) for fu in popen_list if fu.poll() is not None])))
			raise Exception
	"""

	def power_of_two(self, value):
		"""
		Returns the closest power of two value to the input value (from the low side)
		"""
		return int(math.pow(2, math.floor(math.log(value)/math.log(2))))

	def get_best_nbins(self, parf):
		"""
		Calculates the best number of bins for folding based on the period value from the parfile (parf)
		and sampling interval (tsamp)
		"""
		try:
			f = open(parf, 'r')
			parlines = f.read().splitlines()
			f.close()
			res=[ii for ii in parlines if re.search("F0", ii) is not None]
			if len(res) > 0:
				f0=float(re.sub("\s+", " ", res[0]).split(" ")[1])
				nbins=self.power_of_two(int(math.ceil((1000.0/f0)/self.sampling)))
			else:
				res=[ii for ii in parlines if re.search("P0", ii) is not None]
				if len(res) > 0:
					p0=float(re.sub("\s+", " ", res[0]).split(" ")[1])
					nbins=self.power_of_two(int(math.ceil(p0*1000.0/self.sampling)))
				else:
					nbins=1024
			if nbins > 1024: return 1024
			else: return nbins
		except: return 1024

	def get_psr_dm(self, parf):
		"""
		Reads parfile and returns pulsar DM
		"""
		dm=0
		try:
			f = open(parf, 'r')
			parlines = f.read().splitlines()
			f.close()
			res=[ii for ii in parlines if re.search("DM", ii) is not None and re.search("DMEPOCH", ii) is None]
			if len(res) > 0:
				dm=float(re.sub("\s+", " ", res[0]).split(" ")[1])
		except: pass
		return dm

	def lcd(self, low, high):
		"""
		Calculates the lowest common denominator of 'high' value between 'low' and 'high'
		"""
		for ii in range(low, high+1): 
			if high % ii == 0: return ii
		return high

	def hcd(self, low, high, value):
		"""
		Calculates the highest common denominator of 'value' value between 'low' and 'high'
		"""
		for ii in range(high, low-1, -1): 
			if value % ii == 0: return ii
		return 1

	# function that checks all processes in the list and kill them if they are still running
	def kill(self):
		if self.log != None: self.log.info("Killing all open processes for SAP=%d TAB=%d..." % (self.sapid, self.tabid))
		for proc in self.procs:
			if proc.poll() is None: # process is still running
				proc.kill()
				if proc != None: proc.communicate()
				if proc != None: proc.poll()
		self.procs = []

	# refresh NFS mounting of locus node (loc) on hoover node
	# by doing 'ls' command
	def hoover_mounting(self, cep2, firstfile, loc):
		uniqdir="/".join(firstfile.split("/")[0:-1]).split("/data/")[-1]
		input_dir="%s/%s_data/%s" % (cep2.hoover_data_dir, loc, uniqdir)
		process = Popen(shlex.split("ls %s" % (input_dir)), stdout=PIPE, stderr=STDOUT)
		process.communicate()
		return input_dir

	# running prepfold(s) for all specified pulsars
	# returning the Popen instances of prepfold calls
	def start_prepfold(self, cmdline, total_chan):
		if total_chan >= 256:
			prepfold_nsubs = 512
		else: prepfold_nsubs = total_chan
#		prepfold_nsubs = total_chan
		self.log.info("Getting proper value of nsubs in prepfold between %d and %d..." % (prepfold_nsubs, total_chan))
		# calculating the least common denominator of total_chan starting from prepfold_nsubs
		prepfold_nsubs = self.lcd(prepfold_nsubs, total_chan)
		self.log.info("Number of subbands, -nsubs, for prepfold is %d" % (prepfold_nsubs))
		self.log.info("Running prepfolds...")
		prepfold_popens=[]  # list of prepfold Popen objects
		for psr in self.psrs:   # pulsar list is empty if --nofold is used
			psr2=re.sub(r'^[BJ]', '', psr)
			prepfold_nbins=self.get_best_nbins("%s/%s.par" % (self.outdir, psr2))
			# first running prepfold with mask (if --norfi was not set)
			if not cmdline.opts.is_norfi or os.path.exists("%s/%s_rfifind.mask" % (self.curdir, self.output_prefix)):
				# we use ../../../ instead of self.outdir for the full-name of the parfile, because in this case prepfold crashes
				# I suppose it happens because name of the file is TOO long for Tempo
				cmd="prepfold -noscales -nooffsets -noxwin -psr %s -par ../../../%s.par -n %d -nsub %d -fine -nopdsearch -mask %s_rfifind.mask -o %s_%s %s %s.fits" % \
					(psr, psr2, prepfold_nbins, prepfold_nsubs, self.output_prefix, psr, self.output_prefix, cmdline.opts.prepfold_extra_opts, self.output_prefix)
				try: # if prepfold fails right away (sometimes happens with error like this:
					# Read 0 set(s) of polycos for PSR 1022+1001 at 56282.138888888891 (DM = 3.6186e-317)
					# MJD 56282.139 out of range (    0.000 to     0.000)
					# isets = 0
					# I guess something to do with how Tempo treats parfile. When this happens, we try to rerun prepfold with the same
					# command but without using -par option
					prepfold_popen = self.start_and_go(cmd, workdir=self.curdir, immediate_status_check=True)
				except Exception:
					self.log.warning("***** Prepfold failed when using par-file. Will try the same command but without using -par option *****")
					cmd="prepfold -noscales -nooffsets -noxwin -psr %s -n %d -nsub %d -fine -nopdsearch -mask %s_rfifind.mask -o %s_%s %s %s.fits" % \
						(psr, prepfold_nbins, prepfold_nsubs, self.output_prefix, psr, self.output_prefix, cmdline.opts.prepfold_extra_opts, self.output_prefix)
					prepfold_popen = self.start_and_go(cmd, workdir=self.curdir)
				prepfold_popens.append(prepfold_popen)
				time.sleep(5) # will sleep for 5 secs, in order to give prepfold enough time to finish 
						# with temporary files lile resid2.tmp otherwise it can interfere with next prepfold call

			# running prepfold without mask
			if not cmdline.opts.is_norfi or os.path.exists("%s/%s_rfifind.mask" % (self.curdir, self.output_prefix)): 
				output_stem="_nomask"
			else: output_stem=""
			# we use ../../../ instead of self.outdir for the full-name of the parfile, because in this case prepfold crashes
			# I suppose it happens because name of the file is TOO long for Tempo
			cmd="prepfold -noscales -nooffsets -noxwin -psr %s -par ../../../%s.par -n %d -nsub %d -fine -nopdsearch -o %s_%s%s %s %s.fits" % \
				(psr, psr2, prepfold_nbins, prepfold_nsubs, psr, self.output_prefix, output_stem, cmdline.opts.prepfold_extra_opts, self.output_prefix)
			try: # same reasoning as above
				prepfold_popen = self.start_and_go(cmd, workdir=self.curdir, immediate_status_check=True)
			except Exception:
				self.log.warning("***** Prepfold failed when using par-file. Will try the same command but without using -par option *****")
				cmd="prepfold -noscales -nooffsets -noxwin -psr %s -n %d -nsub %d -fine -nopdsearch -o %s_%s%s %s %s.fits" % \
					(psr, prepfold_nbins, prepfold_nsubs, psr, self.output_prefix, output_stem, cmdline.opts.prepfold_extra_opts, self.output_prefix)
				prepfold_popen = self.start_and_go(cmd, workdir=self.curdir)
			prepfold_popens.append(prepfold_popen)
			time.sleep(5) # again will sleep for 5 secs, in order to give prepfold enough time to finish 
					# with temporary files like resid2.tmp otherwise it can interfere with next prepfold call
		return prepfold_popens

	# making prepfold diagnostic plots	
	def make_prepfold_plots(self, obs):

		# running convert on prepfold ps to pdf and png
		self.log.info("Running convert on ps to pdf and png of the plots...")
		prepfold_ps=glob.glob("%s/*.pfd.ps" % (self.curdir))
		for psfile in prepfold_ps:
			base=psfile.split("/")[-1].split(".ps")[0]
			cmd="convert %s.ps %s.pdf" % (base, base)
			# have separate try/except blocks for each convert command to allow to continue processing
			# in case something is wrong with ps-file
			try:
				self.execute(cmd, workdir=self.curdir)
			except: pass
			cmd="convert -rotate 90 %s.ps %s.png" % (base, base)
			try:
				self.execute(cmd, workdir=self.curdir)
			except : pass
			cmd="convert -rotate 90 -crop 200x140-0 %s.ps %s.th.png" % (base, base)
			try:
				self.execute(cmd, workdir=self.curdir)
			except: pass

		# getting the list of *.pfd.bestprof files and read chi-sq values for all folded pulsars
		psr_bestprofs=rglob(self.outdir, "*.pfd.bestprof", 3)
		if len(psr_bestprofs) > 0:
			self.log.info("Reading chi-squared values and adding to chi-squared.txt...")
			# also preparing montage command to create combined plot
			montage_cmd="montage -background none -pointsize 10.2 "

			# reading log file to get RFI fraction from rfifind
			logf = open(self.log.get_logfile(), 'r')
			rfifracs = [ff.split("(")[1].split("%")[0].lstrip() for ff in logf.read().splitlines() if re.search("Number of  bad", ff) is not None]
			logf.close()
			if np.size(rfifracs) == 0: rfifrac="" # in case we did not run rfifind
			elif np.size(rfifracs) > 1: rfifrac = rfifracs[-1]
			else: rfifrac = rfifracs[0]

			chif=open("%s/%s_sap%03d_tab%04d_stokes%d_part%03d_chi-squared.txt" % (self.outdir, obs.id, self.sapid, self.tabid, self.stokes_index, self.part), 'w')
			thumbs=[] # list of thumbnail files
			# check first if all available *bestprof files are those created without mask. If so, then allow to make
			# a diagnostic combined plot using prepfold plots without a mask
			good_bestprofs=[file for file in psr_bestprofs if re.search("_nomask_", file) is None]
			if len(good_bestprofs) == 0:
				good_bestprofs=[file for file in psr_bestprofs]
			for bp in good_bestprofs:
				psr=bp.split("/")[-1].split("_")[0]
				thumbfile=bp.split(self.outdir+"/")[-1].split(".pfd.bestprof")[0] + ".pfd.th.png"
				# getting current number for SAP and TA beam
				try: # we need this try block in case User manually creates sub-directories with some test bestprof files there
				     # which will also be found by rglob function. So, we need to exclude them by catching an Exception
				     # on a wrong-formed string applying int()
					cursapid=int(thumbfile.split("_SAP")[-1].split("_")[0])
				except: continue
				curprocdir=thumbfile.split("_SAP")[-1].split("_")[1]
				if self.sapid != cursapid or self.procdir != curprocdir:
					continue
				thumbs.append(thumbfile)
				# calculating the S/N (profile significance) from the bestprof file
				snr=0.0
				try:
					snrlog="snr-presto.log"
					if not os.path.exists("%s/%s" % (self.curdir, snrlog)):
						cmd="snr.py --presto --snrmethod=Off --auto-off --plot --saveonly %s | tee %s/%s" % (bp, self.curdir, snrlog)    
						self.execute(cmd, workdir=self.curdir, is_os=True)
						tmp = np.genfromtxt("%s/%s" % (self.curdir, snrlog), skip_header=13, skip_footer=2, usecols=(4,4), dtype=float, unpack=True)[0]
						snr = float(tmp[0])
				except:
					if self.sapid == cursapid and self.procdir == curprocdir:
						self.log.warning("Warning: can't read file %s or calculate S/N of the profile" % (bp))

				chif.write("file=%s obs=%s_SAP%d_%s_%s S/N=%g%s\n" % (thumbfile, self.code, cursapid, curprocdir, psr, snr, rfifrac != "" and " RFI=%s" % (rfifrac) or ""))
				montage_cmd += "-label '%s SAP%d %s\n%s\nS/N = %g' %s " % (self.code, cursapid, curprocdir, psr, snr, thumbfile)

			chif.close()
			cmd="mv %s_sap%03d_tab%04d_stokes%d_part%03d_chi-squared.txt chi-squared.txt" % (obs.id, self.sapid, self.tabid, self.stokes_index, self.part)
			self.execute(cmd, workdir=self.outdir)

			# creating combined plots
			# only creating combined plo % (snrtmpdir)ts when ALL corresponding thumbnail files exist. It is possible, when there are 2+ beams on
			# the same node, that bestprof files do exist, but thumbnails were not created yet at the time when chi-squared.txt is
			# getting created for another beam. And this will cause "montage" command to fail
			# At the end combined plot will eventually be created for this node during the procesing of the last beam of this node
			if len(thumbs) > 0 and len([ff for ff in thumbs if os.path.exists(ff)]) == len(thumbs):
				# creating combined plots
				self.log.info("Combining all pfd.th.png files in a single combined plot...")
				montage_cmd += "combined.png"
				self.execute(montage_cmd, workdir=self.outdir)
				# making thumbnail version of the combined plot
				cmd="convert -resize 200x140 -bordercolor none -border 150 -gravity center -crop 200x140-0-0 +repage combined.png combined.th.png"
				self.execute(cmd, workdir=self.outdir)


	# running single-pulse analysis (prepdata + single_pulse_search.py) for all specified pulsars
	def run_prepdata(self, cmdline):
		self.log.info("Running single-pulse analysis (prepdata + single_pulse_search.py for a pulsar(s) DM")
		for psr in self.psrs:   # pulsar list is empty if --nofold is used
			psr2=re.sub(r'^[BJ]', '', psr)
			if cmdline.opts.is_nofold:
				psrdm = 0.0
			else:
				psrdm=self.get_psr_dm("%s/%s.par" % (self.outdir, psr2))
			# running prepdata with mask (if --norfi was not set)
			if not cmdline.opts.is_norfi or os.path.exists("%s/%s_rfifind.mask" % (self.curdir, self.output_prefix)):
				cmd="prepdata -noscales -nooffsets -noclip -nobary -dm %f -mask %s_rfifind.mask -o %s_%s %s %s.fits" % \
					(psrdm, self.output_prefix, psr, self.output_prefix, cmdline.opts.prepdata_extra_opts, self.output_prefix)
			else: # running prepdata without mask
				cmd="prepdata -noscales -nooffsets -noclip -nobary -dm %f -o %s_%s %s %s.fits" % \
					(psrdm, psr, self.output_prefix, cmdline.opts.prepdata_extra_opts, self.output_prefix)
			self.execute(cmd, workdir=self.curdir)
			# running single_pulse_search.py on .dat file (created either with mask or not)
			cmd="single_pulse_search_lotaas.py --no-width-label %s_%s.dat" % (psr, self.output_prefix)
			self.execute(cmd, workdir=self.curdir)

	# Running "RRATs" analysis (prepsubband + single_pulse_search.py for a range of DMs)
	def run_rrats_analysis(self, cmdline, total_chan):
		self.log.info("Running RRATs analysis (prepsubband + single_pulse_search.py for a range of DMs)")
		# calculating the greatest common denominator of number pf channels from 1024 (highest possible for -nsub in prepsubband) down
		maxnsub = self.hcd(1, 1024, total_chan)
		for psr in self.psrs:   # pulsar list is empty if --nofold is used
			psr2=re.sub(r'^[BJ]', '', psr)
			if cmdline.opts.is_nofold:
				psrdm = 0.0
			else:
				psrdm=self.get_psr_dm("%s/%s.par" % (self.outdir, psr2))
			dmstep=0.01
			numdms=1000 # 1000 is the maximum
			# parsing prepsubband extra options for 'numdms', if there are no such option we will use default numdms=1000, if there
			# are several, we will use only the last one
			pattern = "\-numdms\s+(\d+)"
			numdms_list=re.findall(pattern, cmdline.opts.prepdata_extra_opts)
			if len(numdms_list) != 0: numdms = int(numdms_list[-1])
			# parsing prepsubband extra options for 'dmstep', if there are no such option we will use default dmstep=0.01, if there
			# are several, we will use only the last one
			pattern = "\-dmstep\s+([\d\.]+)"
			dmstep_list=re.findall(pattern, cmdline.opts.prepdata_extra_opts)
			if len(dmstep_list) != 0: dmstep = float(dmstep_list[-1])
			# parsing prepsubband extra options for 'lodm', if there are no such option we will calculate it as lodm=psrdm-0.5*dmstep*numdms
			# if there are several, we will only use the last one
			pattern = "\-lodm\s+([\d\.]+)"
			lodm_list=re.findall(pattern, cmdline.opts.prepdata_extra_opts)
			if len(lodm_list) != 0: 
				lodm = float(lodm_list[-1])
			else:
				lodm = psrdm - 0.5*dmstep*numdms # we want to cover +-5 in DM with 0.01 steps (by default), i.e. ~1000 DM trials
			if lodm <= 0.0: lodm = 0.01 # because we will do DM=0 with prepdata anyway
			# running prepdata for DM=0 with mask (if --norfi was not set)
			if not cmdline.opts.is_norfi or os.path.exists("%s/%s_rfifind.mask" % (self.curdir, self.output_prefix)):
				cmd="prepdata -noscales -nooffsets -noclip -nobary -dm 0.00 -mask %s_rfifind.mask -o %s_%s_DM0.00 %s %s.fits" % \
					(self.output_prefix, psr, self.output_prefix, cmdline.opts.prepdata_extra_opts, self.output_prefix)
			else: # running prepdata for DM=0 without mask
				cmd="prepdata -noscales -nooffsets -noclip -nobary -dm 0.00 -o %s_%s_DM0.00 %s %s.fits" % \
					(psr, self.output_prefix, cmdline.opts.prepdata_extra_opts, self.output_prefix)
			self.execute(cmd, workdir=self.curdir)

			# running prepsubband with mask (if --norfi was not set)
			if not cmdline.opts.is_norfi or os.path.exists("%s/%s_rfifind.mask" % (self.curdir, self.output_prefix)):
				cmd="prepsubband -noscales -nooffsets -noclip -nobary -nsub %d -mask %s_rfifind.mask -lodm %f -dmstep %f -numdms %d -o %s_%s %s %s.fits" % \
					(maxnsub, self.output_prefix, lodm, dmstep, numdms, psr, self.output_prefix, cmdline.opts.prepsubband_extra_opts, self.output_prefix)
			else: # running prepsubband without mask
				cmd="prepsubband -noscales -nooffsets -noclip -nobary -nsub %d -lodm %f -dmstep %f -numdms %d -o %s_%s %s %s.fits" % \
					(maxnsub, lodm, dmstep, numdms, psr, self.output_prefix, cmdline.opts.prepsubband_extra_opts, self.output_prefix)
			self.execute(cmd, workdir=self.curdir)

			# running single_pulse_search.py for *_DM*dat files
			cmd='single_pulse_search_lotaas.py -p -g "%s_%s_DM*.dat"' % (psr, self.output_prefix)
			self.execute(cmd, workdir=self.curdir)

			# running single_pulse_search_lotaas.py on .singlepulse file to make a plot (for all DMs)
			cmd='single_pulse_search_lotaas.py -t 5.5 --no-width-label -g "%s_%s_DM*.singlepulse"' % (psr, self.output_prefix)
			self.execute(cmd, workdir=self.curdir)

			# convert single-pulse ps-file (for all DMs) to png
			cmd="convert %s_%s_singlepulse.ps %s_%s_singlepulse.png" % (psr, self.output_prefix, psr, self.output_prefix)
			self.execute(cmd, workdir=self.curdir)

			# running single_pulse_search_lotaas.py on .singlepulse file to make a plot (for all DMs except DM=0)
			cmd='single_pulse_search_lotaas.py -t 5.5 --no-width-label --dms %f -g "%s_%s_DM*.singlepulse"' % (lodm, psr, self.output_prefix)
			self.execute(cmd, workdir=self.curdir)

			# convert single-pulse ps-file (for all DMs except DM=0) to png
			sp_ps=glob.glob("%s/%s_%s_DMs*_singlepulse.ps" % (self.curdir, psr, self.output_prefix))
			if len(sp_ps) != 0:
				cmd="convert %s %s_singlepulse.png" % (sp_ps[0].split("/")[-1], sp_ps[0].split("/")[-1].split("_singlepulse.ps")[0])
				self.execute(cmd, workdir=self.curdir)

			# removing *_DM*.dat files
			cmd="rm -f %s_%s_DM*.dat" % (psr, self.output_prefix)
			self.execute(cmd, workdir=self.curdir)

	# function that does last steps in processing, creating tarball, copyting it, changing ownership, etc...
	def finish_off(self, obs, cep2, cmdline, is_archive_to_summary=True, target_summary_dir=""):

		is_iquv = False
		if (self.tab.is_coherent and obs.stokesCS == "IQUV") or (self.tab.is_coherent == False and obs.stokesIS == "IQUV"):
			is_iquv = True

		# tarball name
		if is_iquv:
			tarname="%s_SAP%03d_B%03d_S%d_P%03d%s" % (cep2.pipeid, self.sapid, self.tabid, self.stokes_index, self.part, self.archive_suffix)
		else:
			tarname="%s_SAP%03d_B%03d_P%03d%s" % (cep2.pipeid, self.sapid, self.tabid, self.part, self.archive_suffix)
		# variables for rsync'ing the data
		verbose=""
		if cmdline.opts.is_debug:
			verbose="-v"
		if target_summary_dir=="":
			if not cmdline.opts.is_globalfs:
				output_dir="%s/%s_%s/%s%s" % \
					(self.processed_dir_root, cep2.processed_dir_prefix, self.summary_node, cmdline.opts.outdir == "" and cep2.pipeid or cmdline.opts.outdir, self.summary_node_dir_suffix)
			else:
				output_dir="%s/%s/%s%s" % \
					(self.processed_dir_root, cep2.processed_dir_prefix, cmdline.opts.outdir == "" and cep2.pipeid or cmdline.opts.outdir, self.summary_node_dir_suffix)
		else:
			output_dir=target_summary_dir
		output_archive="%s/%s" % (output_dir, tarname)

		# function that checks that there are no other stokes files, except the needed one
		# when we deal with IQUV data
		# should be called only when is_iquv is True
		def is_stokes_ok_for_the_tarball(tf, st):
			tfbase=tf.split("/")[-1]
			if "_S0" in tfbase or "_S1" in tfbase or "_S2" in tfbase or "_S3" in tfbase:
				if "_S%d" % (self.stokes_index) in tfbase: return True
				else: return False
			return True

		# function that filters out part-related files that belong to a different part 
		# (if for example there were more than 1 part on the same locus node)
		# returns True if file should stay on the list
		def is_file_ok_for_the_tarball(tf, part):
			tfbase=tf.split("/")[-1]
			if "_part" in tfbase:
				if "_part%d." % (self.part) in tfbase: return True
				else: return False
			if "_PART" in tfbase:
				if "_PART%d." % (self.part) in tfbase: return True
				else: return False
			if "_PL" in tfbase:
				for bl in tfbase.split("_PL"):
					if "_P" in bl:
						if "_PSR_" in bl and "_P" not in bl.split("_PSR_")[0] and "_P" not in bl.split("_PSR_")[1]: continue
						if "_P%03d" % (self.part) in bl: continue
						else: return False
			else:
				if "_P" in tfbase:
					if "_PSR_" in tfbase and "_P" not in tfbase.split("_PSR_")[0] and "_P" not in tfbase.split("_PSR_")[1]:
						return True
					if "_P%03d" % (self.part) in tfbase: return True
					else: return False
			return True

		if not cmdline.opts.is_feedback:
			if not cmdline.opts.is_cobalt:
				# copying parset file to output directory
				self.log.info("Copying original parset file to output directory...")
				cmd="cp -f %s %s" % (obs.parset, self.outdir)
				self.execute(cmd, workdir=self.outdir)
			if is_archive_to_summary and not cmdline.opts.is_globalfs:
				# Make a tarball of all the plots for this beam
				self.log.info("Making a tarball of all the files with extensions: %s" % (", ".join(self.extensions)))
				tar_list=[]
				for ext in self.extensions:
					ext_list=rglob(self.curdir, ext, 3)
					tar_list.extend(ext_list)
				tar_list.extend(glob.glob("%s/*.par" % (self.outdir)))
				tar_list.extend(glob.glob("%s/*.parset" % (self.outdir)))
				# filtering the list
				tar_list = [tf for tf in tar_list if is_file_ok_for_the_tarball(tf, self.part)]
				if is_iquv:
					tar_list = [tf for tf in tar_list if is_stokes_ok_for_the_tarball(tf, self.stokes_index)]
				cmd="tar -cv --ignore-failed-read -f %s %s" % (tarname, " ".join([f.split(self.outdir+"/")[1] for f in tar_list]))
				try: # --ignore-failed-read does not seem to help with tar failing for some beams
        		             # like file was changed during the tar, though tarball seem to be fine
					self.execute(cmd, workdir=self.outdir)
				except: pass

				# copying archive file to summary node
				if not cmdline.opts.is_globalfs:
					self.log.info("Copying archive file to %s:%s" % (self.summary_node, output_dir))
					cmd="rsync %s -axP --rsh='ssh -x' %s %s:%s" % (verbose, tarname, self.summary_node, output_archive)
					self.execute(cmd, workdir=self.outdir)

		# finish
		self.end_time=time.time()
		self.total_time=self.end_time - self.start_time
		self.log.info("UTC stop time is: %s" % (time.asctime(time.gmtime())))
		self.log.info("Total running time: %.1f s (%.2f hrs)" % (self.total_time, self.total_time/3600.))

		# flushing log file and copy it to outdir on local node and summary node
		self.log.flush()
		if not cmdline.opts.is_log_append: cmd="cp -f %s %s" % (cep2.get_logfile(), self.outdir)
		else: cmd="cat %s >> %s/%s" % (cep2.get_logfile(), self.outdir, cep2.get_logfile().split("/")[-1])
		os.system(cmd)
		if not cmdline.opts.is_globalfs:
			cmd="rsync %s -axP --rsh='ssh -x' %s %s:%s" % (verbose, cep2.get_logfile(), self.summary_node, output_dir)
		else:
			cmd="cp %s %s %s" % (verbose, cep2.get_logfile(), output_dir)
		proc = Popen(shlex.split(cmd), stdout=PIPE, stderr=STDOUT, cwd=self.outdir)
		proc.communicate()

		if not cmdline.opts.is_feedback:
			# create a new full archive tarball
			tar_list=[]
			for ext in self.full_archive_exts:
				ext_list=rglob(self.curdir, ext, 3)
				tar_list.extend(ext_list)
			tar_list.extend(glob.glob("%s/*.par" % (self.outdir)))
			tar_list.extend(glob.glob("%s/*.parset" % (self.outdir)))
			tar_list.extend(glob.glob("%s/%s" % (self.outdir, cep2.get_logfile().split("/")[-1])))
			# filtering the list
			tar_list = [tf for tf in tar_list if is_file_ok_for_the_tarball(tf, self.part)]
			if is_iquv:
				tar_list = [tf for tf in tar_list if is_stokes_ok_for_the_tarball(tf, self.stokes_index)]
			cmd="tar -cv --ignore-failed-read -f %s %s" % (tarname, " ".join([f.split(self.outdir+"/")[1] for f in tar_list]))
			try:
				# --ignore-failed-read does not seem to help with tar failing for some beams
				# like file was changed during the tar, though tarball seem to be fine
				proc = Popen(shlex.split(cmd), stdout=PIPE, stderr=STDOUT, cwd=self.outdir)
				proc.communicate()
			except: pass

		# updating tarname filename in the feedback unit
		self.fbunit.filename_update("%s/%s" % (self.outdir, tarname))

		# specific to dragnet
		if cep2.cluster_headnode == "dragnet":
			# changing ownership to 'dragnet' group
			cmd="chgrp -R dragnet %s" % (self.outdir)
			os.system(cmd)
		# changing the file permissions to be re-writable for group
		cmd="chmod -R g+w %s" % (self.outdir)
		os.system(cmd)

	# main processing function
	def run(self, obs, cep2, cmdline, log):
		# if there are no pulsars to fold we set --nofold option to True
		if len(self.psrs) == 0 and cmdline.opts.is_nofold == False:
			cmdline.opts.is_nofold = True

		# if we are not using dspsr directly to read *.h5 files
		if not cmdline.opts.is_with_dal:
			self.run_nodal(obs, cep2, cmdline, log)
		else:
			self.run_dal(obs, cep2, cmdline, log)

	# main run function to use 2bf2fits for conversion
	def run_nodal(self, obs, cep2, cmdline, log):
		try:
			self.log = log
			self.logfile = cep2.get_logfile().split("/")[-1]		
			self.start_time=time.time()	

			# start logging
			self.log.info("%s SAP=%d TAB=%d %s(%s%s Stokes: %s)    UTC start time is: %s  @node: %s" % \
				(obs.id, self.sapid, self.tabid, obs.FE and ", ".join(self.tab.stationList) + " " or "", obs.FE and "FE/" or "", self.code, \
				self.stokes, time.asctime(time.gmtime()), cep2.current_node))

			# Re-creating root output directory
			cmd="mkdir -m 775 -p %s" % (self.outdir)
			self.execute(cmd)

			# creating Par-file in the output directory or copying existing one
			if not cmdline.opts.is_nofold: 
				tmpdir = self.get_parfile(cmdline, cep2)

			# Creating curdir dir
			cmd="mkdir -m 775 -p %s" % (self.curdir)
			self.execute(cmd)

			is_iquv = False
			if (self.tab.is_coherent and obs.stokesCS == "IQUV") or (self.tab.is_coherent == False and obs.stokesIS == "IQUV"):
				is_iquv = True

			# if we run the whole processing and not just plots
			if not cmdline.opts.is_plots_only and not cmdline.opts.is_nodecode and not cmdline.opts.is_feedback:
				if (len(self.tab.location) > 1 or cmdline.opts.is_all_parts_at_once) and \
				not (is_iquv and ((self.tab.is_coherent and obs.nsplitsCS == 1) \
				or (not self.tab.is_coherent and obs.nsplitsIS == 1))): # it means we are on hoover nodes, so dir with the input data is different
				        		        	                # also we need to moint locus nodes that we need
					self.log.info("Re-mounting locus nodes on 'hoover' node %s: %s" % (cep2.current_node, " ".join(self.tab.location)))
					input_files=[]
					for loc in self.tab.location:
						if loc in self.tab.rawfiles.keys():
							# first "mounting" corresponding locus node
							input_dir = self.hoover_mounting(cep2, self.tab.rawfiles[loc][0], loc)
							if is_iquv:
								# in this case we get files only for needed Stokes
								input_file=["%s/%s" % (input_dir, f.split("/" + obs.id + "/")[-1]) for f in self.tab.rawfiles[loc] \
									if re.search("_S%d_" % (self.stokes_index), f)]
								# copy *.h5 files (want to keep them)
								for f in [ff for ff in self.tab.rawfiles[loc] if re.search("_S%d_" % (self.stokes_index), ff)]:
									cmd="cp -f %s/%s.h5 ." % (input_dir, f.split("/" + obs.id + "/")[-1].split(".raw")[0])
									try:
										self.execute(cmd, workdir=self.curdir)
									except: pass
							else:
								# in this case we get all available files
								input_file=["%s/%s" % (input_dir, f.split("/" + obs.id + "/")[-1]) for f in self.tab.rawfiles[loc]]
								# copy *.h5 files (want to keep them)
								for f in self.tab.rawfiles[loc]:
									cmd="cp -f %s/%s.h5 ." % (input_dir, f.split("/" + obs.id + "/")[-1].split(".raw")[0])
									try:
										self.execute(cmd, workdir=self.curdir)
									except: pass

							input_files.extend(input_file)
					input_file=" ".join(input_files)
				else:
					if not cmdline.opts.is_globalfs:	
						if cep2.current_node in self.tab.rawfiles.keys():
							if is_iquv:
								input_file=" ".join([f for f in self.tab.rawfiles[cep2.current_node] if re.search("_S%d_" % (self.stokes_index), f)])
								# copy *.h5 files
								for f in [ff for ff in self.tab.rawfiles[cep2.current_node] if re.search("_S%d_" % (self.stokes_index), ff)]:
									cmd="cp -f %s.h5 ." % (f.split(".raw")[0])
									try:
										self.execute(cmd, workdir=self.curdir)
									except: pass
							else:
								input_file=" ".join(self.tab.rawfiles[cep2.current_node])
								# copy *.h5 files
								for f in self.tab.rawfiles[cep2.current_node]:
									cmd="cp -f %s.h5 ." % (f.split(".raw")[0])
									try:
										self.execute(cmd, workdir=self.curdir)
									except: pass
						else: input_file=""
					else: # global FS
						input_files=[]
						for val in self.tab.rawfiles.values(): input_files.extend(val)
						if is_iquv:
							needed_files=[f for f in input_files if re.search("_S%d_" % (self.stokes_index), f)]
							input_file=" ".join(needed_files)
							# copy *.h5 files
							for f in needed_files:
								cmd="cp -f %s.h5 ." % (f.split(".raw")[0])
								try:
									self.execute(cmd, workdir=self.curdir)
								except: pass
						else:
							input_file=" ".join(input_files)
							# copy *.h5 files
							for f in input_files:
								cmd="cp -f %s.h5 ." % (f.split(".raw")[0])
								try:
									self.execute(cmd, workdir=self.curdir)
								except: pass

				self.log.info("Input data: %s" % (input_file))

			self.output_prefix="%s_SAP%d_%s" % (obs.id, self.sapid, self.procdir)
			if is_iquv: self.output_prefix+="_S%d" % (self.stokes_index)
			self.log.info("Output file(s) prefix: %s" % (self.output_prefix))

			total_chan = self.tab.nrSubbands*self.nrChanPerSub

			# if we run the whole processing and not just plots
			if not cmdline.opts.is_plots_only and not cmdline.opts.is_feedback:
				# converting raw data to 8 bits
				if not cmdline.opts.is_nodecode and cmdline.opts.is_raw_to_8bit:
					# creating directory for 8-bit raw data (raw-8bit)
					cmd="mkdir -m 775 -p %s/%s" % (self.curdir, self.raw_8bit_dir)
					self.execute(cmd)
					verbose=""
					if cmdline.opts.is_debug: verbose="-v"
					self.log.info("Converting raw 32-bit data to 8 bits...")
					cmd="python %s/digitize.py %s -s %g -o %s/%s %s" % (cep2.lofarsoft_bin, verbose, cmdline.opts.digitize_sigma, self.curdir, self.raw_8bit_dir, input_file.replace(".raw", ".h5"))
					self.execute(cmd, workdir="%s/%s" % (self.curdir, self.raw_8bit_dir))

				# running data conversion (2bf2fits)
				if not cmdline.opts.is_nodecode:
					verbose=""
					if cmdline.opts.is_debug: verbose="-v"
					if cmdline.opts.is_cobalt: # Cobalt-specific call without using parset file
						for ss in obs.saps:
							if ss.sapid == self.sapid:
								sap = ss
						nsamples=int(49152 / (self.nrChanPerSub * self.downsample_factor)) # this roughly corresponds to a ~250ms block
						# then we tweak this value to be sure it's common denominator of the total number of samples
						nsamples=self.lcd(nsamples, os.path.getsize(input_file)/(4*self.tab.nrSubbands*self.nrChanPerSub))
						cmd="2bf2fits %s %s -append -nbits 8 -A %d -sigma %d -nsubs %d -sap %d -tab %d -stokes %d -o %s \
-nsamples %d -nchans %d -ra %s -dec %s -psr %s -clock %d -band %s -startdate %s -starttime %s -samptime %g -duration %g -subs %s -obsid %s -observer %s %s %s" % \
(verbose, self.raw2fits_extra_options, cmdline.opts.decode_nblocks, cmdline.opts.decode_sigma, self.tab.nrSubbands, self.sapid, \
self.tabid, self.stokes_index, self.output_prefix, nsamples, self.nrChanPerSub, str(self.tab.rarad), \
str(self.tab.decrad), sap.source, obs.sampleClock, obs.bandFilter, obs.startdate, obs.starttime, \
self.sampling/1000., obs.duration, sap.subbandList, obs.id, obs.projectPI.split()[0].split(',')[0], cmdline.opts.bf2fits_extra_opts, input_file)
					else: # BG/P call using parset file
						cmd="2bf2fits -revend %s %s -parset %s -append -nbits 8 -A %d -sigma %d -nsubs %d -sap %d -tab %d -stokes %d -duration %g -o %s %s %s" % \
(verbose, self.raw2fits_extra_options, obs.parset, cmdline.opts.decode_nblocks, cmdline.opts.decode_sigma, self.tab.nrSubbands, self.sapid, self.tabid, self.stokes_index, obs.duration, \
self.output_prefix, cmdline.opts.bf2fits_extra_opts, input_file)
					self.execute(cmd, workdir=self.curdir)
					# fixing the coordinates
					#cmd="fix_fits_coordinates.py %s.fits" % (self.output_prefix)
					#self.execute(cmd, workdir=self.curdir)

				# running RFI excision, checking
				if not cmdline.opts.is_norfi:
					zapstr=""
					# we should zap 1st chan as after 2nd polyphase it has junk
					if self.nrChanPerSub > 1:
						zapstr="-zapchan 0:%d:%d" % (total_chan-1, self.nrChanPerSub)
					self.log.info("Creating RFI mask...")
					cmd="rfifind -o %s -psrfits -noclip -blocks 16 %s %s %s.fits" % (self.output_prefix, zapstr, cmdline.opts.rfifind_extra_opts, self.output_prefix)
					rfifind_popen = self.start_and_go(cmd, workdir=self.curdir)

					# start subdyn.py
					if not cmdline.opts.is_skip_subdyn:
						self.log.info("Producing RFI report...")
						samples_to_average=int((cmdline.opts.subdyn_time_average * 1000.) / self.sampling)
						cmd="python %s/subdyn.py --psrfits --saveonly -n %d --title \"%s %s Averaging: %g s \" %s.fits" % (cep2.lofarsoft_bin, samples_to_average, obs.id, self.procdir, cmdline.opts.subdyn_time_average, self.output_prefix)
						subdyn_popen = self.start_and_go(cmd, workdir=self.curdir)

					# waiting for rfifind to finish
					self.waiting("rfifind", rfifind_popen)

				if not cmdline.opts.is_nofold:	
					# running prepfold with and without mask
					if not cmdline.opts.is_skip_prepfold:
						prepfold_popens = self.start_prepfold(cmdline, total_chan)

					# now running dspsr stuff...
					if not cmdline.opts.is_skip_dspsr:

						verbose="-q"
						if cmdline.opts.is_debug: verbose="-v"
						dspsr_popens=[] # list of dspsr Popen objects
						for psr in self.psrs: # pulsar list is empty f --nofold is used
							psr2=re.sub(r'^[BJ]', '', psr)
							dspsr_nbins=self.get_best_nbins("%s/%s.par" % (tmpdir, psr2))
							cmd="dspsr -b %d -A %s -E %s/%s.par %s -O %s_%s -t %d %s %s.fits" % \
								(dspsr_nbins, self.dspsr_folding_options, tmpdir, psr2, verbose, psr, \
								self.output_prefix, cmdline.opts.nthreads, cmdline.opts.dspsr_extra_opts, self.output_prefix)
							dspsr_popen = self.start_and_go(cmd, workdir=self.curdir)
							dspsr_popens.append(dspsr_popen)

						# waiting for dspsr to finish
						self.waiting_list("dspsr", dspsr_popens)

						# fixing coordinates in the ar-file
						for psr in self.psrs:  # pulsar list is empty if --nofold is used
							fix_header_coords(self, self, psr, self.curdir, self.output_prefix)

						for psr in self.psrs:  # pulsar list is empty if --nofold is used
							self.log.info("Running post-dspsr processing for pulsar %s..." % (psr))

							# running common DSPSR post-processing
							dspsr_postproc(self, self, cmdline, obs, psr, total_chan, self.tab.nrSubbands, self.curdir, self.output_prefix)

			# running pav
			if not cmdline.opts.is_skip_dspsr and not cmdline.opts.is_feedback:
				make_dspsr_plots(self, self, cmdline, obs, self.tab.nrSubbands, self.curdir, self.output_prefix)

			# Starting pdmp
			if not cmdline.opts.is_plots_only and not cmdline.opts.is_feedback:
				if not cmdline.opts.is_skip_dspsr and not cmdline.opts.is_nofold and not cmdline.opts.is_nopdmp:
					pdmp_popens = start_pdmp(self, self, cmdline, obs, self.tab.nrSubbands, self.curdir, self.output_prefix)
		
				# waiting for prepfold to finish
				try:
					if not cmdline.opts.is_nofold and not cmdline.opts.is_skip_prepfold: 
						self.waiting_list("prepfold", prepfold_popens)
				# we let PULP to continue if prepfold has crashed, as dspsr can run ok, or other instances of prepfold could finish ok as well
				except: pass

			# if we want to run prepdata to create a time-series and make a list of TOAs
			try:
				if not cmdline.opts.is_plots_only and not cmdline.opts.is_feedback and cmdline.opts.is_single_pulse:	
					# for IQUV only running this for Stokes I
					if is_iquv:
						if self.stokes_index == 0: self.run_prepdata(cmdline)
					else:
						self.run_prepdata(cmdline)
			# we let PULP to continue if prepdata (or single_pulse_search.py) has crashed, as the rest of the pipeline can finish ok
			except: pass

			# making prepfold diagnostic plots
			if not cmdline.opts.is_nofold and not cmdline.opts.is_feedback:
				self.make_prepfold_plots(obs)

			# Running "RRATs" analysis (prepsubband + single_pulse_search.py for a range of DMs)
			try:
				if not cmdline.opts.is_plots_only and not cmdline.opts.is_feedback and cmdline.opts.is_rrats:
					# for IQUV only running this for Stokes I
					if is_iquv:
						if self.stokes_index == 0: self.run_rrats_analysis(cmdline, total_chan)
					else:
						self.run_rrats_analysis(cmdline, total_chan)
			except: pass

			# waiting for subdyn to finish
			if not cmdline.opts.is_plots_only and not cmdline.opts.is_feedback and \
			not cmdline.opts.is_norfi and not cmdline.opts.is_skip_subdyn:
				self.waiting("subdyn.py", subdyn_popen)

			# waiting for pdmp to finish
			if not cmdline.opts.is_plots_only and not cmdline.opts.is_feedback:
				if not cmdline.opts.is_skip_dspsr and not cmdline.opts.is_nopdmp and not cmdline.opts.is_nofold: 
					finish_pdmp(self, self, pdmp_popens, cmdline, obs, self.curdir, self.output_prefix)
					# making diagnostic plot after pdmp
					make_dspsr_plots(self, self, cmdline, obs, self.tab.nrSubbands, self.curdir, self.output_prefix, True)

			# Running spectar.py to calculate pulsar spectra
			if not cmdline.opts.is_skip_dspsr and not cmdline.opts.is_feedback:
				calc_psr_spectra(self, self, cmdline, obs, self.sapid, self.tabid, self.curdir, self.output_prefix)

			# finishing off the processing...
			self.finish_off(obs, cep2, cmdline)

			# flushing feedback file to disk (success)
			self.fbunit.flush(100, cep2, False)

			# removing tmpdir
			if not cmdline.opts.is_nofold: 
				cmd="rm -rf %s" % (tmpdir)
				os.system(cmd)

		except Exception:
			# flushing feedback file to disk first (error)
			self.fbunit.flush(0, cep2, False)
			self.log.exception("Oops... 'run_nodal' function for %s%s has crashed!" % (obs.FE and "FE/" or "", self.code))
			self.kill()
			raise

		# kill all open processes
		self.kill()
		# remove reference to PulpLogger class from processing unit
		self.log = None

	# main run function using dspsr to directly read *.h5 files
	def run_dal(self, obs, cep2, cmdline, log):
		try:
			self.log = log
			self.logfile = cep2.get_logfile().split("/")[-1]		
			self.start_time=time.time()	

			# start logging
			self.log.info("%s SAP=%d TAB=%d %s(%s%s Stokes: %s)    UTC start time is: %s  @node: %s" % \
				(obs.id, self.sapid, self.tabid, obs.FE and ", ".join(self.tab.stationList) + " " or "", obs.FE and "FE/" or "", self.code, \
				self.stokes, time.asctime(time.gmtime()), cep2.current_node))

			# Re-creating root output directory
			cmd="mkdir -m 775 -p %s" % (self.outdir)
			self.execute(cmd)

			# creating Par-file in the output directory or copying existing one
			if not cmdline.opts.is_nofold: 
				tmpdir = self.get_parfile(cmdline, cep2)

			# Creating curdir dir
			cmd="mkdir -m 775 -p %s" % (self.curdir)
			self.execute(cmd)

			is_iquv = False
			if (self.tab.is_coherent and obs.stokesCS == "IQUV") or (self.tab.is_coherent == False and obs.stokesIS == "IQUV"):
				is_iquv = True

			# if not just making plots...
			if not cmdline.opts.is_plots_only and not cmdline.opts.is_feedback and not cmdline.opts.is_nodecode:
				if (len(self.tab.location) > 1 or cmdline.opts.is_all_parts_at_once) and \
				not (is_iquv and ((self.tab.is_coherent and obs.nsplitsCS == 1) \
				or (not self.tab.is_coherent and obs.nsplitsIS == 1))): # it means we are on hoover nodes, so dir with the input data is different
				        		        	                # also we need to moint locus nodes that we need
					self.log.info("Re-mounting locus nodes on 'hoover' node %s: %s" % (cep2.current_node, " ".join(self.tab.location)))
					input_files=[]
					for loc in self.tab.location:
						if loc in self.tab.rawfiles.keys():
							# first "mounting" corresponding locus node
							input_dir = self.hoover_mounting(cep2, self.tab.rawfiles[loc][0], loc)
							# dspsr needs all polarizations S* files to be in the current directory together with h5 files,
							# so we have to make soft links to input files
							self.log.info("Making links to input files in the current directory...")
							if is_iquv:
								# in this case we get files only for needed Stokes
								infiles=[ff for ff in self.tab.rawfiles[loc] if re.search("_S%d_" % (self.stokes_index), ff)]
							else:
								infiles=self.tab.rawfiles[loc]
							for f in infiles:
								# links to the *.raw files
								cmd="ln -sf %s/%s ." % (input_dir, f.split("/")[-1])
								try:
									self.execute(cmd, workdir=self.curdir)
								except: pass
								# copy *.h5 files (want to keep them)
								cmd="cp -f %s/%s.h5 ." % (input_dir, f.split("/" + obs.id + "/")[-1].split(".raw")[0])
								try:
									self.execute(cmd, workdir=self.curdir)
								except: pass
							input_file=["%s.h5" % (f.split("/" + obs.id + "/")[-1]).split(".raw")[0] for f in infiles]
							input_files.extend(input_file)
				else:
					if not cmdline.opts.is_globalfs:
						if cep2.current_node in self.tab.rawfiles.keys():
							# make a soft links in the current dir (in order for processing to be consistent with the case when data are in many nodes)
							self.log.info("Making links to input files in the current directory...")
							if is_iquv:
								# in this case we get files only for needed Stokes
								infiles=[ff for ff in self.tab.rawfiles[cep2.current_node] if re.search("_S%d_" % (self.stokes_index), ff)]
							else:
								infiles=self.tab.rawfiles[cep2.current_node]
							for f in infiles:
								# links to the *.raw files
								cmd="ln -sf %s ." % (f)
								try:
									self.execute(cmd, workdir=self.curdir)
								except: pass
								# copy *.h5 files
								cmd="cp -f %s.h5 ." % (f.split(".raw")[0])
								try:
									self.execute(cmd, workdir=self.curdir)
								except: pass
							input_files=["%s.h5" % (f.split("/" + obs.id + "/")[-1]).split(".raw")[0] for f in infiles]
					else: # global FS
						self.log.info("Making links to input files in the current directory...")
						input_files=[]
						for val in self.tab.rawfiles.values(): input_files.extend(val)
						if is_iquv:
							# in this case we get files only for needed Stokes
							infiles=[ff for ff in input_files if re.search("_S%d_" % (self.stokes_index), ff)]
						else:
							infiles=input_files
						for f in infiles:
							# links to the *.raw files
							cmd="ln -sf %s ." % (f)
							try:
								self.execute(cmd, workdir=self.curdir)
							except: pass
							# copy *.h5 files
							cmd="cp -f %s.h5 ." % (f.split(".raw")[0])
							try:
								self.execute(cmd, workdir=self.curdir)
							except: pass
						input_files=["%s.h5" % (f.split("/" + obs.id + "/")[-1]).split(".raw")[0] for f in infiles]

				self.log.info("Input data: %s" % ("\n".join(input_files)))

			self.output_prefix="%s_SAP%d_%s" % (obs.id, self.sapid, self.procdir)
			if is_iquv:
				self.output_prefix+="_S%d" % (self.stokes_index)
			self.log.info("Output file(s) prefix: %s" % (self.output_prefix))

			if self.tab.is_coherent:
				proc_subs = self.nrSubsPerFile * cmdline.opts.nsplitsCS
				if cmdline.opts.first_freq_splitCS * self.nrSubsPerFile + proc_subs > self.tab.nrSubbands:
					proc_subs -= (cmdline.opts.first_freq_splitCS * self.nrSubsPerFile + proc_subs - self.tab.nrSubbands)  
			else:
				proc_subs = self.nrSubsPerFile * cmdline.opts.nsplitsIS
				if cmdline.opts.first_freq_splitIS * self.nrSubsPerFile + proc_subs > self.tab.nrSubbands:
					proc_subs -= (cmdline.opts.first_freq_splitIS * self.nrSubsPerFile + proc_subs - self.tab.nrSubbands)  
			nsubs_eff = min(self.tab.nrSubbands, proc_subs)
			total_chan = nsubs_eff * self.nrChanPerSub

			if not cmdline.opts.is_plots_only and not cmdline.opts.is_feedback:
				# converting raw data to 8 bits
				if not cmdline.opts.is_nodecode and cmdline.opts.is_raw_to_8bit:
					# creating directory for 8-bit raw data (raw-8bit)
					cmd="mkdir -m 775 -p %s/%s" % (self.curdir, self.raw_8bit_dir)
					self.execute(cmd)
					verbose=""
					if cmdline.opts.is_debug: verbose="-v"
					self.log.info("Converting raw 32-bit data to 8 bits...")
					cmd="python %s/digitize.py %s -s %g -o %s/%s %s" % (cep2.lofarsoft_bin, verbose, cmdline.opts.digitize_sigma, self.curdir, self.raw_8bit_dir, " ".join(["%s/%s" % (self.curdir, ff) for ff in input_files]))
					self.execute(cmd, workdir="%s/%s" % (self.curdir, self.raw_8bit_dir))

				# getting the list of S? files, the number of which is how many freq splits we have
				# we also sort this list by split number
				s0_files=sorted(input_files, key=lambda input_file: int(input_file.split("_P")[-1].split("_")[0]))
				if not cmdline.opts.is_nodecode:
					verbose="-q"
					if cmdline.opts.is_debug: verbose="-v"

					# running dspsr for every pulsar for all frequency splits
					self.log.info("Running processing for all frequency splits...")
					# loop on pulsars
					for psr in self.psrs: # pulsar list is empty if --nofold is used
						psr2=re.sub(r'^[BJ]', '', psr)
						dspsr_nbins=self.get_best_nbins("%s/%s.par" % (tmpdir, psr2))
						if not cmdline.opts.is_nofold and not cmdline.opts.is_skip_dspsr:
							self.log.info("Running dspsr for pulsar %s..." % (psr))
						# loop on frequency splits
						for ii in range(len(s0_files)):
							# refreshing NFS mounting of locus nodes
							for loc in self.tab.location:
								if loc in self.tab.rawfiles:
									self.hoover_mounting(cep2, self.tab.rawfiles[loc][0], loc)
							fpart=int(s0_files[ii].split("_P")[-1].split("_")[0])
							if not cmdline.opts.is_nofold and not cmdline.opts.is_skip_dspsr:
								cmd="dspsr -b %d -A %s %s -E %s/%s.par -O %s_%s_P%d -t %d %s %s" % \
									(dspsr_nbins, self.dspsr_folding_options, verbose, tmpdir, psr2, \
									psr, self.output_prefix, fpart, cmdline.opts.nthreads, cmdline.opts.dspsr_extra_opts, s0_files[ii])
								self.execute(cmd, workdir=self.curdir)
							# run digifil with coherent dedispersion for further single-pulse analysis
							# works only for CV data and we need to run it for every pulsar (different DMs)
							if cmdline.opts.is_single_pulse and self.code=="CV":
								try:
									if not cmdline.opts.is_nofold:
										psrdm=self.get_psr_dm("%s/%s.par" % (tmpdir, psr2))
									else: psrdm = 0.0
									self.log.info("Running digifil for pulsar %s... (split %d)" % (psr, ii))
									cmd="digifil %s -B 512 -b 8 -F %d:D -D %f -o %s_%s_P%d.fil %s %s" % \
										(verbose, total_chan > 1 and total_chan or 2 * total_chan, psrdm, psr, self.output_prefix, fpart, \
										cmdline.opts.digifil_extra_opts, s0_files[ii])
									self.execute(cmd, workdir=self.curdir)
								# we let PULP to continue if digifil has crashed, as the rest of the pipeline can finish ok
								except: pass

						# running psradd to add all freq channels together
						if not cmdline.opts.is_nofold and not cmdline.opts.is_skip_dspsr:
							self.log.info("Adding frequency channels together...")
							ar_files=glob.glob("%s/%s_%s_P*.ar" % (self.curdir, psr, self.output_prefix))
							cmd="psradd -R -m time -o %s_%s.ar %s" % (psr, self.output_prefix, " ".join(ar_files))
							self.execute(cmd, workdir=self.curdir)

							# fixing coordinates in the ar-file
							fix_header_coords(self, self, psr, self.curdir, self.output_prefix)

							# running common DSPSR post-processing
							dspsr_postproc(self, self, cmdline, obs, psr, total_chan, nsubs_eff, self.curdir, self.output_prefix)
							# removing ar-files from dspsr for every frequency split
							if not cmdline.opts.is_debug:
								remove_list=glob.glob("%s/%s_%s_P*.ar" % (self.curdir, psr, self.output_prefix))
								cmd="rm -f %s" % (" ".join(remove_list))
								self.execute(cmd, workdir=self.curdir)

						# adding frequency parts together (if more than 1)
						# and running rfifind on fil-file if needed
						if cmdline.opts.is_single_pulse and self.code=="CV":
							try:
								# we sort the list of fil-files in decrease frequency order (needed for sigproc_splice)
								fil_files=sorted(glob.glob("%s/%s_%s_P*.fil" % (self.curdir, psr, self.output_prefix)), key=lambda x: int(x.split("_PART")[-1].split(".fil")[0]), reverse=True)
								if len(fil_files) > 1:
									self.log.info("Adding frequency channels together for single-pulse analysis...")
									cmd="sigproc_splice -a -o %s_%s.fil %s" % (psr, self.output_prefix, " ".join(fil_files))
									self.execute(cmd, workdir=self.curdir)
									if not cmdline.opts.is_debug:
										cmd="rm -f %s" % (" ".join(fil_files))
										self.execute(cmd, workdir=self.curdir)
								else:   # in case there is only one split, just rename it
									if len(fil_files) == 1:
										cmd="mv -f %s %s_%s.fil" % (fil_files[0], psr, self.output_prefix)
										self.execute(cmd, workdir=self.curdir)
								# running rfifind
								if not cmdline.opts.is_norfi:
									self.log.info("Running rfifind for pulsar %s on a filterbank file...")
									cmd="rfifind -o %s_%s -noclip -blocks 16 %s %s_%s.fil" % (psr, self.output_prefix, cmdline.opts.rfifind_extra_opts, psr, self.output_prefix)
									self.execute(cmd, workdir=self.curdir)
							except: pass

					# removing links for input .raw files
					if not cmdline.opts.is_debug:
						cmd="rm -f %s" % (" ".join(["%s.raw" % (f.split(".h5")[0]) for f in input_files]))
						self.execute(cmd, workdir=self.curdir)

			# running pav
			if not cmdline.opts.is_skip_dspsr and not cmdline.opts.is_feedback:
				make_dspsr_plots(self, self, cmdline, obs, nsubs_eff, self.curdir, self.output_prefix)

			# Running pdmp
			if not cmdline.opts.is_plots_only and not cmdline.opts.is_feedback and \
			not cmdline.opts.is_skip_dspsr and not cmdline.opts.is_nopdmp and not cmdline.opts.is_nofold:
				pdmp_popens = start_pdmp(self, self, cmdline, obs, nsubs_eff, self.curdir, self.output_prefix)
				# waiting for pdmp to finish
				finish_pdmp(self, self, pdmp_popens, cmdline, obs, self.curdir, self.output_prefix)
				# making diagnostic plot after pdmp
				make_dspsr_plots(self, self, cmdline, obs, nsubs_eff, self.curdir, self.output_prefix, True)

			# running single-pulse analysis
			if not cmdline.opts.is_plots_only and not cmdline.opts.is_feedback and cmdline.opts.is_single_pulse and self.code=="CV":	
				run_prepdata_CV(self, self, cmdline, self.outdir, self.curdir, self.output_prefix)

			# Running spectar.py to calculate pulsar spectra
			if not cmdline.opts.is_skip_dspsr and not cmdline.opts.is_feedback:
				calc_psr_spectra(self, self, cmdline, obs, self.sapid, self.tabid, self.curdir, self.output_prefix)

			# finishing off the processing...
			self.finish_off(obs, cep2, cmdline)

			# flushing feedback file to disk (success)
			self.fbunit.flush(100, cep2, False)

			# removing tmpdir
			if not cmdline.opts.is_nofold: 
				cmd="rm -rf %s" % (tmpdir)
				os.system(cmd)

		except Exception:
			# flushing feedback file to disk first (error)
			self.fbunit.flush(0, cep2, False)
			self.log.exception("Oops... 'run_dal' function for %s%s has crashed!" % (obs.FE and "FE/" or "", self.code))
			self.kill()
			raise

		# kill all open processes
		self.kill()
		# remove reference to PulpLogger class from processing unit
		self.log = None


class CSUnit(PipeUnit):
	def __init__(self, obs, cep2, cmdline, tab, log, curstokes, part=0):
		PipeUnit.__init__(self, obs, cep2, cmdline, tab, log, curstokes, part)
		self.code = "CS"
		self.stokes = obs.stokesCS[curstokes]
		self.beams_root_dir = "stokes"
		self.raw2fits_extra_options="-CS -H"
		self.nrChanPerSub = obs.nrChanPerSubCS
		self.nrSubsPerFile = obs.nrSubsPerFileCS
		self.sampling = obs.samplingCS
		self.downsample_factor = obs.downsample_factorCS
		self.summary_node = cep2.summary_nodes[self.code]
		if cep2.cluster_headnode[:5] == "head0" or cep2.cluster_headnode[:3] == "cpu": # CEP4
			self.summary_node_dir_suffix = "/pulp/cs" # "_CSplots"
			self.outdir_suffix = "/pulp/cs" # "_red"
		else:
			self.summary_node_dir_suffix = "_CSplots" # "_CSplots"
			self.outdir_suffix = "_red" # "_red"
		self.archive_suffix = "_bf.tar" # "_bf.tar.gz" # "_pulpCS.tar.gz"
		# setting outdir and curdir directories
		self.set_outdir(obs, cep2, cmdline)
		# initialiaze feedback unit
		self.feedback_init(obs, cep2, cmdline)

class ISUnit(PipeUnit):
	def __init__(self, obs, cep2, cmdline, tab, log, curstokes, part=0):
		PipeUnit.__init__(self, obs, cep2, cmdline, tab, log, curstokes, part)
		self.code = "IS"
		self.stokes = obs.stokesIS[curstokes]
		self.beams_root_dir = "incoherentstokes"
		self.raw2fits_extra_options = "-CS -H -IS"
		self.nrChanPerSub = obs.nrChanPerSubIS
		self.nrSubsPerFile = obs.nrSubsPerFileIS
		self.sampling = obs.samplingIS
		self.downsample_factor = obs.downsample_factorIS
		self.summary_node = cep2.summary_nodes[self.code]
		if cep2.cluster_headnode[:5] == "head0" or cep2.cluster_headnode[:3] == "cpu": # CEP4
			self.summary_node_dir_suffix = "/pulp/is" # "_redIS"
			self.outdir_suffix = "/pulp/is" # "_redIS"
		else:
			self.summary_node_dir_suffix = "_redIS" # "_redIS"
			self.outdir_suffix = "_redIS" # "_redIS"
		self.archive_suffix = "_bf.tar" # "_bf.tar.gz" # "_pulpIS.tar.gz"
		# setting outdir and curdir directories
		self.set_outdir(obs, cep2, cmdline)
		# initialiaze feedback unit
		self.feedback_init(obs, cep2, cmdline)

class FE_CSUnit(PipeUnit):
	def __init__(self, obs, cep2, cmdline, tab, log, curstokes, part=0):
		PipeUnit.__init__(self, obs, cep2, cmdline, tab, log, curstokes, part)
		self.code = "CS"
		self.stokes = obs.stokesCS[curstokes]
		self.beams_root_dir = "stokes"
		self.raw2fits_extra_options="-CS -H"
		self.nrChanPerSub = obs.nrChanPerSubCS
		self.nrSubsPerFile = obs.nrSubsPerFileCS
		self.sampling = obs.samplingCS
		self.downsample_factor = obs.downsample_factorCS
		self.summary_node = cep2.summary_nodes[self.code]
		if cep2.cluster_headnode[:5] == "head0" or cep2.cluster_headnode[:3] == "cpu": # CEP4
			self.summary_node_dir_suffix = "/pulp/cs" # "_CSplots"
			self.outdir_suffix = "/pulp/cs" # "_red"
		else:
			self.summary_node_dir_suffix = "_CSplots" # "_CSplots"
			self.outdir_suffix = "_red" # "_red"
		self.archive_suffix = "_bf.tar" # "_bf.tar.gz" # "_pulpCS.tar.gz"
		# re-assigning procdir from BEAMN to station name
		if obs.FE and self.tab.stationList[0] != "": 
			self.procdir = self.tab.stationList[0]
		# setting outdir and curdir directories
		self.set_outdir(obs, cep2, cmdline)
		# initialiaze feedback unit
		self.feedback_init(obs, cep2, cmdline)
