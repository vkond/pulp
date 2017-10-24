###################################################################
#
# Class CVUnit - processsing unit for CV data
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
from pulp_pipeunit import PipeUnit, dspsr_postproc, make_dspsr_plots, start_pdmp, finish_pdmp

class CVUnit(PipeUnit):
	def __init__(self, obs, cep2, cmdline, tab, log, curstokes, part=0):
		PipeUnit.__init__(self, obs, cep2, cmdline, tab, log, curstokes, part)
		self.code = "CV"
		self.stokes = obs.stokesCS
		self.beams_root_dir = "rawvoltages"
		self.raw2fits_extra_options=""
		self.nrChanPerSub = obs.nrChanPerSubCS
		self.nrSubsPerFile = obs.nrSubsPerFileCS
		self.sampling = obs.samplingCS
		self.downsample_factor = obs.downsample_factorCS
		self.summary_node = cep2.summary_nodes[self.code]
		if cep2.cluster_headnode[:5] == "head0" or cep2.cluster_headnode[:3] == "cpu": # CEP4
			self.summary_node_dir_suffix = "/pulp/cs"  # "_CVplots"
			self.outdir_suffix = "/pulp/cs"  # "_red"
		else:
			self.summary_node_dir_suffix = "_CVplots"  # "_CVplots"
			self.outdir_suffix = "_red"  # "_red"
		self.archive_suffix = "_bf.tar" # "_bf.tar.gz" # "_pulpCV.tar.gz"
		# extensions of the files to copy to archive (parfile and parset will be also included)
		self.extensions=["*.pdf", "*.ps", "*png", "*.ar", "*.AR", "*pdmp*", "*.rv", "*.out", "*.h5", "*.per", "*.posn", "*.txt"]
		# setting outdir and curdir directories
		self.set_outdir(obs, cep2, cmdline)
		# initialiaze feedback unit
		self.feedback_init(obs, cep2, cmdline)

	# main CV processing function
	def run(self, obs, cep2, cmdline, log):
		# if there are no pulsars to fold we set --nofold option to True
		if len(self.psrs) == 0 and cmdline.opts.is_nofold == False:
			cmdline.opts.is_nofold = True
		# if we are not using dspsr directly to read *.h5 files
		if cmdline.opts.is_nodal:
			self.run_nodal(obs, cep2, cmdline, log)
		else:
			self.run_dal(obs, cep2, cmdline, log)


	# main run function using bf2puma2 for conversion
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

			# if not just making plots...
			if not cmdline.opts.is_plots_only and not cmdline.opts.is_feedback and not cmdline.opts.is_nodecode:
				if len(self.tab.location) > 1: # it means we are on hoover nodes, so dir with the input data is different
        	        	                          # also we need to moint locus nodes that we need
					self.log.info("Re-mounting locus nodes on 'hoover' node %s: %s" % (cep2.current_node, " ".join(self.tab.location)))
					input_files=[]
					for loc in self.tab.location:
						if loc in self.tab.rawfiles:
							# first "mounting" corresponding locus node
							input_dir = self.hoover_mounting(cep2, self.tab.rawfiles[loc][0], loc)
							# bf2puma2 assumes all polarizations S* files to be in the current directory, so we have to
							# make soft links to input files
							self.log.info("Making links to input files in the current directory...")
							for f in self.tab.rawfiles[loc]:
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
							input_file=["%s" % (f.split("/" + obs.id + "/")[-1]) for f in self.tab.rawfiles[loc]]
							input_files.extend(input_file)
				else:
					if cep2.current_node in self.tab.rawfiles:
						# make a soft links in the current dir (in order for processing to be consistent with the case when data are in many nodes)
						self.log.info("Making links to input files in the current directory...")
						for f in self.tab.rawfiles[cep2.current_node]:
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
						input_files=["%s" % (f.split("/" + obs.id + "/")[-1]) for f in self.tab.rawfiles[cep2.current_node]]

				self.log.info("Input data: %s" % ("\n".join(input_files)))

			self.output_prefix="%s_SAP%d_%s" % (obs.id, self.sapid, self.procdir)
			self.log.info("Output file(s) prefix: %s" % (self.output_prefix))

			#total_chan = self.tab.nrSubbands * self.nrChanPerSub
			proc_subs = self.nrSubsPerFile * cmdline.opts.nsplitsCS
			if cmdline.opts.first_freq_splitCS * self.nrSubsPerFile + proc_subs > self.tab.nrSubbands:
				proc_subs -= (cmdline.opts.first_freq_splitCS * self.nrSubsPerFile + proc_subs - self.tab.nrSubbands)  
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
					cmd="python %s/digitize.py %s -s %g -o %s/%s %s" % (cep2.lofarsoft_bin, verbose, cmdline.opts.digitize_sigma, self.curdir, self.raw_8bit_dir, " ".join(["%s/%s" % (self.curdir, ff.replace(".raw", ".h5")) for ff in input_files]))
					self.execute(cmd, workdir="%s/%s" % (self.curdir, self.raw_8bit_dir))

				# getting the list of "_S0_" files, the number of which is how many freq splits we have
				# we also sort this list by split number
				s0_files=sorted([f for f in input_files if re.search("_S0_", f) is not None], key=lambda input_file: int(input_file.split("_P")[-1].split("_")[0]))
				if not cmdline.opts.is_nodecode:
					# checking if extra options for complex voltages processing were set in the pipeline
					nblocks=""
					if cmdline.opts.nblocks != -1: nblocks = "-b %d" % (cmdline.opts.nblocks)
					is_all_for_scaling=""
					if cmdline.opts.is_all_times: is_all_for_scaling="-all_times"
					is_write_ascii=""
					if cmdline.opts.is_write_ascii: is_write_ascii="-t"
					verbose=""
					if cmdline.opts.is_debug: verbose="-v"

					# running bf2puma2 command for all frequency splits...
					self.log.info("Running bf2puma2 for all frequency splits...")
					# loop on frequency splits
					for ii in range(len(s0_files)):
						# refreshing NFS mounting of locus nodes
						for loc in self.tab.location: 
							if loc in self.tab.rawfiles: 
								self.hoover_mounting(cep2, self.tab.rawfiles[loc][0], loc)
						cmd="bf2puma2 -f %s -h %s -p %s -hist_cutoff %f %s %s %s %s" % (s0_files[ii], cep2.puma2header, obs.parset, cmdline.opts.hist_cutoff, verbose, nblocks, is_all_for_scaling, is_write_ascii)
						self.execute(cmd, workdir=self.curdir)

				# removing links for input files
				if not cmdline.opts.is_debug:
					cmd="rm -f %s" % (" ".join(input_files))
					self.execute(cmd, workdir=self.curdir)

				if not cmdline.opts.is_nofold and not cmdline.opts.is_skip_dspsr:
					# getting the MJD of the observation
					# for some reason dspsr gets wrong MJD without using -m option
					input_prefix=s0_files[0].split("_S0_")[0]
					bf2puma_outfiles=glob.glob("%s/%s_SB*_CH*" % (self.curdir, input_prefix))
					self.log.info("Reading MJD of observation from the header of output bf2puma2 files...")
					cmd="head -25 %s | grep MJD" % (bf2puma_outfiles[0])
					mjdline=os.popen(cmd).readlines()
					if np.size(mjdline)>0:
						obsmjd=mjdline[0][:-1].split(" ")[-1]
						self.log.info("MJD = %s" % (obsmjd))
					else:
						self.log.error("Can't read the header of file '%s' to get MJD" % (bf2puma_outfiles[0]))
						self.kill()
						raise Exception

					verbose="-q"
					if cmdline.opts.is_debug: verbose="-v"

					# running dspsr for every frequency channel. We run it in bunches of number of channels in subband
					# usually it is 16 which is less than number of cores in locus nodes. But even if it is 32, then it should be OK (I think...)
					self.log.info("Running dspsr for each frequency channels...")
					# loop on pulsars
					for psr in self.psrs:
						self.log.info("Running dspsr for pulsar %s..." % (psr))
						psr2=re.sub(r'^[BJ]', '', psr)
						psrdm=self.get_psr_dm("%s/%s.par" % (tmpdir, psr2))
						dspsr_nbins=self.get_best_nbins("%s/%s.par" % (tmpdir, psr2))
						for bb in range(0, len(bf2puma_outfiles), self.nrChanPerSub):
							self.log.info("For %s" % (", ".join(bf2puma_outfiles[bb:bb+self.nrChanPerSub])))
							dspsr_popens=[] # list of dspsr Popen objects
							for cc in range(bb, bb+self.nrChanPerSub):
								input_file=bf2puma_outfiles[cc]
								cmd="dspsr -m %s -b %d -A %s %s -E %s/%s.par -O %s_%s_SB%s -t %d %s %s" % \
									(obsmjd, dspsr_nbins, self.dspsr_folding_options, verbose, tmpdir, psr2, psr, self.output_prefix, \
										input_file.split("_SB")[1], cmdline.opts.nthreads, cmdline.opts.dspsr_extra_opts, input_file)
								dspsr_popen = self.start_and_go(cmd, workdir=self.curdir)
								dspsr_popens.append(dspsr_popen)

								# running the single-pulse analysis
								if cmdline.opts.is_single_pulse:
									try:
										# getting coordinates string of the pulsar
										psr_ra=pu.coord_to_string(*pu.rad_to_hms(self.tab.rarad))
										psr_dec=pu.coord_to_string(*pu.rad_to_dms(self.tab.decrad))
										if self.tab.decrad < 0: psr_coords=psr_ra+psr_dec
										else: psr_coords=psr_ra+"+"+psr_dec
										# run digifil with coherent dedispersion for further single-pulse analysis
										# works only for CV data and we need to run it for every pulsar (different DMs)
										# also contrary to DAL case, here we need to specify correct site and pulsar coordinates when using Puma2-format input files
										cmd="digifil -set site=LOFAR -set name=%s -set coord=%s %s -F 2:D -D %f -b 8 -o %s_%s_SB%s.fil %s %s" % \
											(psr, psr_coords, verbose, psrdm, psr, self.output_prefix, input_file.split("_SB")[1], cmdline.opts.digifil_extra_opts, input_file)
										dspsr_popen = self.start_and_go(cmd, workdir=self.curdir)
										dspsr_popens.append(dspsr_popen)
									except: pass

							# waiting for dspsr to finish
							self.waiting_list("dspsr", dspsr_popens)

						# running psradd to add all freq channels together
						self.log.info("Adding frequency channels together...")
						ar_files=glob.glob("%s/%s_%s_SB*_CH*.ar" % (self.curdir, psr, self.output_prefix))
						cmd="psradd -R -m time -o %s_%s.ar %s" % (psr, self.output_prefix, " ".join(ar_files))
						self.execute(cmd, workdir=self.curdir)

                                                # fixing coordinates in the ar-file
						fix_header_coords(self, self, psr, self.curdir, self.output_prefix)

						# running common DSPSR post-processing
						dspsr_postproc(self, self, cmdline, obs, psr, total_chan, nsubs_eff, self.curdir, self.output_prefix)
						# removing files created by dspsr for each freq channel
						if not cmdline.opts.is_debug:
							remove_list=glob.glob("%s/%s_%s_SB*_CH*.ar" % (self.curdir, psr, self.output_prefix))
							cmd="rm -f %s" % (" ".join(remove_list))
							self.execute(cmd, workdir=self.curdir)

						# adding frequency parts created by digifil together
						# and running rfifind on fil-file if needed
						if cmdline.opts.is_single_pulse:
							try:
								self.log.info("Adding frequency channels together for single-pulse analysis...")
								# we sort the list of fil-files in decrease frequency order (needed for sigproc_splice)
								def filsort(filfile):
									sb=int(filfile.split("_SB")[-1].split("_")[0])
									ch=int(filfile.split("_CH")[-1].split(".fil")[0])
									return sb*self.nrChanPerSub+ch
								fil_files=sorted(glob.glob("%s/%s_%s_SB*_CH*.fil" % (self.curdir, psr, self.output_prefix)), key=filsort, reverse=True)
								cmd="sigproc_splice -a -o %s_%s.fil %s" % (psr, self.output_prefix, " ".join(fil_files))
								self.execute(cmd, workdir=self.curdir)
								# removing files created by digifil for each freq channel
								if not cmdline.opts.is_debug:
									cmd="rm -f %s" % (" ".join(fil_files))
									self.execute(cmd, workdir=self.curdir)
								# running rfifind
								if not cmdline.opts.is_norfi:
									self.log.info("Running rfifind for pulsar %s on a filterbank file...")
									cmd="rfifind -o %s_%s -noclip -blocks 16 %s %s_%s.fil" % (psr, self.output_prefix, cmdline.opts.rfifind_extra_opts, psr, self.output_prefix)
									self.execute(cmd, workdir=self.curdir)
							except: pass

					# removing files created by bf2puma2 for each freq channel
					if not cmdline.opts.is_debug:
						remove_list=glob.glob("%s/%s_SB*_CH*" % (self.curdir, input_prefix))
						cmd="rm -f %s" % (" ".join(remove_list))
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
			if not cmdline.opts.is_plots_only and not cmdline.opts.is_feedback and cmdline.opts.is_single_pulse:
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
			self.log.exception("Oops... 'run_nodal' function for %s%s has crashed!" % (obs.FE and "FE/" or "", self.code))
			self.kill()
			raise

		# kill all open processes
		self.kill()
		# remove reference to PulpLogger class from processing unit
		self.log = None

class FE_CVUnit(CVUnit):
	def __init__(self, obs, cep2, cmdline, tab, log, curstokes, part=0):
		CVUnit.__init__(self, obs, cep2, cmdline, tab, log, curstokes, part)
		# re-assigning procdir from BEAMN to station name
		if obs.FE and self.tab.stationList[0] != "": 
			self.procdir = self.tab.stationList[0]
		# setting outdir and curdir directories
		self.set_outdir(obs, cep2, cmdline)
		# initialiaze feedback unit
		self.feedback_init(obs, cep2, cmdline)
