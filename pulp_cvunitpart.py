###################################################################
#
# Class CVUnitPart - processsing unit for CV data when there are
# several frequency splits and procesing is done without using
# hoover nodes
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
from pulp_pipeunit import PipeUnit
from pulp_pipeunitpart import PipeUnitPart
from pulp_cvunit import CVUnit

# class for processing of part of BW for CV data
class CVUnitPart(PipeUnitPart, CVUnit):
	def __init__(self, obs, cep2, cmdline, tab, log, curstokes, node, part):
                PipeUnitPart.__init__(self, obs, cep2, cmdline, tab, log, curstokes, node, part)
		CVUnit.__init__(self, obs, cep2, cmdline, tab, log, curstokes, part)

        # main CV processing function
        def run(self, obs, cep2, cmdline, log):
		CVUnit.run(self, obs, cep2, cmdline, log)

	# main run function using bf2puma2 for conversion
	def run_nodal(self, obs, cep2, cmdline, log):
		try:
			self.log = log
			self.logfile = cep2.get_logfile().split("/")[-1]		
			self.start_time=time.time()	

			# start logging
			self.log.info("%s SAP=%d TAB=%d PART=%d %s(%s%s Stokes: %s)    UTC start time is: %s  @node: %s" % \
				(obs.id, self.sapid, self.tabid, self.part, obs.FE and ", ".join(self.tab.stationList) + " " or "", obs.FE and "FE/" or "", self.code, \
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
				# first we make links to *.h5 and *.raw files in the current directory for the target node

				# make a soft links in the current dir (in order for processing to be consistent with the case when data are in many nodes)
				self.log.info("Making links to input files in the current directory...")
				for f in [ff for ff in self.tab.rawfiles[cep2.current_node] if re.search("_P%03d_" % self.part, ff)]:
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

				# second, we need to rsync other files for this part from other locus nodes
				# forming the list of dependent locus nodes that have the rest of the data for this part
				dependent_nodes=[key for key in self.tab.rawfiles.keys() if any("_P%03d_" % self.part in ff for ff in self.tab.rawfiles[key])]
				if cep2.current_node in dependent_nodes: dependent_nodes.remove(cep2.current_node)

				# rsyncing the rest of the data (if we are not using global FS)
				if not cmdline.opts.is_globalfs:
					self.log.info("Rsync'ing other *.h5 and *.raw files from other locus nodes to the current directory...")
			                verbose=""
                			if cmdline.opts.is_debug: verbose="-v"
				else:
					self.log.info("Making links to input files in the current directory for the rest of the files...")

				if not cmdline.opts.is_globalfs:
					for node in dependent_nodes:
						# get the list of necessary *.raw files to copy from this node
						target_files = [ff for ff in self.tab.rawfiles[node] if re.search("_P%03d_" % self.part, ff)]
						h5_files = ["%s.h5" % (ff.split(".raw")[0]) for ff in target_files]
						target_files.extend(h5_files)
						self.log.info("Node: %s   Files: %s" % (node, ", ".join(target_files)))
						cmd="rsync %s -rLptgoDxP --rsh='ssh -x' %s:%s ." % (verbose, node, (" %s:" % node).join(target_files))
						self.execute(cmd, workdir=self.curdir)
				else:
					for node in dependent_nodes:
						# get the list of necessary *.raw files to copy from this node
						target_files = [ff for ff in self.tab.rawfiles[node] if re.search("_P%03d_" % self.part, ff)]
						self.log.info("Node: %s   Files: %s" % (node, ", ".join(target_files)))
						for f in target_files:
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

				input_files=[ff.split("/")[-1] for ff in glob.glob("%s/*_P%03d_*.raw" % (self.curdir, self.part))]
				self.log.info("Input data: %s" % ("\n".join(input_files)))

			self.output_prefix="%s_SAP%d_%s_PART%d" % (obs.id, self.sapid, self.procdir, self.part)
			self.log.info("Output file(s) prefix: %s" % (self.output_prefix))

			if not cmdline.opts.is_globalfs:
				target_summary_dir="%s/%s_%s/%s%s/%s/SAP%d/%s" % (self.processed_dir_root, cep2.processed_dir_prefix, self.summary_node, \
					cmdline.opts.outdir == "" and cep2.pipeid or cmdline.opts.outdir, self.summary_node_dir_suffix, self.beams_root_dir, self.sapid, self.procdir)
			else:
				target_summary_dir="%s/%s/%s%s/%s/SAP%d/%s" % (self.processed_dir_root, cep2.processed_dir_prefix, \
					cmdline.opts.outdir == "" and cep2.pipeid or cmdline.opts.outdir, self.summary_node_dir_suffix, self.beams_root_dir, self.sapid, self.procdir)

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

				if not cmdline.opts.is_nodecode:
					# getting the list of "_S0_" files, the number of which is how many freq splits we have
					# we also sort this list by split number
					s0_file=[f for f in input_files if re.search("_S0_", f) is not None][0]
					# checking if extra options for complex voltages processing were set in the pipeline
					nblocks=""
					if cmdline.opts.nblocks != -1: nblocks = "-b %d" % (cmdline.opts.nblocks)
					is_all_for_scaling=""
					if cmdline.opts.is_all_times: is_all_for_scaling="-all_times"
					is_write_ascii=""
					if cmdline.opts.is_write_ascii: is_write_ascii="-t"
					verbose=""
					if cmdline.opts.is_debug: verbose="-v"

					# running bf2puma2 command for this frequency split
					self.log.info("Running bf2puma2 for the frequency split %d ..." % self.part)
					cmd="bf2puma2 -f %s -h %s -p %s -hist_cutoff %f %s %s %s %s" % (s0_file, cep2.puma2header, \
						obs.parset, cmdline.opts.hist_cutoff, verbose, nblocks, is_all_for_scaling, is_write_ascii)
					self.execute(cmd, workdir=self.curdir)

				# removing input .raw files (links or rsynced)
				if not cmdline.opts.is_debug:
					cmd="rm -f %s" % (" ".join(input_files))
					self.execute(cmd, workdir=self.curdir)

				if not cmdline.opts.is_nofold and not cmdline.opts.is_skip_dspsr:
					# getting the MJD of the observation
					# for some reason dspsr gets wrong MJD without using -m option
					input_prefix=s0_file.split("_S0_")[0]
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
					self.log.info("Running processing for each frequency channel for the frequency split %d ..." % self.part)
					# loop on pulsars
					for psr in self.psrs:
						psr2=re.sub(r'^[BJ]', '', psr)
						psrdm=self.get_psr_dm("%s/%s.par" % (tmpdir, psr2))
						dspsr_nbins=self.get_best_nbins("%s/%s.par" % (tmpdir, psr2))
						self.log.info("Running dspsr for pulsar %s..." % (psr))
						for bb in range(0, len(bf2puma_outfiles), self.nrChanPerSub):
							self.log.info("For %s" % (", ".join(bf2puma_outfiles[bb:bb+self.nrChanPerSub])))
							dspsr_popens=[] # list of dspsr Popen objects
							for cc in range(bb, bb+self.nrChanPerSub):
								input_file=bf2puma_outfiles[cc]
								cmd="dspsr -m %s -b %d -A -L %d %s -E %s/%s.par -O %s_%s_SB%s -t %d %s %s" % \
									(obsmjd, dspsr_nbins, cmdline.opts.tsubint, verbose, tmpdir, psr2, psr, self.output_prefix, \
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

						# running psradd to add all freq channels together for this frequency split
						self.log.info("Adding frequency channels together for the frequency split %d ..." % self.part)
						ar_files=glob.glob("%s/%s_%s_SB*_CH*.ar" % (self.curdir, psr, self.output_prefix))
						cmd="psradd -R -m time -o %s_%s.ar %s" % (psr, self.output_prefix, " ".join(ar_files))
						self.execute(cmd, workdir=self.curdir)
						
						# removing files created by dspsr for each freq channel
						if not cmdline.opts.is_debug:
							remove_list=glob.glob("%s/%s_%s_SB*_CH*.ar" % (self.curdir, psr, self.output_prefix))
							cmd="rm -f %s" % (" ".join(remove_list))
							self.execute(cmd, workdir=self.curdir)

						# adding frequency parts created by digifil together
						if cmdline.opts.is_single_pulse:
							try:
								self.log.info("Adding frequency channels together for single-pulse analysis for the frequency split %d..." % self.part)
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
							except: pass

					# removing files created by bf2puma2 for each freq channel
					if not cmdline.opts.is_debug:
						remove_list=glob.glob("%s/%s_SB*_CH*" % (self.curdir, input_prefix))
						cmd="rm -f %s" % (" ".join(remove_list))
						self.execute(cmd, workdir=self.curdir)

				# creating beam directory on summary node (where to copy ar-files)
				self.log.info("Creating directory %s on summary node %s..." % (target_summary_dir, self.summary_node))
				if cmdline.opts.is_slurm:
					docker_cmd_prefix=docker_cmd_suffix=""
					if cmdline.opts.is_docker:
						docker_cmd_prefix=cep2.docker_cmd_prefix
						docker_cmd_suffix="%s -e SLURM_JOB_ID=%s -e SLURM_JOB_NAME=%s %s" % \
							(cep2.docker_common_opts, cep2.slurm_jobid, cep2.slurm_jobname, cmdline.opts.docker_container)
					cmd="%ssrun %s -c 1 --jobid=%s --job-name=%s %s -w %s %s mkdir -m 775 -p %s" % \
						(docker_cmd_prefix, cep2.srun_general_opts, cep2.slurm_jobid, cep2.slurm_jobname, cep2.slurm_extra_opts, self.summary_node, docker_cmd_suffix, target_summary_dir)
				else:
					cmd="ssh -t %s 'mkdir -m 775 -p %s'" % (self.summary_node, target_summary_dir)
				p = Popen(shlex.split(cmd), stdout=PIPE, stderr=STDOUT)
				p.communicate()

				# rsync'ing the ar-files and h5-files to the summary node to be combined and further processed
				# also copy par-file(s) to have it/them in the summary directory
				# also copy .fil files if CV obs and is_single_pulse==True
				verbose=""
				if cmdline.opts.is_debug: verbose="-v"
				if cmdline.opts.is_single_pulse:
					if not cmdline.opts.is_globalfs:
						self.log.info("Rsync'ing *.ar, *.fil files and *.h5 files to %s:%s" % (self.summary_node, target_summary_dir))
					fil_list=[]
					for psr in self.psrs:
						fil_list.append("%s/%s_%s.fil" % (self.curdir, psr, self.output_prefix))
				else:
					if not cmdline.opts.is_globalfs:
						self.log.info("Rsync'ing *.ar files and *.h5 files to %s:%s" % (self.summary_node, target_summary_dir))
				ar_list=[]
				for psr in self.psrs:
					ar_list.append("%s/%s_%s.ar" % (self.curdir, psr, self.output_prefix))
				# making list of h5 files
				h5_list = glob.glob("%s/*.h5" % (self.curdir))
				# also making list of par-files
				par_list = glob.glob("%s/*.par" % (self.outdir))
				if cmdline.opts.is_single_pulse:
					if not cmdline.opts.is_globalfs:
						cmd="rsync %s -axP --rsh='ssh -x' %s %s %s %s %s:%s" % \
							(verbose, " ".join(ar_list), " ".join(fil_list), " ".join(h5_list), " ".join(par_list), self.summary_node, target_summary_dir)
						self.execute(cmd, workdir=self.curdir)
				else:
					if not cmdline.opts.is_globalfs:
						cmd="rsync %s -axP --rsh='ssh -x' %s %s %s %s:%s" % (verbose, " ".join(ar_list), " ".join(h5_list), " ".join(par_list), self.summary_node, target_summary_dir)
						self.execute(cmd, workdir=self.curdir)

                        # finishing off the processing...
                        self.finish_off(obs, cep2, cmdline, False, target_summary_dir)

                        # flushing feedback file to disk (success)
                        self.fbunit.flush(100, cep2, False)

                        # removing tmpdir
                        if not cmdline.opts.is_nofold:
                                cmd="rm -rf %s" % (tmpdir)
                                os.system(cmd)

		except Exception:
			self.log.exception("Oops... 'run_dal' function for %s%s has crashed!" % (obs.FE and "FE/" or "", self.code))
			self.kill()
			raise

		# kill all open processes
		self.kill()
		# remove reference to PulpLogger class from processing unit
		self.log = None

# class for processing of part of BW for CV data in case of FE observation
class FE_CVUnitPart(CVUnitPart):
	def __init__(self, obs, cep2, cmdline, tab, log, curstokes, node, part):
                CVUnitPart.__init__(self, obs, cep2, cmdline, tab, log, curstokes, node, part)
                # re-assigning procdir from BEAMN to station name
                if obs.FE and self.tab.stationList[0] != "":
                        self.procdir = self.tab.stationList[0]
                # setting outdir and curdir directories
                self.set_outdir(obs, cep2, cmdline)
		# initialiaze feedback unit
		self.feedback_init(obs, cep2, cmdline)
