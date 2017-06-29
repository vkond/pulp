###################################################################
#
# Class PipeUnitPart - aka Processing unit per Beam when there are
# several frequency splits and procesing is done without using
# hoover nodes; and other derivative classes
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
from pulp_pipeunit import PipeUnit, CSUnit, ISUnit, FE_CSUnit, rglob

class PipeUnitPart(PipeUnit):
        def __init__(self, obs, cep2, cmdline, tab, log, curstokes, node, part):
		PipeUnit.__init__(self, obs, cep2, cmdline, tab, log, curstokes, part)
                self.location = node  # processing target node of this bandwidth part

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

	# main run function using dspsr to directly read *.h5 files
	def run_dal(self, obs, cep2, cmdline, log):
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
				if (self.tab.is_coherent and obs.stokesCS == "IQUV") or (self.tab.is_coherent == False and obs.stokesIS == "IQUV"):
					# in this case we get files only for needed Stokes
					infiles=[ff for ff in self.tab.rawfiles[cep2.current_node] if re.search("_S%d_P%03d_" % (self.stokes_index, self.part), ff)]
				else: 
					infiles=[ff for ff in self.tab.rawfiles[cep2.current_node] if re.search("_P%03d_" % self.part, ff)]

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

                                # second, we need to rsync other files for this part from other locus nodes
                                # forming the list of dependent locus nodes that have the rest of the data for this part
				if (self.tab.is_coherent and obs.stokesCS == "IQUV") or (self.tab.is_coherent == False and obs.stokesIS == "IQUV"):
                                	dependent_nodes=[key for key in self.tab.rawfiles.keys() if any("_S%d_P%03d_" % (self.stokes_index, self.part) in ff for ff in self.tab.rawfiles[key])]
				else:
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
						if (self.tab.is_coherent and obs.stokesCS == "IQUV") or (self.tab.is_coherent == False and obs.stokesIS == "IQUV"):
							target_files = [ff for ff in self.tab.rawfiles[node] if re.search("_S%d_P%03d_" % (self.stokes_index, self.part), ff)]
						else:
							target_files = [ff for ff in self.tab.rawfiles[node] if re.search("_P%03d_" % self.part, ff)]
						h5_files = ["%s.h5" % (ff.split(".raw")[0]) for ff in target_files]
						target_files.extend(h5_files)
						self.log.info("Node: %s   Files: %s" % (node, ", ".join(target_files)))
						cmd="rsync %s -rLptgoDxP --rsh='ssh -x' %s:%s ." % (verbose, node, (" %s:" % node).join(target_files))
						self.execute(cmd, workdir=self.curdir)
				else:
                                	for node in dependent_nodes:
                                        	# get the list of necessary *.raw files to copy from this node
						if (self.tab.is_coherent and obs.stokesCS == "IQUV") or (self.tab.is_coherent == False and obs.stokesIS == "IQUV"):
        	                                	target_files = [ff for ff in self.tab.rawfiles[node] if re.search("_S%d_P%03d_" % (self.stokes_index, self.part), ff)]
						else:
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

                                input_files=[ff.split("/")[-1] for ff in glob.glob("%s/*_P%03d_*.h5" % (self.curdir, self.part))]
                                self.log.info("Input data: %s" % ("\n".join(input_files)))

			if (self.tab.is_coherent and obs.stokesCS == "IQUV") or (self.tab.is_coherent == False and obs.stokesIS == "IQUV"):
                        	self.output_prefix="%s_SAP%d_%s_S%d_PART%d" % (obs.id, self.sapid, self.procdir, self.stokes_index, self.part)
			else:
                        	self.output_prefix="%s_SAP%d_%s_PART%d" % (obs.id, self.sapid, self.procdir, self.part)
                        self.log.info("Output file(s) prefix: %s" % (self.output_prefix))

			if not cmdline.opts.is_globalfs:
				target_summary_dir="%s/%s_%s/%s%s/%s/SAP%d/%s" % (self.processed_dir_root, cep2.processed_dir_prefix, self.summary_node, \
        	                        cmdline.opts.outdir == "" and cep2.pipeid or cmdline.opts.outdir, self.summary_node_dir_suffix, self.beams_root_dir, self.sapid, self.procdir)
			else:
				target_summary_dir="%s/%s/%s%s/%s/SAP%d/%s" % (self.processed_dir_root, cep2.processed_dir_prefix, \
        	                        cmdline.opts.outdir == "" and cep2.pipeid or cmdline.opts.outdir, self.summary_node_dir_suffix, self.beams_root_dir, self.sapid, self.procdir)

			proc_subs = self.nrSubsPerFile
                        if self.nrSubsPerFile * (self.part + 1) > self.tab.nrSubbands:
                                proc_subs -= (self.nrSubsPerFile * (self.part + 1) - self.tab.nrSubbands)
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
                                        cmd="python %s/release/share/pulsar/bin/digitize.py %s -s %g -o %s/%s %s" % (cep2.lofarsoft, verbose, cmdline.opts.digitize_sigma, self.curdir, self.raw_8bit_dir, " ".join(["%s/%s" % (self.curdir, ff) for ff in input_files])) 
                                        self.execute(cmd, workdir="%s/%s" % (self.curdir, self.raw_8bit_dir))

                                if not cmdline.opts.is_nodecode:
                                        s0_file=input_files[0]
                                        verbose="-q"
                                        if cmdline.opts.is_debug: verbose="-v"

                                        # running dspsr for every pulsar for this frequency split
                                        self.log.info("Running processing for the frequency split %d ..." % self.part)
                                        # loop on pulsars
                                        for psr in self.psrs:
                                                psr2=re.sub(r'^[BJ]', '', psr)
                                                dspsr_nbins=self.get_best_nbins("%s/%s.par" % (tmpdir, psr2))
						if not cmdline.opts.is_nofold and not cmdline.opts.is_skip_dspsr:
	                                                self.log.info("Running dspsr for pulsar %s..." % (psr))
        	                                        cmd="dspsr -b %d -A -L %d %s -E %s/%s.par -O %s_%s -t %d %s %s" % \
								(dspsr_nbins, cmdline.opts.tsubint, verbose, tmpdir, psr2, psr, self.output_prefix, \
								cmdline.opts.nthreads, cmdline.opts.dspsr_extra_opts, s0_file)
                                	                self.execute(cmd, workdir=self.curdir)
						# run digifil with coherent dedispersion for further single-pulse analysis
						# works only for CV data and we need to run it for every pulsar (different DMs)
						if cmdline.opts.is_single_pulse and self.code=="CV":
							try:
								if not cmdline.opts.is_nofold:
									psrdm=self.get_psr_dm("%s/%s.par" % (tmpdir, psr2))
								else: psrdm=0.0
								self.log.info("Running digifil for pulsar %s..." % (psr))
								cmd="digifil %s -B 512 -b 8 -F %d:D -D %f -o %s_%s.fil %s %s" % \
									(verbose, total_chan > 1 and total_chan or 2 * total_chan, psrdm, psr, self.output_prefix, \
									cmdline.opts.digifil_extra_opts, s0_file)
								self.execute(cmd, workdir=self.curdir)
							# we let PULP to continue if digifil has crashed, as the rest of the pipeline can finish ok
							except: pass

                                        # removing input .raw files (links or rsynced)
                                        if not cmdline.opts.is_debug:
                                                cmd="rm -f %s" % (" ".join(["%s.raw" % (f.split(".h5")[0]) for f in input_files]))
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
				if cmdline.opts.is_single_pulse and self.code=="CV":
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
				if cmdline.opts.is_single_pulse and self.code=="CV":
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


        # main run function to use 2bf2fits for conversion
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

                                self.log.info("Copying .h5 file to the current directory...")
                                if not cmdline.opts.is_globalfs:
                                        if cep2.current_node in self.tab.rawfiles.keys():

                				if (self.tab.is_coherent and obs.stokesCS == "IQUV") or (self.tab.is_coherent == False and obs.stokesIS == "IQUV"):
		        		        	# in this case we get files only for needed Stokes
			                		infiles=[ff for ff in self.tab.rawfiles[cep2.current_node] if re.search("_S%d_P%03d_" % (self.stokes_index, self.part), ff)]
                				else: 
		                    			infiles=[ff for ff in self.tab.rawfiles[cep2.current_node] if re.search("_P%03d_" % self.part, ff)]

                                                for f in infiles:
                                                        # forming the input file string
				                	input_file=f  # should be only one file for a part
                                                        # copy *.h5 files
                                                        cmd="cp -f %s.h5 ." % (f.split(".raw")[0])
				                	try:
	                                                        self.execute(cmd, workdir=self.curdir)
                					except: pass
                                else: # global FS
                                        input_files=[]
                                        for val in self.tab.rawfiles.values(): input_files.extend(val)
                			if (self.tab.is_coherent and obs.stokesCS == "IQUV") or (self.tab.is_coherent == False and obs.stokesIS == "IQUV"):
		        		        # in this case we get files only for needed Stokes
                                                infiles=[ff for ff in input_files if re.search("_S%d_P%03d_" % (self.stokes_index, self.part), ff)]
                                        else:
		                    		infiles=[ff for ff in input_files if re.search("_P%03d_" % self.part, ff)]

                                        for f in infiles:
                                                # forming the input file string
				               	input_file=f  # should be only one file for a part
                                                # copy *.h5 files
                                                cmd="cp -f %s.h5 ." % (f.split(".raw")[0])
				               	try:
                                                        self.execute(cmd, workdir=self.curdir)
                				except: pass

                                self.log.info("Input data: %s" % (input_file))

			if (self.tab.is_coherent and obs.stokesCS == "IQUV") or (self.tab.is_coherent == False and obs.stokesIS == "IQUV"):
                        	self.output_prefix="%s_SAP%d_%s_S%d_PART%d" % (obs.id, self.sapid, self.procdir, self.stokes_index, self.part)
			else:
                        	self.output_prefix="%s_SAP%d_%s_PART%d" % (obs.id, self.sapid, self.procdir, self.part)
                        self.log.info("Output file(s) prefix: %s" % (self.output_prefix))

			if not cmdline.opts.is_globalfs:
				target_summary_dir="%s/%s_%s/%s%s/%s/SAP%d/%s" % (self.processed_dir_root, cep2.processed_dir_prefix, self.summary_node, \
        	                        cmdline.opts.outdir == "" and cep2.pipeid or cmdline.opts.outdir, self.summary_node_dir_suffix, self.beams_root_dir, self.sapid, self.procdir)
			else:
				target_summary_dir="%s/%s/%s%s/%s/SAP%d/%s" % (self.processed_dir_root, cep2.processed_dir_prefix, \
        	                        cmdline.opts.outdir == "" and cep2.pipeid or cmdline.opts.outdir, self.summary_node_dir_suffix, self.beams_root_dir, self.sapid, self.procdir)

			proc_subs = self.nrSubsPerFile
                        if self.nrSubsPerFile * (self.part + 1) > self.tab.nrSubbands:
                                proc_subs -= (self.nrSubsPerFile * (self.part + 1) - self.tab.nrSubbands)
                        nsubs_eff = min(self.tab.nrSubbands, proc_subs)
                        total_chan = nsubs_eff * self.nrChanPerSub

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
                                        cmd="python %s/release/share/pulsar/bin/digitize.py %s -s %g -o %s/%s %s" % (cep2.lofarsoft, verbose, cmdline.opts.digitize_sigma, self.curdir, self.raw_8bit_dir, input_file.replace(".raw", ".h5"))
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
                                                nsamples=self.lcd(nsamples, os.path.getsize(input_file)/(4*nsubs_eff*self.nrChanPerSub))
                                                cmd="2bf2fits %s %s -append -nbits 8 -A %d -sigma %d -nsubs %d -sap %d -tab %d -part %d -stokes %d -o %s \
-nsamples %d -nchans %d -ra %s -dec %s -psr %s -clock %d -band %s -startdate %s -starttime %s -samptime %g -duration %g -subs %s -obsid %s -observer %s %s %s" % \
(verbose, self.raw2fits_extra_options, cmdline.opts.decode_nblocks, cmdline.opts.decode_sigma, nsubs_eff, self.sapid, \
self.tabid, self.part, self.stokes_index, self.output_prefix, nsamples, self.nrChanPerSub, str(self.tab.rarad), \
str(self.tab.decrad), sap.source, obs.sampleClock, obs.bandFilter, obs.startdate, obs.starttime, \
self.sampling/1000., obs.duration, sap.subbandList, obs.id, obs.projectPI.split()[0].split(',')[0], cmdline.opts.bf2fits_extra_opts, input_file)
                                        else: # BG/P call using parset file
                                                cmd="2bf2fits -revend %s %s -parset %s -append -nbits 8 -A %d -sigma %d -nsubs %d -sap %d -tab %d -part %d -stokes %d -duration %g -o %s %s %s" % \
(verbose, self.raw2fits_extra_options, obs.parset, cmdline.opts.decode_nblocks, cmdline.opts.decode_sigma, nsubs_eff, self.sapid, self.tabid, self.part, self.stokes_index, obs.duration, \
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
                                                samples_to_average=int((cmdline.opts.subdyn_time_average * 1000.)/ self.sampling)
                                                cmd="python %s/release/share/pulsar/bin/subdyn.py --psrfits --saveonly -n %d --title \"%s %s Averaging: %g s \" %s.fits" % (cep2.lofarsoft, samples_to_average, obs.id, self.procdir, cmdline.opts.subdyn_time_average, self.output_prefix)
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
                                                for psr in self.psrs: # pulsar list is empty if --nofold is used
                                                        psr2=re.sub(r'^[BJ]', '', psr)
                                                        dspsr_nbins=self.get_best_nbins("%s/%s.par" % (tmpdir, psr2))
                                                        cmd="dspsr -b %d -A -L %d -E %s/%s.par %s -O %s_%s -t %d %s %s.fits" % \
                                                                (dspsr_nbins, cmdline.opts.tsubint, tmpdir, psr2, verbose, psr, \
                                                                self.output_prefix, cmdline.opts.nthreads, cmdline.opts.dspsr_extra_opts, self.output_prefix)
                                                        dspsr_popen = self.start_and_go(cmd, workdir=self.curdir)
                                                        dspsr_popens.append(dspsr_popen)

                                                # waiting for dspsr to finish
                                                self.waiting_list("dspsr", dspsr_popens)

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
                                        if (self.tab.is_coherent and obs.stokesCS == "IQUV") or (self.tab.is_coherent == False and obs.stokesIS == "IQUV"):
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
                                        if (self.tab.is_coherent and obs.stokesCS == "IQUV") or (self.tab.is_coherent == False and obs.stokesIS == "IQUV"):
                                                if self.stokes_index == 0: self.run_rrats_analysis(cmdline, total_chan)
                                        else:
                                        	self.run_rrats_analysis(cmdline, total_chan)
                        except: pass

                        # waiting for subdyn to finish
                        if not cmdline.opts.is_plots_only and not cmdline.opts.is_feedback and \
                        not cmdline.opts.is_norfi and not cmdline.opts.is_skip_subdyn:
                                self.waiting("subdyn.py", subdyn_popen)

			if not cmdline.opts.is_plots_only and not cmdline.opts.is_feedback:

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

                                # rsync'ing all relevant files to the summary node to be combined and further processed
                                # also copy par-file(s) to have it/them in the summary directory
                                verbose=""
                                if cmdline.opts.is_debug: verbose="-v"
				if not cmdline.opts.is_globalfs:
	                                self.log.info("Rsync'ing all relevant files to %s:%s" % (self.summary_node, target_summary_dir))
				copy_list=[]
                        	for ext in self.extensions + ["*.h5", "*.ar", "*.AR", "*.pfd"]:
                                	ext_list=rglob(self.curdir, ext, 3)
					copy_list.extend(ext_list)
                     	   	copy_list.extend(glob.glob("%s/*.par" % (self.outdir)))
				if not cmdline.opts.is_globalfs:
                                	cmd="rsync %s -axP --rsh='ssh -x' %s %s:%s" % (verbose, " ".join(copy_list), self.summary_node, target_summary_dir)
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
                        self.log.exception("Oops... 'run_nodal' function for %s%s has crashed!" % (obs.FE and "FE/" or "", self.code))
                        self.kill()
                        raise

                # kill all open processes
                self.kill()
                # remove reference to PulpLogger class from processing unit
                self.log = None


class CSUnitPart(PipeUnitPart, CSUnit):
        def __init__(self, obs, cep2, cmdline, tab, log, curstokes, node, part):
		PipeUnitPart.__init__(self, obs, cep2, cmdline, tab, log, curstokes, node, part)
                CSUnit.__init__(self, obs, cep2, cmdline, tab, log, curstokes, part)

class ISUnitPart(PipeUnitPart, ISUnit):
        def __init__(self, obs, cep2, cmdline, tab, log, curstokes, node, part):
		PipeUnitPart.__init__(self, obs, cep2, cmdline, tab, log, curstokes, node, part)
                ISUnit.__init__(self, obs, cep2, cmdline, tab, log, curstokes, part)

class FE_CSUnitPart(PipeUnitPart, FE_CSUnit):
        def __init__(self, obs, cep2, cmdline, tab, log, curstokes, node, part):
		PipeUnitPart.__init__(self, obs, cep2, cmdline, tab, log, curstokes, node, part)
                FE_CSUnit.__init__(self, obs, cep2, cmdline, tab, log, curstokes, part)
