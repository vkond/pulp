###################################################################
#
# Class Pipeline - main processing class
# Other class, aka Processing Units (per Beam) are defined
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
from pulp_pipeunit import *
from pulp_cvunit import *
from pulp_pipeunitpart import *
from pulp_cvunitpart import *

# The main processing class
class Pipeline:
	def __init__(self, obs, cep2, cmdline, log):
		self.units = []   # list of processing units
		self.feedbacks = [] # list of feedback units
		self.sum_popens = []  # list of Popen make_summary processes
		self.summary_dirs = {}  # dictionary, key - "CS", "IS", or "CV" (same as for summary nodes dict in CEP2Info class, value - summary dir
		# extensions of the files to copy to archive in summary (*_nopfd*.tgz)
		#self.summary_archive_exts=["*.log", "*.txt", "*.pdf", "*.ps", "*.bestprof", "*.inf", "*.rfirep", "*png", "*parset", "*.par", "*_sp_inf.tar.gz", "*_sp_singlepulse.tar.gz"]
		self.summary_archive_exts=["*.log", "*.txt", "*.pdf", "*.ps", "*.pfd", "*.bestprof", "*.polycos", "*.inf", "*.rfirep", "*png", "*.ar", "*.AR", "*pdmp*", "*_rfifind*", "*.dat", "*.singlepulse", "*.h5", "*.fil", "*.rv", "*.out", "*parset", "*.par", "*_sp_inf.tar.gz", "*_sp_singlepulse.tar.gz", "*.per", "*.posn", "*.out", "*.txt"]
		# extensions for summaries which do not require combining parts together
		self.summary_archive_exts_nocombine=["*.log", "*.txt", "*.pdf", "*.ps", "*.bestprof", "*.rfirep", "*png", "*.rv", "*.out", "*parset", "*.par", "*_sp_inf.tar.gz", "*_sp_singlepulse.tar.gz", "*.per", "*.posn", "*.h5", "*.out", "*.txt"]

		# prefix and suffix for summary archive name, in between them there will CS, IS, CV code
		self.summary_archive_prefix="_summary"
		self.summary_archive_suffix=".tar"
		self.number_failed_pipes = 0   # number of failed pipelines
		self.number_failed_summaries = 0  # number of failed summaries
		self.log = None   # to be initialized later when doing summaries (in order to execute, etc. functions to be similar in syntax with those in PipeUnit)

		already_assigned_nodes=[] # list of nodes that already assigned as host nodes for processing of one of the parts
		# initializing Processing Units based on the list of beams to process
                for beam in cmdline.beams:
			tab = 0
                        sapid=int(beam.split(":")[0])
                        tabid=int(beam.split(":")[1])
			for ss in obs.saps:
				if ss.sapid == sapid:
					for tt in ss.tabs:
						if tt.tabid == tabid:
							tab = tt
							break

			if tab == 0: # in case no TABs were found (e.g. data is missing, etc)
				log.error("No info or data available for the beam %d:%d!" % (sapid, tabid))
				raise Exception

			# here we pseudo-randomly choose "active" stokes value (0, 1, 2, or 3) for CV data
			# Active in a sense, that other stokes raw files (for the given part) will be copied to this node. If active stokes is always 0
			# then there are significantly increase in data volume for lower nodes on CEP2
			# This value should depend on ObsID (or PipeID), so it is the same when pipeline is needed to be rerun.
			# It also depends on beam number, and is different for different beams
			active_stokes=((int(cep2.get_pipeid()[1:])%16)/4 + sapid*len(obs.saps) + tabid)%4
			

			# checking if location (locus node) of raw data or at least processed data is known
			# If not, then skip this beam
			if len(tab.location) == 0: continue

			is_iquv_not_split = False
			if (tab.is_coherent and obs.stokesCS == "IQUV" and obs.nsplitsCS == 1) or (tab.is_coherent == False and obs.stokesIS == "IQUV" and obs.nsplitsIS == 1):
				is_iquv_not_split = True

			# if we do not want to use hoover nodes than we should make first a list of nodes
			# where processing will happen for parts of the band
			if (len(tab.location) > 1 or cmdline.opts.is_all_parts_at_once) and cmdline.opts.is_nohoover and not is_iquv_not_split:
				for curstokes in np.arange(4): # loop on stokes files, S0-S3
					if tab.is_coherent and obs.CV and curstokes != active_stokes: continue  # for CV data we process all stokes files alltogether
					if tab.is_coherent and obs.CS and obs.stokesCS == "I" and curstokes > 0: continue # for CS StokesI there are only S0 files
					if tab.is_coherent == False and obs.IS and obs.stokesIS == "I" and curstokes > 0: continue # for IS StokesI there are only S0 files
					if cmdline.opts.stokes != -1 and cmdline.opts.stokes != curstokes: continue  # if we want to process particular Stokes, then we skip others

					if tab.is_coherent:
						ffs = cmdline.opts.first_freq_splitCS
						nsplits = cmdline.opts.nsplitsCS
					else:
						ffs = cmdline.opts.first_freq_splitIS
						nsplits = cmdline.opts.nsplitsIS
					for part in np.arange(ffs, ffs + nsplits):
						part_nodes=[key for key in tab.rawfiles.keys() if any("_S%d_P%03d_" % (curstokes, part) in ff for ff in tab.rawfiles[key])]
						if len(part_nodes) != 0:
							node=part_nodes[0]
							if node in already_assigned_nodes:
								# i.e. this node is already taken to process some other part, so we have to check if one of
								# other nodes that have this part is not taken yet
								# getting the list of nodes from "part_nodes" that are not in the "already_assigned_nodes"
								good_nodes=list(set(part_nodes)-set(part_nodes).intersection(set(already_assigned_nodes)))
								if len(good_nodes) != 0: # means that there is separate node available
									node=good_nodes[0]
									already_assigned_nodes.append(node)
								else: # in case of CV data we still have to check S1-S3 locus nodes
									if tab.is_coherent and obs.CV:
										for cvstokes in list(set(np.arange(4))-set(np.arange(4)).intersection(set([active_stokes]))):
											part_nodes1=[key for key in tab.rawfiles.keys() if any("_S%d_P%03d_" % (cvstokes, part) in ff for ff in tab.rawfiles[key])]
											if len(part_nodes1) != 0:
												node=part_nodes1[0]
												if node in already_assigned_nodes:
													good_nodes=list(set(part_nodes1)-set(part_nodes1).intersection(set(already_assigned_nodes)))
													if len(good_nodes) != 0: # means that there is separate node available
														node=good_nodes[0]
														already_assigned_nodes.append(node)
														break	
												else:
													already_assigned_nodes.append(node)
													break
											else: # this is possible if we want to run only summaries, but one of the locus nodes with the data is not accessible
												if not cmdline.opts.is_summary:
													log.error("No locus node available to process the data for Part=%d for SAP=%d TAB=%d." % (part, sapid, tabid))
													raise Exception
									
								# if there are no good nodes we just process the data on the same node
							else:
								already_assigned_nodes.append(node)
	
#								log.error("No unique separate node available to process the Part=%d for SAP=%d TAB=%d.\nTry process without --no-hoover option!" % \
#									(part, sapid, tabid))
#								raise Exception
						else: # this is possible if we want to run only summaries, but one of the locus nodes with the data is not accessible
							if not cmdline.opts.is_summary:
								log.error("No locus node available to process the data for Part=%d for SAP=%d TAB=%d." % (part, sapid, tabid))
								raise Exception

						if not tab.is_coherent:
							unit = ISUnitPart(obs, cep2, cmdline, tab, log, curstokes, node, part)
						if tab.is_coherent and tab.specificationType != "flyseye":
							if obs.CS:
								unit = CSUnitPart(obs, cep2, cmdline, tab, log, curstokes, node, part)
							elif obs.CV:
								unit = CVUnitPart(obs, cep2, cmdline, tab, log, curstokes, node, part)
							else:
								log.error("Can't initialize processing pipeline unit for SAP=%d TAB=%d PART=%d STOKES=%d on node %s" % (sapid, tabid, part, curstokes, node))
								raise Exception
						if tab.is_coherent and tab.specificationType == "flyseye":
							if obs.stokesCS[0] == "I":  # for I and IQUV
								unit = FE_CSUnitPart(obs, cep2, cmdline, tab, log, curstokes, node, part)
							elif obs.stokesCS[0] == "X": # for XY (old format) or XXYY
								unit = FE_CVUnitPart(obs, cep2, cmdline, tab, log, curstokes, node, part)
							else:
								log.error("Can't initialize processing pipeline FE unit for SAP=%d TAB=%d PART=%d STOKES=%d on node %s" % (sapid, tabid, part, curstokes, node))
								raise Exception

						# adding unit to the list
						self.units.append(unit)

			else: # only one locus node in the beam
				for curstokes in np.arange(4): # loop on stokes files, S0-S3
					if tab.is_coherent and obs.CV and curstokes != active_stokes: continue  # for CV data we process all stokes files alltogether
					if tab.is_coherent and obs.CS and obs.stokesCS == "I" and curstokes > 0: continue # for CS StokesI there are only S0 files
					if tab.is_coherent == False and obs.IS and obs.stokesIS == "I" and curstokes > 0: continue # for IS StokesI there are only S0 files
					if cmdline.opts.stokes != -1 and cmdline.opts.stokes != curstokes: continue  # if we want to process particular Stokes, then we skip others

					if not tab.is_coherent:
						unit = ISUnit(obs, cep2, cmdline, tab, log, curstokes)
					if tab.is_coherent and tab.specificationType != "flyseye":
						if obs.CS:
							unit = CSUnit(obs, cep2, cmdline, tab, log, curstokes)
						elif obs.CV:
							unit = CVUnit(obs, cep2, cmdline, tab, log, curstokes)
						else:
							log.error("Can't initialize processing pipeline unit for SAP=%d TAB=%d STOKES=%d" % (sapid, tabid, curstokes))
							raise Exception
					if tab.is_coherent and tab.specificationType == "flyseye":
						if obs.stokesCS[0] == "I":  # for I and IQUV
							unit = FE_CSUnit(obs, cep2, cmdline, tab, log, curstokes)
						elif obs.stokesCS[0] == "X": # for XY (old format) or XXYY
							unit = FE_CVUnit(obs, cep2, cmdline, tab, log, curstokes)
						else:
							log.error("Can't initialize processing pipeline FE unit for SAP=%d TAB=%d STOKES=%d" % (sapid, tabid, curstokes))
							raise Exception

					# adding unit to the list
					self.units.append(unit)

		if len(self.units) == 0:
			log.info("None beams to process!")
			raise Exception

		# creating main output directory on locus092 for CS data and on locus094 for IS data, and on locus093 for CV data
		# before that we also remove this directory if user flag is_del was set
		if cep2.cluster_headnode[:5] == "head0" or cep2.cluster_headnode[:3] == "cpu": # CEP4
			unique_outdirs=["%s:%s:%s/%s/%s%s" % (unit.code, unit.summary_node, unit.processed_dir_root, cep2.processed_dir_prefix, \
				cmdline.opts.outdir == "" and cep2.pipeid or cmdline.opts.outdir, \
				unit.summary_node_dir_suffix) for unit in self.units if unit.summary_node != "" and unit.summary_node_dir_suffix != ""]
		else:
			unique_outdirs=["%s:%s:%s/%s_%s/%s%s" % (unit.code, unit.summary_node, unit.processed_dir_root, cep2.processed_dir_prefix, unit.summary_node, \
				cmdline.opts.outdir == "" and cep2.pipeid or cmdline.opts.outdir, \
				unit.summary_node_dir_suffix) for unit in self.units if unit.summary_node != "" and unit.summary_node_dir_suffix != ""]
		unique_outdirs=np.unique(unique_outdirs)
		fbindex = self.get_number_fbunits(obs, cmdline) # file index for feedback file
		for uo in unique_outdirs:
			data_code=uo.split(":")[0]
			node=uo.split(":")[1]
			sumdir=uo.split(":")[-1]
			self.summary_dirs[data_code] = sumdir
			fbunit = FeedbackUnit(fbindex, node, sumdir)
			fbunit.update("%s/" % (sumdir), "%s/.%s%s%s.fb" % (cep2.get_logdir(), cep2.pipeid, self.summary_archive_prefix, data_code), data_code, obs, True)
			fbunit.flush(0, cep2, True)
			self.feedbacks.append(fbunit)
			fbindex += 1
			# deleting previous results if option --del was set
			if cmdline.opts.is_delete:
				log.info("Deleting previous summary results on %s: %s" % (node, sumdir))
				if cmdline.opts.is_slurm:
					# for this simple command we do not need to use docker image
					slurm_summaries_extra_opts=cep2.slurm_summaries_extra_opts
					if not cmdline.opts.is_globalfs:
						slurm_summaries_extra_opts="%s -w %s" % (cep2.slurm_summaries_extra_opts, node)
					docker_cmd_prefix=docker_cmd_suffix=""
					if cmdline.opts.is_docker:
						docker_cmd_prefix=cep2.docker_cmd_prefix
						docker_cmd_suffix="%s -e SLURM_JOB_ID=%s -e SLURM_JOB_NAME=%s %s" % \
							(cep2.docker_common_opts, cep2.slurm_jobid, cep2.slurm_jobname, cmdline.opts.docker_container)
					cmd="%ssrun %s -c 1 --jobid=%s --job-name=%s %s %s rm -rf %s" % \
						(docker_cmd_prefix, cep2.srun_general_opts, cep2.slurm_jobid, cep2.slurm_jobname, slurm_summaries_extra_opts, docker_cmd_suffix, sumdir)
				else:
					# for this simple command we do not need to use docker image
					if not cmdline.opts.is_globalfs:
						cmd="ssh -t %s 'rm -rf %s'" % (node, sumdir)
					else:
						cmd="rm -rf %s" % (sumdir)
				p = Popen(shlex.split(cmd), stdout=PIPE, stderr=STDOUT)
				p.communicate()
			log.info("Creating output summary directory on %s: %s" % (node, sumdir))
			if cmdline.opts.is_slurm:
				# for this simple command we do not need to use docker image
				slurm_summaries_extra_opts=cep2.slurm_summaries_extra_opts
				if not cmdline.opts.is_globalfs:
					slurm_summaries_extra_opts="%s -w %s" % (cep2.slurm_summaries_extra_opts, node)
				docker_cmd_prefix=docker_cmd_suffix=""
				if cmdline.opts.is_docker:
					docker_cmd_prefix=cep2.docker_cmd_prefix
					docker_cmd_suffix="%s -e SLURM_JOB_ID=%s -e SLURM_JOB_NAME=%s %s" % \
						(cep2.docker_common_opts, cep2.slurm_jobid, cep2.slurm_jobname, cmdline.opts.docker_container)
				cmd="%ssrun %s -c 1 --jobid=%s --job-name=%s %s %s mkdir -m 775 -p %s" % \
					(docker_cmd_prefix, cep2.srun_general_opts, cep2.slurm_jobid, cep2.slurm_jobname, slurm_summaries_extra_opts, docker_cmd_suffix, sumdir)
			else:
				# for this simple command we do not need to use docker image
				if not cmdline.opts.is_globalfs:
					cmd="ssh -t %s 'mkdir -m 775 -p %s'" % (node, sumdir)
				else:
					cmd="mkdir -m 775 -p %s" % (sumdir)
			p = Popen(shlex.split(cmd), stdout=PIPE, stderr=STDOUT)
			p.communicate()
		
		# Defining output directories for all local locus nodes
		# We have to do it here rather than in 'run' function of PipeUnit because there can be more than 1 beam per node
		# and then there can be a clash when we will delete (if -del is set) the output directory 
		#unit_outdirs = ["%s:%s" % (unit.outdir.split(cep2.processed_dir_prefix + "_")[1].split("/")[0], unit.outdir) for unit in self.units]
		unit_outdirs = ["%s:%s" % (unit.location, unit.outdir) for unit in self.units]
		unit_outdirs=np.unique(unit_outdirs)

		# Deleting all local output directories if --del is set 
		if cmdline.opts.is_delete:
			for uo in unit_outdirs:
				node=uo.split(":")[0]
				outdir=uo.split(":")[1]
				log.info("Deleting previous processed results on %s: %s" % (node, outdir))
				if cmdline.opts.is_slurm:
					# for this simple command we do not need to use docker image
					slurm_extra_opts=cep2.slurm_extra_opts
					if not cmdline.opts.is_globalfs:
						slurm_extra_opts="%s -w %s" % (cep2.slurm_extra_opts, node)
					docker_cmd_prefix=docker_cmd_suffix=""
					if cmdline.opts.is_docker:
						docker_cmd_prefix=cep2.docker_cmd_prefix
						docker_cmd_suffix="%s -e SLURM_JOB_ID=%s -e SLURM_JOB_NAME=%s %s" % \
							(cep2.docker_common_opts, cep2.slurm_jobid, cep2.slurm_jobname, cmdline.opts.docker_container)
					cmd="%ssrun %s -c 1 --jobid=%s --job-name=%s %s %s rm -rf %s" % \
						(docker_cmd_prefix, cep2.srun_general_opts, cep2.slurm_jobid, cep2.slurm_jobname, slurm_extra_opts, docker_cmd_suffix, outdir)
				else:
					# for this simple command we do not need to use docker image
					if not cmdline.opts.is_globalfs:
						cmd="ssh -t %s 'rm -rf %s'" % (node, outdir)
					else:
						cmd="rm -rf %s" % (outdir)
				p = Popen(shlex.split(cmd), stdout=PIPE, stderr=STDOUT)
				p.communicate()
				Popen(shlex.split("stty sane"), stderr=open(os.devnull, 'rb')).wait()
		log.info("")


	# get total number of feedback units (so to know what index to assign for summary feedbacks)
        # function to initialize feedback unit (should be called after self.outdir is set)
        def get_number_fbunits(self, obs, cmdline):
                nunits=0
		for sap in obs.saps:
			for tab in sap.tabs:
				if tab.is_coherent:
					if obs.stokesCS == "IQUV": nunits += 4 * cmdline.opts.nsplitsCS
					else: nunits += cmdline.opts.nsplitsCS
				else:
					if obs.stokesIS == "IQUV": nunits += 4 * cmdline.opts.nsplitsIS
					else: nunits += cmdline.opts.nsplitsIS
		return nunits

	# kicking off the pipeline
	def start(self, obs, cep2, cmdline, log):
		# and here we start...
		log.info("Starting PULP processing for:")

		for unit in self.units:
			# if data for this beam are on several nodes, then we have to log in to hoover node...
			tabpart=""
			locations="[#locations = %d, #files = %d]" % (len(unit.tab.location), unit.tab.numfiles)

			is_iquv_not_split = False
			if (unit.tab.is_coherent and obs.stokesCS == "IQUV") or (unit.tab.is_coherent == False and obs.stokesIS == "IQUV"):
				stoki="STOKES=%s " % (unit.stokes)
				if (unit.tab.is_coherent and obs.nsplitsCS == 1) or (unit.tab.is_coherent == False and obs.nsplitsIS == 1):
					is_iquv_not_split = True
					target_nodes=[key for key in unit.tab.rawfiles.keys() if any("_S%d_" % (unit.stokes_index) in ff for ff in unit.tab.rawfiles[key])]
					nfiles=int(np.sum([len([val for val in unit.tab.rawfiles[key] if re.search("_S%d_" % (unit.stokes_index), val)]) for key in target_nodes]))
					locations="[#locations = %d, #files = %d]" % (len(target_nodes), nfiles)
					locus=target_nodes[0]
			else: stoki=""

			if (len(unit.tab.location) > 1 or cmdline.opts.is_all_parts_at_once) and not is_iquv_not_split:
				if cmdline.opts.is_nohoover: # are not using hoover nodes
					locus=unit.location
					tabpart="PART=%d " % (unit.part)
					if (unit.tab.is_coherent and obs.stokesCS == "IQUV") or (unit.tab.is_coherent == False and obs.stokesIS == "IQUV"):
						target_nodes=[key for key in unit.tab.rawfiles.keys() if any("_S%d_P%03d_" % (unit.stokes_index, unit.part) in ff for ff in unit.tab.rawfiles[key])]
						nfiles=int(np.sum([len([val for val in unit.tab.rawfiles[key] if re.search("_S%d_P%03d_" % (unit.stokes_index, unit.part), val)]) for key in target_nodes]))
					else:
						target_nodes=[key for key in unit.tab.rawfiles.keys() if any("_P%03d_" % unit.part in ff for ff in unit.tab.rawfiles[key])]
						nfiles=int(np.sum([len([val for val in unit.tab.rawfiles[key] if re.search("_P%03d_" % unit.part, val)]) for key in target_nodes]))
					locations="[#locations = %d, #files = %d]" % (len(target_nodes), nfiles)
				else: locus=cep2.hoover_nodes[0]
			else:
				if not is_iquv_not_split:
					locus=unit.tab.location[0]
					
			use_pulp_switch=""
			if cmdline.opts.is_auto and cep2.cluster_headnode[:5] == "lhn00":
				use_pulp_switch="use Pulp; "

			if cmdline.opts.is_slurm:
				slurm_extra_opts=cep2.slurm_extra_opts
				if not cmdline.opts.is_globalfs:
					slurm_extra_opts="%s -w %s" % (cep2.slurm_extra_opts, locus)
				docker_cmd_prefix=docker_cmd_suffix=""
				if cmdline.opts.is_docker:
					docker_cmd_prefix=cep2.docker_cmd_prefix
					docker_cmd_suffix="%s -e SLURM_JOB_ID=%s -e SLURM_JOB_NAME=%s %s" % \
						(cep2.docker_common_opts, cep2.slurm_jobid, cep2.slurm_jobname, cmdline.opts.docker_container)
				cmd="%ssrun %s -c %d --jobid=%s --job-name=%s %s %s %s'%spulp.py %s --noinit --local %s --beams %d:%d%s --confdir %s %s'%s" % \
					(docker_cmd_prefix, cep2.srun_general_opts, cmdline.opts.nthreads, cep2.slurm_jobid, cep2.slurm_jobname, slurm_extra_opts, \
					docker_cmd_suffix, cep2.start_shell, use_pulp_switch, cmdline.opts.is_auto and "--id %s" % obs.id or "", \
					stoki != "" and "--stokes %d" % unit.stokes_index or "", unit.sapid, unit.tabid, tabpart != "" and "/%d" % unit.part or "", \
					cmdline.opts.confdir, " ".join(cmdline.options), cep2.end_shell)
				log.info("CMD = %s" % (cmd))
			else:
				docker_cmd_suffix=""
				if cmdline.opts.is_docker:
					docker_cmd_suffix="%s -e SLURM_JOB_ID=%s -e SLURM_JOB_NAME=%s %s" % \
						(cep2.docker_common_opts, cep2.slurm_jobid, cep2.slurm_jobname, cmdline.opts.docker_container)
				cmd="ssh -t %s %s '%spulp.py %s --noinit --local %s --beams %d:%d%s --confdir %s %s'" % \
					(locus, docker_cmd_suffix, use_pulp_switch, cmdline.opts.is_auto and "--id %s" % obs.id or "", \
					stoki != "" and "--stokes %d" % unit.stokes_index or "", unit.sapid, unit.tabid, tabpart != "" and "/%d" % unit.part or "", \
					cmdline.opts.confdir, " ".join(cmdline.options))
			if cmdline.opts.is_auto or cmdline.opts.is_debug:
				log.info("   Starting processing on locus node: %s..." % (locus))
				log.info("     CMD: %s" % (cmd))
			unit.parent = Popen(shlex.split(cmd), stdout=PIPE, stderr=STDOUT)
			Popen(shlex.split("stty sane"), stderr=open(os.devnull, 'rb')).wait()
			log.info("SAP=%d TAB=%d %s%s%s(%s%s) on %s (pid=%d)  %s" % \
				(unit.sapid, unit.tabid, stoki, tabpart, unit.tab.specificationType == "flyseye" and ", ".join(unit.tab.stationList) + " " or "", \
				unit.tab.specificationType == "flyseye" and "FE/" or "", \
				unit.tab.is_coherent and (obs.CS and "CS" or "CV") or "IS", \
				locus, unit.parent.pid, locations))
			time.sleep(1) # wait 1 sec (otherwise terminal output gets messed up often)


	# here we finish the pipeline, creating extra files, convert ps(pdf) to png, etc...
	# create FE status maps, TA heat maps...
	def finish(self, obs, cep2, cmdline, log):

		try:
			# waiting for processing in individual locus nodes to finish
			# unless we want to do just a summary
			if not cmdline.opts.is_summary:
				run_units = [u.parent.pid for u in self.units if u.parent.poll() is None]
				Popen(shlex.split("stty sane"), stderr=open(os.devnull, 'rb')).wait()
				log.info("Still running [%d]: %s" % (len(run_units), run_units))
				for unit in self.units:
					Popen(shlex.split("stty sane"), stderr=open(os.devnull, 'rb')).wait()
					log.info("waiting...")
					unit.parent.communicate()
					log.info("Process pid=%d has finished, status=%d" % (unit.parent.pid, unit.parent.returncode))
					run_units = [u.parent.pid for u in self.units if u.parent.poll() is None]
					if len(run_units) > 0: log.info("Still running [%d]: %s" % (len(run_units), run_units))

				# loop over finished processes to see if they all finished OK
				failed_units = [u for u in self.units if u.parent.returncode > 0]
				Popen(shlex.split("stty sane"), stderr=open(os.devnull, 'rb')).wait()
				self.number_failed_pipes = len(failed_units)
				log.info("Failed beams [%d]: %s" % (self.number_failed_pipes, ",".join(["%s:%s%s" % (u.sapid, u.tabid, \
					((len(u.tab.location)>1 or cmdline.opts.is_all_parts_at_once) and cmdline.opts.is_nohoover) and " (P%d)" % u.part or "") for u in failed_units])))
				if self.number_failed_pipes > 0:
					log.info("*** Summaries will not be complete! Re-run processing for the failed beams using --beams option. ***")

			if not cmdline.opts.is_nosummary:
				self.sum_popens=[]
				log.info("Starting summaries...")

				use_pulp_switch=""
				if cmdline.opts.is_auto and cep2.cluster_headnode[:5] == "lhn00":
					use_pulp_switch="use Pulp; "


				# when using Slurm we should shrink the job allocation to use only nodes necessary for summaries
				if cmdline.opts.is_slurm:
					# adding 1 more because the main Pulp process is also run under Slurm
					#need_only_nnodes=len(self.summary_dirs.items()) + 1
					need_only_nnodes=2 # maximum we will need 3 job steps (cpus), and original allocation should have 20 tasks per node,
							   # so we should just fit in 1 node
					log.info("Shrinking Slurm job %s (%s) allocation to %d node(s) to finish up with summaries..." % (cep2.slurm_jobid, cep2.slurm_jobname, need_only_nnodes))
					try:
						docker_cmd_prefix=cep2.docker_cmd_prefix
						cmd="%sscontrol -v update jobid=%s NumNodes=%d" % (docker_cmd_prefix, cep2.slurm_jobid, need_only_nnodes)
						##cmd="scontrol update jobid=%s NumTasks=%d NumCPUs=%d" % (cep2.slurm_jobid, need_only_nnodes, need_only_nnodes)
						##cmd="scontrol update jobid=%s NumTasks=%d" % (cep2.slurm_jobid, need_only_nnodes)
						log.info(cmd)
						os.system(cmd)
						#p = Popen(shlex.split(cmd), stdout=PIPE, stderr=STDOUT)
						#p = Popen(shlex.split(cmd), stdout=PIPE, stderr=PIPE, shell=True)
						#(sout, serr) = p.communicate()
						#if sout != NULL and sout != "": log.info(sout)	
						#if serr != NULL and serr != "": log.info(serr)	
					except:
						log.warning("scontrol failed")

				# starting separate pulp.py on summary nodes just to finish up
				for (sumcode, sumdir) in self.summary_dirs.items():
					sumnode = cep2.summary_nodes[sumcode]
					if cmdline.opts.is_slurm:
						slurm_summaries_extra_opts=cep2.slurm_summaries_extra_opts
						if not cmdline.opts.is_globalfs:
							slurm_summaries_extra_opts="%s -w %s" % (cep2.slurm_summaries_extra_opts, sumnode)
						docker_cmd_prefix=docker_cmd_suffix=""
						if cmdline.opts.is_docker:
							docker_cmd_prefix=cep2.docker_cmd_prefix
							docker_cmd_suffix="%s -e SLURM_JOB_ID=%s -e SLURM_JOB_NAME=%s %s" % \
								(cep2.docker_common_opts, cep2.slurm_jobid, cep2.slurm_jobname, cmdline.opts.docker_container)
						cmd="%ssrun %s -c %d --jobid=%s --job-name=%s %s %s %s'%spulp.py %s --noinit --summary --local --beams %s --confdir %s %s'%s" % \
							(docker_cmd_prefix, cep2.srun_general_opts, cmdline.opts.nthreads, cep2.slurm_jobid, cep2.slurm_jobname, slurm_summaries_extra_opts, \
							docker_cmd_suffix, cep2.start_shell, use_pulp_switch, cmdline.opts.is_auto and "--id %s" % (obs.id) or "", \
							sumcode, cmdline.opts.confdir, " ".join(cmdline.options), cep2.end_shell)
					else:
						docker_cmd_suffix=""
						if cmdline.opts.is_docker:
							docker_cmd_suffix="%s -e SLURM_JOB_ID=%s -e SLURM_JOB_NAME=%s %s" % \
								(cep2.docker_common_opts, cep2.slurm_jobid, cep2.slurm_jobname, cmdline.opts.docker_container)
						cmd="ssh -t %s %s '%spulp.py %s --noinit --summary --local --beams %s --confdir %s %s'" % \
							(sumnode, docker_cmd_suffix, use_pulp_switch, cmdline.opts.is_auto and "--id %s" % (obs.id) or "", \
							sumcode, cmdline.opts.confdir, " ".join(cmdline.options))
					sum_popen = Popen(shlex.split(cmd), stdout=PIPE, stderr=STDOUT)
					self.sum_popens.append(sum_popen)
					log.info("Making summaries on %s... (pid=%d)" % (sumnode, sum_popen.pid))

				run_units = [p.pid for p in self.sum_popens if p.poll() is None]
				Popen(shlex.split("stty sane"), stderr=open(os.devnull, 'rb')).wait()
				log.info("Still running [%d]: %s" % (len(run_units), run_units))
				for proc in self.sum_popens:
					Popen(shlex.split("stty sane"), stderr=open(os.devnull, 'rb')).wait()
					log.info("waiting...")
					proc.communicate()
					log.info("Process pid=%d has finished, status=%d" % (proc.pid, proc.returncode))
					run_units = [p.pid for p in self.sum_popens if p.poll() is None]
					finished_units = [p for p in self.sum_popens if p.poll() is not None]
					for fu in finished_units:
						if fu.returncode != 0: raise Exception
					if len(run_units) > 0: log.info("Still running [%d]: %s" % (len(run_units), run_units))

				# loop over finished summaries to see if they all finished OK
				failed_summaries = [s for s in self.sum_popens if s.returncode > 0]
				Popen(shlex.split("stty sane"), stderr=open(os.devnull, 'rb')).wait()
				self.number_failed_summaries = len(failed_summaries)
				log.info("%d failed summaries" % (self.number_failed_summaries))

		except Exception:
			log.exception("Oops... 'finish' function of the pipeline has crashed!")
			self.kill(log)
			raise

	# return number of failed pipelines
	def get_number_failed_pipes(self):
		return self.number_failed_pipes

	# return number of failed summaries
	def get_number_failed_summaries(self):
		return self.number_failed_summaries

	# execute command on local node (similar to execute in PipeUnit)
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
               			proc = Popen(shlex.split(cmd), stdout=PIPE, stderr=STDOUT, cwd=workdir, shell=shell)
       	        		self.log.process2log(proc)
        	       		proc.communicate()
				status = proc.poll()
			job_end = time.time()
			job_total_time = job_end - job_start
        	       	self.log.info("Finished at UTC %s, status=%s, Total running time: %.1f s (%.2f hrs)" % \
					(time.asctime(time.gmtime()), status, job_total_time, job_total_time/3600.))
			self.log.info("")
			# if job is not successful
			if status != 0:
				raise Exception
		except Exception:
			self.log.exception("Oops... job has crashed!\n%s\nStatus=%s" % (cmd, status))
			Popen(shlex.split("stty sane"), stderr=open(os.devnull, 'rb')).wait()
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
                                process = Popen(shlex.split(cmd), cwd=workdir, shell=shell)
                                time.sleep(10)  # waiting 10 secs to see if process crashes right away
                                if process.poll() is not None and process.poll() != 0:
                                        raise Exception
                                else: process.kill()  # if process is still running, it means that cmd is good, so we kill it in order to
                                                        # restart it with proper stdout/stderr and add it to the list
                        process = Popen(shlex.split(cmd), stdout=PIPE, stderr=STDOUT, cwd=workdir, shell=shell)
                        status=process.returncode
                        self.log.info("Job pid=%d, not waiting for it to finish." % (process.pid))
                        return process
                except Exception:
                        self.log.exception("Oops... job has crashed!\n%s\nStatus=%s" % (re.sub("\n", "\\\\n", cmd), status))
			Popen(shlex.split("stty sane"), stderr=open(os.devnull, 'rb')).wait()
                        raise Exception
	
	# waiting for the command to finish
        def waiting(self, prg, popen):
                """
                Waiting for process to finish
                """
                try:
                        job_start = time.time()
                        self.log.info("Waiting for %s to finish, pid=%d" % (prg, popen.pid))
                        (sout, serr) = popen.communicate()
                        # we pipe serr to sout, so no need to log serr
                        self.log.info(sout)
                        job_end = time.time()
                        job_total_time = job_end - job_start
                        self.log.info("Process pid=%d (%s) has finished at UTC %s, status=%d, Waiting time: %.1f s (%.2f hrs)" % \
                                (popen.pid, prg, time.asctime(time.gmtime()), popen.returncode, job_total_time, job_total_time/3600.))
                        self.log.info("")
                        # if job is not successful
                        if popen.poll() != 0:
                                raise Exception
                except Exception:
                        self.log.exception("Oops... %s has crashed!\npid=%d, Status=%s" % (prg, popen.pid, popen.returncode))
			Popen(shlex.split("stty sane"), stderr=open(os.devnull, 'rb')).wait()
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
                        self.log.info("Processes of %s have finished at UTC %s, Waiting time: %.1f s (%.2f hrs)" % \
                                (prg, time.asctime(time.gmtime()), job_total_time, job_total_time/3600.))
                        self.log.info("")
                except Exception:
                        self.log.exception("Oops... %s has crashed!\npids = %s" % (prg, ",".join(["%d" % (fu.pid) for fu in popen_list if fu.poll() is not None])))
			Popen(shlex.split("stty sane"), stderr=open(os.devnull, 'rb')).wait()
                        raise Exception


	# function that checks all processes in the list and kill them if they are still running
	def kill(self, log=None):
		if log != None: log.info("Killing all open processes...")
		for unit in self.units:
			if unit.parent != None and unit.parent.poll() is None:
				unit.parent.kill()
				if unit.parent != None: unit.parent.communicate()
				if unit.parent != None: unit.parent.poll()
		self.units = []
		# killing summary processes if open
		for sum_popen in self.sum_popens:
			if sum_popen != None and sum_popen.poll() is None:
				sum_popen.kill()
				if sum_popen != None: sum_popen.communicate()
				if sum_popen != None: sum_popen.poll()
		self.sum_popens=[]

	# make feedback file
	def make_feedback(self, obs, cep2, cmdline, log=None):
		self.log = log
		sumnode=cep2.get_current_node()
		sumcode = cmdline.opts.beam_str
		sumdir = self.summary_dirs[sumcode]

		# moving log-files to corresponding directory
		if self.log != None: self.log.info("Moving log-files...")
		if cep2.cluster_headnode[:5] == "head0" or cep2.cluster_headnode[:3] == "cpu": # CEP4
			workunits = [u for u in self.units if u.code == sumcode]
		else:
			workunits = [u for u in self.units if u.summary_node == sumnode and u.code == sumcode]
		for unit in workunits:
			if os.path.exists("%s/%s_sap%03d_beam%04d.log" % (sumdir, obs.id, unit.sapid, unit.tabid)):
				if not cmdline.opts.is_log_append:	
					cmd="mv -f %s_sap%03d_beam%04d.log %s/SAP%d/%s" % \
						(obs.id, unit.sapid, unit.tabid, unit.beams_root_dir, unit.sapid, unit.procdir)
					if self.log != None: self.execute(cmd, workdir=sumdir)
					else:
						proc = Popen(shlex.split(cmd), stdout=PIPE, stderr=STDOUT, cwd=sumdir)
						proc.communicate()
				else:
					# appending log from sumdir to the one in corresponing beam directory
					cmd="cat %s/%s_sap%03d_beam%04d.log >> %s/%s/SAP%d/%s/%s_sap%03d_beam%04d.log" % \
						(sumdir, obs.id, unit.sapid, unit.tabid, sumdir, unit.beams_root_dir, unit.sapid, unit.procdir, obs.id, unit.sapid, unit.tabid)
					if self.log != None: self.execute(cmd, is_os=True)
					else: os.system(cmd)
					# removing log from sumdir
					cmd="rm -f %s/%s_sap%03d_beam%04d.log" % (sumdir, obs.id, unit.sapid, unit.tabid)
					if self.log != None: self.execute(cmd, workdir=sumdir)
					else: os.system(cmd)

		tarname="%s%s%s%s" % (cep2.pipeid, self.summary_archive_prefix, sumcode, self.summary_archive_suffix)
		# updating the Feedback unit
		if cep2.cluster_headnode[:5] == "head0" or cep2.cluster_headnode[:3] == "cpu": # CEP4
			fbunit=[u for u in self.feedbacks if u.path == sumdir][0]
		else:
			fbunit=[u for u in self.feedbacks if u.node == sumnode and u.path == sumdir][0]
		fbunit.update("%s/%s" % (sumdir, tarname), "%s/.%s%s%s.fb" % (cep2.get_logdir(), cep2.pipeid, self.summary_archive_prefix, sumcode), sumcode, obs, True)
		fbunit.flush(100, cep2, True)


	# run necessary processes to organize summary info on summary nodes
	# to be run locally on summary node
	def make_summary(self, obs, cep2, cmdline, log):
		try:
			if cmdline.opts.is_feedback: 
				self.make_feedback(obs, cep2, cmdline, log)
			else:
				is_to_combine=False
				if cmdline.opts.is_nohoover: # are not using hoover nodes
					for u in self.units: # checking if we actually have data spread across several nodes
                        			is_iquv = False
			                        if (u.tab.is_coherent and obs.stokesCS == "IQUV") or (u.tab.is_coherent == False and obs.stokesIS == "IQUV"):
	                        		        is_iquv = True
						if (len(u.tab.location) > 1 or cmdline.opts.is_all_parts_at_once) and not (is_iquv and ((u.tab.is_coherent and obs.nsplitsCS == 1) \
						or (u.tab.is_coherent == False and obs.nsplitsIS == 1))):
							is_to_combine=True
							break
					if is_to_combine: # if we do have data across many nodes, we first call combine function before doing summaries
						if is_iquv and (obs.nsplitsCS > 1 or obs.nsplitsIS > 1):
							self.combine_parts_IQUV(obs, cep2, cmdline, log)
						else:
							self.combine_parts(obs, cep2, cmdline, log)
				if obs.CV:
					self.make_summary_CV(obs, cep2, cmdline, log, is_to_combine)
				else:
					self.make_summary_CS_IS(obs, cep2, cmdline, log, is_to_combine)

		except Exception:
			log.exception("Oops... 'make_summary' function on %s has crashed!" % (cep2.get_current_node()))
			raise


	# function that combines the ar-files from separate parts for each of the beams for the ObsID
	# and runs further processing in the combined file before summaries
	def combine_parts(self, obs, cep2, cmdline, log):
		try:
			self.log = log
			start_time=time.time()	
			sumnode=cep2.get_current_node()
			sumcode = cmdline.opts.beam_str
			sumdir = self.summary_dirs[sumcode]

			# start logging
			self.log.info("Combining parts for all relevant %s beams on %s:%s    UTC start time is: %s  @node: %s" % \
					(sumcode, sumnode, sumdir, time.asctime(time.gmtime()), sumnode))	

			# combining parts and do further processing for one beam at a time
			# first, we make a list of beams that require this combining. We can't use self.units list directly 
			# as it has element for each part
			if cep2.cluster_headnode[:5] == "head0" or cep2.cluster_headnode[:3] == "cpu": # CEP4
				workunits = [u for u in self.units if u.code == sumcode]
			else:
				workunits = [u for u in self.units if u.summary_node == sumnode and u.code == sumcode]
			beams_combine=[]
			for unit in workunits:
				if (len(unit.tab.location) > 1 or cmdline.opts.is_all_parts_at_once) and cmdline.opts.is_nohoover:
					beams_combine.append("%d:%d" % (unit.sapid, unit.tabid))
			beams_combine=np.unique(beams_combine)
			self.log.info("The parts of these beams will be combined: %s" % (", ".join(beams_combine)))

			# main loop for combining parts of the beam
			for beam in beams_combine:
				sapid=int(beam.split(":")[0])
				tabid=int(beam.split(":")[1])
				ref_unit = None # we will use it to access necessary info that is same for all parts (e.g. list of pulsars)
				for unit in self.units:
					if unit.sapid == sapid and unit.tabid == tabid:
						ref_unit = unit
						break
				if ref_unit == None:
					self.log.error("Can't find any processing units for the beam %s!" % (beam))
					raise Exception
				self.log.info("Combining parts for the beam %s..." % (beam))

				curdir = "%s/%s/SAP%d/%s" % (sumdir, ref_unit.beams_root_dir, sapid, ref_unit.procdir)
				output_prefix="%s_SAP%d_%s" % (obs.id, sapid, ref_unit.procdir)
				
				if ref_unit.tab.is_coherent:
					proc_subs = ref_unit.nrSubsPerFile * cmdline.opts.nsplitsCS
					if cmdline.opts.first_freq_splitCS * ref_unit.nrSubsPerFile + proc_subs > ref_unit.tab.nrSubbands:
						proc_subs -= (cmdline.opts.first_freq_splitCS * ref_unit.nrSubsPerFile + proc_subs - ref_unit.tab.nrSubbands)  
				else:
					proc_subs = ref_unit.nrSubsPerFile * cmdline.opts.nsplitsIS
					if cmdline.opts.first_freq_splitIS * ref_unit.nrSubsPerFile + proc_subs > ref_unit.tab.nrSubbands:
						proc_subs -= (cmdline.opts.first_freq_splitIS * ref_unit.nrSubsPerFile + proc_subs - ref_unit.tab.nrSubbands)  
				nsubs_eff = min(ref_unit.tab.nrSubbands, proc_subs)
				total_chan = nsubs_eff * ref_unit.nrChanPerSub

				# loop on pulsars
				if not cmdline.opts.is_plots_only and not cmdline.opts.is_feedback:
					for psr in ref_unit.psrs:
						self.log.info("Combining parts for pulsar %s..." % (psr))
	
        	                                # running psradd to add all freq channels together
                	                        self.log.info("Adding frequency splits together...")
                        	                ar_files=glob.glob("%s/%s_%s_P*.ar" % (curdir, psr, output_prefix))
						if len(ar_files) == 0:
							self.log.info("skipped")
						else:
                                	        	cmd="psradd -R -m time -o %s_%s.ar %s" % (psr, output_prefix, " ".join(ar_files))
                                        		self.execute(cmd, workdir=curdir)

							# check if number of ar-files is equal to number of splits
							if len(ar_files) < int(proc_subs/ref_unit.nrSubsPerFile):
								self.log.warning("Number of ar-files is smaller than number of splits!")
								proc_subs -= (int(proc_subs/ref_unit.nrSubsPerFile) - len(ar_files)) * ref_unit.nrSubsPerFile
								nsubs_eff = min(ref_unit.tab.nrSubbands, proc_subs)
								total_chan = nsubs_eff * ref_unit.nrChanPerSub

							# fixing coordinates in the ar-file
							fix_header_coords(self, ref_unit, psr, curdir, output_prefix)

							# running common DSPSR post-processing
							dspsr_postproc(self, ref_unit, cmdline, obs, psr, total_chan, nsubs_eff, curdir, output_prefix)
							
                                	       	        # removing ar-files from dspsr for every frequency split
                                        	       	if not cmdline.opts.is_debug:
                                                	       	remove_list=glob.glob("%s/%s_%s_P*.ar" % (curdir, psr, output_prefix))
                                                        	cmd="rm -f %s" % (" ".join(remove_list))
       	                                                	self.execute(cmd, workdir=curdir)

						# combining fil-files after digifil if we do single-pulse analysis for CV data only
						if cmdline.opts.is_single_pulse and obs.CV:
							try:
								self.log.info("Adding frequency splits together for single-pulse analysis...")
								# we sort the list of fil-files in decrease frequency order (needed for sigproc_splice)
                        	                		fil_files=sorted(glob.glob("%s/%s_%s_P*.fil" % (curdir, psr, output_prefix)), key=lambda x: int(x.split("_PART")[-1].split(".fil")[0]), reverse=True)
								if len(fil_files) == 0:
									self.log.info("skipped")
								else:
									cmd="sigproc_splice -a -o %s_%s.fil %s" % (psr, output_prefix, " ".join(fil_files))
									self.execute(cmd, workdir=curdir)
									if not cmdline.opts.is_debug:
										cmd="rm -f %s" % (" ".join(fil_files))
										self.execute(cmd, workdir=curdir)
									# running rfifind
									if not cmdline.opts.is_norfi:	
										self.log.info("Running rfifind for pulsar %s on a filterbank file...")
										cmd="rfifind -o %s_%s -noclip -blocks 16 %s %s_%s.fil" % (psr, output_prefix, cmdline.opts.rfifind_extra_opts, psr, output_prefix)
										self.execute(cmd, workdir=curdir)

									# running prepdata and single_pulse_search.py
									run_prepdata_CV(self, ref_unit, cmdline, sumdir, curdir, output_prefix)
							except: pass

				# running pav
				if not cmdline.opts.is_skip_dspsr and not cmdline.opts.is_feedback:
					make_dspsr_plots(self, ref_unit, cmdline, obs, nsubs_eff, curdir, output_prefix)

				# Running pdmp
				if not cmdline.opts.is_plots_only and not cmdline.opts.is_feedback and \
				not cmdline.opts.is_skip_dspsr and not cmdline.opts.is_nopdmp and not cmdline.opts.is_nofold:
					pdmp_popens = start_pdmp(self, ref_unit, cmdline, obs, nsubs_eff, curdir, output_prefix)
					# waiting for pdmp to finish
					finish_pdmp(self, ref_unit, pdmp_popens, cmdline, obs, curdir, output_prefix)
					# making diagnostic plot after pdmp
					make_dspsr_plots(self, ref_unit, cmdline, obs, nsubs_eff, curdir, output_prefix, True)

				# Running spectar.py to calculate pulsar spectra
				if not cmdline.opts.is_skip_dspsr and not cmdline.opts.is_feedback:
					calc_psr_spectra(self, ref_unit, cmdline, obs, sapid, tabid, curdir, output_prefix)

		except Exception:
			self.log.exception("Oops... 'combine_parts' function on %s has crashed" % (cep2.get_current_node()))
			raise

	# function that combines the ar-files from separate parts for each of the beams for the ObsID
	# and for each Stokes (I, Q, U, and V) separately
	# and runs further processing in the combined file before summaries
	def combine_parts_IQUV(self, obs, cep2, cmdline, log):
		try:
			self.log = log
			start_time=time.time()	
			sumnode=cep2.get_current_node()
			sumcode = cmdline.opts.beam_str
			sumdir = self.summary_dirs[sumcode]

			# start logging
			self.log.info("Combining parts for all relevant %s beams on %s:%s    UTC start time is: %s  @node: %s" % \
					(sumcode, sumnode, sumdir, time.asctime(time.gmtime()), sumnode))	

			# combining parts and do further processing for one beam at a time
			# first, we make a list of beams that require this combining. We can't use self.units list directly 
			# as it has element for each part
			if cep2.cluster_headnode[:5] == "head0" or cep2.cluster_headnode[:3] == "cpu": # CEP4
				workunits = [u for u in self.units if u.code == sumcode]
			else:
				workunits = [u for u in self.units if u.summary_node == sumnode and u.code == sumcode]
			beams_combine=[]
			for unit in workunits:
				if (len(unit.tab.location) > 1 or cmdline.opts.is_all_parts_at_once) and cmdline.opts.is_nohoover:
					beams_combine.append("%d:%d" % (unit.sapid, unit.tabid))
			beams_combine=np.unique(beams_combine)
			self.log.info("The parts of these beams will be combined: %s" % (", ".join(beams_combine)))

			# main loop for combining parts of the beam
			for beam in beams_combine:
				sapid=int(beam.split(":")[0])
				tabid=int(beam.split(":")[1])
				ref_unit = None # we will use it to access necessary info that is same for all parts (e.g. list of pulsars)
				for unit in self.units:
					if unit.sapid == sapid and unit.tabid == tabid:
						ref_unit = unit
						break
				if ref_unit == None:
					self.log.error("Can't find any processing units for the beam %s!" % (beam))
					raise Exception
				self.log.info("Combining parts for the beam %s..." % (beam))

				curdir = "%s/%s/SAP%d/%s" % (sumdir, ref_unit.beams_root_dir, sapid, ref_unit.procdir)
				if ref_unit.tab.is_coherent:
					proc_subs = ref_unit.nrSubsPerFile * cmdline.opts.nsplitsCS
					if cmdline.opts.first_freq_splitCS * ref_unit.nrSubsPerFile + proc_subs > ref_unit.tab.nrSubbands:
						proc_subs -= (cmdline.opts.first_freq_splitCS * ref_unit.nrSubsPerFile + proc_subs - ref_unit.tab.nrSubbands)  
				else:
					proc_subs = ref_unit.nrSubsPerFile * cmdline.opts.nsplitsIS
					if cmdline.opts.first_freq_splitIS * ref_unit.nrSubsPerFile + proc_subs > ref_unit.tab.nrSubbands:
						proc_subs -= (cmdline.opts.first_freq_splitIS * ref_unit.nrSubsPerFile + proc_subs - ref_unit.tab.nrSubbands)  
				nsubs_eff = min(ref_unit.tab.nrSubbands, proc_subs)
				total_chan = nsubs_eff * ref_unit.nrChanPerSub

				# loop on different stokes
				for pol in np.arange(4):
	
					self.log.info("Combining parts for the Stokes %d..." % (pol))
					output_prefix="%s_SAP%d_%s_S%d" % (obs.id, sapid, ref_unit.procdir, pol)

					if not cmdline.opts.is_plots_only and not cmdline.opts.is_feedback:

						# loop on pulsars
						for psr in ref_unit.psrs:
							self.log.info("Combining parts for pulsar %s..." % (psr))
		
       		        	                        # running psradd to add all freq channels together
               		        	                self.log.info("Adding frequency splits together...")
                       		        	        ar_files=glob.glob("%s/%s_%s_P*.ar" % (curdir, psr, output_prefix))
							if len(ar_files) == 0:
								self.log.info("skipped")
							else:
       		                        	        	cmd="psradd -R -m time -o %s_%s.ar %s" % (psr, output_prefix, " ".join(ar_files))
               		                        		self.execute(cmd, workdir=curdir)

								# check if number of ar-files is equal to number of splits
								if len(ar_files) < int(proc_subs/ref_unit.nrSubsPerFile):
									self.log.warning("Number of ar-files is smaller than number of splits!")
									proc_subs -= (int(proc_subs/ref_unit.nrSubsPerFile) - len(ar_files)) * ref_unit.nrSubsPerFile
									nsubs_eff = min(ref_unit.tab.nrSubbands, proc_subs)
									total_chan = nsubs_eff * ref_unit.nrChanPerSub

								# fixing coordinates in the ar-file
								fix_header_coords(self, ref_unit, psr, curdir, output_prefix)

								# running common DSPSR post-processing
								dspsr_postproc(self, ref_unit, cmdline, obs, psr, total_chan, nsubs_eff, curdir, output_prefix)
							
                               			       	        # removing ar-files from dspsr for every frequency split
                                       			       	if not cmdline.opts.is_debug:
                                               			       	remove_list=glob.glob("%s/%s_%s_P*.ar" % (curdir, psr, output_prefix))
                                                       			cmd="rm -f %s" % (" ".join(remove_list))
                                                       			self.execute(cmd, workdir=curdir)

					# running pav
					if not cmdline.opts.is_skip_dspsr and not cmdline.opts.is_feedback:
						make_dspsr_plots(self, ref_unit, cmdline, obs, nsubs_eff, curdir, output_prefix)

					# Running pdmp
					if not cmdline.opts.is_plots_only and not cmdline.opts.is_feedback and \
					not cmdline.opts.is_skip_dspsr and not cmdline.opts.is_nopdmp and not cmdline.opts.is_nofold:
						pdmp_popens = start_pdmp(self, ref_unit, cmdline, obs, nsubs_eff, curdir, output_prefix)
						# waiting for pdmp to finish
						finish_pdmp(self, ref_unit, pdmp_popens, cmdline, obs, curdir, output_prefix)
						# making diagnostic plot after pdmp
						make_dspsr_plots(self, ref_unit, cmdline, obs, nsubs_eff, curdir, output_prefix, True)

					# Running spectar.py to calculate pulsar spectra
					if not cmdline.opts.is_skip_dspsr and not cmdline.opts.is_feedback:
						calc_psr_spectra(self, ref_unit, cmdline, obs, sapid, tabid, curdir, output_prefix)

		except Exception:
			self.log.exception("Oops... 'combine_parts_IQUV' function on %s has crashed" % (cep2.get_current_node()))
			raise


	# run necessary processes to organize summary info on summary nodes for CV data
	# to be run locally on summary node
	def make_summary_CV(self, obs, cep2, cmdline, log, is_to_combine=False):

		self.log = log
		start_time=time.time()	
		sumnode=cep2.get_current_node()
		data_code = cmdline.opts.beam_str
		sumdir = self.summary_dirs[data_code]

		# start logging
		self.log.info("Summaries on %s:%s    UTC start time is: %s  @node: %s" % \
				(sumnode, sumdir, time.asctime(time.gmtime()), sumnode))

		# extracting files from archive files for all beams
		# and moving log-files to corresponding directory
		if cep2.cluster_headnode[:5] == "head0" or cep2.cluster_headnode[:3] == "cpu": # CEP4
			workunits = [u for u in self.units if u.code == data_code]
		else:
			workunits = [u for u in self.units if u.summary_node == sumnode and u.code == data_code]
		if not cmdline.opts.is_globalfs: # not CEP4
			self.log.info("Extracting archives in summary nodes, removing archives, moving log-files...")
			already_extracted=[] # list of archives already extracted (to prevent the same archive to be extracted many times in case multiple parts)
			for unit in workunits:
				result_archive="%s_SAP%03d_B%03d_P%03d%s" % (cep2.pipeid, unit.sapid, unit.tabid, unit.part, unit.archive_suffix)
	                        if result_archive in already_extracted: continue
        	                else: already_extracted.append(result_archive)
				if os.path.exists("%s/%s" % (sumdir, result_archive)):
					# extracting archive
					cmd="tar xvf %s" % (result_archive)
					self.execute(cmd, workdir=sumdir)
					# removing archive
					if cmdline.opts.nsplitsCS < 2:
						cmd="rm -f %s" % (result_archive)
						self.execute(cmd, workdir=sumdir)
					# moving log-file to corresponding SAP/BEAM directory
					if os.path.exists("%s/%s_sap%03d_beam%04d.log" % (sumdir, obs.id, unit.sapid, unit.tabid)) and \
						os.path.exists("%s/%s/SAP%d/%s" % (sumdir, unit.beams_root_dir, unit.sapid, unit.procdir)):
						if not cmdline.opts.is_log_append:	
							cmd="mv -f %s_sap%03d_beam%04d.log %s/SAP%d/%s" % \
								(obs.id, unit.sapid, unit.tabid, unit.beams_root_dir, unit.sapid, unit.procdir)
							self.execute(cmd, workdir=sumdir)
						else:
							# appending log from sumdir to the one in corresponing beam directory
							cmd="cat %s/%s_sap%03d_beam%04d.log >> %s/%s/SAP%d/%s/%s_sap%03d_beam%04d.log" % \
								(sumdir, obs.id, unit.sapid, unit.tabid, sumdir, unit.beams_root_dir, unit.sapid, unit.procdir, obs.id, unit.sapid, unit.tabid)
							self.execute(cmd, is_os=True)
							# removing log from sumdir
							cmd="rm -f %s_sap%03d_beam%04d.log" % (obs.id, unit.sapid, unit.tabid)
							self.execute(cmd, workdir=sumdir)
				else:
					if not os.path.exists("%s/%s" % (sumdir, unit.curdir.split(unit.outdir + "/")[1])):
						self.log.warning("Warning! Neither archive file %s nor corresponding directory tree exists in: %s. Summary won't be complete" % (result_archive, sumdir))
		else: # CEP4 with GlobalFS (we do not need to extract archives but only moving log files)
			self.log.info("Moving log-files...")
			for unit in workunits:
				# moving log-file to corresponding SAP/BEAM directory
				if os.path.exists("%s/%s_sap%03d_beam%04d.log" % (sumdir, obs.id, unit.sapid, unit.tabid)) and \
					os.path.exists("%s/%s/SAP%d/%s" % (sumdir, unit.beams_root_dir, unit.sapid, unit.procdir)):
					if not cmdline.opts.is_log_append:	
						cmd="mv -f %s_sap%03d_beam%04d.log %s/SAP%d/%s" % \
							(obs.id, unit.sapid, unit.tabid, unit.beams_root_dir, unit.sapid, unit.procdir)
						self.execute(cmd, workdir=sumdir)
					else:
						# appending log from sumdir to the one in corresponing beam directory
						cmd="cat %s/%s_sap%03d_beam%04d.log >> %s/%s/SAP%d/%s/%s_sap%03d_beam%04d.log" % \
							(sumdir, obs.id, unit.sapid, unit.tabid, sumdir, unit.beams_root_dir, unit.sapid, unit.procdir, obs.id, unit.sapid, unit.tabid)
						self.execute(cmd, is_os=True)
						# removing log from sumdir
						cmd="rm -f %s_sap%03d_beam%04d.log" % (obs.id, unit.sapid, unit.tabid)
						self.execute(cmd, workdir=sumdir)

		# create beam_process_node.txt file (this is only if file does not exist or it is empty)
		beam_process_node_file="%s/beam_process_node.txt" % (sumdir)
		if not os.path.exists(beam_process_node_file) or os.path.getsize(beam_process_node_file) == 0:
			self.log.info("Creating the beam_process_node.txt file...")
        	        bpnf=open(beam_process_node_file, 'w')
			for unit in workunits:
				for node in unit.tab.location:
					if node in unit.tab.rawfiles:
						for rf in unit.tab.rawfiles[node]:
							bpnf.write("%s %s%s\n" % \
							(node, rf, unit.tab.specificationType == "flyseye" and " [%s]" % (",".join(unit.tab.stationList)) or ""))
			bpnf.close()

		# creating combined DSPSR plots
		# first check if there are diagnostic plots after pdmp
		dspsr_diags=rglob(sumdir, "*_diag_pdmp.png", 3)
		# if not, then we check for diagnostic plots before pdmp
		if len(dspsr_diags) == 0:
			dspsr_diags=rglob(sumdir, "*_diag.png", 3)
		if len(dspsr_diags) > 0:
			self.log.info("Creating DSPSR summary diagnostic plots...")
			if len(dspsr_diags) > 1: cmd="montage %s -background none -mode concatenate -tile %dx dspsr_status.png" % (" ".join(dspsr_diags), int(math.sqrt(len(dspsr_diags))))
			else: cmd="cp -f %s dspsr_status.png" % (dspsr_diags[0])
			self.execute(cmd, workdir=sumdir)

		if os.path.exists("%s/dspsr_status.png" % (sumdir)):
			self.log.info("Copying dspsr status file to status.png ...")
			cmd="mv dspsr_status.png status.png"
			self.execute(cmd, workdir=sumdir)

		# creating thumbnail version of status.png if it exists
		if os.path.exists("%s/status.png" % (sumdir)):		
			self.log.info("Making a thumbnail version of status.png file...")
			cmd="convert -scale 200x140-0-0 status.png status.th.png"
			self.execute(cmd, workdir=sumdir)
		else:
			self.log.info("No status.png created")

		# creating TA heatmaps 
		# only when folding, and only if pulsars are set from the command line, or 'parset' or 'sapfind' or 'sapfind3' or "tabfind+" keywords are used (or
		# nothing is given for --pulsar option)
		# otherwise, different TA beams will be folded for different pulsars, and TA heatmap does not have sense
		if data_code == "CV" and not cmdline.opts.is_nofold and (len(cmdline.psrs) == 0 or (len(cmdline.psrs) != 0 and cmdline.psrs[0] != "tabfind")):
			for sap in obs.saps:
				# getting number of _coherent_ TABs
				nrTABs = len([1 for tab in sap.tabs if tab.is_coherent])
				if sap.nrRings > 0 or nrTABs > 1:
					if len(cmdline.psrs) != 0 and cmdline.psrs[0] != "parset" and cmdline.psrs[0] != "sapfind" and \
										cmdline.psrs[0] != "sapfind3" and cmdline.psrs[0] != "tabfind+":
						psrs = cmdline.psrs # getting list of pulsars from command line
					else: 
						if len(cmdline.psrs) == 0:
							if sap.source != "" and check_pulsars(sap.source, cmdline, cep2, None): psrs = [sap.source]
							else:
								if len(sap.psrs) > 0: psrs = [sap.psrs[0]]
								else: psrs = []
						else:
							if cmdline.psrs[0] == "parset":
								if sap.source != "" and check_pulsars(sap.source, cmdline, cep2, None): psrs = [sap.source]
								else: psrs = []
							if cmdline.psrs[0] == "sapfind" or cmdline.psrs[0] == "sapfind3":
								if len(sap.psrs) > 0:
									if cmdline.psrs[0] == "sapfind": psrs = [sap.psrs[0]]
									else: psrs = sap.psrs
								else: psrs = []
							if cmdline.psrs[0] == "tabfind+":
								if sap.source != "" and check_pulsars(sap.source, cmdline, cep2, None): psrs = [sap.source]
								else: psrs = []
								if len(sap.psrs) > 0: psrs += sap.psrs
								psrs = list(np.unique(psrs))
					if sap.nrRings > 0: self.log.info("Creating TA heatmap with %d rings for SAP=%d..." % (sap.nrRings, sap.sapid))
					else: self.log.info("Creating TA heatmap with %d TA beams for SAP=%d..." % (nrTABs, sap.sapid))

					# here creating chi-squared.txt file base on .fscr.AR files and using psrstat to get S/N
					# first check for paz files if we do rfi zapping
					arfs=rglob(sumdir, "*.paz.fscr.AR", 3)
					# if not, then we check for diagnostic plots before pdmp
					if len(arfs) == 0:
						arfs=rglob(sumdir, "*.fscr.AR", 3)
					if len(arfs) > 0:
                        			chif=open("%s/chi-squared.txt" % (sumdir), 'w')
						for arf in arfs:
							try:
                                        			psr=arf.split("/")[-1].split("_")[0]
                                        			cursapid=int(arf.split("_SAP")[-1].split("_")[0])
                                				curprocdir=arf.split("_SAP")[-1].split("_")[1].split(".")[0]
								cmd="psrstat -Q -j FTp -c snr %s" % (arf)
								snrline=os.popen(cmd).readlines()
								if np.size(snrline) > 0:
									chif.write("file=%s obs=%s_SAP%d_%s_%s S/N=%s\n" % (arf, data_code, cursapid, curprocdir, psr, snrline[0].split()[-1]))
							except: pass
						chif.close()

					for psr in psrs:
						self.log.info(psr)
						# I need this try/except block here, to avoid situation, when there are more than 1 pulsar in the SAP
						# but processing was done only for 1 - this is the case when pulsar is not specified in the command line
						try:
							cmd="cat %s/chi-squared.txt | grep _SAP%d | grep %s > %s/%s-chi-squared.txt" % (sumdir, sap.sapid, psr, sumdir, psr)
							self.execute(cmd, is_os=True)
							if not cmdline.opts.is_cobalt:	
								cmd="plot_LOFAR_TA_multibeam3.py --sap %d --chi %s-chi-squared.txt --parset %s.parset --out_logscale %s_SAP%d_%s_TA_heatmap_log.png --out_linscale %s_SAP%d_%s_TA_heatmap_linear.png --target %s" % (sap.sapid, psr, obs.id, obs.id, sap.sapid, psr, obs.id, sap.sapid, psr, psr)
							else: # for Cobalt we first need to make pseudo parset file
								obs.pseudo_parset_generator("%s/%s.pseudo.parset" % (sumdir, obs.id))
								cmd="plot_LOFAR_TA_multibeam3.py --sap %d --chi %s-chi-squared.txt --parset %s.pseudo.parset --out_logscale %s_SAP%d_%s_TA_heatmap_log.png --out_linscale %s_SAP%d_%s_TA_heatmap_linear.png --target %s" % (sap.sapid, psr, obs.id, obs.id, sap.sapid, psr, obs.id, sap.sapid, psr, psr)
							self.execute(cmd, workdir=sumdir)
							cmd="rm -f %s-chi-squared.txt" % (psr)
							self.execute(cmd, workdir=sumdir)
							# combining TA heatmap log and linear plots
							cmd="convert %s_SAP%d_%s_TA_heatmap_log.png %s_SAP%d_%s_TA_heatmap_linear.png -append ta_heatmap_sap%d_%s.png" % (obs.id, sap.sapid, psr, obs.id, sap.sapid, psr, sap.sapid, psr)
							self.execute(cmd, workdir=sumdir)
						except: 
							self.log.info("Can't make a heatmap plot for pulsar %s" % (psr))
							# removing temporary PSR-chi-squared.txt file
							cmd="rm -f %s/%s-chi-squared.txt" % (sumdir, psr)
							self.execute(cmd, workdir=sumdir)

						# if we do single-pulse analysis then we try to do RRATs heatmap as well
						if cmdline.opts.is_rrats or cmdline.opts.is_single_pulse:
							try:
								psr2=re.sub(r'^[BJ]', '', psr)
								parfile="%s/%s.par" % (sumdir, psr2)
								ref_unit=workunits[0]
								if cmdline.opts.is_nofold:
									psrdm = 0.0
								else:
									if not os.path.exists(parfile):
										parfiles=glob.glob("%s/*/SAP*/BEAM*/*.par" % (sumdir))
										if len(parfiles) > 0: 
											psrdm=ref_unit.get_psr_dm(parfiles[0])
										else: psrdm=0.0
									else: psrdm=ref_unit.get_psr_dm(parfile)
								if not cmdline.opts.is_cobalt:	
									cmd="RRAT_heatmap2.py --parset %s.parset --target %s --dm %f --sap %d --out_logscale %s_SAP%d_%s_RRAT_heatmap_log.png --out_linscale %s_SAP%d_%s_RRAT_heatmap_linear.png" % (obs.id, psr, psrdm, sap.sapid, obs.id, sap.sapid, psr, obs.id, sap.sapid, psr)
								else: # for Cobalt we first need to make pseudo parset file
									obs.pseudo_parset_generator("%s/%s.pseudo.parset" % (sumdir, obs.id))
									cmd="RRAT_heatmap2.py --parset %s.pseudo.parset --target %s --dm %f --sap %d --out_logscale %s_SAP%d_%s_RRAT_heatmap_log.png --out_linscale %s_SAP%d_%s_RRAT_heatmap_linear.png" % (obs.id, psr, psrdm, sap.sapid, obs.id, sap.sapid, psr, obs.id, sap.sapid, psr)
								self.execute(cmd, workdir=sumdir)
								# combining RRAT heatmap log and linear plots
								cmd="convert %s_SAP%d_%s_RRAT_heatmap_log.png %s_SAP%d_%s_RRAT_heatmap_linear.png -append %s_sap%d_rratmap.png" % (obs.id, sap.sapid, psr, obs.id, sap.sapid, psr, psr, sap.sapid)
								self.execute(cmd, workdir=sumdir)
							except:
								self.log.info("Can't make a RRAT heatmap plot for pulsar %s in SAP%d" % (psr, sap.sapid))

					# combining TA heatmaps for different pulsars
					heatmaps=glob.glob("%s/ta_heatmap_sap%d_*.png" % (sumdir, sap.sapid))
					if len(heatmaps) > 0:
						if len(heatmaps) > 1: cmd="convert %s +append ta_heatmap_sap%d.png" % (" ".join(heatmaps), sap.sapid)
						else: cmd="mv %s ta_heatmap_sap%d.png" % (heatmaps[0], sap.sapid)
						self.execute(cmd, workdir=sumdir)
						# remove temporary png files
						cmd="rm -f %s" % (" ".join(heatmaps))
						self.execute(cmd, workdir=sumdir)

					# combining RRAT heatmap for different pulsars
					if cmdline.opts.is_rrats or cmdline.opts.is_single_pulse:
						rratmaps=glob.glob("%s/*_sap%d_rratmap.png" % (sumdir, sap.sapid))
						if len(rratmaps) > 0:
							if len(rratmaps) > 1: cmd="convert %s +append rratmap_sap%d.png" % (" ".join(rratmaps), sap.sapid)
							else: cmd="mv %s rratmap_sap%d.png" % (rratmaps[0], sap.sapid)
							self.execute(cmd, workdir=sumdir)
							# remove temporary png files
							cmd="rm -f %s" % (" ".join(rratmaps))
							self.execute(cmd, workdir=sumdir)

			# combining TA heatmaps for different SAPs
			heatmaps=glob.glob("%s/ta_heatmap_sap*.png" % (sumdir))
			if len(heatmaps) > 0:
				if len(heatmaps) > 1: cmd="convert %s -append TAheatmap_status.png" % (" ".join(heatmaps))
				else: cmd="mv %s TAheatmap_status.png" % (heatmaps[0])
				self.execute(cmd, workdir=sumdir)
				# remove temporary png files
				cmd="rm -f %s" % (" ".join(heatmaps))
				self.execute(cmd, workdir=sumdir)
				# making a thumbnail version of TA heatmap combined plot
				cmd="convert -scale 200x140-0-0 TAheatmap_status.png TAheatmap_status.th.png"
				self.execute(cmd, workdir=sumdir)

			# combining RRAT heatmap for different SAPs
			if cmdline.opts.is_rrats or cmdline.opts.is_single_pulse:
				rratmaps=glob.glob("%s/rratmap_sap*.png" % (sumdir))
				if len(rratmaps) > 0:
					if len(rratmaps) > 1: cmd="convert %s -append rratmap.png" % (" ".join(rratmaps))
					else: cmd="mv %s rratmap.png" % (rratmaps[0])
					self.execute(cmd, workdir=sumdir)
					# remove temporary png files
					cmd="rm -f %s" % (" ".join(rratmaps))
					self.execute(cmd, workdir=sumdir)

		# Make a tarball of all the plots (summary archive)
		if is_to_combine:
			extensions = self.summary_archive_exts
		else: 
			extensions = self.summary_archive_exts_nocombine
		self.log.info("Making a final summary tarball of all files with extensions: %s" % (", ".join(extensions)))
		tarname="%s%s%s%s" % (cep2.pipeid, self.summary_archive_prefix, data_code, self.summary_archive_suffix)
		tar_list=[]
		for ext in extensions:
			ext_list=rglob(sumdir, ext, 3)
			tar_list.extend(ext_list)
		cmd="tar -cv --ignore-failed-read -f %s %s" % (tarname, " ".join([f.split(sumdir+"/")[1] for f in tar_list]))
		try: # --ignore-failed-read does not seem to help with tar failing for some beams
                     # like file was changed during the tar, though tarball seem to be fine
			self.execute(cmd, workdir=sumdir)
		except: pass

		# finish
		end_time=time.time()
		total_time= end_time- start_time
		self.log.info("UTC stop time is: %s" % (time.asctime(time.gmtime())))
		self.log.info("Total running time: %.1f s (%.2f hrs)" % (total_time, total_time/3600.))

		# flushing log file and copy it to summary node
		self.log.flush()
		if not cmdline.opts.is_log_append: cmd="cp -f %s %s" % (cep2.get_logfile(), sumdir)
		else: cmd="cat %s >> %s/%s" % (cep2.get_logfile(), sumdir, cep2.get_logfile().split("/")[-1])
		os.system(cmd)

		# delete file from the archive first
		cmd="tar --delete -v --ignore-failed-read -f %s %s" % (tarname, cep2.get_logfile().split("/")[-1])
		try:
			proc = Popen(shlex.split(cmd), stdout=PIPE, stderr=STDOUT, cwd=sumdir)
			proc.communicate()
			# adding log file to the archive and gzip it
			cmd="tar -rv --ignore-failed-read -f %s %s" % (tarname, cep2.get_logfile().split("/")[-1])
			# --ignore-failed-read does not seem to help with tar failing for some beams
			# like file was changed during the tar, though tarball seem to be fine
			proc = Popen(shlex.split(cmd), stdout=PIPE, stderr=STDOUT, cwd=sumdir)
			proc.communicate()
		except: pass
		# avoid gzipping now
		#cmd="gzip -f %s" % (tarname)
		#proc = Popen(shlex.split(cmd), stdout=PIPE, stderr=STDOUT, cwd=sumdir)
		#proc.communicate()

		# updating the Feedback unit
		self.make_feedback(obs, cep2, cmdline)

		# specific to dragnet
		if cep2.cluster_headnode == "dragnet":
			# changing ownership to 'dragnet' group
			cmd="chgrp -R dragnet %s" % (sumdir)
			os.system(cmd)
		# changing the file permissions to be re-writable for group
		cmd="chmod -R g+w %s" % (sumdir)
		os.system(cmd)

	# run necessary processes to organize summary info on summary nodes for CS and IS data
	# to be run locally on summary node
	def make_summary_CS_IS(self, obs, cep2, cmdline, log, is_to_combine=False):

		self.log = log
		start_time=time.time()	
		sumnode=cep2.get_current_node()
		data_code = cmdline.opts.beam_str
		sumdir = self.summary_dirs[data_code]

		# start logging
		self.log.info("Summaries on %s:%s    UTC start time is: %s  @node: %s" % \
				(sumnode, sumdir, time.asctime(time.gmtime()), sumnode))

		# extracting files from archive files for all beams
		# and moving log-files to corresponding directory
		if cep2.cluster_headnode[:5] == "head0" or cep2.cluster_headnode[:3] == "cpu": # CEP4
			workunits = [u for u in self.units if u.code == data_code]
		else:
			workunits = [u for u in self.units if u.summary_node == sumnode and u.code == data_code]
		if not cmdline.opts.is_globalfs: # not CEP4
			self.log.info("Extracting archives in summary nodes, removing archives, moving log-files...")
			already_extracted=[] # list of archives already extracted (to prevent the same archive to be extracted many times in case multiple parts)
			for unit in workunits:
                		is_iquv = False
                		if (unit.tab.is_coherent and obs.stokesCS == "IQUV") or (unit.tab.is_coherent == False and obs.stokesIS == "IQUV"):
                        		is_iquv = True
				if is_iquv:
					result_archive="%s_SAP%03d_B%03d_S%d_P%03d%s" % (cep2.pipeid, unit.sapid, unit.tabid, unit.stokes_index, unit.part, unit.archive_suffix)
				else:
					result_archive="%s_SAP%03d_B%03d_P%03d%s" % (cep2.pipeid, unit.sapid, unit.tabid, unit.part, unit.archive_suffix)
				if result_archive in already_extracted: continue
				else: already_extracted.append(result_archive)
				if os.path.exists("%s/%s" % (sumdir, result_archive)):
					# extracting archive
					cmd="tar xvf %s" % (result_archive)
					self.execute(cmd, workdir=sumdir)
					# removing archive
					if data_code == "IS":
						if cmdline.opts.nsplitsIS < 2:
							cmd="rm -f %s" % (result_archive)
							self.execute(cmd, workdir=sumdir)
					else:
						if cmdline.opts.nsplitsCS < 2:
							cmd="rm -f %s" % (result_archive)
							self.execute(cmd, workdir=sumdir)
					# moving log-file to corresponding SAP/BEAM directory
					if is_iquv:
						beamlog="%s_sap%03d_beam%04d_stokes%d.log" % (obs.id, unit.sapid, unit.tabid, unit.stokes_index)
					else:
						beamlog="%s_sap%03d_beam%04d.log" % (obs.id, unit.sapid, unit.tabid)
					if os.path.exists("%s/%s" % (sumdir, beamlog)) and \
						os.path.exists("%s/%s/SAP%d/%s" % (sumdir, unit.beams_root_dir, unit.sapid, unit.procdir)):
						if not cmdline.opts.is_log_append:	
							cmd="mv -f %s %s/SAP%d/%s" % \
								(beamlog, unit.beams_root_dir, unit.sapid, unit.procdir)
							self.execute(cmd, workdir=sumdir)
						else:
							# appending log from sumdir to the one in corresponing beam directory
							cmd="cat %s/%s >> %s/%s/SAP%d/%s/%s" % \
								(sumdir, beamlog, sumdir, unit.beams_root_dir, unit.sapid, unit.procdir, beamlog)
							self.execute(cmd, is_os=True)
							# removing log from sumdir
							cmd="rm -f %s" % (beamlog)
							self.execute(cmd, workdir=sumdir)
				else:
					if not os.path.exists("%s/%s" % (sumdir, unit.curdir.split(unit.outdir + "/")[1])):
						self.log.warning("Warning! Neither archive file %s nor corresponding directory tree exists in: %s. Summary won't be complete" % (result_archive, sumdir))

		else: # CEP4 with GlobalFS (we do not need to extract archives but only moving log files)
			self.log.info("Moving log-files...")
			# moving log-file to corresponding SAP/BEAM directory
			for unit in workunits:
                		is_iquv = False
                		if (unit.tab.is_coherent and obs.stokesCS == "IQUV") or (unit.tab.is_coherent == False and obs.stokesIS == "IQUV"):
                        		is_iquv = True
				if is_iquv:
					beamlog="%s_sap%03d_beam%04d_stokes%d.log" % (obs.id, unit.sapid, unit.tabid, unit.stokes_index)
				else:
					beamlog="%s_sap%03d_beam%04d.log" % (obs.id, unit.sapid, unit.tabid)
				if os.path.exists("%s/%s" % (sumdir, beamlog)) and \
					os.path.exists("%s/%s/SAP%d/%s" % (sumdir, unit.beams_root_dir, unit.sapid, unit.procdir)):
					if not cmdline.opts.is_log_append:	
						cmd="mv -f %s %s/SAP%d/%s" % \
							(beamlog, unit.beams_root_dir, unit.sapid, unit.procdir)
						self.execute(cmd, workdir=sumdir)
					else:
						# appending log from sumdir to the one in corresponing beam directory
						cmd="cat %s/%s >> %s/%s/SAP%d/%s/%s" % \
							(sumdir, beamlog, sumdir, unit.beams_root_dir, unit.sapid, unit.procdir, beamlog)
						self.execute(cmd, is_os=True)
						# removing log from sumdir
						cmd="rm -f %s" % (beamlog)
						self.execute(cmd, workdir=sumdir)

		# extracting *_sp_singlepulse.tgz in case it exists
		# for example, when we re-run only summaries, then there are no already corresponding *.singlepulse files in the  BEAM* directories
		if cmdline.opts.is_rrats or cmdline.opts.is_single_pulse:
			tarname="%s_sp_singlepulse.tar.gz" % (obs.id)
			if os.path.exists("%s/%s" % (sumdir, tarname)):
				self.log.info("Extracting *.singlepulse files to make a RRAT heatmap...")
				cmd="tar xvfz %s" % (tarname)
				self.execute(cmd, workdir=sumdir)

		# dictionary to have RFI fractions (extracted from the log files)
		rfis={}
		beamlogs=[log for log in rglob(sumdir, "*.log", 3)]
		for beamlog in beamlogs:
                        # excluding log-files in root directory (*_summary.log, *_pulp.log, pipeline.log, etc)
                        if beamlog.find("/SAP") < 0: continue
			# reading log file to get RFI fraction from rfifind
                        logf = open(beamlog, 'r')
      			rfifracs = [ff.split("(")[1].split("%")[0].lstrip() for ff in logf.read().splitlines() if re.search("Number of  bad", ff) is not None]
                        logf.close()
	  		cursapid=int(beamlog.split("/SAP")[-1].split("/")[0])
			curprocdir=beamlog.split("/SAP")[-1].split("/")[1]
                        if np.size(rfifracs) == 0: rfifrac="" # in case we did not run rfifind
                        elif np.size(rfifracs) > 1: rfifrac = rfifracs[-1]
                        else: rfifrac = rfifracs[0]
			rfis["%d_%s" % (cursapid, curprocdir)] = rfifrac

		# getting the list of *.pfd.bestprof files and read chi-sq values for all folded pulsars
		if not cmdline.opts.is_nofold:
			self.log.info("Reading chi-squared values and adding to chi-squared.txt...")
                        # also preparing montage command to create combined plot
	       	        montage_cmd="montage -background none -pointsize 10.2 "
	       	        montage_cmd_pdf="montage -geometry 100% -adjoin -tile 1x1 -pointsize 12 "
                	chif=open("%s/chi-squared.txt" % (sumdir), 'w')
     	       	        psr_bestprofs=sorted(rglob(sumdir, "*.pfd.bestprof", 3))
			if len(psr_bestprofs) > 0:
				# check first if all available *bestprof files are those created without mask. If so, then allow to make
				# a diagnostic combined plot using prepfold plots without a mask
				good_bestprofs=sorted([file for file in psr_bestprofs if re.search("_nomask_", file) is None])
				if len(good_bestprofs) == 0:
					good_bestprofs=[file for file in psr_bestprofs]
				else: # now we need to check that we have pfd files for all beams/stations
					# the reason to check is that for some of the beams we can only have pfd files with _nomask_
					# then len(good_bestprofs) != 0, but we miss these beams in good_bestprofs
					nomask_bestprofs=sorted([file for file in psr_bestprofs if re.search("_nomask_", file) is not None])
					if len(good_bestprofs) != len(nomask_bestprofs):
						good_procdirs=[bp.split(sumdir+"/")[-1].split("_SAP")[-1].split("_")[1] for bp in good_bestprofs]
						nomask_procdirs=[bp.split(sumdir+"/")[-1].split("_SAP")[-1].split("_")[1] for bp in nomask_bestprofs]
						only_nomask_procdirs=list(set(nomask_procdirs)-set(nomask_procdirs).intersection(set(good_procdirs)))
						only_nomask_bestprofs=[ff for ff in psr_bestprofs if ff.split(sumdir+"/")[-1].split("_SAP")[-1].split("_")[1] in only_nomask_procdirs]
						# adding bestprof files that were only made with _nomask_ to the list of good bestprofs
						good_bestprofs.extend(only_nomask_bestprofs)
					
                       		for bp in good_bestprofs:
               	               		psr=bp.split("/")[-1].split("_")[0]
                	        	thumbfile=bp.split(sumdir+"/")[-1].split(".pfd.bestprof")[0] + ".pfd.th.png"
					# we need it for combined.pdf
					bigfile=thumbfile.split(".th.png")[0] + ".png"

               	        		# getting current number for SAP and TA beam (or station name for FE)
					try: # we need this try block in case User manually creates sub-directories with some test bestprof files there
					     # which will also be found by rglob function. So, we need to exclude them by catching an Exception
					     # on a wrong-formed string applying int()
                      		  		cursapid=int(thumbfile.split("_SAP")[-1].split("_")[0])
					except: continue
               	                	curprocdir=thumbfile.split("_SAP")[-1].split("_")[1]

					# getting Stokes index if this is IQUV data
					stokes=""
					if is_iquv:
						stokes = bigfile.split("%s_S" % (curprocdir))[-1].split("_")[0]

                                	# calculating the peak S/N from the bestprof file
	                                snr=0.0
        	                        try:
                	                        prof = np.loadtxt(bp, comments='#', usecols=(1,1), dtype=float, unpack=True)[0]
                        	                rms=np.std(np.sort(prof)[0:int(np.size(prof)/2.)])
                                	        mean=np.mean(np.sort(prof)[0:int(np.size(prof)/2.)])
                                        	snr=(np.max(prof)-mean)/rms
	                                except:
                	                        self.log.warning("Warning: can't read file %s or calculate peak S/N of the profile" % (bp))

					# getting the part number if there are many splits
					curpart=""
					if thumbfile != thumbfile.split("_PART")[0]:
						curpart=thumbfile.split("_PART")[-1].split("_")[0]

					# check if key exists in a rfis dictionary
					if "%d_%s" % (cursapid, curprocdir) in rfis:
						rfifrac = rfis["%d_%s" % (cursapid, curprocdir)]
					else: rfifrac = ""
	                        	chif.write("file=%s obs=%s_SAP%d_%s%s%s_%s S/N=%g%s\n" % (thumbfile, data_code, \
								cursapid, curprocdir, stokes != "" and "_S%s" % (stokes) or "", curpart != "" and "_PART%s" % (curpart) or "", psr, snr, rfifrac != "" and " RFI=%s" % (rfifrac) or ""))
					if os.path.exists("%s/%s" % (sumdir, thumbfile)):
        		                	montage_cmd += "-label '%s SAP%d %s%s%s\n%s\nS/N = %g' %s " % (data_code, \
								cursapid, curprocdir, stokes != "" and "\nSTOKES_%s" % ("IQUV"[int(stokes)]) or "", \
								curpart != "" and "%sPART%s" % (stokes != "" and " " or "\n", curpart) or "", psr, snr, thumbfile)
					if os.path.exists("%s/%s" % (sumdir, bigfile)):
						montage_cmd_pdf += "-label '%s SAP%d %s%s%s\n%s\nS/N = %g' %s " % (data_code, \
								cursapid, curprocdir, stokes != "" and "\nSTOKES_%s" % ("IQUV"[int(stokes)]) or "", \
								curpart != "" and "%sPART%s" % (stokes != "" and " " or "\n", curpart) or "", psr, snr, bigfile)
              	        chif.close()

			# creating combined plots
			if len(psr_bestprofs) > 0:
                        	self.log.info("Combining all pfd.th.png files in a single combined plot...")
        	        	montage_cmd += "combined.png"
				self.execute(montage_cmd, workdir=sumdir)
				# making thumbnail version of the combined plot
				self.log.info("Making a thumbnail version of combined plot...")
               	        	cmd="convert -resize 200x140 -bordercolor none -border 150 -gravity center -crop 200x140-0-0 +repage combined.png combined.th.png"
               			self.execute(cmd, workdir=sumdir)

				# creating combined PDF plot with all prepfold plots - ONLY for FE 
				if data_code == "CS" and obs.FE:
                        		self.log.info("Combining all pfd.pdf files in a single combined multi-page PDF file...")
        	        		montage_cmd_pdf += "combined.pdf"
					self.execute(montage_cmd_pdf, workdir=sumdir)

		# create beam_process_node.txt file (this is only if file does not exist or it is empty)
		beam_process_node_file="%s/beam_process_node.txt" % (sumdir)
		if not os.path.exists(beam_process_node_file) or os.path.getsize(beam_process_node_file) == 0:
			self.log.info("Creating the beam_process_node.txt file...")
        	        bpnf=open(beam_process_node_file, 'w')
			for unit in workunits:
				for node in unit.tab.location:
					if node in unit.tab.rawfiles:
						for rf in unit.tab.rawfiles[node]:
							bpnf.write("%s %s%s\n" % \
							(node, rf, unit.tab.specificationType == "flyseye" and " [%s]" % (",".join(unit.tab.stationList)) or ""))
			bpnf.close()

		# removing old version of all status png files (if exist)
		self.log.info("Removing previous status png files (if any) ...")
		cmd="rm -f %s %s %s" % (" ".join(glob.glob("%s/*status*.png" % (sumdir))), " ".join(glob.glob("%s/ta_heatmap_sap*.png" % (sumdir))), " ".join(glob.glob("%s/rratmap_sap*.png" % (sumdir))))
		self.execute(cmd, workdir=sumdir)

		# creating TA heatmaps 
		# only when folding, and only if pulsars are set from the command line, or 'parset' or 'sapfind' or 'sapfind3' or "tabfind+" keywords are used (or
		# nothing is given for --pulsar option)
		# otherwise, different TA beams will be folded for different pulsars, and TA heatmap does not have sense
		if data_code == "CS" and not cmdline.opts.is_nofold and (len(cmdline.psrs) == 0 or (len(cmdline.psrs) != 0 and cmdline.psrs[0] != "tabfind")):
			for sap in obs.saps:
				# getting number of _coherent_ TABs
				nrTABs = len([1 for tab in sap.tabs if tab.is_coherent])
				if sap.nrRings > 0 or nrTABs > 1:
					if len(cmdline.psrs) != 0 and cmdline.psrs[0] != "parset" and cmdline.psrs[0] != "sapfind" and \
										cmdline.psrs[0] != "sapfind3" and cmdline.psrs[0] != "tabfind+":
						psrs = cmdline.psrs # getting list of pulsars from command line
					else: 
						if len(cmdline.psrs) == 0:
							if sap.source != "" and check_pulsars(sap.source, cmdline, cep2, None): psrs = [sap.source]
							else:
								if len(sap.psrs) > 0: psrs = [sap.psrs[0]]
								else: psrs = []
						else:
							if cmdline.psrs[0] == "parset":
								if sap.source != "" and check_pulsars(sap.source, cmdline, cep2, None): psrs = [sap.source]
								else: psrs = []
							if cmdline.psrs[0] == "sapfind" or cmdline.psrs[0] == "sapfind3":
								if len(sap.psrs) > 0:
									if cmdline.psrs[0] == "sapfind": psrs = [sap.psrs[0]]
									else: psrs = sap.psrs
								else: psrs = []
							if cmdline.psrs[0] == "tabfind+":
								if sap.source != "" and check_pulsars(sap.source, cmdline, cep2, None): psrs = [sap.source]
								else: psrs = []
								if len(sap.psrs) > 0: psrs += sap.psrs
								psrs = list(np.unique(psrs))
					if sap.nrRings > 0: self.log.info("Creating TA heatmap with %d rings for SAP=%d..." % (sap.nrRings, sap.sapid))
					else: self.log.info("Creating TA heatmap with %d TA beams for SAP=%d..." % (nrTABs, sap.sapid))
					for psr in psrs:
						self.log.info(psr)
						# I need this try/except block here, to avoid situation, when there are more than 1 pulsar in the SAP
						# but processing was done only for 1 - this is the case when pulsar is not specified in the command line
						try:
							cmd="cat %s/chi-squared.txt | grep _SAP%d | grep %s > %s/%s-chi-squared.txt" % (sumdir, sap.sapid, psr, sumdir, psr)
							self.execute(cmd, is_os=True)
							if not cmdline.opts.is_cobalt:	
								cmd="plot_LOFAR_TA_multibeam3.py --sap %d --chi %s-chi-squared.txt --parset %s.parset --out_logscale %s_SAP%d_%s_TA_heatmap_log.png --out_linscale %s_SAP%d_%s_TA_heatmap_linear.png --target %s" % (sap.sapid, psr, obs.id, obs.id, sap.sapid, psr, obs.id, sap.sapid, psr, psr)
							else: # for Cobalt we first need to make pseudo parset file
								obs.pseudo_parset_generator("%s/%s.pseudo.parset" % (sumdir, obs.id))
								cmd="plot_LOFAR_TA_multibeam3.py --sap %d --chi %s-chi-squared.txt --parset %s.pseudo.parset --out_logscale %s_SAP%d_%s_TA_heatmap_log.png --out_linscale %s_SAP%d_%s_TA_heatmap_linear.png --target %s" % (sap.sapid, psr, obs.id, obs.id, sap.sapid, psr, obs.id, sap.sapid, psr, psr)
							self.execute(cmd, workdir=sumdir)
							cmd="rm -f %s-chi-squared.txt" % (psr)
							self.execute(cmd, workdir=sumdir)
							# combining TA heatmap log and linear plots
							cmd="convert %s_SAP%d_%s_TA_heatmap_log.png %s_SAP%d_%s_TA_heatmap_linear.png -append ta_heatmap_sap%d_%s.png" % (obs.id, sap.sapid, psr, obs.id, sap.sapid, psr, sap.sapid, psr)
							self.execute(cmd, workdir=sumdir)
						except: 
							self.log.info("Can't make a heatmap plot for pulsar %s" % (psr))
						# if we do single-pulse analysis then we try to do RRATs heatmap as well
						if cmdline.opts.is_rrats or cmdline.opts.is_single_pulse:
							try:
								psr2=re.sub(r'^[BJ]', '', psr)
								parfile="%s/%s.par" % (sumdir, psr2)
								ref_unit=workunits[0]
								if cmdline.opts.is_nofold:
									psrdm = 0.0
								else:
									if not os.path.exists(parfile):
										parfiles=glob.glob("%s/*/SAP*/BEAM*/*.par" % (sumdir))
										if len(parfiles) > 0: 
											psrdm=ref_unit.get_psr_dm(parfiles[0])
										else: psrdm=0.0
									else: psrdm=ref_unit.get_psr_dm(parfile)
								if not cmdline.opts.is_cobalt:	
									cmd="RRAT_heatmap2.py --parset %s.parset --target %s --dm %f --sap %d --out_logscale %s_SAP%d_%s_RRAT_heatmap_log.png --out_linscale %s_SAP%d_%s_RRAT_heatmap_linear.png" % (obs.id, psr, psrdm, sap.sapid, obs.id, sap.sapid, psr, obs.id, sap.sapid, psr)
								else: # for Cobalt we first need to make pseudo parset file
									obs.pseudo_parset_generator("%s/%s.pseudo.parset" % (sumdir, obs.id))
									cmd="RRAT_heatmap2.py --parset %s.pseudo.parset --target %s --dm %f --sap %d --out_logscale %s_SAP%d_%s_RRAT_heatmap_log.png --out_linscale %s_SAP%d_%s_RRAT_heatmap_linear.png" % (obs.id, psr, psrdm, sap.sapid, obs.id, sap.sapid, psr, obs.id, sap.sapid, psr)
								self.execute(cmd, workdir=sumdir)
								# combining RRAT heatmap log and linear plots
								cmd="convert %s_SAP%d_%s_RRAT_heatmap_log.png %s_SAP%d_%s_RRAT_heatmap_linear.png -append %s_sap%d_rratmap.png" % (obs.id, sap.sapid, psr, obs.id, sap.sapid, psr, psr, sap.sapid)
								self.execute(cmd, workdir=sumdir)
							except:
								self.log.info("Can't make a RRAT heatmap plot for pulsar %s in SAP%d" % (psr, sap.sapid))

					# combining TA heatmaps for different pulsars
					heatmaps=glob.glob("%s/ta_heatmap_sap%d_*.png" % (sumdir, sap.sapid))
					if len(heatmaps) > 0:
						if len(heatmaps) > 1: cmd="convert %s +append ta_heatmap_sap%d.png" % (" ".join(heatmaps), sap.sapid)
						else: cmd="mv %s ta_heatmap_sap%d.png" % (heatmaps[0], sap.sapid)
						self.execute(cmd, workdir=sumdir)
						# remove temporary png files
						cmd="rm -f %s" % (" ".join(heatmaps))
						self.execute(cmd, workdir=sumdir)

					# combining RRAT heatmap for different pulsars
					if cmdline.opts.is_rrats or cmdline.opts.is_single_pulse:
						rratmaps=glob.glob("%s/*_sap%d_rratmap.png" % (sumdir, sap.sapid))
						if len(rratmaps) > 0:
							if len(rratmaps) > 1: cmd="convert %s +append rratmap_sap%d.png" % (" ".join(rratmaps), sap.sapid)
							else: cmd="mv %s rratmap_sap%d.png" % (rratmaps[0], sap.sapid)
							self.execute(cmd, workdir=sumdir)
							# remove temporary png files
							cmd="rm -f %s" % (" ".join(rratmaps))
							self.execute(cmd, workdir=sumdir)

			# combining TA heatmaps for different SAPs
			heatmaps=glob.glob("%s/ta_heatmap_sap*.png" % (sumdir))
			if len(heatmaps) > 0:
				if len(heatmaps) > 1: cmd="convert %s -append TAheatmap_status.png" % (" ".join(heatmaps))
				else: cmd="mv %s TAheatmap_status.png" % (heatmaps[0])
				self.execute(cmd, workdir=sumdir)
				# remove temporary png files
				cmd="rm -f %s" % (" ".join(heatmaps))
				self.execute(cmd, workdir=sumdir)
				# making a thumbnail version of TA heatmap combined plot
				cmd="convert -scale 200x140-0-0 TAheatmap_status.png TAheatmap_status.th.png"
				self.execute(cmd, workdir=sumdir)

			# combining RRAT heatmap for different SAPs
			if cmdline.opts.is_rrats or cmdline.opts.is_single_pulse:
				rratmaps=glob.glob("%s/rratmap_sap*.png" % (sumdir))
				if len(rratmaps) > 0:
					if len(rratmaps) > 1: cmd="convert %s -append rratmap.png" % (" ".join(rratmaps))
					else: cmd="mv %s rratmap.png" % (rratmaps[0])
					self.execute(cmd, workdir=sumdir)
					# remove temporary png files
					cmd="rm -f %s" % (" ".join(rratmaps))
					self.execute(cmd, workdir=sumdir)

		# creating combined DSPSR plots
		if not cmdline.opts.is_skip_dspsr:
			# first check if there are diagnostic plots after pdmp
			dspsr_diags=rglob(sumdir, "*_diag_pdmp.png", 3)
			# if not, then we check for diagnostic plots before pdmp
			if len(dspsr_diags) == 0:
				dspsr_diags=rglob(sumdir, "*_diag.png", 3)
			if len(dspsr_diags) > 0:
				self.log.info("Creating DSPSR summary diagnostic plots...")
				if len(dspsr_diags) > 1: cmd="montage %s -background none -mode concatenate -tile %dx dspsr_status.png" % (" ".join(sorted(dspsr_diags)), int(math.sqrt(len(dspsr_diags))))
				else: cmd="cp -f %s dspsr_status.png" % (dspsr_diags[0])
				self.execute(cmd, workdir=sumdir)
				# making a thumbnail version of combined DSPSR plot
				cmd="convert -scale 200x140-0-0 dspsr_status.png dspsr_status.th.png"
				self.execute(cmd, workdir=sumdir)

		# creating FE status maps
		# FE combined map if exist - should be called FE_status.png
		if data_code == "CS" and not cmdline.opts.is_nofold and obs.FE:
			self.log.info("Creating FE status maps...")
			try:
				cmd="lofar_status_map.py"
				self.execute(cmd, workdir=sumdir)
			except Exception:
				self.log.info("lofar_status_map.py failed")
			femaps=glob.glob("%s/*_core_status.png" % (sumdir))
			femaps.extend(glob.glob("%s/*_remote_status.png" % (sumdir)))
			femaps.extend(glob.glob("%s/*_intl_status.png" % (sumdir)))
			# creating combined plots
			if len(femaps) > 0:
				cmd="convert -append %s FE_status.png" % (" ".join(femaps))
				self.execute(cmd, workdir=sumdir)
				# removing individual maps
				cmd="rm -f %s" % (" ".join(femaps))
				self.execute(cmd, workdir=sumdir)
				# making a thumbnail version of FE status map
				cmd="convert -scale 200x140-0-0 FE_status.png FE_status.th.png"
				self.execute(cmd, workdir=sumdir)

		# Combining different status maps into one 'status.png' to be shown in web-summary page 
		# combining FE maps to general status.png
		if os.path.exists("%s/FE_status.png" % (sumdir)):
			self.log.info("Copying FE status map file to status.png ...")
			cmd="cp -f FE_status.png status.png"
			self.execute(cmd, workdir=sumdir)
			cmd="mv FE_status.th.png status.th.png"
			self.execute(cmd, workdir=sumdir)

		# combining TA heatmaps to general status.png
		if os.path.exists("%s/TAheatmap_status.png" % (sumdir)):
			if os.path.exists("%s/status.png" % (sumdir)): # means that FE maps were created
				self.log.info("Appending TA heatmap map file to status.png ...")
				cmd="montage status.png TAheatmap_status.png -background none -mode concatenate -tile 2x .temp_status.png"
				self.execute(cmd, workdir=sumdir)
				cmd="montage status.th.png TAheatmap_status.th.png -background none -mode concatenate -tile 2x .temp_status.th.png"
				self.execute(cmd, workdir=sumdir)
				cmd="mv .temp_status.png status.png"
				self.execute(cmd, workdir=sumdir)
				cmd="mv .temp_status.th.png status.th.png"
				self.execute(cmd, workdir=sumdir)
				cmd="rm -f TAheatmap_status.th.png"
				self.execute(cmd, workdir=sumdir)
			else:  # if status.png does not exist yet
				self.log.info("Copying TA heatmap map file to status.png ...")
				cmd="cp -f TAheatmap_status.png status.png"
				self.execute(cmd, workdir=sumdir)
				cmd="mv TAheatmap_status.th.png status.th.png"
				self.execute(cmd, workdir=sumdir)

		# combining dspsr summary plots to general status.png
		if os.path.exists("%s/dspsr_status.png" % (sumdir)):
			if os.path.exists("%s/status.png" % (sumdir)): # means that either FE maps or TA heatmap(s) were created
				self.log.info("Appending dspsr status file to status.png ...")
				cmd="montage status.png dspsr_status.png -background none -mode concatenate -tile 2x .temp_status.png"
				self.execute(cmd, workdir=sumdir)
				cmd="montage status.th.png dspsr_status.th.png -background none -mode concatenate -tile 2x .temp_status.th.png"
				self.execute(cmd, workdir=sumdir)
				cmd="mv .temp_status.png status.png"
				self.execute(cmd, workdir=sumdir)
				cmd="mv .temp_status.th.png status.th.png"
				self.execute(cmd, workdir=sumdir)
				cmd="rm -f dspsr_status.th.png"
				self.execute(cmd, workdir=sumdir)
			else:  # if status.png does not exist yet
				self.log.info("Copying dspsr status file to status.png ...")
				cmd="cp -f dspsr_status.png status.png"
				self.execute(cmd, workdir=sumdir)
				cmd="mv dspsr_status.th.png status.th.png"
				self.execute(cmd, workdir=sumdir)

		# creating thumbnail version of status.png if it exists
		if os.path.exists("%s/status.png" % (sumdir)):		
			self.log.info("Making a thumbnail version of status.png file...")
			cmd="convert -scale 200x140-0-0 status.th.png .temp_status.th.png"
			self.execute(cmd, workdir=sumdir)
			cmd="mv .temp_status.th.png status.th.png"
			self.execute(cmd, workdir=sumdir)
		else:
			self.log.info("No status.png created")

		# making complete tarballs for singlepulse *.inf and *.singlepulse files
		if cmdline.opts.is_rrats or cmdline.opts.is_single_pulse:
			self.log.info("Making tarballs for *.inf and *.singlepulse files...")
			# making tarball for *_DM*.inf files
			tarname="%s_sp_inf.tar" % (obs.id)
			tarlist=glob.glob("%s/*/SAP*/BEAM*/*_DM*.inf" % (sumdir))
			if len(tarlist) > 0:
				# make temporary file that lists all inf-files to be used as input for tar command 
				# this is to avoid very long argument list in the command line
				ff=open("%s/inflist" % (sumdir), 'w')
				ff.writelines("%s\n" % item.split(sumdir+"/")[1] for item in tarlist)
				ff.close()
				cmd="tar -cv --ignore-failed-read -f %s --files-from inflist" % (tarname)
				try: 
					self.execute(cmd, workdir=sumdir)
					# gzipping...
					cmd="gzip -f %s" % (tarname)
					self.execute(cmd, workdir=sumdir)
				except: pass
				# remove individual *_DM*.inf files
				cmd="find ./ -name '*_DM*.inf' -exec rm -f {} \;"  # 'find' command here to avoid very long argument list
				self.execute(cmd, workdir=sumdir)
				# remove inflist
				cmd="rm -f %s/inflist" % (sumdir)
				os.system(cmd)

			# making tarball for *_DM*.singlepulse files
			tarname="%s_sp_singlepulse.tar" % (obs.id)
			tarlist=glob.glob("%s/*/SAP*/BEAM*/*_DM*.singlepulse" % (sumdir))
			if len(tarlist) > 0:
				# make temporary file that lists all .singlepulse files to be used as input for tar command 
				# this is to avoid very long argument list in the command line
				ff=open("%s/splist" % (sumdir), 'w')
				ff.writelines("%s\n" % item.split(sumdir+"/")[1] for item in tarlist)
				ff.close()
				cmd="tar -cv --ignore-failed-read -f %s --files-from splist" % (tarname)
				try: 
					self.execute(cmd, workdir=sumdir)
					# gzipping...
					cmd="gzip -f %s" % (tarname)
					self.execute(cmd, workdir=sumdir)
				except: pass
				# remove individual *_DM*.singlepulse files
				cmd="find ./ -name '*_DM*.singlepulse' -exec rm -f {} \;"  # 'find' command here to avoid very long argument list
				self.execute(cmd, workdir=sumdir)
				# remove splist
				cmd="rm -f %s/splist" % (sumdir)
				os.system(cmd)
			
		# Make a tarball of all the plots (summary archive)
		if is_to_combine:
			extensions = self.summary_archive_exts
		else:
			extensions = self.summary_archive_exts_nocombine
		self.log.info("Making a final summary tarball of all files with extensions: %s" % (", ".join(extensions)))
		tarname="%s%s%s%s" % (cep2.pipeid, self.summary_archive_prefix, data_code, self.summary_archive_suffix)
		tar_list=[]
		for ext in extensions:
			ext_list=rglob(sumdir, ext, 3)
			tar_list.extend(ext_list)
		cmd="tar -cv --ignore-failed-read -f %s %s" % (tarname, " ".join([f.split(sumdir+"/")[1] for f in tar_list]))
		try: # --ignore-failed-read does not seem to help with tar failing for some beams
                     # like file was changed during the tar, though tarball seem to be fine
			self.execute(cmd, workdir=sumdir)
		except: pass

		# finish
		end_time=time.time()
		total_time= end_time- start_time
		self.log.info("UTC stop time is: %s" % (time.asctime(time.gmtime())))
		self.log.info("Total running time: %.1f s (%.2f hrs)" % (total_time, total_time/3600.))

		# flushing log file and copy it to summary node
		self.log.flush()
		if not cmdline.opts.is_log_append: cmd="cp -f %s %s" % (cep2.get_logfile(), sumdir)
		else: cmd="cat %s >> %s/%s" % (cep2.get_logfile(), sumdir, cep2.get_logfile().split("/")[-1])
		os.system(cmd)

		# delete file from the archive first
		cmd="tar --delete -v --ignore-failed-read -f %s %s" % (tarname, cep2.get_logfile().split("/")[-1])
		try:
			proc = Popen(shlex.split(cmd), stdout=PIPE, stderr=STDOUT, cwd=sumdir)
			proc.communicate()
			# adding log file to the archive and gzip it
			cmd="tar -rv --ignore-failed-read -f %s %s" % (tarname, cep2.get_logfile().split("/")[-1])
			# --ignore-failed-read does not seem to help with tar failing for some beams
			# like file was changed during the tar, though tarball seem to be fine
			proc = Popen(shlex.split(cmd), stdout=PIPE, stderr=STDOUT, cwd=sumdir)
			proc.communicate()
		except: pass
		# avoid gzipping now
		#cmd="gzip -f %s" % (tarname)
		#proc = Popen(shlex.split(cmd), stdout=PIPE, stderr=STDOUT, cwd=sumdir)
		#proc.communicate()

		# updating the Feedback unit
		self.make_feedback(obs, cep2, cmdline)

		# specific to dragnet
		if cep2.cluster_headnode == "dragnet":
			# changing ownership to 'dragnet' group
			cmd="chgrp -R dragnet %s" % (sumdir)
			os.system(cmd)
		# changing the file permissions to be re-writable for group
		cmd="chmod -R g+w %s" % (sumdir)
		os.system(cmd)
