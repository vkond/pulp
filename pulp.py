#!/usr/bin/env python
#
# This is Python implementation of PULP, Pulsar standard pipeline
# written by A2 in shell
#
# Vlad, Oct 25, 2011 (c)
#
# Major update - Vlad, Mar-Apr, 2014
#  - upgrades related to Cobalt
#  - upgrades related to incorporate to Automatic framework
###################################################################
import os, sys, glob, time
import getopt
import numpy as np
import cPickle
import re
import logging
import signal
import subprocess, shlex
from subprocess import PIPE, STDOUT, Popen
from pulp_parset import Observation
from pulp_usercmd import CMDLine
from pulp_sysinfo import CEP2Info
from pulp_logging import PulpLogger
from pulp_pipeline import Pipeline
from pulp_feedback import FeedbackUnit

# exit function that cleans terminal before exiting
# after Ctrl-C and when using "ssh -t" terminal gets messed up, so one has to reset it
# the command "stty sane" allows to reset terminal without clearing it (it puts all esc sequences
# to its default values)
def quit (status):
	Popen(shlex.split("stty sane"), stderr=open(os.devnull, 'rb')).wait()
	sys.exit(status)

# Main function
def main(wrapper=None):

	# Check first that all Python modules that we need are present in the system
	try:
		modules_to_check=['os', 'sys', 'glob', 'time', 'getopt', 'numpy', 'cPickle', 're', 'logging', 'signal', 'subprocess', 'shlex',
                                  'h5py', 'tarfile', 'optparse', 'random', 'string', 'math', 'psr_utils',
                                  'pulp_logging', 'pulp_feedback', 'pulp_sysinfo', 'pulp_usercmd', 'pulp_parset',
                                  'pulp_pipeunit', 'pulp_pipeunitpart', 'pulp_pipeline', 'pulp_cvunitpart', 'pulp_cvunit']
		for module in modules_to_check:
			__import__('imp').find_module(module)
	except ImportError as ex:
		print "One or more Python modules are missing: [%s]" % (", ".join(ex.args))
		raise

#	# first define signal handler for SIGTERM (if system or user sends kill signal), "kill -9" can _not_ be caught 
	def sigterm_handler(signum, frame):
		raise Exception("SIGTERM signal got caught!")
	signal.signal(signal.SIGTERM, sigterm_handler)

        # parsing command line
	cmdline = CMDLine("", wrapper)

	if not cmdline.opts.is_auto:
		# getting the version of the pulp.py (svn revision number)
		try:
			cmd="svn info"
			proc = Popen(shlex.split(cmd), stdout=PIPE, stderr=STDOUT, cwd=os.environ["LOFARSOFT"] + "/src/Pulsar/scripts/pypulp")
			(out, err) = proc.communicate()
			VERSION="r" + out.split("Revision: ")[1].split("\n")[0]
		except Exception:
			VERSION="unknown version"
		cmdline.set_version(VERSION)

	if not cmdline.opts.is_auto or cmdline.opts.is_local:
		# creating name of the Logger
		logger_name = "PULP"
		# start forming name of the logfile
		logfile = cmdline.opts.obsid
		# set local variable obsid
		obsid = cmdline.opts.obsid

	if cmdline.opts.is_local and cmdline.opts.beam_str != "":
		if not cmdline.opts.is_summary:
			beam=cmdline.opts.beam_str.split(",")[0]
                        if re.search(r'[^\:\d\/]+', beam) is not None:
                                print "Option --beams can only has digits, colons, commas, and in some cases / for parts!\nCan't start local processing. Exiting..."
                                return 1
                        if re.search(r'[\:]+', beam) is None:
                                print "Option --beams should have at least one colon!\nCan't start local processing. Exiting..."
                                return 1
			(sap, tab) = beam.split(":")
			if sap == "" or tab == "":
                                print "Option --beams has at least one empty SAP or TAB value!\nCan't start local processing. Exiting..."
                                return 1
			# getting SAP and TAB ids for local processing
			sapid = int(sap)
			tabid = int(tab.split("/")[0])
			partid = -1
			stokesid = cmdline.opts.stokes

			logger_name += beam
			# if there is no "part" part, then form the name of the log-file in a usual way
			if re.search(r'[\/]+', tab) is None:
				# forming finally the name of the local log-file
				if stokesid == -1:
					logfile = "%s_sap%03d_beam%04d.log" % (obsid, sapid, tabid)
				else:
					logfile = "%s_sap%03d_beam%04d_stokes%d.log" % (obsid, sapid, tabid, stokesid)
			else: # if there are several BW splits for one beam, then we also add "part" number in the name of the log-file
				part = tab.split("/")[-1]
				if part == "":
                                	print "Option --beams has empty PART value!\nCan't start local processing. Exiting..."
                                	return 1
				partid = int(part)
				if stokesid == -1:
					logfile = "%s_sap%03d_beam%04d_part%02d.log" % (obsid, sapid, tabid, partid)
				else:
					logfile = "%s_sap%03d_beam%04d_stokes%d_part%02d.log" % (obsid, sapid, tabid, stokesid, partid)
		else:   
			# when --summary then name of the type of summary should be given in --beams (CS, IS or CV)
			# when run locally
			logger_name += "_summary%s" % (cmdline.opts.beam_str)
			# in this case summary locus node is given in --beams option
			logfile += "_summary%s.log" % (cmdline.opts.beam_str)
	else:
		if not cmdline.opts.is_auto: 
			logfile += "_pulp.log"

	# initializing the Logger
	if not cmdline.opts.is_auto or cmdline.opts.is_local:
		log = PulpLogger(logger_name)
	else:
		try:
			log = wrapper.logger
		except Exception as ex:
			print "Pulsar auto-pipeline wrapper object passed to the pulp is bad: [%s]" % (", ".join(ex.args))
			raise

	try:
		# in AUTO mode first get PipeID, ObsID and locations/filenames of input/output data
		if cmdline.opts.is_auto and not cmdline.opts.is_local:
			try:	
				if len(wrapper.input_data['coherent'].data) != 0:
					infile=wrapper.input_data['coherent'].data[0].file
				else:
					if len(wrapper.input_data['incoherent'].data) != 0:
						infile=wrapper.input_data['incoherent'].data[0].file
					else:
						log.error("No input files are given in the wrapper object!")
						return 1
				# ObsID
				obsid=infile.split("/")[-1].split("_")[0]

				if len(wrapper.output_data['data'].data) != 0:
					outfile=wrapper.output_data['data'].data[0].file
				else:
                                        log.error("No output files are given in the wrapper object!")
                                        return 1
				# PipeID
				# this is done specifically to handle different data structure on CEP4
				cluster_headnode=os.popen('hostname').readlines()[0].strip().split(".")[0]
				if cluster_headnode[:5] == "head0" or cluster_headnode[:3] == "cpu":
					pipeid=outfile.split("/pulp")[0].split("/")[-1]
					processed_dir=outfile.split("/pulp")[0]
				else: # for CEP2
					pipeid=outfile.split("/")[-1].split("_")[0]
					processed_dir="/".join(outfile.split("/")[:-1])
			except Exception as ex:
				log.error("Error accessing wrapper mapfile attributes: [%s]" % (", ".join(ex.args)))
				raise

		sysinfo_file = "sysinfo.b"
		if not cmdline.opts.is_noinit:
			# getting info about the CEP2
			cep2 = CEP2Info()
			if not cmdline.opts.is_auto:
		        	# forming the name of the log-file as OBSID_pulp.log in $HOME/.pulp/<OBSID> dir
		        	# if obsid is not given than the log-file will be _pulp.log
				cluster_headnode=os.popen('hostname').readlines()[0].strip().split(".")[0]
				if cluster_headnode[:5] == "head0" or cluster_headnode[:3] == "cpu":	
					cep2.home="/data/LOFAR_PULSAR_ARCHIVE"
				cep2.set_logdir("%s/.pulp/%s" % (cep2.home, obsid))
				cep2.set_pipeid(obsid)
                                cep2.set_feedbackfile("%s/.pulp/Observation%s_feedback" % (cep2.home, cep2.get_pipeid()[1:])) # assume feedback filename and exclude leading "L" from the ObsID
			else:
				cep2.logdir = wrapper.job_dir
				cep2.set_pipeid(pipeid)
				cep2.set_processed_dir(processed_dir)
				cluster_headnode=os.popen('hostname').readlines()[0].strip().split(".")[0]
				if cluster_headnode[:5] == "head0" or cluster_headnode[:3] == "cpu":	
					cep2.processed_dir_prefix = "/".join(processed_dir.split("/data")[-1].split("/")[:-1])
				cep2.set_feedbackfile(wrapper.parset_feedback_file) # use the feedback filename as supplied by the wrapper
			# checking if we are using SLURM. If so, then setting JobID and Jobname
			if cmdline.opts.is_slurm:
				try:
					cep2.set_slurm_jobid()
					cep2.set_slurm_jobname()
					if cep2.slurm_jobid == "" or cep2.slurm_jobname == "":
						raise
				except Exception as ex:
					log.error("Error in getting SLURM_JOB_ID and/or SLURM_JOB_NAME! Either not use --slurm with Pulp, or do resource allocation with salloc: [%s]" % (", ".join(ex.args)))
					raise
				
			# saving sysinfo configuration to file
                	sysfd = open (cep2.get_logdir() + "/" + sysinfo_file, "wb")
	                cPickle.dump(cep2, sysfd, True)  # when False, the dumpfile is ascii
        	        sysfd.close()
		else:
			confdir=cmdline.get_confdir()
			# loading sysinfo setup from the file
			sysfd = open(confdir + "/" + sysinfo_file, "rb")
			cep2=cPickle.load(sysfd)
			sysfd.close()

		if not cmdline.opts.is_auto or cmdline.opts.is_local:
			if cmdline.opts.is_auto:
				cep2.set_local_logdir(cep2.get_pipeid()[1:])
			cep2.set_logfile("%s/%s" % (cep2.get_local_logdir(cmdline), logfile))
			# adding file handler to our Logger
			logfh = logging.FileHandler(cep2.get_logfile(), mode='%s' % (cmdline.opts.is_log_append and "a" or "w"))
			log.addHandler(logfh)
		else:
			cep2.set_logfile(wrapper.config.get("logging", "log_file"))
			# cleaning potential left-over feedback files...
			cmd="rm -f %s/.*.fb" % (cep2.get_logdir())
			os.system(cmd)

		# setting current node and dir
		cep2.set_current_node()
		cep2.set_current_dir()

		# starting logging...
		start_pipe_time=time.time()
		log.info("UTC time is:  %s" % (time.asctime(time.gmtime())))

		obsconf_file = cep2.get_logdir() + "/obs.b"
		pipeline_file = cep2.get_logdir() + "/pipeline.b"

		# checking validity of the options
		cmdline.check_options(cep2, log)

		if not cmdline.opts.is_noinit and not cmdline.opts.is_auto:
			# checking connections to processig nodes
			# we do this in all cases except when we use Slurm _and_ GlobalFS
			if not (cmdline.opts.is_slurm and cmdline.opts.is_globalfs):
				cep2.check_connection(log)
				# print down nodes
				cep2.print_down_nodes(log)

		log.info("\nInitializing...")
		if not cmdline.opts.is_noinit:
			# initializing our Observation, collecting info from parset, etc.
			obs = Observation(obsid, cep2, cmdline, log, wrapper)

			# checking if rawdata available on the cluster for user-specified beams
			cmdline.is_rawdata_available(cep2, obs, log)

			# saving obs configuration to file
                	obsfd = open (obsconf_file, "wb")
	                cPickle.dump(obs, obsfd, True)  # when False, the dumpfile is ascii
        	        obsfd.close()
		else:
			# loading obs setup from the file
			obsfd = open(obsconf_file, "rb")
			obs=cPickle.load(obsfd)
			obsfd.close()

		# updating cmdline default parameters based on obtained info about Observation
		# such as, number of frequency splits
		# and about FWHMs of CS and IS beams (depends on what observation is, HBA or LBA)
		cmdline.update_default_values(obs, cep2, log)
		# printing info about observation
		obs.print_info(cmdline, log)

		if not cmdline.opts.is_local: 
			# print summary
			cmdline.print_summary(cep2, obs, log)

		# printing info both to STDOUT and logfile
		cep2.print_info(cmdline, log)

		# if --beam option is not set, it means that we start the pipeline from main node
		if not cmdline.opts.is_local:	
			# initializing pulsar pipeline
			psrpipe = Pipeline(obs, cep2, cmdline, log)
			# saving pipeline config to file
                	pipefd = open (pipeline_file, "wb")
              		cPickle.dump(psrpipe, pipefd, True)
			pipefd.close()
			# kick off the pipeline
			if not cmdline.opts.is_summary:
				psrpipe.start(obs, cep2, cmdline, log)
			# wait for all childs to finish and prepare logs, all files in order
			# convert, FE maps, etc.
			try:
				psrpipe.finish(obs, cep2, cmdline, log)
			except KeyboardInterrupt as ex:
				log.exception("User interruption...")
				if psrpipe != None:
					psrpipe.kill(log)
				raise
			except Exception as ex:
				if psrpipe != None:
					psrpipe.kill(log)
				raise
#			except SystemExit as ex:
#				log.exception("Caught SystemExit...")
#				if psrpipe != None:
#					psrpipe.kill(log)
#				raise

			# get failures
			nfailures = psrpipe.get_number_failed_pipes()
			nfailures += psrpipe.get_number_failed_summaries()

			# making feedback file (always even if there are failures)
			log.info("Making final feedback file %s..." % (cep2.get_feedbackfile()))
			# open feedback file
			fbunit = FeedbackUnit(-1, cep2.current_node, "/".join(cep2.get_logfile().split("/")[:-1]))
               		fbunit.update(cep2.get_logfile(), cep2.get_feedbackfile(), "Log", None, True)
               		fbunit.flush(100, cep2, False, True)
			# get list of all other feedback files, open them, append to the main feedback file and delete them
			fbunits=glob.glob(cep2.get_logdir() + "/" + ".%s*.fb" % (cep2.pipeid))
			for fb in fbunits:
				cmd="cat %s >> %s" % (fb, cep2.get_feedbackfile())
				os.system(cmd)
				cmd="rm -f %s" % (fb)
				os.system(cmd)
			# appending info about the number of data products to the feedback file
	                fb=open(cep2.get_feedbackfile(), 'a')
               		fb.write("%s = %d\n" % (cep2.feedback_nrofoutputs_keyword, len(fbunits)))
               		fb.close()
			# specific to dragnet
			if cep2.cluster_headnode == "dragnet":
				# changing group ownership to dragnet
				cmd="chgrp dragnet %s" % (cep2.get_feedbackfile())
				os.system(cmd)
				# making feedback file writabel to group
				cmd="chmod g+w %s" % (cep2.get_feedbackfile())
				os.system(cmd)

                        # end of the pipeline...
                        end_pipe_time=time.time()
                        pipe_total_time = end_pipe_time - start_pipe_time
                        log.info("Finished")
                        log.info("UTC time is:  %s" % (time.asctime(time.gmtime())))
                        log.info("Total wall time:  %.1f s (%.2f hrs)" % (pipe_total_time, pipe_total_time/3600.))
                        
                        # flushing log file and copy it to summary nodes
                        if not cmdline.opts.is_auto:
                                log.flush()
			try:
				for (sumcode, sumdir) in psrpipe.summary_dirs.items():
					sumnode = cep2.summary_nodes[sumcode]
					if not cmdline.opts.is_globalfs:
						cmd="rsync -avxP --rsh='ssh -x' %s %s:%s" % (cep2.get_logfile(), sumnode, sumdir)
					else:
						cmd="cp %s %s" % (cep2.get_logfile(), sumdir)
					proc = Popen(shlex.split(cmd), stdout=PIPE, stderr=STDOUT)
                                       	proc.communicate()
                        except Exception:
       	                        log.error("Copying of main log-file '%s' to %s:%s failed" % (cep2.get_logfile(), sumnode, sumdir))

			# raise Exception if there were failures in the beams processing or summaries processing
			if nfailures != 0:
				log.error("One or more beams/summaries have failed!")
				raise Exception
                        			
		else:
			# loading pipeline config from the file
			if os.path.exists(pipeline_file):
				pipefd = open(pipeline_file, "rb")
				psrpipe=cPickle.load(pipefd)
				pipefd.close()
			else: psrpipe = Pipeline(obs, cep2, cmdline, log)
			if not cmdline.opts.is_summary:
				# running processing for particular beam
				for unit in psrpipe.units:
					try:
						if partid != -1:
							if unit.sapid == sapid and unit.tabid == tabid and unit.part == partid:
								if stokesid == -1 or (stokesid != -1 and unit.stokes_index == stokesid):
									unit.run(obs, cep2, cmdline, log)
						else:
							if unit.sapid == sapid and unit.tabid == tabid:
								if stokesid == -1 or (stokesid != -1 and unit.stokes_index == stokesid):
									unit.run(obs, cep2, cmdline, log)
					except KeyboardInterrupt as ex:
						log.exception("User interruption...")
						unit.kill() # killing all open processes
						raise
					except Exception as ex:
						unit.kill() # killing all open processes
						raise
#					except SystemExit as ex:
#						log.exception("Caught SystemExit...")
#						unit.kill() # killing all open processes
#						raise
					
			else:   # running local pulp to make summary actions
				try:
					psrpipe.make_summary(obs, cep2, cmdline, log)
				except KeyboardInterrupt as ex:
					log.exception("User interruption...")
					if psrpipe != None:
						psrpipe.kill(log)
					raise
				except Exception as ex:
					if psrpipe != None:
						psrpipe.kill(log)
					raise
#				except SystemExit as ex:
#					log.exception("Caught SystemExit...")
#					if psrpipe != None:
#						psrpipe.kill(log)
#					raise

	except Exception as ex:
		log.exception("Oops... pulp has crashed! [%s]" % (", ".join(ex.args)))
		raise

	if not cmdline.opts.is_auto:
                log.flush()

	if not cmdline.opts.is_auto or cmdline.opts.is_local:
		logfh.close()
		log.removeHandler(logfh)

	# removing log file from ~/.pulp/obsid dir if it is for pulp on local node
	if cmdline.opts.is_local:
		cmd="rm -rf %s" % (cep2.get_logfile())
		os.system(cmd)

	if not cmdline.opts.is_auto or cmdline.opts.is_local:
		logging.shutdown()

	# At this point there were no errors, so we return 0
	return 0

###  M A I N ###
if __name__ == "__main__":
	try:
		status = main()
		if status != 0: raise Exception
		else: quit(0)
	except Exception as ex:
		print "PULP, Standalone mode, error in main(): [%s]" % (", ".join(ex.args))
		quit(1)
	except KeyboardInterrupt as ex:
		print "User interruption in Standalone mode: [%s]"  % (", ".join(ex.args))
		quit(1)
#	except SystemExit as ex:
#		print "SystemExit in Standalone mode: [%s]"  % (", ".join(ex.args))
#		quit(1)


# pulp Class for the Automatic framework
class pulp:
	"""
	The main class of the pulsar pipeline
	"""

	def __init__(self, wrapper):
		self.wrapper = wrapper

	def go(self):
		return main(self.wrapper)
