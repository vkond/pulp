###################################################################
#
# Class Observation description similar to obsinfo class used to gather
# info about observation
# (rewritten to handle only new format of parset files)
# Additional classes SAPBeam and TABeam are defined
#

import os, sys, glob
import numpy as np
import h5py
import time
import re
import math
import subprocess
from pulp_logging import PulpLogger

# Global function that calculates the distance between two points on sphere
def radial_distance(rarad1, decrad1, rarad2, decrad2):
	"""
	radial_distance: prog calculates radial distance between
	given (rarad1, decrad1) and (rarad2, decrad2) position (in radians). 
	It uses the formula used in ATNF pulsar catalogue to calculate distances 
	between the pulsars (Haversine formula).
	Return value: distance (in deg)
	"""
	dist = 2.*math.asin(math.sqrt((math.sin((decrad1-decrad2)/2.))**2 + math.cos(decrad1)*math.cos(decrad2)*(math.sin((rarad1-rarad2)/2.))**2))
	dist = (180. * dist) / math.pi
	return dist

# find all pulsars within max distance from selected center, sort them by flux
# and return 3 brightest ones
# coordinates are in radians
# max distance in degrees
def find_pulsars(rarad, decrad, cmdline, max_distance):
	radeg = (180. * rarad) / math.pi
	decdeg = (180. * decrad) / math.pi
	psrras=np.abs(cmdline.ras - radeg)
	psrdecs=np.abs(cmdline.decs - decdeg)
	# first, select only pulsars _roughly_ within required distance
	crit=(psrras<1.5*max_distance)&(psrdecs<1.5*max_distance)
	psrbs=cmdline.catpsrs[crit]
	psrras=cmdline.ras[crit]
	psrdecs=cmdline.decs[crit]
	psrs400=cmdline.s400[crit]
	# for selected pulsars calculate the distance precisely
	psrdist=np.array([radial_distance(rarad, decrad, (psrras[ii]/180.)*math.pi, (psrdecs[ii]/180.)*math.pi) for ii in np.arange(np.size(psrras))])
	psrdist-=max_distance
	crit=(psrdist<0)
	psrs400=psrs400[crit]
	psrbs=psrbs[crit]
	reo=re.compile(r'^\*$')
	psrs400=[float(reo.sub('0.0', s)) for s in psrs400]
	# sort by flux in reverse order (getting the reversed list of indices)
	ind=[ii for ii in reversed(np.argsort(psrs400))]
	psrbs=psrbs[ind]
	if np.size(psrbs) > 3: psrbs=psrbs[:3]
	return psrbs

# Class SAPBeam or "Station" beam describes the parameters specific for this particular 
# station beam
class SAPBeam:
	# si - system info object
	# root - "parent" class of the sap beam
	def __init__(self, id):
		self.sapid = id            # ID number of the SAP beam

		self.rarad=self.decrad=0
		self.source=""

		self.nrTiedArrayBeams = 0  # number of TA beams including possible IS beam
		self.tabs=[]               # list of TA Beams (including possible IS beam)

		self.nrRings = 0           # number of TA rings (if used)
		self.ringSize = 0          # size of TA ring (in deg)
		self.nrSubbands = 0        # number of subbands (can be different for different SAPs)
		self.subbandList=""        # range of subbands, e.g. 77..320
		self.subbands=[]           # list of subbands
		self.psrs = []             # up to 3 brightest pulsars in the SAP

	# updating attributes using parset file
	def update(self, plines, si, cmdline, root, log=None, wrapper=None):

		# Getting the list of subbands
		res=[ii for ii in plines if "Observation.Beam[%d].subbandList" % (self.sapid) in ii]
		try:
			# getting range of subbands
			self.subbandList=res[0].split("=", 1)[-1].strip().lstrip("[").rstrip("]")
			# parsing the string with subband ranges to get list of subbands
			self.subbands = root.getSubbands(self.subbandList)
			# getting total number of Subbands
			self.nrSubbands = len(self.subbands)
		except: pass

	        # getting info about the pointing
		res=[ii for ii in plines if "Observation.Beam[%d].angle1" % (self.sapid) in ii and "AnaBeam" not in ii]
		try:
			self.rarad=np.float64(res[0].split("=", 1)[-1].strip())
		except: pass

		res=[ii for ii in plines if "Observation.Beam[%d].angle2" % (self.sapid) in ii and  "AnaBeam" not in ii]
		try:
               		self.decrad=np.float64(res[0].split("=", 1)[-1].strip())
		except: pass

	        # getting info about Source name
		res=[ii for ii in plines if "Observation.Beam[%d].target" % (self.sapid) in ii]
		try:
                	self.source=res[0].split("=", 1)[-1].strip().strip("'\"").replace(' ', '_')
		except: pass

		# Getting number of TA Beams in station beam including possible IS beam
		res=[ii for ii in plines if "Observation.Beam[%d].nrTiedArrayBeams" % (self.sapid) in ii]
		try:
			self.nrTiedArrayBeams=int(res[0].split("=", 1)[-1].strip())
		except: pass

		# Getting number of TA rings
		res=[ii for ii in plines if "Observation.Beam[%d].nrTabRings" % (self.sapid) in ii]
		try:
			self.nrRings=int(res[0].split("=", 1)[-1].strip())
		except: pass

		# Getting the size of the TA ring (in deg)
		if self.nrRings != 0:
			res=[ii for ii in plines if "Observation.Beam[%d].tabRingSize" % (self.sapid) in ii]
			try:
				self.ringSize=np.float64(res[0].split("=", 1)[-1].strip())
				self.ringSize = self.ringSize * (180./math.pi)
			except: pass

		# find up to 3 brightest pulsars in the SAP
		if not cmdline.opts.is_nofold:
			if len(cmdline.psrs) == 0 or (len(cmdline.psrs) != 0 and cmdline.psrs[0] == "sapfind") or \
				(len(cmdline.psrs) != 0 and cmdline.psrs[0] == "sapfind3") or \
				(len(cmdline.psrs) != 0 and cmdline.psrs[0] == "tabfind+"):
				self.psrs = find_pulsars(self.rarad, self.decrad, cmdline, cmdline.opts.fwhm_IS/2.)
				if len(cmdline.psrs) != 0 and (cmdline.psrs[0] == "sapfind" or cmdline.psrs[0] == "tabfind+") and len(self.psrs) > 0: 
					self.psrs = self.psrs[:1]

		# initializing SAP beams objects and making the list of SAP beams
		for tid in xrange(self.nrTiedArrayBeams):
			tab=TABeam(tid, self.sapid, self.nrSubbands)
			tab.update(plines, si, cmdline, root, self.rarad, self.decrad, log, wrapper)
			self.tabs.append(tab)

	# updating attributes using .h5 files
	def update_from_H5(self, si, cmdline, root, log=None, wrapper=None):
		# get list of .h5 files for the current SAP
		saph5list = glob.glob("%s/%s_SAP%03d_B*_S0_P000_bf.h5" % (si.get_logdir(), root.id, self.sapid))
		# first, we read any 1 file to collect general info from the SUB_ARRAY_POINTING group
		s5 = h5py.File(saph5list[0], 'r')
		self.rarad = np.float64((s5['SUB_ARRAY_POINTING_%03d' % (self.sapid)].attrs['POINT_RA'] * math.pi) / 180.)
		self.decrad = np.float64((s5['SUB_ARRAY_POINTING_%03d' % (self.sapid)].attrs['POINT_DEC'] * math.pi) / 180.)
		self.nrTiedArrayBeams = s5['SUB_ARRAY_POINTING_%03d' % (self.sapid)].attrs['OBSERVATION_NOF_BEAMS']
		s5.close()
		# getting the target at the center of the SAP (or at least closest to the center)
		srcmap={}
		for ff in saph5list:
			s5 = h5py.File(ff, 'r')
			beams=[int(i.split("BEAM_")[1]) for i in s5['SUB_ARRAY_POINTING_%03d' % (self.sapid)].keys() if 'BEAM' in i]
			for b in beams:
				raoff=s5['SUB_ARRAY_POINTING_%03d' % (self.sapid)]['BEAM_%03d' % (b)].attrs['POINT_OFFSET_RA']
				decoff=s5['SUB_ARRAY_POINTING_%03d' % (self.sapid)]['BEAM_%03d' % (b)].attrs['POINT_OFFSET_DEC']
				offset=math.sqrt(raoff*raoff+decoff*decoff)
				srcmap[(ff,b)]=offset
			s5.close()
		(centerfile, centerbeam)=min(srcmap, key=lambda k: srcmap[k])
		s5 = h5py.File(centerfile, 'r')
		self.source = s5['SUB_ARRAY_POINTING_%03d' % (self.sapid)]['BEAM_%03d' % (centerbeam)].attrs['TARGETS'][0].replace(' ', '_')
		# also getting the Chan/sub and subband_width for this beam
		nchanpersub = s5['SUB_ARRAY_POINTING_%03d' % (self.sapid)]['BEAM_%03d' % (centerbeam)].attrs['CHANNELS_PER_SUBBAND']
		subwidth = s5['SUB_ARRAY_POINTING_%03d' % (self.sapid)]['BEAM_%03d' % (centerbeam)].attrs['SUBBAND_WIDTH'] # in Hz
		s5.close()
		# getting the total number of subbands in this SAP
		# and also getting the list of freqs of the channels!
		cfreqs=[]
		saph5list = glob.glob("%s/%s_SAP%03d_B%03d_S0_P*_bf.h5" % (si.get_logdir(), root.id, self.sapid, centerbeam)) # list of all parts for one of the beams
		for ff in saph5list:
			coord = 1
			s5 = h5py.File(ff, 'r')
			self.nrSubbands += s5['SUB_ARRAY_POINTING_%03d' % (self.sapid)]['BEAM_%03d' % (centerbeam)]['STOKES_0'].attrs['NOF_SUBBANDS']
			if s5['SUB_ARRAY_POINTING_%03d' % (self.sapid)]['BEAM_%03d' % (centerbeam)]['COORDINATES']['COORDINATE_0'].attrs['AXIS_NAMES'][0] == "Frequency": coord = 0 
			cfreqs.extend(list(s5['SUB_ARRAY_POINTING_%03d' % (self.sapid)]['BEAM_%03d' % (centerbeam)]['COORDINATES']['COORDINATE_%d' % (coord)].attrs['AXIS_VALUES_WORLD']))
			s5.close()
		cfreqs.sort() # sorting, min freq is first
		if nchanpersub > 1:
			cfreqs=np.array(cfreqs)
			cfreqs.shape=(np.size(cfreqs)/nchanpersub, nchanpersub)
			cfreqs=np.mean(cfreqs, axis=1)
		self.subbands=[int((fr - root.lower_band_edge*1e6 + 0.5*subwidth)/subwidth) for fr in cfreqs] # we add "0.5*subwidth" to avoid problems with rounding
		self.subbandList="%d..%d" % (self.subbands[0], self.subbands[-1])

		# find up to 3 brightest pulsars in the SAP
		if not cmdline.opts.is_nofold:
			if len(cmdline.psrs) == 0 or (len(cmdline.psrs) != 0 and cmdline.psrs[0] == "sapfind") or \
				(len(cmdline.psrs) != 0 and cmdline.psrs[0] == "sapfind3") or \
				(len(cmdline.psrs) != 0 and cmdline.psrs[0] == "tabfind+"):
				self.psrs = find_pulsars(self.rarad, self.decrad, cmdline, cmdline.opts.fwhm_IS/2.)
				if len(cmdline.psrs) != 0 and (cmdline.psrs[0] == "sapfind" or cmdline.psrs[0] == "tabfind+") and len(self.psrs) > 0: 
					self.psrs = self.psrs[:1]

		# initializing SAP beams objects and making the list of SAP beams
		for tid in xrange(self.nrTiedArrayBeams):
			try:
				tab=TABeam(tid, self.sapid, self.nrSubbands)
				tab.update_from_H5(si, cmdline, root, log, wrapper)
				self.tabs.append(tab)
			except Exception:
				msg="Error in getting info for the beam %d:%d. Probably data are missing." % (self.sapid, tid)
				if log != None: log.warning(msg)
				else: print msg
		# select FE tabs
		fetabs=[ii for ii in xrange(len(self.tabs)) if self.tabs[ii].is_coherent and len(self.tabs[ii].stationList) == 1]
		if len(fetabs) > 0:
			fesrc=list(set([self.tabs[ii].specificationType for ii in fetabs]))
			if len(fesrc) == 1: # it means all stations pointed to the same src
				root.FE = True
				for ii in fetabs: self.tabs[ii].specificationType = "flyseye"
			else: # not FE
				for ii in fetabs: self.tabs[ii].specificationType = "manual"

# Class TABeam describes the parameters specific for particular TA beam, which can be truly
# TA beam or beam from individual station in FE mode
class TABeam:
	def __init__(self, id, sapid, nrSubbands):
		self.tabid = id            # ID number of the TA beam
		self.parent_sapid = sapid  # ID of the parent SAP beam

		self.rarad=self.decrad=0
		self.raoffset=self.decoffset=0   # offsets in radians (used for pseudo parset generator)
		self.is_coherent=False     # coherent or incoherent beam
		self.DM=0                  # dispersion measure
		self.stationList=[]        # stations that form this beam, only used to FE now (?) to indicate one station
		self.specificationType=""  # "flyseye" for FE, "manual" (or "ring"?) for coherent 
		self.location=[]           # list of locus nodes with the data for this beam
		self.rawfiles={}           # dictionary that keeps info about all raw files
                                           # key - locus node, value - list of rawfiles on this node (with full path)
                                           # self.location - is just a list of keys
		self.assigned_files=[]     # list of assigned files from the parset file (for this particular beam)
					   # we use this to check its number against the actual number of present files
		self.is_data_available = True # flag that tells whether rawdata available or not for this beam
		self.nrSubbands = nrSubbands  # duplicating number of subbands from parent SAP
		self.numfiles=0            # number of all files for this beam (sum of rawfiles for each node)

	# updating attributes using the parset file
	def update(self, plines, si, cmdline, root, sapRA, sapDEC, log=None, wrapper=None):

	        # getting info about the pointing (offsets from the SAP center)
		res=[ii for ii in plines if "Observation.Beam[%d].TiedArrayBeam[%d].angle1" % (self.parent_sapid, self.tabid) in ii and "AnaBeam" not in ii]
		try:
                	self.raoffset=np.float64(res[0].split("=", 1)[-1].strip())
                	self.rarad=self.raoffset + sapRA
			if self.rarad < 0:           self.rarad += 2.*math.pi
			if self.rarad >= 2.*math.pi: self.rarad -= 2.*math.pi
		except: pass

		res=[ii for ii in plines if "Observation.Beam[%d].TiedArrayBeam[%d].angle2" % (self.parent_sapid, self.tabid) in ii and "AnaBeam" not in ii]
		try:
                	self.decoffset=np.float64(res[0].split("=", 1)[-1].strip())
                	self.decrad=self.decoffset + sapDEC
			if self.decrad > 0.5*math.pi:
				self.decrad = math.pi - self.decrad
				self.rarad += math.pi
				if self.rarad < 0:           self.rarad += 2.*math.pi
				if self.rarad >= 2.*math.pi: self.rarad -= 2.*math.pi
			if self.decrad < -0.5*math.pi:
				self.decrad = -math.pi - self.decrad
				self.rarad += math.pi
				if self.rarad < 0:           self.rarad += 2.*math.pi
				if self.rarad >= 2.*math.pi: self.rarad -= 2.*math.pi
		except: pass

		# getting if beam is coherent or not
		res=[ii.split("=", 1)[-1].strip() for ii in plines if "Observation.Beam[%d].TiedArrayBeam[%d].coherent" % (self.parent_sapid, self.tabid) in ii]
		try:
			tf=re.compile(r"True|true|T|t|1")
			res0=[ii for ii in res if tf.search(ii) is not None]
			if len(res0) > 0: self.is_coherent=True
		except: pass

	        # getting DM
		res=[ii for ii in plines if "Observation.Beam[%d].TiedArrayBeam[%d].dispersionMeasure" % (self.parent_sapid, self.tabid) in ii]
		try:
                	self.DM=np.float64(res[0].split("=", 1)[-1].strip())
		except: pass

	        # getting specification type of the beam
		res=[ii for ii in plines if "Observation.Beam[%d].TiedArrayBeam[%d].specificationType" % (self.parent_sapid, self.tabid) in ii]
		try:
	       		self.specificationType=res[0].split("=", 1)[-1].strip()
		except: pass

		# getting station list
		res=[ii for ii in plines if "Observation.Beam[%d].TiedArrayBeam[%d].stationList" % (self.parent_sapid, self.tabid) in ii]
		try:
                        self.stationList = res[0].split("=", 1)[-1].strip().lstrip("[").rstrip("]").split(",")
		except: pass

		# checking where the raw data are
		missing_nodes=[]
		if not cmdline.opts.is_auto:
			if not cmdline.opts.is_globalfs:
				missing_nodes=self.get_tab_rawdata(si, cmdline, root, log)
			else:
				self.rawfiles_manual_mapping(si, cmdline, root, log)
		else: # forming 'rawfiles' dictionary from the given input data
			# if we are not processing all splits then we skip files that we do not need
			#
			# if we are on CEP4 in auto mode, then in the pipeline parset, all "host" values are the same = "CEP4"
			# so we still need to use manual mapping
			if si.cluster_headnode[:5] == "head0" or si.cluster_headnode[:3] == "cpu": # CEP4
				self.rawfiles_manual_mapping(si, cmdline, root, log)	
			else:
	                        hosts = [x.host for x in wrapper.input_data['coherent'].data]
        	                hosts.extend([x.host for x in wrapper.input_data['incoherent'].data])
                	        files = [x.file for x in wrapper.input_data['coherent'].data]
                        	files.extend([x.file for x in wrapper.input_data['incoherent'].data])
				for (loc, infile) in zip(hosts, [ff.replace(".h5", ".raw") for ff in files]):
					if self.is_coherent:
						if cmdline.opts.nsplitsCS != root.nsplitsCS:
							fpart=int(infile.split("/")[-1].split("_P")[-1].split("_")[0])
							if fpart < cmdline.opts.first_freq_splitCS or fpart >= cmdline.opts.first_freq_splitCS + cmdline.opts.nsplitsCS:
								continue
					else:
						if cmdline.opts.nsplitsIS != root.nsplitsIS:
							fpart=int(infile.split("/")[-1].split("_P")[-1].split("_")[0])
							if fpart < cmdline.opts.first_freq_splitIS or fpart >= cmdline.opts.first_freq_splitIS + cmdline.opts.nsplitsIS:
								continue
					if loc in self.rawfiles: self.rawfiles[loc].append(infile)
					else: self.rawfiles[loc]=[infile]
				
			self.location=self.rawfiles.keys()
			# getting the total number of files available
			self.numfiles = sum([len(self.rawfiles[loc]) for loc in self.location])

		# Now getting the list of assigned files for this beam from the Parset file
		res=[ii.split("=", 1)[-1].strip().lstrip("[").rstrip("]") for ii in plines if "%s_SAP%03d_B%03d" % (root.id, self.parent_sapid, self.tabid) in ii]
		try:
			self.assigned_files=[el for el in res[0].split(",") if "%s_SAP%03d_B%03d" % (root.id, self.parent_sapid, self.tabid) in el]
			# if we are processing not all the frequency splits, then we include only files from those splits that we need	
			if self.is_coherent:
				if cmdline.opts.nsplitsCS != root.nsplitsCS:
					tmp_assigned_files=[]
					for sp in xrange(cmdline.opts.first_freq_splitCS, cmdline.opts.first_freq_splitCS + cmdline.opts.nsplitsCS):
						tmp_assigned_files.extend([el for el in self.assigned_files if "_P%03d" % (sp) in el])
					self.assigned_files=tmp_assigned_files
			else:
				if cmdline.opts.nsplitsIS != root.nsplitsIS:
					tmp_assigned_files=[]
					for sp in xrange(cmdline.opts.first_freq_splitIS, cmdline.opts.first_freq_splitIS + cmdline.opts.nsplitsIS):
						tmp_assigned_files.extend([el for el in self.assigned_files if "_P%03d" % (sp) in el])
					self.assigned_files=tmp_assigned_files
		except: pass

		# Now checking that the number of available files is the same as assigned number
		log.info("%s" % (",".join(self.rawfiles)))
		if self.numfiles < len(self.assigned_files):
			self.is_data_available = False
			missing_files=self.assigned_files[:]
			# checking first what files are missing
			for loc in self.location:
				missing_files=list(set(missing_files)-set([s.split("/")[-1] for s in self.rawfiles[loc]]))
			if len(missing_files) == 0:
				msg="This is weird... There are should be at least one file missing..."	
				if log != None: log.warning(msg)
				else: print msg
			msg="Warning: The number of available files (%d) is less than assigned (%d) for the beam %d:%d!" % (self.numfiles, np.size(self.assigned_files), self.parent_sapid, self.tabid)
			if len(missing_files) > 0:
				msg += "\n[Missing files]: %s" % (",".join(missing_files))
			if len(missing_nodes) > 0:
				msg += "\n[Missing nodes]: %s" % (",".join(missing_nodes))
			if log != None: log.warning(msg)
			else: print msg

		if not cmdline.opts.is_auto:
			# if tab.location is empty we still have to fill it based on where the processed data are, because Pipeine _needs_ to know
			# on what nodes to start processing even if it is only for re-doing plots...
			# So, we try to look for log-file for the particular beam to determine the locus node
			if len(self.location) == 0:
				if not cmdline.opts.is_globalfs:
					missing_nodes=self.get_processed_data(si, cmdline, root, missing_nodes, log)


	# checking where the raw data are
	def get_processed_data(self, si, cmdline, root, missnodes, log=None, nsplits=0):

		if nsplits == 0:
			if self.is_coherent: nsplits = root.nsplitsCS
			else: nsplits = root.nsplitsIS
		if len(si.alive_nodes) > 0:
       	        	cexeclocus=si.cexec_nodes[si.alive_nodes[0]]
               		if len(si.alive_nodes) > 1:
                       		for s in si.alive_nodes[1:]:
					if s in si.hoover_nodes: pass
					else: cexeclocus += ",%s" % (si.cexec_nodes[s].split(":")[1])
                       		for s in si.alive_nodes[1:]:
					if s in si.hoover_nodes: cexeclocus += " %s" % (si.cexec_nodes[s])
					else: pass
		if len(self.stationList) != 0 and self.stationList[0] != "":
			cmd="%s %s 'find /*/%s*/%s_red*/*/SAP%d \( -name \"%s\" -or -name \"BEAM%d\" \) -type d' | grep -v denied | grep -v such | grep -v match | grep -v xauth | grep -v connect | grep -v closed | egrep -v \'\\*\\*\\*\\*\\*\'" % (si.cexeccmd, cexeclocus, si.processed_dir_prefix, root.id, self.parent_sapid, self.stationList[0], self.tabid)
		else:
			cmd="%s %s 'find /*/%s*/%s_red*/*/SAP%d -name \"BEAM%d\" -type d' | grep -v denied | grep -v such | grep -v match | grep -v xauth | grep -v connect | grep -v closed | egrep -v \'\\*\\*\\*\\*\\*\'" % (si.cexeccmd, cexeclocus, si.processed_dir_prefix, root.id, self.parent_sapid, self.tabid)
       	        cexec_output=[line.strip() for line in os.popen(cmd).readlines()]
		#log.info("CMD = %s" % (cmd))
		#log.info("OUTPUT = %s" % (", ".join(cexec_output)))
		loc=""
		for l in xrange(len(cexec_output)):
			if re.match("^-----", cexec_output[l]) is not None:
				loc=cexec_output[l].split(" ")[1].split("-")[0].split(".")[0]
			else: # it means that we found the file and loc now has the locus node name where processed data are
			      # but we continue loop, as there can be several locus nodes with processed data for 1 beam if there were several splits
				self.location.append(loc)

			if len(self.location) > 1 and not self.is_coherent: # it means that our beam is Incoherent and in the location list there are likely summary directory is also included
				if si.summary_nodes["IS"] in self.location: # if summary node is in the list then remove the first occurrence from the list
					self.location.remove(si.summary_nodes["IS"])
		if loc=="":
			msg="Warning: Neither raw or even processed data available for beam %d:%d" % (self.parent_sapid, self.tabid)
			if log != None: log.warning(msg)
			else: print msg
		else: self.location = np.unique(self.location)

		#for node in self.rawfiles.keys():
		#	log.info("NODE = %s" % (node))
		#	log.info("FILES = %s" % (", ".join(self.rawfiles[node])))

		# Now we are getting the list of potential rawfiles (even if they were deleted)
		if len(self.location) != 0:
               		cexeclocus=si.cexec_nodes[self.location[0]]
               		if len(self.location) > 1:
	      			for s in self.location[1:]:
					cexeclocus += ",%s" % (si.cexec_nodes[s].split(":")[1])
			if len(self.stationList) != 0 and self.stationList[0] != "":
				#cmd="%s %s 'ls -1L $(find /*/%s*/%s_red*/*/SAP%d -name \"%s_SAP%03d_B%03d_S*_bf.h5\")' | grep -v denied | grep -v such | grep -v match | grep -v xauth | grep -v connect | grep -v closed | egrep -v \'\\*\\*\\*\\*\\*\'" % (si.cexeccmd, cexeclocus, si.processed_dir_prefix, root.id, self.parent_sapid, root.id, self.parent_sapid, self.tabid)
				cmd="%s %s 'find -L /*/%s*/%s_red*/*/SAP%d -name \"%s_SAP%03d_B%03d_S*_bf.h5\"' | grep -v denied | grep -v such | grep -v match | grep -v xauth | grep -v connect | grep -v closed | egrep -v \'\\*\\*\\*\\*\\*\'" % (si.cexeccmd, cexeclocus, si.processed_dir_prefix, root.id, self.parent_sapid, root.id, self.parent_sapid, self.tabid)
			else:
				#cmd="%s %s 'ls -1L $(find /*/%s*/%s_red*/*/SAP%d/BEAM%d -name \"%s_SAP%03d_B%03d_S*_bf.h5\")' | grep -v denied | grep -v such | grep -v match | grep -v xauth | grep -v connect | grep -v closed | egrep -v \'\\*\\*\\*\\*\\*\'" % (si.cexeccmd, cexeclocus, si.processed_dir_prefix, root.id, self.parent_sapid, self.tabid, root.id, self.parent_sapid, self.tabid)
				cmd="%s %s 'find -L /*/%s*/%s_red*/*/SAP%d/BEAM%d -name \"%s_SAP%03d_B%03d_S*_bf.h5\"' | grep -v denied | grep -v such | grep -v match | grep -v xauth | grep -v connect | grep -v closed | egrep -v \'\\*\\*\\*\\*\\*\'" % (si.cexeccmd, cexeclocus, si.processed_dir_prefix, root.id, self.parent_sapid, self.tabid, root.id, self.parent_sapid, self.tabid)
               		cexec_output=[line.strip() for line in os.popen(cmd).readlines()]
			#log.info("CMD = %s" % (cmd))
			#log.info("OUTPUT = %s" % (", ".join(cexec_output)))
			for l in xrange(len(cexec_output)):
				if re.match("^-----", cexec_output[l]) is not None:
					loc=cexec_output[l].split(" ")[1].split("-")[0].split(".")[0]
				else:
					# first we are checking that this file (corresponding raw-file in the /data/ObsId) DO NOT already present
					# in the current dictionary self.rawfiles. This is to deal with the case when rawdata are in /home area
					# that is visible from all locus nodes. And also to exclude situation when both raw data and processed data
					# are present
					rawfile_exist=False
					#rfile="%s%s/%s%s/%s%s" % (si.rawdir, cexec_output[l].split(si.rawdir)[-1].split("/")[0], si.rawdir_prefix_specificator, root.id, si.rawdir_suffix_specificator, cexec_output[l].replace(".h5", ".raw").split("/")[-1])
                                        rfile_cmd="%s %s 'ls -d %s%s/%s%s/%s' | grep -v denied | grep -v such | grep -v match | grep -v xauth | grep -v connect | grep -v closed | egrep -v \'\\*\\*\\*\\*\\*\'" % (si.cex     eccmd, cexeclocus, si.rawdir, cexec_output[l].split(si.rawdir)[-1].split("/")[0], si.rawdir_prefix_specificator, root.id, si.rawdir_suffix_specificator)
                                        rfile_dir=[line.strip() for line in os.popen(rfile_cmd).readlines()][1]
                                        rfile="%s%s" % (rfile_dir, cexec_output[l].replace(".h5", ".raw").split("/")[-1])
					for cloc in self.rawfiles:
						if rfile in self.rawfiles[cloc]:
							rawfile_exist=True
							break
					# if we are processing not all the frequency splits, then we include only files from those splits that we need
					if self.is_coherent:
						if cmdline.opts.nsplitsCS != nsplits:
							fpart=int(cexec_output[l].split("/")[-1].split("_P")[-1].split("_")[0])
							if fpart < cmdline.opts.first_freq_splitCS or fpart >= cmdline.opts.first_freq_splitCS + cmdline.opts.nsplitsCS:
								rawfile_exist=True
					else:
						if cmdline.opts.nsplitsIS != nsplits:
							fpart=int(cexec_output[l].split("/")[-1].split("_P")[-1].split("_")[0])
							if fpart < cmdline.opts.first_freq_splitIS or fpart >= cmdline.opts.first_freq_splitIS + cmdline.opts.nsplitsIS:
								rawfile_exist=True
					# adding the file if it's not already added
					if not rawfile_exist:
						# converting the full filename in the processing location to file name in the /data/ObsId
						#rfile="%s%s/%s%s/%s%s" % (si.rawdir, cexec_output[l].split(si.rawdir)[-1].split("/")[0], si.rawdir_prefix_specificator, root.id, si.rawdir_suffix_specificator, cexec_output[l].replace(".h5", ".raw").split("/")[-1])
						if loc in self.rawfiles: self.rawfiles[loc].append(rfile)
						else: self.rawfiles[loc]=[rfile]	
					# excluding this locus node from the list of missing nodes
					#if loc in missnodes:
					#	missnodes=list(set(missnodes)-set([loc]))
		#log.info("+++++++ END ++++++++")
		#for node in self.rawfiles.keys():
		#	log.info("NODE = %s" % (node))
		#	log.info("FILES = %s" % (", ".join(self.rawfiles[node])))
		if len(self.rawfiles.keys()) != 0:
			missnodes=self.rawfiles.keys()[:]
		return missnodes

	# filling in self.rawfiles and self.location by manually assigning nodes for all raw files 
	# for the case when we use GlobalFS. If Slurm is used, then in anyway any available node will be used based on the Slurm scheduling
	# but if we use SSH, then we should have this mapping
	def rawfiles_manual_mapping(self, si, cmdline, root, log=None, nsplits=0):
		if nsplits == 0:
			if self.is_coherent: nsplits = root.nsplitsCS
			else: nsplits = root.nsplitsIS
		cmd="find %s*/%s%s -name \"%s_SAP%03d_B%03d_S*_bf.raw\"" % (si.rawdir, si.rawdir_prefix_specificator, root.id, root.id, self.parent_sapid, self.tabid)
       	        fileslist=[line.strip() for line in os.popen(cmd).readlines()]
		avnodes=len(root.nodeslist)
		for ff in xrange(len(fileslist)):
			# first we are checking that this file (exactly the same with full path) DO NOT already present
			# in the current dictionary self.rawfiles. This is to deal with the case when rawdata are in /home area
			# that is visible from all locus nodes
			rawfile_exist=False
			for cloc in self.rawfiles:
				if fileslist[ff] in self.rawfiles[cloc]:
					rawfile_exist=True
					break
			# if we are processing not all the frequency splits, then we include only files from those splits that we need
			if self.is_coherent:
				if cmdline.opts.nsplitsCS != nsplits:
					fpart=int(fileslist[ff].split("/")[-1].split("_P")[-1].split("_")[0])
					if fpart < cmdline.opts.first_freq_splitCS or fpart >= cmdline.opts.first_freq_splitCS + cmdline.opts.nsplitsCS:
						rawfile_exist=True
			else:
				if cmdline.opts.nsplitsIS != nsplits:
					fpart=int(fileslist[ff].split("/")[-1].split("_P")[-1].split("_")[0])
					if fpart < cmdline.opts.first_freq_splitIS or fpart >= cmdline.opts.first_freq_splitIS + cmdline.opts.nsplitsIS:
						rawfile_exist=True
			# adding the file if it's not already added
			if not rawfile_exist:
				loc=root.nodeslist[ff%avnodes]
				if loc in self.rawfiles: self.rawfiles[loc].append(fileslist[ff])
				else: self.rawfiles[loc]=[fileslist[ff]]	

		# list of all nodes
		self.location=self.rawfiles.keys()
		if len(self.location) == 0:
			self.is_data_available = False
			msg="No data available for beam %d:%d" % (self.parent_sapid, self.tabid)
			if log != None: log.warning(msg)
			else: print msg
		else:
			# getting the total number of files available
			self.numfiles = sum([len(self.rawfiles[loc]) for loc in self.location])


	# checking where the raw data are
	def get_tab_rawdata(self, si, cmdline, root, log=None, nsplits=0):
		# Determining where the raw data are....
                # forming string with all locus nodes needed to check in one cexec command
		# here we are using only nodes that are alive
		#log.info("++++++++ get_tab_rawdata ++++++++")
		missing_nodes=root.assigned_nodeslist[:]
		if cmdline.opts.is_locate_rawdata: nodeslist=si.alive_nodes
		else: nodeslist=root.nodeslist
		if nsplits == 0:
			if self.is_coherent: nsplits = root.nsplitsCS
			else: nsplits = root.nsplitsIS
		if len(nodeslist) > 0:
       	        	cexeclocus=si.cexec_nodes[nodeslist[0]]
               		if len(nodeslist) > 1:
                       		for s in nodeslist[1:]:
					cexeclocus += ",%s" % (si.cexec_nodes[s].split(":")[1])
			cmd="%s %s 'find %s*/%s%s -name \"%s_SAP%03d_B%03d_S*_bf.raw\"' | grep -v denied | grep -v such | grep -v match | grep -v xauth | grep -v connect | grep -v closed | egrep -v \'\\*\\*\\*\\*\\*\'" % (si.cexeccmd, cexeclocus, si.rawdir, si.rawdir_prefix_specificator, root.id, root.id, self.parent_sapid, self.tabid)
       	        	cexec_output=[line.strip() for line in os.popen(cmd).readlines()]
			#log.info("CMD = %s" % (cmd))
			#log.info("OUTPUT = %s" % (", ".join(cexec_output)))
			for l in xrange(len(cexec_output)):
				#log.info("l = %d (%d), line = %s" % (l, len(cexec_output), cexec_output[l]))
				if re.match("^-----", cexec_output[l]) is not None:
					loc=cexec_output[l].split(" ")[1].split("-")[0].split(".")[0]
				else:
					# first we are checking that this file (exactly the same with full path) DO NOT already present
					# in the current dictionary self.rawfiles. This is to deal with the case when rawdata are in /home area
					# that is visible from all locus nodes
					rawfile_exist=False
					for cloc in self.rawfiles:
						if cexec_output[l] in self.rawfiles[cloc]:
							rawfile_exist=True
							break
					# if we are processing not all the frequency splits, then we include only files from those splits that we need
					if self.is_coherent:
						if cmdline.opts.nsplitsCS != nsplits:
							fpart=int(cexec_output[l].split("/")[-1].split("_P")[-1].split("_")[0])
							if fpart < cmdline.opts.first_freq_splitCS or fpart >= cmdline.opts.first_freq_splitCS + cmdline.opts.nsplitsCS:
								rawfile_exist=True
					else:
						if cmdline.opts.nsplitsIS != nsplits:
							fpart=int(cexec_output[l].split("/")[-1].split("_P")[-1].split("_")[0])
							if fpart < cmdline.opts.first_freq_splitIS or fpart >= cmdline.opts.first_freq_splitIS + cmdline.opts.nsplitsIS:
								rawfile_exist=True
					# adding the file if it's not already added
					if not rawfile_exist:
						if loc in self.rawfiles: self.rawfiles[loc].append(cexec_output[l])
						else: self.rawfiles[loc]=[cexec_output[l]]	
					# excluding this locus node from the list of missing nodes
					if loc in missing_nodes:
						missing_nodes=list(set(missing_nodes)-set([loc]))

			# list of all nodes
			self.location=self.rawfiles.keys()
			#log.info("LOCATIONS: %s" % (", ".join(self.location)))
			#log.info("MISSING NODES = %s" % (", ".join(missing_nodes)))
			#log.info("ALL FILES=%s" % (", ".join([", ".join(self.rawfiles[loc]) for loc in self.location])))
			#for loc in self.location:
			#	log.info("LOCATION = %s   FILES = %d" % (loc, len(self.rawfiles[loc])))
			#	for ff in self.rawfiles[loc]:
			#		log.info("%s" % (ff))
			if len(self.location) == 0:
				self.is_data_available = False
				msg="No data available for beam %d:%d" % (self.parent_sapid, self.tabid)
				if log != None: log.warning(msg)
				else: print msg
			else:
				# getting the total number of files available
				self.numfiles = sum([len(self.rawfiles[loc]) for loc in self.location])
				#log.info("NUMFILES = %d" % (self.numfiles))
		return missing_nodes


	# updating attributes using .h5 files
	def update_from_H5(self, si, cmdline, root, log=None, wrapper=None):

		# get list of .h5 files for the current TAB
		tabh5list = glob.glob("%s/%s_SAP%03d_B%03d_S*_bf.h5" % (si.get_logdir(), root.id, self.parent_sapid, self.tabid))
		# first, we read any 1 file to collect general info from the BEAM group
		t5 = h5py.File(tabh5list[0], 'r')
		self.rarad = np.float64((t5['SUB_ARRAY_POINTING_%03d' % (self.parent_sapid)]['BEAM_%03d' % (self.tabid)].attrs['POINT_RA'] * math.pi) / 180.)
		self.decrad = np.float64((t5['SUB_ARRAY_POINTING_%03d' % (self.parent_sapid)]['BEAM_%03d' % (self.tabid)].attrs['POINT_DEC'] * math.pi) / 180.)
		self.raoffset = np.float64((t5['SUB_ARRAY_POINTING_%03d' % (self.parent_sapid)]['BEAM_%03d' % (self.tabid)].attrs['POINT_OFFSET_RA'] * math.pi) / 180.)
		self.decoffset = np.float64((t5['SUB_ARRAY_POINTING_%03d' % (self.parent_sapid)]['BEAM_%03d' % (self.tabid)].attrs['POINT_OFFSET_DEC'] * math.pi) / 180.)
		self.DM = np.float64(t5['SUB_ARRAY_POINTING_%03d' % (self.parent_sapid)]['BEAM_%03d' % (self.tabid)].attrs['DISPERSION_MEASURE'])
		signal_sum = str(t5['SUB_ARRAY_POINTING_%03d' % (self.parent_sapid)]['BEAM_%03d' % (self.tabid)].attrs['SIGNAL_SUM']).lower()
		nof_stokes = int(t5['SUB_ARRAY_POINTING_%03d' % (self.parent_sapid)]['BEAM_%03d' % (self.tabid)].attrs['OBSERVATION_NOF_STOKES'])
		if signal_sum == "coherent": self.is_coherent = True
		self.stationList = [s for s in t5['SUB_ARRAY_POINTING_%03d' % (self.parent_sapid)]['BEAM_%03d' % (self.tabid)].attrs['STATIONS_LIST']]
		self.specificationType = "manual"
		if self.is_coherent and len(self.stationList) == 1:
			# for now we just put the target name in this field. Later in the Observation class we will check if targets
			# are the same for all these tabs, and then it means it is FE (flyseye)
			self.specificationType = t5['SUB_ARRAY_POINTING_%03d' % (self.parent_sapid)]['BEAM_%03d' % (self.tabid)].attrs['TARGETS'][0]
		t5.close()

		# open part0 file to get nsubs/file and to determine number of splits as it is not yet determined at time when this function is called
		t5p0 = h5py.File("%s/%s_SAP%03d_B%03d_S0_P000_bf.h5" % (si.get_logdir(), root.id, self.parent_sapid, self.tabid), 'r')
		nsubs_file = t5p0['SUB_ARRAY_POINTING_%03d' % (self.parent_sapid)]['BEAM_%03d' % (self.tabid)]['STOKES_0'].attrs['NOF_SUBBANDS']
		t5p0.close()
		nsplits = int(self.nrSubbands / nsubs_file)
		if self.nrSubbands % nsubs_file != 0: nsplits += 1
		# changing cmdline split-related options in the _copy_ of Cmdline class
		if self.is_coherent:
			# for CS/CV
			if cmdline.opts.first_freq_splitCS >= nsplits: cmdline.opts.first_freq_splitCS = 0
			if cmdline.opts.nsplitsCS == -1: cmdline.opts.nsplitsCS = nsplits
			if cmdline.opts.first_freq_splitCS + cmdline.opts.nsplitsCS > nsplits:
				cmdline.opts.nsplitsCS -= (cmdline.opts.first_freq_splitCS + cmdline.opts.nsplitsCS - nsplits)
		else:
			# for IS
			if cmdline.opts.first_freq_splitIS >= nsplits: cmdline.opts.first_freq_splitIS = 0
			if cmdline.opts.nsplitsIS == -1: cmdline.opts.nsplitsIS = nsplits
			if cmdline.opts.first_freq_splitIS + cmdline.opts.nsplitsIS > nsplits:
				cmdline.opts.nsplitsIS -= (cmdline.opts.first_freq_splitIS + cmdline.opts.nsplitsIS - nsplits)

		# checking where the raw data are
		missing_nodes=[]
		if not cmdline.opts.is_auto:
			if not cmdline.opts.is_globalfs:
				missing_nodes=self.get_tab_rawdata(si, cmdline, root, log, nsplits)
			else:
				self.rawfiles_manual_mapping(si, cmdline, root, log, nsplits)
		else: # forming 'rawfiles' dictionary from the given input data
			# if we are not processing all splits then we skip files that we do not need
			#
			# if we are on CEP4 in auto mode, then in the pipeline parset, all "host" values are the same = "CEP4"
			# so we still need to use manual mapping
			if si.cluster_headnode[:5] == "head0" or si.cluster_headnode[:3] == "cpu": # CEP4
				self.rawfiles_manual_mapping(si, cmdline, root, log)	
			else:
	                        hosts = [x.host for x in wrapper.input_data['coherent'].data if "_SAP%03d_B%03d_" % (self.parent_sapid, self.tabid) in x.file]
        	                hosts.extend([x.host for x in wrapper.input_data['incoherent'].data if "_SAP%03d_B%03d_" % (self.parent_sapid, self.tabid) in x.file])
                	        files = [x.file for x in wrapper.input_data['coherent'].data if "_SAP%03d_B%03d_" % (self.parent_sapid, self.tabid) in x.file]
                        	files.extend([x.file for x in wrapper.input_data['incoherent'].data if "_SAP%03d_B%03d_" % (self.parent_sapid, self.tabid) in x.file])
				for (loc, infile) in zip(hosts, [ff.replace(".h5", ".raw") for ff in files]):
					if self.is_coherent:
						if cmdline.opts.nsplitsCS != nsplits:
							fpart=int(infile.split("/")[-1].split("_P")[-1].split("_")[0])
							if fpart < cmdline.opts.first_freq_splitCS or fpart >= cmdline.opts.first_freq_splitCS + cmdline.opts.nsplitsCS:
								continue
					else:
						if cmdline.opts.nsplitsIS != nsplits:
							fpart=int(infile.split("/")[-1].split("_P")[-1].split("_")[0])
							if fpart < cmdline.opts.first_freq_splitIS or fpart >= cmdline.opts.first_freq_splitIS + cmdline.opts.nsplitsIS:
								continue
					if loc in self.rawfiles: self.rawfiles[loc].append(infile)
					else: self.rawfiles[loc]=[infile]
				
			self.location=self.rawfiles.keys()
			# getting the total number of files available
			self.numfiles = sum([len(self.rawfiles[loc]) for loc in self.location])

		# assigned files
		self.assigned_files=["%s.raw" % (f.split("/")[-1].split(".h5")[0]) for f in tabh5list]
		# if we are processing not all the frequency splits, then we include only files from those splits that we need	
		if self.is_coherent:
			if cmdline.opts.nsplitsCS != nsplits:
				tmp_assigned_files=[]
				for sp in xrange(cmdline.opts.first_freq_splitCS, cmdline.opts.first_freq_splitCS + cmdline.opts.nsplitsCS):
					tmp_assigned_files.extend([el for el in self.assigned_files if "_P%03d" % (sp) in el])
				self.assigned_files=tmp_assigned_files
		else:
			if cmdline.opts.nsplitsIS != nsplits:
				tmp_assigned_files=[]
				for sp in xrange(cmdline.opts.first_freq_splitIS, cmdline.opts.first_freq_splitIS + cmdline.opts.nsplitsIS):
					tmp_assigned_files.extend([el for el in self.assigned_files if "_P%03d" % (sp) in el])
				self.assigned_files=tmp_assigned_files
		
		# Now checking that the number of available files is the same as assigned number
		if self.numfiles < len(self.assigned_files):
			self.is_data_available = False
			missing_files=self.assigned_files[:]
			# checking first what files are missing
			for loc in self.location:
				missing_files=list(set(missing_files)-set([s.split("/")[-1] for s in self.rawfiles[loc]]))
			if len(missing_files) == 0:
				msg="This is weird... There are should be at least one file missing..."	
				if log != None: log.warning(msg)
				else: print msg
			if cmdline.opts.is_auto:
				msg="Warning: The number of available files (%d) is less than assigned (%d) for the beam %d:%d!" % (self.numfiles, np.size(self.assigned_files), self.parent_sapid, self.tabid)
				if len(missing_files) > 0:
					msg += "\n[Missing files]: %s" % (",".join(missing_files))
				if len(missing_nodes) > 0:
					msg += "\n[Missing nodes]: %s" % (",".join(missing_nodes))
				if log != None: log.warning(msg)
				else: print msg

		if not cmdline.opts.is_auto:
			# if tab.location is empty we still have to fill it based on where the processed data are, because Pipeline _needs_ to know
			# on what nodes to start processing even if it is only for re-doing plots...
			# So, we try to look for log-file for the particular beam to determine the locus node
			if len(self.location) == 0 or self.numfiles < len(self.assigned_files):
				if not cmdline.opts.is_globalfs:
					missing_nodes=self.get_processed_data(si, cmdline, root, missing_nodes, log, nsplits)
			if self.numfiles < len(self.assigned_files):
				msg="Warning: The number of available files (%d) is less than assigned (%d) for the beam %d:%d!" % (self.numfiles, np.size(self.assigned_files), self.parent_sapid, self.tabid)
				if len(missing_files) > 0:
					msg += "\n[Missing files]: %s" % (",".join(missing_files))
				if len(missing_nodes) > 0:
					msg += "\n[Missing nodes]: %s" % (",".join(missing_nodes))
				if log != None: log.warning(msg)
				else: print msg

		# if one (or more) of the locus nodes is down, then our number of assigned files will be the same, so we won't notice if anything is missing
		# so, we have to calculate how many files we expect from knowing number of subs per file, or nsplits
		if nof_stokes == 1: nexpected_files = nsplits
		else: nexpected_files = 4 * nsplits

		if self.numfiles < nexpected_files:
			self.is_data_available = False
			msg="Warning: The number of available files (%d) is less than expected (%d) for the beam %d:%d!" % (self.numfiles, nexpected_files, self.parent_sapid, self.tabid)
			if log != None: log.warning(msg)
			else: print msg


# Class Observation with info from the parset file
class Observation:
	# si - system info
	# log - logger class
	def __init__(self, id, si, cmdline, log=None, wrapper=None):
		self.id = id
		self.parset = cmdline.opts.parset
		self.parset_dir = si.parset_dir

		self.project = "" # the name of the campaign
		self.projectPI = "" # PI of the project

		self.nrTiedArrayBeams = 0 # total number of TABs in all SAPs
		self.nrBeams = 0  # number of station beams
		self.saps=[]      # list of SAP beams (objects of SAPBeam class)

		self.startdate=self.starttime=""
                self.duration=0  # duration of obs in seconds

                self.antenna=self.antenna_config=self.band=self.bandFilter=""
		self.nstations=self.ncorestations=self.nremotestations=0
		self.stations=[]
		self.assigned_nodeslist=[]
		self.nodeslist=[]
		self.IM=self.IS=self.CS=self.CV=self.FE=self.OCD=False    # False until it's checked to be otherwise
		self.stokesIS = ""         # Stokes parameters for IS and CS
		self.stokesCS = ""

		self.nrSubbands = 0        # number of subbands
		self.subbandList=""        # range of subbands, e.g. 77..320
		self.subbands=[]           # list of subbands
		self.subbandWidth = 0      # width of subband in kHz
		self.nrSubsPerFileIS = 0   # number of subbands per file
		self.nrSubsPerFileCS = 0    
		self.nrChanPerSubIS = 0    # number of channels per subband (for IS)
		self.nrChanPerSubCS = 0    # number of channels per subband (for CS)
		self.nsplitsCS = 0         # number of frequency splits for CS/CV
		self.nsplitsIS = 0         # number of frequency splits for IS
		self.sampleClock = 0       # clock in MHz (200 or 160)
		self.downsample_factorIS = 0  # time integration (downsampling) factor, can be different for IS and CS 
		self.downsample_factorCS = 0  
		self.samplingIS = 0        # sampling interval in ms (depends on on integration steps, clock, number of channels)
		self.samplingCS = 0
		self.lower_band_edge = 0   # lowest frequency for the band in MHz
		self.bw = 0                # bandwidth (in MHz) - recorded bandwidth
		self.freq_extent = 0       # freq separation between lowest and highest channels recorded (if there are gaps in frequency)
		self.cfreq = 0             # central freq (in MHz)

		if cmdline.opts.is_cobalt:
			# collect metadata from .h5 files rather than from the parset file
			self.update_from_H5(si, cmdline, log, wrapper)
		else:
			# check if parset file exists and if so, then update the info
			if self.is_parset(): self.update(si, cmdline, log, wrapper)
			else:
				msg="Can't find the parset file '%s' for ObsID = %s" % (self.parset, self.id)
				if log != None: log.error(msg)
				else: print msg
				raise Exception

	# return True if parset file was found, and False otherwise
	def is_parset (self):
		if self.parset != "": 
			if os.path.exists(self.parset):	return True
			else: return False
		else:   # checking old parset location
			self.parset = "%s/%s/%s.parset" % (self.parset_dir, self.id, self.id)
			if os.path.exists(self.parset):	return True
			else:   # checking new parset location
				self.parset = "%s/%s.parset" % (self.parset_dir, self.id)
				if os.path.exists(self.parset):	return True
				else: return False

	# parsing the string with ranges of subbands recorded to get list of subbands
	def getSubbands(self, sblist):
		subs=[]
		sbparts=sblist.split(",")
		for ss in sbparts:
			sedges=ss.split("..")
			if len(sedges) == 1: subs.append(int(sedges[0]))
			else: subs.extend(range(int(sedges[0]), int(sedges[1])+1))
		subs.sort() # sorting with smallest being the first
		return subs

	# update info based on parset file
	def update (self, si, cmdline, log=None, wrapper=None):

                # reading first parset file into the list of lines
                f = open(self.parset, 'r')
                # ignoring also comments and empty lines
                comments_and_empty=re.compile(r"(^\s*#+.*$)|(^\s*$)")
                plines = [ff for ff in f.read().splitlines() if comments_and_empty.search(ff) is None]
                f.close()

                # Getting the Date of observation
                res=[ii for ii in plines if "Observation.startTime" in ii]
		try:
			self.startdate = res[0].split("=", 1)[-1].strip().strip("'\"")
			self.starttime = self.startdate.split(" ")[1]
			# checking if startdate is after Jan 27, 2012 or not. If not throw the error
			c1 = time.strptime(self.startdate, "%Y-%m-%d %H:%M:%S")
			c2 = time.strptime("2012-01-27 00:00:00", "%Y-%m-%d %H:%M:%S")
		except: pass
		else:
			if time.mktime(c1) < time.mktime(c2):
				msg="The observation was before Jan 27, 2012. Python version of pulsar pipeline can be run\n\
only for observations that were taken after this date"
				if log!=None: log.error(msg)
				else: print msg
				raise Exception

                # Getting the Duration
                res=[ii for ii in plines if "Observation.stopTime" in ii]
		try:
                        stoptime = res[0].split("=", 1)[-1].strip().strip("'\"")
			c1 = time.strptime(self.startdate, "%Y-%m-%d %H:%M:%S")
			c2 = time.strptime(stoptime, "%Y-%m-%d %H:%M:%S")
			self.duration=time.mktime(c2)-time.mktime(c1)  # difference in seconds
			self.startdate  = self.startdate.split(" ")[0]
		except: pass

                # Getting the Antenna info (HBA or LBA) and Filter setting
                res=[ii for ii in plines if "Observation.bandFilter" in ii]
		try:
			self.bandFilter = res[0].split("=", 1)[-1].strip()
			self.antenna = res[0].split("=", 1)[-1].strip().split("_")[0]
			self.band = res[0].split("=", 1)[-1].strip().split("A_")[-1]
		except: pass

		# changing cmdline FWHM-related options in the _copy_ of Cmdline class
                if cmdline.opts.fwhm_CS < 0.0:
                        if self.antenna == "HBA": cmdline.opts.fwhm_CS = si.fwhm_hba
                        if self.antenna == "LBA": cmdline.opts.fwhm_CS = si.fwhm_lba
                if cmdline.opts.fwhm_IS < 0.0:
                        if self.antenna == "HBA": cmdline.opts.fwhm_IS = si.fov_hba
                        if self.antenna == "LBA": cmdline.opts.fwhm_IS = si.fov_lba

                # Getting the Antenna config (LBA_OUTER, LBA_INNER, HBA_JOINED, etc.)
                res=[ii for ii in plines if "Observation.antennaSet" in ii]
		try:
			self.antenna_config = res[0].split("=", 1)[-1].strip()
		except: pass

                # Getting the name of the Campaign
                res=[ii for ii in plines if "Observation.Campaign.name" in ii]
		try:
			self.project = res[0].split("=", 1)[-1].strip().strip("'\"")
		except: pass

		# Getting the name of the Campaign's PI
		res=[ii for ii in plines if "Observation.Campaign.PI" in ii]
		try:
                	self.projectPI=res[0].split("=", 1)[-1].strip().strip("'\"")
		except: pass

                # Getting the stations and their number (including separately the number of CS and RS)
                res=[ii for ii in plines if "OLAP.storageStationNames" in ii]
		try:
			stations = res[0].split("=", 1)[-1].strip().lstrip("[").rstrip("]")
			# removing LBA and HBA from station names, replacing HBA ears HBA0 to /0 and HBA1 to /1
                        stations=stations.replace("HBA0", "/0")
                        stations=stations.replace("HBA1", "/1")
                        stations=stations.replace("HBA", "")
                        stations=stations.replace("LBA", "")
			self.stations = stations.split(",")
			self.nstations = len(self.stations)
			self.ncorestations = stations.count("CS")
			self.nremotestations = stations.count("RS")
		except: pass

		if not cmdline.opts.is_auto:
			# checking "locations" keywords first as in the new parset files (as from Jan 27, 2012) "mountpoints" can give wrong values
			res=[ii for ii in plines if "Output_Beamformed.locations" in ii]
			try:
				self.assigned_nodeslist=res[0].split("=", 1)[-1].strip().lstrip("[").rstrip("]").split(",")
				self.assigned_nodeslist=[n.split(":")[0] for n in self.assigned_nodeslist]
			except: pass
		else:
			#
			# if we are on CEP4 in auto mode, then in the pipeline parset, all "host" values are the same = "CEP4"
			# so we still need to use manual mapping
			if si.cluster_headnode[:5] == "head0" or si.cluster_headnode[:3] == "cpu": # CEP4
				self.assigned_nodeslist = si.locus_nodes[:]
			else:
				self.assigned_nodeslist = [x.host for x in self.input_data['incoherent'].data]
				self.assigned_nodeslist.extend([x.host for x in self.input_data['coherent'].data])

		self.assigned_nodeslist = list(set(self.assigned_nodeslist))

		# checking if all nodes in assigned_nodeslist are in alive
		if not cmdline.opts.is_locate_rawdata and len(self.assigned_nodeslist) > 0:
			if not (cmdline.opts.is_slurm and cmdline.opts.is_globalfs):
				nodes_unavail=list(set(self.assigned_nodeslist)-set(si.alive_nodes).intersection(set(self.assigned_nodeslist)))
				self.nodeslist=list(set(si.alive_nodes).intersection(set(self.assigned_nodeslist)))
				if len(nodes_unavail) > 0:
					msg="Warning! Some raw data are on nodes that are not available [%d]: %s" % (len(nodes_unavail), ", ".join(nodes_unavail))
					if log != None: log.warning(msg)
					else: print msg
			else:
				self.nodeslist=self.assigned_nodeslist

                # check if online coherent dedispersion (OCD) was used
                res=[ii for ii in plines if "OLAP.coherentDedisperseChannels" in ii]
		try:
			if res[0].split("=", 1)[-1].strip().lower()[0] == 't':
				self.OCD=True
		except: pass

		# getting info about the Type of the data (CV, CS, IS, FE, Imaging, etc.)
                res=[ii for ii in plines if "Output_CoherentStokes.enabled" in ii]
		try:
			if res[0].split("=", 1)[-1].strip().lower()[0] == 't':
				self.CS = True
                                self.CV = False
                                res=[ii for ii in plines if "OLAP.CNProc_CoherentStokes.which" in ii]
				self.stokesCS=res[0].split("=", 1)[-1].strip()
                                # in the transition phase there were some parset with just XY
                                # this means just 2 files, one for X, another Y
                                # now is always XXYY, i.e. 4 files get written
                                if self.stokesCS == "XXYY" or self.stokesCS == "XY":
					self.CV = True
                                        self.CS = False
		except: pass

	        # check if data are incoherent stokes data
		res=[ii for ii in plines if "Output_IncoherentStokes.enabled" in ii]
		try:
                	if res[0].split("=", 1)[-1].strip().lower()[0] == 't':
                        	self.IS = True
				res=[ii for ii in plines if "OLAP.CNProc_IncoherentStokes.which" in ii]
				self.stokesIS=res[0].split("=", 1)[-1].strip()
		except: pass

		# at ~05.07.2012 the logic in the Parset files has changed, and the flag Output_CoherentStokes.enabled = True
		# ONLY for CS data and not for CV data.  For CV data one needs to check Output_Beamformed.enabled flag
		if self.IS == False and self.CS == False:
			# checking the Output_Beamformed flag
			res=[ii for ii in plines if "Output_Beamformed.enabled" in ii]
			try:
                		if res[0].split("=", 1)[-1].strip().lower()[0] == 't':
                        		self.CV = True
					res=[ii for ii in plines if "OLAP.CNProc_CoherentStokes.which" in ii]
					self.stokesCS=res[0].split("=", 1)[-1].strip()
			except: pass
			else:
					if self.stokesCS != "XXYY" and self.stokesCS != "XY":
						msg="Wrong Stokes setting '%s' for CV observation for ObsID = %s" % (self.stokesCS, self.id)
						if log != None: log.error(msg)
						else: print msg
						raise Exception

	        # check if data are imaging
		res=[ii for ii in plines if "Output_Correlated.enabled" in ii]
		try:
                	if res[0].split("=", 1)[-1].strip().lower()[0] == 't':
				self.IM = True
		except: pass

	        # check if data are fly's eye mode data
		res=[ii for ii in plines if "PencilInfo.flysEye" in ii]
		try:
                	if res[0].split("=", 1)[-1].strip().lower()[0] == 't':
				self.FE = True
		except: pass

                # if Stokes are still undetermined (e.g. if obs is IM), then rereading default stokes for CS
                if self.stokesCS == "":
			res=[ii for ii in plines if "OLAP.CNProc_CoherentStokes.which" in ii]
			try:
				self.stokesCS=res[0].split("=", 1)[-1].strip()
			except: pass
                if self.stokesIS == "":
			res=[ii for ii in plines if "OLAP.CNProc_IncoherentStokes.which" in ii]
			try:
                                self.stokesIS=res[0].split("=", 1)[-1].strip()
			except: pass

		# Getting the list of subbands
		res=[ii for ii in plines if "Observation.subbandList" in ii]
		try:
			# getting range of subbands
			self.subbandList=res[0].split("=", 1)[-1].strip().lstrip("[").rstrip("]")
			# parsing the string with subband ranges to get list of subbands
			self.subbands = self.getSubbands(self.subbandList)
			# getting total number of Subbands
			self.nrSubbands = len(self.subbands)
		except: pass

		# in new parset files (after Jan 27, 2012) there are new keywords for number of
		# chans per subband and this number can be different for IS and CS
		res=[ii for ii in plines if "OLAP.CNProc_IncoherentStokes.channelsPerSubband" in ii]
		try:
			# getting number of IS channels
			self.nrChanPerSubIS=int(res[0].split("=", 1)[-1].strip())
		except: pass

		res=[ii for ii in plines if "OLAP.CNProc_CoherentStokes.channelsPerSubband" in ii]
		try:
			# getting number of CS channels
			self.nrChanPerSubCS=int(res[0].split("=", 1)[-1].strip())
		except: pass

		# getting the number of subbands per File for IS
		res=[ii for ii in plines if "OLAP.CNProc_IncoherentStokes.subbandsPerFile" in ii]
		try:
			# getting number of subbands per file for IS
			self.nrSubsPerFileIS=int(res[0].split("=", 1)[-1].strip())
		except: pass

		# getting the number of subbands per File for CS
		res=[ii for ii in plines if "OLAP.CNProc_CoherentStokes.subbandsPerFile" in ii]
		try:
			# getting number of subbands per file for CS
			self.nrSubsPerFileCS=int(res[0].split("=", 1)[-1].strip())
		except: pass

		# calculate the number of frequency splits for CS/CV and IS
		self.nsplitsCS = int(self.nrSubbands / self.nrSubsPerFileCS)
		if self.nrSubbands % self.nrSubsPerFileCS != 0: self.nsplitsCS += 1
		self.nsplitsIS = int(self.nrSubbands / self.nrSubsPerFileIS)
		if self.nrSubbands % self.nrSubsPerFileIS != 0: self.nsplitsIS += 1

		# changing cmdline split-related options in the _copy_ of Cmdline class
		# for CS/CV
		if cmdline.opts.first_freq_splitCS >= self.nsplitsCS: cmdline.opts.first_freq_splitCS = 0
		if cmdline.opts.nsplitsCS == -1: cmdline.opts.nsplitsCS = self.nsplitsCS
		if cmdline.opts.first_freq_splitCS + cmdline.opts.nsplitsCS > self.nsplitsCS:
			cmdline.opts.nsplitsCS -= (cmdline.opts.first_freq_splitCS + cmdline.opts.nsplitsCS - self.nsplitsCS)
		# for IS
		if cmdline.opts.first_freq_splitIS >= self.nsplitsIS: cmdline.opts.first_freq_splitIS = 0
		if cmdline.opts.nsplitsIS == -1: cmdline.opts.nsplitsIS = self.nsplitsIS
		if cmdline.opts.first_freq_splitIS + cmdline.opts.nsplitsIS > self.nsplitsIS:
			cmdline.opts.nsplitsIS -= (cmdline.opts.first_freq_splitIS + cmdline.opts.nsplitsIS - self.nsplitsIS)

		# Getting the sample clock
		res=[ii for ii in plines if "Observation.sampleClock" in ii]
		try:
			# getting the clock
			self.sampleClock=int(res[0].split("=", 1)[-1].strip())
		except: pass

		if self.sampleClock == 0: # if keyword 'Observation.sampleClock' is missing in the parset file
			res=[ii for ii in plines if "Observation.clockMode" in ii]
			try:
				# getting the clock
				self.sampleClock=int(res[0].split("=", 1)[-1].strip().split("Clock")[1])
			except: pass

		# Getting width of the subband (in kHz)
		res=[ii for ii in plines if "Observation.subbandWidth" in ii]
		try:
			self.subbandWidth=float(res[0].split("=", 1)[-1].strip())
		except: pass

		if self.subbandWidth == 0 and self.sampleClock != 0:
			self.subbandWidth = ( ( self.sampleClock / 2. ) / 512. ) * 1000.
		
		# getting timeIntegrationFactors for IS and calculating sampling intervals
		res=[ii for ii in plines if "OLAP.CNProc_IncoherentStokes.timeIntegrationFactor" in ii]
		try:
			self.downsample_factorIS=int(res[0].split("=", 1)[-1].strip())
			if self.downsample_factorIS != 0 and self.sampleClock != 0 and self.nrChanPerSubIS != 0:
				self.samplingIS = self.downsample_factorIS / ((self.sampleClock * 1000. * 1000. / 1024.) / self.nrChanPerSubIS) * 1000.
		except: pass

		# getting timeIntegrationFactors for CS and calculating sampling intervals
		res=[ii for ii in plines if "OLAP.CNProc_CoherentStokes.timeIntegrationFactor" in ii]
		try:
			self.downsample_factorCS=int(res[0].split("=", 1)[-1].strip())
			if self.downsample_factorCS != 0 and self.sampleClock != 0 and self.nrChanPerSubCS != 0:
				self.samplingCS = self.downsample_factorCS / ((self.sampleClock * 1000. * 1000. / 1024.) / self.nrChanPerSubCS) * 1000.
		except: pass

		# Calculating the total BW (in MHz)
		if self.nrSubbands != 0 and self.subbandWidth != 0:
			self.bw = self.subbandWidth * self.nrSubbands / 1000.

		# Calculating the central freq (in MHz)
		if self.band != "":
			try:
				lower_band_freq = int(self.band.split("_")[0])
				if self.sampleClock == 200:
					if lower_band_freq > 200:
						self.lower_band_edge = 200.0
					elif lower_band_freq < 200 and lower_band_freq > 100:
						self.lower_band_edge = 100.0
					else: self.lower_band_edge = 0.0
				if self.sampleClock == 160:
					if lower_band_freq >= 160:
						self.lower_band_edge = 160
					elif lower_band_freq < 160 and lower_band_freq >= 80:
						self.lower_band_edge = 80
					else: self.lower_band_edge = 0

				if len(self.subbands) > 0 and self.subbandWidth != 0 and (self.nrChanPerSubIS != 0 or self.nrChanPerSubCS != 0):
					try:
						# CS has a priority
						if self.nrChanPerSubCS != 0: nchanpersub = self.nrChanPerSubCS
						else: nchanpersub = self.nrChanPerSubIS
						if nchanpersub > 1: # it means that we ran 2nd polyphase
							lofreq = self.lower_band_edge + (self.subbandWidth / 1000.) * self.subbands[0] - 0.5 * \
									(self.subbandWidth / 1000.) - 0.5 * (self.subbandWidth / 1000. / nchanpersub)
							hifreq = self.lower_band_edge + (self.subbandWidth / 1000.) * (self.subbands[-1] + 1) - 0.5 * \
									(self.subbandWidth / 1000.) - 0.5 * (self.subbandWidth / 1000. / nchanpersub)
						else:
							lofreq = self.lower_band_edge + (self.subbandWidth / 1000.) * self.subbands[0] - 0.5 * (self.subbandWidth / 1000.)
							hifreq = self.lower_band_edge + (self.subbandWidth / 1000.) * (self.subbands[-1] + 1) - 0.5 * (self.subbandWidth / 1000.)
						self.freq_extent = hifreq - lofreq
						self.cfreq = lofreq + self.freq_extent/2.
					except: pass
			except: pass

		# Getting number of Station beams
		res=[ii for ii in plines if "Observation.nrBeams" in ii]
		try:
			self.nrBeams=int(res[0].split("=", 1)[-1].strip())
		except: pass

		# initializing SAP beams objects and making the list of SAP beams
		for sid in xrange(self.nrBeams):
			sap=SAPBeam(sid)
			sap.update(plines, si, cmdline, self, log, wrapper)
			self.saps.append(sap)

		# calculating the total number of TABs in all SAPs
		self.nrTiedArrayBeams = sum([sap.nrTiedArrayBeams for sap in self.saps])

		# updating the nodeslist if we are looking up for raw data in all alive nodes rather than using the nodeslist from Parset
		if cmdline.opts.is_locate_rawdata:
			self.nodeslist=[]
			for sap in self.saps:
				for tab in sap.tabs:
					self.nodeslist.extend(tab.location)
			self.nodeslist=list(set(self.nodeslist))


	# prints info about the observation
	def print_info(self, cmdline, log=None):
		"""
		prints info about the observation
		"""
		if log != None:
			log.info("\n================================================================")
			log.info("ObsID: %s   Project: %s   PI: %s" % (self.id, self.project, self.projectPI))	
			if not cmdline.opts.is_cobalt:
				log.info("Parset: %s" % (self.parset))
			log.info("Start UTC: %s %s  Duration: %s" % \
				(self.startdate, self.starttime, self.duration>=3600. and "%.1fh" % \
				(float(self.duration/3600.)) or "%.1fm" % (float(self.duration/60.))))
			mode=""
			if self.FE: 
				if not self.CS and not self.CV:
					mode+="FE (" + self.stokesCS + ") "
				else:
					mode+="FE "
			if self.IM: mode+="Im "
			if self.IS: mode+="IS (" + self.stokesIS + ") "
			if self.CS: mode+="CS (" + self.stokesCS + ") "
			if self.CV: mode+="CV "
			log.info("%s   Band: %s   Mode: %s   OCD: %s" % (self.antenna_config, self.band, mode, self.OCD and "yes" or "no"))
			log.info("#stations: %d [%dCS, %dRS]" % (self.nstations, self.ncorestations, self.nremotestations))
			log.info("Clock: %d MHz   Fc: %g MHz   BW: %g MHz" % (self.sampleClock, self.cfreq, self.bw))
			log.info("#subbands: %d [%s]   SubWidth: %g kHz" % (self.nrSubbands, self.subbandList, self.subbandWidth))
			if self.nrChanPerSubIS == self.nrChanPerSubCS or self.nrChanPerSubIS == 0 or self.nrChanPerSubCS == 0:
				nchanspersub = (self.nrChanPerSubIS != 0 and str(self.nrChanPerSubIS) or str(self.nrChanPerSubCS))
			elif self.IS and not self.CS and not self.CV:
				nchanspersub = str(self.nrChanPerSubIS)
			elif not self.IS and (self.CS or self.CV):
				nchanspersub = str(self.nrChanPerSubCS)
			else: nchanspersub = "%d (IS), %d (%s)" % (self.nrChanPerSubIS, self.nrChanPerSubCS, self.CS and "CS" or "CV")
			if self.downsample_factorIS == self.downsample_factorCS or self.downsample_factorIS == 0 or self.downsample_factorCS == 0:
				dfactor = (self.downsample_factorIS != 0 and str(self.downsample_factorIS) or str(self.downsample_factorCS))
			elif self.IS and not self.CS and not self.CV:
				dfactor = str(self.downsample_factorIS)
			elif not self.IS and (self.CS or self.CV):
				dfactor = str(self.downsample_factorCS)
			else: dfactor = "%d (IS), %d (%s)" % (self.downsample_factorIS, self.downsample_factorCS, self.CS and "CS" or "CV")
			if self.samplingIS == self.samplingCS or self.samplingIS == 0.0 or self.samplingCS == 0.0:
				sampling = (self.samplingIS != 0.0 and str(self.samplingIS) or str(self.samplingCS)) + " ms"
			elif self.IS and not self.CS and not self.CV:
				sampling = str(self.samplingIS) + " ms"
			elif not self.IS and (self.CS or self.CV):
				sampling = str(self.samplingCS) + " ms"
			else: sampling = "%g ms (IS), %g ms (%s)" % (self.samplingIS, self.samplingCS, self.CS and "CS" or "CV")
			log.info("#chans/sub: %s   Downsample Factor: %s" % (nchanspersub, dfactor))
			log.info("Sampling: %s" % (sampling))
			if self.nrBeams > 1:
				log.info("#SAPs: %d" % (self.nrBeams))
				for sap in self.saps:
					log.info("%d Target: %s   #TABs: %d%s" % (sap.sapid, sap.source, sap.nrTiedArrayBeams, \
						sap.nrRings > 0 and " #Rings: %d RingSize: %g deg" % (sap.nrRings, sap.ringSize) or ""))
			else:
				log.info("#SAPs: %d   Target: %s   #TABs: %d%s" % (self.nrBeams, self.saps[0].source, self.saps[0].nrTiedArrayBeams, \
					self.saps[0].nrRings > 0 and " #Rings: %d RingSize: %g deg" % (self.saps[0].nrRings, self.saps[0].ringSize) or ""))
			log.info("================================================================")

	# generates pseudo parset file to be used by other scripts, like
	# plot_LOFAR_TA_multibeam*.py, RRAT_heatmap*.py
	# file name is given by pseudoparset parameters (full name)
	def pseudo_parset_generator(self, pseudoparset):
		# if this file already exists, do nothing
		if os.path.exists(pseudoparset): return
		output_line="#"
		output_line+="\nObservation.startTime = '%s %s'" % (self.startdate, self.starttime) 
		output_line+="\nOLAP.storageStationNames = [%s]" % (",".join(["%s%s%s" % (s.split("/")[0], self.antenna, s.split("/")[-1]) for s in self.stations]))
		for sap in self.saps:
			output_line+="\nObservation.Beam[%d].subbandList = [%s]" % (sap.sapid, sap.subbandList)
			output_line+="\nObservation.Beam[%d].angle1 = %s" % (sap.sapid, str(sap.rarad))
			output_line+="\nObservation.Beam[%d].angle2 = %s" % (sap.sapid, str(sap.decrad))
			for tab in sap.tabs:
				output_line+="\nObservation.Beam[%d].TiedArrayBeam[%d].angle1 = %s" % (sap.sapid, tab.tabid, str(tab.raoffset))
				output_line+="\nObservation.Beam[%d].TiedArrayBeam[%d].angle2 = %s" % (sap.sapid, tab.tabid, str(tab.decoffset))
				output_line+="\nObservation.Beam[%d].TiedArrayBeam[%d].coherent = %s" % (sap.sapid, tab.tabid, str(tab.is_coherent).lower())
		with open(pseudoparset, "w") as pseudo:
			pseudo.write(output_line)

	# update info based on .h5 files
	def update_from_H5 (self, si, cmdline, log=None, wrapper=None):

                # copying all .h5 files to log-directory
		if not cmdline.opts.is_auto:
			if not (cmdline.opts.is_slurm and cmdline.opts.is_globalfs):
				self.assigned_nodeslist = self.nodeslist = si.alive_nodes # for Cobalt we will always search
			else:
				self.assigned_nodeslist = self.nodeslist = si.cexec_nodes.keys()	
			if not cmdline.opts.is_globalfs:
				if len(self.nodeslist) > 0:
        	                	cexeclocus=si.cexec_nodes[self.nodeslist[0]]
	        	                if len(self.nodeslist) > 1:
        	        	                for s in self.nodeslist[1:]:
                	        	                cexeclocus += ",%s" % (si.cexec_nodes[s].split(":")[1])
					# first try to copy all *.h5 files from the processed directories if they exist. This is necessary in case when _not_all_ raw data were deleted
					# Then the pipeline would fail, because list of .h5 files is not empty but not full
					# We copy first from the processed directories in case some of the .h5 files have zero size (due to problems, no disk space, etc. - I had
					# this problem before). Then, these .h5 files will be overwritten by the original ones, if they still exist
                                        if si.cluster_headnode[:4] == "drag" or si.cluster_headnode[:3] == "drg":
                                            cmd="%s %s 'find /*/%s* -name \"%s_SAP*_B*_S*_P*_bf.h5\" -exec cp -f {} %s \;' | grep -v denied | grep -v such | grep -v match | grep -v xauth | grep -v connect | grep -v closed | grep -v information | grep -v missing | grep -v overwrite | egrep -v \'\\-\\-\\-\\-\\-\' | egrep -v \'\\*\\*\\*\\*\\*\'" % (si.cexeccmd, cexeclocus, si.processed_dir_prefix, self.id, si.get_logdir())
                                        else:
					    cmd="%s %s 'cp -f -L $(find /*/%s* -name \"%s_SAP*_B*_S*_P*_bf.h5\") %s' | grep -v denied | grep -v such | grep -v match | grep -v xauth | grep -v connect | grep -v closed | grep -v information | grep -v missing | grep -v overwrite | egrep -v \'\\-\\-\\-\\-\\-\' | egrep -v \'\\*\\*\\*\\*\\*\'" % (si.cexeccmd, cexeclocus, si.processed_dir_prefix, self.id, si.get_logdir())
        		                if cmdline.opts.is_debug: log.debug("cmd: %s" % cmd)
                		        os.system(cmd)
					# copy .h5 from the original directories
                                        if si.cluster_headnode[:4] == "drag" or si.cluster_headnode[:3] == "drg":
                                            cmd="%s %s 'find %s*/%s%s -name \"%s_SAP*_B*_S*_P*_bf.h5\" -exec cp -f {} %s \;' | grep -v denied | grep -v such | grep -v match | grep -v xauth | grep -v connect | grep -v closed | grep -v information | grep -v missing | grep -v overwrite | egrep -v \'\\-\\-\\-\\-\\-\' | egrep -v \'\\*\\*\\*\\*\\*\'" % (si.cexeccmd, cexeclocus, si.rawdir, si.rawdir_prefix_specificator, self.id, self.id, si.get_logdir())
                                        else:
					    cmd="%s %s 'cp -f -L $(find %s*/%s%s -name \"%s_SAP*_B*_S*_P*_bf.h5\") %s' | grep -v denied | grep -v such | grep -v match | grep -v xauth | grep -v connect | grep -v closed | grep -v information | grep -v missing | grep -v overwrite | egrep -v \'\\-\\-\\-\\-\\-\' | egrep -v \'\\*\\*\\*\\*\\*\'" % (si.cexeccmd, cexeclocus, si.rawdir, si.rawdir_prefix_specificator, self.id, self.id, si.get_logdir())
	        	                if cmdline.opts.is_debug: log.debug("cmd: %s" % cmd)
        	        	        os.system(cmd)
			else: # here we use global FS
				cmd="cp -f -L $(find %s*/%s%s -name \"%s_SAP*_B*_S*_P*_bf.h5\") %s > /dev/null" % (si.rawdir, si.rawdir_prefix_specificator, self.id, self.id, si.get_logdir())
        	                if cmdline.opts.is_debug: log.debug("cmd: %s" % cmd)
       	        	        os.system(cmd)
                else:
			#
			# if we are on CEP4 in auto mode, then in the pipeline parset, all "host" values are the same = "CEP4"
			# so we still need to use manual mapping
			if si.cluster_headnode[:5] == "head0" or si.cluster_headnode[:3] == "cpu": # CEP4
				self.nodeslist = si.locus_nodes[:]
				if cmdline.opts.is_debug: log.debug("copy CEP4:*.h5 %s" % (si.get_logdir()))
				infiles=[f.file for f in wrapper.input_data['incoherent'].data]
				infiles.extend([f.file for f in wrapper.input_data['coherent'].data])
				for source in infiles:
					try:
						subprocess.check_call(['cp', "%s" % (source), si.get_logdir()])
					except Exception:
						if cmdline.opts.is_debug: log.debug("failed")
						log.warning("Can't copy %s to %s. Probably data are missing." % (source, si.get_logdir()))
			else:
                        	self.nodeslist = []
				for source in wrapper.input_data['incoherent'].data:
					if cmdline.opts.is_debug: log.debug("copy %s:%s %s" % (source.host, source.file, si.get_logdir()))
					try:
						if not cmdline.opts.is_globalfs:
							subprocess.check_call(['scp', "%s:%s" % (source.host, source.file), si.get_logdir()])
						else:
							subprocess.check_call(['cp', "%s" % (source.file), si.get_logdir()])
						self.nodeslist.append(source.host)
					except Exception:
						if cmdline.opts.is_debug: log.debug("failed")
						log.warning("Can't copy %s:%s to %s. Probably data are missing." % (source.host, source.file, si.get_logdir()))
				for source in wrapper.input_data['coherent'].data:
					if cmdline.opts.is_debug: log.debug("copy %s:%s %s" % (source.host, source.file, si.get_logdir()))
					try:
						if not cmdline.opts.is_globalfs:
							subprocess.check_call(['scp', "%s:%s" % (source.host, source.file), si.get_logdir()])
						else:
							subprocess.check_call(['cp', "%s" % (source.file), si.get_logdir()])
						self.nodeslist.append(source.host)
					except Exception:
						if cmdline.opts.is_debug: log.debug("failed")
						log.warning("Can't copy %s:%s to %s. Probably data are missing." % (source.host, source.file, si.get_logdir()))
					
			self.nodeslist = list(set(self.nodeslist))
			self.assigned_nodeslist = self.nodeslist

		if len(self.nodeslist) > 0:
			# get list of .h5 files
			h5list=glob.glob("%s/%s*_SAP*_B*_S*_P*_bf.h5" % (si.get_logdir(), self.id))
			if len(h5list) == 0:
				msg="No .h5 files are found. Exiting..."
				if log != None: log.error(msg)
				else: print msg
				raise Exception
			# first, we read any 1 file to collect general info from the ROOT group
			sapid=int(h5list[0].split("/")[-1].split("_SAP")[-1].split("_")[0])
			tabid=int(h5list[0].split("/")[-1].split("_B")[-1].split("_")[0])
			f5 = h5py.File(h5list[0], 'r')
			self.project = f5.attrs['PROJECT_ID']
			self.projectPI = f5.attrs['PROJECT_PI']
			self.startdate = f5.attrs['OBSERVATION_START_UTC']
			self.starttime = self.startdate.split("T")[1]
			if not self.starttime[-1].isdigit(): # for whatever reason there can be "Z" in the end of the time string
				self.starttime = self.starttime[:-1]
			self.startdate  = self.startdate.split("T")[0]
			self.duration=f5.attrs['TOTAL_INTEGRATION_TIME']
			self.bandFilter = f5.attrs['FILTER_SELECTION']
			self.antenna = self.bandFilter.split("_")[0]
			self.band = self.bandFilter.split("A_")[-1]
			self.antenna_config = f5.attrs['ANTENNA_SET']

			# changing cmdline FWHM-related options in the _copy_ of Cmdline class
        	        if cmdline.opts.fwhm_CS < 0.0:
                	        if self.antenna == "HBA": cmdline.opts.fwhm_CS = si.fwhm_hba
                        	if self.antenna == "LBA": cmdline.opts.fwhm_CS = si.fwhm_lba
	                if cmdline.opts.fwhm_IS < 0.0:
        	                if self.antenna == "HBA": cmdline.opts.fwhm_IS = si.fov_hba
                	        if self.antenna == "LBA": cmdline.opts.fwhm_IS = si.fov_lba

			stations=f5.attrs['OBSERVATION_STATIONS_LIST']
			self.stations = [s[0:5] for s in stations]
			self.nstations = f5.attrs['OBSERVATION_NOF_STATIONS']
			self.ncorestations = len([s for s in stations if s[0:2] == "CS"])
			self.nremotestations = len([s for s in stations if s[0:2] == "RS"])
			self.sampleClock = int(f5.attrs['CLOCK_FREQUENCY'])
			# based on the clock and used filter, we calculate the lowest frequency of the band
			lower_band_freq = int(self.band.split("_")[0])
			if self.sampleClock == 200:
				if lower_band_freq > 200:
					self.lower_band_edge = 200.0
				elif lower_band_freq < 200 and lower_band_freq > 100:
					self.lower_band_edge = 100.0
				else: self.lower_band_edge = 0.0
			if self.sampleClock == 160:
				if lower_band_freq >= 160:
					self.lower_band_edge = 160
				elif lower_band_freq < 160 and lower_band_freq >= 80:
					self.lower_band_edge = 80
				else: self.lower_band_edge = 0

			self.bw = f5.attrs['BANDWIDTH']
			self.subbandWidth = f5['SUB_ARRAY_POINTING_%03d' % (sapid)]['BEAM_%03d' % (tabid)].attrs['SUBBAND_WIDTH'] / 1000. # in kHz
			self.nrBeams = f5.attrs['OBSERVATION_NOF_SUB_ARRAY_POINTINGS']

			# initializing SAP beams objects and making the list of SAP beams
			for sid in xrange(self.nrBeams):
				try:
					sap=SAPBeam(sid)
					sap.update_from_H5(si, cmdline, self, log, wrapper)
					self.saps.append(sap)
				except Exception: 
					msg="Error in getting info for the SAP=%d. Probably data are missing." % (sid)
					if log != None: log.warning(msg)
					else: print msg

			# calculating the total number of TABs in all SAPs
			self.nrTiedArrayBeams = sum([sap.nrTiedArrayBeams for sap in self.saps])

			f5.close()

			# updating the nodeslist as it's possible we are only processing few splits
			if not cmdline.opts.is_auto:
				self.nodeslist=[]
				for sap in self.saps:
					for tab in sap.tabs:
						self.nodeslist.extend(tab.location)

			# getting number of subs/file for coherent beams
			for sap in self.saps:
				for tab in sap.tabs:
					if tab.is_coherent:
						f5 = h5py.File("%s/%s_SAP%03d_B%03d_S0_P000_bf.h5" % (si.get_logdir(), self.id, sap.sapid, tab.tabid), 'r')
						self.nrSubsPerFileCS = f5['SUB_ARRAY_POINTING_%03d' % (sap.sapid)]['BEAM_%03d' % (tab.tabid)]['STOKES_0'].attrs['NOF_SUBBANDS']
						self.nrChanPerSubCS = f5['SUB_ARRAY_POINTING_%03d' % (sap.sapid)]['BEAM_%03d' % (tab.tabid)].attrs['CHANNELS_PER_SUBBAND']
						if str(f5['SUB_ARRAY_POINTING_%03d' % (sap.sapid)]['BEAM_%03d' % (tab.tabid)].attrs['DEDISPERSION']).lower() == "coherent":
							self.OCD = True
						nof_stokes = f5['SUB_ARRAY_POINTING_%03d' % (sap.sapid)]['BEAM_%03d' % (tab.tabid)].attrs['OBSERVATION_NOF_STOKES']
						if nof_stokes == 1: 
							self.stokesCS = "I"
							self.CS = True
						else:
							if bool(f5['SUB_ARRAY_POINTING_%03d' % (sap.sapid)]['BEAM_%03d' % (tab.tabid)].attrs['COMPLEX_VOLTAGE']):
								self.CV = True
								self.stokesCS = "XXYY"
							else:
								self.CS = True
								self.stokesCS = "IQUV"
						self.samplingCS = f5['SUB_ARRAY_POINTING_%03d' % (sap.sapid)]['BEAM_%03d' % (tab.tabid)].attrs['SAMPLING_TIME'] * 1000. # in ms
						self.downsample_factorCS = int(self.samplingCS * ((self.sampleClock * 1000. / 1024.) / self.nrChanPerSubCS) + 0.5)
						f5.close()
						break
				else: continue
				break

			# getting number of subs/file for incoherent beams
			for sap in self.saps:
				for tab in sap.tabs:
					if not tab.is_coherent:
						self.IS = True
						f5 = h5py.File("%s/%s_SAP%03d_B%03d_S0_P000_bf.h5" % (si.get_logdir(), self.id, sap.sapid, tab.tabid), 'r')
						self.nrSubsPerFileIS = f5['SUB_ARRAY_POINTING_%03d' % (sap.sapid)]['BEAM_%03d' % (tab.tabid)]['STOKES_0'].attrs['NOF_SUBBANDS']
						self.nrChanPerSubIS = f5['SUB_ARRAY_POINTING_%03d' % (sap.sapid)]['BEAM_%03d' % (tab.tabid)].attrs['CHANNELS_PER_SUBBAND']
						nof_stokes = f5['SUB_ARRAY_POINTING_%03d' % (sap.sapid)]['BEAM_%03d' % (tab.tabid)].attrs['OBSERVATION_NOF_STOKES']
						if nof_stokes == 1: self.stokesIS = "I"
						else: self.stokesIS = "IQUV"
						self.samplingIS = f5['SUB_ARRAY_POINTING_%03d' % (sap.sapid)]['BEAM_%03d' % (tab.tabid)].attrs['SAMPLING_TIME'] * 1000. # in ms
						self.downsample_factorIS = int(self.samplingIS * ((self.sampleClock * 1000. / 1024.) / self.nrChanPerSubIS) + 0.5)
						f5.close()
						break
				else: continue
				break
			if self.nrSubsPerFileIS == 0: self.nrSubsPerFileIS = self.nrSubsPerFileCS
			if self.nrSubsPerFileCS == 0: self.nrSubsPerFileCS = self.nrSubsPerFileIS

			# merging subbands lists from all SAPs
			for sap in self.saps:
				self.subbands.extend(sap.subbands)
			self.subbands=list(set(self.subbands)) # getting list of unique subbands
			self.subbands.sort()
			self.subbandList="%d..%d" % (self.subbands[0], self.subbands[-1])
			self.nrSubbands = len(self.subbands)
			if self.nrSubbands > 0 and self.subbandWidth != 0 and (self.nrChanPerSubIS != 0 or self.nrChanPerSubCS != 0):
				try:
					# CS has a priority
					if self.nrChanPerSubCS != 0: nchanpersub = self.nrChanPerSubCS
					else: nchanpersub = self.nrChanPerSubIS
					if nchanpersub > 1: # it means that we ran 2nd polyphase
						lofreq = self.lower_band_edge + (self.subbandWidth / 1000.) * self.subbands[0] - 0.5 * \
								(self.subbandWidth / 1000.) - 0.5 * (self.subbandWidth / 1000. / nchanpersub)
						hifreq = self.lower_band_edge + (self.subbandWidth / 1000.) * (self.subbands[-1] + 1) - 0.5 * \
								(self.subbandWidth / 1000.) - 0.5 * (self.subbandWidth / 1000. / nchanpersub)
					else:
						lofreq = self.lower_band_edge + (self.subbandWidth / 1000.) * self.subbands[0] - 0.5 * (self.subbandWidth / 1000.)
						hifreq = self.lower_band_edge + (self.subbandWidth / 1000.) * (self.subbands[-1] + 1) - 0.5 * (self.subbandWidth / 1000.)
					self.freq_extent = hifreq - lofreq
					self.cfreq = lofreq + self.freq_extent/2.
				except: pass

			# changing cmdline split-related options in the _copy_ of Cmdline class
			# calculate the number of frequency splits for CS/CV and IS
			self.nsplitsCS = int(self.nrSubbands / self.nrSubsPerFileCS)
			if self.nrSubbands % self.nrSubsPerFileCS != 0: self.nsplitsCS += 1
			self.nsplitsIS = int(self.nrSubbands / self.nrSubsPerFileIS)
			if self.nrSubbands % self.nrSubsPerFileIS != 0: self.nsplitsIS += 1

			# changing cmdline split-related options in the _copy_ of Cmdline class
			# for CS/CV
			if cmdline.opts.first_freq_splitCS >= self.nsplitsCS: cmdline.opts.first_freq_splitCS = 0
			if cmdline.opts.nsplitsCS == -1: cmdline.opts.nsplitsCS = self.nsplitsCS
			if cmdline.opts.first_freq_splitCS + cmdline.opts.nsplitsCS > self.nsplitsCS:
				cmdline.opts.nsplitsCS -= (cmdline.opts.first_freq_splitCS + cmdline.opts.nsplitsCS - self.nsplitsCS)
			# for IS
			if cmdline.opts.first_freq_splitIS >= self.nsplitsIS: cmdline.opts.first_freq_splitIS = 0
			if cmdline.opts.nsplitsIS == -1: cmdline.opts.nsplitsIS = self.nsplitsIS
			if cmdline.opts.first_freq_splitIS + cmdline.opts.nsplitsIS > self.nsplitsIS:
				cmdline.opts.nsplitsIS -= (cmdline.opts.first_freq_splitIS + cmdline.opts.nsplitsIS - self.nsplitsIS)

			# removing .h5 files in the end from the log directory
			cmd="rm -f %s" % (" ".join(h5list))
			os.system(cmd)
		else:
			msg="All locus nodes are down! Can't get any of .h5 files to get a metadata!"
			if log != None: log.error(msg)
			else: print msg
			raise Exception
