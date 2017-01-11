###################################################################
#
# Class FeedbackUnit - collects feedback info for every output
#                      tarball file from the processing + log file
#

import os, sys, os.path, tarfile

# class for keeping feedback info for processed data type (CS, IS, or CV)
class FeedbackUnit:
	def __init__(self, index, node, path, sapid=-1, tabid=-1, stokes_index=0, part=0):
		self.index = index
		self.node = node
		self.path = path
		self.sapid = sapid
		self.tabid = tabid
		self.stokes_index = stokes_index
		self.part = part
		self.fbfile = "" # output feedback filename for the current unit

		self._data = {
		 'fileFormat' : "PULP",
		 'filename' : "",
		 'location' : node + ":" + self.path,
		 'percentageWritten' : 0,
		 'size' : 0,
		 'datatype' : "", # CoherentStokes, IncoherentStoks, SummaryCoherentStokes, ...
		 'fileContent' : []
		}

	# function to update tarball filename and related fields
	def filename_update(self, filename):
		self._data['filename'] = filename.split("/")[-1]
		# file size in bytes
		try:
			self._data['size'] = os.path.getsize(self.path + "/" + self._data['filename'])
		except Exception: self._data['size'] = 0
		# getting contents of the tarball file
		try:
			tar=tarfile.open(self.path + "/" + self._data['filename'])
			self._data['fileContent']=[ii.name for ii in tar.getmembers()]
		except Exception: pass

	# function to update the feedback unit when the actual file exists
	def update(self, filename, fbfile, code, obs, is_summary=False):
		self.fbfile = fbfile
		self._data['filename'] = filename.split("/")[-1]
		# file size in bytes
		try:
			self._data['size'] = os.path.getsize(self.path + "/" + self._data['filename'])
		except Exception: self._data['size'] = 0

		# updating the data type
		if code == "CS": self._data['datatype'] = "CoherentStokes"
		elif code == "IS": self._data['datatype'] = "IncoherentStokes"
		elif code == "CV": self._data['datatype'] = "ComplexVoltages"
		else: self._data['datatype'] == "Unknown"
		if is_summary and self._data['datatype'] != "Unknown":
			self._data['datatype'] = "Summary" + self._data['datatype']

		# getting contents of the tarball file
		try:
			tar=tarfile.open(self.path + "/" + self._data['filename'])
			self._data['fileContent']=[ii.name for ii in tar.getmembers()]
		except Exception: pass

		# filling the Beam info
		if not is_summary:
			for ss in obs.saps:
				if ss.sapid == self.sapid:
					sap = ss
					for tt in ss.tabs:
						if tt.tabid == self.tabid:
							tab = tt
							break
			if tab.specificationType == "flyseye":
				self._beam = "FlysEyeBeam"
				self._sampling = obs.samplingCS / 1000.
				self._nchans = obs.nrChanPerSubCS
				if obs.stokesCS[0] == 'X':
					self._stokes = ['Xre', 'Xim', 'Yre', 'Yim']
				else:
					self._stokes = list(list(obs.stokesCS)[self.stokes_index])
			else:
				if code == "CV" or code == "CS": 
					self._beam = "CoherentStokesBeam"
					self._sampling = obs.samplingCS / 1000.
					self._nchans = obs.nrChanPerSubCS
					if obs.stokesCS[0] == 'X':
						self._stokes = ['Xre', 'Xim', 'Yre', 'Yim']
					else:
						self._stokes = list(list(obs.stokesCS)[self.stokes_index])
				if code == "IS": 
					self._beam = "IncoherentStokesBeam"
					self._sampling = obs.samplingIS / 1000.
					self._nchans = obs.nrChanPerSubIS
					if obs.stokesIS[0] == 'X':
						self._stokes = ['Xre', 'Xim', 'Yre', 'Yim']
					else:
						self._stokes = list(list(obs.stokesIS)[self.stokes_index])
			self._data["%s.SAP" % self._beam] = self.sapid
			self._data["%s.TAB" % self._beam] = self.tabid
			self._data["%s.samplingTime" % self._beam] = self._sampling
			self._data["%s.dispersionMeasure" % self._beam] = tab.DM
			if code == "CV" or code == "CS": subs_file = obs.nrSubsPerFileCS
			if code == "IS": subs_file = obs.nrSubsPerFileIS
			nsplits = int(len(sap.subbands) / subs_file)
               		if len(sap.subbands) % subs_file != 0:
				nsplits += 1
			if nsplits == 1:
				subs = sap.subbands[:]
			else:
				if self.part == nsplits - 1:
					subs = sap.subbands[self.part*subs_file:]
				else:	
					subs = sap.subbands[self.part*subs_file:(self.part+1)*subs_file]
			self._data["%s.nrSubbands" % self._beam] = len(subs)
			self._data["%s.stationSubbands" % self._beam] = "[" + ",".join(["%d" % (sub) for sub in subs]) + "]"
			self._data["%s.centralFrequencies" % self._beam] = \
				"[" + ",".join(["%f" % ((obs.lower_band_edge * 1000. + sub * obs.subbandWidth) * 1000.) for sub in subs]) + "]" # in Hz
			self._data["%s.channelWidth" % self._beam] = obs.subbandWidth / self._nchans
			self._data["%s.channelsPerSubband" % self._beam] = self._nchans
			self._data["%s.stokes" % self._beam] = "[" + ",".join(["'%s'" % (s) for s in self._stokes]) + "]"
			# writing Station info for FE
			if tab.specificationType == "flyseye":
				self._data["%s.stationName" % self._beam] = tab.stationList[0][0:5]
				self._data["%s.antennaFieldName" % self._beam] = tab.stationList[0][5:]
			else: # writing coordinates for TABs (for CS/CV only)
				if code == "CV" or code == "CS":
					self._data["%s.Pointing.equinox" % self._beam] = "J2000"
					self._data["%s.Pointing.coordType" % self._beam] = "RA-DEC"
					self._data["%s.Pointing.angle1" % self._beam] = tab.rarad
					self._data["%s.Pointing.angle2" % self._beam] = tab.decrad
					self._data["%s.Offset.equinox" % self._beam] = "J2000"
					self._data["%s.Offset.coordType" % self._beam] = "RA-DEC"
					self._data["%s.Offset.angle1" % self._beam] = tab.rarad - sap.rarad
					self._data["%s.Offset.angle2" % self._beam] = tab.decrad - sap.decrad

	# function to flush the feedback info to the file
	def flush(self, how_much_is_written, cep2, is_summary, is_log=False):
		self._data['percentageWritten'] = how_much_is_written
		with open(self.fbfile, 'w') as fb:
			# if it is PULP log-file, then just write 
			if is_log: fb.write("%s=%s/%s\n" % (cep2.feedback_log_prefix, self._data['location'], self._data['filename']))
			else:
				for key, value in self._data.items():
					fb.write("%s[%d].%s=%s\n" % (cep2.feedback_prefix, self.index, key, value))
