import sys
import math
#import matplotlib.pyplot as plt
from scipy import stats
import numpy as np
#import statsmodels.api as sm
import argparse
import pysam


class Arguments:
	def __init__(self):
		i = ""
		c = ""
		t = ""
		o = ""

	def parse(self):
		parser = argparse.ArgumentParser()
		parser.add_argument("-c", help = "bam file - treated sample", required = False, type = str)
		parser.add_argument("-t", help = "bam file - control sample", required = False, type = str)
		parser.add_argument("-i", help = "input file - counts", required = False, type = str)
		parser.add_argument("-o", help = "output file", required = True, type = str)
		parser.add_argument("-id", help = "id list", required = False, type = str)
		args = parser.parse_args()
		self.i = args.i
		self.c = args.c
		self.t = args.t
		self.o = args.o
		self.id = args.id

class Transcript:
	def __init__(self):
		self.stops_control = []
		self.stops_modification = []
		self.stops_control_norm = []
		self.norm_c = 0.0
		self.FC = []
		self.reactivity = []
		self.idt = ""
		self.FC_lim = 0
		self.length = 0

	def norm_mean(self):
		self.FC = []
		self.FC_pos = []
		#print self.idt
		#log2(FC)
		for i in range(0, len(self.stops_modification)):
			try:
				if self.stops_control[i] != 0 and self.stops_modification[i] != 0: 
					self.FC.append(math.log(float(self.stops_modification[i]) / self.stops_control[i], 2))
					self.FC_pos.append(i)
			except KeyError:
				pass
		#self.plot_scatter(self.stops_modification, self.stops_control, "_65_71_before_normalization")

	def norm_kernel_densities(self):
		density = stats.gaussian_kde(self.FC)
		den = density(self.FC)
		index = np.argmax(den)
		#print density
		self.norm_c = math.pow(2, self.FC[index])
		max_y = den[index]
		f_o = open("maxy.txt", 'w')
		f_o.write(str(max_y))
		print max_y
		#print self.norm_c
		for i in range(0, self.length):
			try:
				self.stops_control_norm.append(self.stops_control[i] * self.norm_c)
			except KeyError:
				pass
		#self.plot_scatter(self.stops_modification, self.stops_control_norm, "_65_71_normalized")
		#plt.figure()
		# plot histgram of sample
		#plt.hist(self.FC, bins=20, normed=1)
		# plot estimated density
		#plt.legend()
		#plt.show()

	def plot_scatter(self, a, b, t):
		plt.scatter(a, b)
		plt.savefig(str(self.idt) + t, facecolor='w', edgecolor='w')
		plt.close()

	def plot_hist(self):
		#histogram FC
		plt.hist(self.FC, bins = 30, normed=1, facecolor='green', alpha=0.5)
		plt.xlabel('log2(FC)')
		plt.title(self.idt)
		plt.savefig(str(self.idt) + "_FC_hist", facecolor='w', edgecolor='w')
		plt.close()
    	#plt.show()

	def reactivity_c(self):
		self.reactivity = []
		for i in range(0, len(self.stops_modification)):
			try:
				self.reactivity.append(max(self.stops_modification[i] - self.stops_control_norm[i], 0))
			except KeyError:
				self.reactivity.append(0)

	def norm2_8(self):
		FC = []
		FC_filtr = []
		for i in range(0, len(self.stops_modification)):
			try:
				if self.stops_control_norm[i] != 0 and self.stops_modification[i] != 0: 
					FC.append(math.log(float(self.stops_modification[i]) / self.stops_control_norm[i], 2))
				else:
					FC.append(0)
			except KeyError:
				FC.append(0)
		#print FC
		self.FC_lim = np.std(FC)

		sort_react = sorted(self.reactivity)
		n = int(math.floor(len(sort_react) * 0.9))
		m = int(math.floor(len(sort_react) * 0.98))
		c = sum(sort_react[n:m]) / (m - n)
		up = sort_react[m]
		
		file_out = open(arg.o, "w")

		#print len(FC), len(self.reactivity)

		for key in range(0, len(self.reactivity)):
			if FC[key] >= self.FC_lim or (FC[key] == 0 and self.reactivity[key] > 0):
				FC_filtr.append("P")
				if self.reactivity[key] > up:
					self.reactivity[key] = up / c
				else:
					self.reactivity[key] = self.reactivity[key] / c
			else:
				FC_filtr.append("F")
				if self.reactivity[key] > up:
					self.reactivity[key] = up / c
				else:
					self.reactivity[key] = self.reactivity[key] / c

			#print self.idt , key + 1, self.stops_modification[key], self.stops_control[key], self.stops_control_norm[key]  ,self.reactivity[key] , FC[key], FC_filtr[key] 
			
			out = "\t".join([self.idt , str(key + 1) , str(self.stops_modification[key]) ,  str(self.stops_control[key]) , str(self.stops_control_norm[key]) , str(self.reactivity[key]) , str(FC[key]) , str(FC_filtr[key])])
			file_out.write(out + "\n")

	def out_empty(self):
		file_out = open(arg.o, "w")
		#print arg.o
		for key in range(0, len(self.stops_modification)):
			out = "\t".join([self.idt , str(key + 1) , str(self.stops_modification[key]) ,  str(self.stops_control[key]) , str(0) , str(0) , str(0) , str(0)])
			file_out.write(out + "\n")
		f_o = open("maxy.txt", 'w')
		f_o.write(str(0))
class Input:
	def get_stops(self, bam,idt): # counting stops for each position in transcript
		reads={}
		samfile = pysam.AlignmentFile(bam, "rb")
		iter = samfile.fetch(idt)
		for read in iter:
			try:
				pos = read.reference_start
				if pos != 0:
					if read.is_reverse:
						pass
					#if read.is_read1 and read.is_proper_pair:
					else:
						reads[read.reference_start] += 1
					#else:
					#	print read.is_read2
			except KeyError:
				reads[read.reference_start] = 1
		for i in range(0, max(reads.keys())):
			try:
				reads[i]
			except KeyError:
				reads[i] = 0
		return reads

	def input_f(self):
		id_tp = ""
		start = False
		id_p = ""
		if arg.i:
			with open(arg.i) as inp:
				for i in inp:
					#print i
					line = i.split()
					#print id_p, line[0]
					if id_p == line[0]:
						t.stops_modification.append(int(line[3]))
						t.stops_control.append(int(line[2]))
						id_p = line[0]
					elif start == True:
						t.idt = id_p
						t.length = len(t.stops_control)
						t.norm_mean()
						try:
							t.norm_kernel_densities()
							t.reactivity_c()
							t.norm2_8()
						except ValueError:
							t.out_empty()

						start == False
						id_p = line[0]		
					if start == False:
						t = Transcript()
						id_p = line[0]
						t.stops_modification.append(int(line[3]))
						t.stops_control.append(int(line[2]))
						start = True
				try:
					t.idt = id_p
					t.length = len(t.stops_control)
					t.norm_mean()
				except UnboundLocalError:
					pass
				try:
					t.norm_kernel_densities()
					t.reactivity_c()
					t.norm2_8()
				except ValueError:
					t.out_empty()
		elif arg.c and arg.t:
			#get list of transcripts from bam file
			samfile = pysam.AlignmentFile(arg.c, "rb")
			ids = {}
			for i in samfile.header["SQ"]:
				if not arg.id:
					ids[i['SN']] = i['LN']
				else:
					if i['SN'] in arg.id.split(","):
						ids[i['SN']] = i['LN']
			#for each transcript
			for idt in ids.keys():
				t = Transcript()
				t.idt = idt
				t.length = ids[idt]
				t.stops_modification = self.get_stops(arg.t, t.idt)
				t.stops_control = self.get_stops(arg.c, t.idt)
				if len(t.stops_modification) > 0.4 * ids[idt] and len(t.stops_control) > 0.4 * ids[idt]:
					t.idt = id_p
					t.norm_mean()
					if sum(t.FC):
						t.norm_kernel_densities()
						t.reactivity_c()
						t.norm2_8()

arg = Arguments()
arg.parse()
a = Input().input_f()
