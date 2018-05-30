
import sys
import math
#import matplotlib.pyplot as plt
#from scipy import stats
#import numpy as np
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
		parser.add_argument("-o", help = "output file", required = True, type = str)
		args = parser.parse_args()
		self.c = args.c
		self.t = args.t
		self.o = args.o

class Transcript:
	def __init__(self):
		self.stops_control = []
		self.stops_modification = []
		self.length = 0
		self.id = ""
	def output(self):

		with open(arg.o, 'a') as out:
			for i in (0, self.length):
				out.write(self.id + "\t" + str(i + 1)+ "\t" + self.stops_control[i] + "\t" + self.stops_modification[i] + "\n")


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
		id_p = ""
		outp = open(arg.o, 'w').close() #clear file
		if arg.c and arg.t:
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
				t.output()


arg = Arguments()
arg.parse()
a = Input().input_f()
