import os
import sys
import subprocess
import matplotlib.pyplot as plt
import matplotlib
import pickle
import locale
import re
import numpy as np
import math
import itertools
from multiprocessing import Pool


curPath = "."
parentPath = ".."
defaultBinary = parentPath + "/bin/Shot"
synPrefix = curPath + "/synData/"
modelPrefix = curPath + "/model/"
outPrefix = curPath + "/out/"
figPrefix = curPath + "/fig/"

optLossModel = "--use-modelloss"
optLossCore = "--use-coreloss"
optLossIter = "--use-iteration"
optIteration = "--iteration"
optEpsilon = "--epsilon"
optRow = "--plain"
optTrans = "--scan"
optSzMode = "--size-mode"
optSzThr = "--size-thread"
optSzRank = "--size-rank"

defaultParam = {
	'Loss': optLossIter,
	'Thread': 4,
	'Iteration': 1,
	'Epsilon': 0.0,
	'Mode': [4],
	'NNZ': [100000],
	'NNZDim': [],
	'Dim': [10000],
	'DimRank': [20000],
	'Rank': [8]
}

rangeParam = {
	'Mode': [3, 4, 5, 6], 
	'NNZ': [10000, 50000, 100000, 500000, 1000000, 500000, 1000000, 5000000, 10000000],
	'Dim': [1000, 2000, 5000, 10000, 100000, 1000000, 10000000],
	'Rank': [3, 4, 5, 6, 7, 8, 10, 12, 14, 16],
	'Method': ["TuckerALS", "MET1", "plain", "scan"], 
	'OTrial': [1, 2, 3],
	'ITrial': [1, 2, 3]
}

defaultLoss = defaultParam['Loss']
defaultThread = defaultParam['Thread']
defaultIteration = defaultParam['Iteration']
defaultEpsilon = defaultParam['Epsilon']
defaultMode = defaultParam['Mode'][0]
defaultRank = defaultParam['Rank'][0]

# Customizing matplotlib
figSize = (6, 4.5) 
matplotlib.rcParams.update({
		'font.size': 18,
		'font.family': 'serif',
		'lines.linewidth': 3, 
		'lines.markersize': 7, 
		'figure.dpi': 300,
		'nbagg.scanparent': True,
		'legend.fontsize': 'small',
		'legend.scatterpoints': 1,
		'legend.fancybox': True
		})

labelSet = {
	"plain": "S-HOT",
	"scan": "S-HOT$_{scan}$",
	"TuckerALS": "NaiveTucker",
	"MET1": "TuckerMatlab",
	"Haten2": "Haten2"
}

colorSet = {
	"plain": 'b',
	"scan": 'r',
	"TuckerALS": "g",
	"MET1": "k"
}

lineSet = {
	"plain": '--',
	"scan": '-',
	"TuckerALS": '-.', 
	"MET1": ':'
}

faceSet = {
	"plain": 's',
	"scan": 'o',
	"TuckerALS": '^',
	"MET1": '*',
	"Haten2": 'x'
}

def unitMatlabCall(method, filename, outfilename, mode, rank, **kwargs):

	szIteration = kwargs["Iteration"] if kwargs.has_key("Iteration") else defaultIteration
	szEpsilon = kwargs["Epsilon"] if kwargs.has_key("Epsilon") else defaultEpsilon

	cmd = ["matlab", "-nodisplay", "-r", "\"cd %s/baselines; try run%s(\'.%s\', %d, %d, %f, %d); catch display(lasterr); end; quit force\""%(curPath, method, filename, mode, rank, szEpsilon, szIteration), "|", "tee", outfilename]
	cmd = " ".join(cmd)
	print "Processing:\t" + cmd

	outMessage = [];

	proc = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE);
	while True:
		out = proc.stdout.read(1)
		if out == '' and proc.poll() != None:
			break
		if out != '':
			sys.stdout.write(out)
			sys.stdout.flush()
			outMessage.append(out)
	
	print "Done:\t" + cmd
	text = ''.join(outMessage)
	outMessage = text.split('\n')
	outMessage = outMessage[:-1]
	return outMessage


def unitLearningCall(dataFileName, modelFileName, outputFileName, **kwargs):

	method = "--" + kwargs["Method"]
	szMode = kwargs["Mode"] if kwargs.has_key("Mode") else defaultMode
	szThread = str(kwargs["Thread"] if kwargs.has_key("Thread") else defaultThread)
	szRankBase = str(kwargs["Rank"] if kwargs.has_key("Rank") else defaultRank)
	szIteration = str(kwargs["Iteration"] if kwargs.has_key("Iteration") else defaultIteration)
	szEpsilon = str(kwargs["Epsilon"] if kwargs.has_key("Epsilon") else defaultEpsilon)
	termCond = kwargs["TermCond"] if kwargs.has_key("TermCond") else defaultLoss
	szRank = szRankBase
	for i in range(szMode-1):
		szRank = szRank + ":" + szRankBase
	szMode = str(szMode);

	outFile = ["2>&1", "|", "tee", outputFileName]

	cmd = [defaultBinary, method, termCond, optSzMode, szMode, optSzRank, szRank, optSzThr, szThread,
				optIteration, szIteration, optEpsilon, szEpsilon, dataFileName, modelFileName] + outFile

	cmd = " ".join(cmd)

	print "Processing:\t" + cmd

	outMessage = [];

	proc = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE);
	while True:
		out = proc.stdout.read(1)
		if out == '' and proc.poll() != None:
			break
		if out != '':
			sys.stdout.write(out)
			sys.stdout.flush()
			outMessage.append(out)
	
	print "Done:\t" + cmd
	text = ''.join(outMessage)
	outMessage = text.split('\n')
	outMessage = outMessage[:-1]
	return outMessage

def parseMatlabOutput(outMsg):
	stepList = []
	lineAcc = []
	for line in outMsg[:-1]:
		if not line.startswith("Iter "):
			continue

		if line.split()[2] == "mode":
			lineAcc.append(-1)
		else:
			num = line.split()[1].split(":")[0]
			lineAcc = [int(num)] + map(lambda x: float(x), line.split()[2:5]) + lineAcc
			stepList.append(lineAcc)
			lineAcc = []
	if len(stepList) > 0:
		line = stepList[-1] 
		time = float(line[1])
		err = float(line[3]);
		loss = float(line[2])
		print time, 0, err, loss, stepList
		print time, 0, err, [stepList[0], stepList[-1]]
		return time, 0, err, loss, stepList
	else:
		return False;

def parseOutput(outMsg):
	stepList = []
	lineAcc = []
	loss = []
	for line in outMsg[:-2]:
		if line.startswith("Iter#"):
		    num = line.split()[0].split("#")[1]
		    lineAcc = [int(num)] + map(lambda x: float(x), line.split()[1:]) + lineAcc
		    stepList.append(lineAcc)
		    lineAcc = []
		else:
		    pass
	line = outMsg[-2].split()
	time = float(line[-1])
        err = 0
	line = outMsg[-1].split()
	evalTime = float(line[-1])
	return time, evalTime, err, [stepList[0], stepList[-1]]

def synTensor(params):
	(mode, dim, nnz, ot) = params
	fileID = "syn%d_d%.0E_n%.0E_ot%d" % (mode, dim, nnz, ot)
	fullsynpath = synPrefix + fileID
	if (not os.path.isfile(fullsynpath)):
		print "Synthesizing:\t" + fullsynpath
		filePtr = open(fullsynpath, "w")
		indexList = np.random.choice(dim, (nnz, mode+1), replace=True)
		for line in indexList:
			filePtr.write(' '.join(map(lambda x: str(x), line[:-1]))+" "+str.format("{0:<.6f}", float(line[-1])/dim) + "\n")
	
		filePtr.close();
		print "Done:\t" + fullsynpath
	
	return fileID

def plot_NNZ_WCK_Conv(data, idx):

	for mode, dim, rank in itertools.product(defaultParam['Mode'], defaultParam['Dim'], defaultParam['Rank']):
		fileName = figPrefix + ("NNZ_WCK_m%d_d%.0E_r%d_i%d" % (mode, dim, rank, idx + 1)) + ".pdf"
		if (os.path.isfile(fileName)):
			continue

		resLine = {} 
                errLine = {}
		for method, nnz in itertools.product(rangeParam['Method'], rangeParam['NNZ']):
			line = resLine.get(method, [])
                        err = errLine.get(method, [])
			try:
				line.append(data[mode, dim, nnz, method, rank][idx])
				err.append(data[mode, dim, nnz, method, rank][idx+2])
				resLine[method] = line
                                errLine[method] = err
			except:
				pass

		fig = plt.figure(figsize=figSize)
		print fileName
		for key in rangeParam['Method']:
			try:
				print key, resLine[key]
				plt.plot(rangeParam['NNZ'], resLine[key],
						label=labelSet[key],
						linestyle=lineSet[key],
						color=colorSet[key],
						marker=faceSet[key]
						)
			except:
				pass

		plt.xlabel('The number of Non-Zero')
		if (idx == 0):
			plt.ylabel('Elapsed time for\n'+str(idx+1)+' iteration (sec)')
		else:
			plt.ylabel('Elapsed time\nuntil convergence (sec)')
		plt.tight_layout(rect=[0, 0, 1, 0.95])
		plt.xscale('log');
		plt.yscale('log');
		plt.draw()
		plt.savefig(fileName, scanparent=True)
		plt.close(fig)

def plot_Dim_WCK_Conv(data, idx):

	for mode, nnz, rank in itertools.product(defaultParam['Mode'], defaultParam['NNZ'] + defaultParam['NNZDim'], defaultParam['Rank']):
		fileName = figPrefix + ("Dim_WCK_m%d_n%.0E_r%d_i%d" % (mode, nnz, rank, idx + 1)) + ".pdf"
		if (os.path.isfile(fileName)):
			continue

		resLine = {} 
                errLine = {}
		for method, dim in itertools.product(rangeParam['Method'], rangeParam['Dim']):
			line = resLine.get(method, [])
                        err = errLine.get(method, [])
			try:
				line.append(data[mode, dim, nnz, method, rank][idx])
				err.append(data[mode, dim, nnz, method, rank][idx+2])
				resLine[method] = line
                                errLine[method] = err
			except:
				pass

		fig = plt.figure(figsize=figSize)
		print fileName
		for key in rangeParam['Method']:
			try:
				print key, resLine[key]
				plt.plot(rangeParam['Dim'][:len(resLine[key])], resLine[key],
						label=labelSet[key],
						linestyle=lineSet[key],
						color=colorSet[key],
						marker=faceSet[key]
						)
			except:
				pass

		plt.xlabel('Tensor dimensionality')
		if (idx == 0):
			plt.ylabel('Elapsed time for\n'+str(idx+1)+' iteration (sec)')
		else:
			plt.ylabel('Elapsed time\nuntil convergence (sec)')
		plt.tight_layout(rect=[0, 0, 1, 0.95])
		plt.xscale('log');
		plt.yscale('log');
		plt.draw()
		plt.savefig(fileName, scanparent=True)
		plt.close(fig)

def plot_Mode_WCK_Conv(data, idx):

	for nnz, dim, rank in itertools.product(defaultParam['NNZ'], defaultParam['Dim'], defaultParam['Rank']):
		fileName = figPrefix + ("Mode_WCK_d%.0E_n%.0E_r%d_i%d" % (dim, nnz, rank, idx + 1)) + ".pdf"
		if (os.path.isfile(fileName)):
			continue

		resLine = {}
                errLine = {}
		for method, mode in itertools.product(rangeParam['Method'], rangeParam['Mode']):
			line = resLine.get(method, [])
                        err = errLine.get(method, [])
			try:
				line.append(data[mode, dim, nnz, method, rank][idx])
				err.append(data[mode, dim, nnz, method, rank][idx+2])
				resLine[method] = line
                                errLine[method] = err
			except:
				pass

		fig = plt.figure(figsize=figSize)
		print fileName
		for key in rangeParam['Method']:
			try:
				print key, resLine[key]
				plt.plot(rangeParam['Mode'][:len(resLine[key])], resLine[key], 
						label=labelSet[key],
						linestyle=lineSet[key],
						color=colorSet[key],
						marker=faceSet[key]
						)
			except:
				pass

		plt.xlabel('Tensor order')
		if (idx == 0):
			plt.ylabel('Elapsed time for\n'+str(idx+1)+' iteration (sec)')
		else:
			plt.ylabel('Elapsed time\nuntil convergence (sec)')
		plt.xticks(range(3, rangeParam['Mode'][-1]+1))
		plt.xlim(np.array([2, rangeParam['Mode'][-1]])+0.5)
                plt.ylim([1, 1E+3])
		plt.yscale('log');
		plt.tight_layout(rect=[0, 0, 1, 0.95])
		plt.draw()
		plt.savefig(fileName, scanparent=True)
		plt.close(fig)

def plot_Rank_WCK_Conv(data, idx):

	for mode, nnz, dim in itertools.product(defaultParam['Mode'], defaultParam['NNZ'], defaultParam['Dim'] + defaultParam['DimRank']):
		fileName = figPrefix + ("Rank_WCK_m%d_d%.0E_n%.0E_i%d" % (mode, dim, nnz, idx + 1)) + ".pdf"
		if (os.path.isfile(fileName)):
			continue



		resLine = {}
                errLine = {}
		for method, rank in itertools.product(rangeParam['Method'], rangeParam['Rank']):
			line = resLine.get(method, [])
                        err = errLine.get(method, [])
			try:
				line.append(data[mode, dim, nnz, method, rank][idx])
				err.append(data[mode, dim, nnz, method, rank][idx+2])
				resLine[method] = line
                                errLine[method] = err
			except:
				pass
		
		fig = plt.figure(figsize=figSize)
		print fileName
		for key in rangeParam['Method']:
			try:
				print key, resLine[key]
				plt.plot(rangeParam['Rank'][:len(resLine[key])], resLine[key],
						label=labelSet[key],
						linestyle=lineSet[key],
						color=colorSet[key],
						marker=faceSet[key]
						)
			except:
				pass

		plt.xlabel('Rank')
		if (idx == 0):
			plt.ylabel('Elapsed time for\n'+str(idx+1)+' iteration (sec)')
		else:
			plt.ylabel('Elapsed time\nuntil convergence (sec)')
		plt.tight_layout(rect=[0, 0, 1, 0.95])
		plt.yscale('log');
		plt.draw()
		plt.savefig(fileName, scanparent=True)
		plt.close(fig)

if __name__ == "__main__":

	resHash = {}

	expSet = []
	expSet = expSet + list(itertools.product(rangeParam['Mode'], defaultParam['NNZ'], defaultParam['Dim'], rangeParam['Method'], defaultParam['Rank'], rangeParam['OTrial'], rangeParam['ITrial']))
	expSet = expSet + list(itertools.product(defaultParam['Mode'], rangeParam['NNZ'], defaultParam['Dim'], rangeParam['Method'], defaultParam['Rank'], rangeParam['OTrial'], rangeParam['ITrial']))
	expSet = expSet + list(itertools.product(defaultParam['Mode'], defaultParam['NNZ'], rangeParam['Dim'], rangeParam['Method'], defaultParam['Rank'], rangeParam['OTrial'], rangeParam['ITrial']))
	expSet = expSet + list(itertools.product(defaultParam['Mode'], defaultParam['NNZDim'], rangeParam['Dim'], rangeParam['Method'], defaultParam['Rank'], rangeParam['OTrial'], rangeParam['ITrial']))
	expSet = expSet + list(itertools.product(defaultParam['Mode'], defaultParam['NNZ'], defaultParam['Dim'] + defaultParam['DimRank'], rangeParam['Method'], rangeParam['Rank'], rangeParam['OTrial'], rangeParam['ITrial']))

	for mode, nnz, dim, method, rank, ot, it in expSet:
		fileID = "syn%d_d%.0E_n%.0E_ot%d" % (mode, dim, nnz, ot)
		datafilename = synPrefix + fileID
		if (not os.path.isfile(datafilename)):
			fileID = synTensor((mode, dim, nnz, ot))
		
		postfix = ("_%s_r%d_t%d"%(method, rank, it))
		modelfilename = modelPrefix + fileID + postfix
		outfilename = outPrefix + fileID + postfix
		if (not os.path.isfile(outfilename)):
			if method in ["plain", "scan"]:
				outMsg = unitLearningCall(datafilename, modelfilename, outfilename, Method=method, Mode=mode, Rank=rank)
				os.remove(modelfilename)
			else:
				unitMatlabCall(method, datafilename, outfilename, mode, rank)

		if method in ["plain", "scan"]:
			outMsg = open(outfilename, "r").readlines()
			retVal = parseOutput(outMsg)
		else:
			outMsg = open(outfilename, "r").readlines()
                        print method, 
			retVal = parseMatlabOutput(outMsg)
		if not retVal == False:
			resHash[mode, dim, nnz, method, rank, ot, it] = retVal

	expSet = []
	expSet = expSet + list(itertools.product(defaultParam['Mode'], defaultParam['NNZDim'], rangeParam['Dim'], rangeParam['Method'], defaultParam['Rank']))
	expSet = expSet + list(itertools.product(defaultParam['Mode'], defaultParam['NNZ'], rangeParam['Dim'], rangeParam['Method'], defaultParam['Rank']))
	expSet = expSet + list(itertools.product(defaultParam['Mode'], rangeParam['NNZ'], defaultParam['Dim'], rangeParam['Method'], defaultParam['Rank']))
	expSet = expSet + list(itertools.product(defaultParam['Mode'], defaultParam['NNZ'], defaultParam['Dim'], rangeParam['Method'], rangeParam['Rank']))
	expSet = expSet + list(itertools.product(rangeParam['Mode'], defaultParam['NNZ'], defaultParam['Dim'], rangeParam['Method'], defaultParam['Rank']))
	expSet = expSet + list(itertools.product(defaultParam['Mode'], defaultParam['NNZ'], defaultParam['Dim'] + defaultParam['DimRank'], rangeParam['Method'], rangeParam['Rank']))

	# Data post processing: report average value
	for mode, nnz, dim, method, rank in expSet:
		otSumVal = [0, 0, 0, 0]
		for ot in rangeParam['OTrial']:
			itSumVal = [0, 0, 0, 0]
			for it in rangeParam['ITrial']:
				try:
					retVal = resHash[mode, dim, nnz, method, rank, ot, it]
					# loss print
					# print mode, dim, nnz, method, rank, retVal[-2]
					# 1. only collects w.c.k for 1iter, until convergence
					itSumVal[0] = itSumVal[0] + retVal[-1][0][1]
					itSumVal[1] = itSumVal[1] + retVal[-1][-1][1]
					itSumVal[2] = itSumVal[2] + retVal[-1][0][1]**2
					itSumVal[3] = itSumVal[3] + retVal[-1][-1][1]**2
					#print mode, nnz, dim, method, rank, "IT", it, itSumVal
				except:
					pass

			if itSumVal != [0, 0, 0, 0]:
				itSumVal[0] = itSumVal[0]/len(rangeParam['ITrial'])
				itSumVal[1] = itSumVal[1]/len(rangeParam['ITrial'])
				itSumVal[2] = itSumVal[2]/len(rangeParam['ITrial'])
				itSumVal[3] = itSumVal[3]/len(rangeParam['ITrial'])
				resHash[mode, dim, nnz, method, rank, ot] = itSumVal
				otSumVal[0] = otSumVal[0] + itSumVal[0]
				otSumVal[1] = otSumVal[1] + itSumVal[1]
				otSumVal[2] = otSumVal[2] + itSumVal[2]
				otSumVal[3] = otSumVal[3] + itSumVal[3]

		if otSumVal != [0, 0, 0, 0]:
			#print mode, nnz, dim, method, rank, "OT", ot, otSumVal
			otSumVal[0] = otSumVal[0]/len(rangeParam['OTrial'])
			otSumVal[1] = otSumVal[1]/len(rangeParam['OTrial'])
			otSumVal[2] = otSumVal[2]/len(rangeParam['OTrial'])
			otSumVal[3] = otSumVal[3]/len(rangeParam['OTrial'])
			resHash[mode, dim, nnz, method, rank] = otSumVal

	for mode, nnz, dim, method, rank in expSet:
            try:
                otSumVal = resHash[mode, dim, nnz, method, rank]
                otSumVal[2] = math.sqrt(otSumVal[2] - otSumVal[0]**2)
                otSumVal[3] = math.sqrt(otSumVal[3] - otSumVal[1]**2)
                print mode, nnz, dim, method, rank, otSumVal
            except:
                pass
	
	plot_Rank_WCK_Conv(resHash, 0)
	plot_Mode_WCK_Conv(resHash, 0)
	plot_Dim_WCK_Conv(resHash, 0)
	plot_NNZ_WCK_Conv(resHash, 0)
	print "Completed!"
