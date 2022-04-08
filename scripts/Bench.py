from datetime import datetime
from itertools import *

import os, sys, getopt, errno
import re

import shutil
import subprocess

import json

from sklearn.model_selection import ParameterGrid


class Experiment:
    Restarts = 1
    def __init__(self, exeFilePath, inputFilePath, inputFiles, outputFilePath, parameterFilePath, paramterFileName, environmentPath):
        self.exeFilePath = exeFilePath
        self.inputFilePath = inputFilePath
        self.inputFiles = inputFiles
        self.outputFilePath = outputFilePath
        self.parameterFilePath = parameterFilePath
        self.paramterFileName = paramterFileName
        self.environmentPath = environmentPath

    def Run(self):        
        # creates the directory in which to store the log files in
        timestampForFiles = datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
        logFolderName = os.path.join(self.outputFilePath, timestampForFiles)
        
        CreateDirIfNecessary(logFolderName)

        problemInstanceDir = os.path.join(logFolderName, 'ProblemInstances', '')
        CreateDirIfNecessary(problemInstanceDir)

        for inputFileName in self.inputFiles:
            shutil.copy2(os.path.join(self.inputFilePath, inputFileName), problemInstanceDir) # copy files into benchmark dir
        
        parameterPermutations = self.CreateParameterPermutations(logFolderName)

        # create subfolder for each permutation
        for shorthand, parameterPermutation in parameterPermutations.items():
            permutationDir = os.path.join(logFolderName, shorthand, '')

            # TODO: scan subfolder and check membership of every file with self.inputFiles (which can contain a regex expression)
            for inputFileName in self.inputFiles:
                instanceDir = os.path.join(permutationDir, inputFileName[:-5])
                CreateDirIfNecessary(instanceDir)

                self.ExecuteRepeatedly(instanceDir, inputFileName, permutationDir, self.Restarts)
        

    def CreateParameterPermutations(self, benchRootDir):
        parameterPermutations = {}
        
        paramPool = self.CreateParameterPool() # returns ParameterGrid from sklearn

        # Full factorial analysis
        count = 0
        for permutation in paramPool:
            shorthand = 'P-' + str(count)
            
            parameterPermutations[shorthand] = permutation

            permutationDir = os.path.join(benchRootDir, shorthand, '')
            CreateDirIfNecessary(permutationDir)

            with open(os.path.join(permutationDir, 'parameters.json'), 'w') as paramFile:
                json.dump(permutation, paramFile, indent = 4)
                
            count += 1

        return parameterPermutations

    def CreateParameterPool(self):  
        if self.paramterFileName != '':
            parameterPermutations = []
            parameters = {}
            with open(os.path.join(self.parameterFilePath, self.paramterFileName)) as f:
                parameters = json.load(f)

            parameterPermutations.append(parameters)

            return parameterPermutations
            
        # TODO: Modify parameter file
        # raise NotImplementedError

        parameterRanges = {}

        # sklearn.ParameterGrid creates full factorial parameter combinations
        return ParameterGrid(parameterRanges)

    def ExecuteRepeatedly(self, logFolderName, inputFileName, permutationDir, restarts):
        currentRound = 0
        for i in range(restarts):
            # TODO: instead of creating a subfolder for each run, consider writing the files to a tmp directory, then rename and move into logFolderName
            runDir = os.path.join(logFolderName, 'Run-' + str(i), '')
            CreateDirIfNecessary(runDir)
            
            cmd = [self.exeFilePath.replace('\\', '/'), 
            "-i", self.inputFilePath.replace('\\', '/'), 
            "-f", inputFileName, 
            "-o", runDir.replace('\\', '/'), 
            "-p", os.path.join(permutationDir, 'parameters.json').replace('\\', '/')]

            print(cmd)
            self.StartProcess(cmd)
            
            currentRound += 1
            print('{} - {} with {}: {} / {} restarts done'.format(
                datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
                inputFileName[:-5],
                permutationDir[-4:-1],
                currentRound, 
                restarts))

            
    def StartProcess(self, cmd):
        osenv = os.environ.copy()
        osenv["PATH"] = self.environmentPath + osenv["PATH"]
        
        p = subprocess.check_call(cmd, env=osenv) #, shell=True, stdout=open(os.devnull, "w")
        
        """
        p = subprocess.Popen(cmd, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.STDOUT) # close_fds=True
        
        p.stdin.close()
        output, err = p.communicate()
            
        p.terminate()
        
        p.wait()
        """

def CreateDirIfNecessary(outputDir):
    if not os.path.exists(os.path.dirname(outputDir)):
        try:
            os.makedirs(os.path.dirname(outputDir))
        except OSError as exc: # Guard against race condition
            if exc.errno != errno.EEXIST:
                raise


def main():
    workingDirectory = os.getcwd()

    exeFilePath = os.path.join(workingDirectory, '../../../x64/Release/', 'ContainerLoadingApplication.exe')
    inputFilePath = os.path.join(workingDirectory, 'Converted/CJCM_Unibo/', '')
    outputFilePath = os.path.join(workingDirectory, 'Benchmarks', '') 
    parameterFilePath = os.path.join(workingDirectory, 'Parameters') 
    paramterFileName = 'Parameters-TSBP.json'
    environmentPath = '$(GUROBI_HOME)/include;C:/Program Files/or-tools/include;$(GUROBI_HOME)/lib;C:/Program Files/or-tools/lib;C:/Program Files/SCIPOptSuite 7.0.2/lib;C:/Program Files/SCIPOptSuite 7.0.1/include;'

    # specify input files for benchmark and copy into benchmark subfolder.
    inputFiles = [
    'E00N10.json',
    #'E00N15.json',
    #'E00N23.json',
    #'E00X23.json',
    #'E02F17.json',
    #'E02F20.json',
    #'E02F22.json',
    #'E02N20.json',
    #'E03N10.json',
    #'E03N15.json',
    #'E03N16.json',
    #'E03N17.json',
    #'E03X18.json',
    #'E04F15.json',
    #'E04F17.json',
    #'E04F19.json',
    #'E04F20.json',
    #'E04N15.json',
    #'E04N17.json',
    #'E04N18.json',
    #'E05F15.json',
    #'E05F18.json',
    #'E05F20.json',
    #'E05N15.json',
    #'E05N17.json',
    #'E05X15.json',
    #'E07F15.json',
    #'E07N10.json',
    #'E07N15.json',
    #'E07X15.json',
    #'E08F15.json',
    #'E08N15.json',
    #'E10N10.json',
    #'E10N15.json',
    #'E10X15.json',
    #'E13N10.json',
    #'E13N15.json',
    #'E13X15.json',
    #'E15N10.json',
    #'E15N15.json',
    #'E20F15.json',
    #'E20X15.json'
    ]

    experiment = Experiment(exeFilePath, inputFilePath, inputFiles, outputFilePath, parameterFilePath, paramterFileName, environmentPath)
    experiment.Run()
    

if __name__ == "__main__":
   main()