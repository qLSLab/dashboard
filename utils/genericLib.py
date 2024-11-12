import os
import sys
import time
import re
import pandas as pd
import base64
import datetime
import io
import xml.etree.ElementTree as ET
from dash import html

def setWorkingDirs(wrkDir=None, dizDirNames=None):
    """Set the working directories as recived from input dictionary.

    If arguments are omitted, current working directory is set as
    workingDirectory and no directory is nested.
    If wrkDir does not exists, the script create it.
    wrkDir/dir1
          /dir2
          /dir3

    Keyword arguments:
     wrkDir: The working directory path. if None it is set to cwd
     dizDirNames: The dictionary containing the names of the nested
                  directories, in the form: {'dirName', 'dirPath'}

    Return:
    (wrkDirPath, dizDirPaths) where:
    wrkDirPath: The working dir path
    dizDirPaths: The dictionary containing the paths of the nested directories,
     in the form: {'dirName', 'dirPath'}
    """
    cwDir = os.getcwd()
    dizDirPaths = {}
    if wrkDir is None:
        wrkDirPath = cwDir + os.sep
    else:
        if not os.path.exists(wrkDir):
            print("the working directory", wrkDir, "does not exist. Fix it.")
            sys.exit()
        else:
            wrkDirPath = os.path.abspath(wrkDir)
            if wrkDirPath[-1] != os.sep:
                wrkDirPath += os.sep
    if dizDirNames is not None:
        for el in dizDirNames:
            dizDirPaths[el] = pathJoinCheck(dizDirNames[el], wrkDirPath)
    return (wrkDirPath, dizDirPaths)

def pathJoinCheck(dir2add, rootPath="."):
    """Check the existence of a path and build it if necessary."""
    path = os.path.join(rootPath, dir2add)
    if not os.path.exists(path):
        os.makedirs(path)
    return path

# setting working dirs
CWDIR, dDirs = setWorkingDirs(
    dizDirNames={
        "raw": "rawData",
        "out": "outputs",
        "mods": "models",
        "logs": "logs",
        "figs": "figures",
        "utils": "utils",
        "map": "maps",
        "scripts": "scripts"
    },
)

def logFileOpen(logDIR=None, timeStamp=None, aim=None):
    cwDir = os.getcwd()
    if logDIR is None:
        logDIR = cwDir + os.sep
    if timeStamp is None:
        timeStamp = getTimeStamp()
    logFileName = os.path.join(logDIR, timeStamp + "_" + aim + ".log")
    logStream = open(logFileName, mode="w")
    return logStream

def getBaseName(path):
    """The program filename without extension"""
    return os.path.basename(path).split(".")[0]

def getTimeStamp():
    return time.strftime("%Y%m%d%H%M%S", time.localtime())

def toLog(logStream, string):
    logStream.write(string + "\n")

def writeLineByLineToFile(stream, dataToWrite, spaziatore):
    stream.write(spaziatore.join(str(s) for s in dataToWrite) + "\n")

def beingOrExit(prgPath):
    if not os.path.exists(prgPath):
        print("The file ", prgPath, "does not exist, check the path")
        sys.exit()

def dizReaProd(s):
    termini = s.strip().split(" + ")
    diz = {}
    for termine in termini:

        coeffMetabo = termine.split()
        coeff = coeffMetabo[0]
        if isCoeff(coeff) is True:
            coeff = float(coeffMetabo[0])
            metabolita = " ".join(coeffMetabo[1:])
        else:
            metabolita = " ".join(coeffMetabo)
            coeff = 1.0
        diz[metabolita] = coeff
    return diz


def isCoeff(s):
    """Determine if a string splitted on the spaces the first element is the
    stoichiometric coefficient or not.
    Example: if string is "2 A" return True; if string is "A" it returns False"""
    answer = re.match("((\d+(\.\d*)?|\.\d+)([eE][-+]?\d+)?)$", s)
    if answer is not None:
        return True
    else:
        return False

def getPrimeNumberInRange(start, stop):
    l = []
    for num in range(start,stop):
        if all(num%i!=0 for i in range(2,num)):
            l.append(num)
    return l

def unique(a):
    """return the list with duplicate elements removed"""
    return list(set(a))

def intersect(a, b):
    """return the intersection of two lists"""
    return list(set(a) & set(b))


def union(a, b):
    """return the union of two lists"""
    return list(set(a) | set(b))


def difference(a, b):
    """return the difference of two lists"""
    return list(set(a) - set(b))

def extractRegexFromItem(item, regex):
    sItem = pd.Series(data=item)
    if sItem.empty == False:
        dfItem = sItem.str.extractall(regex)
    else:
        dfItem = []
    return dfItem

def parse_contents(contents, filename):
    MODELDIR = dDirs["mods"]
    RAWDIR = dDirs["raw"]
    MAPDIR = dDirs["map"]


    try:
        content_type, content_string = contents.split(',')
        decoded = base64.b64decode(content_string)
        fileWExt = os.path.splitext(filename)
        try:
            if "xml" in fileWExt[1]:
                tree = ET.ElementTree(ET.fromstring(decoded))
                tree.write(os.path.join(MODELDIR, filename))
                return filename

            elif "svg" in fileWExt[1]:
                tree = ET.ElementTree(ET.fromstring(decoded))
                tree.write(os.path.join(MAPDIR, filename))
                return filename

            elif "tsv" in fileWExt[1]:
                # print("content_type\n", content_type, "\n\ncontent_string\n", content_string)
                df = pd.read_csv(
                    io.StringIO(decoded.decode('utf-8')), sep = "\t")
                df.to_csv(os.path.join(RAWDIR, filename), sep = "\t", index = False)
                return filename
        except Exception as e:
            print(e)
            return html.Div([
                'There was an error processing this file.'
            ])

    except:
        filename = filename[0]
        content_type, content_string = contents[0].split(',')
        decoded = base64.b64decode(content_string)
        fileWExt = os.path.splitext(filename)
        try:
            if "xml" in fileWExt[1]:
                tree = ET.ElementTree(ET.fromstring(decoded))
                tree.write(os.path.join(MODELDIR, filename))
                return filename

            elif "svg" in fileWExt[1]:
                tree = ET.ElementTree(ET.fromstring(decoded))
                tree.write(os.path.join(MAPDIR, filename))
                return filename

            elif "tsv" in fileWExt[1]:
                # print("content_type\n", content_type, "\n\ncontent_string\n", content_string)
                df = pd.read_csv(
                    io.StringIO(decoded.decode('utf-8')), sep = "\t")
                df.to_csv(os.path.join(RAWDIR, filename), sep = "\t", index = False)
                return filename
        except Exception as e:
            print(e)
            return html.Div([
                'There was an error processing this file.'
            ])



def getRowdataList(lParams):
    dParams = {}
    for p in lParams:
        dParams[p] = None
    return [dParams] * 5
