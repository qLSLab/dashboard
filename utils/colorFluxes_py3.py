#/usr/bin/python2.7

#### colorFluxes.py

# given that the fileVaues lines have the following structure:
# idString value1 value2 ... valueN
# program calling:
# colorAndThickInkMap.py map.svg fileValues valueCol

import sys
from collections import defaultdict
from lxml import etree as ET
import matplotlib as mpl
import matplotlib.cm as cmx
import matplotlib.pyplot as plt
import math



### Parameters

maxThikc = 8.0
minThick = 1.0
epsilon = .00001
gray = (0.839216, 0.839216, 0.839216)

### Functions

### is_number --
### Function to manage matlab "dirty"-strings

def is_number(s):
    try:
        float(s)
        return True
    except ValueError:
        return False


def thicknessMapping(x, maxFlux):
    x = (x * maxThikc) / maxFlux
    x = max(x, minThick)
    return x


def findMapBound(valuesList):
    ## finding the min and max values excluding INFs
    ## print "-----\n" , valuesList

    while max(valuesList) == float('Inf'):
        valuesList.remove(float('Inf'))
    maxVal = max(valuesList)

    ## print "max Val = ", maxVal

    while min(valuesList) == -float('Inf'):
        valuesList.remove(-float('Inf'))
    minVal = min(valuesList)
    low = min(x for x in valuesList if x != 0)

    ## print "min Val = ", minVal
    ## print "low =", low
    ## fixing the boundaries for the color map to be simmetric around 0
    mapBoundary = max(abs(minVal), maxVal)
    ## print "mapBoundaryValue = ", mapBoundary
    return [mapBoundary, minVal, maxVal, low]


def fixStrokeColor(wordsListStr, toColorStr):
    wordsList = wordsListStr.split(';')
    match = [s for s in wordsList if "stroke:" in s]
    ## print wordsList
    iStroke = wordsList.index(match[0])
    # colorStr = match[0].translate(None, 'stroke:')
    # In python 3 it changes: Basically it says 'translate nothing to nothing' (first two parameters) and translate any punctuation or numbers to None (i.e. remove them).
    colorStr = match[0].translate(str.maketrans('','','stroke:'))
    ## print "          match = ", match, "colorStr= ", colorStr, "iStroke = ", iStroke
    wordsList[iStroke] = "stroke:" + toColorStr
    s = ";"
    ## print s.join(wordsList)
    return s.join(wordsList)


def fixStrokewidth(wordsListStr, toWidthStr):
    wordsList = wordsListStr.split(';')
    match = [s for s in wordsList if "stroke-width:" in s]
    ## print wordsList
    iStroke = wordsList.index(match[0])
    # widthStr = match[0].translate(None, 'stroke:')
    widthStr = match[0].translate(str.maketrans('','','stroke:'))
    ## print "          match = ", match, "colorStr= ", colorStr, "iStroke = ", iStroke
    wordsList[iStroke] = "stroke-width:" + toWidthStr
    s = ";"
    ## print s.join(wordsList)
    return s.join(wordsList)


def fixStrokeDasharray(wordsListStr, toDashStr):
    wordsList = wordsListStr.split(';')
    match = [s for s in wordsList if "stroke-dasharray:" in s]
    ## print wordsList
    if match != []:
        iStroke = wordsList.index(match[0])
        # print('match[0]\t', match[0])
        # print('iStroke\t', iStroke)
        # dashStr = match[0].translate(None, 'stroke:')
        dashStr = match[0].translate(str.maketrans('','','stroke:'))
        # print('dashStr\t', dashStr)
        ## print "          match = ", match, "colorStr= ", colorStr, "iStroke = ", iStroke
        wordsList[iStroke] = "stroke-dasharray:" + toDashStr
        ## print s.join(wordsList)
    s = ";"
    ## print s.join(wordsList)
    return s.join(wordsList)


# def sign_0(k):                            # Aaaaaargh!
#     if k == 0.0:
#         return 0.0
#     elif k > 0.0:
#         return 1.0
#     else:
#         return -1.0

def sign_0(k):
    if k == 0.0:
        return 0.0
    else:
        return math.copysign(1.0, k)



# INPUT

fileMapName = sys.argv[1]
fileValuesName = sys.argv[2]
idColValColumn = int(sys.argv[3])
idThickValColumn = int(sys.argv[4])
filename = sys.argv[5]

SODIPODI = 'http://sodipodi.sourceforge.net/DTD/sodipodi-0.dtd'
CC = 'http://creativecommons.org/ns#'
SVG = 'http://www.w3.org/2000/svg'
DC = 'http://purl.org/dc/elements/1.1/'
XLINK = 'http://www.w3.org/1999/xlink'
RDF = 'http://www.w3.org/1999/02/22-rdf-syntax-ns#'
INKSCAPE = 'http://www.inkscape.org/namespaces/inkscape'

NSMAP = {
    None: SVG,
    'sodipodi' :  SODIPODI,
    'cc' :  CC,
    'svg' :  SVG,
    'dc' :  DC,
    'xlink' :  XLINK,
    'rdf' :  RDF,
    'inkscape' :  INKSCAPE
    }

# print "parameters:\t", fileMapName

### Reading the values file,  starting point to map in color fluxes.

with open(fileValuesName, 'r') as fileValues:
    linesValues = fileValues.readlines()
fileValues.close

# print('linesValues\n', linesValues)


linesInFileVal = len(linesValues)

# print "#lines in dict", linesInFileVal

# for i in range(0, linesInFileVal):
#   print linesValues[i]


### Building the color dictionary

colDict = {}
tipDict = {}
colThickDict = defaultdict(list) # To read it colThickDict[row][0][element]

for lines in linesValues    :
    line = lines.split('\t')
    # print(line)
    # print('select\t', line[idColValColumn])
    # print('is_number\t', is_number(line[idColValColumn]))
    if is_number(line[idColValColumn]) :
        tipDict[line[0]] = float(line[idColValColumn])
        colDict[line[0]] = abs(float(line[idColValColumn]))
        colThickDict[line[0]] = abs(float(line[idThickValColumn]))

# print(colDict)
# print(tipDict)


[mapBoundaryColValue, minColVal, maxColVal, low] = findMapBound(colDict.values())

# print(mapBoundaryColValue, minColVal, maxColVal, low)
[mapBoundaryColValue, minColVal, maxColVal, low] = findMapBound(tipDict.values())
# print(mapBoundaryColValue, minColVal, maxColVal, low)

[mapBoundaryThickValue, minThickVal, maxThickVal, low] = findMapBound(colThickDict.values())
# print(mapBoundaryThickValue, minThickVal, maxThickVal, low)

minColVal = float(sys.argv[6])
maxColVal = float(sys.argv[7])
low = float(sys.argv[8])

### Setting up the color mapping.

# palette = Cubehelix.make(start=1, rotation=1, n=256, gamma=0.5)

colorbar = str(sys.argv[9])
#jet = cm = plt.get_cmap('viridis_r')
#jet = cm = plt.get_cmap('gist_rainbow_r')
jet = cm = plt.get_cmap(colorbar)

norm = mpl.colors.LogNorm(vmin = 0.1, vmax = (abs(maxColVal)), clip = False)
m = cmx.ScalarMappable(norm=norm, cmap=jet)

# a more convinient way to import the svg file
root = ET.parse(open(fileMapName, 'r'))


for element in root.iter():
    if element.tag.split("}")[1] == 'title':
        element.text = "Heterologous pathways Map" # element.text contains the text among <>text<\>
        ## print element.text

### For every value in the dictionary look for any occurence in the
### map and change its color.

thickness = minThick
count = 0
for key, value in colDict.items():
    count = count + 1
    i = 0
    ## print "<-----", key, value
    if value == float('Inf'):
        color = mpl.colors.rgb2hex(m.to_rgba(maxColVal))
    elif value == -float('Inf'):
        color = mpl.colors.rgb2hex(m.to_rgba(minColVal))
    elif value == 0:
        color = mpl.colors.rgb2hex(gray)
    else:
        color = mpl.colors.rgb2hex(m.to_rgba(value))

    for element in root.iter():
        eTag = element.tag.split("}")[1]
        ## print "-Tag- ", eTag
        if (eTag == 'line') or (eTag == 'path') :
            eID = element.get("id")
            ## print "      eID- ", eID
            R_ID = "R_"+key
            ## print R_ID
            if (eID != None) and (eID == R_ID): # and (eID[len(eID)-1] != "_"):
                                                # extract those elements whose attribute id is ...
                ## print "-->R_ID ", R_ID
                ## print "      eID- ", eID
                styleStr =  element.get("style")
                ## print "      Style- ", styleStr
                element.set("style", fixStrokeColor(styleStr, color))


## Add FORWARD color

count = 0
for key, value in tipDict.items():
    count = count + 1
    i = 0
    ## print "<-----", key, value
    if sign_0(value) == 1.0:
        color = mpl.colors.rgb2hex(m.to_rgba(abs(value)))
    elif value == float('Inf'):
        color = mpl.colors.rgb2hex(m.to_rgba(maxColVal))
    elif value == -float('Inf'):
        color = mpl.colors.rgb2hex(m.to_rgba(minColVal))
    else:
        color = mpl.colors.rgb2hex(gray)

    for element in root.iter():
        eTag = element.tag.split("}")[1]
        ## print "-Tag- ", eTag
        if (eTag == 'line') or (eTag == 'path') :
            eID = element.get("id")
            # print "      eID- ", eID
            R_ID = "F_" + key
            ## print R_ID
            if (eID != None) and (eID == R_ID): # and (eID[len(eID)-1] != "_"):
                                        # extract those elements whose attribute id is ...
                ## print "-->R_ID ", R_ID
                ## print "      eID- ", eID
                styleStr =  element.get("style")
                ## print "      Style- ", styleStr
                element.set("style", fixStrokeColor(styleStr, color))


# aggiungo colore BACKWARD
count = 0
for key, value in tipDict.items():
    count = count + 1
    i = 0
    ## print "<-----", key, value
    if sign_0(value) == -1.0:
        color = mpl.colors.rgb2hex(m.to_rgba(abs(value)))
    elif value == float('Inf'):
        color = mpl.colors.rgb2hex(m.to_rgba(maxColVal))
    elif value == -float('Inf'):
        color = mpl.colors.rgb2hex(m.to_rgba(minColVal))
    else:
        color = mpl.colors.rgb2hex(gray)

    for element in root.iter():
        eTag = element.tag.split("}")[1]
        ## print "-Tag- ", eTag
        if (eTag == 'line') or (eTag == 'path') :
            eID = element.get("id")
            ## print "      eID- ", eID
            R_ID = "B_" + key
            ## print R_ID
            if (eID != None) and (eID == R_ID): # and (eID[len(eID)-1] != "_"):
                                        # extract those elements whose attribute id is ...
                ## print "-->R_ID ", R_ID
                ## print "      eID- ", eID
                styleStr =  element.get("style")
                ## print "      Style- ", styleStr
                element.set("style", fixStrokeColor(styleStr, color))


# Thickness and colour are in two distinct dictionary because the
# dictionaries can have different reactions.

count = 0
for key, value in colThickDict.items():
    count = count + 1
    i = 0
    ## print "<-----", key, value
    if value > epsilon:
        if value == float('Inf'):
            thickness = maxThikc
        else:
            thickness = thicknessMapping(value, maxThickVal)
    elif value < -epsilon:
        if value == -float('Inf'):
            thickness = maxThikc
        else:
            thickness = thicknessMapping(abs(value), maxThickVal)
    else:
        thickness = minThick

    for element in root.iter():
        eTag = element.tag.split("}")[1]
        if (eTag == 'line') or (eTag == 'path') :
            eID = element.get("id")
            R_ID = "R_"+key
            if (eID != None) and (eID == R_ID): # and (eID[len(eID)-1] != "_"):
                                        # extract those elements whose attribute id is ...
                ## print "-->R_ID ", R_ID
                ## print "      eID- ", eID
                styleStr =  element.get("style")
                element.set("style", fixStrokewidth(styleStr, str(thickness)))


# Stroke

count = 0
for key, value in colThickDict.items():
    count = count + 1
    i = 0
    ## print "<-----", key, value
    if value == 0:
        dash = "5,5"
    else:
        dash = "none"

    for element in root.iter():
        eTag = element.tag.split("}")[1]
        if (eTag == 'line') or (eTag == 'path') :
            eID = element.get("id")
            R_ID = "R_"+key
            if (eID != None) and (eID == R_ID): # and (eID[len(eID)-1] != "_"):
                                        # extract those elements whose attribute id is ...
                ## print "-->R_ID ", R_ID
                ## print "      eID- ", eID
                styleStr =  element.get("style")
                element.set("style", fixStrokeDasharray(styleStr, str(dash)))


# Backward stroke.

count = 0
for key, value in tipDict.items():
    count = count + 1
    i = 0
    ## print "<-----", key, value
    if sign_0(value) == 1.0:
        dash = "1,2"
    elif sign_0(value) == 0.0:
        dash = "1,2"
    else:
        dash = "none"

    for element in root.iter():
        eTag = element.tag.split("}")[1]
        ## print "-Tag- ", eTag
        if (eTag == 'line') or (eTag == 'path') :
            eID = element.get("id")
            ## print "      eID- ", eID
            R_ID = "B_"+key
            ## print R_ID
            if (eID != None) and (eID == R_ID): #and (eID[len(eID)-1] != "_"):
                                        # extract those elements whose attribute id is ...
                ## print "-->R_ID ", R_ID
                ## print "      eID- ", eID
                styleStr =  element.get("style")
                element.set("style", fixStrokeDasharray(styleStr, str(dash)))


# Forward stroke.

count = 0
for key, value in tipDict.items():
    count = count + 1
    i = 0
    # print "<-----", key, value
    if sign_0(value) == -1.0:
        dash = "1,2"
    elif sign_0(value) == 0.0:
        dash = "1,2"
    else:
        dash = "none"

    for element in root.iter():
        eTag = element.tag.split("}")[1]
        ## print "-Tag- ", eTag
        if (eTag == 'line') or (eTag == 'path') :
            eID = element.get("id")
            ## print "      eID- ", eID
            R_ID = "F_" + key
            ## print R_ID
            if (eID != None) and (eID == R_ID): # and (eID[len(eID)-1] != "_"):
                                        # extract those elements whose attribute id is ...
                # print "-->R_ID ", R_ID
                # print "      eID- ", eID
                styleStr =  element.get("style")
                element.set("style", fixStrokeDasharray(styleStr, str(dash)))


defsFlag = 0
for element in root.iter():
    eTag = element.tag.split("}")[1]
    if (eTag == 'defs'): # The legend already exists then modify stops values.
        defsFlag = 1
        iIndex = element

# print "defsFlag = ", defsFlag

# LG = 'linearGradient'
# linearGradient = ET.SubElement(iIndex, ET.QName(SVG, LG) , id="Gradient", x1="0", x2="1", y1="1", y2="1", nsmap=NSMAP)
# color = mpl.colors.rgb2hex(m.to_rgba(minColVal))
# lgStop = ET.SubElement(linearGradient, ET.QName(SVG, "stop"), offset="0", nsmap=NSMAP)
# lgStop.set("stop-color", color)
# color = mpl.colors.rgb2hex(m.to_rgba(0))
# lgStop = ET.SubElement(linearGradient, ET.QName(SVG, "stop"), offset="0.5", nsmap=NSMAP)
# lgStop.set("stop-color", color)
# color = mpl.colors.rgb2hex(m.to_rgba(maxColVal))
# lgStop = ET.SubElement(linearGradient, ET.QName(SVG, "stop"), offset="1", nsmap=NSMAP)
# lgStop.set("stop-color", color)


# for element in root.iter():
#     eTag = element.tag.split("}")[1]
#     print eTag
#     if(eTag == "g"):
#         eID = element.get("id")
#         if(eID == "layer1"):
#             layer = element


for element in root.iter():
    eTag = element.tag.split("}")[1]
    # print eTag
    if(eTag == "svg"):
        svgEl = element


perc50 = round(10**(math.log10(maxColVal - low) * 0.5), 3)
perc25 = round(10**(math.log10(maxColVal - low) * 0.25), 3)
perc75 = round(10**(math.log10(maxColVal - low) * 0.75), 3)

# COLORBAR ELEMENTS

# legenda = ET.SubElement(svgEl, ET.QName(SVG, "rect") , id="legenda", x="10", y="40", rx="15", ry="15", width="1000", height="20", fill="url(#Gradient)", nsmap=NSMAP)


legTextStart = ET.SubElement(svgEl, ET.QName(SVG, "text") ,
                                 transform = "matrix(1 0 0 1 5 5)",
                                 nsmap = NSMAP)
legTextStart.text = str(round(low, 3))
legTextStart.set("font-family", "'MyriadPro-Regular'")
legTextStart.set("font-size", "4")

legTextStart = ET.SubElement(svgEl, ET.QName(SVG, "text"),
                                 transform = "matrix(1 0 0 1 10 10)",
                                 nsmap = NSMAP)
legTextStart.text = "|" #str(round(low, 3))
legTextStart.set("font-family", "'MyriadPro-Regular'")
legTextStart.set("font-size", "4")

legTextMiddle = ET.SubElement(svgEl, ET.QName(SVG, "text"),
                                  transform = "matrix(1 0 0 1 45 5)",
                                  nsmap = NSMAP)
legTextMiddle.text = str(perc25)
legTextMiddle.set("font-family", "'MyriadPro-Regular'")
legTextMiddle.set("font-size", "4")

legTextMiddle = ET.SubElement(svgEl, ET.QName(SVG, "text"),
                                  transform = "matrix(1 0 0 1 50 10)",
                                  nsmap = NSMAP)
legTextMiddle.text = "|" #str(perc25)
legTextMiddle.set("font-family", "'MyriadPro-Regular'")
legTextMiddle.set("font-size", "4")

legTextMiddle = ET.SubElement(svgEl, ET.QName(SVG, "text"),
                                  transform = "matrix(1 0 0 1 89 5)",
                                  nsmap = NSMAP)
legTextMiddle.text = str(perc50)
legTextMiddle.set("font-family", "'MyriadPro-Regular'")
legTextMiddle.set("font-size", "4")

legTextMiddle = ET.SubElement(svgEl, ET.QName(SVG, "text"),
                                  transform = "matrix(1 0 0 1 94 10)",
                                  nsmap = NSMAP)
legTextMiddle.text = "|" #str(perc50)
legTextMiddle.set("font-family", "'MyriadPro-Regular'")
legTextMiddle.set("font-size", "4")

legTextMiddle = ET.SubElement(svgEl, ET.QName(SVG, "text"),
                                  transform = "matrix(1 0 0 1 135 5)",
                                  nsmap = NSMAP)
legTextMiddle.text = str(perc75)
legTextMiddle.set("font-family", "'MyriadPro-Regular'")
legTextMiddle.set("font-size", "4")

legTextEnd = ET.SubElement(svgEl, ET.QName(SVG, "text"),
                               transform = "matrix(1 0 0 1 142 10)",
                               nsmap = NSMAP)
legTextEnd.text = "|" #str(maxColVal)
legTextEnd.set("font-family", "'MyriadPro-Regular'")
legTextEnd.set("font-size", "4")

legTextEnd = ET.SubElement(svgEl, ET.QName(SVG, "text"),
                               transform = "matrix(1 0 0 1 182 5)",
                               nsmap = NSMAP)
legTextEnd.text = str(round(maxColVal, 3))
legTextEnd.set("font-family", "'MyriadPro-Regular'")
legTextEnd.set("font-size", "4")

legTextEnd = ET.SubElement(svgEl, ET.QName(SVG, "text"),
                               transform = "matrix(1 0 0 1 190 10)",
                               nsmap = NSMAP)
legTextEnd.text = "|" #str(maxColVal)
legTextEnd.set("font-family", "'MyriadPro-Regular'")
legTextEnd.set("font-size", "4")

## l'argomento transform indica le coordinate dell'elemento che stiamo aggiungendo
legTextFlux = ET.SubElement(svgEl, ET.QName(SVG, "text"),
                                transform = "matrix(1 0 0 1 220 10)",
                                nsmap = NSMAP)

legTextFlux.text = "---- Dashed line: flux = 0"
legTextFlux.set("font-family", "'MyriadPro-Regular'")
legTextFlux.set("font-size", "4")
legTextFlux.set("fill", "gray")

titleMap =  sys.argv[10] + " - " +  sys.argv[11] + " - " +  sys.argv[12] + " + " +  sys.argv[13] + " - Time point: " +  sys.argv[14] + "h - sf: " +  sys.argv[15]

## l'argomento transform indica le coordinate dell'elemento che stiamo aggiungendo
legTextFlux = ET.SubElement(svgEl, ET.QName(SVG, "text"),
                                transform = "matrix(1 0 0 1 220 15)",
                                nsmap = NSMAP)

legTextFlux.text = titleMap
legTextFlux.set("font-family", "'MyriadPro-Regular'")
legTextFlux.set("font-size", "3")


# <rect id="rect2" x="10" y="20" rx="15" ry="15" width="500" height="20" fill="url(#Gradient)"/>
# <text transform = "matrix(1 0 0 1 10 15)" font-family="'MyriadPro-Regular'" font-size="12">start</text>
# <text transform = "matrix(1 0 0 1 260 15)" font-family="'MyriadPro-Regular'" font-size="12">middle</text>
# <text transform = "matrix(1 0 0 1 510 15)" font-family="'MyriadPro-Regular'" font-size="12">end</text>

# legTextMiddle = ET.SubElement(svgEl, ET.QName(SVG, "text"), transform = "matrix(1 0 0 1 500 30)", nsmap = NSMAP)
# legTextMiddle.text = str(0)
# legTextMiddle.set("font-family", "'MyriadPro-Regular'")
# legTextMiddle.set("font-size", "24")

# print "-----------------"
# for child in iIndex:
#     print child.tag
# print "-----------------"


lista = list(root.iter())
# print lista

# scrittura output


with open(filename, 'wb') as newMapFile:
    newMapFile.write(ET.tostring(root, pretty_print=True))
newMapFile.close

#### end of file -- colorFluxes.py
