import matplotlib.pyplot as plt
import os
import numpy as np
import utils.genericLib as gL
from lxml import etree as ET

FIGUREDIR = gL.dDirs["figs"]

def plotODvsBiomasslux(lX, lOD, lBIOM, outputFigName):
    fig = plt.figure(figsize=(13.3, 10))
    ax = fig.add_subplot(111)

    ax.scatter(
        lX, lOD, alpha=0.7, marker="o", s=300, label="OD / OD(t0)"
    )

    ax.scatter(
        lX, lBIOM,
        alpha=0.7,
        marker="x",
        s=300,
        label="Biomass flux / Biomass flux(t0)",
    )


    plt.xlabel("Time (h)", size=18)
    plt.xticks(fontsize=18)
    plt.yticks(fontsize=18)
    plt.legend(fontsize=18)
    plt.savefig(os.path.join(FIGUREDIR, outputFigName + ".png"))
    plt.close()

def colorRxns(tree, lCandidateRxns, color):
    root = tree.getroot()
    for element in root.getiterator():
        if element.tag.endswith('path'):
            pathTag = element.tag
    for element in root.getiterator(pathTag):
        if element.attrib["id"][2:] in lCandidateRxns:
            lStyle = element.attrib["style"].split(";")
            lStyleUpdate = []
            for styleElement in lStyle:
                if styleElement.split(":")[0] == "stroke":
                    lStyleUpdate.append(styleElement.split(":")[0] + ":" + color)
                elif styleElement.split(":")[0] == "stroke-width":
                    lStyleUpdate.append(styleElement.split(":")[0] + ":3px")
                else:
                    lStyleUpdate.append(styleElement.split(":")[0] + ":" + styleElement.split(":")[1])
            element.attrib["style"] = ";".join(lStyleUpdate)
    return tree

def plotHeatmap(df, outputFileName, titleFig):
    masked_data = np.ma.masked_invalid(df.values) # np.ma.masked_invalid(data) creates a mask for NaN values
    #in the data. When we pass masked_data to plt.imshow(), the NaN values are excluded from the heatmap.

    # By default, Matplotlib includes NaN values in the colorbar legend. To exclude them,
    # we can use the set_bad function of the colormap.
    # Create a colormap and set NaN values to be transparent
    cmap = plt.cm.viridis
    cmap.set_bad(color='none')
    plt.imshow(masked_data, cmap=cmap)
    # Displaying a color bar to understand which color represents which range of data
    plt.colorbar()
    # Assigning labels of x-axis according to dataframe
    plt.xticks(range(df.shape[1]), df.columns, fontsize=6, rotation = 45)
    # Assigning labels of y-axis according to dataframe
    plt.yticks(range(df.shape[0]), df.index, fontsize=2)
    plt.title(titleFig, fontsize = 10)
    # Displaying the figure
    plt.savefig(outputFileName, format = "pdf")
    plt.close("all")

def heatmap_deletions(data, outputFileName, titleFig):
    data.fillna(0, inplace = True)
    cmap = plt.cm.viridis
    cmap.set_bad(color='none')

    fig,ax = plt.subplots(figsize=(20,20))
    im = ax.imshow(data,cmap = cmap)
    ax.set_aspect(aspect="auto")
    cbar = ax.figure.colorbar(im, ax=ax)
    cbar.ax.tick_params(labelsize=18)
    cbar.ax.set_ylabel("% Biomass reduction", fontsize = 20, rotation=-90, va="bottom")
    ax.set_xticks(range(data.shape[1]))
    ax.set_xticklabels(data.columns, fontsize=15, rotation = 45)
    ax.set_yticks(range(data.shape[0]))
    if data.shape[0] < 10:
        ax.set_yticklabels(data.index, fontsize=10)
    elif data.shape[0] >= 10 and  data.shape[0] < 15:
        ax.set_yticklabels(data.index, fontsize=10)
    else:
        ax.set_yticklabels(data.index, fontsize=6)
    ax.set_title(titleFig, fontsize = 18)
    plt.savefig(outputFileName, format = "pdf")
    plt.close("all")

def heatmap(data, outputFileName, titleFig):
    masked_data = np.ma.masked_invalid(data.values)
    cmap = plt.cm.viridis
    cmap.set_bad(color='none')

    fig,ax = plt.subplots(figsize=(20,20))
    # x0 = data.shape[1] - 0.5
    # x1 = x0 + data.shape[1]
    # ax.imshow(data,cmap = cmap, extent=[x0, x1, -0.5, a.shape[0] - 0.5])
    im = ax.imshow(data,cmap = cmap)
    ax.set_aspect(aspect="auto")
    cbar = ax.figure.colorbar(im, ax=ax)
    cbar.ax.set_ylabel("Percentage variation", fontsize = 16, rotation=-90, va="bottom")
    ax.set_xticks(range(data.shape[1]))
    ax.set_xticklabels(data.columns, fontsize=15, rotation = 45)
    ax.set_yticks(range(data.shape[0]))
    if data.shape[0] < 10:
        ax.set_yticklabels(data.index, fontsize=12)
    elif data.shape[0] >= 10 and  data.shape[0] < 15:
        ax.set_yticklabels(data.index, fontsize=8)
    else:
        ax.set_yticklabels(data.index, fontsize=4)
    ax.set_title(titleFig, fontsize = 18)
    plt.savefig(outputFileName, format = "pdf")
    plt.close("all")
