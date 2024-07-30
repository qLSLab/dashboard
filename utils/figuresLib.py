import matplotlib.pyplot as plt
import os
import utils.genericLib as gL
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
