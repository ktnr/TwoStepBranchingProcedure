import json

import math

import os, sys
import re

import matplotlib
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

def ReadJsonAsDict(fileName):
    data = {}
    with open(fileName) as f:
        data = json.load(f)

    return data


def PlotSolution(solution):
    fig = plt.figure()
    ax = fig.add_subplot(111)
    
    W = max(solution["PlacedAreaVector"])
    H = len(solution["DeactivatedAreaVector"]) - 1

    ax.set_xlim((0, W))
    ax.set_ylim((0, H))
    ax.set_aspect('equal')

    items = solution["Items"]

    itemArea = 0.0
    for i, item in enumerate(items):
        x1 = item["X"]
        y1 = item["Y"]

        dx = item["Dx"]
        dy = item["Dy"]
        x2 = x1 + dx
        y2 = y1 + dy

        itemArea += dx * dy

        w = (x2 - x1)
        h = (y2 - y1)

        color = matplotlib.colors.to_hex([ (x2 - x1) / W, (y2 - y1) / H, (w * h)/ (W * H) ])

        rectPlot = matplotlib.patches.Rectangle((x1,y1), x2 - x1, y2 - y1, color=color)
        ax.add_patch(rectPlot)

        ax.annotate(i, (x1 + ((x2 - x1) / 2.0), y1 + ((y2 - y1) / 2.0)), color='b', weight='bold', 
                fontsize=6, ha='center', va='center')

    if itemArea > W * H:
        print(f"Item area {itemArea} > {W * H}")
        raise ValueError("Item area is larger than container area.")

    plt.show()

def main():
    path = "Output/2022-01-21_19-23-48/"
    fileName = "Solution.json"

    jsonDict = ReadJsonAsDict(path + fileName)

    PlotSolution(jsonDict)

if __name__ == "__main__":
    main()
