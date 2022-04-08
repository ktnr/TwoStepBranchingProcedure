from asyncio.windows_events import NULL
import json
import math

import os, sys
import re

def ReadJsonAsDict(fileName):
    data = {}
    with open(fileName) as f:
        data = json.load(f)

    return data

def ReadUniboFormatAsArray(fileName):
    with open(fileName) as f:
        lines = f.readlines()

        return lines

class Bin:
    def __init__(self, dx, dy):
        self.Dx = dx
        self.Dy = dy
        self.Area = dx * dy

class Item:
    def __init__(self, externId, dx, dy):
        self.X = 0
        self.Y = 0
        self.Dx = dx
        self.Dy = dy
        self.Area = dx * dy
        self.ExternId = externId

def ConvertFromUniboFormat(lines, fileName):
    name = fileName.split(".")[0]

    numberOfItems = int(lines[0].replace("\n", ""))
    """
    feasibility = ""
    if "F" in name:
        feasibility = "F"
    elif "N" in name:
        feasibility = "N"
    elif "X" in name:
        feasibility = "U"
    else:
        raise ValueError("Problem instance status missing in file name.")

    epsilonFileName = int(name[1:3])
    numberOfItemsFileName = int(name[-2:])

    if numberOfItemsFileName != numberOfItems:
        raise ValueError(f"Input file name declares {numberOfItemsFileName} but the file itself {numberOfItems}.")
    """

    secondLineSplitted = lines[1].replace("\n", "").split(" ")
    containerDx, containerDy = int(secondLineSplitted[0]), int(secondLineSplitted[1])

    newContainer = Bin(containerDx, containerDy, 1)

    areaSum = 0
    newItems = []
    for i, itemLine in enumerate(lines[2:]):
        splittedLine = itemLine.replace("\n", "").split(" ")

        dx = int(splittedLine[1])
        dy = int(splittedLine[2])

        demand = int(splittedLine[3])
        copies = int(splittedLine[4])
        profit = int(splittedLine[5])

        if demand > 1 or copies > 1 or profit > 0:
            raise ValueError("Item has demand {demand}, copies {copies}, and proft {profit}.")

        areaSum += dx * dy
        newItems.append(Item(i, dx, dy, 1, 1))

    numberOfItemTypes = len(newItems)
    if numberOfItemTypes != numberOfItems:
        raise ValueError(f"Input file has {numberOfItemTypes} but {numberOfItems} were expected.")

    """
    epsilonCalculated = 100 * (1.0 - (float(areaSum) / float(containerDx * containerDy)))
    epsilonFloor = math.floor(epsilonCalculated)
    epsilonCeiling = math.ceil(epsilonCalculated)
    if epsilonFileName not in range(epsilonFloor, epsilonCeiling + 1):
        print(f"Stated epsilon = {epsilonFileName} not in [{epsilonFloor}, {epsilonCeiling}] = calculated epsilon")
    """

    newJsonDict = {}

    newJsonDict["Name"] = name
    newJsonDict["InstanceType"] = "2D-OPP"
    newJsonDict["NumberItemTypes"] = numberOfItemTypes

    newJsonDict["Container"] = newContainer
    newJsonDict["ItemTypes"] = newItems

    return newJsonDict

def ConvertFromDictionary(jsonDict, fileName):
    name = jsonDict["Name"]
    items = jsonDict["ItemTypes"]
    container = jsonDict["Container"]

    containerDx = container["Length"]
    containerDy = container["Width"]
    #containerArea = container["Cost"]

    #if containerArea != containerDx * containerDy:
    #    raise ValueError("Wrong container area in input data.")

    newContainer = Bin(containerDx, containerDy)

    newItems = []
    for i, item in enumerate(items):
        dx = int(item["Length"])
        dy = int(item["Width"])

        """
        demand = item["DemandMax"]
        if demand != None and int(demand) > 1:
            raise ValueError("Item has demand > 1, (weakly) heterogeneous instance.")

        area = int(item["Value"])
        if area != dx * dy:
            #print(f"{fileName} ({name}): Wrong item area in input data: {area} != {dx} * {dy} = {dx * dy}")
            raise ValueError(f"{fileName} ({name}): Wrong item area in input data: {area} != {dx} * {dy} = {dx * dy}")
        """

        newItems.append(Item(i, dx, dy))

    numberOfItemTypes = len(newItems)

    newJsonDict = {}

    newJsonDict["Name"] = name
    #newJsonDict["InstanceType"] = "2D-OPP"
    #newJsonDict["NumberItemTypes"] = numberOfItemTypes

    newJsonDict["Container"] = newContainer
    newJsonDict["Items"] = newItems

    return newJsonDict

def main():
    sourcePath = "data/inputCLP/MSB-630_Unibo"
    targetPath = "data/input/MSB-630_Unibo"
    for root, dirs, files in os.walk(sourcePath):
        print(root)
        for fileName in files:
            if fileName.endswith(".json"):
                print(f"Processing file {fileName}")

                data = ReadJsonAsDict(os.path.join(sourcePath, fileName))
                #data = ReadUniboFormatAsArray(os.path.join(sourcePath, fileName))

                convertedFile = {}
                try:
                    convertedFile = ConvertFromDictionary(data, fileName)
                    #convertedFile = ConvertFromUniboFormat(data, fileName)
                except ValueError as ex:
                    print(ex)

                newFileName = fileName.split(".")[0] + ".json"
                with open(os.path.join(targetPath, newFileName), "w") as fp:
                    json.dump(convertedFile, fp, indent=4, default=lambda x: x.__dict__)

if __name__ == "__main__":
    main()