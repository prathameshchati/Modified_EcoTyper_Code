import pandas as pd
import numpy as np
import seaborn as sns
import subprocess
import os
import scanpy as sc
import random
import sys

arguments = sys.argv

metadataPath = arguments[1]

rootFolder = "../Generate_Cell_State_Abundances_Predefined_States_Test_Outputs/"

fullMetadataPath = rootFolder + metadataPath

metadata = pd.read_csv(fullMetadataPath, delimiter = "\t")
print("Metadata read successfully")

sampleColumn = arguments[2]
CellTypeColumn = arguments[3]
cellStateColumn = arguments[4]

rootOutput = "/duo4/users/pchati/ecotyper/Generate_Cell_State_Abundances_Predefined_States_Test_Outputs/"

outputPath = arguments[5]

outputPath = rootOutput + outputPath

def generateH(metadata, sampleColumn, CellTypeColumn, cellStateColumn, outputPath):
    
    cellTypes = set(metadata[CellTypeColumn])
    sampleNames = set(metadata[sampleColumn])
    
    for cellType in cellTypes:
        
        print("Generating state abundance matrix for cell type:", cellType)
        
        cellTypeMetadata = metadata[metadata[CellTypeColumn] == cellType]
        
        cellTypeStates = set(cellTypeMetadata[cellStateColumn])

        cellTypeStatesRows = []
        for state in cellTypeStates:
#             cellTypeStateName = cellType + "_" + "S" + str(state)
#             cellTypeStatesRows.append(cellTypeStateName)
              cellTypeStatesRows.append(state)

        cellTypeH = pd.DataFrame(columns = sampleNames, index = cellTypeStatesRows)
    
        for sample in sampleNames:
            cellTypeSampleMetadata = cellTypeMetadata[cellTypeMetadata[sampleColumn] == sample]
            cellTypeStateAbudnaceList = []
            numCellType = len(cellTypeSampleMetadata)
            for state in cellTypeStates:
                cellTypeSampleStateMetadata = cellTypeSampleMetadata[cellTypeSampleMetadata[cellStateColumn] == state]
                numCellTypeState = len(cellTypeSampleStateMetadata)
                if (numCellTypeState == 0):
                    cellTypeStateAbudnaceList.append(0)
                else: 
                    stateAbundance = numCellTypeState/numCellType
                    cellTypeStateAbudnaceList.append(stateAbundance)

            cellTypeH[sample] = cellTypeStateAbudnaceList

        cellTypeH.fillna(0)
        
        outputCellTypePath = outputPath + "/" + cellType
        
        directoryExists = os.path.exists(outputCellTypePath)
        
        if (directoryExists):
            print(directoryExists)
        else:
            mkdirCellType = "mkdir " + outputCellTypePath
            subprocess.call(mkdirCellType, shell=True)
        
        saveFile = outputCellTypePath + "/" + "state_abundances.txt"
        
        cellTypeH.to_csv(saveFile, sep = "\t")
        

generateH(metadata, sampleColumn, CellTypeColumn, cellStateColumn, outputPath)