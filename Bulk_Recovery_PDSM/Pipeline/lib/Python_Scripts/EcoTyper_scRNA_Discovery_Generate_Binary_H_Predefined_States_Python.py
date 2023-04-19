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

metadata = pd.read_csv(fullMetadataPath, delimiter = "\t", index_col = 0)
print("Metadata read successfully")

sampleColumn = arguments[2]
CellTypeColumn = arguments[3]
cellStateColumn = arguments[4]

rootOutput = "/duo4/users/pchati/ecotyper/Generate_Cell_State_Abundances_Predefined_States_Test_Outputs/"

outputPath = arguments[5]

outputPath = rootOutput + outputPath


def generateBinaryH(metadata, sampleColumn, CellTypeColumn, cellStateColumn, outputPath):
    
    cellTypes = set(metadata[CellTypeColumn])
#     cellTypes = ["Malignant"]
    
    for cellType in cellTypes:
        
        cellTypeMetadata = metadata[metadata[CellTypeColumn] == cellType]
        
        cellTypeStates = set(cellTypeMetadata[cellStateColumn])

        cellTypeStatesRows = []
        for state in cellTypeStates:
            # cellTypeStateName = cellType + "_" + "S" + str(state)
            # cellTypeStatesRows.append(cellTypeStateName)
            cellTypeStatesRows.append(state)

            
        cellIDs = list(cellTypeMetadata.index.values)
#         cellIDs = list(cellTypeMetadata["ID"])


#         cellTypeBinaryH = pd.DataFrame(columns = cellIDs, index = cellTypeStatesRows)
#         cellTypeBinaryH = cellTypeBinaryH.T

        cellTypeBinaryH = pd.DataFrame(index = cellIDs)
    
        cellTypeBinaryH.index.name = "Index"

        for state, stateName in zip(cellTypeStates, cellTypeStatesRows):
            cellTypeStateMetadata = cellTypeMetadata[cellTypeMetadata[cellStateColumn] == state]
            
            listOfOnes = [1]*len(cellTypeStateMetadata)
            
            cellTypeStateMetadata[stateName] = listOfOnes
            
#             cellTypeStateCellIDs = list(cellTypeStateMetadata.index.values)
            cellTypeStateCellIDs = cellTypeStateMetadata[[stateName]]
            
            cellTypeBinaryH = pd.concat([cellTypeBinaryH, cellTypeStateCellIDs], axis=1).fillna(0)

            
#             binaryAssignment = []
            
#             for cID in cellIDs:
#                 if (cID in cellTypeStateCellIDs):
#                     binaryAssignment.append(1)
#                 else:
#                     binaryAssignment.append(0)
            
    
        cellTypeBinaryH = cellTypeBinaryH.T
        
#         return cellTypeBinaryH

        outputCellTypePath = outputPath + "/" + cellType
        
        directoryExists = os.path.exists(outputCellTypePath)
        
        if (directoryExists):
            print(directoryExists)
        else:
            mkdirCellType = "mkdir " + outputCellTypePath
            subprocess.call(mkdirCellType, shell=True)
        
        saveFile = outputCellTypePath + "/" + "Binary_H.txt"
        
        cellTypeBinaryH.to_csv(saveFile, sep = "\t")
        

generateBinaryH(metadata, sampleColumn, CellTypeColumn, cellStateColumn, outputPath)