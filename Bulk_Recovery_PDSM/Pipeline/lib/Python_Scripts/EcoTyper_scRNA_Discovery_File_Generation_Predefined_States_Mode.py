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

def generateMappingToInitialState(metadata, sampleColumn, CellTypeColumn, cellStateColumn, outputPath):

    cellTypes = set(metadata[CellTypeColumn])
#     cellTypes = ["MonocyteDC"]
    
    for cType in cellTypes:
        cellTypeMetadata = metadata[metadata[CellTypeColumn] == cType]
        cellStates = list(set(cellTypeMetadata[cellStateColumn]))
        mappingInitialStatesDF = pd.DataFrame()
        mappingInitialStatesDF["State"] = cellStates
        mappingInitialStatesDF["InitialState"] = cellStates
        
        # SAVE OUTPUTS
        cellTypeOutputPath = outputPath + "/" + cType + "/mapping_to_initial_states.txt"
        mappingInitialStatesDF.to_csv(cellTypeOutputPath, sep = "\t")

def generateInitialStateAssignment(metadata, sampleColumn, CellTypeColumn, cellStateColumn, outputPath):
    
    cellTypes = set(metadata[CellTypeColumn])
#     cellTypes = ["MonocyteDC"]
    
    for cType in cellTypes:
        cellTypeMetadata = metadata[metadata[CellTypeColumn] == cType]
        cellTypeMetadata = cellTypeMetadata.set_index("ID")
#         cellTypeMetadata.index.name = "ID"
        cellTypeInitialStateAssignmentDF = cellTypeMetadata[[cellStateColumn]] 
        cellTypeInitialStateAssignmentDF = cellTypeInitialStateAssignmentDF.rename(columns = {cellStateColumn: "State"})
        cellTypeInitialStateAssignmentDF = cellTypeInitialStateAssignmentDF.sort_values(by = "State")

        # SAVE OUTPUTS
        cellTypeOutputPath = outputPath + "/" + cType + "/initial_state_assignment.txt"
        cellTypeInitialStateAssignmentDF.to_csv(cellTypeOutputPath, sep = "\t")
        
        # GENERATE STATE ASSIGNMENT
        initalStateList = list(cellTypeInitialStateAssignmentDF["State"])
        cellTypeInitialStateAssignmentDF["InitialState"] = initalStateList
        
        # SAVE OUTPUTS
        cellTypeOutputPath = outputPath + "/" + cType + "/state_assignment.txt"
        cellTypeInitialStateAssignmentDF.to_csv(cellTypeOutputPath, sep = "\t")


        
        
# def generateStateAssignment(metadata, outputPath):
#     cellTypes = set(metadata[CellTypeColumn])
    
#     for cType in cellTypes:
# #         print(cType)
#         initialStateAssignmentCellTypePath = outputPath + "/" + cType + "/initial_state_assignment.txt"
#         cellTypeInitialStateAssignmentDF = pd.read_csv(initialStateAssignmentCellTypePath, delimiter = "\t", index_col = 0)
#         initalStateList = list(cellTypeInitialStateAssignmentDF["State"])
#         cellTypeInitialStateAssignmentDF["InitialState"] = initalStateList
        
#         # SAVE OUTPUTS
#         cellTypeOutputPath = outputPath + "/" + cType + "/state_assignment.txt"
#         cellTypeInitialStateAssignmentDF.to_csv(cellTypeOutputPath, sep = "\t")
        
def generateRankData(metadata, sampleColumn, CellTypeColumn, cellStateColumn, outputPath):

    cellTypes = list(set(metadata[CellTypeColumn]))
#     cellTypes = ["MonocyteDC"]
    
    rankData = pd.DataFrame(index = cellTypes)
    
    rankList = []
    falseList = []
    for cType in cellTypes:
        cellTypeMetadata = metadata[metadata[CellTypeColumn] == cType]
        cellTypeStates = set(cellTypeMetadata[cellStateColumn])
        cellTypeRank = len(cellTypeStates)
        rankList.append(cellTypeRank)
        falseList.append('FALSE')
        
    rankData.index.name = "CellType"
    rankData["Chosen_Rank"] = rankList
    rankData["Redefined"] = falseList
    
    rankDataOuputPath = outputPath + "/" + "rank_data.txt"
    rankData.to_csv(rankDataOuputPath, sep = "\t")
    
    
generateMappingToInitialState(metadata, sampleColumn, CellTypeColumn, cellStateColumn, outputPath)

generateInitialStateAssignment(metadata, sampleColumn, CellTypeColumn, cellStateColumn, outputPath)

# generateStateAssignment(metadata, outputPath)

generateRankData(metadata, sampleColumn, CellTypeColumn, cellStateColumn, outputPath)
    
