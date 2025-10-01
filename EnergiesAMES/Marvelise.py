import pandas as pd
from pandarallel import pandarallel
pandarallel.initialize(progress_bar=True)

pd.set_option('display.float_format', '{:.6f}'.format)

marvelColumns = ["nu1", "nu2", "l2", "nu3", "e/f", "J", "Em", "Uncertainty", "Transitions"]

# marvelEnergies = pd.read_csv("14NH3-NewEnergies-MoreNbAssignments.txt", delim_whitespace=True, names=marvelColumns, dtype=str)
marvelEnergies = pd.read_csv("MarvelEnergies.txt", delim_whitespace=True, names=marvelColumns, dtype=str)
# marvelEnergies = marvelEnergies[marvelEnergies["J"] == "0"]
print(marvelEnergies.head(20).to_string(index=False))
# marvelEnergies = marvelEnergies[marvelEnergies["Transitions"].astype(int) > 1]

# symmetryMap = {
#     "A1'": "1",
#     "A2'": "2",
#     "E'": "3",
#     "A1\"": "4",
#     "A2\"": "5",
#     "E\"": "6"
# }

# marvelEnergies["Gamma"] = marvelEnergies["Gamma"].map(symmetryMap)
# marvelEnergies["inv"] = marvelEnergies["inv"].map(inversionMap)
marvelEnergies["J"] = marvelEnergies["J"].astype(int)
JMax = marvelEnergies["J"].max()
maxEnergy = marvelEnergies["Em"].astype(float).max()

# statesFileColumns = ["Gamma", "i", "E", "GammaRot", "J", "K", "Kc", "GammaVib", "n1", "n2", "n3", "Ci", "m1", "m2", "m3", "m4", "i2"]
# states = pd.read_csv("greppedEnergiesKrot.dat", delim_whitespace=True, names=statesFileColumns, dtype=str)
statesFileColumns = ["Gamma", "i", "E", "(1", "GammaRot", ";1", "J", "K", "Kc", ")1", "(2", "GammaVib", ";2", "n1", "n2", "n3", ")2",  "Ci", "(3", "m1", "m2", "m3", "m4", ")3", "(4", "i2", ")4"] # For them pesky (;)
states = pd.read_csv("CS2-Refit2.dat", delim_whitespace=True, names=statesFileColumns, dtype=str)
# statesFileColumns = ["i", "E", "g", "J", "Gamma", "n1", "n2", "n3", "GammaVib", "J'", "K'", "GammaRot", "Ci", "::", "i2", "m1", "m2", "m3", "i3"]
# states = pd.read_csv("/scratch/vol1/asmola/CS2/EnergiesAMES/CS2_AMES1_P48_02_step4_states_0-159.states", delim_whitespace=True, names=statesFileColumns, dtype=str)

print(states.head(20).to_string(index=False))
states["J"] = states["J"].astype(int)
states["E"] = states["E"].astype(float)
states = states[states["J"] <= JMax]
states = states[states["E"] <= maxEnergy + 100]
print(states.head(20).to_string(index=False))

def findParity(row):
    irrepParityMap = {"A1": 1, "A2": -1, "B1": -1, "B2": 1}
    if int(row["J"]) % 2 == 0:
        isEven = True
    else:
        isEven = False
    efParityMap = {False: {1: "f", -1: "e"}, True: {1: "e", -1: "f"}}
    row["e/f"] = efParityMap[isEven][irrepParityMap[row["Gamma"]]]
    return row

states = states.parallel_apply(lambda x:findParity(x), result_type="expand", axis=1)


print(states.head(20).to_string(index=False))

def findMatchingLevel(row, states):
    matchingStates = states[states["J"] == row["J"]]
    matchingStates["Obs-Calc"] = abs(float(row["Em"]) - states["E"])
    matchingStates = matchingStates.sort_values(by="Obs-Calc")
    matchingState = matchingStates.head(1).squeeze()
    row["TroveID"] = matchingState["i"]
    row["Etrove"] = matchingState["E"]
    row["Obs-Calc"] = matchingState["Obs-Calc"]
    return row

marvelEnergies = marvelEnergies.parallel_apply(lambda x:findMatchingLevel(x, states), result_type="expand", axis=1)
print(marvelEnergies.head(20).to_string(index=False))
marvelEnergiesReduced = marvelEnergies.dropna()
marvelEnergies = marvelEnergies.to_string(index=False, header=False)
marvelisedFile = "CS2-Matched-KRot5-REFIT2-AMES.out"
with open(marvelisedFile, "w+") as FileToWriteTo:
    FileToWriteTo.write(marvelEnergies)
print("New states file ready!")

marvelEnergiesReduced = marvelEnergiesReduced.to_string(index=False, header=False)
marvelisedFileReduced = "CS2-Matched-KRot5-REFIT2-AMES-Reduced.out"
with open(marvelisedFileReduced, "w+") as FileToWriteTo:
    FileToWriteTo.write(marvelEnergiesReduced)
print("New states file ready!")


marvelEnergies = generateTagColumn(marvelEnergies)
states = generateTagColumn(states)

columnsToConvertToInteger = ["J", "n1", "n2", "n3", "n4", "inversion"]
for column in columnsToConvertToInteger:
    states[column] = states[column].astype(int)

states = pd.merge(states, marvelEnergies.drop(["J", "Gamma", "Nb"], axis=1), on="Tag", how="left")

columnsToConvertToFloat = ["Em", "E", "Calc", "weight", "Uncertainty"]
for column in columnsToConvertToFloat:
    states[column] = states[column].astype(float)
    
pd.set_option('display.float_format', '{:.6f}'.format)

def marvelise(row):
    # Parameters that get used for the calculated uncertainty
    deltaB = 0.002
    deltaOmega = 0.3
    calculatedEnergy = row["Calc"]
    if row["Em"] >= 0 and abs(row["E"] - row["Em"]) < 1:
        row["Marvel"] = "Ma"
        marvelUncertainty = row["Uncertainty"]
        row["n1"] = row["nu1"]
        row["n2"] = row["nu2"]
        row["n3"] = row["nu3"]
        row["n4"] = row["nu4"] 
        row["l3"] = row["L3"]
        row["l4"] = row["L4"]
        row["inversion"] = row["inv"]
        row["K'"] = row["K"]
        row["weight"] = f"{marvelUncertainty:.{6}f}"
        marvelEnergy = row["Em"] 
        row["E"] = f"{marvelEnergy:.{6}f}" 
    else:
        row["Marvel"] = "Ca"
        estimatedUncertainty = deltaB*row["J"]*(row["J"] + 1) + deltaOmega*(row["n1"] + row["n2"] + row["n3"] + row["n4"])
        row["weight"] = f"{estimatedUncertainty:.{6}f}"
        row["E"] = f"{calculatedEnergy:.{6}f}"
    row["Calc"] = f"{calculatedEnergy:.{6}f}"
    return row

print("\n")
print("Beginning Marvelisation...")
states = states.parallel_apply(lambda x:marvelise(x), result_type="expand", axis=1)
states["inversion"] = states["inversion"].astype(int)
print("\n")
print("Marvelisation complete!")
statesFileColumns = ["i", "E", "g", "J", "weight", "p", "Gamma", "Nb", "n1", "n2", "n3", "n4", "l3", "l4", "inversion", "J'", "K'", "pRot", "v1", "v2", "v3", "v4", "v5", "v6", "GammaVib", "Marvel", "Calc"]
states = states[statesFileColumns]

columnWidth = 1
def formatColumns(value):
    global columnWidth
    return f'{value: >{columnWidth}}'
def reformatColumns(dataFrame, columnReformattingOptions):
    global columnWidth 
    collectedReformattedColumns = []
    for key in columnReformattingOptions.keys():
        columnWidth = key
        dataFrameToReformat = dataFrame[columnReformattingOptions[key]]
        collectedReformattedColumns += [dataFrameToReformat.applymap(formatColumns)]
    dataFrame = pd.concat(collectedReformattedColumns, axis=1, join="inner")
    return dataFrame

columnReformattingOptions = {
    12: ['i', 'E', 'weight', 'Calc'],
    6: ["g", "n1", "J'"],
    2: ["Gamma", "pRot"],
    10: ["Nb"],
    3: ["n2", "n3", "n4", "l4", "K'", "v2", "v3", "v4", "v5", "v6", "Marvel"],
    7: ["J", "v1"],
    4: ["inversion"],
    5: ["l3", "GammaVib"],
    1: ["p"]
}
print("\n")
print("Reformatting columns...")
states = reformatColumns(states, columnReformattingOptions)
print("\n")
print("Reformatting complete!")

pd.set_option('display.float_format', '{:.6f}'.format)

states = states[statesFileColumns]
print("\n")
print("Concationation complete!")
print(states.head(20).to_string(index=False))
print("\n")
print("Printing states file...")
states = states.to_string(index=False, header=False)
statesFile = "14N-1H3__CoYuTe-Marvelised-2024.states"
with open(statesFile, "w+") as FileToWriteTo:
    FileToWriteTo.write(states)
print("New states file ready!")