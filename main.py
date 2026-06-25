import difflib
from difflib import SequenceMatcher
import math

radii = {}
def makeRadii():
    radii['H']  = 0.31
    radii['C']  = 0.76
    radii['N']  = 0.71
    radii['O']  = 0.66
    radii['F']  = 0.57
    radii['CL'] = 1.02
    radii['BR'] = 1.20
    radii['I']  = 1.39
    radii['P']  = 1.07
    radii['S']  = 1.05
    radii['SI'] = 1.11
    radii['NA'] = 1.66
    radii['B']  = 0.84
    radii['AL'] = 1.21
    radii['LI'] = 1.28
    radii['K']  = 2.03
    radii['SE'] = 1.20
    radii['SN'] = 1.39
    radii['RB'] = 2.20
    radii['AS'] = 1.19
    radii['TE'] = 1.38
    radii['BA'] = 2.15
    radii['CA'] = 1.76
    radii['D']  = 0.31
    radii['MG'] = 1.41
    radii['BE'] = 0.96
    radii['CS'] = 2.44
    radii['SR'] = 1.95
makeRadii()

#ps the xyz might be zyx because thats what matched the files change when implimenting probably
def parse(file1):
    output = ""
    with open("pdbParser/output/output.pdb", "w") as file:
        file.write(output)

    with open(file1, "r") as f:
        lines = f.readlines()

    loopNumber = 0
    atomNum = 0
    connections = {}
    positions = {}
    elements = {}

    #header
    #HEADER 
    #CRYST1

    #check to see if it ends with the cif file
    if(file1.endswith(".cif")):
        print("Ran CIF file")
        i = 0
        lineIt = 0
        labelIt = symbolIt = xIt = yIt = zIt = None
        while i < len(lines):
            line = lines[i].strip()

            if line.startswith("_cell_length_a"):
                a = float(line.split()[1])
            if line.startswith("_cell_length_b"):
                b = float(line.split()[1])
            if line.startswith("_cell_length_c"):
                c = float(line.split()[1])

            id = {}

            if line == "loop_":
                i += 1
                while i < len(lines) and lines[i].strip().startswith("_atom_site"):
                    if lines[i].strip() == "_atom_site_label":
                        labelIt = lineIt
                    if lines[i].strip() == "_atom_site_type_symbol":
                        symbolIt = lineIt
                    if lines[i].strip() == "_atom_site_fract_y":
                        yIt = lineIt
                    if lines[i].strip() == "_atom_site_fract_x":
                        xIt = lineIt
                    if lines[i].strip() == "_atom_site_fract_z":
                        zIt = lineIt
                    lineIt += 1
                    i += 1
                if None in (labelIt, symbolIt, yIt, xIt, zIt):
                    print("Sorry your file didn't include enough data. Please input another file and try again")
                    return
                while i < len(lines):
                    if(lines[i].strip() == ""): 
                        i += 1
                        break
                    element = lines[i].strip().split()[symbolIt]
                    x = float(lines[i].strip().split()[xIt]) * a
                    y = float(lines[i].strip().split()[yIt]) * b
                    z = float(lines[i].strip().split()[zIt]) * c
                    atomNum += 1
                    connections[atomNum] = []
                    positions[atomNum] = [x, y, z]
                    elements[atomNum] = element
                    output += "ATOM      " + str(atomNum) + "  " + lines[i].strip().split()[0][0] + "   MOL A   " + lines[i].strip().split()[0][1] + "      " + f"{x: .3f}" + "  " + f"{y: .3f}" + "  " + f"{z: .3f}" + "  1.00  0.00           " + element + "\n" 
                    i += 1  
            i += 1
    #check for other file types(mol/pdb/sdf/xyz/cor) res
    if(file1.endswith(".xyz")):
        print("Ran XYZ file")
        i = 0
        while i < len(lines):
            line = lines[i].strip()
            i += 1
            while i < len(lines) and i > 1:
                element = lines[i].strip().split()[0]
                x = float(lines[i].strip().split()[1])
                y = float(lines[i].strip().split()[2])
                z = float(lines[i].strip().split()[3])
                atomNum += 1
                connections[atomNum] = []
                positions[atomNum] = [x, y, z]
                elements[atomNum] = element
                output += "ATOM      " + str(atomNum) + "  " + element + "   MOL A   " + "1" + "      " + f"{x: .3f}" + "  " + f"{y: .3f}" + "  " + f"{z: .3f}" + "  1.00  0.00           " + element + "\n"
                i += 1
        i += 1
    
    if(file1.endswith(".mol") or file1.endswith(".sdf")):
        if(file1.endswith(".mol")): print("Ran MOL file")
        if(file1.endswith(".mol2")): print("Ran MOL2 file")
        if(file1.endswith(".sdf")): print("Ran SDF file")
        i = 0
        atomLen = float(lines[4].strip().split()[0])
        while i < len(lines):
            line = lines[i].strip()
            i += 1
            while i <= atomLen + 4 and i > 4:
                element = lines[i].strip().split()[3]
                x = float(lines[i].strip().split()[0])
                y = float(lines[i].strip().split()[1])
                z = float(lines[i].strip().split()[2])
                atomNum += 1
                connections[atomNum] = []
                positions[atomNum] = [x, y, z]
                elements[atomNum] = element
                output += "ATOM      " + str(atomNum) + "  " + element + "   MOL A   " + "1" + "      " + f"{x: .3f}" + "  " + f"{y: .3f}" + "  " + f"{z: .3f}" + "  1.00  0.00           " + element + "\n"
                i += 1
        i += 1

    if(file1.endswith(".cor")):
        print("Ran COR file")
        i = 0
        while i < len(lines):
            line = lines[i].strip()
            while "*" in lines[i].strip().split():
                i += 1
            i += 1
            while i < len(lines):
                element = lines[i].strip().split()[3]
                x = float(lines[i].strip().split()[4])
                y = float(lines[i].strip().split()[5])
                z = float(lines[i].strip().split()[6])
                atomNum = lines[i].strip().split()[0]
                connections[atomNum] = []
                positions[atomNum] = [x, y, z]
                elements[atomNum] = element
                output += "ATOM      " + str(atomNum) + "  " + element + "   MOL A   " + lines[i].strip().split()[1] + "      " + f"{x: .3f}" + "  " + f"{y: .3f}" + "  " + f"{z: .3f}" + "  1.00  0.00           " + element + "\n"
                i += 1
        i += 1

    #iterate through each atom and compare it to every other atom
    for center in positions:
        #print(f"Center: {center}")
        cx = positions[center][0]
        cy = positions[center][1]
        cz = positions[center][2]
        output += "CONECT    " + str(center)
        for outer in positions:
            if outer != center:
                #print(f"Outer: {outer}")
    #get the distance, get the dir by subtracting x y and z and add them up than square root it
                ox = positions[outer][0]
                oy = positions[outer][1]
                oz = positions[outer][2]
                distanceVector = [ox - cx, oy - cy, oz - cz]
                distance = math.sqrt(distanceVector[0]**2 + distanceVector[1]**2 + distanceVector[2]**2)
    #if the distance of the outer from the center is less than the two radii added together than multiplied by the tolerance factor
                if distance < (radii[elements[center]] + radii[elements[outer]]) * 1.2:
                    output += "    " + str(outer)

        output += "\n"
    output += "END"
    with open("pdbParser/output/output.pdb", "w") as file:
        file.write(output)

def test(file1, file2):
# Calculate the similarity ratio
    with open(file1, "r") as f:
        content1 = f.read()
    with open(file2, "r") as f:
        content2 = f.read()

    ratio = SequenceMatcher(None, content1, content2).ratio()
    print(f"Similarity: {ratio * 100:.2f}%")

    if(ratio == 1.0): return
    diff = difflib.ndiff(content1, content2)
    clean_diff = []

    for item in diff:
        if item.startswith("?"):
            continue
            
        prefix = item[0]
        char = item[2]
        
        if prefix == "-":
            # Red text for deletions
            clean_diff.append(f"\033[91m{char}\033[0m")
        elif prefix == "+":
            # Green text for additions
            clean_diff.append(f"\033[92m{char}\033[0m")
        elif prefix == " ":
            # Normal text for matches
            clean_diff.append(char)

    print("".join(clean_diff))

def inputFun():
    includedFiles = ["mol", "xyz", "sdf", "cor", "cif", "all"]

    inputType = input("Please enter file type: ")
    if(inputType.lower() not in includedFiles):
        print("Sorry that file is not supported.")
        inputType = input("Please enter file type: ")
    elif(inputType.lower() == "all"):
        for file in includedFiles:
            if file != "all":
                parse(f"pdbParser/Data/{file.upper()}/data1.{file.lower()}")
                file1 = "pdbParser/output/output.pdb"
                file2 = "pdbParser/Data/PDB/data1.pdb"
                test(file1, file2)
    else:
        parse(f"pdbParser/Data/{inputType.upper()}/data1.{inputType.lower()}")
        file1 = "pdbParser/output/output.pdb"
        file2 = "pdbParser/Data/PDB/data1.pdb"
        test(file1, file2)
    
inputFun()