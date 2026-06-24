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

def parse(file1):
    with open(file1, "r") as f:
        lines = f.readlines()

    i = 0
    loopNumber = 0
    atomNum = 0
    output = ""
    connections = {}
    positions = {}
    elements = {}
    lineIt = 0

    #header
    #HEADER 
    #CRYST1
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

        if line == "loop_" and loopNumber == 0:
            i += 1
            while i < len(lines) and lines[i].strip().startswith("_atom_site"):
                if lines[i].strip() == "_atom_site_label":
                    labelIt = lineIt
                if lines[i].strip() == "_atom_site_type_symbol":
                    symbolIt = lineIt
                if lines[i].strip() == "_atom_site_fract_x":
                    xIt = lineIt
                if lines[i].strip() == "_atom_site_fract_y":
                    yIt = lineIt
                if lines[i].strip() == "_atom_site_fract_z":
                    zIt = lineIt
                lineIt += 1
                i += 1
            if None in (labelIt, symbolIt, xIt, yIt, zIt):
                return
            while i < len(lines) and lines[i].strip() != "" and lines[i].strip() != "loop_":

                element = lines[i].strip().split()[symbolIt]
                x = float(lines[i].strip().split()[xIt]) * a
                y = float(lines[i].strip().split()[yIt]) * b
                z = float(lines[i].strip().split()[zIt]) * c
                name = lines[i].strip().split()[labelIt]
                atomNum += 1
                id[name] = atomNum
                connections[atomNum] = []
                positions[atomNum] = [x, y, z]
                elements[atomNum] = element
                output += "ATOM      " + str(atomNum) + "  " + lines[i].strip().split()[0][0] + "   MOL A   " + lines[i].strip().split()[0][1] + "      " + f"{y: .3f}" + "  " + f"{x: .3f}" + "  " + f"{z: .3f}" + "  1.00  0.00           " + element + "\n" 
                i += 1

            loopNumber += 1
                    
                    
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
    with open("./output/output.pdb", "w") as file:
        file.write(output)

def test(file1, file2):
# Calculate the similarity ratio

    
    with open(file1, "r") as f:
        content1 = f.read()
    with open(file2, "r") as f:
        content2 = f.read()

    ratio = SequenceMatcher(None, content1, content2).ratio()
    print(f"Similarity: {ratio * 100:.2f}%")

parse("./Data/CIF/data1.cif")
file1 = "./output/output.pdb"
file2 = "./Data/PDB/data1.pdb"

test(file1, file2)