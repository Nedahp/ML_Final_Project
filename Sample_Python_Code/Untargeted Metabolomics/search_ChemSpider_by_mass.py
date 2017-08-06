


# search ChemSpider database for compounds whose mass is equal to the input mass in a predefined error margine
# input: a mass measurement, an error margine(such as 10ppm)
# output: a list of compounds in ChemSpider whose masses fall into the interval (mass-margine,mass+margine)
def search_by_mass(mass,margine):
    # pip install chemspipy
    from chemspipy import ChemSpider
    # register to generate a sequrity code
    cs = ChemSpider('dfdc677d-e7d3-435b-a74e-bfe6167a3899')
    # search the data base
    CSIDs = []
    for result in cs.simple_search_by_mass(mass,margine):
        CSIDs.append(result.csid)
    return CSIDs



# read an excel file to extract list of mass measurements
# input: a file name
# output: a list for mass measurements
def readFile(fileName):
    import pandas as pd
    from pandas import Series, DataFrame
    productMasses = pd.read_excel(fileName)
    mass_list = productMasses['ProductMass']
    return mass_list



# for each product, find a list of compound matches in ChemSpider based on mass
# input: a list for masses, an error margine(such as 10ppm)
# output: a dictionary with each mass as the key and a list of matched compounds as the value
def extract_compounds_with_matchedMasses_from_ChemSpider(productMasses,massMargine):
    matched_in_ChemSpider = {}
    for i in range(len(productMasses)):
        print i
        compounds = []
        mass = productMasses[i]
        compounds = search_by_mass(mass,massMargine)
        if compounds != []:
            matched_in_ChemSpider.update({i:compounds})
    return matched_in_ChemSpider



# retrive a molecule's mol file from hard drive
# input: an excel file with information of molecules, index of the query molecule
# output: the query molecul's mol file
def read_product_molFile(massFile_Name, index):
    import pandas as pd
    from pandas import Series, DataFrame
    path_to_product_MolFile_folder = '/Users/Neda/Desktop/Operators_Project/Prof_Lee_project/project2+cho/code/'
    productMasses = pd.read_excel(massFile_Name)
    parentMetkeggID = productMasses['parentMetkeggID'][index]
    parentMetkeggID = str(parentMetkeggID)
    productID = productMasses['productID'][index]
    productID = str(productID)
    local_path_to_molFile = parentMetkeggID +'/product_' + productID + '.mol'
    global_path_to_molFile= path_to_product_MolFile_folder + local_path_to_molFile
    f = open(global_path_to_molFile,'r')
    product_molFile = f.read()
    return product_molFile



# extract number of atoms and bonds for a molecule from its molfile
# input: a mol file
# output: number of atoms, number of bonds
def extract_bondNum_atomNum(molFile):
    lines = molFile.splitlines()
    count_nonEmpty_lines = 0;
    for line in lines:
        if (line != ' ') and (line != ''):
            a = line.split() 
            b = filter(None, a)
            if b[len(b)-1] == 'V2000':
                count_nonEmpty_lines = count_nonEmpty_lines + 1
                if count_nonEmpty_lines == 1:
                    first_line = line
    first_line = first_line.split(' ')
    first_line = filter(None, first_line)
    atom_num = int(first_line[0])
    bond_num = int(first_line[1])
    return atom_num,bond_num



# extract atoms types and their indexes from a mol file
# input: number of atoms, number of bonds, mol file
# output: a dictionary(atomsTable) with index of an atom as the key and the atom type as the value
def extract_atomTable(atomNum,bondNum,molFile):
    count_nonEmpty_lines = -1;
    atomsTable = {}
    firstLine_is_seen = False
    lines = molFile.splitlines()
    for line in lines:
        if (line != ' ') and (line != ''):
            a = line.split() 
            b = filter(None, a)
            if b[len(b)-1] == 'V2000':
                firstLine_is_seen = True
            if firstLine_is_seen:
                count_nonEmpty_lines = count_nonEmpty_lines + 1
            if (count_nonEmpty_lines > 0) and (count_nonEmpty_lines <= atomNum):
                split_line = line.split(' ')
                split_line = filter(None, split_line)
                atomType = split_line[3]
                atomsTable.update({count_nonEmpty_lines:atomType})
    return atomsTable



# extract the connection of atoms and the number of bonds between each pair from a mol file
# input: a dictionary with atom types(atomsTable), number of atoms, number of bonds, mol file
# output: a dictionary(atoms_bonds_Table) with an index as the key and a list of [atom1, atom2, num_bind] as the value
def extract_atoms_bonds_table(atomTable,atomNum,bondNum,molFile):
    count_nonEmpty_lines = -1 * (1 + atomNum);
    atoms_bonds_Table = {}
    firstLine_is_seen = False
    lines = molFile.splitlines()
    for line in lines:
        if (line != ' ') and (line != ''):
            a = line.split() 
            b = filter(None, a)
            if b[len(b)-1] == 'V2000':
                firstLine_is_seen = True
            if firstLine_is_seen:
                count_nonEmpty_lines = count_nonEmpty_lines + 1
            if (count_nonEmpty_lines > 0) and (count_nonEmpty_lines <= bondNum):
                split_line = line.split(' ')
                split_line = filter(None, split_line)
                atom1 = atomTable[int(split_line[0])]
                atom2 = atomTable[int(split_line[1])]
                bonds = int(split_line[2])
                atoms_bonds_Table.update({count_nonEmpty_lines:[atom1,atom2,bonds]})
    return atoms_bonds_Table



# compare two molecule based on number of atoms, bonds and connections, and return if they are the same
# input: two atoms_bonds_Tables
# output: a boolean
def compare_atoms_bonds_Tables(atoms_bonds_Table1,atoms_bonds_Table2):
    for i in range(len(atoms_bonds_Table2)):
        atom1_1 = atoms_bonds_Table1[i+1][0]
        atom1_2 = atoms_bonds_Table1[i+1][1]
        bond1 = atoms_bonds_Table1[i+1][2]
    
        atom2_1 = atoms_bonds_Table2[i+1][0]
        atom2_2 = atoms_bonds_Table2[i+1][1]
        bond2 = atoms_bonds_Table2[i+1][2]
    
        if (((atom1_1 == atom2_1) and (atom1_2 == atom2_2)) or (atom1_1 == atom2_2) and (atom1_2 == atom2_1)):
            if bond1 == bond2:
                is_the_same = True
            else:
                is_the_same = False
        else:
            is_the_same = False
    return is_the_same



# compare two mol files and return if they are the same
# input: two mol files
# output: a boolean
def compare_two_molFiles(mol_1,mol_2):
    atomNum1,bondNum1 = extract_bondNum_atomNum(mol_1)
    atomTable1 = extract_atomTable(atomNum1,bondNum1,mol_1)
    atoms_bonds_Table1 = extract_atoms_bonds_table(atomTable1,atomNum1,bondNum1,mol_1)

    atomNum2,bondNum2 = extract_bondNum_atomNum(mol_2)
    atomTable2 = extract_atomTable(atomNum2,bondNum2,mol_2)
    atoms_bonds_Table2 = extract_atoms_bonds_table(atomTable2,atomNum2,bondNum2,mol_2)

    if (atomNum1 == atomNum2) and (bondNum1 == bondNum2):
        is_the_same = compare_atoms_bonds_Tables(atoms_bonds_Table1,atoms_bonds_Table2)
    else:
        is_the_same = False
    return is_the_same



# compare a molecule with the compounds matched its mass from ChemSpider database, use mol files for the comparision
# input: a list of matched ChemSpider compunds, information of the input molecule
# output: a dictionary: key = input molecule, value = list of mathed compounds based on mass and mol file 
def find_matches(matched_in_ChemSpider,massFile_Name):
    from chemspipy import ChemSpider
    cs = ChemSpider('dfdc677d-e7d3-435b-a74e-bfe6167a3899')
    for i in matched_in_ChemSpider.keys():
        print i
        # intialiaztion
        matched_compounds = []
        matches = {}
        # load mol file info of the product
        product_molFile = read_product_molFile(massFile_Name, i)
        # for each compound in data base with almost the same mass
        for CSID in matched_in_ChemSpider[i]:
            # extract the compound's mol file
            c = cs.get_compound(CSID)
            ChemSpider_compound_mol_info = c.mol_2d
            # compare the product's and compound's mol files
            is_the_same = compare_two_molFiles(product_molFile,ChemSpider_compound_mol_info)
            # add the compound to the list if it's molfile is the same as the product's
            if is_the_same:
                matched_compounds.append(CSID)
        # if at least one compound found as a match
        if matched_compounds != []:
            matches.update({i:matched_compounds})
    # return the whole matches for products
    return matches



# main function
def mainFunction(massFile_Name,massMargine):
    # read the products'masses
    productMasses = readFile(massFile_Name)
    # for each product mass find compounds in chemspider with the same mass
    matched_in_ChemSpider_by_mass = extract_compounds_with_matchedMasses_from_ChemSpider(productMasses,massMargine)
    # extract the mol information for found chemspider compounds
    # compare the mol information to the product's mol file
    matched_in_ChemSpider_by_molFile = find_matches(matched_in_ChemSpider_by_mass,massFile_Name)
    return matched_in_ChemSpider_by_molFile 



# call main on HilNeg dataset mass measurments with error a predefined margine 
HilNeg_result = mainFunction('HilNeg_productMass_noCHO_noKEGG.xlsx',0.0001)

