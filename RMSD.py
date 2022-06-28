import Bio.PDB

def RMSD (structure_1, structure_2):
	''' Calculate the RMSD between two protein structures '''
	type1 = Bio.PDB.Polypeptide.PPBuilder().build_peptides(Bio.PDB.PDBParser(QUIET=True).get_structure('X', structure_1), aa_only=True)
	length1 = type1[-1][-1].get_full_id()[3][1]
	fixed = [atom['CA'] for atom in type1[0]]
	#Second structure
	type2 = Bio.PDB.Polypeptide.PPBuilder().build_peptides(Bio.PDB.PDBParser(QUIET=True).get_structure('X', structure_2), aa_only=True)
	length2 = type2[-1][-1].get_full_id()[3][1]
	moving = [atom['CA'] for atom in type2[0]]
	#Choose the length of the smallest structure
	lengths = [length1, length2]
	smallest = min(int(item) for item in lengths)
	#Find RMSD
	sup = Bio.PDB.Superimposer()
	sup.set_atoms(fixed[:smallest], moving[:smallest])
	sup.apply(Bio.PDB.PDBParser(QUIET=True).get_structure('X', structure_2)[0].get_atoms())
	RMSD = round(sup.rms, 4)
	print(RMSD)

RMSD('6ooy.pdb', '6ooy.pdb')
