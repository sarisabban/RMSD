import Bio.PDB

def RMSD(structure_1, structure_2):
	''' Calculate the RMSD between two protein structures '''
	builder = Bio.PDB.Polypeptide.PPBuilder()
	STR1 = builder.build_peptides(Bio.PDB.PDBParser(QUIET=True)\
		.get_structure('Structure 1', structure_1), aa_only=True)
	STR2 = builder.build_peptides(Bio.PDB.PDBParser(QUIET=True)\
		.get_structure('Structure 2', structure_2), aa_only=True)
	fixed  = [atom['CA'] for poly in STR1 for atom in poly]
	moving = [atom['CA'] for poly in STR2 for atom in poly]
	lengths = [len(fixed), len(moving)]
	print(lengths)
	smallest = min(lengths)
	sup = Bio.PDB.Superimposer()
	sup.set_atoms(fixed[:smallest], moving[:smallest])
	sup.apply(Bio.PDB.PDBParser(QUIET=True)\
		.get_structure('Structure 2', structure_2)[0].get_atoms())
	RMSD = round(sup.rms, 4)
	print(RMSD)

RMSD('6ooy.pdb', '6ooz.pdb')
