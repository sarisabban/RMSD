import sys
from pymol import cmd

cmd.load(sys.argv[1], 'F1')
cmd.load(sys.argv[2], 'F2')
RMSD = round(cmd.align('F1','F2')[0], 3)
print(RMSD)




def RMSD(structure_1, structure_2):
	'''
	Calculate the RMSD between two protein structures using Biopython
	The Biopython algorithm is poorly designed and only aligns local motifs
	rather than full protein structures/complexes.
	'''
	import Bio.PDB
	builder = Bio.PDB.Polypeptide.PPBuilder()
	STR1 = builder.build_peptides(Bio.PDB.PDBParser(QUIET=True)\
		.get_structure('Structure 1', structure_1), aa_only=True)
	STR2 = builder.build_peptides(Bio.PDB.PDBParser(QUIET=True)\
		.get_structure('Structure 2', structure_2), aa_only=True)
	fixed  = [atom for poly in STR1 for atom in poly]
	moving = [atom for poly in STR2 for atom in poly]
	lengths = [len(fixed), len(moving)]
	smallest = min(lengths)
	sup = Bio.PDB.Superimposer()
	sup.set_atoms(fixed[:smallest], moving[:smallest])
	sup.apply(Bio.PDB.PDBParser(QUIET=True)\
		.get_structure('Structure 2', structure_2)[0].get_atoms())
	RMSD = round(sup.rms, 4)
	print(RMSD)
