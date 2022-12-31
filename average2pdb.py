from Bio.PDB import PDBParser

# Parse the PDB files and extract the coordinates of the atoms
pdb_parser = PDBParser()
structure_1 = pdb_parser.get_structure('structure_1', 'pdb_file_1.pdb')
structure_2 = pdb_parser.get_structure('structure_2', 'pdb_file_2.pdb')

# Create a new PDB file to write the average coordinates to
out_file = open('average_structure.pdb', 'w')

# Iterate over the atoms in both structures and calculate the average coordinates
for model_1, model_2 in zip(structure_1, structure_2):
    for chain_1, chain_2 in zip(model_1, model_2):
        for residue_1, residue_2 in zip(chain_1, chain_2):
            for atom_1, atom_2 in zip(residue_1, residue_2):
                # Calculate the average coordinates
                avg_x = (atom_1.get_coord()[0] + atom_2.get_coord()[0]) / 2
                avg_y = (atom_1.get_coord()[1] + atom_2.get_coord()[1]) / 2
                avg_z = (atom_1.get_coord()[2] + atom_2.get_coord()[2]) / 2
                avg_coord = (avg_x, avg_y, avg_z)
                
                # Write the average coordinates to the output PDB file
                out_file.write('ATOM  {}{}{}{}{}{}{}{}{}\n'.format(
                    atom_1.get_serial(),
                    atom_1.get_name(),
                    atom_1.get_altloc(),
                    atom_1.get_resname(),
                    atom_1.get_chain(),
                    atom_1.get_resid(),
                    atom_1.get_icode(),
                    '{:8.3f}'.format(avg_x),
                    '{:8.3f}'.format(avg_y),
                    '{:8.3f}'.format(avg_z)
                ))
                
# Close the output file
out_file.close()
