from ase.io.cube import read_cube_data
data, atoms = read_cube_data('beta_2_potential.cube')
print(data[1][1][1])
print(atoms)
