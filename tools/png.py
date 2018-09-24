from ase.io import read, write

data = read("CONTCAR")
write("test.png", data, rotation="-60z, 120x")
