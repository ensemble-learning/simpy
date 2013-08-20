"""packmol input template
"""

INP = """#
tolerance 2.0

filetype pdb

output sim.pdb

structure %pdbfile% 
  number %n% 
  inside box 0. 0. 0. %l% %l% %l%
end structure

"""
