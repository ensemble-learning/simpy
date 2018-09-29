o = open("layer.inp", 'w')
l_slab = 20
l_tail = 3
l_body = 16
l_head = 3
l_water = 60

#infor

water = [0,0,0]
water[0] = 2000
water[1] = l_slab + l_tail + l_body + l_head
water[2] = water[1] + l_water

mol1 = [0,0,0,0,0,0,0]
mol1[0] = 48
mol1[1] = l_slab
mol1[2] = l_slab + l_tail + l_body + l_head
mol1[3] = 15
mol1[4] = l_slab + l_tail + l_body
mol1[5] = 3
mol1[6] = l_slab + l_tail

mol2 = [0,0,0,0,0,0,0]
mol2[0] = 48
mol2[1] = l_slab + l_tail + l_body + l_head + l_water
mol2[2] = l_slab + l_tail + l_body + l_head + l_water + l_head + l_body + l_tail
mol2[3] = 15
mol2[4] = l_slab + l_tail + l_body + l_head + l_water + l_head 
mol2[5] = 3
mol2[6] = l_slab + l_tail + l_body + l_head + l_water + l_head + l_body  

o.write("""#
# build a interface with water in the middle of the box
#

tolerance 2.0 
filetype pdb
output interface.pdb

structure water.pdb
  number %d 
  inside box 0. 0. %d. 32. 32. %d.
end structure

structure e1.pdb 
  number %d 
  inside box 0. 0. %d. 32. 32. %d.
  atoms %d 
    over plane 0. 0. 1. %d.
  end atoms
  atoms %d 
    below plane 0. 0. 1. %d.
  end atoms
end structure 

structure e1.pdb 
  number %d
  inside box 0. 0. %d. 32. 32. %d.
  atoms %d
    below plane 0. 0. 1. %d
  end atoms
  atoms %d
    over plane 0. 0. 1. %d.
  end atoms
end structure

"""%( water[0], water[1], water[2], mol1[0],  mol1[1], mol1[2], mol1[3], mol1[4], mol1[5], mol1[6],mol2[0],  mol2[1], mol2[2], mol2[3], mol2[4], mol2[5], mol2[6]))
