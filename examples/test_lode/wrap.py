from ase.io import read,write

c = read("water_100angs.xyz",":")
for i in range(2):
    c[i].wrap()
write("water_100angs.xyz",c)
