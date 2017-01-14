#!/usr/bin/env python3

#############################
# A quick script to generate
# some molecule files used
# in testing. This make the
# C++ parsing easier
#############################

import argparse

parser = argparse.ArgumentParser()
parser.add_argument("xyz", type=str, help="Path to xyz file")
parser.add_argument("bas", type=str, help="Path to basis file")
parser.add_argument("out", type=str, help="Path to output file")
args = parser.parse_args()

ang_to_bohr = 1.0/0.52917721067


print("");
print("-----------------------------------")
print("- OPTIONS                         -")
print("-----------------------------------")
print("XYZ file: {}".format(args.xyz))
print("BAS file: {}".format(args.bas))
print("Out file: {}".format(args.out))

# Read geometry / XYZ
geo = [ l.strip() for l in open(args.xyz, 'r').readlines()[2:] ]
geo = [ l.split() for l in geo if l ]

# Read basis set
flines = [ l.strip() for l in open(args.bas, 'r').readlines() ]
flines = [ l for l in flines if l and not l.startswith('!') ]
 
atombas = { }

i = 1  # skip initial ****
while i < len(flines):
  atom = { }
  atomsym = flines[i].split()[0]
  atom['Shells'] = [] 

  # Nprim2 counts the general contractions as multiple primitives
  nprim1 = 0
  nprim2 = 0
 
  i += 1
  while flines[i] != '****':
    lsplt = flines[i].split()
    styp = lsplt[0].upper()
    nshell = int(lsplt[1])

    shell = {}
    shell['Type'] = styp


    alpha = []
    coef = []
    i += 1
    for j in range(0, nshell):
      lsplt = flines[i].split()
      alpha.append(lsplt[0])
      coef.append(lsplt[1:])
      i += 1

    shell['Alpha'] = alpha
    shell['Coef'] = coef  

    # add shell to atom 
    atom['Shells'].append(shell)

    nprim1 += len(shell['Alpha'])
    nprim2 += len(shell) * len(shell['Coef'][0])

  atom['NPRIM1'] = nprim1
  atom['NPRIM2'] = nprim2

  # add atom to basis set
  atombas[atomsym] = atom

  i += 1


# Apply basis set to molecule
mol = []

for a in geo:
  mol.append({
               'Sym' : a[0],
               'XYZ' : [ float(a[1])*ang_to_bohr, float(a[2])*ang_to_bohr, float(a[3])*ang_to_bohr ],
               'BAS' : atombas[a[0]]
             })
             

print("-----------------------------------")

# Create file
with open(args.out, 'w') as f:
  f.write(str(len(mol)) + "\n")
  for a in mol:
    f.write("{} {} {} {}\n".format(a['Sym'], len(a['BAS']['Shells']), a['BAS']['NPRIM1'], a['BAS']['NPRIM2']))
    f.write("{} {} {}\n".format(a['XYZ'][0], a['XYZ'][1], a['XYZ'][2]))

    for b in a['BAS']['Shells']:
      nprim = len(b['Alpha'])
      ngen = len(b['Coef'][0])
      f.write('{} {} {}\n'.format(b['Type'], nprim, ngen))

      for i in range(0, nprim):
        f.write('{}'.format(b['Alpha'][i]))

        for j in range(0, ngen):
          f.write('    {}'.format(b['Coef'][i][j]))

        f.write("\n")

