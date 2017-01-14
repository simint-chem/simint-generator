for I in Al2Cl6.xyz  Anthracene.xyz  Benzene.xyz  C4H4Cl2F2.xyz  Ethane.xyz  R-camphor.xyz  Water.xyz
do
  for J in 6-31gss.gbs  aug-cc-pvqz.gbs  aug-cc-pvtz.gbs  dzp.gbs  roos-ano-tz.gbs  sto-3g.gbs
  do
    ./gen_mol.py xyz/$I bas/$J ../$I.$J.mol
  done
done
