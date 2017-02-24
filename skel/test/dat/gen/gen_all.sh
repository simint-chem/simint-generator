for I in Al2Cl6  anthracene  benzene  C4H4Cl2F2  ethane  R-camphor  water
do
  for J in 6-31gss  aug-cc-pvqz  aug-cc-pvtz  dzp  roos-ano-tz  sto-3g
  do
    ./gen_mol.py xyz/$I.xyz bas/$J.gbs ../$I.$J.mol
  done
done
