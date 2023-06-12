
USALIGN_HOME=~/Apps/USalign/

$USALIGN_HOME/USalign ./3cna_rcsb.pdb ./2pel_rcsb.pdb -ter 2 -mm 0 -outfmt 0 > 3cnaA_to_2pelA_USalign_SQ_noTransRotMat.log

$USALIGN_HOME/USalign ./3cna_rcsb.pdb ./2pel_rcsb.pdb -ter 2 -mm 0 -outfmt 0 -m out.dat > 3cnaA_to_2pelA_USalign_SQ.log
cat out.dat >> 3cnaA_to_2pelA_USalign_SQ.log
rm out.dat

$USALIGN_HOME/USalign ./3cna_rcsb.pdb ./2pel_rcsb.pdb -ter 2 -mm 3 -outfmt 0 -m out.dat > 3cnaA_to_2pelA_USalign_CP.log
cat out.dat >> 3cnaA_to_2pelA_USalign_CP.log
rm out.dat

$USALIGN_HOME/USalign ./3cna_rcsb.pdb ./2pel_rcsb.pdb -ter 2 -mm 5 -outfmt 0 -m out.dat > 3cnaA_to_2pelA_USalign_fNS.log
cat out.dat >> 3cnaA_to_2pelA_USalign_fNS.log
rm out.dat

$USALIGN_HOME/USalign ./3cna_rcsb.pdb ./2pel_rcsb.pdb -ter 2 -mm 6 -outfmt 0 -m out.dat > 3cnaA_to_2pelA_USalign_sNS.log
cat out.dat >> 3cnaA_to_2pelA_USalign_sNS.log
rm out.dat

