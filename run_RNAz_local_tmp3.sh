#! /bin/bash
rm -rf RNAz_out_temp
mkdir RNAz_out_temp
for f in temp/*.aln; do
    b=`basename $f .aln`
    echo Processing alignment $b
    rnazWindow.pl --window=120 --slide=40 --min-seqs=4 temp/$b.aln | RNAz --both-strands -p 0.9 | rnazCluster.pl --html  --html-dir ${b} > RNAz_out_temp/${b}.Window.rnaz.dat    
done
