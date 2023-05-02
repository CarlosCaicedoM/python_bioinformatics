#! /bin/bash
rm -rf results_sRNA/RNAz_out
mkdir results_sRNA/RNAz_out
for f in results_sRNA/IGRs_alignments/*.aln; do
    b=`basename $f .aln`
    echo Processing alignment $b
    rnazSelectSeqs.pl -i 70 -n 4 results_sRNA/IGRs_alignments/$b.aln | RNAz --both-strands -p 0.9 > results_sRNA/RNAz_out/${b}.out
done


