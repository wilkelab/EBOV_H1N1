#!/bin/sh
printf "type\tfile\t\tnucl_ident_pi\n" > align_stats.txt

echo "l2014"
printf "EBOV_2014\tl2014\t" >> align_stats.txt
../scripts/calc_pi.py EBOV_2014/l2014.aln >> align_stats.txt
echo "np2014"
printf "EBOV_2014\tnp2014\t" >> align_stats.txt
../scripts/calc_pi.py EBOV_2014/np2014.aln >> align_stats.txt

echo "l_all_outbreaks"
printf "EBOV_all_outbreaks\tl_nuc\t" >> align_stats.txt
../scripts/calc_pi.py EBOV_all_outbreaks/l_nuc.aln >> align_stats.txt
echo "np_all_outbreaks"
printf "EBOV_all_outbreaks\tnp_nuc\t" >> align_stats.txt
../scripts/calc_pi.py EBOV_all_outbreaks/np_nuc.aln >> align_stats.txt
echo "genome_all_outbreaks"
printf "EBOV_all_outbreaks\tfull_nuc\t" >> align_stats.txt
../scripts/calc_pi.py EBOV_all_outbreaks/full_nuc.aln >> align_stats.txt

echo "HA09_April_noDupStrains"
printf "H1N1\tHA09_April_noDupStrains\t" >> align_stats.txt
../scripts/calc_pi.py H1N1/HA09_April_noDupStrains.aln >> align_stats.txt
echo "NP09_April_noDupStrains"
printf "H1N1\tNP09_April_noDupStrains\t" >> align_stats.txt
../scripts/calc_pi.py H1N1/NP09_April_noDupStrains.aln >> align_stats.txt



