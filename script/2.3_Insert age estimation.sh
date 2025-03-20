## Insert age estimation 

/home/dell/raid/pub_software/repeat/RepeatMasker/util/buildSummary.pl afa_chr_ngap1.fa.out > summary.out

/home/dell/raid/pub_software/repeat/RepeatMasker/util/calcDivergenceFromAlign.pl -s sesame.divsum afa_chr_ngap1.fa.cat

Perl /home/dell/raid/pub_software/repeat/RepeatMasker/util/createRepeatLandscape.pl -div sesame.divsum -g 2780000000 > sesame.html
