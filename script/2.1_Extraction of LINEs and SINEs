# African Hedgehog
## Step 1: Constructing a Non-Redundant Library

# Search and export repeat families for Eulipotyphla, its ancestor nodes, and all its descendant groups:
~/raid/pub_software/repeat/RepeatMasker/famdb.py -i ~/raid/pub_software/repeat/RepeatMasker/Libraries/RepeatMaskerLib.h5 families -f embl -ad Eulipotyphla > Eulipotyphla_ad.embl

# Convert EMBL format to FASTA for later merging with RepeatModeler2 results
~/raid/pub_software/repeat/RepeatMasker/util/buildRMLibFromEMBL.pl Eulipotyphla_ad.embl > Eulipotyphla_ad.fasta

~/raid/pub_software/repeat/RepeatMasker/util/buildRMLibFromEMBL.pl ~/raid/pub_software/repeat/RepeatMasker/Libraries/RMRBSeqs.embl > RMRBSeqs.embl.fasta

# Merge with de novo predicted sequences
cat afa-families.fa Eulipotyphla_ad.fasta RMRBSeqs.embl.fasta > all_final.fasta

# Remove redundancy using cd-hit
cd-hit -i all_final.fasta -o new.fasta -n 5 -d 0 -aL 0.99 -c 0.8 -s 0.8 -T 60 -M 30000
mv new.fasta cdhit.fa

# Convert all bases to uppercase
# seqkit seq -u afa_chr_ngap.fa > afa_chr_ngap1.fa

## Step 2: Identifying Repetitive Sequences

getorf -sequence family.fasta -outseq family.orf -minsize 300
pfam_scan.pl -cpu 100 -fasta family.orf -dir ~/raid/pub_software/repeat/Pfam_db/ > pfam.results
awk '{if ($6~/^PF/) {print $1}}' pfam.results | sed 's/\#/ /1' | awk '{print $1}' > pf.domains.names

cat pf.domains.names | sort | uniq -c | awk '{print $2,$1-1}' | sort -k 2 > pf.domains.count

## Step 3: Genome Alignment and Multiple Sequence Comparison

# Split repeat-masked sequences
sed "s/\// /g" family.fasta > family_namechg.fasta
faSplit byname family_namechg.fasta family_ind_se
mkdir family_ind_blast
cd family_ind_seq

i=1
for j in $(ls *.fa)
do
   kk[i]=$j
   i=$(($i+1))
done

cd ../family_ind_blast
for ((i=1401;i<=1680;i++))
do
nohup bash ~/raid/pub_software/repeat/TE_ManAnnot/bin/make_fasta_from_blast.sh ../afa_chr_ngap1.fa ../family_ind_seq/${kk[$i]} 0 500 &
done

# Remove empty alignment files
find . -type f -size 0 -delete

# Improve efficiency by reducing sequence numbers
mkdir family_ind_msa
cd family_ind_msa
for i in $(ls ../family_ind_blast/)
do
   samtools faidx $i
   j=$(cat $i.fai | wc -l)
   if [ $j -gt 120 ]
   then
       bash ../ready_for_MSA1.sh $i 100 25
   else
       cp $i $i.rdmSubset
   fi
done

# Perform multiple sequence alignment
cd family_ind_msa
i=1
for j in $(ls *.rdmSubset)
do
   kk[i]=$j
   i=$(($i+1))
done

mkdir family_ind_maf
cd family_ind_maf

for ((i=401;i<=633;i++))
do
mafft --reorder --thread -1 ../family_ind_msa/${kk[$i]} > ${kk[$i]}.maf &
done

# Remove families with only one aligned sequence
find . -type f -size 0 -delete

# Correct multiple sequence alignments
mkdir family_ind_coffee
cd family_ind_coffee
for i in $(ls ../family_ind_maf)
do
   t_coffee -other_pg seq_reformat -in ../family_ind_maf/$i -action +rm_gap 80 | head -n -12 > $i.t
done

## Step 4: Obtain Consensus Sequences and Annotate

mkdir family_ind_cons
cd family_ind_cons
for i in $(ls ../family_ind_coffee)
do
   cons -plurality 0.1 -sequence ../family_ind_coffee/$i -outseq $i.cons
done

# Annotate using HMM
mkdir family_ind_hmm
cd family_ind_hmm
for i in $(ls ../family_ind_maf)
do
   hmmbuild --cpu 86 --amino $i.hmm ../family_ind_cons/$i.t.cons
done

for i in $(ls *.hmm)
do
   hmmsearch --cpu 86 $i ../uniprot_sprot.fasta > $i.out
done

# TE-AID analysis
for i in $(ls ../family_ind_cons)
do
   ~/raid/pub_software/repeat/TE-Aid/TE-Aid -q ../family_ind_cons/$i -m 10 -g ../afa_chr_ngap1.fa -o $i.1
done
