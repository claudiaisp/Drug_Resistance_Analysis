#!/bin/sh
#for a in /media/claudia/Backup/claudia/Sequences/Assemblies/*
#do	
	#n1="/media/claudia/Backup/claudia/Sequences/Assemblies/"
	#name1=${a#"$n1"}
	#amrfinder -n "$a" -o "/media/claudia/Backup/claudia/Sequences/Finder/${name1}"
	
#done

for b in /media/claudia/Backup/claudia/Sequences/Assemblies/*
do
	n2="/media/claudia/Backup/claudia/Sequences/Assemblies/"
	name2=${b#"$n2"}
	blastn -query $b -db megares_drugs_database -outfmt "6 qseqid sseqid slen length pident mismatch gapopen qstart qend sstart send evalue bitscore" -perc_identity 95.00 -out /media/claudia/Backup/claudia/Sequences/Megares/Drugs/${name2}
	blastn -query $b -db megares_modified_database -outfmt "6 qseqid sseqid slen length pident mismatch gapopen qstart qend sstart send evalue bitscore" -perc_identity 95.00 -out /media/claudia/Backup/claudia/Sequences/Megares/Modified/${name2}
done


