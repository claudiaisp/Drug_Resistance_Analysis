#!/bin/sh
for a in /"path were the files are"/*
do	
	n1="/"path were the files are"/"
	name1=${a#"$n1"}
	amrfinder -n "$a" -o "/media/claudia/Backup/claudia/Sequences/Finder/${name1}"
	
done

for b in /"path were the files are"/*
do
	n2="/"path were the files are"/"
	name2=${b#"$n2"}
	blastn -query $b -db megares_drugs_database -outfmt "6 qseqid sseqid slen length pident mismatch gapopen qstart qend sstart send evalue bitscore" -perc_identity 95.00 -out /"path where the files will be stored"/${name2}
	blastn -query $b -db megares_modified_database -outfmt "6 qseqid sseqid slen length pident mismatch gapopen qstart qend sstart send evalue bitscore" -perc_identity 95.00 -out /"path where the files will be stored"/${name2}
done


