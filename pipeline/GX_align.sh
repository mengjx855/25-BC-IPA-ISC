# grep ref seq and rename id
seqkit grep -r -f Csp.key_gene.tsv ref_seqs/Csp.GCF_000155085.1.ffn > Csp.key_gene.ffn
seqkit grep -r -f Csp.key_gene.tsv ref_seqs/Csp.GCF_000155085.1.faa > Csp.key_gene.faa

# blast
makeblastdb -in Csp.key_gene.faa -out Csp.key_gene.faa -dbtype prot
blastp -query ref_seqs/Bco.GCF_034355335.1.faa -db Csp.key_gene.faa -evalue 1e-5 -max_target_seqs 5 -num_threads 10 -outfmt 6 -out Bco_faa.btp
blastx -query ref_seqs/Bco.GCF_034355335.1.ffn -db Csp.key_gene.faa -evalue 1e-5 -max_target_seqs 5 -num_threads 10 -outfmt 6 -out Bco_ffn.btp
blastx -query ref_seqs/Bco.GCF_034355335.1.fna -db Csp.key_gene.faa -evalue 1e-5 -max_target_seqs 5 -num_threads 10 -outfmt 6 -out Bco_fna.btp

# fldH gene
seqkit grep -r -p 'fldH' Csp.key_gene.faa | seqkit replace -p "^" -r "C.sporogenes|" > key_fldH.faa
cat Bco_faa.btp | grep fldH | sort -rk 3 | cut -f1 | head -1 | seqkit grep -f - ref_seqs/Bco.GCF_034355335.1.faa | seqkit replace -p " .*" -r "|fldH" | seqkit replace -p "^" -r "B.coccoides|" >> key_fldH.faa
cat Eco_faa.btp | grep fldH | sort -rk 3 | cut -f1 | head -1 | seqkit grep -f - ref_seqs/Eco.GCF_000005845.2.faa | seqkit replace -p " .*" -r "|fldH" | seqkit replace -p "^" -r "E.coli|" >> key_fldH.faa
muscle -align key_fldH.faa -output key_fldH.msa
trimal -in key_fldH.msa -out key_fldH.msa.trim