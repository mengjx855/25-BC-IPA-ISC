# Figure 1
#download raw_fq
cat sample_name | parallel -j 4 run_dwn_sra.sh {} raw_fq

#qc
cat sample_name | parallel -j 4 run_fastp_rmhost.sh raw_fq/{}_1.fastq.gz human clean_fq/{}

#mapping
seqkit replace -p "\s.*" -r "" genomes/GCF_900461125.1_49699_D02_genomic.fna > Blautia_coccoides.fna
minimap2 -d Blautia_coccoides.mmi Blautia_coccoides.fna

for i in HeQ_2017.PRJEB15371 YanQ_2023c.PRJEB67456 SchirmerM_2018.PRJNA389280 WengY_2019.PRJNA429990; do
    cat $i.clean_fq.filepath | parallel -j 10 --colsep="\t" run_minimap2.sh {2} Blautia_coccoides.mmi $i/{1}
    ls $i/*bam | parallel -j 10 --plus samtools coverage -d 0 -o {..}.cvm {}
    ls $i/*cvm | parallel -j 20 -q perl -e 'open I, "$ARGV[0]"; $x=0;while(<I>){chomp; next if /numreads/; @s=split/\s+/; $x+=$s[3]}; print "$ARGV[0]\t$x\n"' {} | sed 's/.cvm//g' | csvtk join --left-join -t -f 1 - /share/data1/mjx/meta/20240803_IBD_pub_data_metagenome/SchirmerM_2018.PRJNA389280/clean_fq.lib_size | csvtk add-header -t -n "sample,rc,libs" > $i.profile;
done

# Figure 5
seqkit replace -p '.*' -r 'Bcoccoides.{nr}' GCF_034355335.1_ASM3435533v1_genomic.fna -o Bcoccoides.fa
seqkit replace -p '.*' -r 'Csporogenes.{nr}' GCF_000155085.1_ASM15508v1_genomic.fna -o Csporogenes.fa
seqkit replace -p '.*' -r 'Prussellii.{nr}' GCF_003012055.1_ASM301205v1_genomic.fna -o Prussellii.fa

minimap2 -d ref_genome.mmi ref_genome.fa # combine three genome
cat clean_fq.filelist | parallel -j 8 --colsep="\t" run_minimap2.sh {5} ref_genome.mmi minimap2/{4}
ls *sort.bam | parallel -j 20 --plus samtools coverage -o {/..}.cvg {}
combine_file_zy_folder_allsample.py -D minimap2/ -suffix .cvg -t 1 -n 1 -v 4 -o minimap2.rc