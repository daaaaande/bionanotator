

[![DOI](https://zenodo.org/badge/257535671.svg)](https://zenodo.org/badge/latestdoi/257535671)

# bionanotator - annotate bionano .csv with annovar and more
### what it does
- takes the raw bionano output and splits it into vcf and tsv file: .vcf is used for annovar, .tsv includes useful columns later pasted into the final outfile
- detects for transversions- if start coordinate is bigger than end, it flips those to preven bedtools errors in further analysis
- annotates with annovar for included databases (hg38, table_annovar)
- scans included genes for cancer hallmark genes + leukemia predisposition genes
- uses bedtools to count the overlapping regions in the dvg_merged database with at least 50% overlap (query vs db match  AND db match vs query)
- if allele frequencies are found, it leaves the data in the output file. if all allele frequencies columns are empty , these will be removed
- summarizes results in one .tsv file ($sample.annotated.finished.tsv)

>>  built for use on Linux

### requirements
- PERL installed : sudo apt install perl -y
- annovar installation with all needed hg38 databases:
- place of table_annovar.pl  in ../../ of bionanotator.pl
- hg38 annovar databases in ../../humanhg38:  refGene,cytoBand,dgvMerged,gwasCatalog,wgEncodeRegDnaseClustered,genomicSuperDups,wgRna
- for the hit counts to work, remove first column from dgv_merged.txt annovar database file to turn it into a valid .bed file: ../../humanhg38/dgv_merged_bedfile.bed
 ( can be done with `cat dgv_merged.txt | cut -f 2- > dgv_merged_bedfile.bed`)

- bedtools intersect installed
# HOWTO
`perl bionanotator.pl --i ../test1_bnout_raw.csv --del 0 --sample test_name --case g `
### parameters
`--in in.csv ` (--i); input bionano .csv file, can be somatic or germline case

`--hallmark_file hallmark_genes.tsv` (--h); Cancer hallmark .tsv mapping file, part of this repo.

`--logfile logfile.log` (--l); Logfile for current run

`--sample testsample` (--s); samplename for this run

`--delete_in_between_files 0 ` (--d ); [0/1] delete annovar in-betwen files when annotation is done

`--case g `(--c ); [g/s] somatic (s) vs germline (g) case, somatic excludes parent info columns and adapts to different input file


### input file
- the input file for bionanotator is a .csv
- septarator is ";"
- header is present
- example for germline case :
```bash
test@test$ head bionanotator_infile.csv
#hSmapEntryID;QryContigID;RefcontigID1;RefcontigID2;QryStartPos;QryEndPos;RefStartPos;RefEndPos;Confidence;Type;XmapID1;XmapID2;LinkID;QryStartIdx;QryEndIdx;RefStartIdx;RefEndIdx;Zygosity;Genotype;GenotypeGroup;RawConfidence;RawConfidenceLeft;RawConfidenceRight;RawConfidenceCenter;SVsize;SVfreq;orientation;Sample;Algorithm;Size;Present_in_%_of_BNG_control_samples;Present_in_%_of_BNG_control_samples_with_the_same_enzyme;Fail_assembly_chimeric_score;OverlapGenes;NearestNonOverlapGene;NearestNonOverlapGeneDistance;PutativeGeneFusion;Found_in_parents_assemblies;Found_in_parents_molecules;Found_in_self_molecules;Mother_molecule_count;Father_molecule_count;Self_molecule_count
104;502;1;1;7557703.6;7580046.6;10180177.0;10204152.0;0.99;deletion;6;6;-1;967;968;1314;1315;heterozygous;1;-1;17. Feb;1240.51;1249.24;17. Feb;1631.9;0.564;NA ;KB0060_C___De_novo;assembly_comparison;1632;0.0;0.0;not_applicable;UBE4B;KIF1B;6554.0;-;mother;mother;yes;24;2;29
404;312;1;1;15708806.2;15728287.7;98931448.0;98953857.0;1.00;deletion;31;31;-1;3203;3205;18688;18691;heterozygous;1;-1;64.45;4036.55;2630.41;64.45;2927.5;0.469;NA ;KB0060_C___De_novo;assembly_comparison;2928;0.5;0.0;not_applicable;LPPR5;JC244945;38532.0;-;mother;mother;yes;25;0;32
505;201;1;1;38542692.1;38558056.6;167788315.9;167805848.0;1.00;deletion;68;68;-1;7658;7659;28268;28271;heterozygous;1;-1;42.86;4261.15;9506.92;42.86;2167.5;0.436;NA ;KB0060_C___De_novo;assembly_comparison;2169;0.5;0.0;not_applicable;MPZL1;ADCY10;3540.0;-;mother;mother;yes;26;1;23
887;91;2;2;64016582.5;64025409.6;64072911.0;64075906.5;1.00;insertion;76;76;-1;11913;11915;12807;12808;heterozygous;1;-1;1748.29;15074.26;5698.48;1748.29;5831.7;0.420;NA ;KB0060_C___De_novo;assembly_comparison;5832;0.0;0.0;not_applicable;DQ600650;PELI1;16745.5;-;mother;mother;yes;32;0;30
947;92;2;2;4229734.6;4234987.2;4178883.5;4186522.0;1.00;deletion;77;77;-1;620;621;643;645;heterozygous;1;-1;184.75;673.29;20095.23;184.75;2385.9;0.567;NA ;KB0060_C___De_novo;assembly_comparison;2386;0.0;0.0;not_applicable;-;JA429818;37109.0;-;mother;mother;yes;38;2;26
1141;472;2;2;8879965.1;8882824.8;122289954.0;122300947.0;1.00;deletion;109;109;-1;1652;1653;23605;23606;heterozygous;1;-1;1479.19;2094.37;2179.48;1479.19;8133.3;0.457;NA ;KB0060_C___De_novo;assembly_comparison;8134;0.0;0.0;not_applicable;DQ589229;DQ571825;203560.0;-;mother;mother;yes;19;0;23

```

### contribute
 Feedback is always welcome!
