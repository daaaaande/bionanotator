# bionanotator - annotate bionano .vcf with annovar and more

### built for use on Linux
### requirements
- annovar installation with all needed hg38 databases:
- place of table_annovar.pl  in ../../ of bionanotator.pl
- hg38 annovar databases in ../../humanhg38:  refGene,cytoBand,dgvMerged,gwasCatalog,wgEncodeRegDnaseClustered,genomicSuperDups,wgRna,gnomad30_genome

# HOWTO
`perl bionanotator.pl --i ../bionano_out_new.vcf --del 0`
### parameters
`--in in.vcf ` (--i); input bionano .vcf file

`--hallmark_file hallmark_genes.tsv` (--h); Cancer hallmark .tsv mapping file  

`--logfile logfile.log` (--l); Logfile for current run

`--sample testsample` (--s); samplename for this run

`--delete_in_between_files 1 ` (--d ); [0/1] delete annovar in-betwen files when annotation is done


### example input .vcf file
```bash
test@test$ head bionano_out_new.vcf
1       207552951       207571434       0       0
2       1512267 1538271 0       0
3       41687791        41694625        0       0
4       99054229        99079600        0       0
6       21814479        21847958        0       0
7       51195881        51195881        0       0
test@test$
```
- note that there is no header, ref and alt columns are both 0
- more input is not needed

### contributing
 Feedback is always welcome!


### roadmap
- add bigwig value column
- change empty column remover from static to flexible based on file 
