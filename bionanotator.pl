#/usr/bin/perl -w
use strict;
use Getopt::Long qw(GetOptions);

# defaults for file parameters
my$hallmark_mapping_file="hallmark_genes.tsv"; # unusual mapping file, not one gene per line
my$infile="in.vcf";
my$logfile="logfile.log";
my$sample="testsample";
# true/false params
my$delete_empty_cols=1;
my$delete_in_between_files=1;

GetOptions('in|i=s'=>\$infile,'sample|s=s'=>\$sample,'hallmark_file|h=s'=>\$hallmark_mapping_file,'delete_empty_cols|d=i'=>\$delete_empty_cols,'logfile|log=s'=>\$logfile,'delete_in_between_files|del=i'=>\$delete_in_between_files)|| warn "Using default parameters:\nin=in.vcf\ndelete_empty_cols=1\nhallmark_file=hallmark_genes.tsv!\nlogfile=logfile.log\ndelete_in_between_files=1\n";
chomp($infile,$hallmark_mapping_file,$sample,$delete_empty_cols,$delete_in_between_files,$logfile);

open(ER,'>>',$logfile)||die "$!\n logfile not found. use --log path/to/your_logfile.log\n";

# hallmark mapping hash filling
print ER "reading hallmark mapping file $hallmark_mapping_file...\n";
open(MA,$hallmark_mapping_file) || die "$!\n hallmark mapping file not found\n";
my@allemappings= <MA>;
my%mapping_hash=();  # mapping hash: gene name is key, hallmark mechanism is value
foreach my $mapline (@allemappings){
  chomp $mapline;
  my@mappingline_parts=split(/\s+/,$mapline);
  my$hallmarktype_full=shift @mappingline_parts;# getting hallmark description properly with some regex cleaning by shifting first array index
  $hallmarktype_full=~s/HALLMARK//;# cleanup of the hallmark name
  $hallmarktype_full=~s/\s+//;
  $hallmarktype_full=~s/_//;
  shift @mappingline_parts; # remove web adress
  # rest of the line is all hallmark genes, need to be cleaned up before saving
  foreach my $hallmg (@mappingline_parts){
    $hallmg=~s/\s+//;
    if($hallmg=~/[A-Z]/){  # checking for empty lines and gene names
      $mapping_hash{"$hallmg"}="$hallmarktype_full";
    }
  }
}
my@allehallmarkg=keys %mapping_hash;
my@all_hm_genes= values %mapping_hash;
close MA;


# check for infile
my$infilecheck=`ls -1f $infile`;
if($infilecheck=~/$infile/gix){
  print ER "found infile. creating rundir...\n";
}
else{
  die "could not find infile\ntry --in yourfile.vcf\n";
}

# prepare the annovar run
mkdir "run_$sample";
print ER "created run dir run_$sample\n";

chdir "run_$sample";
print ER "now in run_$sample/.\n";



# find infile
# create rundir
# go to rundir
# execute annovar - maybe place and DB place as param?
my$annovar_ex=`perl ../../table_annovar.pl  ../$infile ../../humanhg38/ -buildver hg38 -out $sample.annotated -polish -protocol refGene,cytoBand,dgvMerged,gwasCatalog,wgEncodeRegDnaseClustered,genomicSuperDups,wgRna,gnomad30_genome -operation gx,r,r,r,r,r,r,f -nastring . >$logfile`;
print ER "executed annovar: $annovar_ex\n";


# delete all non-csv files in run dir
if($delete_in_between_files){
  my$del_between_files=`find . -type f -not -name '*.txt' -print0 | xargs -0 rm --`;
  print ER "deleted all in-between files per users request: $del_between_files\n";
}


my$outfile_found=`ls -1f *.txt`;
chomp $outfile_found;
print ER "found annovar outfile: $outfile_found\n";


# the final outfile will be : run_$sample/$sample_annotated.tsv
open(FIN,">$sample.annotated.tsv") || die "$!\ncould not write the final outfile\n";
open(OU,$outfile_found)|| die "$!\ncould not open annovar outfile in run_$sample\n";

my@all_out_lines=<OU>;

my$line_count=0;
my$new_header_line="Chr	Start	End	Ref	Alt	Func.refGene\tGene.refGene\tCancer.hallmark	GeneDetail.refGene	cytoBand	dgvMerged	gwasCatalog	wgEncodeRegDnaseClustered	genomicSuperDups	wgRna\n";
print FIN "$new_header_line"; # print the header into final output file
my$hallm="none";
foreach my $outline(@all_out_lines){
  $line_count=$line_count +1;
  chomp $outline;
  if($line_count == 1){# header
  # save header, edit as needed

  }
  else{
    my@all_line_parts=split(/\t/,$outline);
    # find cancer hallmarks in gene names

    my@genes_all=split(/\;/,$all_line_parts[6]);# find each gene
    $hallm="none";
    foreach my $single_gene(@genes_all){
      # check hash for fitting hits
      $single_gene=~s/\"//g;
      if(grep(/$single_gene?/,@allehallmarkg)){
        if($mapping_hash{$single_gene}=~/[A-Z]/){
        $hallm=$mapping_hash{$single_gene};# one hallmark max per SV for now
        }
      }
    }

    # all other file lines
    # split by column
    my$new_line="$all_line_parts[0]\t$all_line_parts[1]\t$all_line_parts[2]\t$all_line_parts[3]\t$all_line_parts[4]\t$all_line_parts[5]\t$all_line_parts[6]\t$hallm\t$all_line_parts[7]\t$all_line_parts[8]\t$all_line_parts[11]\t$all_line_parts[12]\t$all_line_parts[13]\t$all_line_parts[14]\t$all_line_parts[15]\n";
    # cleanup
    $new_line=~s/\t\s+/\t/g;
    print FIN $new_line;

  }
}
#
print "done creating $sample.annotated.tsv.\n";
1;
