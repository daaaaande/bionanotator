#/usr/bin/perl -w
use strict;
use Getopt::Long qw(GetOptions);
# get the .csv file raw output from bionano
# extract all information for bionanotator .vcf file
# create .vcf file for bionanotator
# extract other useful columns and prepare them for later re-attachment to final output of bionanotator

my$infile="in.csv";# the bionano out
my$vcf_file="vcf.out"; # the vcf with the SV in it
my$tsv_file="tsv.out"; # all extra info tsv for the final output merge
my$case="g";# germline, with mother/father information= "g", somatic, without father/mother information="s"
# TODO: add case and parameter for somatic vs germline- the input columns are different in that case, maybe only for the later pasted into file.
my$include="";# columns to extra add to the output , list of numbers splitted by a ,


GetOptions('in|i=s'=>\$infile,'vcf|v=s'=>\$vcf_file,'tsv|t=s'=>\$tsv_file,'case|c=s'=>\$case,'include|i_c=s'=>\$include)|| warn "perl prepare_bionano_output.pl --in IN --vcf VCF --tsv TSV --case g\n";
chomp ($infile,$vcf_file,$tsv_file,$case,$include);

# read the infile, split into useful and not useful parts
open(IN,$infile)|| die "could not open infile : $!\n";
my@in_csv_file=<IN>;

# ucsc genome browser : url building parts
my$url_start="https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg38&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=";
my$url_end="&hgsid=856936115_4aA1WK7ohwZycxkESPOCzhHUjzEm";


open(VCF,">",$vcf_file)|| die "could not open VCF outfile, exiting... \n";
open(TSV,">",$tsv_file)|| die "could not open TSV outfile, exiting... \n";

############### the from user desired extra columns ###############
my$need_to_get_header=0;
my@header_names=();# here we collect the names of the infile header
my$include_add_cols="";
my@splitted_numcols_by_comma=split(",",$include);
#my$test=join(/;/,@splitted_numcols_by_comma);
#print "all parts of infile wanted detected::$test::\n";
if(scalar(@splitted_numcols_by_comma)>0){
      # get the column name of the original infile
      $need_to_get_header=1;
}
#################### germline case ################################
if($case eq "g"){# germluine case, the raw infile has the parents information
  # print header into TSV file, vcf file does not need headers
  foreach my $infile_line (@in_csv_file){
        chomp $infile_line;
        if($infile_line=~/#/){
              # header
              if($need_to_get_header){
                    print "collecting extra header columns for @splitted_numcols_by_comma\n...";
                    my@lineparts=split(/\;/,$infile_line);
                    foreach my $header_name_num (@splitted_numcols_by_comma){
                          push(@header_names,"$lineparts[$header_name_num]")
                    }
              }
              my$user_added_head=join("\t",@header_names);
              # print header line only once we have a # starting line after we collected the extra column names aswell
              print TSV "coords\talt_chrom_trans\ttype\tconfidence\tzygosity\tSVFreq\tSVSize\tOverlapGenes\tnonOverlapgene\tNonOveGenedist\tfound_in_parents\tm_mole_found\tf_mole_found\tself_mole_found\tucsc_link\tQryContigID\t$user_added_head\n";

        }
        else{
              my$coords_done="";
              # non-header
              my@lineparts=split(/\;/,$infile_line);
              # has parts : SmapEntryID	QryContigID	RefcontigID1	RefcontigID2	QryStartPos	QryEndPos	RefStartPos	RefEndPos	Confidence	Type	XmapID1	XmapID2	LinkID	QryStartIdx	QryEndIdx	RefStartIdx	RefEndIdx	Zygosity	Genotype	GenotypeGroup	RawConfidence	RawConfidenceLeft	RawConfidenceRight	RawConfidenceCenter	SVsize	SVfreq	orientation	Sample	Algorithm	Size	Present_in_%_of_BNG_control_samples	Present_in_%_of_BNG_control_samples_with_the_same_enzyme	Fail_assembly_chimeric_score	OverlapGenes	NearestNonOverlapGene	NearestNonOverlapGeneDistance	PutativeGeneFusion	Found_in_parents_assemblies	Found_in_parents_molecules	Found_in_self_molecules	Mother_molecule_count	Father_molecule_count	Self_molecule_count
              # get coordinates- for vcf and for coords format that will be attached later
              my$chr=$lineparts[2];# RefcontigID1
              my$start_raw=$lineparts[6];# RefStartPos
              my$end_raw=$lineparts[7];
              # now regex out the .?? at the end of start and end
              $start_raw=~s/\.[0-9]//;
              $end_raw=~s/\.[0-9]//;

              # check for end >start for all coordinates. if not, reverse the order

              if($end_raw < $start_raw){
                # here we have to swap the direction
                $coords_done="chr"."$chr".":"."$end_raw"."-"."$start_raw";
                print VCF "$chr\t$end_raw\t$start_raw\t0\t0\n";
              }
              elsif($end_raw >= $start_raw){
                # stay in old order
                $coords_done="chr"."$chr".":"."$start_raw"."-"."$end_raw";
                print VCF "$chr\t$start_raw\t$end_raw\t0\t0\n";
              }
              else{
                print "line : $infile_line\n error: could not decide for what is bigger, start or end \n";
              }

              my$full_url="$url_start"."$coords_done"."$url_end";


              #print "check: $coords_done\n";
              my$extra_cols;
              foreach my $col_num_to_incl (@splitted_numcols_by_comma){
                    # attach the numbered column into the
                    $extra_cols=$extra_cols."\t"."$lineparts[$col_num_to_incl]";
              }

              # for the .vcf file we are done here. need to print header and attach 00 to it, then we are done with that
              # now the other tsv file: here we collect all usefulk columns and hand them later over to bionanotator.pl for final merging
              # column numbers we want besides the coordinates:
                                # edited coords raw | alt chrom trans |type           |confidence|      zygocity|         SVFreq|     SVSize|        OverlapGenes|nonOverlapgene|   NonOveGenedist|found_in_parents|m_mole_found|f_mole_found|self_mole_found|      aded ucsc url | ID of contig| user added things
              my$full_line_things_needed="$coords_done\t$lineparts[3]\t$lineparts[9]\t$lineparts[8]\t$lineparts[17]\t$lineparts[25]\t$lineparts[29]\t$lineparts[33]\t$lineparts[34]\t$lineparts[35]\t$lineparts[37]\t$lineparts[40]\t$lineparts[41]\t$lineparts[42]\t$full_url\t$lineparts[1]$extra_cols\n";
              #$full_line_things_needed=~s/\t\w/\t/;
              $full_line_things_needed=~s/\n\s/\n/;
              $full_line_things_needed=~s/\r//;
              $full_line_things_needed=~s/\n.*/\n/;
              print TSV $full_line_things_needed;


        }
  }
}
if($case eq "s"){# somatic case, the raw infile has no parents information
  # print header into TSV file, vcf file does not need headers
  foreach my $infile_line (@in_csv_file){
        chomp $infile_line;
        if($infile_line=~/#/){
              # header
              # header
              if($need_to_get_header){
                    print "collecting extra header columns for @splitted_numcols_by_comma\n...";
                    my@lineparts=split(/\;/,$infile_line);
                    foreach my $header_name_num (@splitted_numcols_by_comma){
                          push(@header_names,"$lineparts[$header_name_num]")
                    }
              }
              my$user_added_head=join("\t",@header_names);
              # print header line only once we have a # starting line after we collected the extra column names aswell
              print TSV "coords\talt_chrom_trans\ttype\tconfidence\tSVFreq\tSVSize\tOverlapGenes\tnonOverlapgene\tNonOveGenedist\tself_mole_found\tucsc_link\tQryContigID\t$user_added_head\n";

        }
        else{
              my$coords_done="";# start with empty- then fillup
              # non-header
              my@lineparts=split(/\;/,$infile_line);
              # has parts : SmapEntryID	QryContigID	RefcontigID1	RefcontigID2	QryStartPos	QryEndPos	RefStartPos	RefEndPos	Confidence	Type	XmapID1	XmapID2	LinkID	QryStartIdx	QryEndIdx	RefStartIdx	RefEndIdx	Zygosity	Genotype	GenotypeGroup	RawConfidence	RawConfidenceLeft	RawConfidenceRight	RawConfidenceCenter	SVsize	SVfreq	orientation	Sample	Algorithm	Size	Present_in_%_of_BNG_control_samples	Present_in_%_of_BNG_control_samples_with_the_same_enzyme	Fail_assembly_chimeric_score	OverlapGenes	NearestNonOverlapGene	NearestNonOverlapGeneDistance	PutativeGeneFusion	Found_in_parents_assemblies	Found_in_parents_molecules	Found_in_self_molecules	Mother_molecule_count	Father_molecule_count	Self_molecule_count
              # get coordinates- for vcf and for coords format that will be attached later
              my$chr=$lineparts[2];# RefcontigID1
              my$start_raw=$lineparts[6];# RefStartPos
              my$end_raw=$lineparts[7];
              # now regex out the .?? at the end of start and end
              $start_raw=~s/\.[0-9]//;
              $end_raw=~s/\.[0-9]//;
              # check for transversions- here we have to swap start and end in vcf
              if($end_raw < $start_raw){
                # here we have to swap the direction
                $coords_done="chr"."$chr".":"."$end_raw"."-"."$start_raw";
                print VCF "$chr\t$end_raw\t$start_raw\t0\t0\n";
              }
              elsif($end_raw >= $start_raw){
                # stay in old order
                $coords_done="chr"."$chr".":"."$start_raw"."-"."$end_raw";
                print VCF "$chr\t$start_raw\t$end_raw\t0\t0\n";
              }
              else{
                print "line : $infile_line\n error: could not decide for what is bigger, start or end of SV \n";
              }
              # for the .vcf file we are done here. need to print header and attach 00 to it, then we are done with that
              my$full_url="$url_start"."$coords_done"."$url_end";

              my$extra_cols;
              foreach my $col_num_to_incl (@splitted_numcols_by_comma){
                    # attach the numbered column into the
                    $extra_cols=$extra_cols."\t"."$lineparts[$col_num_to_incl]";
              }
              # now the other tsv file: here we collect all usefulk columns and hand them later over to bionanotator.pl for final merging
              # column numbers we want besides the coordinates:
                                # edited coords raw | alt chrom trans |type           |confidence|             SVFreq|     SVSize|        OverlapGenes|nonOverlapgene|   NonOveGenedist|found_in_parents|m_mole_found|f_mole_found|self_mole_found|
              my$full_line_things_needed="$coords_done\t$lineparts[3]\t$lineparts[9]\t$lineparts[8]\t$lineparts[22]\t$lineparts[26]\t$lineparts[30]\t$lineparts[31]\t$lineparts[32]\t$lineparts[35]\t$full_url\t$lineparts[1]$extra_cols\n";
              #$full_line_things_needed=~s/\t\w/\t/;
              $full_line_things_needed=~s/\n\s/\n/;
              $full_line_things_needed=~s/\r//;
              $full_line_things_needed=~s/\n.*/\n/;
              print TSV $full_line_things_needed;


        }
  }
}
print "splitted $infile into $vcf_file and $tsv_file\n";
1;

# TODO: add ucsc genome browser clickable link
# add following things to the output back again:
#svfreq
#qrycontigID
# thats the second column in both cases to simply add, maybe at the end?
# we need to simply insert the coords column into a last column with the url we wil build
