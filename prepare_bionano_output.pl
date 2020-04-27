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
GetOptions('in|i=s'=>\$infile,'vcf|v=s'=>\$vcf_file,'tsv|t=s'=>\$tsv_file)|| warn "perl prepare_bionano_output.pl --in IN --vcf VCF --tsv TSV\n";
chomp ($infile,$vcf_file,$tsv_file);

# read the infile, split into useful and not useful parts
open(IN,$infile)|| die "could not open infile : $!\n";
my@in_csv_file=<IN>;


open(VCF,">",$vcf_file)|| die "could not open VCF outfile, exiting... \n";
open(TSV,">",$tsv_file)|| die "could not open TSV outfile, exiting... \n";


# print header into TSV file, vcf file does not need headers
print TSV "coords\talt_chrom_trans\ttype\tconfidence\tzygocity\tSVFreq\tSVSize\tOverlapGenes\tnonOverlapgene\tNonOveGenedist\tfound_in_parents\tm_mole_found\tf_mole_found\tself_mole_found\n";
foreach my $infile_line (@in_csv_file){
      chomp $infile_line;
      if($infile_line=~/#/){
            # header
      }
      else{
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
            my$coords_done="chr"."$chr".":"."$start_raw"."-"."$end_raw";
            #print "check: $coords_done\n";
            print VCF "$chr\t$start_raw\t$end_raw\t0\t0\n";
            # for the .vcf file we are done here. need to print header and attach 00 to it, then we are done with that
            # now the other tsv file: here we collect all usefulk columns and hand them later over to bionanotator.pl for final merging
            # column numbers we want besides the coordinates:
                              # edited coords raw | alt chrom trans |type           |confidence|      zygocity|         SVFreq|     SVSize|        OverlapGenes|nonOverlapgene|   NonOveGenedist|found_in_parents|m_mole_found|f_mole_found|self_mole_found|
            my$full_line_things_needed="$coords_done\t$lineparts[3]\t$lineparts[9]\t$lineparts[8]\t$lineparts[17]\t$lineparts[25]\t$lineparts[29]\t$lineparts[33]\t$lineparts[34]\t$lineparts[35]\t$lineparts[37]\t$lineparts[40]\t$lineparts[41]\t$lineparts[42]\n";
            #$full_line_things_needed=~s/\t\w/\t/;
            $full_line_things_needed=~s/\n\s/\n/;
            $full_line_things_needed=~s/\r//;
            $full_line_things_needed=~s/\n.*/\n/;
            print TSV $full_line_things_needed;


      }
}


print "splitted $infile into $vcf_file and $tsv_file\n";
1;
# TODO:
#- integrate into bionaotator
# paste into final outfile
# check
