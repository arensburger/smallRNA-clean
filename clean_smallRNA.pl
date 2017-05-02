#!/usr/bin/perl
# August 2015 Modification of the scrip clean_single.pl modified for small RNAs
use strict;
require File::Temp;
use File::Temp ();
use File::Basename;
use Getopt::Long;
use File::Path;

# CONSTANTS
# adapter trimming parameters
my $TRIMMOMATIC_PATH = "./Trimmomatic-0.32"; # java program location, must be in same folder as this script for java reasons I don't get
my $SORTMELOC = "/home/peter/smallRNA-clean/sortmerna"; # must keep full path name for sortme program
my $MINLEN = 19;
my $threeprimeadapter = "AGATCGGAAGAGCACACGTCT";

# other parameters
my $readspair1; #filename of first pair fastq
my $outputname; #base name for output files
my $outputdir; #outputdirectory
my $threads = `grep -c processor /proc/cpuinfo`; #number of threads to use
$threads =~ s/\s//g
;
#return date and time
sub datetime {
	use POSIX qw/strftime/;
	return (strftime('%D %T',localtime));
}

#####read and check the inputs
GetOptions(
	'f:s'   => \$readspair1,
	'o:s'	=> \$outputname,
	'd:s'	=> \$outputdir,
	't:s'   => \$threads,
);

unless ($readspair1) {
	die ("usage: perl clean_smallRNA.pl -f <REQUIRED: FASTQ file, must end in .fq> -o <OPTIONAL output name> -d <OPIONAL output directory (default current directory)> -t <OPTIONAL: number of threads to use, default $threads>\n");
}
unless ($outputname) {
	$outputname = basename($readspair1, ".fq"); #take the name from first member of the pair
}

#test if output directory has been specified, if not set output to current directoy
unless ($outputdir) {
	$outputdir = `pwd`; #set output to current directoy
	chomp $outputdir;
}

#test if output directory exists if not, create it
unless ( -e $outputdir )
{
	my $dir_test = mkdir($outputdir, 0777);
	unless ($dir_test) { 
		die "cannot create directory $outputdir\n";
	}
}

#create temporary files remove paired
my $unpaired_output = File::Temp->new( UNLINK => 1, SUFFIX => '.fastq' ); # temporary file with chastity results pair1

#create log file
if ( -e $outputdir) {
	open (LOG, ">$outputdir/Log.txt") or die ("cannot create $outputdir/Log.txt");
}
else {
	open (LOG, ">Log.txt") or die ("cannot create Log.txt");
}

########### start cleanning ######################
print LOG datetime, " Initial count\n";
print LOG datetime, " File $readspair1, total FASTQ reads: ", count_fastq($readspair1, "fq"), "\n"; # do a basic count
print LOG "\n";

print LOG datetime, " MD5sum: ", md5sum($readspair1), "\n";

# clip adapters (needs to be done before collapsing because it works on fastq files
print LOG datetime, " Clipping adapters and quality filter with Trimmomatic\n";
my $noadapter_name = clipadapters($readspair1);
print LOG datetime, " File with unpaired reads, FASTQ reads: ", count_fastq($noadapter_name, "fq"), "\n";
print LOG "\n";

# collapse the reads, removing duplicates
print LOG datetime, " Collapsing reads, removing duplicates\n";
my $collapsed_name = $outputname . "-collapsed.fa";
#`cp $noadapter_name noadapterfile.fq`;
`fastx_collapser -i $noadapter_name -o $collapsed_name -Q33`;
if ($? < 0) {
	print LOG "fastx_collapser did not run correctly, aborting\n";
	die "fastx_collapser did not run correctly\n";
}
print LOG datetime, " File with collapsed reads, FASTA reads: ", count_fastq($collapsed_name, "fa"), "\n";
print LOG "\n";

# removing 3' adapters
print LOG datetime, " Removing 3' adapter\n";
my $nothreeprime_name = $outputname . "-nothreeprime.fa";
`fastx_clipper -a $threeprimeadapter -c -l $MINLEN -i $collapsed_name -o $nothreeprime_name`;
if ($? < 0) {
        print LOG "fastx_clipper did not run correctly, aborting\n";
        die "fastx_clipper did not run correctly\n";
}
print LOG datetime, " File with 3' adapters removed, FASTA reads: ", count_fastq($nothreeprime_name, "fa"), "\n";
print LOG "\n";

#artifact filter
print LOG datetime, " Removing artifacts from unpaired file\n";
filter_artifact($nothreeprime_name);
print LOG datetime, " File with unpaired reads, FASTA reads: ", count_fastq($nothreeprime_name, "fa"), "\n";
print LOG "\n";

#remove ribosomal sequence
print LOG datetime, " Removing ribosomal sequences using SortMeRNA\n";
my $noribo_name = $outputname . "-noribo";
my $noribo_name2 = $noribo_name . ".fa"; # needed because sortmRNA adds the .fa extension
ribosome_removal($nothreeprime_name, $noribo_name);
print LOG datetime, " File with ribosomes removed, FASTA reads: ", count_fastq($noribo_name, "fa"), "\n";
print LOG "\n";

#print the data files 
my $unpairedoutname = $outputdir . "/" . $outputname .  ".fa";
`mv $noribo_name $unpairedoutname`;
`mv $collapsed_name $outputdir`;
`mv $nothreeprime_name $outputdir`;
print LOG datetime, " Data file are written in $unpairedoutname\n";
print LOG "\n";
print STDERR "Done\n";

######## subroutines ###################
sub count_fastq { # now can count fasta too
#	print STDERR "count fastq...\n";
	my ($title, $type) = @_;
	my $rawcount;

	my $txtcount = `wc -l $title`;
	if ($txtcount =~ /^(\d+)\s/) {
		my $count;
		my $value = $1;
		if ($type eq "fq") {
			$count = $value/4;
		}
		elsif ($type eq "fa") {
			$count = $value/2;
		}
		else {
			die "did not recognize type $type in subroutine count_fastq\n";
		}
		return (commify($count));
	}
	else {
		die "could not count lines in file $title using wc -l\n";
	}
}

sub clipadapters {
	my ($inputfile1) = @_;
	my $forward_unpaired = File::Temp->new( UNLINK => 1, SUFFIX => '.fastq' );
	my $trim_run = `java -jar $TRIMMOMATIC_PATH/trimmomatic-0.32.jar SE -threads $threads -phred33 $inputfile1 $forward_unpaired ILLUMINACLIP:$TRIMMOMATIC_PATH/illuminaClipping.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:$MINLEN`;
	if ($? < 0) {
		print LOG "Trimmomatic did not run correctly, aborting\n";
		die "Trimmomatic did not run correctly\n";
	}
	return($forward_unpaired);
}

sub ribosome_removal{
	my ($infile, $outfile) = @_;
	open (INPUT1, $infile) or die "cannot open file $infile\n";
        my $ribooutput = File::Temp->new( UNLINK => 1);
	my $ribooutput2 = File::Temp->new( UNLINK => 1);
	`$SORTMELOC/sortmerna --ref $SORTMELOC/rRNA_databases/silva-arc-23s-id98.fasta,$SORTMELOC/index/silva-arc-23s-db:$SORTMELOC/rRNA_databases/misc_rRNA.fasta,$SORTMELOC/index/misc_rRNA-db:$SORTMELOC/rRNA_databases/silva-bac-16s-id90.fasta,$SORTMELOC/index/silva-bac-16s-db:$SORTMELOC/rRNA_databases/rfam-5.8s-database-id98.fasta,$SORTMELOC/index/rfam-5.8s-database-db:$SORTMELOC/rRNA_databases/silva-bac-23s-id98.fasta,$SORTMELOC/index/silva-bac-23s-db:$SORTMELOC/rRNA_databases/rfam-5s-database-id98.fasta,$SORTMELOC/index/rfam-5s-db:$SORTMELOC/rRNA_databases/silva-euk-18s-id95.fasta,$SORTMELOC/index/silva-euk-18s-db:$SORTMELOC/rRNA_databases/silva-euk-28s-id98.fasta,$SORTMELOC/index/silva-euk-28s-db:$SORTMELOC/rRNA_databases/silva-arc-16s-id95.fasta,$SORTMELOC/index/silva-arc-16s-db --reads $infile --other $ribooutput2 -a $threads --fastx --aligned $ribooutput --sam`;
        if ($? < 0) {
                print LOG "sortmeRNA did not run correctly, aborting\n";
                die "sortmeRNA did not run correctly\n";
        }

	my $tempname = "$ribooutput2" . ".fa"; #necessary because sortmeRNA adds .fa to the output file without asking
	`mv $tempname $outfile`;	
}

sub filter_artifact {
	my ($inputfile) = @_;
	my $tx = File::Temp->new( UNLINK => 1, SUFFIX => '.fastq' ); # temporary file with hits
	`fastx_artifacts_filter -i $inputfile -o $tx -Q 33`;

	my $wctext = `wc -m $tx`; #word count of the ouptput used to check that fastx ran
	if ($wctext =~ /^0\s/) {
		die "There was an error, either fastx_artifacts_filter is not installed or all the sequences were removed as artifacts\n";
	}
	else {
		`mv $tx $inputfile`;
	}
}

# from perl cookbook, to put commas on big numbers
sub commify {
    my $text = reverse $_[0];
    $text =~ s/(\d\d\d)(?=\d)(?!\d*\.)/$1,/g;
    return scalar reverse $text
}

sub md5sum {
	my ($filename) = @_;
	my $textout = `md5sum $filename`;
	return($textout);
}

