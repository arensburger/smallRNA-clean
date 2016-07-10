#!/usr/bin/perl

# July 2016 Takes the output of cleanning script where md5 are reported and checks which ones are identical

use strict;
use Getopt::Long;
use Cwd;

my $dir = cwd(); # by default the directory is the currrent workding directory
my $help;

#####read and check the inputs
GetOptions(
	'd:s'   => \$dir,
	'h:s'   => \$help
);

if ($help) {
	die "usage: perl checkmd5 -d <directory to look for log file, current directory by default>, -h\n";
}

# check the logs
my $logs = `find $dir -name 'Log.txt' | xargs grep -h 'MD5sum'`;
if ($logs eq '') {
	die "did not find any Log.txt files in directory $dir\n";
}

my %md5freq; # md5 number as key and occurence as value

# parse the log output
my @data = split ("\n", $logs);
my $files; #names of the files checked
foreach my $line (@data) {
	if ($line =~ /MD5sum:\s+(\S+)\s+(\S+)/) {
		$md5freq{$1} += 1;
		$files .= $2 . "\t" . $1 . "\n";
		
	}
	else {
		die "cannot read line\n$line";	
	}
}

# print the results
print "checked files\n$files";
my $allclear=1;
foreach my $key (keys %md5freq) {
	if ($md5freq{$key} > 1) {
		print "\n$key occurs $md5freq{$key} times\n";
		$allclear = 0;
	}
}
if ($allclear) {
	print "\nAll clear\n";
}



