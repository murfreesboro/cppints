#! /usr/bin/perl -w

use strict;
use warnings;
use diagnostics;

#
# after using the eval.py, we will get a data file;
# here this file is used to collect the data for 
# clear compare
#
#
my $file1 = $ARGV[0];
my $file2 = $ARGV[1];

my $pattern = "total K4 FLOPS";

# get the file list
my @filelist;
open (READ, "$file1")  || die "can not open the $file1 for reading!\n";
while(<READ>) {
	if ($_ =~ /for file/) {
		my $line = $_;
		$line =~ s/^\s+//;
		my @tmp = split /\s+/, $line;
		push(@filelist,$tmp[2]);
	}
}
close(READ);

# now print out stuff
print "for the first file\n";
open (READ, "$file1")  || die "can not open the $file1 for reading!\n";
my $iFile;
foreach $iFile(@filelist) {
	seek(READ,0,0);
	my $getit = 0;
	while(<READ>) {
		if ($_ =~ /$iFile/) {
			$getit = 1;
			while(<READ>) {
				if ($_ =~ /$pattern/) {
					print "$iFile\t", $_;
					last;
				}
			}
			last;
		}
	}
	if ($getit == 0) {
		die "did not find the data section for the file: $iFile\n";
	}	
}
close(READ);

# now print out stuff
print "for the second file\n";
open (READ, "$file2")  || die "can not open the $file2 for reading!\n";
foreach $iFile(@filelist) {
	seek(READ,0,0);
	my $getit = 0;
	while(<READ>) {
		if ($_ =~ /$iFile/) {
			$getit = 1;
			while(<READ>) {
				if ($_ =~ /$pattern/) {
					print "$iFile\t", $_;
					last;
				}
			}
			last;
		}
	}
	if ($getit == 0) {
		die "did not find the data section for the file: $iFile\n";
	}	
}
close(READ);

