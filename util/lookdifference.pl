#! /usr/bin/perl -w

use strict;
use warnings;
use diagnostics;

#
# take a look for the difference of comparison made
# file is the log file to read the comparison difference
# data type is either double or float
#

my $file;
my $data_type = "double";
if ("$#ARGV" == 0) {
	$file = $ARGV[0];
}elsif("$#ARGV" == 1) {
	$file = $ARGV[0];
	$data_type = $ARGV[1];
} else {
	die "input parameter is wrong!\n";
}

# check data type
if ($data_type ne "double" && $data_type ne "float") {
	die "input data type is invalid!\n";
}

# set up the thresh value
# for double we set it as 10^-7
# for float we set it as 10^-3
my $thresh = 0.0000001;
if ($data_type eq "float") {
	$thresh = 0.001;
}

# prepare file writing
open (READ, "$file") || die "can not open the $file for reading to see comparison difference!!!\n";
while(<READ>) {
	if ($_ =~ /difference/) {
		my @data = split /\s+/, $_;
		my $d = $data[2];
		if ($d > $thresh) {
			print $_, "\n";
		}
	}
}
close(READ);


