#! /usr/bin/perl -w

use strict;
use warnings;
use diagnostics;

#
# re-organize the generate entry cpp file by dividing 
# the long switch list into several sections according 
# to the angular momentum
#
sub findL; 

# set the hgp_os_eri.cpp
my $file = "hgp_os_eri.cpp";
my $out = "$file".".out";

# prepare
open (READ, "$file")  || die "can not open the $file for reading!\n";

# see what is the maxL?
my $maxL = 0;
while(<READ>) {
	if ($_ =~ /case/ && $_ !~ /default/) {
		my $line = <READ>;

		# truncate the ( part
		my @tmp = split /\(/, $line;

		# get the function name
		my $i = $tmp[0];
		@tmp = split /\s+/, $i;
		$i = $tmp[2];
		print $i, "\n";

		# now get the angular momentum
		@tmp = split /_/, $i;
		my $len = @tmp;
		my @am;
		for($i=3; $i<$len; $i=$i+1) {
			push(@am,$tmp[$i]);
		}

		# now find L
		foreach $i(@am) {
			my $L = findL($i);
			if ($L>$maxL) {
				$maxL = $L;
			}
		}
	}
}
#print $maxL, "\n";
if ($maxL == 1) {
	die "maxL is only SP/S/P, we do not need to do this\n";
}

# now let's output the result
# firstly it's the head part
# we stop by the { part
open (WRITE, ">$out") || die "can not open the $out for writing!\n";
seek(READ,0,0);
my $pos = 0;
while(<READ>) {
	print(WRITE $_);
	if ($_ =~ /^{/) {
		$pos = tell(READ);
		last;
	}
}

# S/P/SP is one section
# then D,F,G,H etc.
my $L;
my $start = 1;
for($L=$start; $L<=$maxL; $L = $L + 1) {

	# print the switch in output
	if ($L == 1) {
		my $line = "  if (maxL<=1) {\n";
		print WRITE $line;
		$line = "    switch(LCode) {\n";
		print WRITE $line;
	}else{
		my $line = "  }else if (maxL==$L) {\n";
		print WRITE $line;
		$line = "    switch(LCode) {\n";
		print WRITE $line;
	}

	# rewind the file handle
	seek(READ,$pos,0);
	while(<READ>) {
		if ($_ =~ /case/ && $_ !~ /default/) {

			# set the reading in lines
			my $line0 = $_;
			my $line1 = <READ>;
			my $line2 = <READ>;
			my $maxl = 0;

			# truncate the ( part
			my @tmp = split /\(/, $line1;

			# get the function name
			my $i = $tmp[0];
			@tmp = split /\s+/, $i;
			$i = $tmp[2];

			# now get the angular momentum
			@tmp = split /_/, $i;
			my $len = @tmp;
			my @am;
			for($i=3; $i<$len; $i=$i+1) {
				push(@am,$tmp[$i]);
			}

			# now find L
			foreach $i(@am) {
				my $l = findL($i);
				if ($l>$maxl) {
					$maxl = $l;
				}
			}

			# now let's see whether to print out lines
			if ($L == 1 && $maxl <= 1) {
				print WRITE $line0;
				print WRITE $line1;
				print WRITE $line2;
			}elsif ($maxl == $L) {
				print WRITE $line0;
				print WRITE $line1;
				print WRITE $line2;
			}
			
		}
	}

	#end part of switch code
	print WRITE "    case default:\n";
   print WRITE "#ifdef DEBUG\n";
   print WRITE "      printf(\"%s %lld\\n\",\"Un-recognized LCode in the integrals calculation \", LCode);\n";
   print WRITE "assert(0);\n";
   print WRITE "#endif\n";
   print WRITE "      return true;\n";
   print WRITE "      break;\n";
	print WRITE "    }\n";
	print WRITE "\n";
}
print WRITE "  }\n";
print WRITE "}\n";

close(READ);
close(WRITE);

sub findL {

	my $s = $_[0];
	if ($s eq "S" || $s eq "s") {
		return 0;
	}elsif ($s eq "SP" || $s eq "sp") {
		return 1;
	}elsif ($s eq "P" || $s eq "p") {
		return 1;
	}elsif ($s eq "D" || $s eq "d") {
		return 2;
	}elsif ($s eq "F" || $s eq "f") {
		return 3;
	}elsif ($s eq "G" || $s eq "g") {
		return 4;
	}elsif ($s eq "H" || $s eq "h") {
		return 5;
	}elsif ($s eq "I" || $s eq "i") {
		return 6;
	}elsif ($s eq "K" || $s eq "k") {
		return 7;
	}elsif ($s eq "L" || $s eq "l") {
		return 8;
	}elsif ($s eq "M" || $s eq "m") {
		return 9;
	}elsif ($s eq "N" || $s eq "n") {
		return 10;
	}else{
		die "highest L support in findL is 10, please add more\n";
	}
	return -1;
}
