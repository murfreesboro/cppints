#! /usr/bin/perl -w

use strict;
use warnings;
use diagnostics;
use File::Basename;

#
# this is used to generate the cmake file: CMakeLists.txt
# for the given integral module
#

my $num_args = $#ARGV + 1;
if ($num_args != 2) {
	print "\nUsage: cmakegen.pl folder_name order \n";
	print "folder_name is the top integral folder, for example; /home/mike/project/gints_engine/hgp_os/ \n";
	print "inside this folder, you should have three sub folders, namely energy, first_deriv and second_deriv \n";
	print "order is between 0-2; which chooses the energy, 1st deriv or 2ed deriv sub folder\n";
	exit;
}

# let's try to get the work order
my $folder = $ARGV[0];
if (! -d "$folder") {
	print "input folder is $folder\n";
	die "the input folder does not exist!!\n";
}

# set the work order
my $order = $ARGV[1];
my $src_folder;
my $work_folder = "$folder";
if ($order == 0) {
	$src_folder = "\${PROJECT_SOURCE_DIR}/src/gints_engine/hgp_os/energy";
	$work_folder = "$work_folder"."/energy";
}elsif ($order == 1) {
	$src_folder = "\${PROJECT_SOURCE_DIR}/src/gints_engine/hgp_os/first_deriv";
	$work_folder = "$work_folder"."/first_deriv";
}elsif ($order == 2) {
	$src_folder = "\${PROJECT_SOURCE_DIR}/src/gints_engine/hgp_os/second_deriv";
	$work_folder = "$work_folder"."/second_deriv";
}else{
	print "input order is $order\n";
	die "the input order is incorrect, only between 0 to 2!!\n";
}

# check whether work folder exist?
if (! -d "$work_folder") {
	print "input folder is $work_folder\n";
	die "the work folder does not exist for order $order!!\n";
}

# now let's start the file
my $file = "CMakeLists.txt";
my $project = "SRC_GINTS_ENGINE_CXX";
open (WRITE, ">$file")  || die "can not open the $file for writinging!\n";
print (WRITE "#\n");
if ($order == 0) {
	print (WRITE "# build the analytical integral engine for energy\n");
	$project = "$project"."_D0";
}elsif ($order == 1) {
	print (WRITE "# build the analytical integral engine for 1st derivatives\n");
	$project = "$project"."_D1";
}elsif ($order == 2) {
	print (WRITE "# build the analytical integral engine for 2ed derivatives\n");
	$project = "$project"."_D2";
}
print (WRITE "#\n");
my $s = "set($project";
print (WRITE "$s\n");

# loop over the the top work folder
my @cppfiles = glob("$work_folder/*.cpp");
my $f;
foreach $f(@cppfiles) {
	my $f1 = basename($f);
	print (WRITE  "   $src_folder/$f1\n");
}

# now go into each sub folder to get it
my @folders = glob("$work_folder/*");
my $d;
foreach $d(@folders) {
	if ( -d "$d") {
		my @cppfiles = glob("$d/*.cpp");
		my $f;
		foreach $f(@cppfiles) {
			my $f1 = basename($f);
			print (WRITE  "   $src_folder/$f1\n");
		}
	}
}
$s = "   )";
print (WRITE "$s\n");

# add in library
$s = "add_library(gints_engine_d0 \${$project})";
if ($order == 1) {
	$s = "add_library(gints_engine_d1 \${$project})";
}elsif ($order == 2) {
	$s = "add_library(gints_engine_d2 \${$project})";
}
print (WRITE "$s\n");

# close the file
close(WRITE);

