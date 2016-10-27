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
if ($num_args != 3) {
	print "\nUsage: cmakegen.pl folder_name order module \n";
	print "folder_name is the top integral folder, for example; /home/mike/project/gints_engine/hgp_os/ \n";
	print "in the folder you should have three sub folders, namely energy, first_deriv and second_deriv \n";
	print "order is between 0-2; which chooses the energy, 1st deriv or 2ed deriv sub folder\n";
	print "module is the name of integral module, for example; twobodyoverlap\n";
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

# let's see whether entry cpp file exist?
my $module = $ARGV[2];
my $entry_cpp_file = "hgp_os_"."$module";
if ($order == 0) {
	$entry_cpp_file = "$entry_cpp_file".".cpp";
}elsif ($order == 1) {
	$entry_cpp_file = "$entry_cpp_file"."_d1.cpp";
}elsif ($order == 2) {
	$entry_cpp_file = "$entry_cpp_file"."_d2.cpp";
}
my $cpp_file = "$work_folder"."/"."$entry_cpp_file";
if (! -e "$cpp_file") {
	print "input module is $module\n";
	print "top cpp file is $cpp_file\n";
	die "however, the top cpp file for running the integral module does not exist!!\n";
}
$entry_cpp_file = "$src_folder"."/"."$entry_cpp_file";

# append the module to the work folder
# check whether work folder exist?
$work_folder = "$work_folder"."/"."$module";
if (! -d "$work_folder") {
	print "input module is $work_folder\n";
	die "the work folder does not exist for order $order!!\n";
}

# now let's start the file
my $file = "CMakeLists.txt";
my $project  = "SRC_GINTS_ENGINE_";
my $libname  = "gints_engine_"."$module"."_";
open (WRITE, ">$file")  || die "can not open the $file for writinging!\n";
print (WRITE "#\n");
my $cap_module = $module;
$cap_module =~ tr/a-z/A-Z/;
if ($order == 0) {
	print (WRITE "# build the analytical integral engine for energy of $module\n");
	$project = "$project"."$cap_module"."_D0";
	$libname = "$libname"."d0";
}elsif ($order == 1) {
	print (WRITE "# build the analytical integral engine for 1st derivatives of $module\n");
	$project = "$project"."$cap_module"."_D1";
	$libname = "$libname"."d1";
}elsif ($order == 2) {
	print (WRITE "# build the analytical integral engine for 2ed derivatives of $module\n");
	$project = "$project"."$cap_module"."_D2";
	$libname = "$libname"."d2";
}
print (WRITE "#\n");
my $s = "set($project";
print (WRITE "$s\n");

# now list all of files
print (WRITE  "   $entry_cpp_file\n");
my @cppfiles = glob("$work_folder/*.cpp");
my $f;
foreach $f(@cppfiles) {
	my $f1 = basename($f);
	print (WRITE  "   $src_folder/$module/$f1\n");
}
$s = "   )";
print (WRITE "$s\n");

# add in library
$s = "add_library($libname \${$project})";
print (WRITE "$s\n");

# close the file
close(WRITE);

