//
// CPPINTS: A C++ Program to Generate Analytical Integrals Based on Gaussian
// Form Primitive Functions
//
// Copyright (C) 2012-2014 Fenglai Liu
// This softare uses the MIT license as below:
//
//	Permission is hereby granted, free of charge, to any person obtaining 
//	a copy of this software and associated documentation files (the "Software"), 
//	to deal in the Software without restriction, including without limitation 
//	the rights to use, copy, modify, merge, publish, distribute, sublicense, 
//	and/or sell copies of the Software, and to permit persons to whom the Software 
//	is furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in all 
// copies or substantial portions of the Software.
//						    
//	THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, 
//	INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR 
//	PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE 
//	FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, 
//	ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
//
//
#include <boost/filesystem.hpp>
#include <boost/algorithm/string.hpp>   // string handling
#include <boost/algorithm/string/split.hpp>
#include <boost/lexical_cast.hpp>
#include "shellsymbol.h"
#include "inttype.h"
#include "infor.h"
using namespace boost::filesystem;
using namespace boost;
using namespace inttype;
using namespace infor;

LineParse::LineParse(string input) 
{
	line = input;
	trim(line);
	number = 0;
	isEmpty = false;
	isComment = false;
	if(line.empty()){ 
		isEmpty = true;
	} else if(line.at(0) == '#' || line.at(0) == '!') {
		isComment = true;
	} else { 
		// two steps: 
		// first one is to check whether it contains some comment 
		// and get rid of the comment part
		// second, split the substr without a comment
		vector<string>sstr;
		split(sstr, line, is_any_of("#") || is_any_of("!"));
		string s=sstr[0];
		trim(s);
		split(pieces, s, is_any_of(" ") || is_any_of("="), token_compress_on);
		number = pieces.size();
	}
}

bool WordConvert::toInt(string s, int& x) {
	try {
		x = lexical_cast<int>(s);
	} catch(boost::bad_lexical_cast& error) {
		return false;
	}
	return true;
}

bool WordConvert::compare(string s1, string s2){
	to_upper(s1);
	to_upper(s2);
	trim(s1);
	trim(s2);
	if (s1 == s2) {
		return true;
	} else {
		return false;
	}
}

void WordConvert::capitalize(string& s){
	to_upper(s);
}

bool WordConvert::isInt(string s) {
	if (s[0] == '+' || s[0] == '-') {
		for(int i=1;i<(int)s.size();i++)
			if (!isdigit(s[i])) return false;
	}else{
		for(int i=0;i<(int)s.size();i++)
			if (!isdigit(s[i])) return false;
	}
	return true;
}

Infor::Infor(const string& input):enforceVRR(NO_ENFORCE_ON_RR),
	enforceHRR(NO_ENFORCE_ON_RR),maxL(5),auxMaxL(6),derivOrder(0),
	maxL_vrrPrinting(12),maxL_hrrPrinting(10),maxL_singleFile(12),
	vec_form(USE_SCR_VEC),M_limit(10),fmt_error(12),hrr_method("hgp"),
	vrr_method("os")
{ 
	// open the input file
	ifstream inf;
	const char * input_file = input.c_str();
	inf.open(input_file,ios::in);
	if (!inf) {
		crash(true, "Can not open input in infor class");
	}

	// now let's process each line of information
	string line;
	WordConvert w;
	inf.seekg(0,ios::beg);
	while(getline(inf,line)) {

		// read in this line
		LineParse l(line);

		// shall we bypass?
		if (l.isCom() || l.isEmp()) continue;

		// vrr printing option
		if (w.compare(l.findValue(0), "vrr_print_option")) {
			string value = l.findValue(1);
			w.capitalize(value);
			if (value == "ARRAY") {
				enforceVRR = ENFORCE_RR_WITH_ARRAY;
			}else if (value == "VAR") {
				enforceVRR = ENFORCE_RR_WITH_VAR;
			}else{
				crash(true, "Invalid option given in processing vrr_print_option");
			}
		}

		// hrr printing
		if (w.compare(l.findValue(0), "hrr_print_option")) {
			string value = l.findValue(1);
			w.capitalize(value);
			if (value == "ARRAY") {
				enforceHRR = ENFORCE_RR_WITH_ARRAY;
			}else if (value == "VAR") {
				enforceHRR = ENFORCE_RR_WITH_VAR;
			}else{
				crash(true, "Invalid option given in processing hrr_print_option");
			}
		}

		// check the HRR algorithm
		if (w.compare(l.findValue(0), "hrr_method")) {
			string value   = l.findValue(1);
			w.capitalize(value);
			if (value == "HGP") {
				hrr_method = "hgp";
			}else if (value == "NONE") {
				hrr_method = "none";
			}else{
				crash(true, "Invalid name given for HRR method.");
			}
		}

		// check the VRR algorithm
		if (w.compare(l.findValue(0), "vrr_method")) {
			string value   = l.findValue(1);
			w.capitalize(value);
			if (value == "OS") {
				vrr_method = "os";
			}else{
				crash(true, "Invalid name given for VRR method.");
			}
		}

		// highest L to create integral files
		if (w.compare(l.findValue(0), "maxl")) {
			string value   = l.findValue(1);
			int tmp = 0;
			if (!w.toInt(value,tmp)) {
				crash(true, "In Infor we can not process maxl. not an integer");
			}
			if (tmp >= S) {
				maxL = tmp;
			}else{
				crash(true, "Invalid maximum L value given for infor class, negative integer.");
			}
		}

		// highest L for aux shell to create integral files
		// this is used in 2/3 body of ERI
		if (w.compare(l.findValue(0), "aux_max_l")) {
			string value   = l.findValue(1);
			int tmp = 0;
			if (!w.toInt(value,tmp)) {
				crash(true, "In Infor we can not process aux_max_l. not an integer");
			}
			if (tmp >= S) {
				auxMaxL = tmp;
			}else{
				crash(true, "Invalid maximum L value given for aux_max_l, negative integer.");
			}
		}

		// to see what kind of job we are going to perform
		if (w.compare(l.findValue(0), "job_list")) {
			int nJobs = l.getNPieces()-1; 
			for(int i=1; i<=nJobs; i++) {
				string value = l.findValue(i);
				w.capitalize(value);
				int job = getOperIntName(value);
				joblist.push_back(job);
			}
		}

		// the highest derivatives order
		if (w.compare(l.findValue(0), "deriv_order")) {
			string value = l.findValue(1);
			int tmp = 0;
			if (!w.toInt(value,tmp)) {
				crash(true, "In Infor we can not process deriv_order. not an integer");
			}
			if (tmp <= 2 && tmp >= 0) {
				derivOrder = tmp;
			}else{
				crash(true, "Invalid derivatives order given for infor class.");
			}
		}

		// the highest angular momentum that we use variable rather 
		// than vector in the result cpp file for vrr part
		if (w.compare(l.findValue(0), "vrr_print_max_l")) {
			string value   = l.findValue(1);
			int tmp = 0;
			if (!w.toInt(value,tmp)) {
				crash(true, "In Infor we can not process vrr_print_max_l. not an integer");
			}
			if (tmp >= S) {
				maxL_vrrPrinting = tmp;
			}else{
				crash(true, "Invalid negative vrr_print_max_l value given for infor class.");
			}
		}

		// the highest angular momentum that we use variable rather 
		// than vector in the result cpp file for hrr part
		if (w.compare(l.findValue(0), "hrr_print_max_l")) {
			string value   = l.findValue(1);
			int tmp = 0;
			if (!w.toInt(value,tmp)) {
				crash(true, "In Infor we can not process hrr_print_max_l. not an integer");
			}
			if (tmp >= S) {
				maxL_hrrPrinting = tmp;
			}else{
				crash(true, "Invalid negative hrr_print_max_l value given for infor class.");
			}
		}

		// the highest angular momentum that we will keep everything
		// in a single file
		// else we need to split the cpp file into the vrr.cpp and hrr.cpp
		if (w.compare(l.findValue(0), "single_cpp_max_l")) {
			string value   = l.findValue(1);
			int tmp = 0;
			if (!w.toInt(value,tmp)) {
				crash(true, "In Infor we can not process single_cpp_max_l. not an integer");
			}
			if (tmp >= S) {
				maxL_singleFile = tmp;
			}else{
				crash(true, "Invalid negative single_cpp_max_l value given for infor class.");
			}
		}

		// what kind of vector the program will use?
		if (w.compare(l.findValue(0), "vector_form")) {
			string value = l.findValue(1);
			w.capitalize(value);
			if (value == "TBB_VEC") {
				vec_form = TBB_VEC;
			}else if (value == "STD_VEC") {
				vec_form = STD_VEC;
			}else if (value == "SCR") {
				vec_form = USE_SCR_VEC;
			}else{
				crash(true, "Invalid option given in processing vector_form");
			}
		}

		// set M limit
		if (w.compare(l.findValue(0), "m_limit")) {
			string value = l.findValue(1);
			int tmp = 0;
			if (!w.toInt(value,tmp)) {
				crash(true, "In Infor we can not process m_limit. not an integer");
			}
			if (tmp != 8 && tmp != 9 && tmp != 10) {
				crash(true, "Invalid M limit given for infor class, only 8, 9, 10 allowed right now.");
			}else{
				M_limit = tmp;
			}

			// set fmt error
			// this is closely related to the setting of M limit
			if (M_limit == 8 || M_limit == 9) {
				fmt_error = 13;
			}else if (M_limit == 10) {
				fmt_error = 12;
			}
		}
	}

	// close the input file
	inf.close();

	// we need to check the job list,
	// it should not be empty
	if (joblist.size() == 0) {
		crash(true,"Empty job list given in infor class");
	}

	// also we need to check that whether the maxl for single cpp file
	// is larger than the HRR maxl
	// this is used to make sure that HRR part is always in vector form
	if (maxL_singleFile<maxL_hrrPrinting) {
		crash(true,"max L for spliting cpp file should be at least same "
				"with max L used in HRR variable/vector setting");
	};

	// finally, create the peoject folder
	string project = getProjectName();
	path p(project.c_str());
	if (exists(p)) {
		remove_all(p);
	}
	create_directory(p);

	// create energy - first derivatives - second derivatives folder
	// for each derivative order, we also create the operator folder
	for(int i=0; i<=derivOrder; i++) {

		// get the folder name
		string folder = "energy";
		if (i == 1) folder = "first_deriv";
		if (i == 2) folder = "second_deriv";
		path p1(folder.c_str());

		// create the folder
		path p0(p);
		p0 /= p1;
		create_directory(p0);

		// now let's create each job folder
		for(int iJob = 0; iJob<(int)joblist.size(); iJob++) {
			int job = joblist[iJob];
			string jobname = getOperStringName(job);

			// for the two body and three body eri, the codes 
			// will be in the eri folder, too
			if (jobname == "TWOBODYERI" || jobname == "THREEBODYERI") {
				jobname = "ERI";
			};

			// now create the folder
			// if folder is already there, we bypass it
			to_lower(jobname);
			path p2(p0);
			path p3(jobname.c_str());
			p2 /= p3;
			if (exists(p2)) continue;
			create_directory(p2);
		}
	}
}

string Infor::getProjectName() const 
{
	string project;
	if (hasHRR()) {
		project = hrr_method + "_" + vrr_method;
	}else{
		project = vrr_method;
	}
	return project;
}

string Infor::getProjectFileDir(const int& oper, const string& fileName, 
		int order, bool withTmpWorkDir) const 
{
	// project name
	string project = getProjectName();
	path p0(project.c_str());

	// get the deriv name
	string folder = "energy";
	if (order == 1) folder = "first_deriv";
	if (order == 2) folder = "second_deriv";
	path p1(folder.c_str());

	// integral folder name
	string jobname = getOperStringName(oper);
	to_lower(jobname);
	path p2(jobname.c_str());

	// we note that it could be in the tmp work dir
	string workDir = "tmp_work_dir";
	path p3(workDir.c_str());

	// finally it's the integral file
	path p4(fileName.c_str());

	// now let's connect them together
	if (withTmpWorkDir) {
		p0 /= p1;
		p0 /= p2;
		p0 /= p3;
		p0 /= p4;
	}else{
		p0 /= p1;
		p0 /= p2;
		p0 /= p4;
	}
	return p0.string();
}

string Infor::getProjectTmpFileDir(const int& oper, int order) const 
{
	// project name
	string project = getProjectName();
	path p0(project.c_str());

	// get the deriv name
	string folder = "energy";
	if (order == 1) folder = "first_deriv";
	if (order == 2) folder = "second_deriv";
	path p1(folder.c_str());

	// integral folder name
	string jobname = getOperStringName(oper);
	to_lower(jobname);
	path p2(jobname.c_str());

	// now create the whole path
	string workDir = "tmp_work_dir";
	path p3(workDir.c_str());
	p0 /= p1;
	p0 /= p2;
	p0 /= p3;
	return p0.string();
}

int Infor::getVRRMethod() const 
{
	if (vrr_method == "os") {
		return OS;
	}else{
		crash(true, "Invalid VRR method in the getVRRMethod");
	}
	return -1;
}

