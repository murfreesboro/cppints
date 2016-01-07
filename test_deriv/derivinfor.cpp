#include "libgen.h"
#include "textread.h"
#include "derivinfor.h"
#include "boost/lexical_cast.hpp"
using namespace textread;
using namespace derivinfor;


void DerivInforArray::readInformation(const string& fileName)
{

	// let me open the file to locate the code section
	ifstream inf;
	const char * input_file = fileName.c_str();
	inf.open(input_file,ios::in);
	if (!inf) {
		crash(true,"can not open the file in readInformation of derivInforArray");
	}

	// locate the code
	string line;
	inf.seekg(0,ios::beg);
	bool getIt = false;
	string scode = boost::lexical_cast<string>(code);
	while(getline(inf,line)) {
		LineParse l(line);
		if (l.getNPieces() > 0) {
			string val = l.findValue(0);
			if (val == scode) {
				getIt = true;
				break;
			}
		}
	}
	if (! getIt) {
		crash(true, "fail to get the cpp code in the given data file for extracting deriv array");
	}

	// now counting how many data fields
	WordConvert w;
	while(getline(inf,line)) {
		LineParse l(line);
		if (l.isEmp()) break;

		// whether this is line contains the position information?
		string p1 = l.findValue(0);
		w.capitalize(p1);
		if (p1 == "BRA1" || p1 == "BRA2" || p1 == "KET1" || p1 == "KET2") {

			// set up an empty deriv
			DerivInfor derivInfor(code,jobOrder);

			// get the pos1
			Int pos1 = -1;
			if (p1 == "BRA1") {
				pos1 = BRA1;
			}else if (p1 == "BRA2") {
				pos1 = BRA2;
			}else if (p1 == "KET1") {
				pos1 = KET1;
			}else {
				pos1 = KET2;
			}

			// let's see whether we have pos2
			Int pos2 = -1;
			if (l.getNPieces() >= 2) {
				string p2 = l.findValue(1);
				w.capitalize(p2);
				if (p2 == "BRA1") {
					pos2 = BRA1;
				}else if (p2 == "BRA2") {
					pos2 = BRA2;
				}else if (p2 == "KET1") {
					pos2 = KET1;
				}else {
					pos2 = KET2;
				}
			}

			// now push in pos infor
			derivInfor.updatePos(pos1,pos2);

			// now let's read in dir
			string line2;
			while(getline(inf,line2)) {

				// we break when we read in the end sign
				if (line2.find("##") != std::string::npos) {
					break;
				}

				// now prepare for read in
				LineParse l(line2);

				// read in first dir
				Int dir1 = -1;
				string p1 = l.findValue(0);
				if (p1 == "X" || p1 == "x") {
					dir1 = DERIV_X;
				}else if (p1 == "Y" || p1 == "y") {
					dir1 = DERIV_Y;
				}else {
					dir1 = DERIV_Z;
				}

				// read in second dir
				Int dir2 = -1;
				if (l.getNPieces() >= 2) {
					string p2 = l.findValue(1);
					if (p2 == "X" || p2 == "x") {
						dir2 = DERIV_X;
					}else if (p2 == "Y" || p2 == "y") {
						dir2 = DERIV_Y;
					}else {
						dir2 = DERIV_Z;
					}
				}

				// finally update infor
				derivInfor.updateDir(dir1,dir2);
			}

			// now push in 
			derivInforVec.push_back(derivInfor);
		}
	}
	
	// finally close the file
	inf.close();
}

Int DerivInforArray::getTotalNumDeriv() const
{
	Int n=0;
	Int nDeriv = getLenDerivInforArray();
	for(Int i=0; i<nDeriv; i++) {
		const DerivInfor& deriv = getDerivInfor(i);
		Int len = deriv.getDerivDirLen();
		n += len;
	}
	return n;
}

void DerivInforArray::print() const 
{
	cout << "job order is " << jobOrder << endl;
	cout << "LCode is     " << code << endl;
	Int nDeriv = getLenDerivInforArray();
	for(Int i=0; i<nDeriv; i++) {
		const DerivInfor& deriv = getDerivInfor(i);
		Int len = deriv.getDerivDirLen();
		if (jobOrder == 1) {
			cout << "position " << deriv.getDerivPos() << endl;
			cout << "dir len " << len << endl;
			for(Int j=0; j<len; j++) {
				cout << "dir " << deriv.getDerivDirection(j) << endl;
			}
		}else if (jobOrder == 2) {
			Int p1, p2;
			deriv.getDerivPos(p1,p2);
			cout << "position 1 " << p1 << " position 2 " << p2 << endl;
			for(Int j=0; j<len; j++) {
				Int d1, d2;
				deriv.getDerivDirection(j,d1,d2);
				cout << "dir 1 " << d1 << " dir 2 " << d2 << endl;
			}
		}
	}
}
