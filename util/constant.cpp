//
// This simple program is used to generate the constants 
// for the cppints program
// all of constants are defined as macro
// fenglai liu
//

#include<cstdlib>
#include<cstdio>
#include<string>
#include<cmath>
#include "boost/lexical_cast.hpp"
using boost::lexical_cast;
using namespace std;
int main() 
{
	cout << "/*******************************************************" << endl;
	cout << " *THIS FILE PROVIDES CONSTANTS FOR THE CPPINTS PROGRAM *" << endl;
	cout << " *ALL CONSTANTS HAVE 18 DECIMALS ACCURACY              *" << endl;
	cout << " *******************************************************/" << endl;
	cout << endl << endl << endl;

	// generate the number of 2m+1
	cout << "/*******************************************************" << endl;
	cout << " *          CONSTANTS 1/(2M+1) M FROM 0 TO 100         *" << endl;
	cout << " *******************************************************/" << endl;
	string mainpart = "ONEOVER";
	for(int i=0; i<=100; i++) {
		long double r      = i; 
		long double result = 1.0E0/(2.0E0*r+1.0E0);
		string n = lexical_cast<string>(2*i+1);
		string name = mainpart + n;
		printf("#define   %-20s   %-20.18Lf\n", name.c_str(), result);
	}
	cout << endl << endl << endl;

	// generate other constants
	cout << "/*******************************************************" << endl;
	cout << " *                    OTHER CONSTANTS                  *" << endl;
	cout << " *******************************************************/" << endl;
	long double c  = 1.0;
	long double pi = 4.0*atan(c);
	printf("#ifndef PI\n");
	printf("#define   %-20s   %-20.18Lf\n", "PI", pi);
	printf("#endif\n");
	long double result = 2.0E0/(sqrt(pi));
	printf("#define   %-20s   %-20.18Lf\n", "TWOOVERSQRTPI", result);

	return 0;
}
