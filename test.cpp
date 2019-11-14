#include "polygamma.hpp"
#include <iostream>

using namespace std;

int main() {
  cout << "-0.5 " << diGamma(complex<double>(-0.5, 0)) << endl;
	cout << "0.5 " << diGamma(complex<double>(0.5, 0)) << endl;
	cout << "0.5i " << diGamma(complex<double>(0, 0.5)) << endl;
	cout << "20 " << diGamma(complex<double>(20, 0)) << endl;
	
	cout << endl;
	
	cout << "-0.5 1 " << pGamma(1, complex<double>(-0.5, 0)) << endl;
	cout << "0.5 1 " << pGamma(1, complex<double>(0.5, 0)) << endl;
	cout << "0.5i 1 " << pGamma(1, complex<double>(0, 0.5)) << endl;
	cout << "18 1 " << pGamma(1, complex<double>(18,0)) << endl;
	cout << "20 1 " << pGamma(1, complex<double>(20,0)) << endl;
	
	cout << endl;
	
	cout << "-0.5 2 " << pGamma(2, complex<double>(-0.5, 0)) << endl;
	cout << "0.5 2 " << pGamma(2, complex<double>(0.5, 0)) << endl;
	cout << "0.5i 2 " << pGamma(2, complex<double>(0,0.5)) << endl;
	cout << "18 2 " << pGamma(2, complex<double>(18,0)) << endl;
	cout << "20 2 " << pGamma(2, complex<double>(20,0)) << endl;
}
