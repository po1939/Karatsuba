#include <iostream>
#include "posint.h"
using namespace std;

int main() {

  // Set base to 16^3
  // This causes numbers to display in hex.  
  // However, base must be no greater than 2^15 when
  // using 32-bit integers, so we can only fit 3 hex 
  // digits (or 12 bits) per 32-bit int.  It is more
  // space-efficient to use the default, but then the
  // numbers will print in binary.
  
  PosInt::setBase(16,3);
  
  // PosInt::setBase(10, 4);

  // x = 2^128
	double time1;
	double time2;
	int threshold[3];
	
	PosInt p(0);
	PosInt increment(1);
	PosInt a(2147483647);
        PosInt b(1312311235);
	PosInt c(2147483647);
        PosInt d(1312311235);
	bool flag =0;
	
	for (int j =0;j<3;j++){
		p.set(0);
		flag =0;
		for (int i = 0;i<1000;i++){
			p.add(increment);		
	
			a.set(2147483647);
        		b.set(1312311235);
        		a.pow(p);
			b.pow(p);
		
        		int start_s=clock();
        		a.fastMul(b);
        		int stop_s=clock();
        		time1 = (stop_s-start_s)/double(CLOCKS_PER_SEC)*1000;

        		c.set(2147483647);
        		d.set(1312311235);
			c.pow(p);
			d.pow(p);       

        		start_s=clock();
        		c.mul(d);
        		stop_s=clock();
        		time2 = (stop_s-start_s)/double(CLOCKS_PER_SEC)*1000;
 
			if ((i+1)%100 == 0){
				cout <<"p is "<<i+1<<". fastMul time: "<<time1<<"  mul time: "<<time2<<endl;
			}
			if (flag == 0 && time1 < time2) {
				cout <<"\n--------fastMul is faster. p is: ";
				p.print_array(cout);
				cout <<"--------\n";
				flag =1;
				threshold[j] = i+1;
			}
		}
	}
	int tTotal=0;
	int tAvg=0;
	for (int i =0;i<3;i++){
		tTotal += threshold[i];
	}
	tAvg = tTotal/3;
	cout << "\nThreshold average is: "<<tAvg<<endl;



 // z = random number between 0 and 2^128 - 1
  //  srand(0);
 
/*
  srand(time(NULL));
  PosInt z;
  z.rand(x);

  cout << "x = 2^128 = " << x << endl;
  cout << "x         = ";
  x.print_array(cout);
  cout << endl << endl;

  cout << "z = rand(x) = " << z << endl;
  cout << "z           = ";
  z.print_array(cout);
  cout << endl << endl;

  // z = z^2
  z.mul(z);

  cout << "z^2 = " << z << endl;

  // z = z mod x
  z.mod(x);

  cout << "z^2 mod x = " << z << endl;
*/
}
  
