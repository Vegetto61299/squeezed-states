#include "class.h"
#include <fstream>
int main()
{
	ofstream file("data.txt");

	state s;
	s.squeeze(0,1,0);
	
	cout << real(s.scs);//% is elementwise multiplication c=a%b means c(0)=a(0)*b(0)... so it's modulo squared
	file.close();
	return 0;
}