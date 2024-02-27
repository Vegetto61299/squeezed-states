#include "class.h"
#include <fstream>
int main()
{
	ofstream file("data.txt");

	state s;
	s.squeeze(0,1,0);
	cout << real(s.scs);
	file.close();
	return 0;
}