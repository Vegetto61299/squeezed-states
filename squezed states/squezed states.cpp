#include "class.h"
#include <fstream>
int main()
{

	state s;
	s.plotstate(0,1);
	s.plotstatems(0,1);
	s.plotstate(1, 0);
	s.plotstatems(1, 0);
	s.plotstate(1, 1);
	s.plotstatems(1, 1);
	s.plotstate(1, 2);
	s.plotstatems(1, 2);
	return 0;
}