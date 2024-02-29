#include "class.h"
#include <fstream>
int main()
{
	state s;
	std::thread one(&state::plotstate, s, 0, 1);
	std::thread two(&state::plotstatems, s, 0, 1);
	std::thread three(&state::plotstate, s, 1, 0);
	std::thread four(&state::plotstatems, s, 1, 0);
	std::thread five(&state::plotstate, s, 1, 1);
	std::thread six(&state::plotstatems, s, 1, 1);
	std::thread seven(&state::plotstate, s, 1, 2);
	std::thread eight(&state::plotstatems, s, 1, 2);
	one.join();
	two.join();
	three.join();
	four.join();
	five.join();
	six.join();
	seven.join();
	eight.join(); 
	return 0;
}