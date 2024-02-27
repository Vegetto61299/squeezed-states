#include "class.h"
#include <fstream>
int main()
{
	ofstream file("data.txt");

	squeeze state;
	file << real(state.phi);
	file.close();
}

