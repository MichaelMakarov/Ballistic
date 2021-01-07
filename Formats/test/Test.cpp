#include <iostream>
#include "TleFormat.h"
#include <PZ90.h>

void TestTleFormat();

int main()
{
	std::cout << "...Testing TleFormat...\n";
	TestTleFormat();
	return 0;
}

void TestTleFormat()
{
	const char* filepath = "E:/Users/projects/TleFile.tle";
	ball::formats::TleFormat tle1;
	if (!ball::formats::TleFromFile(filepath, tle1))
		std::cout << "Failed to load tle from file!\n";
	ball::types::OsculParams op = tle1.ToOsculParams(std::make_unique<ball::tasks::PZ90>()->Mu());
}