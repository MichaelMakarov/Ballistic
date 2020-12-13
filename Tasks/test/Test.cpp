#include "GeoPotential.h"
#include "EGM96.h"
#include "Constants.h"
#include "AstroValues.h"
#include "MathFunctions.h"
#include <iostream>
#include <iomanip>

using namespace ball::tasks;
using namespace ball;

void TestConversions();
void TestGeoPotential();

int main()
{
	std::cout << "...Conversion test...\n";
	TestConversions();
	std::cout << "...GeoPotential tests...\n";
	TestGeoPotential();
	return 0;
}

void TestConversions()
{
	geometry::XYZ p(-4688980.289, -11060428.914, 238914.750);
	std::cout << "Initial position in orthogonal: " << p << std::endl;
	auto rbl = GCS_OrthoToSpher(p);
	std::cout << "Converted to spherical: " << rbl << std::endl;
	auto xyz = GCS_SpherToOrtho(rbl);
	std::cout << "Converted to orthogonal: " << xyz << std::endl;

}

void TestGeoPotential()
{
	GeoPotential gp(
		std::move(std::unique_ptr<PotentialEGM96>{ new PotentialEGM96() }),
		6378137,
		3.98600441800E+14,
		/*egm96::R,
		egm96::Mu,*/
		16);
	geometry::RBL p(6378137.000 + 5.63755000000E+06, 0.0, math::DegToRad(349.45));
	double a = gp(p);// / c.R;
	std::cout << std::setprecision(16);
	std::cout << "The point is " << p << std::endl;
	std::cout << "The potential is " << a << std::endl;
}