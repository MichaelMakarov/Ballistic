#include "GeoPotential.h"
#include "EGM96.h"
#include "GeneralConstants.h"
#include "Conversions.h"
#include "StaticAtmosphere.h"
#include <iostream>
#include <iomanip>

using namespace ball::space;
using namespace ball::geometry;

void test_conversions();
void test_geopotential();
void test_atmosphere();

int main()
{
	test_conversions();
	test_geopotential();
	test_atmosphere();
	return 0;
}

void test_conversions()
{
	std::cout << "\n...Test conversions...\n";

	auto xyzPosition{ XYZ(-4688980.289, -11060428.914, 238914.750) };
	std::cout << "Initial position in orthogonal: " << xyzPosition << std::endl;
	auto rblPosition = GCS_ortho_to_spher(xyzPosition);
	std::cout << "Converted to spherical: " << rblPosition << std::endl;
	auto xyzCheck = GCS_spher_to_ortho(rblPosition);
	std::cout << "Converted to orthogonal: " << xyzCheck << std::endl;

}

void test_geopotential()
{
	std::cout << "\n...Test geopotential...\n";

	std::cout << std::setprecision(16);
	auto gp{ GeoPotential(std::move(std::unique_ptr<EGM96>{ new EGM96() }), 16) };
	double delta = ball::math::deg_to_rad(5);
	double latitude{ 0 };
	double longitude{ ball::math::deg_to_rad(349.45) };
	double radius{ 6378137.000 + 5.63755000000E+06 };
	for (size_t i = 0; i < 60; ++i)
	{
		auto position{ RBL(radius, latitude, longitude) };
		std::cout << "The point is " << position << std::endl;
		std::cout << "The potential is " << gp(position) << std::endl;
		latitude += 0.1 * delta;
		longitude += delta;
	}
}

void test_atmosphere()
{
	std::cout << "\n...Test atmosphere...\n";

	EGM96 gm;
	auto atmosphere{ StaticAtmosphere81(gm.R(), gm.Fl()) };
	double delta = 1000.0;
	double height{ 10 };
	auto rblPosition{ RBL(gm.R() + height, 0, ball::math::deg_to_rad(349.45)) };
	std::cout << "Position: " << rblPosition << std::endl;
	for (size_t i = 0; i < 60; ++i)
	{
		std::cout << "h = " << height << "; rho = " << 
			atmosphere.density(GCS_spher_to_ortho(rblPosition), ball::time::JD2000) << std::endl;
		height += delta;
		rblPosition.R += delta;
	}
}