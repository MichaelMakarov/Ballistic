#include <iostream>
#include "Conversions.h"
using namespace ball;

int main()
{
	using namespace general::time;
	using namespace general::math;

	std::cout << "\n...conversion between spherical and orthogonal coordinate systems...\n";

	auto xyz_pos{ general::math::Vec3(-4688980.289, -11060428.914, 238914.750) };
	std::cout << "point in orthogonal CS: " << xyz_pos << std::endl;
	auto rbl_pos = ort_to_sph(xyz_pos);
	std::cout << "point converted to spherical CS: " << rbl_pos << std::endl;
	auto xyz_back = sph_to_ort(rbl_pos);
	std::cout << "point converted back to orthogonal CS: " << xyz_back << std::endl;

	std::cout << "\n...calculating solar coordinates...\n";

	auto time{ JD(DateTime(2019, 1, 4, 0, 0, 0)) };
	int minutes{ 10 };
	for (size_t k{ 0 }; k < 100; ++k) {
		//ACS_solar_position()
	}

	return 0;
}