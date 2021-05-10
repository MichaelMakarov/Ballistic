#include <iostream>
#include <iomanip>
#include <filesystem>
#include "Conversions.h"
#include "Atmosphere1981.h"
#include "Atmosphere2004.h"
#include "SolarSystem.h"
#include "EGM96.h"
#include "general/Mathematics.h"
#include "DataProviders.h"

using namespace ball;
using namespace general::time;
using namespace general::math;

int main()
{
	std::cout << std::setprecision(15);

	auto path = std::filesystem::current_path().parent_path().parent_path().parent_path();
	std::string celestrackfile = path.string() + "/resources/celestrack/spaceweather.txt";
	auto spaceweather = load_spaceweather_data(celestrackfile);
	auto datetime{ DateTime(2021, 1, 20, 16, 33, 0) };
	auto jd{ JD(datetime) };
	const auto& data = spaceweather[datetime.get_date()];

	auto stat_atm{ Atmosphere1981(EGM96::R(), EGM96::Fl()) };
	auto dynm_atm{ Atmosphere2004(EGM96::R(), EGM96::Fl(), data.F10_7, data.F81, data.Kpsum / 8) };
	const double delta = 10e3;
	double height{ 10e3 };
	auto [sunsph, sunort] = Sun::positionACS(jd);
	sunort = ACS_to_GCS(sunort, sidereal_time_true(jd));
	auto sphpos{ Vec3(EGM96::R() + height, 0, deg_to_rad(349.45)) };
	std::cout << "F10,7 = " << data.F10_7 << "; F81 = " << data.F81 << "; " << data.Kpsum / 8 << std::endl;
	std::cout << "Position: " << sphpos << std::endl;
	for (size_t i = 0; i < 60; ++i) {
		auto pos{ sph_to_ort(sphpos) };
		std::cout << "h = " << height
			<< "; msa rho = " << stat_atm.density(pos, jd)
			<< "; mda rho = " << dynm_atm.density(pos, jd, sunort, sunsph.y())
			<< std::endl;
		height += delta;
		sphpos.x() += delta;
	}
	return 0;
}