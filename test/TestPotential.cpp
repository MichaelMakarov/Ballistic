#include <iostream>
#include <iomanip>
#include <filesystem>
#include "Atmosphere1981.h"
#include "Atmosphere2004.h"
#include "GeoPotential.h"
#include "EGM96.h"
#include "Structures.h"
#include "DataProviders.h"
#include "SolarSystem.h"
#include <Conversions.h>

using namespace general::math;
using namespace general::time;
using namespace ball;

Vec3 atm_deceleration(const Vec3& vel, const double sb, const double density)
{
	return -density * vel.length() * sb * vel;
}


int main()
{
	std::cout << std::setprecision(16);
	std::cout << "Comparison of accelerations caused by geopotential and atmosphere...\n";

	auto path = std::filesystem::current_path().parent_path().parent_path().parent_path();
	std::string celestrackfile = path.string() + "/resources/celestrack/spaceweather.txt";
	auto spaceweather = load_spaceweather_data(celestrackfile);

	auto gm = EGM96("resources/gravity/egm96.txt");

	auto x{ State(
		5480653.6, 3721134.68, 1460430.37, 
		-3526.85815, 2909.53662, 5781.30371, 
		0.0204, 
		JD(datetime_from_str("2018-06-07 12:27:10.000")), 
		1) 
	};
	const Vec3 pos{ {x.vec[0], x.vec[1], x.vec[2]} };
	const Vec3 vel{ {x.vec[3], x.vec[4], x.vec[5]} };
	const auto& data = spaceweather[x.T.to_datetime().get_date()];
	auto [sunsph, sunort] = Sun::positionACS(x.T);
	sunort = ACS_to_GCS(sunort, sidereal_time_true(x.T));

	for (size_t i = 0; i < gm.count(); ++i) {
		auto gpt = GeoPotential(EGM96::Mu(), EGM96::R(), gm, 2 + i);
		std::cout << "gpt harmonics " << 2 + i << ": " << gpt.derivatives(pos) << std::endl;
	}

	auto atm84 = Atmosphere1981(EGM96::R(), EGM96::Fl());
	auto atm04 = Atmosphere2004(EGM96::R(), EGM96::Fl(), data.F10_7, data.F81, data.Kpsum / 8);

	std::cout << "stat atm: " << atm_deceleration(vel, x.Sb, atm84.density(pos, x.T)) << std::endl;
	std::cout << "dynm atm: " << atm_deceleration(vel, x.Sb, atm04.density(pos, x.T, sunort, sunsph[1])) << std::endl;
	std::cout << "sun: " << Sun::acceleration(pos, x.T) << std::endl;
	std::cout << "moon: " << Moon::acceleration(pos, x.T) << std::endl;

	return 0;
}