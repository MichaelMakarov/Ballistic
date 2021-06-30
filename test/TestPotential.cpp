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

void test_siderealtime() {
	constexpr auto jd = JD(DateTime(2021, 10, 29, 17, 6, 18));
	constexpr double meansid = sidereal_time_mean(jd);
	const double truesid = sidereal_time_true(jd);
	std::cout << "mean sidereal time: " << decimaltime_to_hms(meansid * 12 / PI) << std::endl;
	std::cout << "true sidereal time: " << decimaltime_to_hms(truesid * 12 / PI) << std::endl;
	std::cout << "correct time value: 19:39:17.2599\n\n";
}

int main()
{
	test_siderealtime();

	std::cout << std::setprecision(16);
	std::cout << "Comparison of accelerations caused by geopotential and atmosphere...\n";

	auto path = std::filesystem::current_path().parent_path().parent_path().parent_path();
	std::string celestrackfile = path.string() + "/resources/celestrack/spaceweather.txt";
	auto spaceweather = load_spaceweather_data(celestrackfile);

	auto gm = EGM96("resources/gravity/egm96.txt");

	Params<> x{
		{
			5480653.6, 3721134.68, 1460430.37,
			-3526.85815, 2909.53662, 5781.30371 
		},
		JD(datetime_from_str("2018-06-07 12:27:10.000")), 
		0.0204,
		1
	};

	const Vec3 pos = slice<0, 2>(x.vec);
	const Vec3 vel = slice<3, 5>(x.vec);
	const auto& data = spaceweather[jd_to_datetime(x.T).get_date()];
	auto [sunsph, sunort] = Sun::positionACS(x.T);
	sunort = ACS_to_GCS(sunort, sidereal_time_true(x.T));

	for (size_t i = 0; i < gm.count(); ++i) {
		auto gpt = GeoPotential(EGM96::Mu(), EGM96::R(), gm, 2 + i);
		std::cout << "gpt harmonics " << 2 + i << ": " << gpt.acceleration(pos) << std::endl;
	}

	auto atm84 = Atmosphere1981(EGM96::R(), EGM96::Fl());
	auto atm04 = Atmosphere2004(EGM96::R(), EGM96::Fl(), data.F10_7, data.F81, data.Kpsum / 8);

	std::cout << "stat atm: " << atm_deceleration(vel, x.vec[6], atm84.density(pos, x.T)) << std::endl;
	std::cout << "dynm atm: " << atm_deceleration(vel, x.vec[6], atm04.density(pos, x.T, sunort, sunsph[1])) << std::endl;
	std::cout << "sun: " << Sun::acceleration(pos, x.T) << std::endl;
	std::cout << "moon: " << Moon::acceleration(pos, x.T) << std::endl;

	return 0;
}