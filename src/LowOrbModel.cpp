#include "..\include\ballistic\models\LowOrbModel.h"
#include "..\include\ballistic\models\LowOrbModel.h"
#include "..\include\ballistic\models\LowOrbModel.h"
#include "..\include\ballistic\models\LowOrbModel.h"
#include "..\include\ballistic\models\LowOrbModel.h"
#include "..\include\ballistic\models\LowOrbModel.h"
#include "..\include\ballistic\models\LowOrbModel.h"
#include "..\include\ballistic\models\LowOrbModel.h"
#include "..\include\ballistic\models\LowOrbModel.h"
#include "..\include\ballistic\models\LowOrbModel.h"
#include "LowOrbModel.h"

namespace ball
{
	Vec6 MSA::function(const Vec6& vec, const general::time::JD& t)
	{
		double w2{ _eW * _eW };
		auto pos = general::math::Vec3{ { vec[0], vec[1], vec[2] } };
		double h = height_from_gcsposition(pos, _eR, _eFl);
		if (h > this->MinHeight && h < this->MaxHeight) {
			// geopotential aceleration with a centrifugal and a coriolis force
			auto acpot{ _geopotential.acceleration(pos) };
			double density{ _atmosphere.density(pos, t) };
			// atmosphere deceleration a = v * s * density, 
			// s - a ballistic coefficient,
			// v - a velocity of the vehicle,
			// density - a density of the atmosphere
			double acatm = density * this->sBall *
				std::sqrt(vec[3] * vec[3] + vec[4] * vec[4] + vec[5] * vec[5]);
			// the addition all the components
			return Vec6{
				{
					vec[3], vec[4], vec[5],
					acpot[0] + w2 * vec[0] + 2 * _eW * vec[4] - acatm * vec[3],
					acpot[1] + w2 * vec[1] - 2 * _eW * vec[3] - acatm * vec[4],
					acpot[2] - acatm * vec[5] 
				}
			};
		}
		throw std::runtime_error("Height is out of bounds!");
	}

	Vec6 MDASM::function(const Vec6& vec, const general::time::JD& t)
	{
		double w2{ _eW * _eW };
		auto pos = general::math::Vec3{ vec[0], vec[1], vec[2] };
		double h = height_from_gcsposition(pos, _eR, _eFl);
		if (h > this->MinHeight && h < this->MaxHeight) {
			auto [sunsph, sunort] = Sun::positionACS(t);
			sunort = ACS_to_GCS(sunort, sidereal_time_mean(t));
			double density = _atmosphere.density(pos, t, sunort, sunsph[1]);
			// deceleration by atmosphere
			double atmv = density * this->sBall *
				std::sqrt(vec[3] * vec[3] + vec[4] * vec[4] + vec[5] * vec[5]);
			// acelerations by geopotential, sun and moon
			auto potv = _geopotential.acceleration(pos);
			//auto pots = acceleration_by_masspoint(pos, sunort, Sun::Mu());
			//auto potm = Moon::acceleration(pos, t);
			//potv += pots + potm;
			// the addition all the components
			return Vec6
			{
				vec[3], vec[4], vec[5],
				potv[0] + w2 * vec[0] + 2 * _eW * vec[4] - atmv * vec[3],
				potv[1] + w2 * vec[1] - 2 * _eW * vec[3] - atmv * vec[4],
				potv[2] - atmv * vec[5]
			};
		}
		throw std::runtime_error("Height is out of bounds!");
	}

	general::math::Vec<55> MDASM::extfunction(const general::math::Vec<55>& vec, const general::time::JD& t)
	{
		auto pos = slice<0, 2>(vec);
		double h = height_from_gcsposition(pos, _eR, _eFl);
		if (h > this->MinHeight && h < this->MaxHeight) {
			const double w2{ _eW * _eW };
			const double vel2{ vec[3] * vec[3] + vec[4] * vec[4] + vec[5] * vec[5] };
			const double vel = std::sqrt(vel2);
			const double sidt = sidereal_time_true(t);
			// solar position
			auto [sunsph, sunort] = Sun::positionACS(t);
			sunort = ACS_to_GCS(sunort, sidt);
			// lunar position
			auto [_, moonort] = Moon::positionACS(t);
			moonort = ACS_to_GCS(moonort, sidt);
			// density of atmosphere
			const double rhov = _atmosphere.density(pos, t, sunort, sunsph[1]) * vel;
			// deceleration by atmosphere
			const double atmv = rhov * this->sBall;
			// acelerations by geopotential, sun and moon
			auto [potv, potm] = _geopotential.derivatives(pos);
			potv += acceleration_by_masspoint(pos, sunort, Sun::Mu());
			potv += acceleration_by_masspoint(pos, moonort, Moon::Mu());
			auto dacc = accelerationdiff_by_masspoint(pos, sunort, Sun::Mu());
			dacc += accelerationdiff_by_masspoint(pos, moonort, Moon::Mu());
			potm(0, 0) += dacc[0]; potm(1, 1) += dacc[1]; potm(2, 2) += dacc[2];
			// the addition all the components
			general::math::Vec<55> res{
				vec[3], vec[4], vec[5],
				potv[0] + w2 * vec[0] + 2 * _eW * vec[4] - atmv * vec[3],
				potv[1] + w2 * vec[1] - 2 * _eW * vec[3] - atmv * vec[4],
				potv[2] - atmv * vec[5]
			};
			size_t i{ 6 };
			//equation of variations
			for (size_t k = 0; k < 7; ++k) {
				res[i + 0] = vec[i + 3];
				res[i + 1] = vec[i + 4];
				res[i + 2] = vec[i + 5];
				res[i + 3] = vec[i] * (w2 + potm(0, 0)) + vec[i + 1] * potm(0, 1) + vec[i + 2] * potm(0, 2) + 2 * _eW * vec[i + 4] -
					atmv * (vec[i + 3] + vec[3] * (vec[i + 3] * vec[3] + vec[i + 4] * vec[4] + vec[i + 5] * vec[5]) / vel2) - rhov * vec[3] * vec[i + 6];
				res[i + 4] = vec[i] * potm(1, 0) + vec[i + 1] * (w2 + potm(1, 1)) + vec[i + 2] * potm(1, 2) - 2 * _eW * vec[i + 3] -
					atmv * (vec[i + 4] + vec[4] * (vec[i + 3] * vec[3] + vec[i + 4] * vec[4] + vec[i + 5] * vec[5]) / vel2) - rhov * vec[4] * vec[i + 6];
				res[i + 5] = vec[i] * potm(2, 0) + vec[i + 1] * potm(2, 1) + vec[i + 2] * (w2 + potm(2, 2)) -
					atmv * (vec[i + 5] + vec[5] * (vec[i + 3] * vec[3] + vec[i + 4] * vec[4] + vec[i + 5] * vec[5]) / vel2) - rhov * vec[5] * vec[i + 6];
				res[i + 6] = 0.0;
				i += 7;
			}
			return res;
			//return Vec13{
			//	// equations for vector of state
			//	vec[3], vec[4], vec[5],
			//	potv[0] + w2 * vec[0] + 2 * _eW * vec[4] - atmv * vec[3],
			//	potv[1] + w2 * vec[1] - 2 * _eW * vec[3] - atmv * vec[4],
			//	potv[2] - atmv * vec[5],
			//	// equations for variations
			//	vec[9], vec[10], vec[11],
			//	vec[6] * (w2 + potm(0, 0)) + vec[7] * potm(0, 1) + vec[8] * potm(0, 2) + 2 * _eW * vec[10] -
			//		atmv * (vec[9] + vec[3] * (vec[9] * vec[3] + vec[10] * vec[4] + vec[11] * vec[5]) / vel2) - rhov * vec[3] * vec[12],
			//	vec[6] * potm(0, 0) + vec[7] * (w2 + potm(0, 1)) + vec[8] * potm(0, 2) - 2 * _eW * vec[9] -
			//		atmv * (vec[10] + vec[4] * (vec[9] * vec[3] + vec[10] * vec[4] + vec[11] * vec[5]) / vel2) - rhov * vec[4] * vec[12],
			//	vec[6] * potm(0, 0) + vec[7] * potm(0, 1) + vec[8] * (w2 + potm(0, 2)) -
			//		atmv * (vec[11] + vec[5] * (vec[9] * vec[3] + vec[10] * vec[4] + vec[11] * vec[5]) / vel2) - rhov * vec[5] * vec[12],
			//	0.0
			//};
		}
		throw std::runtime_error("Height is out of bounds!");
	}

	Vec6 MDASMVb::function(const Vec6& vec, const general::time::JD& t)
	{
		using namespace general::math;
		const double w2{ _eW * _eW };
		const Vec3 pos{ vec[0], vec[1], vec[2] };
		const Vec3 vel{ vec[3], vec[4], vec[5] };
		double h = height_from_gcsposition(pos, _eR, _eFl);
		if (h > this->MinHeight && h < this->MaxHeight) {
			auto [sunsph, sunort] = Sun::positionACS(t);
			sunort = ACS_to_GCS(sunort, sidereal_time_mean(t));
			double density = _atmosphere.density(pos, t, sunort, sunsph[1]);
			// deceleration by atmosphere
			double acatm = density * ball_coefficient(pos, vel, t) * vel.length();
			// acelerations by geopotential, sun and moon
			auto acpot = _geopotential.acceleration(pos);
			auto acsun = acceleration_by_masspoint(pos, sunort, Sun::Mu());
			auto acmoon = Moon::acceleration(pos, t);
			acpot += acsun + acmoon;
			// the addition all the components
			return Vec6 {
				vec[3], vec[4], vec[5],
				acpot[0] + w2 * vec[0] + 2 * _eW * vec[4] - acatm * vec[3],
				acpot[1] + w2 * vec[1] - 2 * _eW * vec[3] - acatm * vec[4],
				acpot[2] - acatm * vec[5]
			};
		}
		throw std::runtime_error("Height is out of bounds!");
	}

	Vec16 MDASMVb::extfunction(const Vec16& vec, const general::time::JD& t)
	{
		using namespace general::math;
		double w2{ _eW * _eW };
		Vec3 pos{ vec[0], vec[1], vec[2] };
		Vec3 vel{ vec[3], vec[4], vec[5] };
		double h = height_from_gcsposition(pos, _eR, _eFl);
		if (h > this->MinHeight && h < this->MaxHeight) {
			const double sidt = sidereal_time_mean(t);
			// solar position
			auto [sunsph, sunort] = Sun::positionACS(t);
			sunort = ACS_to_GCS(sunort, sidt);
			// lunar position
			auto [_, moonort] = Moon::positionACS(t);
			moonort = ACS_to_GCS(moonort, sidt);
			// density of atmosphere
			double rho = _atmosphere.density(pos, t, sunort, sunsph[1]);
			rho *= ball_coefficient(pos, vel, t);
			const double v = vel.length();
			// deceleration by atmosphere
			double atmv = rho * v;
			// acelerations by geopotential, sun and moon
			auto [potv, potm] = _geopotential.derivatives(pos);
			potv += acceleration_by_masspoint(pos, sunort, Sun::Mu());
			potv += acceleration_by_masspoint(pos, moonort, Moon::Mu());
			auto dacc = accelerationdiff_by_masspoint(pos, sunort, Sun::Mu());
			dacc += accelerationdiff_by_masspoint(pos, moonort, Moon::Mu());
			potm(0, 0) += dacc[0]; potm(1, 1) += dacc[1]; potm(2, 2) += dacc[2];
			// the addition all the components
			return Vec16
			{
				// equations for vector of state
				vec[3], vec[4], vec[5],
				potv[0] + w2 * vec[0] + 2 * _eW * vec[4] - atmv * vec[3],
				potv[1] + w2 * vec[1] - 2 * _eW * vec[3] - atmv * vec[4],
				potv[2] - atmv * vec[5],
				0,//this->sBall,
				// equations for variations
				vec[9], vec[10], vec[11],
				vec[6] * (w2 + potm(0, 0)) + vec[7] * potm(0, 1) + vec[8] * potm(0, 2) -
				vec[9] * atmv + 0
			};
		}
		throw std::runtime_error("Height is out of bounds!");
	}
}