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
	Vec6 StatAtmModel::function(const Vec6& vec, const general::time::JD& t)
	{
		double w2{ _eW * _eW };
		auto pos = general::math::Vec3{ { vec[0], vec[1], vec[2] } };
		double h = GCS_height_from_position(pos, _eR, _eFl);
		if (h > this->MinHeight && h < this->MaxHeight) {
			// geopotential aceleration with a centrifugal and a coriolis force
			auto acpot{ _geopotential.derivatives(pos) };
			double density{ _atmosphere.density(pos, t) };
			// atmosphere deceleration a = v * s * rho, 
			// s - a ballistic coefficient,
			// v - a velocity of the vehicle,
			// rho - a rho of the atmosphere
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

	Vec6 DynAtmModel::function(const Vec6& vec, const general::time::JD& t)
	{
		double w2{ _eW * _eW };
		auto pos = general::math::Vec3{ { vec[0], vec[1], vec[2] } };
		double h = GCS_height_from_position(pos, _eR, _eFl);
		if (h > this->MinHeight && h < this->MaxHeight) {
			auto [sunsph, sunort] = Sun::positionACS(t);
			sunort = ACS_to_GCS(sunort, sidereal_time_avr(t));
			double density = _atmosphere.density(pos, t, sunort, sunsph[1]);
			// deceleration by atmosphere
			double acatm = density * this->sBall *
				std::sqrt(vec[3] * vec[3] + vec[4] * vec[4] + vec[5] * vec[5]);
			// acelerations by geopotential, sun and moon
			auto acpot = _geopotential.derivatives(pos);
			auto acsun = acceleration_by_masspoint(pos, sunort, Sun::Mu());
			auto acmoon = Moon::acceleration(pos, t);
			acpot += acsun + acmoon;
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

	DynAtmModel::EV DynAtmModel::extfunction(const DynAtmModel::EV& vec, const general::time::JD& t)
	{
		auto pos = slice<0, 2>(vec);
		double h = GCS_height_from_position(pos, _eR, _eFl);
		if (h > this->MinHeight && h < this->MaxHeight) {
			const double w2{ _eW * _eW };
			const double vel{ std::sqrt(vec[3] * vec[3] + vec[4] * vec[4] + vec[5] * vec[5]) };
			const double _vel2{ 1 / vel / vel };
			const double sidt = sidereal_time_avr(t);
			// solar position
			auto [sunsph, sunort] = Sun::positionACS(t);
			sunort = ACS_to_GCS(sunort, sidt);
			// lunar position
			auto [_, moonort] = Moon::positionACS(t);
			moonort = ACS_to_GCS(moonort, sidt);
			// density of atmosphere
			const double rho = _atmosphere.density(pos, t, sunort, sunsph[1]) * vel;
			// deceleration by atmosphere
			const double atmv = rho * this->sBall * vel;
			// acelerations by geopotential, sun and moon
			auto [potv, potm] = _geopotential.fullderivatives(pos);
			potv += acceleration_by_masspoint(pos, sunort, Sun::Mu());
			potv += acceleration_by_masspoint(pos, moonort, Moon::Mu());
			auto dacc = accelerationdiff_by_masspoint(pos, sunort, Sun::Mu());
			dacc += accelerationdiff_by_masspoint(pos, moonort, Moon::Mu());
			potm(0, 0) += dacc[0]; potm(1, 1) += dacc[1]; potm(2, 2) += dacc[2];
			// the addition all the components
			return EV{
				{
					// equations for vector of state
					vec[3], vec[4], vec[5],
					potv[0] + w2 * vec[0] + 2 * _eW * vec[4] - atmv * vec[3],
					potv[1] + w2 * vec[1] - 2 * _eW * vec[3] - atmv * vec[4],
					potv[2] - atmv * vec[5],
					// equations for variations
					vec[9], vec[10], vec[11],
					vec[6] * (w2 + potm(0, 0)) + vec[7] * potm(0, 1) + vec[8] * potm(0, 2) + 2 * _eW * vec[10] -
						atmv * (vec[9] + vec[3] * (vec[9] * vec[3] + vec[10] * vec[4] + vec[11] * vec[5]) / _vel2) - rho * vec[3] * vec[12],
					vec[6] * potm(0, 0) + vec[7] * (w2 + potm(0, 1)) + vec[8] * potm(0, 2) - 2 * _eW * vec[9] -
						atmv * (vec[10] + vec[4] * (vec[9] * vec[3] + vec[10] * vec[4] + vec[11] * vec[5]) / _vel2) - rho * vec[4] * vec[12],
					vec[6] * (w2 + potm(0, 0)) + vec[7] * potm(0, 1) + vec[8] * potm(0, 2) -
						atmv * (vec[11] + vec[5] * (vec[9] * vec[3] + vec[10] * vec[4] + vec[11] * vec[5]) / _vel2) - rho * vec[5] * vec[12],
					0.0
				}
			};
		}
		throw std::runtime_error("Height is out of bounds!");
	}

	Vec6 VarBallModel::function(const Vec6& vec, const general::time::JD& t)
	{
		double w2{ _eW * _eW };
		auto pos = general::math::Vec3{ { vec[0], vec[1], vec[2] } };
		double h = GCS_height_from_position(pos, _eR, _eFl);
		if (h > this->MinHeight && h < this->MaxHeight) {
			auto [sunsph, sunort] = Sun::positionACS(t);
			sunort = ACS_to_GCS(sunort, sidereal_time_avr(t));
			double density = _atmosphere.density(pos, t, sunort, sunsph[1]);
			// deceleration by atmosphere
			double acatm = density * ball_coefficient(t) *
				std::sqrt(vec[3] * vec[3] + vec[4] * vec[4] + vec[5] * vec[5]);
			// acelerations by geopotential, sun and moon
			auto acpot = _geopotential.derivatives(pos);
			auto acsun = acceleration_by_masspoint(pos, sunort, Sun::Mu());
			auto acmoon = Moon::acceleration(pos, t);
			acpot += acsun + acmoon;
			// the addition all the components
			return Vec6{
				{ vec[3], vec[4], vec[5],
				acpot[0] + w2 * vec[0] + 2 * _eW * vec[4] - acatm * vec[3],
				acpot[1] + w2 * vec[1] - 2 * _eW * vec[3] - acatm * vec[4],
				acpot[2] - acatm * vec[5] }
			};
		}
		throw std::runtime_error("Height is out of bounds!");
	}

	VarBallModel::EV VarBallModel::extfunction(const VarBallModel::EV& vec, const general::time::JD& t)
	{
		double w2{ _eW * _eW };
		auto pos = general::math::Vec3{ { vec[0], vec[1], vec[2] } };
		double h = GCS_height_from_position(pos, _eR, _eFl);
		if (h > this->MinHeight && h < this->MaxHeight) {
			const double sidt = sidereal_time_avr(t);
			// solar position
			auto [sunsph, sunort] = Sun::positionACS(t);
			sunort = ACS_to_GCS(sunort, sidt);
			// lunar position
			auto [_, moonort] = Moon::positionACS(t);
			moonort = ACS_to_GCS(moonort, sidt);
			// density of atmosphere
			double rho = _atmosphere.density(pos, t, sunort, sunsph[1]);
			rho *= ball_coefficient(t);
			const double vel{ std::sqrt(vec[3] * vec[3] + vec[4] * vec[4] + vec[5] * vec[5]) };
			// deceleration by atmosphere
			double atmv = rho * vel;
			// acelerations by geopotential, sun and moon
			auto [potv, potm] = _geopotential.fullderivatives(pos);
			potv += acceleration_by_masspoint(pos, sunort, Sun::Mu());
			potv += acceleration_by_masspoint(pos, moonort, Moon::Mu());
			auto dacc = accelerationdiff_by_masspoint(pos, sunort, Sun::Mu());
			dacc += accelerationdiff_by_masspoint(pos, moonort, Moon::Mu());
			potm(0, 0) += dacc[0]; potm(1, 1) += dacc[1]; potm(2, 2) += dacc[2];
			// the addition all the components
			return EV{
				{ 
					// equations for vector of state
					vec[3], vec[4], vec[5],
					potv[0] + w2 * vec[0] + 2 * _eW * vec[4] - atmv * vec[3],
					potv[1] + w2 * vec[1] - 2 * _eW * vec[3] - atmv * vec[4],
					potv[2] - atmv * vec[5],
					// equations for variations
					vec[9], vec[10], vec[11],
					vec[6] * (w2 + potm(0, 0)) + vec[7] * potm(0, 1) + vec[8] * potm(0, 2) - 
					vec[9] * atmv + 0
				}
			};
		}
		throw std::runtime_error("Height is out of bounds!");
	}
}