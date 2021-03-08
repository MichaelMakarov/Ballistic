#pragma once
#include "TranslationModel.h"
#include "gravity/GeoPotential.h"
#include "atmosphere/Atmosphere.h"
#include "Conversions.h"

namespace ball
{
	// A forecast implements atmosphere standard and geopotential model
	class PotAtmModel : public TranslationModel<PotAtmModel>
	{
	private:
		GeoPotential _geopotential;
		StaticAtmosphere81 _atmosphere;
		double _eW, _eFl, _eR;

	public:
		PotAtmModel(
			const IEarth& gravity,
			const size_t harmonics) :
			_atmosphere{ StaticAtmosphere81(gravity.R(), gravity.Fl()) },
			_geopotential{ GeoPotential(gravity, harmonics) },
			_eW{ gravity.W() },
			_eFl{ gravity.Fl() },
			_eR{ gravity.R() }
		{}
		~PotAtmModel() = default;

		// Acelerations calculation using current vector in GCS and time
		general::math::PV function(const general::math::PV& vec, const general::time::JD& t)
		{
			auto ac = general::math::PV(vec.Vel, general::math::Vec3());
			double r = vec.Pos.length();
			double r_2{ r * r };
			double w_2{ _eW * _eW };
			double h = GCS_height_from_position(vec.Pos, _eR, _eFl);
			if (h > this->MinHeight && h < this->MaxHeight) {
				// geopotential aceleration with a centrifugal and a coriolis force
				auto xyzAcPot{ _geopotential.acceleration(vec.Pos) };
				double density = _atmosphere.density(vec.Pos, t);
				double v = vec.Vel.length();
				double acAtm = v * density * this->sBall;
				// atmosphere aceleration a = v * s * rho, 
				// s - a ballistic coefficient,
				// v - a velocity of the vehicle,
				// rho - a density of the atmosphere
				// the addition all the components
				ac.Vel.X = xyzAcPot.X + w_2 * vec.Pos.X + 2 * _eW * vec.Vel.Y - acAtm * vec.Vel.X;
				ac.Vel.Y = xyzAcPot.Y + w_2 * vec.Pos.Y - 2 * _eW * vec.Vel.X - acAtm * vec.Vel.Y;
				ac.Vel.Z = xyzAcPot.Z - acAtm * vec.Vel.Z;
				return ac;
			}
			throw std::runtime_error("Height is out of bounds!");
		}
	};
}