#pragma once
#include "TranslationModel.h"
#include "gravity/GeoPotential.h"
#include "atmosphere/Atmosphere.h"
#include "Conversions.h"

namespace ball
{
	// A forecast implements atmosphere standard and geopotential model
	template<class AtmType>
	class PotAtmModel : public TranslationModel<PotAtmModel<AtmType>>
	{
	private:
		GeoPotential _geopotential;
		std::shared_ptr<IAtmosphere<AtmType>> _pAtmosphere;
		double _eW, _eFl, _eR;

	public:
		PotAtmModel(
			const std::shared_ptr<IEarth>& pGravity,
			const size_t harmonics,
			const std::shared_ptr<IAtmosphere<AtmType>>& pAtmosphere) :
			_pAtmosphere{ pAtmosphere },
			_geopotential{ GeoPotential(pGravity, harmonics) },
			_eW{ pGravity->W() },
			_eFl{ pGravity->Fl() },
			_eR{ pGravity->R() }
		{}
		/*PotAtmModel(const PotAtmModel& f) = default;
		PotAtmModel(PotAtmModel&& f) noexcept : TranslationModel<PotAtmModel<AtmType>>(f),
			_geopotential{ std::move(f._geopotential) },
			_pAtmosphere{ std::move(f._pAtmosphere) },
			_eW{ f._eW }, _eFl{ f._eFl }, _eR{ f._eR }
		{}*/
		~PotAtmModel() = default;

		/*PotAtmModel& operator = (const PotAtmModel& f) = default;
		PotAtmModel& operator = (PotAtmModel&& f) noexcept
		{
			_geopotential = std::move(f._geopotential);
			_pAtmosphere = std::move(f._pAtmosphere);
			_eW = f._eW;
			_eFl = f._eFl;
			_eR = f._eR;
			_eW = _eR = _eFl = 0;
			return *this;
		}*/

		// Acelerations calculation using current vector in GCS and time
		general::math::PV function(const general::math::PV& vec, const general::time::JD& t)
		{
			auto ac = general::math::PV(vec.Vel, general::math::Vec3());
			double r = vec.Pos.length();
			double r_2{ r * r };
			double w_2{ _eW * _eW };
			double h = GCS_height_from_position(vec.Pos, _eR, _eFl);
			if (h < this->MinHeight || h > this->MaxHeight)
				throw std::runtime_error("Height is out of bounds!");
			// geopotential aceleration with a centrifugal and a coriolis force
			auto xyzAcPot{ _geopotential.acceleration(vec.Pos) };
			// atmosphere aceleration a = v * s * rho, 
			// s - a ballistic coefficient,
			// v - a velocity of the vehicle,
			// rho - a density of the atmosphere
			double density = _pAtmosphere->density(vec.Pos, t);
			double v = vec.Vel.length();
			double acAtm = v * density * this->sBall;
			// the addition all the components
			ac.Vel.X = xyzAcPot.X + w_2 * vec.Pos.X + 2 * _eW * vec.Vel.Y - acAtm * vec.Vel.X;
			ac.Vel.Y = xyzAcPot.Y + w_2 * vec.Pos.Y - 2 * _eW * vec.Vel.X - acAtm * vec.Vel.Y;
			ac.Vel.Z = xyzAcPot.Z - acAtm * vec.Vel.Z;
			return ac;
		}
	};
}