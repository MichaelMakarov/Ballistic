#pragma once
#include "Forecast.h"
#include "gravity/GeoPotential.h"
#include "atmosphere/Atmosphere.h"

namespace ball
{
	using namespace space;
	using namespace geometry;
	using namespace time;

	// A forecast implements static atmosphere 81 standard and PZ90 geopotential model
	class PzForecast : public Forecast<PzForecast>
	{
	private:
		GeoPotential _geopotential;
		std::shared_ptr<IAtmosphere> _pAtmosphere;
		double _eW, _eFl, _eR;

	public:
		PzForecast(
			const std::shared_ptr<IGravity>& pGravity,
			const size_t harmonics,
			const std::shared_ptr<IAtmosphere>& pAtmosphere) :
			_pAtmosphere{ pAtmosphere },
			_geopotential{ GeoPotential(pGravity, harmonics) },
			_eW{ pGravity->W() },
			_eFl{ pGravity->Fl() },
			_eR{ pGravity->R() }
		{}
		PzForecast(const PzForecast& f) = default;
		PzForecast(PzForecast&& f) noexcept : Forecast(f),
			_geopotential{ std::move(f._geopotential) },
			_pAtmosphere{ std::move(f._pAtmosphere) },
			_eW{ f._eW }, _eFl{ f._eFl }, _eR{ f._eR }
		{}
		~PzForecast() = default;

		PzForecast& operator = (const PzForecast& f) = default;
		PzForecast& operator = (PzForecast&& f) noexcept
		{
			_geopotential = std::move(f._geopotential);
			_pAtmosphere = std::move(f._pAtmosphere);
			_eW = f._eW;
			_eFl = f._eFl;
			_eR = f._eR;
			_eW = _eR = _eFl = 0;
			return *this;
		}

		// Acelerations calculation using current vector in GCS and time
		PV function(const PV& vec, const JD& t);
	};
}