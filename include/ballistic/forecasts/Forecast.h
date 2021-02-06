#pragma once
#include "general/Geometry.h"
#include "general/Times.h"

namespace ball
{
	using namespace general::geometry;
	using namespace general::time;

	// A forecast interface provides a function for the aceleration calculation.
	// ForecastType is a type of the class inherits the interface is used for static polymorphism
	template<class ForecastType>
	class Forecast
	{
	public:
		Forecast() = default;
		Forecast(const Forecast& f) noexcept : sBall{ f.sBall }, MinHeight{ f.MinHeight }, MaxHeight{ f.MaxHeight } {}
		Forecast(Forecast&& f) noexcept : sBall{ f.sBall }, MinHeight{ f.MinHeight }, MaxHeight{ f.MaxHeight }
		{
			sBall = 0.0;
		}
		~Forecast() = default;

		Forecast& operator = (const Forecast& f) noexcept
		{
			sBall = f.sBall;
			MinHeight = f.MinHeight;
			MaxHeight = f.MaxHeight;
			return *this;
		}
		Forecast& operator = (Forecast&& f) noexcept
		{
			sBall = f.sBall;
			MinHeight = f.MinHeight;
			MaxHeight = f.MaxHeight;
			sBall = 0.0;
			return *this;
		}
		// Acelerations calculation using current vector in GCS and time
		PV function(const PV& vec, const JD& t)
		{
			return static_cast<ForecastType*>(this)->function(vec, t);
		}

	public:
		double MinHeight = 1e5, MaxHeight = 1e8;
		double sBall = 0.0;
	};
}