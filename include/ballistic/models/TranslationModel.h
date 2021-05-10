#pragma once
#include "general/Geometry.h"
#include "general/Times.h"

namespace ball
{
	// A forecast interface provides a function that returns the vector of size 6.
	// ForecastType is a type of the class that inherits the interface (implements static polymorphism)
	template<class Model>
	class TranslationModel
	{
	public:
		TranslationModel() = default;
		~TranslationModel() = default;
		// Acelerations calculation using current vector in GCS and time
		general::math::PV function(
			const general::math::PV& vec,
			const general::time::JD& t)
		{
			return static_cast<Model*>(this)->function(vec, t);
		}

	public:
		double MinHeight = 1e4, MaxHeight = 1e8;
		double sBall = 0.0;
	};
}