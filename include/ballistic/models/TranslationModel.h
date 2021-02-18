#pragma once
#include "general/Geometry.h"
#include "general/Times.h"

namespace ball
{
	// A forecast interface provides a function for the aceleration calculation.
	// ForecastType is a type of the class inherits the interface is used for static polymorphism
	template<class ModelType>
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
			return static_cast<ModelType*>(this)->function(vec, t);
		}

	public:
		double MinHeight = 1e5, MaxHeight = 1e8;
		double sBall = 0.0;
	};
}