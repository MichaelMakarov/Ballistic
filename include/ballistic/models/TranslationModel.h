#pragma once
#include "Structures.h"

namespace ball
{
	// A forecast interface provides a function that returns the vector of dim 6.
	// ForecastType is a type of the class that inherits the interface (implements static polymorphism)
	template<class Model>
	class TranslationModel
	{
	public:
		// Acelerations calculation using current vector in GCS and time
		Vec6 function(const Vec6& vec, const general::time::JD& t)
		{
			return static_cast<Model*>(this)->function(vec, t);
		}

	public:
		double MinHeight = 1e4, MaxHeight = 1e8;
	};
}