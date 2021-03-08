#pragma once
#include "general/Quaternion.h"
#include "general/Times.h"

namespace ball
{
	// A forecast implements solar radiation pressure
	template<class Model>
	class RotationModel
	{
	public:
		// Calculating accelerations
		general::math::Quaternion function(
			const general::math::Quaternion& rot,
			const general::time::JD& t)
		{
			return static_cast<Model*>(this)->function(rot, t);
		}
	};
}