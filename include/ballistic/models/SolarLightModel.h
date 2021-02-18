#pragma once
#include "RotationModel.h"

namespace ball
{
	
	class SolarLightModel : public RotationModel<SolarLightModel>
	{
		SolarLightModel() = default;
		~SolarLightModel() = default;
	};
}