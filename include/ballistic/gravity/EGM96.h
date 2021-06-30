﻿#pragma once
#include "EarthModel.h"

namespace ball
{
	class EGM96 : public IEarth
	{
	public:
		inline size_t count() const override { return _count; }
		const std::vector<std::pair<double, double>>& harmonics() const override { return _harmonics; }

		constexpr static inline double Mu() { return 0.3986004415E+15; }
		constexpr static inline double R() { return 0.6378136300E+07; }
		constexpr static inline double Ec2() { return 6.69437999014e-3; }
		constexpr static inline double W() { return 72.92115e-6; }
		constexpr static inline double Fl() { return 1.0 / 298.257223563; }

		EGM96() = default;
		EGM96(const char* filepath);
		~EGM96() = default;

	private:
		static std::vector<std::pair<double, double>> _harmonics;
		static size_t _count;
	};
}