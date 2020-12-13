namespace ball
{
	namespace tasks
	{
		namespace pz90
		{
			// gravitational constant (PZ-90.11)
			inline constexpr double G{ 6.67259e-11 };
			// light velocity
			inline constexpr double C{ 299792258 };

			// semimajor axis of Earth's ellipsoid m
			inline constexpr double R{ 6378136 };
			// Earth's compression
			inline constexpr double Fl{ 1 / 298.2564151 };
			// a geocentric gravitational constant m/sec
			inline constexpr double Mu{ 398600.4418e9 };
			// Earth's eccentricity
			inline constexpr double Ec{ 0.08181910657053156 };
			// a velocity of Earth's rotation rad/sec
			inline constexpr double W{ 7.292115e-5 };
		}
	}
}