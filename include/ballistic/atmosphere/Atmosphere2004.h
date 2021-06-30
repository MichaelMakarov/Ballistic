#include "Atmosphere.h"
#include "general/Vector.h"

namespace ball
{
	/// <summary>
	/// A model of atmosphere according GOST2004
	/// </summary>
	class Atmosphere2004 : public IAtmosphere<Atmosphere2004, general::math::Vec3, double>
	{
	private:
		double dynamic_density(
			const double h, 
			const general::math::Vec3& position,
			const general::time::JD& time,
			const general::math::Vec3& sun,
			const double inclination) const;
		double static_density(const double h) const;

	public:
		/// <summary>
		/// creating an atmosphere object according to GOST 2004
		/// </summary>
		/// <param name="eR"> - Earth's radius</param>
		/// <param name="eFl"> - Earth's flattening</param>
		/// <param name="F10_7"> - daily averaged index of solar activity (ISA)</param>
		/// <param name="F81"> - ISA averaged during last 81 days</param>
		/// <param name="Kp"> - quasilogarithmic daily averaged index of geomagnetic perturbations (IGP)</param>
		Atmosphere2004(
			const double eR,
			const double eFl,
			const double F10_7,
			const double F81, 
			const double Kp) : _f81{ F81 }, _f10_7{ F10_7 }, _kp{ Kp }, _eR{ eR }, _eFl{ eFl }
		{ }
		~Atmosphere2004() = default;

		/// <summary>
		/// calculating the density of atmosphere
		/// </summary>
		/// <param name="position"> - a point</param>
		/// <param name="time"> - JD</param>
		/// <param name="sun"> - solar position</param>
		/// <param name="inclination"> - solar inclination</param>
		/// <returns>the density</returns>
		double density(
			const general::math::Vec3& position,
			const general::time::JD& time,
			const general::math::Vec3& sun,
			const double inclination) const;

	private:
		const double _eR;
		const double _eFl;
		const double _f81;
		const double _f10_7;
		const double _kp;

		constexpr inline const static double ISA_STEP = 25.0;
		constexpr inline const static double LOWER_H = 120.0;
		constexpr inline const static double UPPER_H = 1500.0;
		constexpr inline const static double NIGHT_DENSITY = 1.58868e-8;

		const static double Isa[6];

		const static double Heights[4];
		const static double C[4][3];

		const static double Aheights[7];
		const static double Ah[7][7];
		const static double Al[7][7];

		const static double Bheights[7];
		const static double Bh[7][5];
		const static double Bl[7][5];

		const static double Cheights[7];
		const static double Ch[7][5];
		const static double Cl[7][5];

		const static double D[7][5];

		const static double Eheights[7];
		const static double Eh[7][9];
		const static double El[7][9];

		const static double Et[7][4];

		const static double Lheights[7];
		const static double Lh[7][5];
		const static double Ll[7][5];

		const static double N[3];
		const static double A[9];
		const static double Phi[7];
	};
}