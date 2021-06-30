#pragma once
#include "TranslationModel.h"
#include "GeoPotential.h"
#include "Atmosphere1981.h"
#include "Atmosphere2004.h"
#include "Conversions.h"
#include "SolarSystem.h"
#include "general/Polynomial.h"
#include <iostream>

namespace ball
{
	/// <summary>
	/// vector of 13 dim (x, y, z, vx, vy, vz, dx, dy, dz, dvx, dvy, dvz, ds)
	/// </summary>
	using Vec13 = general::math::Vec<13>;
	/// <summary>
	/// vector if 16 dim (x, y, z, vx, vy, vz, dx, dy, dz, dvx, dvy, dvz, ds)
	/// </summary>
	using Vec16 = general::math::Vec<16>;
	/// <summary>
	/// A model of forecast implements static atmosphere standard and geopotential model
	/// </summary>
	class MSA : public TranslationModel<MSA>
	{
	private:
		GeoPotential _geopotential;
		Atmosphere1981 _atmosphere;
		double _eW, _eFl, _eR;
	public:
		double sBall;

	public:
		MSA(
			const IEarth& gravity,
			const size_t harmonics,
			const double eMu, const double eR, 
			const double eW, const double eFl,
			const double sb = 0.0) :
			_atmosphere{ Atmosphere1981(eR, eFl) },
			_geopotential{ GeoPotential(eMu, eR, gravity, harmonics) },
			_eW{ eW },
			_eFl{ eFl },
			_eR{ eR },
			sBall{ sb }
		{}
		~MSA() = default;

		/// <summary>
		/// acelerations calculation using current  and time
		/// </summary>
		/// <param name="vec"> - vector in GCS</param>
		/// <param name="t"> - relevant time</param>
		/// <returns>an acceleration</returns>
		Vec6 function(const Vec6& vec, const general::time::JD& t);
	};
	/// <summary>
	/// A model of forecast implements:
	/// geopotential model,
	/// dynamic atmosphere standard (GOST 2004),
	/// simple solar and lunar potential models
	/// 
	/// </summary>
	class MDASM : public TranslationModel<MDASM>
	{
		using Vec7 = general::math::Vec<7>;
	public:
		double sBall;
	private:
		GeoPotential _geopotential;
		Atmosphere2004 _atmosphere;
		double _eW, _eFl, _eR;
	public:
		/// <summary>
		/// Creating a sample of a class
		/// </summary>
		/// <param name="eMu"> - gravitational constant</param>
		/// <param name="eR"> - Earth's equator</param>
		/// <param name="eW"> - angular velocity of rotation</param>
		/// <param name="eFl"> - Earth;s flatening</param>
		/// <param name="f10_7"> - daily averaged ISA</param>
		/// <param name="f81"> - 81 daily averaged ISA</param>
		/// <param name="kp"> - daily averaged IGP</param>
		/// <param name="sb"> - ballistic coefficient</param>
		/// <param name="gravity"> - a gravity model</param>
		/// <param name="harmonics"> - a number of harmonics</param>
		MDASM(
			const IEarth& gravity,
			const size_t harmonics,
			const double eMu, const double eR,
			const double eW, const double eFl,
			const double f10_7, const double f81, const double kp,
			const double sb = 0.0) :
			_atmosphere{ Atmosphere2004(eR, eFl, f10_7, f81, kp) },
			_geopotential{ GeoPotential(eMu, eR, gravity, harmonics) },
			_eW{ eW },
			_eFl{ eFl },
			_eR{ eR },
			sBall{ sb }
		{}
		~MDASM() = default;

		/// <summary>
		/// Calculating the acceleration using position vector and time
		/// </summary>
		/// <param name="vec"> - vector in GCS</param>
		/// <param name="t"> - relevant time</param>
		/// <returns>an aceleration</returns>
		Vec6 function(const Vec6& vec, const general::time::JD& t);

		general::math::Vec<55> extfunction(const general::math::Vec<55>& vec, const general::time::JD& t);
	};

	/// <summary>
	/// A model of forecast implements:
	/// geopotential model,
	/// dynamic atmosphere standard (GOST 2004) with variable ballistic parameter,
	/// simple solar and lunar potential models
	/// </summary>
	class MDASMVb : public TranslationModel<MDASMVb>
	{
	private:
		GeoPotential _geopotential;
		Atmosphere2004 _atmosphere;
		double _eW, _eFl, _eR, _eMu;
		general::math::PowerPolynomial<4> _poly;

	private:
		double ball_coefficient(
			const general::math::Vec3& pos, 
			const general::math::Vec3& vel,
			const general::time::JD& t) noexcept {
			double sidt = sidereal_time_mean(t);
			const auto oscul = oscul_from_ACS(GCS_to_ACS(pos, sidt), GCS_to_ACS(vel, sidt), _eMu);
			return _poly(oscul.latitudearg);
		}

	public:
		/// <summary>
		/// Creating a sample of a class
		/// </summary>
		/// <param name="eMu"> - gravitational constant</param>
		/// <param name="eR"> - Earth's equator</param>
		/// <param name="eW"> - angular velocity of rotation</param>
		/// <param name="eFl"> - Earth;s flatening</param>
		/// <param name="f10_7"> - daily averaged ISA</param>
		/// <param name="f81"> - 81 daily averaged ISA</param>
		/// <param name="kp"> - daily averaged IGP</param>
		/// <param name="gravity"> - a gravity model</param>
		/// <param name="harmonics"> - a number of harmonics</param>
		MDASMVb(
			const IEarth& gravity,
			const size_t harmonics,
			const double eMu, const double eR,
			const double eW, const double eFl,
			const double f10_7, const double f81, const double kp,
			const std::array<double, 5>& sb) :
			_atmosphere{ Atmosphere2004(eR, eFl, f10_7, f81, kp) },
			_geopotential{ GeoPotential(eMu, eR, gravity, harmonics) },
			_eW{ eW },
			_eFl{ eFl },
			_eR{ eR },
			_eMu{ eMu }
		{
			std::memcpy(_poly.data, sb.data(), sizeof(sb));
		}
		~MDASMVb() = default;

		/// <summary>
		/// Calculating the right part of the differential equation using vector and time
		/// </summary>
		/// <param name="vec"> - vector in GCS</param>
		/// <param name="t"> - relevant time</param>
		/// <returns>a vector of position and velocity</returns>
		Vec6 function(const Vec6& vec, const general::time::JD& t);
		/// <summary>
		/// Calculating the extended right part of the differential equations with variations
		/// </summary>
		/// <param name="vec">is a vector in GCS</param>
		/// <param name="t">is relevant time</param>
		/// <returns>an extended vector</returns>
		Vec16 extfunction(const Vec16& vec, const general::time::JD& t);
	};
}