#pragma once
#include "TranslationModel.h"
#include "GeoPotential.h"
#include "Atmosphere1981.h"
#include "Atmosphere2004.h"
#include "Conversions.h"
#include "SolarSystem.h"

namespace ball
{
	/// <summary>
	/// A model of forecast implements static atmosphere standard and geopotential model
	/// </summary>
	class StatAtmModel : public TranslationModel<StatAtmModel>
	{
	private:
		GeoPotential _geopotential;
		Atmosphere1981 _atmosphere;
		double _eW, _eFl, _eR;

	public:
		StatAtmModel(
			const double eMu, const double eR, 
			const double eW, const double eFl,
			const IEarth& gravity,
			const size_t harmonics) :
			_atmosphere{ Atmosphere1981(eR, eFl) },
			_geopotential{ GeoPotential(eMu, eR, gravity, harmonics) },
			_eW{ eW },
			_eFl{ eFl },
			_eR{ eR }
		{}
		~StatAtmModel() = default;

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
	/// </summary>
	class DynAtmModel : public TranslationModel<DynAtmModel>
	{
		using EV = general::math::Vec<13>;
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
		/// <param name="gravity"> - a gravity model</param>
		/// <param name="harmonics"> - a number of harmonics</param>
		DynAtmModel(
			const double eMu, const double eR,
			const double eW, const double eFl,
			const double f10_7, const double f81, const double kp,
			const IEarth& gravity,
			const size_t harmonics) :
			_atmosphere{ Atmosphere2004(eR, eFl, f10_7, f81, kp) },
			_geopotential{ GeoPotential(eMu, eR, gravity, harmonics) },
			_eW{ eW },
			_eFl{ eFl },
			_eR{ eR }
		{}
		~DynAtmModel() = default;

		/// <summary>
		/// Calculating the acceleration using position vector and time
		/// </summary>
		/// <param name="vec"> - vector in GCS</param>
		/// <param name="t"> - relevant time</param>
		/// <returns>an aceleration</returns>
		Vec6 function(const Vec6& vec, const general::time::JD& t);

		EV extfunction(const EV& vec, const general::time::JD& t);
	};
	/// <summary>
	/// A model of forecast implements:
	/// geopotential model,
	/// dynamic atmosphere standard (GOST 2004) with variable ballistic parameter,
	/// simple solar and lunar potential models
	/// </summary>
	class VarBallModel : public TranslationModel<VarBallModel>
	{
		/// <summary>
		/// 16d vector (x, y, z, vx, vy, vz, dx, dy, dz, dvx, dvy, dvz, sb, sam, sfr, sph)
		/// </summary>
		using EV = general::math::Vec<16>;
	private:
		GeoPotential _geopotential;
		Atmosphere2004 _atmosphere;
		double _eW, _eFl, _eR;
	public:
		double sAmpl = 0.0;	// an amplitude of ballistic coefficient variation
		double sFreq = 1.0;	// a frequency of ballistic coefficient variation
		double sPhas = 0.0;	// a phase of ballisitc coefficient variation
	private:
		double ball_coefficient(const general::time::JD& t) noexcept
		{
			return this->sBall + this->sAmpl * std::sin(this->sFreq * t.to_double() + sPhas);
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
		VarBallModel(
			const double eMu, const double eR,
			const double eW, const double eFl,
			const double f10_7, const double f81, const double kp,
			const IEarth& gravity,
			const size_t harmonics) :
			_atmosphere{ Atmosphere2004(eR, eFl, f10_7, f81, kp) },
			_geopotential{ GeoPotential(eMu, eR, gravity, harmonics) },
			_eW{ eW },
			_eFl{ eFl },
			_eR{ eR }
		{}
		~VarBallModel() = default;

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
		EV extfunction(const EV& vec, const general::time::JD& t);
	};
}