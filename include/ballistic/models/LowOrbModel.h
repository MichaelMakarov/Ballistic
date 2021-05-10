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
		general::math::PV function(const general::math::PV& vec, const general::time::JD& t)
		{
			double w_2{ _eW * _eW };
			double h = GCS_height_from_position(vec.pos, _eR, _eFl);
			if (h > this->MinHeight && h < this->MaxHeight) {
				// geopotential aceleration with a centrifugal and a coriolis force
				auto acpot{ _geopotential.acceleration(vec.pos) };
				double density{ _atmosphere.density(vec.pos, t) };
				double acatm = vec.vel.length() * density * this->sBall;
				// atmosphere aceleration a = v * s * rho, 
				// s - a ballistic coefficient,
				// v - a velocity of the vehicle,
				// rho - a density of the atmosphere
				// the addition all the components
				return general::math::PV(
					vec.vel.x(), vec.vel.y(), vec.vel.z(),
					acpot.x() + w_2 * vec.pos.x() + 2 * _eW * vec.vel.y() - acatm * vec.vel.x(),
					acpot.y() + w_2 * vec.pos.y() - 2 * _eW * vec.vel.x() - acatm * vec.vel.y(),
					acpot.z() - acatm * vec.vel.z()
				);
			}
			throw std::runtime_error("Height is out of bounds!");
		}
	};
	/// <summary>
	/// A model of forecast implements dynamic atmosphere standard (GOST 2004) and geopotential model
	/// </summary>
	class DynAtmModel : public TranslationModel<DynAtmModel>
	{
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
		general::math::PV function(const general::math::PV& vec, const general::time::JD& t)
		{
			double w_2{ _eW * _eW };
			double h = GCS_height_from_position(vec.pos, _eR, _eFl);
			if (h > this->MinHeight && h < this->MaxHeight) {
				auto [sunsph, sunort] = Sun::positionACS(t);
				sunort = ACS_to_GCS(sunort, sidereal_time_avr(t));
				double density = _atmosphere.density(vec.pos, t, sunort, sunsph.y());
				double acatm = vec.vel.length() * density * this->sBall;
				// acelerations by geopotential, sun and moon
				auto acpot = _geopotential.acceleration(vec.pos);
				auto acsun = acceleration_by_masspoint(vec.pos, sunort, Sun::Mu());
				auto acmoon =  Moon::acceleration(vec.pos, t);
				acpot += acsun + acmoon;
				// atmosphere aceleration a = v * s * rho, 
				// s - a ballistic coefficient,
				// v - a velocity of the vehicle,
				// rho - a density of the atmosphere
				// the addition all the components
				return general::math::PV(
					vec.vel.x(), vec.vel.y(), vec.vel.z(),
					acpot.x() + w_2 * vec.pos.x() + 2 * _eW * vec.vel.y() - acatm * vec.vel.x(),
					acpot.y() + w_2 * vec.pos.y() - 2 * _eW * vec.vel.x() - acatm * vec.vel.y(),
					acpot.z() - acatm * vec.vel.z()
				);
			}
			throw std::runtime_error("Height is out of bounds!");
		}
	};
}