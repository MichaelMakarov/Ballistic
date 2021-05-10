#include "GeoPotential.h"
#include "Conversions.h"
#include <algorithm>

namespace ball
{
	void GeoPotential::calc_trigonometric(const double coslambda, const double sinlambda)
	{
		_cs[0].first = 1.0; _cs[0].second = 0.0;
		for (size_t i = 1; i <= _count; ++i) {
			_cs[i].first = _cs[i - 1].first * coslambda - _cs[i - 1].second * sinlambda;
			_cs[i].second = _cs[i - 1].second * coslambda + _cs[i - 1].first * sinlambda;
		}
	}

	void GeoPotential::calc_polynoms(const double cosphi, const double sinphi)
	{
		_pnm[0] = 1;
		_pnm[1] = sinphi * std::sqrt(3);
		_pnm[2] = cosphi * std::sqrt(3);
		size_t k{ 3 };
		for (size_t n = 2; n <= _count; ++n) {
			for (size_t m = 0; m < n; ++m) {
				_pnm[k] = std::sqrt(2 * n - 1) * sinphi * _pnm[k - n] - std::sqrt((n - 1 - m) * (n - 1 + m) / (2.0 * n - 3)) * _pnm[(k + 1) - n - n];
				_pnm[k++] *= std::sqrt((2.0 * n + 1) / ((n - m) * (n + m)));
			}
			_pnm[k++] = std::sqrt(1 + 0.5 / n) * cosphi * _pnm[k - n - 1];
		}
	}

	GeoPotential::GeoPotential(
		const double eMu, const double eR,
		const IEarth& model,
		const size_t count) : _harmonics{model.harmonics()}
	{
		_count = std::min(count, model.count());
		_eR = eR;
		_eMu = eMu;
		size_t dim = ((_count + 1) * (_count + 2)) / 2;
		_cs.resize(_count + 1);
		_pnm.resize(dim + 1);
		_pnm[dim] = 0.0;
	}


	GeoPotential::GeoPotential(GeoPotential&& gp) noexcept : _harmonics{ gp._harmonics }
	{
		_count = gp._count;
		_eR = gp._eR;
		_eMu = gp._eMu;
		gp._eR = gp._eMu = 0;
		gp._count = 0;
	}

	double GeoPotential::operator () (const general::math::Vec3& point)
	{
		double result{ 0 };
		const double r{ point.length() };
		const double xy{ std::sqrt(point.x() * point.x() + point.y() * point.y()) };
		const double sinphi{ point.z() / r };
		const double cosphi{ xy / r };
		const double coslambda{ point.x() / xy };
		const double sinlambda{ point.y() / xy };
		const double R_r{ _eR / r };
		double mult{ 1 };
		size_t k{ 0 };
		calc_trigonometric(coslambda, sinlambda);
		calc_polynoms(cosphi, sinphi);
		for (size_t n = 0; n <= _count; ++n) {
			for (size_t m = 0; m <= n; ++m) {
				result += mult * _pnm[k] * (_harmonics[k].first * _cs[m].first + _harmonics[k].second * _cs[m].second);
				k++;
			}
			mult *= R_r;
		}
		return _eMu / r * result;
	}

	constexpr double delta(const size_t m) { return 1.0 - 0.5 * (m == 0); }

	general::math::Vec3 GeoPotential::acceleration(const general::math::Vec3& point)
	{
		using namespace general::math;
		
		const double r{ point.length() };
		const double xy{ std::sqrt(point.x() * point.x() + point.y() * point.y()) };
		const double sinphi{ point.z() / r };
		const double cosphi{ xy / r };
		const double tgphi{ sinphi / cosphi };
		const double coslambda{ point.x() / xy };
		const double sinlambda{ point.y() / xy };
		const double mu_r2{ _eMu / r / r };
		const double R_r{ _eR / r };
		double mult{ R_r * R_r };
		const double zxyr{ point.z() / xy / r };

		const Vec3 dR{ point / r };
		const auto dPhi{ Vec3(-point.x() * zxyr, -point.y() * zxyr, xy / r) };
		const auto dLambda{ Vec3(-point.y() / xy, point.x() / xy, 0) };

		size_t k{ 3 };
		double poly, dpoly;
		double kcs, ksc;
		Vec3 dUn, dUsum;
		
		calc_trigonometric(coslambda, sinlambda);
		calc_polynoms(cosphi, sinphi);

		// calculating the accelerations
		for (size_t n = 2; n <= _count; ++n) {
			dUn.x() = dUn.y() = dUn.z() = 0.0;
			for (size_t m = 0; m <= n; ++m) {
				// current Legendre function
				poly = _pnm[k];
				// a derivative of the current Legendre function
				dpoly = -poly * m * tgphi + _pnm[k + 1] * std::sqrt((n - m) * (n + m + 1) * delta(m));
				// Cnm * cos(m * L) + Snm * sin(m * L)
				kcs = _harmonics[k].first * _cs[m].first + _harmonics[k].second * _cs[m].second;
				// Snm * cos(m * L) - Cnm * sin(m * L)
				ksc = _harmonics[k].second * _cs[m].first - _harmonics[k].first * _cs[m].second;
				dUn.x() -= poly * kcs;		// a derivative by r
				dUn.y() += dpoly * kcs;		// a derivative by phi
				dUn.z() += poly * ksc * m;	// a derivative by lambda
				++k;
			}
			dUsum.x() += (n + 1) * mult * dUn.x();
			dUsum.y() += mult * dUn.y();
			dUsum.z() += mult * dUn.z();
			mult *= R_r;
		}
		dUsum.z() /= cosphi;
		dUsum *= mu_r2;

		return Vec3(
			-mu_r2 * dR.x() + dUsum.x() * dR.x() + dUsum.y() * dPhi.x() + dUsum.z() * dLambda.x(),
			-mu_r2 * dR.y() + dUsum.x() * dR.y() + dUsum.y() * dPhi.y() + dUsum.z() * dLambda.y(),
			-mu_r2 * dR.z() + dUsum.x() * dR.z() + dUsum.y() * dPhi.z()
		);
	}
}