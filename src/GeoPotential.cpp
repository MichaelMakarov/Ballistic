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
		_pnm.resize(dim + 2);
		_pnm[dim] = _pnm[dim + 1] = 0.0;
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
		const double xy{ std::sqrt(point[0] * point[0] + point[1] * point[1]) };
		const double sinphi{ point[2] / r };
		const double cosphi{ xy / r };
		const double coslambda{ point[0] / xy };
		const double sinlambda{ point[1] / xy };
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

	constexpr inline double delta(const size_t m) noexcept { return 1.0 - 0.5 * (m == 0); }
	/// <summary>
	/// a derivative of the legendre function
	/// </summary>
	/// <param name="pnm"> - a legendre function</param>
	/// <param name="pnm1"> - a legendre function of the next derivative degree</param>
	/// <param name="n"> - a degree of the legendre function</param>
	/// <param name="m"> - derivative degree of the legendre function</param>
	/// <param name="tgphi">is z / sqrt(1 - z^2)</param>
	/// <returns>a derivative of the legendre function</returns>
	double dpnm(
		const double pnm, const double pnm1, 
		const size_t n, const size_t m, 
		const double tgphi) noexcept
	{
		return -pnm * tgphi * m + pnm1 * (m < n) * std::sqrt((n - m) * (n + m + 1) * delta(m));
	}

	general::math::Vec3 GeoPotential::derivatives(const general::math::Vec3& point)
	{
		using namespace general::math;
		
		const double r{ point.length() };
		const double xy{ std::sqrt(point[0] * point[0] + point[1] * point[1]) };
		const double sinphi{ point[2] / r };
		const double cosphi{ xy / r };
		const double tgphi{ sinphi / cosphi };
		const double coslambda{ point[0] / xy };
		const double sinlambda{ point[1] / xy };
		const double mu_r2{ _eMu / r / r };
		const double R_r{ _eR / r };
		double mult{ 1 };
		const Matrix3x3 CT = Matrix3x3({
			{ cosphi * coslambda, -sinphi * coslambda, -sinlambda },
			{ cosphi * sinlambda, -sinphi * sinlambda, coslambda },
			{ sinphi, cosphi, 0 }
		});

		size_t k{ 0 };
		double kcs, ksc, poly;
		Vec3 dUn, dUsum;
		
		calc_trigonometric(coslambda, sinlambda);
		calc_polynoms(cosphi, sinphi);

		// calculating the potential derivatives
		for (size_t n = 0; n <= _count; ++n) {
			dUn[0] = dUn[1] = dUn[2] = 0.0;
			for (size_t m = 0; m <= n; ++m) {
				poly = _pnm[k];
				// Cnm * cos(m * L) + Snm * sin(m * L)
				kcs = _harmonics[k].first * _cs[m].first + _harmonics[k].second * _cs[m].second;
				// Snm * cos(m * L) - Cnm * sin(m * L)
				ksc = _harmonics[k].second * _cs[m].first - _harmonics[k].first * _cs[m].second;
				dUn[0] -= poly * kcs;										// a derivative by r
				dUn[1] += dpnm(_pnm[k], _pnm[k + 1], n, m, tgphi) * kcs;	// a derivative by phi
				dUn[2] += poly * ksc * m;									// a derivative by lambda
				++k;
			}
			dUsum[0] += (n + 1) * mult * dUn[0];
			dUsum[1] += mult * dUn[1];
			dUsum[2] += mult * dUn[2];
			mult *= R_r;
		}
		dUsum[2] /= cosphi;
		dUsum *= mu_r2;
		return CT * dUsum;
	}

	std::pair<general::math::Vec3, general::math::Matrix3x3> GeoPotential::fullderivatives(const general::math::Vec3& point)
	{
		using namespace general::math;

		const double r{ point.length() };
		const double xy{ std::sqrt(point[0] * point[0] + point[1] * point[1]) };
		const double sinphi{ point[2] / r };
		const double cosphi{ xy / r }, cosphi2{ cosphi * cosphi };
		const double tgphi{ sinphi / cosphi };
		const double coslambda{ point[0] / xy };
		const double sinlambda{ point[1] / xy };
		const double mu_r2{ _eMu / r / r };
		const double R_r{ _eR / r };
		double mult{ 1 };
		const auto CT = Matrix3x3({
			{ cosphi * coslambda, -sinphi * coslambda, -sinlambda },
			{ cosphi * sinlambda, -sinphi * sinlambda, coslambda },
			{ sinphi, cosphi, 0 }
		});
		const auto C = Matrix3x3({
			{ cosphi * coslambda, cosphi * sinlambda, sinphi },
			{ -sinphi * coslambda, -sinphi * sinlambda, cosphi },
			{ -sinlambda, coslambda, 0 }
		});
		Matrix3x3 G;
		size_t k{ 0 };
		double kcs, ksc, dpoly, poly;
		Vec<3> dUn, dUsum, ddUn;
		Vec<6> ddUsum;

		calc_trigonometric(coslambda, sinlambda);
		calc_polynoms(cosphi, sinphi);

		// calculating the potential derivatives
		for (size_t n = 0; n <= _count; ++n) {
			dUn[0] = dUn[1] = dUn[2] = 0.0;
			ddUn[0] = ddUn[1] = ddUn[2] = 0.0;
			for (size_t m = 0; m <= n; ++m) {
				poly = _pnm[k];
				dpoly = dpnm(_pnm[k], _pnm[k + 1], n, m, tgphi);
				// Cnm * cos(m * L) + Snm * sin(m * L)
				kcs = _harmonics[k].first * _cs[m].first + _harmonics[k].second * _cs[m].second;
				// Snm * cos(m * L) - Cnm * sin(m * L)
				ksc = _harmonics[k].second * _cs[m].first - _harmonics[k].first * _cs[m].second;
				dUn[0] -= poly * kcs;		// a derivative by r
				dUn[1] += dpoly * kcs;		// a derivative by phi
				dUn[2] += poly * ksc * m;	// a derivative by lambda
				// a double derivative by phi
				ddUn[0] += dpnm(_pnm[k + 1], _pnm[k + 2], n, m + 1, tgphi) - m * (poly / cosphi2 + dpoly * tgphi);
				ddUn[1] += dpoly * ksc * m;	// a derivative by phi and lambda
				ddUn[2] = dUn[0] * m * m;	// a double derivative by lambda
				++k;
			}
			dUsum[0] += (n + 1) * mult * dUn[0];
			dUsum[1] += mult * dUn[1];
			dUsum[2] += mult * dUn[2];
			ddUsum[0] -= (n + 2) * (n + 1) * mult * dUn[0];	// a double derivative by r
			ddUsum[1] -= (n + 1) * mult * dUn[1];			// a derivative by r and phi
			ddUsum[2] -= (n + 1) * mult * dUn[2];			// a derivative by r and lambda
			ddUsum[3] += mult * ddUn[0];					// a double derivative by phi
			ddUsum[4] += mult * ddUn[1];					// a derivative by phi and lambda
			ddUsum[5] += mult * ddUn[2];					// a double derivative by lambda
			mult *= R_r;
		}
		dUsum[2] /= cosphi;
		dUsum *= mu_r2;
		ddUsum *= mu_r2 / r;
		ddUsum[2] /= cosphi;
		ddUsum[4] /= cosphi;
		ddUsum[5] /= cosphi * cosphi;
		// matrix of partial derivatives filling
		G(0, 0) = ddUsum[0];
		G(0, 1) = G(1, 0) = ddUsum[1] - dUsum[1] / r;
		G(0, 2) = G(2, 0) = ddUsum[2] - dUsum[2] / r;
		G(1, 1) = dUsum[0] + ddUsum[3];
		G(1, 2) = G(2, 1) = tgphi * dUsum[2] / r + ddUsum[4];
		G(2, 2) = (dUsum[0] - tgphi * dUsum[1]) / r + ddUsum[5];

		return std::make_pair(CT * dUsum, CT * G * C);
	}
}