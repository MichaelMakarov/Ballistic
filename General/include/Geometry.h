#pragma once
#include <ostream>

namespace ball
{
	namespace geometry
	{
		// Struct represents point in the orthogonal coordinate system.
		// Implements adding, substraction, multiplying, can be used for vector implementing.
		struct XYZ
		{
			double X, Y, Z;
			XYZ() : X(0), Y(0), Z(0) {}
			XYZ(
				const double x,
				const double y,
				const double z) : X(x), Y(y), Z(z)
			{}
			~XYZ() {}

			double length() const;

			XYZ& operator += (const XYZ& v);
			XYZ& operator -= (const XYZ& v);
			double operator *= (const XYZ& v);
			XYZ& operator /= (const double n);
			XYZ& operator *= (const double n);

			friend std::ostream& operator << (std::ostream& o, const XYZ& v);

			friend XYZ operator + (const XYZ& f, const XYZ& s);
			friend XYZ operator - (const XYZ& f, const XYZ& s);
			friend XYZ operator * (const XYZ& v, const double n);
			friend XYZ operator / (const XYZ& v, const double n);
			friend XYZ operator * (const double n, const XYZ& v);
		};

		struct RBL
		{
			double R, B, L;

			RBL() : R(0), B(0), L(0) {}
			RBL(
				const double r,
				const double b,
				const double l) : R(r), B(b), L(l)
			{}
			~RBL() {}

			RBL& operator += (const RBL& v);
			RBL& operator -= (const RBL& v);
			RBL& operator /= (const double n);
			RBL& operator *= (const double n);

			friend std::ostream& operator << (std::ostream& o, const RBL& v);

			friend RBL operator + (const RBL& f, const RBL& s);
			friend RBL operator - (const RBL& f, const RBL& s);
			friend RBL operator * (const RBL& v, const double n);
			friend RBL operator / (const RBL& v, const double n);
			friend RBL operator * (const double n, const RBL& v);

		};

		// Struct implements the vector of 6th dimension.
		// Contains 3 position and 3 velocity values.
		struct PV
		{
			double P1, P2, P3, V1, V2, V3;

			PV() : P1(0), P2(0), P3(0), V1(0), V2(0), V3(0)
			{}
			PV(
				const double px,
				const double py,
				const double pz,
				const double vx,
				const double vy,
				const double vz) :
				P1{ px }, P2{ py }, P3{ pz },
				V1{ vx }, V2{ vy }, V3{ vz }
			{}
			PV(
				const XYZ& position,
				const XYZ& velocity) :
				P1{ position.X },
				P2{ position.Y },
				P3{ position.Z },
				V1{ velocity.X },
				V2{ velocity.Y },
				V3{ velocity.Z }
			{}
			PV(
				const RBL& position,
				const RBL& velocity) :
				P1{ position.R },
				P2{ position.B },
				P3{ position.L },
				V1{ velocity.R },
				V2{ velocity.B },
				V3{ velocity.L }
			{}
			PV(const PV& pv) :
				P1{ pv.P1 }, P2{ pv.P2 }, P3{ pv.P3 },
				V1{ pv.V1 }, V2{ pv.V2 }, V3{ pv.V3 }
			{}

			PV& operator += (const PV& pv);
			PV& operator -= (const PV& pv);
			PV& operator *= (const double m);

			friend PV operator + (const PV& f, const PV& s);
			friend PV operator - (const PV& f, const PV& s);
			friend PV operator * (const double m, const PV& pv);
			friend PV operator * (const PV& pv, const double m);

			friend std::ostream& operator << (std::ostream& o, const PV& pv);
		};
	}
}