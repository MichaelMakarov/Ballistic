#pragma once
#include "XYZ.h"
#include "RBL.h"

namespace ball
{
	namespace types
	{
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
				const geometry::XYZ& position, 
				const geometry::XYZ& velocity) :
					P1{ position.X }, 
					P2{ position.Y }, 
					P3{ position.Z },
					V1{ velocity.X }, 
					V2{ velocity.Y }, 
					V3{ velocity.Z }
			{}
			PV(
				const geometry::RBL& position,
				const geometry::RBL& velocity) :
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