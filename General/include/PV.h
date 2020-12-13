#pragma once
#include "XYZ.h"
#include "RBL.h"

namespace ball
{
	namespace types
	{
		struct PV
		{
			double Pos[3], Vel[3];

			PV(
				const double px,
				const double py, 
				const double pz,
				const double vx,
				const double vy,
				const double vz) : 
					Pos{ px, py, pz },
					Vel{ vx, vy, vz }
			{}
			PV(
				const geometry::XYZ& position, 
				const geometry::XYZ& velocity) :
					Pos{ position.X, position.Y, position.Z },
					Vel{ velocity.X, velocity.Y, velocity.Z }
			{}
			PV(
				const double r,
				const double b,
				const double l,
				const double dr,
				const double db,
				const double dl) :
					Pos{ r, b, l },
					Vel{ dr, db, dl }
			{}
			PV(
				const geometry::RBL& position,
				const geometry::RBL& velocity) :
					Pos{ position.R, position.B, position.L },
					Vel{ velocity.R, velocity.B, velocity.L }
			{}
			PV(const PV& pv)
			{
				Pos[0] = pv.Pos[0];
				Pos[1] = pv.Pos[1];
				Pos[2] = pv.Pos[2];
				Vel[0] = pv.Vel[0];
				Vel[1] = pv.Vel[1];
				Vel[2] = pv.Vel[2];
			}

		};
	}
}