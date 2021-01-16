#include "Structures.h"
#include "GeneralConstants.h"
#include "Conversions.h"

namespace ball
{
	namespace geometry
	{
        std::ostream& operator<<(std::ostream& os, const State& p)
        {
            os << "T: " << p.T.to_datetime() << "; Vec: " << p.Vec << "; s = " <<
                p.Sb << "; loop = " << p.Loop;
            return os;
        }

        std::ostream& operator<<(std::ostream& os, const Oscul& p)
        {
            os << "T: " << p.T << "; a = " << p.A << "; e = " << p.E << 
                "; i = " << p.I << "; w = " << p.W << "; o = " << p.O <<
                "; M = " << p.M << "; s = " << p.Sb << "; loop = " << p.Loop;
            return os;
        }

        bool load_tle(std::istream& stream, TLE& f)
        {
            if (!stream) return false;
            char buf[16] = {};
            size_t i;
            stream.read(buf, 2);
            if (std::strcmp(buf, "1 ") != 0)
            {
                f.Header.push_back(buf[0]);
                f.Header.push_back(buf[1]);
                while (buf[0] != '\n')
                {
                    stream.read(buf, 1);
                    f.Header.push_back(buf[0]);
                }
            }
            stream.read(buf, 2);
            stream.read(buf, 5);
            f.Satellite = std::atol(buf);
            for (i = 0; i < 5; ++i) buf[i] = 0;
            stream.read(&f.Classification, 1);
            stream.read(buf, 3);
            f.LaunchYear = std::atol(buf);
            stream.read(buf, 3);
            f.LaunchNumber = std::atol(buf);
            stream.read(f.LaunchPlace, 3);
            stream.read(buf, 3);
            f.EpochYear = std::atol(buf);
            stream.read(buf, 12);
            f.EpochTime = std::atof(buf);
            for (i = 0; i < 12; ++i) buf[i] = 0;
            stream >> f.FirstAcc >> f.SecondAcc >> f.Deceleration;
            while (!stream.eof()) {
                stream.read(buf, 1);
                if (buf[0] == '\n') break;
            }
            stream.read(buf, 8);
            stream >> f.Inclination >> f.AscendLong >>
                f.Eccentricity >> f.PeriapsArg >> f.AvrAnomaly;
            f.Eccentricity *= 1e-7;
            stream.read(buf, 12);
            f.Frequency = std::atof(buf);
            for (i = 0; i < 12; ++i) buf[i] = 0;
            stream.read(buf, 5);
            f.Loop = std::atol(buf);
            return true;
        }

        Oscul TLE::to_oscul(const double mu) const
        {
            auto jd{ time::JD(time::JD2000 + EpochTime) };
            for (size_t i = 0; i < EpochYear; ++i)
                if (i % 4 == 0) jd.add_days(366);
                else jd.add_days(365);
            const double a = space::semimajoraxis_from_period(
                static_cast<double>(time::SEC_PER_DAY) / Frequency,
                mu);
            return geometry::Oscul(
                a,
                Eccentricity,
                Inclination,
                PeriapsArg,
                AscendLong,
                AvrAnomaly,
                0.0,
                jd,
                Loop);
        }
	}
}