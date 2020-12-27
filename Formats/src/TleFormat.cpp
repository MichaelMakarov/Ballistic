#include "TleFormat.h"
#include <fstream>
#include <cstring>
#include "Constants.h"
#include <exception>
#include "AstroValues.h"
#include "PZ90.h"

namespace ball
{
    namespace formats
    {
        bool TleFromFile(const std::string& filepath, TleFormat& f)
        {
            std::ifstream fread(filepath);
            if (!fread) return false;
            char buf[16] = {};
            size_t i;
            fread.read(buf, 2);
            if (std::strcmp(buf, "1 ") != 0)
            {
                f.Header.push_back(buf[0]);
                f.Header.push_back(buf[1]);
                while (buf[0] != '\n')
                {
                    fread.read(buf, 1);
                    f.Header.push_back(buf[0]);
                }
            }
            fread.read(buf, 2);
            fread.read(buf, 5);
            f.Satellite = std::atol(buf);
            for (i = 0; i < 5; ++i) buf[i] = 0;
            fread.read(&f.Classification, 1);
            fread.read(buf, 3);
            f.LaunchYear = std::atol(buf);
            fread.read(buf, 3);
            f.LaunchNumber = std::atol(buf);
            fread.read(f.LaunchPlace, 3);
            fread.read(buf, 3);
            f.EpochYear = std::atol(buf);
            fread.read(buf, 12);
            f.EpochTime = std::atof(buf);
            for (i = 0; i < 12; ++i) buf[i] = 0;
            fread >> f.FirstAcc >> f.SecondAcc >> f.Deceleration;
            while (!fread.eof()) {
                fread.read(buf, 1);
                if (buf[0] == '\n') break;
            }
            fread.read(buf, 8);
            fread >> f.Inclination >> f.AscendLong >>
                f.Eccentricity >> f.PeriapsArg >> f.AvrAnomaly;
            f.Eccentricity *= 1e-7;
            fread.read(buf, 12);
            f.Frequency = std::atof(buf);
            for (i = 0; i < 12; ++i) buf[i] = 0;
            fread.read(buf, 5);
            f.LoopNumber = std::atol(buf);
            fread.close();
            return true;
        }

        types::OsculParams TleFormat::ToOsculParams() const
        {
            time::JD jd = time::JD2000 + EpochTime;
            for (size_t i = 0; i < EpochYear; ++i)
                if (i % 4 == 0) jd.AddDays(366);
                else jd.AddDays(365);
            return types::OsculParams(
                tasks::SemimajorAxisFromPeriod(
                    static_cast<double>(time::SEC_PER_DAY) / Frequency, 
                    tasks::pz90::Mu),
                Eccentricity,
                Inclination,
                PeriapsArg,
                AscendLong,
                AvrAnomaly,
                0.0,
                jd,
                LoopNumber);
        }
    }
}

