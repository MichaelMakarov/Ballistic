#pragma once

namespace ball
{
    namespace tasks
    {
        namespace physical
        {
            // air molar mass at zero height (kg / kmole)
            inline constexpr double Mc{ 28.964420 };
            // Avogadro number (1 / kmole)
            inline constexpr double Na{ 602.257e24};
            // universal gas constant (Joule / (kmole * K))
            inline constexpr double Rs{ 8314.32 };
            // specific gas constant (Joule / (kg * K))
            inline constexpr double R{ 287.05287 };
            // empiric Saterpend's coefficient (K)
            inline  constexpr double S{ 110.4 };
            // empiric Saterpend's coefficient (kgs / (s * m * K^(1/2)))
            inline constexpr double Bs{ 1.458e-6 };
            // adiabate coefficient (Cp / Cv)
            inline constexpr double X{ 1.4 };
            // effective diameter of the molecule (m)
            inline constexpr double SIGMA{ 0.365e-9 };
            // average density of Earth
            inline constexpr double Rho_E{ 5510 };
        }
    }
}