#-------------------------------------------------------------------------------
# Created 19.12.22
# Author: e.s.oyre@astro.uio.no
#-------------------------------------------------------------------------------
#
#                Constants.jl
#
#-------------------------------------------------------------------------------
# Module containing physical constants
#-------------------------------------------------------------------------------

# The defining constants for SI units
const k_B = 1.380649e-23 # J/K............................. Boltzmann's constant
const e = 1.602176634e-19 # C................................. Elementary charge
const c = 299792458 # m/s.............................. Speed of light in vacuum
const h = 6.62607015e-34 # J Hz^-1............................ Planck's constant
# Other constants
const G = 6.67430e-11 # m^3 kg^-1 s^-2.................. Constant of gravitation
const μ_0 = 1.25663706212e-6 # N A^-2...............Vacuum magnetic permeability
const ϵ_0 = 8.8541878128e-12 # F m^-1...............Vacuum electric permittivity
# Masses
const m_e = 9.1093837015e-31 # kg ................................ Electron mass
const m_p = 1.67262192369e-27 # kg ................................. Proton mass

# CGS units
const c_cgs = 1e2c # cm/s
const k_B_cgs = 1.380649e-16 # erg/K
const e_cgs = 4.8032047e-10 # statcoulombs
const m_e_cgs = m_e * 1e3 # grams
const m_p_cgs = m_p * 1e3 # grams

# Derived quantities
const csqrd = c^2 # s^2 m^-2 ........................... The light speed squared
const csqrdinv = 1 / csqrd # s^2 m^-2 ....... Inverse of the light speed squared

# Unit conversion
const J2eV = 6.24150907e18 # 1 J * 6.24e12 eV/J...........................energy
const si2cgs_mass = 1e3 # 1 kg * 1000 g/kg
const si2cgs_numberdensity = 1e-6 # 1/m^3 * 1/(100 cm/m)^3
const si2cgs_velocity = 1e2 # 1 m/s * 100 cm/(m*s)
const si2cgs_charge = 1e-1c_cgs # = 2.99792458e9.  C => statC
const si2cgs_electriccurrent = 1e-1c_cgs # - 2.99792458e9. A => StatC/s = StatA
const si2cgs_currentdensity = 1e-5c_cgs # = 2.99792458e5. A/m^2 => StatA/cm^2
const si2cgs_electricfield = 1e6 / c_cgs # 3.33564095e-5. V/m => StatV/cm
const cgs2si_electricfield = 1/si2cgs_electricfield # = 2.99792458e4. To V/m
const cgs2si_velocity = 1e-2 # 1 cm * 0.01 m/cm

# Particles
const electron_mass = m_e
const electron_charge = -e
const proton_mass = m_p
const proton_charge = e
