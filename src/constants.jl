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
k_B = 1.380649e-23 # J/K..................... Boltzmann's constant
e = 1.602176634e-19 # C....................... Elementary charge
c = 299792458 # m/s.................... Speed of light in vacuum
h = 6.62607015e-34 # J Hz^-1.................. Planck's constant
# Other constants
G = 6.67430e-11 # m^3 kg^-1 s^-2........ Constant of gravitation
μ_0 = 1.25663706212e-6 # N A^-2.......Vacuum magnetic permeability
ϵ_0 = 8.8541878128e-12 # F m^-1.......Vacuum electric permittivity
# Masses
m_e = 9.1093837015e-31 # Kg ........................ Electron mass
m_p = 1.67262192369e-27 # Kg ......................... Proton mass

# CGS units
k_B_cgs = 1.380649e-16 # erg/K
e_cgs = 4.8032047e-10 # statcoulombs
m_e_cgs = m_e * 1e3 # grams
m_p_cgs = m_p * 1e3 # grams

# Derived quantities
csqrd = c^2 # s^2 m^-2 ........ Inverse of the light speed squared
csqrdinv = 1 / csqrd # s^2 m^-2 ........ Inverse of the light speed squared

# Unit conversion
J2eV = 6.24150907e18 # 1 J * 6.24e12 eV/J........................energy
