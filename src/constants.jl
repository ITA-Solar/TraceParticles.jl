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



k_B      = 1.380649e-23 # J/K..................... Boltzmann's constant
e        = 1.602176634e-19 # C....................... Elementary charge
c        = 299792458 # m/s.................... Speed of light in vacuum
G        = 6.67430e−11 # m^3 kg^-1 s^-2........ Constant of gravitation
h        = 6.62607015e−34 # J Hz^-1.................. Planck's constant
μ_0      = 1.25663706212e−6 # N A^-2.......Vacuum magnetic permeability
ϵ_0      = 8.8541878128e−12 # F m^-1.......Vacuum electric permittivity
m_e      = 9.1093837015e-31 # Kg ........................ Electron mass
m_p      = 1.67262192369e-27 # Kg ......................... Proton mass
cSqrdInv = 1/c^2 # s^2 m^-2 ........ Inverse of the light speed squared

# CGS-units
k_B_cgs = 1.380649e-16 # erg/K
e_cgs_esu = 4.80320427e-10 # Fr
m_e_cgs = 9.10938370e-28 # g

# Unit conversion
cgs2SI_t = 1e0  # s * 1 s/s ..........................cgs-seconds to SI
cgs2SI_tg = 1e0 # K * 1 K/K ................................temperature
cgs2SI_u = 1e-2 # cm/s * 1e-2 m/cm = m/s.......................velocity
cgs2SI_b = 1e-4 # G * 1e-4 T/G = T.......................magnetic field
cgs2SI_l = 1e-2 # cm * 1e-2 m/cm = m.............................length
cgs2SI_e = 1e-6 # g cm s^-2 Fr^-1 * 1e-3 kg/g * 1e-2 m/cm * 1e-1 Fr/C
                # = kg m /s^2/C..........................electric field
cgs2SI_r = 1e3  # g/cm^3 * 1e-3 kg/g * 1e6 cm^3/m^3 = kg/m^3 ...density

J2eV = 6.24150907e18 # 1 J * 6.24e12 eV/J........................energy

# Other convertion factors
density2numberdensity_MHD = 1/(m_e + m_p)


