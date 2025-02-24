#-------------------------------------------------------------------------------
# Created 02.08.23
# Author: e.s.oyre@astro.uio.no
#-------------------------------------------------------------------------------

module TestParticles

#________/\_____________________________________________________________________
#
#  Dependencies
#
#____/\______/\_________________________________________________________________

# Mathematics
using LinearAlgebra
using DifferentialEquations
using ForwardDiff
using Interpolations
# Statistics
using Random
using StatsBase
# Data representation
using DataFrames
# I/O
using CSV
using HDF5
using Printf
using JLD2
using BifrostTools
# Parallel computing
using Distributed
# Logging
using Dates


#_______________________________________________________________________________

include("utils.jl")
include("constants.jl")
include("mathematics.jl")
include("statistics.jl")
include("physics.jl")
include("equations_of_motion.jl")
include("callbacks.jl")
include("problem_functions.jl")
include("output_functions.jl")
include("reduction_functions.jl")
include("parameterstructs.jl")
include("differentiations.jl")
include("bifrost.jl")
include("io.jl")
include("dataprocessing.jl")
include("gcastate.jl")
include("get.jl")
include("runs.jl")


#________/\_____________________________________________________________________
#
# Exports
#
#____/\______/\_________________________________________________________________


# runs.jl
export
    rerun_nonthermals,
    rerun_maxiters,
    rerun_idxs
# equations_of_motion.jl
export
    guidingcentreapproximation!,
    gca_2Dxz!,
    lorentzforce!
# statistics.jl
export maxwellianvelocitysample
# callbacks.jl
export
    OutOfDomainCondition_2Dxz,
    outofdomainaffect!,
    GCABreakDownCB
# problem_functions.jl
export
    UposMBvel,
    DposMBvel,
    PposPvel
# output_functions.jl
export
    output_func_max_lightweight
# reduction_functions.jl
export
    SaveBatchAsCSV,
    SaveBatchAsHDF5,
    get_filename
# io.jl
export
    h5_getall,
    h5_getbatch,
    h5_getdataset,
    h5_getenergies,
    h5_getinitialstate,
    save_gcastates
# gcastate.jl
export
    GCAState,
    get_drifts,
    get_x,
    get_y,
    get_z,
    get_vparal,
    get_mu,
    get_time,
    get_energy,
    get_vperp,
    get_larmorradius,
    get_gyrofrequency,
    get_bfield,
    get_efield,
    get_fieldlength,
    get_scalesratio
# get.jl
export
    get_observable,
    get_stateidx,
    get_exbdrift,
    get_fermi,
    get_betatron,
    get_polarisationacc,
    get_driftenergy,
    get_perpenergy,
    get_parallelenergy
# statistics.jl
export
    binmap

#-------------------------------------------------------------------------------
end
