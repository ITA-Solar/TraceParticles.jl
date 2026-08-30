#-------------------------------------------------------------------------------
# Created 02.08.23
# Author: e.s.oyre@astro.uio.no
#-------------------------------------------------------------------------------

module TraceParticles

#________/\_____________________________________________________________________
#
#  Dependencies
#
#____/\______/\_________________________________________________________________

using BifrostTools
using DataFrames
using Dates
using DiffEqCallbacks
using Distributed
using EnumX
using ForwardDiff
using HDF5
using Interpolations
using JLD2
using LinearAlgebra
using Logging
using LoggingExtras
using OrdinaryDiffEq
using Printf
using Random
using SciMLBase
using StaticArrays
using Statistics
using StatsBase


#_______________________________________________________________________________

include("enums.jl")
include("utils.jl")
include("constants.jl")
include("statistics.jl")
include("physics.jl")
include("equations_of_motion.jl")
include("callbacks.jl")
include("problem_functions.jl")
include("callback_specifiers.jl")
include("output_functions.jl")
include("reduction_functions.jl")
include("odeparameters.jl")
include("bifrost.jl")
include("io.jl")
include("dataprocessing.jl")
include("gcastate.jl")
include("get.jl")
include("traceparticlesparameters.jl")
include("traceparticlesproblem.jl")
include("solve.jl")
include("rerun.jl")


#________/\_____________________________________________________________________
#
# Exports
#
#____/\______/\_________________________________________________________________


# equations_of_motion.jl
export
    guidingcentreapproximation!,
    lorentzforce!,
    hybridgcafo!
# enums.jl
export EoMID
# problem_functions.jl
export
    SampleICsFromMHD,
    ICsFromFile,
    PredefinedICs,
    SampleFullOrbit,
    SampleGCA,
    SampleHybrid
# output_functions.jl
export
    output_func_lightweight_gca,
    output_func_lightweight_fo,
    output_func_lightweight_hybrid,
    output_func_max_lightweight
# reduction_functions.jl
export
    SaveBatchAsHDF5,
    get_filename
# odeparameters.jl
export
    FullOrbitParams,
    GCAParams,
    HybridParams,
    HybridParamsWithDetection
# io.jl
export
    h5_getbatchnames,
    h5_getall,
    h5_getbatch,
    h5_getdataset,
    h5_getenergies
# rerun.jl
export
    rerun
# callback_specifiers.jl
export
    OutOfBounds,
    Relativistic,
    HighMagneticGradient
# traceparticlesparameters.jl
export
    TraceParticlesParameters
# traceparticlesproblem.jl
export
    TraceParticlesProblem
# dataprocessing.jl
export
    GCAState,
    get_observable
# physics
export
    get_fullorbit,
    get_guidingcentre
# constants
export
    electron_mass,
    electron_charge,
    proton_mass,
    proton_charge
#-------------------------------------------------------------------------------
end
