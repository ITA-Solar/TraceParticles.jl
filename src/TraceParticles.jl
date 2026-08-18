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


# utils.jl
export
    ElectromagneticFieldInterpolator,
    StaticInterpolation,
    XZInterpolation,
    StaticXZInterpolation
# equations_of_motion.jl
export
    guidingcentreapproximation!,
    lorentzforce!,
    hybridgcafo!
# statistics.jl
export
    binmap
# enums.jl
export
    TerminationCode,
    EoMID
# problem_functions.jl
export
    SampleICsFromMHD,
    ICsFromFile,
    PredefinedICs,
    SampleFullOrbit,
    SampleGCA,
    SampleHybrid,
    initialconditions_mhdsampling
# output_functions.jl
export
    output_func_lightweight,
    output_func_lightweight_hybrid,
    output_func_max_lightweight
# reduction_functions.jl
export
    SaveBatchAsHDF5,
    get_filename,
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
    h5_getenergies,
    h5_getenergies_gca,
    h5_getinitialstate,
    h5_getfinalstate,
    h5_getobservables,
    create_bifrost_itps,
    create_diffeq_ic
# dataprocessing.jl
export
    initialstate_terminationcode,
    finalstate_terminationcode,
    initialstate_nonthermals,
    initialstate_maxiters,
    initialstate_idxs,
    finalstate_idxs,
    particlemaskfunction,
    particlemask,
    replaceparticles,
    generate_probdistr
# rerun.jl
export
    rerun,
# gcastate.jl
export
    GCAState
# get.jl
export
    get_observable
# callback_specifiers.jl
export
    CallbackSpec,
    create_callback,
    OutOfBounds,
    Relativistic,
    HighMagneticGradient
# traceparticlesparameters.jl
export
    TraceParticlesParameters
# traceparticlesproblem.jl
export
    TraceParticlesProblem
#-------------------------------------------------------------------------------
end
