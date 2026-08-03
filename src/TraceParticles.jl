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
using Parameters
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
# constants.jl
export
    J2eV
# equations_of_motion.jl
export
    guidingcentreapproximation!,
    lorentzforce!,
    hybridgcafo!
# statistics.jl
export
    rejectionsample,
    maxwellianvelocitysample,
    binmap
# physics.jl
export
    get_fullorbit,
    get_fullorbit!,
    get_guidingcentre,
    get_guidingcentre!,
    kineticenergy,
    perpendicular_velocity,
    larmorradius,
    characteristicfieldlength,
    scalesratio,
    magneticmoment
# callbacks.jl
export
    OutOfBoundsCondition,
    MagneticGradientCondition,
    MagneticCurvatureCondition,
    ParallelElectricFieldCondition,
    RelativisticConditionGCA,
    HybridSwitchCondition,
    outofboundsaffect!,
    magneticgradientaffect!,
    magneticcurvatureaffect!,
    parallelelectricfieldaffect!,
    relativisticaffect!,
    hybridswitchaffect!,
    hybridswitchaffect_withdetection!,
    TerminationCode
# problem_functions.jl
export
    MHDSamplingGCA,
    MHDSamplingHybrid,
    PredefinedICs,
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
    print_batch_statistics
# odeparameters.jl
export
    GCAParams,
    ParticleParams,
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
    save_gcastates,
    save_energy,
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
# gcastate.jl
export
    GCAState
# get.jl
export
    get_observable
# traceparticlesparameters.jl
export
    TraceParticlesParameters
# traceparticlesproblem.jl
export
    TraceParticlesProblem
#-------------------------------------------------------------------------------
end
