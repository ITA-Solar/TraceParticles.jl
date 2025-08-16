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
using CSV
using DataFrames
using Dates
using DifferentialEquations
using Distributed
using EnumX
using ForwardDiff
using HDF5
using Interpolations
using JLD2
using LinearAlgebra
using Printf
using Random
using StaticArrays
using StatsBase


#_______________________________________________________________________________

include("utils.jl")
include("constants.jl")
include("statistics.jl")
include("physics.jl")
include("equations_of_motion.jl")
include("callbacks.jl")
include("problem_functions.jl")
include("output_functions.jl")
include("reduction_functions.jl")
include("parameterstructs.jl")
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
    switchaffect!,
    TerminationCode
# problem_functions.jl
export
    DposMBvel,
    PposPvel
# output_functions.jl
export
    output_func_lightweight,
    output_func_max_lightweight
# reduction_functions.jl
export
    SaveBatchAsHDF5,
    get_filename
# parameterstructs.jl
export
    GCAParams,
    ParticleParams,
    HybridParams
# io.jl
export
    h5_getall,
    h5_getbatch,
    h5_getdataset,
    h5_getenergies,
    h5_getinitialstate,
    h5_getfinalstate,
    save_gcastates,
    create_bifrost_itps
# dataprosessing.jl
export
    initialstate_terminationcode,
    finalstate_terminationcode,
    initialstate_nonthermal,
    initialstate_maxiters,
    initialstate_idxs,
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
# runs.jl
export
    rerun
#-------------------------------------------------------------------------------
end
