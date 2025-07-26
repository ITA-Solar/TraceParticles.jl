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


# utils.jl
export
    ElectromagneticFieldInterpolator
# runs.jl
export
    rerun
# equations_of_motion.jl
export
    guidingcentreapproximation!,
    gca_2Dxz!,
    lorentzforce!
# statistics.jl
export
    maxwellianvelocitysample,
    binmap
# physics.jl
export
    get_fullorbit,
    get_fullorbit!,
    get_guidingcentre!,
    kineticenergy,
    perpendicular_velocity,
    larmorradius,
    characteristicfieldlength,
    scalesratio,
    magneticmoment
# callbacks.jl
export
    OutOfDomainCondition_3D,
    OutOfDomainCondition_2Dxz,
    RelativisticConditionGCA,
    RelativisticConditionGCA_2Dxz,
    GCABreakDownCondition,
    GCABreakDownCondition_2Dxz,
    outofdomainaffect!,
    relativisticaffect!,
    gcabreakdownaffect!
# problem_functions.jl
export
    UposMBvel,
    DposMBvel,
    DposMBvel_2Dxz,
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
    h5_getfinalstate,
    save_gcastates,
    create_bifrost_itps,
    XZInterpolation
# dataprosessing.jl
export
    initialstate_retmsg,
    finalstate_retmsg,
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
    get_observable,
    get_stateidx,
    get_exbdrift,
    get_fermi,
    get_betatron,
    get_polarisationacc,
    get_driftenergy,
    get_perpenergy,
    get_parallelenergy
# constants.jl
export
    J2eV

#-------------------------------------------------------------------------------
end
