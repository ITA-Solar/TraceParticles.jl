#-------------------------------------------------------------------------------
# Created 02.08.23
# Author: e.s.oyre@astro.uio.no
#-------------------------------------------------------------------------------
#
#                 tp.jl
#
#-------------------------------------------------------------------------------
#
#-------------------------------------------------------------------------------

module TestParticles

#---------------#
#  Use packages #
#---------------#---------------------------------------------------------------
using LinearAlgebra
using Random
using Dates

using DifferentialEquations
using ForwardDiff
using FiniteDifferences
using Interpolations
using Optimization, OptimizationNLopt
using DataFrames
using CSV
using StatsBase
using HDF5

using Printf
using JLD2

using Distributed

using BifrostTools


#-------------------------------------------------------------------------------
"""
        Abstract types
"""
abstract type AbstractMesh end


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
include("io.jl")
include("bifrost.jl")
include("get.jl")
include("dataprocessing.jl")
include("gcastate.jl")
include("runs.jl")


#---------#
# Exports #
#---------#---------------------------------------------------------------------
# physics.jl
export
    get_guidingcentre,
    get_guidingcentre!,
    get_fullorbit,
    get_fullorbit!,
    cosineof_pitchangle,
    kineticenergy,
# equations_of_motion.jl
export
    guidingcentreapproximation!,
    gca_2Dxz!,
    lorentzforce!,
    lorentzforce,
    hybridgcafo!
# statistics.jl
export maxwellianvelocitysample
# callbacks.jl
export
    GCAInvalidCondition,
    GCAValidCondition,
    switch2fo_affect!,
    switch2gca_affect!
# bifrost.jl
export
    get_br_emfield_interpolator,
    get_br_emfield_vecof_interpolators,
    get_br_emfield_numdensity_gastemp_interpolator,
    get_br_var_interpolator
export
    h5_getall,
    h5_getbatch,
    h5_getdataset,
    h5_getenergies
# utils.jl
export
    dropdims,
    create3Daxes,
    discretise,
    UserData,
    MaxValue,
    MinValue
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
    get_exbdrift,
    get_fermi,
    get_betatron,
    get_polarisationacc,
    get_driftenergy,
    get_perpenergy,
    get_parallelenergy
# statistics.jl
export
    binmap,
    rejectionsample
# paramerterstructs.jl
export HybridParams
# runs.jl
export
    rerun_nonthermals,
    rerun_maxiters,
    rerun_idxs

#-------------------------------------------------------------------------------
end
