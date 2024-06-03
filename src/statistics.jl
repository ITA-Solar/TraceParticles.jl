# Created 22.01.24
# Author: e.s.oyre@astro.uio.no
#
#           statistics.jl
#-------------------------------------------------------------------------------
# Contains statistical distributions, random number generators and samplers.
#-------------------------------------------------------------------------------


#---------------#
# Distributions #
#---------------#---------------------------------------------------------------
"""
    normaldistr(
        x::Vector{T} where {T<:Real},
        μ::Real,
        σ::Real
        )

Function returning the values of `x` on a 1D normalised normal distribution with
expectation value `μ` and standard deviation `σ`.

See also [`Utilities.uniformdistr`](@ref).
"""
function normaldistr(
    x::Array{T} where {T<:Real},
    μ::Real,
    σ::Real
    )
    return @.  1/(σ*√(2π))*exp(-0.5((x - μ)/σ)^2)
end # normaldistr
function normaldistr(
    x::Real,
    μ::Real,
    σ::Real
    )
    return 1/(σ*√(2π))*exp(-0.5((x - μ)/σ)^2)
end

function bivariate_normaldistr(
    x::Real,
    y::Real,
    μx::Real,
    μy::Real,
    σx::Real,
    σy::Real,
    corr::Real # correlation
    )
    normfactor = 1/(2π*σx*σy*√(1 - corr^2))
    return normfactor * exp(
        -1/(2*(1 - corr^2)) * (
            (x - μx)^2/σx^2 +
            (y - μy)^2/σy^2 -
            2*corr*(x - μx)*(y - μy)/(σx*σy)
            )
        )
end

"""
    uniformdistr(
        x::Vector{T} where {T<:Real},
        a::Real,
        b::Real
        )
Function returning the values of `x` on a 1D normalised unifrom distribution on
the interval [`a, `b`].

See also [`Utilities.normaldistr`](@ref).
"""
function uniformdistr(
    x::Array{T} where {T<:Real},
    a::Real,
    b::Real
    )
    mask = a .<= x .<= b
    stepheight = 1.0/(b - a)
    prob = zeros(typeof(x[1]), size(x))
    prob[mask] .= stepheight
    return prob
end # function uniformdistr


"""
    bivarate_uniformdistr(x,y,a,b,c,d)
Return the probability of (`x`,`y`) in the bivariate uniform distribution
f(x,y) = U(a,b) x U(c,d)
"""
function bivariate_uniformdistr(
    x::Real,
    y::Real
    ;
    a::Real=0.0,
    b::Real=1.0,
    c::Real=0.0,
    d::Real=1.0
    )
    Δx = b - a
    Δy = d - c
    probability = abs(1/(Δx*Δy))
    return a < x <= b && c < y <= d ? probability : 0.0
end


#----------------#
# Random numbers #
#----------------#--------------------------------------------------------------
"""
    Base.randn(μ, σ, dims)
    
Method which return random variables from a normal distribution with 
expectation value `μ` and standard deviation `σ`.
"""
function Base.randn(
    μ   ::AbstractFloat,
    σ   ::AbstractFloat,
    dims...
    )
    return μ .+ σ.*randn(dims)
end
function Base.randn(
    precision::DataType,
    μ        ::AbstractFloat,
    σ        ::AbstractFloat,
    dims...
    )
    μ = convert.(precision, μ)
    σ = convert.(precision, σ)
    return μ .+ σ.*randn(precision, dims)
end 
function Base.randn(
    rng      ::AbstractRNG,
    precision::DataType,
    μ        ::AbstractFloat,
    σ        ::AbstractFloat,
    dims...
    )
    μ = convert.(precision, μ)
    σ = convert.(precision, σ)
    return μ .+ σ.*randn(rng, precision, dims)
end


"""
    Base.rand(rng, precision, a, b, dims)
Method which return random variables from a uniform distribution on the interval
[a, b] with DataType `precision` and random seed generator `rng`.
"""
function Base.rand(
    a   ::AbstractFloat,
    b   ::AbstractFloat,
    dims...
    )
    return a .+ rand(dims...) .* (b - a)
end 
function Base.rand(
    precision::DataType,
    a        ::AbstractFloat,
    b        ::AbstractFloat,
    dims...
    )
    a = convert(precision, a)
    b = convert(precision, b)
    return a .+ rand(precision, dims...) .* (b - a)
end 
function Base.rand(
    rng ::AbstractRNG,
    a   ::AbstractFloat,
    b   ::AbstractFloat,
    dims...
    )
    return a .+ rand(rng, dims...) .* (b - a)
end
function Base.rand(
    rng      ::AbstractRNG,
    precision::DataType,
    a        ::AbstractFloat,
    b        ::AbstractFloat,
    dims...
    )
    a = convert(precision, a)
    b = convert(precision, b)
    return a .+ rand(rng, precision, dims...) .* (b - a)
end
function Base.rand(
    a   ::AbstractFloat,
    b   ::AbstractFloat,
    dims::Tuple{Vararg{Int}}
    )
    return a .+ rand(dims...) .* (b - a)
end 
function Base.rand(
    precision::DataType,
    a        ::AbstractFloat,
    b        ::AbstractFloat,
    dims::Tuple{Vararg{Int}}
    )
    a = convert(precision, a)
    b = convert(precision, b)
    return a .+ rand(precision, dims...) .* (b - a)
end 
function Base.rand(
    rng ::AbstractRNG,
    a   ::AbstractFloat,
    b   ::AbstractFloat,
    dims::Tuple{Vararg{Int}}
    )
    return a .+ rand(rng, dims...) .* (b - a)
end 
function Base.rand(
    rng      ::AbstractRNG,
    precision::DataType,
    a        ::AbstractFloat,
    b        ::AbstractFloat,
    dims::Tuple{Vararg{Int}}
    )
    a = convert(precision, a)
    b = convert(precision, b)
    return a .+ rand(rng, precision, dims...) .* (b - a)
end
function Base.rand(
    domain::Matrix{T} where {T<:Real}
    )
    numaxes = size(domain)[1]
    r = zeros(numaxes)
    for i = 1:numaxes
        a = domain[i,1]
        b = domain[i,2]
        r[i] = a .+ rand() .* (b - a)
    end
    return  r
end # function rand


#----------#
# Sampling #
#----------#--------------------------------------------------------------------
function maxwellianvelocitysample(
    rng            ::AbstractRNG,
    temperature_itp::Any,
    mass           ::Real,
    args...
    ;
    precision::DataType=Float64,
    )
    T = temperature_itp.(args...)
    μ = 0.0
    σ = sqrt.(tp.k_B*T/mass) # Standard deviation of the Maxwell
    # distribution at this temperature
    return randn.(rng, precision, μ, σ, 3)
end


"""
    importancesampling(
        target  ::Function, 
        proposal::Function, 
        randgen ::Function, 
        N       ::Integer,    
        )
Sample `N` points from the `proposal`-distribution and compute the importance
weights with respect to the `target`-distribution.
"""
function importancesampling(
    target  ::Function, # Target distribution
    proposal::Function, # Proposal distruv
    randgen ::Function, # Random variable generator. Following proposal pdf.
    dims    ::Tuple{Vararg{Integer}}, # Number of samples
    )
    samples = randgen(dims)
    weights = target(samples) ./ proposal(samples)
    return samples, weights
end # function importancesampling


function rejectionsampling(
    target    ::Function,
    maxvalue  ::Real,
    numsamples::Integer,
    domain    ::Matrix{T} where {T<:Real}
    )
    numdims = size(domain)[1]
    positions = zeros((numdims, numsamples))
    accepted = 0
    rejected = 0
    while accepted < numsamples
        pos     = rand(domain)
        yguess  = rand(0.0, maxvalue, 1)
        ytarget = target(pos)[1] # Target should return a single float, but in
        # case this float is in a 1-element Vector we specify the first index.
        if yguess < ytarget
            accepted += 1
            positions[:,accepted] .= pos
        else
            rejected += 1
        end
    end
    acceptencerate = accepted/rejected
    if numdims == 1
        return positions[1,:], acceptencerate
    else
        return positions, acceptencerate
    end
end # function rejection sampling


