# Created 05.01.24 by eilfso
# Author: e.s.oyre@astro.uio.no
#
#           mathematics.jl
#-------------------------------------------------------------------------------
# Contains calculus
#-------------------------------------------------------------------------------


"""
    skewsymmetric_matrix(vector::Vector)
Returns the skew-symmetric matrix of a 3D vector. The matrix is used to 
represent a cross product between two 3D vectors as a matrix product between 
the skew-symmetric matrix and the other vector.
"""
function skewsymmetric_matrix(vector::Vector{Float64})
    return [0.0 -vector[3] vector[2]
        vector[3] 0.0 -vector[1]
        -vector[2] vector[1] 0.0]
end
function skewsymmetric_matrix(vector::Vector{Float32})
    return [0.0f0 -vector[3] vector[2]
        vector[3] 0.0f0 -vector[1]
        -vector[2] vector[1] 0.0f0]
end
