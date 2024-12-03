include("../Utilities.jl")
include("../LoRa_Encoding/ZeroPadding.jl")
#include("../config/parameters.jl")
using Base.Iterators: take, cycle


#= Debugging
bit_M = header_M
CR = 4
DE =  true
isheader = true
SF
=#




function diag_deinterleaving(bit_M::BitMatrix, CR::Int8, DE::Bool)

    if(DE == true)
        SF_LDRO = SF - 2;
    else 
        SF_LDRO = SF;
    end

    # create matrix with size CR_payload + 4 x SF_LDRO
    deinterleaved_matrix = BitMatrix(falses(SF_LDRO, 8)) 

    for x in 1:(CR + 4)
        deinterleaved_matrix[:, 8 - (CR + 4) + x] = circshift(reverse(bit_M[x,1:SF_LDRO]) , x-1)
    end

    return deinterleaved_matrix
end



