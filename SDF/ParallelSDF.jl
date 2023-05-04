module ParallelSDF
using StaticArrays

include("./SDFStructs.jl")
include("./SerialSDF.jl")
using .SDFStructs
using .SerialSDF

"""
Compute the signed distance field for a 2D matrix of booleans with a brute force approach.

True values are considered to be inside of the region, while false values are considered outside.
"""
function bruteSDF2D(img:: Union{Matrix{Bool}, BitMatrix}):: Matrix{Float64}
    result = zeros(size(img))
    @sync Threads.@threads for x in 1:size(img,1)
        for y in 1:size(img, 2)
            # if black, determine distance to closest white pixel.
            # if white, determine distance to closest black pixel.
            # white pixels are "interior", black pixels are "exterior"
            # interior is negative distance. exterior is positive distance
            bestDistance = Inf

            # determine best distance by brute force
            for x1 in 1:size(img,1)
                for y1 in 1:size(img, 2)
                    if img[x1, y1] != img[x, y]
                        # compute distance to boundary of nearest opposite pixel, not the center
                        dist = 0
                        if x1 == x || y1 == y
                            # pixels in same row/col, nearest boundary point is on the center of the edge
                            dist = abs(x1 - x) + abs(y1 - y) - 0.5
                        else
                            # pixels in different row and col, nearest boundary point is on a corner
                            dist = sqrt((abs(x1-x)-0.5)^2 + (abs(y1-y)-0.5)^2)
                        end
                        if dist < bestDistance
                            bestDistance = dist
                        end
                    end
                end
            end
            
            # invert if inside the region
            if img[x, y]
                result[x,y] = -bestDistance
            else
                result[x,y] = bestDistance
            end
        end
    end
    return result
end

"""
Compute the signed distance field for a 2D matrix of booleans with Dijkstra's algorithm.

True values are considered to be inside of the region, while false values are considered outside.
"""
function dijkstraSDF2DSerialUDF(img:: Union{Matrix{Bool}, BitMatrix}):: Matrix{Float64}
    # compute unsigned distance from black pixels, then subtract unsigned distance from white
    # to get SDF
    pos = undef
    neg = undef
    @sync begin
        Threads.@spawn pos = SerialSDF.dijkstraUDF2D(img) 
        Threads.@spawn neg = SerialSDF.dijkstraUDF2D(.!img)
    end
    return pos - neg
end

export bruteSDF2D
export dijkstraSDF2DSerialUDF
end #ParallelSDF