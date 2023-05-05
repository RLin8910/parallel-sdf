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
Compute the signed distance field for a 2D matrix of booleans with Dijkstra'closestX algorithm.

True values are considered to be inside of the region, while false values are considered outside.
"""
function dijkstraSDF2DSerialUDF(img:: Union{Matrix{Bool}, BitMatrix}):: Matrix{Float64}
    # compute unsigned distance from black pixels, then subtract unsigned distance from white
    # to get SDF
    pos = undef
    neg = undef
    @sync begin
        Threads.@spawn pos = SerialSDF.dijkstraUDF2D(img) 
        Threads.@spawn neg = SerialSDF.dijkstraUDF2D(img, true)
    end
    return pos - neg
end

"""
Compute the unsigned distance field for a 2D matrix of booleans with a linear-time distance transform.

True values are considered to be inside of the region, while false values are considered outside.
"""
function linearUDF2D(img:: Union{Matrix{Bool}, BitMatrix}, invert::Bool = false):: Matrix{Float64}
    max_val = size(img,1) + size(img,2) + 2 # "Infinity" value which is guaranteed to be larger than all distances
    result = zeros(size(img))
    g = zeros(size(img))

    # determine the horizontal UDF in the first pass, i.e. the distance to the closest white pixel
    # which is on the same row
    @sync Threads.@threads for x in 1:size(img,1)
        if img[x, 1] != invert
            g[x, 1] = 0
        else
            g[x, 1] = max_val
        end

        # nearest distance pass from below
        for y in 2:size(img,2)
            if img[x,y] != invert
                g[x,y] = 0
            else
                g[x,y] = 1 + g[x, y-1]
            end
        end
        
        # nearest distance pass from above
        for y in (size(img,2)-1):-1:1
            if g[x, y+1] < g[x,y]
                g[x,y] = 1 + g[x, y+1]
            end
        end
    end

    """
    Compute the best squared distance from the pixel at [x,y] to any white pixel on the ith column

    Use squared distance because the actual value is not relevant, only what is closer and further
    """
    function pix_dist(x, y, i)
        return (x-i)^2 + g[i,y]^2
    end

    """
    Compute the segment endpoint. 

    As the scan partitions the region horizontally based on closest x coordinate, `sep` computes the endpoint of the 
    segment. 
    """
    function sep(i, x, y)
        return floor(Int, (x^2 - i^2 + g[x,y]^2 - g[i,y]^2) / (2*(x-i)))
    end

    # Preallocate closestX and endpts arrays to number of threads
    closestX = zeros(Int64, (size(img,1), Threads.nthreads()))
    endpts = zeros(Int64, (size(img,1), Threads.nthreads()))

    # vertical pass - make use of previous best computed horizontal distances to pair with best vertical distance in similar
    # 2-pass approach
    @sync Threads.@threads for y in 1:size(img,2)
        seg = 1
        id = Threads.threadid()
        closestX[1, id] = 1
        endpts[1, id] = 1
        # nearest distance pass from left

        # partition into regions using closestX and endpts which indicate the closest set of coordinates
        # endpts indicates the endpoints of these partitions, while closestX indicates the x coordinates of the point
        
        for x in 2:size(img, 1)
            # create partition
            while seg > 0 && pix_dist(endpts[seg,id], y, closestX[seg,id]) > pix_dist(endpts[seg,id], y, x)
                seg -= 1
            end
            if seg < 1
                # first partition closest coordinate
                seg = 1
                closestX[1,id] = x
            else
                # consequtive partition endpoints and closest coordinates
                endpt = 1+sep(closestX[seg,id], x, y)
                if endpt <= size(img, 1)
                    seg += 1
                    closestX[seg,id] = x
                    endpts[seg,id] = endpt
                end
            end
        end

        # nearest distance pass from right - make use of endpoints
        for x in size(img,1):-1:1
            # Use different distance metric here - calculate based on nearest edge
            diffx = x == closestX[seg,id] ? 0 : abs(x-closestX[seg,id]) - 0.5
            diffy = g[closestX[seg,id], y] == 0 ? 0 : g[closestX[seg,id],y] - 0.5

            result[x,y] = sqrt(diffx^2 + diffy^2)
            # Decrement the segment number
            if x == endpts[seg,id]
                seg -= 1
            end
        end
    end
    
    return result
end

"""
Compute the signed distance field for a 2D matrix of booleans using a linear-time distance transform algorithm.

True values are considered to be inside of the region, while false values are considered outside.
"""
function linearSDF2D(img:: Union{Matrix{Bool}, BitMatrix}):: Matrix{Float64}
    return linearUDF2D(img) - linearUDF2D(img, true)
end

export bruteSDF2D
export dijkstraSDF2DSerialUDF
export linearUDF2D
export linearSDF2D
end #ParallelSDF