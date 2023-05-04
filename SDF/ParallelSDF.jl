module ParallelSDF
using StaticArrays
using DataStructures
using Distributed

include("./SDFStructs.jl")
include("./SerialSDF.jl")
using .SDFStructs
using .SerialSDF

@views function makechunks(X::AbstractVector, n::Integer)
    c = length(X) ÷ n
    return [X[1+c*k:(k == n-1 ? end : c*k+c)] for k = 0:n-1]
end

"""
Compute the signed distance field for a 2D matrix of booleans with a brute force parallelized approach.

True values are considered to be inside of the region, while false values are considered outside.
"""
function bruteSDF2D(img:: Union{Matrix{Bool}, BitMatrix}):: Matrix{Float64}
    result = zeros(size(img))

    chunks = makechunks(collect(1:size(img,1)), Threads.nthreads())

    @sync @Threads.threads for chunk in chunks
        for x in chunk
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
    end
    return result
end

"""
Compute the unsigned distance field for a 2D matrix of booleans with Dijkstra's algorithm.

True values are considered to be inside of the region, while false values are considered outside.
"""
function dijkstraUDF2DParallel(img:: Union{Matrix{Bool}, BitMatrix}):: Matrix{Float64}
    dirs = SVector{8, SVector{2, Int}}([[-1;-1], [-1; 0], [-1;1], [0;1], [1;1], [1;0], [1;-1], [0,-1]])
    
    # initiale data structures
    result = zeros(size(img))
    closed = zeros(Bool, size(img))

    # create partitions
    # priority queues for each partition
    threadCount = Threads.nthreads()
    pqs = Array{PriorityQueue}(undef, threadCount)

	for i in 1:size(pqs,1)
		pqs[i] = PriorityQueue{DijkstraPixel, Float64}(Base.Order.Forward)
	end

    centerX = size(img,1) / 2
    centerY = size(img,2) / 2

    partition = zeros(Int64, size(img))

    # seeding
    for x in 1:size(img,1)
        for y in 1:size(img,2)
            # determine which priority queue to place this pixel in
            # partition image by angle radially about center, scaled as if image were square
            angle = atan(y / centerY - 1, x / centerX - 1)
            partition[x,y] = floor(Int, threadCount * mod2pi(angle) / (2 * pi)) + 1

            if img[x,y]
                for dir in dirs
                    x1 = x + dir[1]
                    y1 = y + dir[2]
                    if x1 < 1 || x1 > size(img,1)
                        continue
                    end
                    if y1 < 1 || y1 > size(img,2)
                        continue
                    end
                    if img[x1, y1]
                        continue
                    end
                    # distance to center of edge
                    dist = 0.5
                    if abs(dir[1]) + abs(dir[2]) > 1
                        # distance to corner
                        dist = sqrt(0.5)
                    end
                    
                    pqs[partition[x,y]][DijkstraPixel(x1, y1, (x-x1) / 2, (y-y1) / 2)] = dist
                end
            end
        end
    end
    
    #propagation
    while true
        # determine closest node in parallel
        minDist = typemax(Float64)
        current = undef
        part = 0

        Threads.@threads for i in 1:threadCount
            if length(pqs[i]) == 0
                continue
            end
            pix, dist = peek(pqs[i])
            if dist < minDist
                minDist = dist
                current = pix
                part = i
            end
        end 

        if part == 0
            return result
        end

        delete!(pqs[part], current)

        if closed[current.x, current.y]
            continue
        end
        
        closed[current.x, current.y] = true
        result[current.x, current.y] = minDist
        
        for dir in dirs
            x1 = current.x + dir[1]
            y1 = current.y + dir[2]
            if x1 < 1 || x1 > size(img,1)
                continue
            end
            if y1 < 1 || y1 > size(img,2)
                continue
            end
            if closed[x1, y1] || img[x1, y1]
                continue
            end
            dx1 = current.dx - dir[1]
            dy1 = current.dy - dir[2]
            
            pixel = DijkstraPixel(x1, y1, dx1, dy1)
            
            pqs[partition[x1,y1]][DijkstraPixel(x1, y1, dx1, dy1)] = sqrt(dx1^2 + dy1^2)
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

"""
Compute the signed distance field for a 2D matrix of booleans with Dijkstra's algorithm.

True values are considered to be inside of the region, while false values are considered outside.
"""
function dijkstraSDF2DParallelUDF(img:: Union{Matrix{Bool}, BitMatrix}):: Matrix{Float64}
    # compute unsigned distance from black pixels, then subtract unsigned distance from white
    # to get SDF
    return dijkstraUDF2DParallel(img) - dijkstraUDF2DParallel(.!img)
end

export bruteSDF2D
export dijkstraUDF2DParallel
export dijkstraSDF2DSerialUDF
end #ParallelSDF