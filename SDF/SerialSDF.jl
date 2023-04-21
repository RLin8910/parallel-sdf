module SerialSDF
using StaticArrays
using DataStructures
using Colors 
include("./SDFStructs.jl")
using .SDFStructs

"""
Convert a full-color image to a thresholded image of booleans only.
"""
function threshold(img:: Matrix{<:Colorant}, threshold::Float64=0.5, channel::Function = ColorTypes.red):: BitMatrix
	result = zeros(Bool, size(img))
	for x in 1:size(img,1)
        for y in 1:size(img, 2)
			result[x,y] = channel(img[x,y]) >= threshold
		end
	end
	return result
end


"""
Compute the signed distance field for a 2D matrix of booleans with a brute force approach.

True values are considered to be inside of the region, while false values are considered outside.
"""
function bruteSDF2D(img:: Union{Matrix{Bool}, BitMatrix}):: Matrix{Float64}
    result = zeros(size(img))
    for x in 1:size(img,1)
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
Compute the unsigned distance field for a 2D matrix of booleans with Dijkstra's algorithm.

True values are considered to be inside of the region, while false values are considered outside.
"""
function dijkstraUDF2D(img:: Union{Matrix{Bool}, BitMatrix}):: Matrix{Float64}
    dirs = SVector{8, SVector{2, Int}}([[-1;-1], [-1; 0], [-1;1], [0;1], [1;1], [1;0], [1;-1], [0,-1]])
    
    # initiale data structures
    result = zeros(size(img))
    closed = zeros(Bool, size(img))
    pq = PriorityQueue{DijkstraPixel, Float64}(Base.Order.Forward)
    
    # seeding
    for x in 1:size(img,1)
        for y in 1:size(img,2)
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
                    pq[DijkstraPixel(x1, y1, (x-x1) / 2, (y-y1) / 2)] = dist
                end
            end
        end
    end
    
    #propagation
    while length(pq) > 0
        pair = dequeue_pair!(pq)
        current = pair.first
        if closed[current.x, current.y]
            continue
        end
        
        closed[current.x, current.y] = true
        result[current.x, current.y] = pair.second
        
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
            
            pq[DijkstraPixel(x1, y1, dx1, dy1)] = sqrt(dx1^2 + dy1^2)
        end
    end

    return result
end

"""
Compute the signed distance field for a 2D matrix of booleans with Dijkstra's algorithm.

True values are considered to be inside of the region, while false values are considered outside.
"""
function dijkstraSDF2D(img:: Union{Matrix{Bool}, BitMatrix}):: Matrix{Float64}
    # compute unsigned distance from black pixels, then subtract unsigned distance from white
    # to get SDF
    return dijkstraUDF2D(img) - dijkstraUDF2D(.!img)
end

export bruteSDF2D
export dijkstraUDF2D
export dijkstraSDF2D
end # SerialSDF