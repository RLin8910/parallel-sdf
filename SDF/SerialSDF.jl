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
Compute the unsigned distance field for a 2D matrix of booleans with Dijkstra's shortest path algorithm.

True values are considered to be inside of the region, while false values are considered outside.
"""
function dijkstraUDF2D(img:: Union{Matrix{Bool}, BitMatrix}, invert:: Bool = false):: Matrix{Float64}
    dirs = SVector{8, SVector{2, Int}}([[-1;-1], [-1; 0], [-1;1], [0;1], [1;1], [1;0], [1;-1], [0,-1]])
    
    # initiale data structures
    result = zeros(size(img))
    closed = zeros(Bool, size(img))
    pq = PriorityQueue{DijkstraPixel, Float64}(Base.Order.Forward)
    
    # seeding
    for x in 1:size(img,1)
        for y in 1:size(img,2)
            if img[x,y] != invert
                for dir in dirs
                    x1 = x + dir[1]
                    y1 = y + dir[2]
                    if x1 < 1 || x1 > size(img,1)
                        continue
                    end
                    if y1 < 1 || y1 > size(img,2)
                        continue
                    end
                    if img[x1, y1] != invert
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
            if closed[x1, y1] || img[x1, y1] != invert
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
Compute the signed distance field for a 2D matrix of booleans with Dijkstra's shortest path algorithm.

True values are considered to be inside of the region, while false values are considered outside.
"""
function dijkstraSDF2D(img:: Union{Matrix{Bool}, BitMatrix}):: Matrix{Float64}
    # compute unsigned distance from black pixels, then subtract unsigned distance from white
    # to get SDF
    return dijkstraUDF2D(img) - dijkstraUDF2D(img, true)
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
    for x in 1:size(img,1)
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

    closestX = Array{Int64}(undef, size(img,1))
    endpts = Array{Int64}(undef, size(img,1))

    # vertical pass - make use of previous best computed horizontal distances to pair with best vertical distance in similar
    # 2-pass approach
    for y in 1:size(img,2)
        seg = 1
        closestX[1] = 1
        endpts[1] = 1
        # nearest distance pass from left

        # partition into regions using closestX and endpts which indicate the closest set of coordinates
        # endpts indicates the endpoints of these partitions, while closestX indicates the x coordinates of the point
        for x in 2:size(img, 1)
            # create partition
            while seg > 0 && pix_dist(endpts[seg], y, closestX[seg]) > pix_dist(endpts[seg], y, x)
                seg -= 1
            end
            if seg < 1
                # first partition closest coordinate
                seg = 1
                closestX[1] = x
            else
                # consequtive partition endpoints and closest coordinates
                endpt = 1+sep(closestX[seg], x, y)
                if endpt <= size(img, 1)
                    seg += 1
                    closestX[seg] = x
                    endpts[seg] = endpt
                end
            end
        end

        # nearest distance pass from right - make use of endpoints
        for x in size(img,1):-1:1
            # Use different distance metric here - calculate based on nearest edge, and take square root
            diffx = x == closestX[seg] ? 0 : abs(x-closestX[seg]) - 0.5
            diffy = g[closestX[seg], y] == 0 ? 0 : g[closestX[seg],y] - 0.5

            result[x,y] = sqrt(diffx^2 + diffy^2)
            # Decrement the segment number
            if x == endpts[seg]
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
export dijkstraUDF2D
export dijkstraSDF2D
export linearUDF2D
export linearSDF2D
end # SerialSDF