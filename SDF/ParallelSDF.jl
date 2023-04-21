module ParallelSDF
using StaticArrays
using DataStructures
include("./SDFStructs.jl")
using .SDFStructs

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

end #ParallelSDF