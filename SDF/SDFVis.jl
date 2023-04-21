module SDFVis
using Colors 
using Images

function normalizeSDF(img:: Matrix{Float64}):: Matrix{Float64}
    maxElt = maximum(img)
    minElt = minimum(img)
    if maxElt != minElt
        result = img ./ max(abs(maxElt), abs(minElt))
        return result
    end
    return copy(img)
end

function toImageSDF(img:: Matrix{Float64}, boundaryThreshold:: Float64 = 0.71):: Matrix{RGB{Float32}}
    normalized = normalizeSDF(img)
    result = zeros(RGB{Float32}, size(normalized))
    for x in 1:size(normalized, 1)
        for y in 1:size(normalized, 2)
            result[x,y] = RGB(max(0, normalized[x,y]), max(0, -normalized[x,y]), abs(img[x,y]) < boundaryThreshold)
        end
    end
    return result
end

export normalizeSDF
export toImageSDF
end # SDFVis