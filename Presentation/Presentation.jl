### A Pluto.jl notebook ###
# v0.19.22

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
end

# ╔═╡ fbec5760-dfb8-11ed-19ce-8753ae969cde
begin
	import Pkg; 
	Pkg.add("StaticArrays")
	Pkg.add("FileIO")
	using Images, FileIO
end

# ╔═╡ 3ef1d920-6a8b-45cd-afc0-4ac1a450a952
begin
	# Parallel SDF dependencies
	Pkg.add("MPI")
	using MPI
end

# ╔═╡ b8ff7564-a124-4605-8506-56b7ceeba8f8
begin
	# my code
	include("../SDF/SerialSDF.jl");
	include("../SDF/SDFVis.jl");
	using .SerialSDF
	using .SDFVis
end

# ╔═╡ 39f9509e-cca3-401e-b0c0-15dcea99d03e
begin
	include("../SDF/ParallelSDF.jl")
	using .ParallelSDF
end

# ╔═╡ e02b0513-471d-43e4-8e3d-510fac112c4d
begin
	# Generated Images
	# create just a couple of random
	dots = rand(Float64,256,256) .< 0.001
	# circle shape
	circle = zeros(Bool, 256,256)
	
	for x in 1:256
	    for y in 1:256
	        dist = sqrt((x-128.5)^2 + (y-128.5)^2)
	        circle[x,y] = dist ≤ 50
	    end
	end
end

# ╔═╡ 5baa54db-95d8-4c71-a991-1b145128d788
begin
	using Colors
	Gray.(dots)
end

# ╔═╡ b3b9540e-fc92-43b6-b826-a74dbc6e325a
Gray.(circle)

# ╔═╡ 6bf774a8-3990-4278-a011-34c7c450fc55
petty = load("./Images/43.png")

# ╔═╡ 6ffc2c09-39a4-4a0e-875d-8224e1ff119f
amogus = load("./Images/amogus.png")

# ╔═╡ e0f93c0e-28f3-4be5-815b-0180cb87c60a
begin
	petty_threshold = SerialSDF.threshold(petty)
	Gray.(petty_threshold)
end

# ╔═╡ f64b95bd-1fbd-42cd-87c2-cd42dba63ca6
begin
	amogus_threshold = SerialSDF.threshold(amogus, 0.5, ColorTypes.alpha)
	Gray.(amogus_threshold)
end

# ╔═╡ 84303126-2191-4bcc-8aa5-7b05ef62d57c
begin
	@time dotsSDFBrute = SerialSDF.bruteSDF2D(dots)
	SDFVis.toImageSDF(dotsSDFBrute)
end

# ╔═╡ 54c6a8ed-03d3-4f71-b240-498782fb0bda
begin
	@time dotsSDFDijkstra = SerialSDF.dijkstraSDF2D(dots)
	SDFVis.toImageSDF(dotsSDFDijkstra)
end

# ╔═╡ a2a42489-278e-45f1-a872-9af4068adbd9
begin
	@time circleSDFBrute = SerialSDF.bruteSDF2D(circle)
	SDFVis.toImageSDF(circleSDFBrute)
end

# ╔═╡ 32c138ab-f116-4127-824b-527520a206fc
begin
	@time circleSDFDijkstra = SerialSDF.dijkstraSDF2D(circle)
	SDFVis.toImageSDF(circleSDFDijkstra)
end

# ╔═╡ db3af62f-2e1a-4e54-ad9e-57bef579e609
begin
	@time amogusSDF = SerialSDF.dijkstraSDF2D(amogus_threshold)
	SDFVis.toImageSDF(amogusSDF, 2.0)
end

# ╔═╡ aa9bf6fe-a58e-4537-82f0-0eee4509c591
begin
	amogus_thres_slider = @bind amogus_thres_val html"<input type=range min=-500 max=500>"
	md"Threshold: $amogus_thres_slider "
end

# ╔═╡ cc3a2409-bdc7-4eb3-b0c3-a9dec01cf69d
Gray.(amogusSDF .< amogus_thres_val)

# ╔═╡ 58a66500-7f07-4cd7-b8a2-020a83177fa6
begin
	@time pettySDF = SerialSDF.dijkstraSDF2D(petty_threshold)
	SDFVis.toImageSDF(pettySDF)
end

# ╔═╡ 08e09dea-9f38-447b-9cd7-dc99a8dfc887
begin
	petty_slider = @bind petty_slider_val html"<input type=range value=0 min=0 max=50>"
	md"Threshold: $petty_slider "
end

# ╔═╡ fb412c5d-f2f4-4126-812c-fa41b2925ce8
begin
	petty_outline = copy(petty)
	for x in 1:size(petty,1)
		for y in 1:size(petty,2)
			if abs(pettySDF[x,y]) < petty_slider_val
				petty_outline[x,y] = RGB(0.1,0.1,0.1)
			end
		end
	end
petty_outline
end

# ╔═╡ c0c48ebd-afb8-4df7-868e-46efbddd5484
@time SerialSDF.dijkstraSDF2D(petty_threshold)

# ╔═╡ 68510fbf-c135-4157-9835-58f2d8950d85
@time ParallelSDF.dijkstraSDF2DSerialUDF(petty_threshold)

# ╔═╡ dc99c890-3401-488d-b551-2f75a6f472a2
@time ParallelSDF.dijkstraSDF2DParallelUDF(petty_threshold)

# ╔═╡ 53c81bbe-305e-4c69-8b85-2c3f56f5bbdc
@time SerialSDF.dijkstraSDF2D(amogus_threshold)

# ╔═╡ f575b998-3cf0-46e4-8340-43545ad58b5f
@time ParallelSDF.dijkstraSDF2DSerialUDF(amogus_threshold)

# ╔═╡ 655b3f4b-e900-4e85-a992-78421549b5cd
# ╠═╡ disabled = true
#=╠═╡
@time ParallelSDF.dijkstraSDF2DParallelUDF(amogus_threshold)
  ╠═╡ =#

# ╔═╡ Cell order:
# ╠═fbec5760-dfb8-11ed-19ce-8753ae969cde
# ╠═b8ff7564-a124-4605-8506-56b7ceeba8f8
# ╠═e02b0513-471d-43e4-8e3d-510fac112c4d
# ╠═5baa54db-95d8-4c71-a991-1b145128d788
# ╠═b3b9540e-fc92-43b6-b826-a74dbc6e325a
# ╠═6bf774a8-3990-4278-a011-34c7c450fc55
# ╠═6ffc2c09-39a4-4a0e-875d-8224e1ff119f
# ╠═e0f93c0e-28f3-4be5-815b-0180cb87c60a
# ╠═f64b95bd-1fbd-42cd-87c2-cd42dba63ca6
# ╠═84303126-2191-4bcc-8aa5-7b05ef62d57c
# ╠═54c6a8ed-03d3-4f71-b240-498782fb0bda
# ╠═a2a42489-278e-45f1-a872-9af4068adbd9
# ╠═32c138ab-f116-4127-824b-527520a206fc
# ╠═db3af62f-2e1a-4e54-ad9e-57bef579e609
# ╠═aa9bf6fe-a58e-4537-82f0-0eee4509c591
# ╠═cc3a2409-bdc7-4eb3-b0c3-a9dec01cf69d
# ╠═58a66500-7f07-4cd7-b8a2-020a83177fa6
# ╠═08e09dea-9f38-447b-9cd7-dc99a8dfc887
# ╠═fb412c5d-f2f4-4126-812c-fa41b2925ce8
# ╠═3ef1d920-6a8b-45cd-afc0-4ac1a450a952
# ╠═39f9509e-cca3-401e-b0c0-15dcea99d03e
# ╠═c0c48ebd-afb8-4df7-868e-46efbddd5484
# ╠═68510fbf-c135-4157-9835-58f2d8950d85
# ╠═dc99c890-3401-488d-b551-2f75a6f472a2
# ╠═53c81bbe-305e-4c69-8b85-2c3f56f5bbdc
# ╠═f575b998-3cf0-46e4-8340-43545ad58b5f
# ╠═655b3f4b-e900-4e85-a992-78421549b5cd
