# Parallel Signed Distance Field Computation
Signed distance functions are functions which compute the Euclidean distance to the boundary points of a set and negates values for inputs which lie within the set. They have many applications, particularly in graphics where they may be used to manipulate surfaces and images. In these use cases, it is often useful to use a signed distance field which stores the result of the signed distance function at every point in space. 

Signed distance fields for 2D images can be computed such that each pixel stores the signed distance from its center to the nearest boundary point. This can also be extended to 3D voxel spaces and higher dimensions, but as the resolution and dimension count grow, the number of points to evaluate the signed distance function at increases quickly.

This project contains several algorithms for computing a signed distance field in 2-dimensional pixel spaces quickly using parallel programming in the Julia programming language. Although higher dimensions are outside of the scope of this project, many of these techniques can easily generalize, which can be useful for 3D graphics or n-dimensional optimization problems.

In order to optimize performance of the signed distance field computation, both parallel and serial considerations were made when comparing the three different algorithms. In particular, serial performance can be significantly optimized compared to brute-force algorithms before even considering parallel optimizations, but this may come at the cost of reduced ability to parallelize the algorithm. 

## Usage
To compute the signed distance field, an image must first be converted into a bit matrix through thresholding. This can be done with `SerialSDF.threshold()`. The signed distance field can then be computed using the appropriate functions in the `SerialSDF` or `ParallelSDF` modules.

(`SerialSDF.`/`ParallelSDF.`)`bruteSDF2D()` is a brute force implementation of the signed distance field computation. Although this is inefficient, it provides a good baseline for comparison.

(`SerialSDF.`/`ParallelSDF.`)`dijkstraSDF2D()` uses a modified version of Dijkstra's graph search, known also as brushfire algorithm, to compute the signed distance field via a graph search. This algorithm uses two passes of `SerialSDF.dijkstraUDF2D()` to first compute the unsigned distance field on the original and inverted versions of the input, then subtracts the inverted image's unsigned distance field to yield the signed distance field.

(`SerialSDF.`/`ParallelSDF.`)`linearSDF2D()` makes use of a linear-time algorithm to compute the signed distance field. Like brushfire algorithm, this algorithm uses two passes to compute the unsigned distance fields of the original and inverted versions of the inputs and subtracts them. Unlike the previous algorithm, computation of the unsigned distance field can be parallelized here, as shown with `ParallelSDF.linearUDF2D()`.

The [Presentation](/Presentation/) folder contains a Pluto notebook with some examples, and the [Benchmarking](/Benchmarking/) folder contains a Jupyter notebook for benchmarking the algorithms.