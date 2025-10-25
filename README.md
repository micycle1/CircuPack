[![](https://jitpack.io/v/micycle1/CircuPack.svg)](https://jitpack.io/#micycle1/CircuPack)
# CircuPack

**Define tangencies. Get precise circle packings.**

A robust Java library that computes Euclidean circle packings from triangulations. Inspired by Gerald Orick’s [GOPack](https://github.com/kensmath/GOPack).

<p align="center">
    <img src="resources/poisson_pack.png" width="31%" height="auto" >
    <img src="resources/hex_pack.png" width="31%" height="auto">
    <img src="resources/10k.png" width="31%" height="auto">
</p>

## Core Idea

Supply a `Triangulation` that defines *which circles are tangent to which*. CircuPack handles all the geometry: it calculates exact radii and center positions so every specified tangency is visually perfect.


## Key Features

*   **Simple Input:** Implement the `Triangulation` interface to define your graph of tangencies.
*   **Automatic Geometry:** No need to guess positions or sizes. The engine computes everything.
*   **Multiple Layouts:** Pack into a disc (MAX_PACK), a custom polygon, or a rectangle [WIP🚧].
*   **Production-Ready:** Stable, reliable and fast (10k circles in ~1s).


## Quick Start

```java
// 1. Define your triangulation (e.g., from Tinfour)
IIncrementalTin tin = new IncrementalTin();
// ... add vertices ...
Triangulation tri = new TinfourTriangulation(tin);

// 2. Initialize the packer
CircuPacker packer = new CircuPacker(tri);
packer.initialize();

// 3. Compute the packing
double maxError = 0.01; // max visual error (overlap, as a fraction of circle radius)
packer.riffle(maxError); // Iterate until max error is below threshold

// 4. Retrieve results
double[] radii = packer.getRadii();
double[] centersX = packer.getCentersX();
double[] centersY = packer.getCentersY();
```
