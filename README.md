# TensorDiagrams.jl

> **⚠️ This package is under development. There could be bugs and APIs may change without notice. We also don't provide decent documentation at this stage replacing it with a minimal example in the example folder.**

A Julia package for constructing, manipulating, and processing tensor network diagrams. In this package, a diagram is a standard tensor network diagram with optional label decorations on the legs, specifying which vector subspace each leg lives in. TensorDiagrams provides a representation of diagrams with support for geometric transformations, diagram composition, and code generation.

## Overview

TensorDiagrams.jl allows you to:

- **Define tensor network diagrams** with nodes (tensors) connected by labeled legs
- **Apply geometric transformations** including reflections and rotations
- **Compose diagrams** by gluing them along shared boundaries
- **Generate contraction code** compatible with `@tensoropt` from TensorOperations.jl
- **Visualize diagrams** using CairoMakie
- **Check consistency** of diagram structures
- **Complete partial diagrams** by filling in valid label combinations

## Core Structures

### TensorNode

A `TensorNode` represents a single tensor with:
- A name identifier
- Legs distributed across four sides: `left`, `right`, `top`, `bottom`
- Allowed labels for each leg (for diagram completion)
- TensorKit properties: dual spaces, flipped arrows, adjoint status

### TensorDiagram

A `TensorDiagram` represents a complete tensor network with:
- A collection of `TensorNode`s
- A contraction pattern (ncon-style indexing)
- Boundary legs connecting to external indices
- Labels indicating vector subspaces
- A scalar factor by which the contraction result should be multiplied


## Installation

```julia
using Pkg
Pkg.add(url="https://github.com/ebelnikola/TensorDiagrams")
```

## Configuration

The package uses YAML configuration files to define names of subspaces:

```yaml
INFINITE_LABELS:
  - x1
  - x2
FINITE_LABELS:
  - z1
  - z2
```

The package differentiates between two classes of subspaces:
- **Infinite subspaces** - rendered with red colored legs in visualizations
- **Finite subspaces** - rendered with green colored legs in visualizations

This distinction is useful for rigorous analysis of tensor networks. 

Load configuration with:
```julia
load_config("path/to/config.yml")
```