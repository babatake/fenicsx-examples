🎯 Overview of the Problem Setup
📌 Target Structure:

    A rectangular waveguide cross-section with width Lx=1.0Lx​=1.0 and height Ly=0.5Ly​=0.5

    The waveguide is modeled as a 2D problem

    All boundaries are perfect electric conductors (PEC) → enforced boundary condition: Hz=0Hz​=0

📌 Objective:

    Eigenvalues: Transverse wave numbers squared kt2kt2​ for each mode

    Eigenvectors: Spatial distributions of the z-component of the magnetic field Hz(x,y)Hz​(x,y)

🔧 Mathematical Formulation
Derived from Source-Free Maxwell's Equations:

Assuming a time-harmonic field dependence of the form ejωt−jkzzejωt−jkz​z, the TE mode satisfies the following Helmholtz equation for Hz(x,y)Hz​(x,y):
∇t2Hz+kt2Hz=0in Ω
∇t2​Hz​+kt2​Hz​=0in Ω

    ∇t2∇t2​: Transverse Laplacian in the (x,y)(x,y)-plane

    kt2=ω2με−kz2kt2​=ω2με−kz2​: Transverse wave number squared

Boundary Condition:

    PEC boundaries → Dirichlet condition:
    Hz=0on ∂Ω
    Hz​=0on ∂Ω

🔁 Weak Formulation via Finite Element Method

The weak (variational) form of the equation becomes:
∫Ω∇Hz⋅∇v dΩ=kt2∫ΩHz⋅v dΩ∀v
∫Ω​∇Hz​⋅∇vdΩ=kt2​∫Ω​Hz​⋅vdΩ∀v

This is formulated as a generalized eigenvalue problem:
[A]{Hz}=λ[B]{Hz},λ=kt2
[A]{Hz​}=λ[B]{Hz​},λ=kt2​

Where:

    [A][A]: Stiffness matrix

    [B][B]: Mass matrix

    λλ: Eigenvalue (related to propagation constant)

📌 What the Script Does
Step	Description
Mesh generation	Discretizes the domain [0,Lx]×[0,Ly][0,Lx​]×[0,Ly​] into triangular elements
Function space definition	Hz∈V⊂H01(Ω)Hz​∈V⊂H01​(Ω) using second-order continuous elements
Matrix assembly	Defines bilinear forms: a(u,v)=∇u⋅∇va(u,v)=∇u⋅∇v, b(u,v)=u⋅vb(u,v)=u⋅v
Apply PEC boundaries	Imposes Hz=0Hz​=0 on all boundaries (Dirichlet condition)
Solve eigenvalue problem	Uses SLEPc to compute eigenvalues λi=kt,i2λi​=kt,i2​ and eigenfunctions Hz,iHz,i​
Visualization	Displays the fundamental mode (lowest eigenvalue) using PyVista
🎯 Key Results

    kt2kt2​: Transverse wave number squared for each mode

    Hz(x,y)Hz​(x,y): Spatial distribution of the magnetic field (z-component)

→ These results form the basis for analyzing dispersion characteristics vs. frequency or propagation constant kzkz​
🔍 Application Examples

    Mode analysis in microwave waveguides

    Determination of cutoff frequencies in shielded or optical waveguides

    Resonant mode analysis in RF structures or cavities
