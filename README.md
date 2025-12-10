# Tomographic reconstruction of femtosecond laser written waveguides [ARCHIVED]

Research code from my postdoc at Medical University of Innsbruck for waveguide tomography using proximal gradient methods (TV-norm regularization).

Published in: 
Nicolas Barré, Ravi Shivaraman, Lisa Ackermann, Simon Moser, Michael Schmidt, Patrick Salter, Martin Booth, and Alexander Jesacher, *"Tomographic refractive index profiling of direct laser written waveguides,"* **[Opt. Express 29](https://doi.org/10.1364/OE.434846)**, 35414-35425 (2021)

**Status**: This proof-of-concept has been reimplemented with improvements in [FluxOptics.jl Tutorial 4](https://anscoil.github.io/FluxOptics.jl/stable/tutorials/04_waveguide_tomography/), which demonstrates that simpler shrinkage (ISTA) performs better than TV-norm for this problem, avoiding staircase artifacts while maintaining reconstruction quality.

For modern waveguide inverse design, see [FluxOptics.jl](https://github.com/anscoil/FluxOptics.jl).
