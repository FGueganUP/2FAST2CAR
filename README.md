# 2FAST2CAR 
This is the latest reasonably stable development version of 2FAST2CAR, the updated version of FASTCAR (dev. Oscar Gayraud and Hugo Geindre). Note it may still be misbehaving, so be cautious :)

# Principle
FASTCAR is designed as a complementary tool to CREST (https://crest-lab.github.io/crest-docs/), enabling a fully automated exploration of the conformational degrees of freedom within reaction profiles at DFT level (with Gaussian) and with a reasonable computational effort.
Starting from a Gaussian output file (opt+freq) and a user-defined parameters list, FASTCAR proceeds accordingly:
- run a CREST calculation, producing a conformer ensemble within an adjustable energy window (DFTB);
- reject rotamers by the use of a internal tool from CREST (cregen);
- reject structures by the computation of a similarity index (invariant RMSD, threshold given by user);
- reoptimise all structures at a first DFT level using Gaussian (+freq if transition state);
- reject duplicate, unconverged structures ; in the case of transition states, also reject structures with imaginary normal mode to dissimilar to the one from the input structure;
- if required, recompute at a higher level of theory (not fully tested yet), and/or perform Intrinsic Reaction Profile calculations, computation of conceptual DFT descriptors or NBO analyses.

See our publication in PCCP: https://doi.org/10.1039/D4CP01721H

Python 3 program, relying on the following dependencies: SpyRMSD, CREST (including XTB), Gaussian, Slurm. 
Required Python libraries: glob, math, numpy, scipy, os, re, shutil, subprocess, sys, unicodedata, matplotlib, datetime

# 2FAST2CAR - additions
In the augmented version, the code was rewritten to ease inclusion of other functions & modification. Extension of the program's capabilities were also included (see the PDF manual for more details). 
