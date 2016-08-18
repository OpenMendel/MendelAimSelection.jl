# MendelAimSelection

This [Julia](http://julialang.org/) package selects the SNPs that are most informative at predicting ancestry for your data — the best Ancestry Informative Markers (AIMs). MendelAimSelection is one component of the umbrella [OpenMendel](https://openmendel.github.io) project.

[![](https://img.shields.io/badge/docs-current-blue.svg)](https://OpenMendel.github.io/MendelAimSelection.jl)

## Installation

*Note: Three OpenMendel packages - [SnpArrays](https://github.com/OpenMendel/SnpArrays.jl), [Search](https://github.com/OpenMendel/Search.jl), and [MendelBase](https://github.com/OpenMendel/MendelBase.jl) must be installed before any Mendel analysis packages will run.*

Within Julia, use the package manager to install MendelAimSelection:

    Pkg.clone("git@github.com:OpenMendel/MendelAimSelection.jl.git")

This package supports Julia v0.4.

## Data Files

To run this analysis package you will need to prepare a Control file and have your data files available. The Control file holds the names of your data files and any optional parameters for the analysis. Details on the general format and contents of the Control and data files can be found on the MendelBase [documentation page](https://openmendel.github.io/MendelBase.jl). Descriptions of the specific options available within the MendelAimSelection analysis package are in its [documentation page](https://openmendel.github.io/MendelAimSelection.jl).

There are example data files in the "docs" subfolder of each Mendel package, for example, ~/.julia/v0.4/MendelAimSelection/docs.

## Running the Analysis

To run this analysis package, first launch Julia. Then load the package with the command:

     julia> using MendelAimSelection

Next, if necessary, change to the directory containing your files, for example,

     julia> cd("~/path/to/data/files/")

Finally, to run the analysis using the parameters in the control file Control_file.txt use the command:

     julia> AimSelection("Control_file.txt")

*Note: The package is called* MendelAimSelection *but the analysis function is called simply* AimSelection.

## Citation

If you use this analysis package in your research, please cite the following reference in the resulting publications:

*Lange K, Papp JC, Sinsheimer JS, Sripracha R, Zhou H, Sobel EM (2013) Mendel: The Swiss army knife of genetic analysis programs. Bioinformatics 29:1568-1570.*

<!--- ## Contributing
We welcome contributions to this Open Source project. To contribute, follow this procedure ... --->

## Acknowledgments

This project is supported by the National Institutes of Health under NIGMS awards R01GM053275 and R25GM103774 and NHGRI award R01HG006139.
