### Overview
Mendel AIM Selection is a component of the umbrella [OpenMendel](https://openmendel.github.io) project. This analysis option selects the SNPs that are most informative at predicting ancestry for your data — the best Ancestry Informative Markers (AIMs).

<!--- ### Appropriate Problems and Data Sets
 ... --->

### Installation
*Note: The three OpenMendel packages (1) [SnpArrays](https://openmendel.github.io/SnpArrays.jl/latest/), (2) [MendelSearch](https://openmendel.github.io/MendelSearch.jl), and (3) [MendelBase](https://openmendel.github.io/MendelBase.jl) must be installed before any other OpenMendel package will run. It is easiest if these three packages are installed in the above order and before any other OpenMendel package.*

Within Julia, use the package manager to install MendelAimSelection:

    (v1.1) pkg> add https://github.com/OpenMendel/MendelAimSelection.jl.git

This package supports Julia v1.0+

### Input Files
The MendelAimSelection analysis package uses the following input files. Example input files can be found in the [docs](https://github.com/OpenMendel/MendelAimSelection.jl/tree/master/docs) subfolder of the MendelAimSelection project. (An analysis won't always need every file type below.)

* [Control File](#control-file): Specifies the names of your data input and output files and any optional parameters (*keywords*) for the analysis. (For a list of common keywords, see [Keywords Table](https://openmendel.github.io/MendelBase.jl/#keywords-table)).
* [Locus File](https://openmendel.github.io/MendelBase.jl/#locus-file): Names and describes the genetic loci in your data.
* [Pedigree File](https://openmendel.github.io/MendelBase.jl/#pedigree-file): Gives information about your individuals, such as name, sex, family structure, and ancestry.
* [Phenotype File](https://openmendel.github.io/MendelBase.jl/#phenotype-file): Lists the available phenotypes.
* [SNP Definition File](https://openmendel.github.io/MendelBase.jl/#snp-definition-file): Defines your SNPs with information such as SNP name, chromosome, position, allele names, allele frequencies.
* [SNP Data File](https://openmendel.github.io/MendelBase.jl/#snp-data-file): Holds the genotypes for your data set. Must be a standard binary PLINK BED file in SNP major format. If you have a SNP data file you must have a SNP definition file.

<a id="control-file"></a>
### Control file
The Control file is a text file consisting of keywords and their assigned values. The format of the Control file is:

	Keyword = Keyword_Value(s)

Below is an example of a simple Control file to run AIM Selection:

	#
	# Input and Output files.
	#
	snpdata_file = 1000genomes_chr1_eas.bed
	snpdefinition_file = 1000genomes_chr1_eas.snp
	pedigree_file = 1000genomes_chr1_eas.ped
	output_file = 1000genomes_chr1_eas Output.txt
	#
	# Analysis parameters for AIM Selection option.
	#

In the example above, the four keywords specify the input and output files: *1000genomes_chr1_eas.bed*, *1000genomes_chr1_eas.snp*, *1000genomes_chr1_eas.ped* and *1000genomes_chr1_eas Output.txt*. The text after the '=' are the keyword values. The names of keywords are *not* case sensitive. (The keyword values *may* be case sensitive.) A list of OpenMendel keywords common to most analysis package can be found [here](https://openmendel.github.io/MendelBase.jl/#keywords-table).

### Data Files
AIM Selection requires a [Control file](https://openmendel.github.io/MendelBase.jl/#control-file), and a [Pedigree file](https://openmendel.github.io/MendelBase.jl/#pedigree-file). Genotype data can be included in the Pedigree file, in which case a [Locus file](https://openmendel.github.io/MendelBase.jl/#locus-file) is required. Alternatively, genotype data can be provided in a [SNP data file](https://openmendel.github.io/MendelBase.jl/#snp-data-file), in which case a [SNP Definition File](https://openmendel.github.io/MendelBase.jl/#snp-definition-file) is required. OpenMendel will also accept [PLINK format](http://zzz.bwh.harvard.edu/plink) FAM and BIM files. Details on the format and contents of the Control and data files can be found on the [MendelBase](https://openmendel.github.io/MendelBase.jl) documentation page. There are example data files in the AIM Selection [docs](https://github.com/OpenMendel/MendelAIMSelection.jl/tree/master/docs) folder.

### Running the Analysis

To run this analysis package, first launch Julia. Then load the package with the command:

     julia> using MendelAimSelection

Next, if necessary, change to the directory containing your files, for example,

     julia> cd("~/path/to/data/files/")

Finally, to run the analysis using the parameters in your Control file, for example, Control_file.txt, use the command:

     julia> AimSelection("Control_file.txt")

*Note: The package is called* MendelAimSelection *but the analysis function is called simply* AimSelection.

<!--- ### Interpreting the results
      ... --->

### Citation

If you use this analysis package in your research, please cite the following reference in the resulting publications:

*Lange K, Papp JC, Sinsheimer JS, Sripracha R, Zhou H, Sobel EM (2013) Mendel: The Swiss army knife of genetic analysis programs. Bioinformatics 29:1568-1570.*

<!--- ### Contributing
We welcome contributions to this Open Source project. To contribute, follow this procedure ... --->

### Acknowledgments

This project is supported by the National Institutes of Health under NIGMS awards R01GM053275 and R25GM103774 and NHGRI award R01HG006139.
