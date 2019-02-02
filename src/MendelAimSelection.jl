__precompile__()

"""
This module orchestrates the selection of the most informative AIMs.
"""
module MendelAimSelection
#
# Required OpenMendel packages and modules.
#
using MendelBase
using SnpArrays
#
# Required external modules.
#
using CSV
using DataFrames
using Distributions
using Missings
using StatsBase

export AimSelection

"""
This is the wrapper function for the AIM Selection analysis option.
"""
function AimSelection(control_file = ""; args...)

  AIM_SELECTION_VERSION :: VersionNumber = v"0.5.0"
  #
  # Print the logo. Store the initial directory.
  #
  print(" \n \n")
  println("     Welcome to OpenMendel's")
  println("   AIM Selection analysis option")
  println("        version ", AIM_SELECTION_VERSION)
  print(" \n \n")
  println("Reading the data.\n")
  initial_directory = pwd()
  #
  # The user specifies the analysis to perform via a set of keywords.
  # Start the keywords at their default values.
  #
  keyword = set_keyword_defaults!(Dict{AbstractString, Any}())
  #
  # Keywords unique to this analysis should be first defined here
  # by setting their default values using the format:
  # keyword["some_keyword_name"] = default_value
  #
  #
  # Process the run-time user-specified keywords that will control the analysis.
  # This will also initialize the random number generator.
  #
  process_keywords!(keyword, control_file, args)
  #
  # Check that the correct analysis option was specified.
  #
  lc_analysis_option = lowercase(keyword["analysis_option"])
  if (lc_analysis_option != "" &&
      lc_analysis_option != "aimselection")
     throw(ArgumentError(
       "An incorrect analysis option was specified.\n \n"))
  end
  keyword["analysis_option"] = "AimSelection"
  #
  # Read the genetic data from the external files named in the keywords.
  #
  (pedigree, person, nuclear_family, locus, snpdata,
    locus_frame, phenotype_frame, pedigree_frame, snp_definition_frame) =
    read_external_data_files(keyword)
  #
  # Execute the specified analysis.
  #
  println(" \nAnalyzing the data.\n")
  execution_error = aim_selection_option(person, snpdata,
    pedigree_frame, snp_definition_frame, keyword)
  if execution_error
    println(" \n \nERROR: Mendel terminated prematurely!\n")
  else
    println(" \n \nMendel's analysis is finished.\n")
  end
  #
  # Finish up by closing, and thus flushing, any output files.
  # Return to the initial directory.
  #
  close(keyword["output_unit"])
  cd(initial_directory)
  return nothing

end # function AimSelection

"""
Ranks SNPs by their ancestry information content. All people should 
be assigned ancestry fractions and be fully typed. Ranks are assigned
by a likelihood ratio heterogeneity test. 
"""
function aim_selection_option(person::Person, snpdata::SnpDataStruct, 
  pedigree_frame::DataFrame, snp_definition_frame::DataFrame, 
  keyword::Dict{AbstractString, Any})
  #
  # Define scalar constants.
  #
  populations = person.populations
  people = person.people
  snps = snpdata.snps
  #
  # Allocate arrays and catalogue the ethnic groups.
  #
  ethnic = blanks(people)
  pedigree_field = names(pedigree_frame)
  ethnic_in_ped = (:Ethnic in pedigree_field)
  if ethnic_in_ped
    copyto!(ethnic, pedigree_frame[:Ethnic])
  else
    throw(ArgumentError("The Ethnic field is not included " *
      "in the pedigree data.\n \n"))
  end
  population = unique(ethnic)
  populations = length(population)
  alleles = zeros(populations)
  genes = zeros(populations)
  dosage = zeros(people)
  pvalue = ones(snps)
  #
  # Loop over the SNPs.
  #
  for snp = 1:snps
    if snpdata.maf[snp] <= 0.01; continue; end
    #
    # Copy the current SNP genotypes into a dosage vector.
    #
    copyto!(dosage, view(snpdata.snpmatrix, :, snp); impute = false)
    #
    # Tally reference alleles and genes in each population.
    #
    xlinked = uppercase(snpdata.chromosome[1]) == "X"
    fill!(alleles, 0.0)
    fill!(genes, 0.0)
    for i = 1:people
      if isnan(dosage[i]); continue; end
      if ethnic[i] == "" || ismissing(ethnic[i]); continue; end
      j = something(findfirst(isequal(ethnic[i]), population), 0)
      if xlinked && person.male[i]
        alleles[j] = alleles[j] + 0.5 * dosage[i]
        genes[j] = genes[j] + 1.0       
      else
        alleles[j] = alleles[j] + dosage[i]
        genes[j] = genes[j] + 2.0
      end
    end
    #
    # Add the maximum loglikelihoods for the different populations.
    #
    lrt = 0.0
    for j = 1:populations
      p = 0.0
      if genes[j] > 0.0
        p = alleles[j] / genes[j]
      end
      if 0.0 < p < 1.0
        n = round(Int, genes[j])
        x = round(Int, alleles[j])
        lrt = lrt + x * log(p) + (n - x)*log(1.0 - p)
      end
    end
    #
    # Subtract the maximum loglikelihood for the entire sample.
    # Based on this, compute the likelihood ratio p-value.
    #
    p = sum(alleles)/sum(genes)
    if 0.0 < p < 1.0
      n = round(Int, sum(genes))
      x = round(Int, sum(alleles))
      lrt = 2.0*(lrt - x * log(p) - (n - x)*log(1.0 - p))
      pvalue[snp] = ccdf(Chisq(1), lrt) 
    end
  end
  #
  # Rank SNPs by their p-values and deposit the rank of each SNP
  # in the SNP definition frame.
  #
  aim_rank = StatsBase.ordinalrank(pvalue)
  snp_definition_frame[:AIMRank] = aim_rank
  #
  # Display the SNP's AIM rank on the screen and in a file, if requested.
  #
  show(snp_definition_frame)
  output_file = string(keyword["output_file"])
  if output_file != ""
    CSV.write(output_file, snp_definition_frame;
      writeheader = true, delim = keyword["output_field_separator"],
      missingstring = keyword["output_missing_value"])
  end

  return execution_error = false
end # function aim_selection_option
#
# Method to obtain path to this package's data files
# so they can be used in the documentation and testing routines.
# For example, datadir("Control file.txt") will return
# "/path/to/package/data/Control file.txt"
#
datadir(parts...) = joinpath(@__DIR__, "..", "data", parts...)

end # module MendelAimSelection

