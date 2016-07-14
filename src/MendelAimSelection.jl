"""
This module orchestrates the selection of the most informative AIMs.
"""
module MendelAimSelection
#
# Other OpenMendel modules.
#
using MendelBase
using SnpArrays
#
# External modules.
#
using DataFrames    # From package DataFrames.
using Distributions # From package Distributions.

export AimSelection

"""
This is the wrapper function for the AIM Selection analysis option.
"""
function AimSelection(control_file = ""; args...)

  const AIM_SELECTION_VERSION :: VersionNumber = v"0.1.0"
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
  keyword = set_keyword_defaults!(Dict{ASCIIString, Any}())
  #
  # Keywords unique to this analysis may be defined here
  # by setting their default values using the format:
  # keyword["some_keyword_name"] = value
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
  # Execute the specifed analysis.
  #
  println(" \nAnalyzing the data.\n")
  execution_error = aim_selection_option(pedigree, person, nuclear_family,
    locus, snpdata, locus_frame, phenotype_frame, pedigree_frame,
    snp_definition_frame, keyword)
  if execution_error
    println(" \n \nERROR: Mendel terminated prematurely!\n")
  else
    println(" \n \nMendel's analysis is finished.\n")
  end
  #
  # Finish up by closing and thus flushing, any output files.
  # Return to the initial directory.
  #
  close(keyword["output_unit"])
  cd(initial_directory)
  return nothing

end # function AimSelection

"""
This function uses a likelihood ratio test to rank SNPs
by their ancestry information content. All people should be
assigned ancestry fractions and fully SNP-typed.
"""
function aim_selection_option(pedigree::Pedigree, person::Person,
  nuclear_family::NuclearFamily, locus::Locus, snpdata::SnpData, 
  locus_frame::DataFrame, phenotype_frame::DataFrame, 
  pedigree_frame::DataFrame, snp_definition_frame::DataFrame,
  keyword::Dict{ASCIIString, Any})
  #
  # Define scalar constants.
  #
  populations = person.populations
  people = person.people
  snps = snpdata.snps
  #
  # Allocate arrays.
  #
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
    copy!(dosage, slice(snpdata.snpmatrix, :, snp); impute = false)
    #
    # Tally reference alleles and genes in each population.
    #
    xlinked = uppercase(snpdata.chromosome[1]) == "X"
    fill!(alleles, 0.0)
    fill!(genes, 0.0)
    for i = 1:people
      if isnan(dosage[i]); continue; end
      if xlinked && person.male[i]
        for j = 1:populations
          alleles[j] = alleles[j] + 0.5 * person.admixture[i, j] * dosage[i]
          genes[j] = genes[j] + person.admixture[i, j]
        end
      else
        for j = 1:populations
          alleles[j] = alleles[j] + person.admixture[i, j] * dosage[i]
          genes[j] = genes[j] + 2.0 * person.admixture[i, j]
        end
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
      n = round(Int, genes[j])
      x = round(Int, alleles[j])
      lrt = lrt + logpdf(Binomial(n, p), x)
    end
    #
    # Subtract the maximum loglikelihood for the entire sample.
    #
    p = sum(alleles)/sum(genes)
    n = round(Int, sum(genes))
    x = round(Int, sum(alleles))
    #
    # Compute the likelihood ratio pvalue.
    #
    lrt = lrt - logpdf(Binomial(n, p), x)
    pvalue[snp] = ccdf(Chisq(1), lrt)
  end
  #
  # Rank snps by their pvalues and deposit the rank of each snp in the
  # snp definition frame.
  #
  aim_rank = ordinalrank(pvalue)
  snp_definition_frame[:AIMRank] = aim_rank
  show(snp_definition_frame)
  return execution_error = false
end # function aim_selection_option

end # module MendelAimSelection

