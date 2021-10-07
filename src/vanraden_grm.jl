"""
    vanraden_grm(snp_mat; diag_weight=1e-4, test_posdef::Bool=false)

Calculate a genomic relationship matrix (GRM) from an input matrix of SNP data,
using the method of Van Raden, 2008 (https://doi.org/10.3168/jds.2007-0980).
Output consists of a named tuple with elements:

- GRM: The estimated GRM
- posdef: Boolean indicating whether the GRM is positive definite. `nothing` is
    returned if `test_posdef` argument is set to false

# Arguments

- `snp_mat::Matrix{Union{Float64, Int64}}`: A matrix of SNP data, with samples
    in rows and SNPs in columns. This matrix should be on the scale of 0, 1, 2
    where 0 indicates one homozygous state, 1 indicates the heterozygous state,
    and 2 indicates the opposite homozygous state. Additionally, the input matrix
    cannot contain any missing data. This matrix will often be in minor allele
    dosage format, such that 2 represents the minor allele, though
    this isn't strictly necessary. Note that the input format is only well-defined
    for biallelic variants.

# Keyword arguments

- `diag_weight=1e-4`: A small positive constant to add to diagonal elements
    of the GRM. Can help ensure that the GRM is positive definite in some cases.
    Set to 0 to disable this feature.
- `test_posdef::Bool=false`: Boolean indicating whether to test if the GRM is
    positive definite. Utilizes the LinearAlgebra module's isposdef() function,
    which attempts to perform a Cholesky decomposition. May add considerable runtime
    for a very large GRM
"""
function vanraden_grm(snp_mat; diag_weight=1e-4, test_posdef::Bool=false)

    if !isa(diag_weight, Number)
        error("diag_weight argument must be numeric - set to 0 to disable")
    end

    println("Calculating GRM from SNP data matrix")
    n_snps = size(snp_mat)[2]
    n_alleles = 2 * size(snp_mat)[1]

    p = zeros(Float16, n_snps)
    for i in 1:n_snps
        p[i] = sum(snp_mat[:, i]) / n_alleles
    end

    ## SNP matrix is recoded to -1, 0, 1
    snp_mat = snp_mat .- 1.
    for i in 1:n_snps
        snp_mat[:, i] = snp_mat[:, i] .- (2 * (p[i] - 0.5))
    end

    GRM = Symmetric(snp_mat * snp_mat' ./ (2 * sum(p .* (1 .- p))))

    ## Add in diagonal weighting if specified
    if diag_weight > 0
        GRM = GRM + diag_weight * I
    end

    ## Test that GRM is positive definite if specified
    if test_posdef
        return (GRM = GRM, posdef = isposdef(GRM))
    else
        return (GRM, posdef = nothing)
    end
end
