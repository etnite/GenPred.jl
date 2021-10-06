
"""
    simple_imp(snp_mat; method = "major")

Return a copy of an input matrix of minor-allele dosage SNP data where all
missing data (NaN) for each SNP has been imputed with either:

- The mean value of the SNP
- The median value of the SNP
- The major allele

This function performs simple, non-positional imputation of missing data in a
SNP matrix. Note that imputation in this sense only refers to the estimation of
missing marker states from present non-missing data, not the imputation of
ungenotyped markers from a reference panel. In many cases a more advanced imputation
method taking biological information into account might be preferred, such as
Beagle, MaCH, etc.

# Arguments

- `snp_mat::Matrix{Float64}`: A matrix of SNP data in minor-allele dosage format (i.e. 0 = homozygous
    major allele, 2 = homozygous minor allele, and 1 = heterozygous) with samples
    in rows, and SNPs in columns. Note that minor-allele dosage format is only
    well-definedfor biallelic variants. Missing data within this matrix should be
    identified with `NaN`, rather than `missing`. Since NaN is part of the
    Float64 type, this matrix does NOT require a Union{Float64, Missing} type

# Keyword arguments

- `method::String=major`: String of either 'mean', 'median', or 'major' specifying what value
    should be used to impute missing values within the SNP matrix
"""
function simple_imp(snp_mat::Matrix{Float64}; method::String = "major")

    out_mat = deepcopy(snp_mat)

    ## You can't define a function inside an if/else block because it will be
    ## limited to that block's scope. What we can do is create a function which
    ## returns another function based on its input (method)
    function ff_sel(meth)
        if meth == "mean"
            ff = x -> mean(x[.!isnan.(x)])
        elseif meth == "median"
            ff = x -> median(x[.!isnan.(x)])
        elseif meth == "major"
            ff = x -> 0.
        else
            error("Please select one of mean, median, or major for method")
        end
    end
    fill_func = ff_sel(method)

        ## Perform replacement of NaN values, one SNP (column) at a time
    for i in axes(snp_mat, 2)
        fill_val = fill_func(out_mat[:, i])
        out_mat[:, i] = replace(out_mat[:, i], NaN => fill_val)
    end

    return(out_mat)
end


"""
    simple_imp!(snp_mat; method = "major")

Identical to simple_imp(), except modifies the input SNP matrix rather than
returning a copy.
"""
function simple_imp!(snp_mat::Matrix{Float64}; method::String = "major")

    ## You can't define a function inside an if/else block because it will be
    ## limited to that block's scope. What we can do is create a function which
    ## returns another function based on its input (method)
    function ff_sel(meth)
        if meth == "mean"
            ff = x -> mean(x[.!isnan.(x)])
        elseif meth == "median"
            ff = x -> median(x[.!isnan.(x)])
        elseif meth == "major"
            ff = x -> 0.
        else
            error("Please select one of mean, median, or major for method")
        end
    end
    fill_func = ff_sel(method)

        ## Perform replacement of NaN values, one SNP (column) at a time
    for i in axes(snp_mat, 2)
        fill_val = fill_func(snp_mat[:, i])
        snp_mat[:, i] = replace(snp_mat[:, i], NaN => fill_val)
    end

    return(snp_mat)
end
