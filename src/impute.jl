
@doc doc"""
Perform simple imputation of missing SNP data

This function performs simple, non-positional imputation of missing data in a
SNP matrix. Note that imputation in this sense only refers to the estimation of
missing marker states from present non-missing data, not the imputation of
ungenotyped markers from a reference panel. In many cases a more advanced imputation
method taking biological information into account might be preferred, such as
Beagle, MaCH, etc.

##### Arguments

- snp_mat: A matrix of floats in minor-allele dosage format (i.e. 0 = homozygous
    major allele, 2 = homozygous minor allele, and 1 = heterozygous) with samples
    in rows, and SNPs in columns. Missing data within this matrix should be
    identified with NaN, rather than the missing type
- method: String of either 'mean', 'median', or 'major' specifying what value
    should be used to impute missing values within the SNP matrix

##### Returns

- snp_mat: The input SNP matrix, with missing data imputed
"""
function simple_impute(snp_mat::Matrix{AbstractFloat}, method::String)

    ## Define functions to calculate fill values
    if method == "mean"
        println("Imputing missing data for each SNP with mean of non-missing calls")
        fill_func = function(x); mean(snp_mat[:, i][.!isnan.(snp_mat[:, i])]); end
    else if method == "median"
        println("Imputing missing data for each SNP with median of non-missing calls")
        fill_func = function(x); median(snp_mat[:, i][.!isnan.(snp_mat[:, i])]); end
    else if method == "major"
        println("Imputing missing data for each SNP with major allele")
        fill_func = function(x); return(0.); end
    else
        error("Please select one of mean, median, or major for method")
    end

    ## Perform replacement of NaN values, one SNP (column) at a time
    for i in axes(snp_mat, 2)
        if NaN ∈ snp_mat[:, i]
            fill_val = fill_func(snp_mat[:, i])
            replace!(snp_mat[:, i], NaN => fill_val)
        else
            continue
        end
    end

    return(snp_mat)
end
