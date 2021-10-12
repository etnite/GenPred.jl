"""
    subset_symmat(M, axes_entries=nothing, retain)

Subset a symmetric matrix using a set of row/column entries to retain. Output
consists of a named tuple with elements:

- `M_sub`: The subsetted matrix
- 'sub_entries': The labels of the rows/columns of the subsetted matrix

# Arguments

- `M`: A symmetric matrix to subset rows/columns out of. Note that the input
    matrix is merely assumed to be symmetric, and does not have to pass
    issymmetric() from the LinearAlgebra module
- `axes_entries=nothing`: A vector labelling the rows and columns of the input
    matrix. If nothing, then the function assumes that the matrix has no labels for
    row/column entries. In this case a vector of integers can be passed for
    `retain` in order to subset
- `retain`: Either a vector of integers specifying which rows/columns of the input
    matrix to retain, or else a vector of strings specifying which entries from
    `
"""
function subset_symmat(M; axes_entries = nothing, retain)

    if size(M, 1) != size(M, 2)
        error("M is not square")
    end

    if isnothing(axes_entries)
        axes_entries = 1:size(M, 1)
    else
        if length(axes_entries) != size(M, 1)
            error("Length of axes_entries must match number of rows and columns of M")
        end
    end

    ind = sort(findall(axes_entries .∈ Ref(retain)))

    return (M_sub = deepcopy(M)[ind, ind], sub_entries = axes_entries[ind])
end
