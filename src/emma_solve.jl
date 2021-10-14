"""
Perform a spectral decomposition of a genomic relationship matrix (K), using
either a Cholesky decomposition or an Eigen decomposition

Cholesky decomposition is used if the number of observations (n) is less than
the sum of the number of random effects (m) and fixed effects (p). Otherwise
an Eigen decomposition is performed
"""
function decompose(;n::Int64, p::Int64, m::Int64,
    Z, K::Union{Nothing, Symmetric{Float64}}, S::Matrix{Float64})

    n <= m + p ? spectral_meth = "eigen" : spectral_meth = "cholesky"

    ## Attempt Cholesky decomposition
    if spectral_meth == "cholesky"
        B = cholesky(K)
        ## R only returns the upper triangle for Cholesky decomp., which is then transposed
        ## We have the lower triangle handy so we can use it directly
        isnothing(K) ? ZBᵀ = Z : ZBᵀ = Z * B.L
        svd_ZBᵀ = svd(ZBᵀ)
        U = svd_ZBᵀ.U
        ϕ = vcat(svd_ZBᵀ.S .^ 2, zeros(n - m))
        SZBᵀ = S * ZBᵀ
        svd_SZBᵀ = try
            svd(SZBᵀ)
        catch
            svd(SZBᵀ + fill(1e-10, shape(SZBᵀ)))
        end

        QR_fact = qr(hcat(X, svd_SZBᵀ.U))
        Q = QR_fact.Q[:, (p + 1):n]
        R = QR_fact.R[p .+ 1:m, p .+ 1:m]
        θ = try
            ans = transpose(R .^ 2) / (svd_SZBᵀ.S .^ 2)
            vcat(ans, zeros(n - p - m))
        catch
            @warn "Cholesky decomposition failed - falling back to Eigen decomposition"
            missing
        end
    end

    ## Perform Eigen decomposition
    if spectral_meth == "eigen" || ismissing(θ)
        offset = sqrt(n)
        if isnothing(K)
            Hb = Z * Z' + offset * I
        else
            Hb = Z * K * Z' + offset * I
        end

        Hb_system = eigen(Symmetric(Hb), sortby = (x -> -x))  # By default Julia sorts eigvals/eigvecs in opposite order of R
        ϕ = Hb_system.values .- offset
        U = Hb_system.vectors

        SHbS_system = eigen(Symmetric(S * Hb * S), sortby = (x -> -x))
        θ = SHbS_system.values[1:(n - p)] .- offset
        Q = SHbS_system.vectors[:, 1:(n - p)]
    end

    return (U = U, Q = Q, ϕ = ϕ, θ = θ)
end


"""
Optimize Likelihood function

Finds minimum of likelihood function for either maximum likelihood (ML) or
restricted maximum likelihood (REML) methods.
"""
function opt_likelihood(;y, n, p, method, low_bound, hi_bound, Q, θ, ϕ)

    ω = Q' * y
    if method == "ML"
        df = n

        ## Here we define an anonymous function that is supplied to optimize(),
        ## which minimizes λ
        like_func = function(λ, df, θ, ω, ϕ)
            df * log(sum(ω.^2 ./ (θ .+ λ))) .+ sum(log.(ϕ .+ λ))
        end
        soln = optimize(x -> like_func(x, df, θ, ω, ϕ), low_bound, hi_bound, Brent())

    else  ## REML method
        df = n - p
        like_func = function(λ, df, θ, ω)
            df * log(sum(ω.^2 ./ (θ .+ λ))) .+ sum(log.(θ .+ λ))
        end
        soln = optimize(x -> like_func(x, df, θ, ω), low_bound, hi_bound, Brent())
    end

    ## Calculate the log-likelihood
    LL = -0.5 * (soln.minimum + df + df * log(2 * π / df))

    return (ω = ω, λ_opt = soln.minimizer, df = df, LL = LL)
end


"""
Generate outputs for emma_solve()
"""
function make_output(;ω, θ, ϕ, df, λ_opt, U, X, y, K, Z, return_SE, return_Hinv, LL)

    ## Calculate genotyping variance (Vu_opt) and error variance (Ve_opt)
    Vu_opt = sum(ω.^2 ./ (θ .+ λ_opt)) / df
    Ve_opt = λ_opt * Vu_opt

    H⁻¹ = U * (U' ./ (ϕ .+ λ_opt))
    W = X' * (H⁻¹ * X)

    ## Calculate fixed effects (BLUEs)
    β = vec((X' * (H⁻¹ * y)) / W)

    if isnothing(K)
        KZᵀ = Z'
    else
        KZᵀ = K * Z'
    end

    ## Calculate random effects (BLUPs)
    KZᵀH⁻¹ = KZᵀ * H⁻¹
    u =  vec(KZᵀH⁻¹ * (y .- (X * β)))   # in R: u <- array(KZt.H⁻¹ %*% (y - X %*% β))

    if return_SE
        W⁻¹ = inv(W)
        β_SE = sqrt.(Vu_opt .* diag(W⁻¹))
        WW = KZᵀH⁻¹ * KZᵀ'
        WWW = KZᵀH⁻¹ * X

        if isnothing(K)
            u_SE = vec(sqrt.(Vu_opt .* (ones(m) .- diag(WW) .+ diag(WWW * W⁻¹ * WWW'))))
        else
            u_SE = vec(sqrt.(Vu_opt .* (diag(K) .- diag(WW) .+ diag(WWW * W⁻¹ * WWW'))))
        end
    else
        u_SE = nothing
        β_SE = nothing
    end

    if !return_Hinv; H⁻¹ = nothing; end

    return (Vu = Vu_opt, Ve = Ve_opt, beta = β, beta_SE = β_SE, u = u, u_SE = u_SE, LL = LL, Hinv = H⁻¹)
end


"""
    emma_solve(y::Vector{Union{Float64, Missing}}, K::Matrix{Float64}=nothing)

Solve a system of mixed linear model equations using the Efficient Mixed Model
Association (EMMA) method.

This function is a port of the `mixed.solve()` function of the rrBLUP R package
(https://doi.org/10.3835/plantgenome2011.08.0024) which utilizes a modified version
of the EMMA algorithm described in Kang et al., 2008 (https://doi.org/10.1534/genetics.107.080101).
EMMA is a fast and efficient method for performing genomic best linear unbiased
prediction (GBLUP) but is limited by its ability to only incorporate one random
effect (genotypes) in addition to the error.

This function requires the use of design matrices to specify the relationships
between the phenotypic vector (y) and the fixed and random effects. In the case
of the random effects, it is usually easiest to simply ensure that the y vector
is sorted identically to the genomic relationship matrix (K) rows/columns and then
rely on the default identity matrix for the random effects design matrix (Z). However,
custom creation of a Z matrix is required if genotypes are replicated across
multiple fixed effect levels.

The code for this function is sparsely commented. Unfortunately the original
R mixed.solve() function code does not contain any comments, and while I understand
the gist of what is happening in each step, I don't understand the logic behind
everything.

# Keyword arguments

- `y::Vector{Union{Missing, Float64}}`: Vector (n × 1) of phenotypic observations.
    Missing values (NA) are omitted, along with the corresponding rows of X and Z.
- `wts::Union{Nothing, Vector{Float64}} = nothing`: An optional vector of weights
    for the phenotypic observations
- `X::Union{Nothing, Matrix{Float64}} = nothing`: Design matrix (n × p) for the
    fixed effects. If not passed, a vector of 1’s is used to model the intercept.
    X must be full column rank (implies β is estimable).
- `Z::Union{Nothing, Matrix{Float64}} = nothing`: Design matrix (n × m) for the
    random effects. If not passed, assumed to be the identity matrix.
- `K::Union{Nothing, Matrix{Float64}} = nothing`: Covariance matrix (m × m) for
    random effects; must be positive semi-definite. If not passed, assumed to be
    the identity matrix.
- `method::String="REML"`: Specifies whether the full ("ML") or restricted
    ("REML") maximum-likelihood method is used.
- `bounds::Vector{Float64} = [1e-09, 1e+09]`: Vector with two elements specifying
    the lower and upper bound for the ridge parameter
- `return_SE::Bool = false`: If true, standard errors are calculated
- `return_Hinv::Bool = false`: If true, the function returns the inverse of
    H = ZKZ' + λI. This is useful for GWAS.
"""
function emma_solve(;y,
    wts::Union{Nothing, Vector{Float64}} = nothing,
    X::Union{Nothing, Matrix{Float64}} = nothing,
    Z::Union{Nothing, Matrix{Float64}} = nothing,
    K::Union{Nothing, Matrix{Float64}, Symmetric{Float64}} = nothing,
    method::String = "REML",
    bounds::Vector{Float64} = [1e-09, 1e+09],
    return_SE::Bool = false,
    return_Hinv::Bool = false)

    n = length(y)
    notNA = findall(.!(ismissing.(y)))
    y = reshape(y, n, 1)

    ## Sanity checks on input arguments
    if !(method ∈ ("ML", "REML"))
        ArgumentError("method should be one of 'ML' or 'REML'")
    end
    if length(bounds) != 2
        ArgumentError("bounds should be a vector of floats with length 2")
    end

    ## Create X as column of ones if it wasn't provided
    if isnothing(X); X = ones(n, 1); end
    p = size(X)[2]

    ## Same thing for the weights vector
    if isnothing(wts); wts = ones(n); end

    ## Set Z to identity matrix if it wasn't provided
    if isnothing(Z)
        Z = Matrix(1.0I, n, n)  # Initialize Float identity matrix of size n x n
    end
    m = size(Z)[2]

    ## Verify compatible dimensions between y, X, Z, and wts
    if size(X)[1] != n
        error("number of rows of X matrix does not equal length of y vector")
    end
    if size(Z)[1] != n
        error("number of rows of Z matrix does not equal length of y vector")
    end
    if length(wts) != n
        error("length of weights vector does not equal length of y vector")
    end

    ## Check if supplied K is positive semidefinite and exit otherwise
    if !isnothing(K)
        K = Symmetric(K)
        diag(K) = diag(K) + 1e6
        if !isposdef(K)
            error("supplied K matrix is not positive semidefinite")
        end
    end

    ## Subset inputs to only include non-missing phenotypic observations
    X = X[notNA, :]
    n = length(notNA)
    y = y[notNA, :]
    Z = Z[notNA, :]
    wts = wts[notNA]

    ## Weight the inputs
    y = y .* wts
    X = X .* wts
    Z = Z .* wts

    ## Test that X is full rank and calculate S matrix
    if rank(X) < p
        error("X is not full rank")
    else
        S = I - (X * inv(X' * X)) * X'
    end

    ## Perform spectral decomposition
    decomp = decompose(n = n, p = p, m = m, Z = Z, K = K, S = S)

    ## Optimize the likelihood function
    optim = opt_likelihood(y = y, n = n, p = p, method = method,
        low_bound = bounds[1], hi_bound = bounds[2], Q = decomp.Q,
        θ = decomp.θ, ϕ = decomp.ϕ)

    ## Construct the output data structures
    answers = make_output(ω = optim.ω, θ = decomp.θ, ϕ = decomp.ϕ, df = optim.df,
        λ_opt = optim.λ_opt, U = decomp.U, X = X, y = y, K = K, Z = Z,
        return_SE = return_SE, return_Hinv = return_Hinv, LL = optim.LL)

    return answers

end
