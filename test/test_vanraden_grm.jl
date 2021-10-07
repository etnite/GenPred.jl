@testset "test vanraden_grm.jl" begin

    ## Input data of 3 samples x 6 SNPs
    ## and corresponding expected output GRM
    snp_data = [2 2 2 2 2 2; 1 1 1 1 1 1; 0 0 0 0 0 0]
    grm = [2. 0. -2.; 0. 0. 0.; -2. 0. 2.]

    ## Run without constant added to diagonal - results in GRM that is not
    ## positive definite
    vr_grm = vanraden_grm(snp_data, diag_weight = 0, test_posdef = true)
    @test vr_grm.GRM == grm
    @test !vr_grm.posdef

    ## Run again with 1e-5 added to diagonal
    vr_grm = vanraden_grm(snp_data, diag_weight = 1e-5, test_posdef = true)
    @test vr_grm.GRM == grm + 1e-5 * I
    @test vr_grm.posdef

    ## Run with float input matrix instead of Int
    vr_grm = vanraden_grm(convert(Matrix{Float64}, snp_data), diag_weight = 0, test_posdef = true)
    @test vr_grm.GRM == grm
    @test !vr_grm.posdef
end
