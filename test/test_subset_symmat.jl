@testset "test subset_symmat.jl" begin

    ## Generate 5x5 symmetric input matrix
    A = [1 2 3 4 5;
        6 7 8 9 10;
        11 12 13 14 15;
        16 17 18 19 20]
    test_mat = A' * A

    ## Labels for rows/cols, and selection of which ones to retain
    ents = ["A", "B", "C", "D", "E"]
    ret = ["D", "B"]

    ## Subset matrix based on matching of strings
    sub_out = subset_symmat(test_mat, axes_entries = ents, retain = ret)
    @test sub_out.sub_entries == ["B", "D"]
    @test sub_out.M_sub == [486 562; 562 654]

    ## Should be same result leaving axes_entries as default nothing and passing
    ## a vector of integers for retain
    ## However, the returned sub_entries should be a vector of integers
    ret = [2, 4]
    sub_out = subset_symmat(test_mat, retain = ret)
    @test sub_out.sub_entries == [2, 4]
    @test sub_out.M_sub == [486 562; 562 654]

end
