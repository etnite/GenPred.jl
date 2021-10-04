using GenPred
using Documenter

DocMeta.setdocmeta!(GenPred, :DocTestSetup, :(using GenPred); recursive=true)

makedocs(;
    modules=[GenPred],
    authors="Brian Ward <brian@brianpward.net>",
    repo="https://github.com/etnite/GenPred.jl/blob/{commit}{path}#{line}",
    sitename="GenPred.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://etnite.github.io/GenPred.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/etnite/GenPred.jl",
)
