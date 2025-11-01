using BuchbergerLearning
using Documenter

DocMeta.setdocmeta!(BuchbergerLearning, :DocTestSetup, :(using BuchbergerLearning); recursive=true)

makedocs(;
    modules=[BuchbergerLearning],
    authors="Nicolaus Jacobsen <jacobsen@rptu.de>",
    sitename="BuchbergerLearning.jl",
    format=Documenter.HTML(;
        canonical="https://Ktrompfl.github.io/BuchbergerLearning.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/Ktrompfl/BuchbergerLearning.jl",
    devbranch="main",
)
