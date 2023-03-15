using NumericalTools
using Documenter

DocMeta.setdocmeta!(NumericalTools, :DocTestSetup, :(using NumericalTools); recursive=true)

makedocs(;
    modules=[NumericalTools],
    authors="Chen Xia <xiachen1996@outlook.com> and contributors",
    repo="https://github.com/physcxia/NumericalTools.jl/blob/{commit}{path}#{line}",
    sitename="NumericalTools.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://physcxia.github.io/NumericalTools.jl",
        edit_link="dev",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/physcxia/NumericalTools.jl",
    devbranch="dev",
)
