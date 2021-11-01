using ElectronGas
using Documenter

DocMeta.setdocmeta!(ElectronGas, :DocTestSetup, :(using ElectronGas); recursive=true)

makedocs(;
    modules=[ElectronGas],
    authors="Kun Chen",
    repo="https://github.com/numericalEFT/ElectronGas.jl/blob/{commit}{path}#{line}",
    sitename="ElectronGas.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://numericalEFT.github.io/ElectronGas.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/numericalEFT/ElectronGas.jl",
    devbranch="master",
)
