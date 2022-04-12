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
        assets=String[]
    ),
    pages=[
        "Home" => "index.md",
        "Manual" => [
            "manual/fock.md",
            "manual/polarization.md",
            "manual/polarization_2D.md",
            "manual/legendreinteraction.md",
            "manual/G_plus_Moroni.md",
            "manual/quasiparticle.md",
        ],
    ]
)

deploydocs(;
    repo="github.com/numericalEFT/ElectronGas.jl",
    devbranch="master"
)
