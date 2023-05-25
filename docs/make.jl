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
            "manual/polarization_3D.md",
            "manual/polarization_2D.md",
            "manual/ladder_3D.md",
            "manual/polarization_approx.md",
            "manual/legendreinteraction.md",
            "manual/G_plus_Moroni.md",
            "manual/quasiparticle.md",
        ],
        "Reference" => Any[
            "lib/convention.md",
            "lib/parameter.md",
            "lib/selfenergy.md",
            "lib/polarization.md",
            "lib/interaction.md",
            "lib/legendreinteraction.md",
            "lib/twopoint.md",
            "lib/BSeq.md",
        ]
    ]
)

deploydocs(;
    repo="github.com/numericalEFT/ElectronGas.jl",
    devbranch="master"
)
