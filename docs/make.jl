using RayTracing
using Documenter

DocMeta.setdocmeta!(RayTracing, :DocTestSetup, :(using RayTracing); recursive=true)

makedocs(;
    modules=[RayTracing],
    authors="Abrudsky Martín Ezequiel <mabrudsky07@email.com>",
    sitename="RayTracing.jl",
    format=Documenter.HTML(;
        canonical="https://mabrudsky.github.io/RayTracing.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/mabrudsky/RayTracing.jl",
    devbranch="main",
)
