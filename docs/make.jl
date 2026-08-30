using Documenter, TraceParticles

# Makie needs a headless-friendly backend and a bitmap format that Documenter
# can embed directly in the generated HTML.
using CairoMakie
CairoMakie.activate!(type="png")

makedocs(
    sitename="TraceParticles.jl",
    modules=[TraceParticles],
    checkdocs=:exports,
    authors="eilifso <eilif.oyre@gmail.com>",
    format=Documenter.HTML(
        prettyurls=get(ENV, "CI", "false") == "true",
        edit_link="develop",
        size_threshold=1_000_000,
        size_threshold_warn=500_000,
    ),
    pages=[
        "Home" => "index.md",
        "Concepts" => "concepts.md",
        "Examples" => "examples.md",
        "Verification gallery" => "verification.md",
        "API" => "api.md",
        "Index" => "index_list.md",
    ]
)

deploydocs(
    repo="github.com/ITA-Solar/TraceParticles.jl.git",
    devbranch="develop",
)
