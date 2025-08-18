using Documenter, TraceParticles

makedocs(
    sitename="TraceParticles.jl",
    modules=[TraceParticles],
    checkdocs=:exports,
    authors="eilifso <eilif.oyre@gmail.com>",
    pages=[
        "Usage" => "index.md",
        "Index" => "index_list.md",
        "Examples" => "examples.md",
        "API" => "api.md",
    ]
)

deploydocs(
    repo="github.com/ITA-Solar/TraceParticles.jl.git",
)
