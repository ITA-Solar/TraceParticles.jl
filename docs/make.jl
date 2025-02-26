using Documenter, TestParticles

makedocs(
    sitename="TestParticles.jl",
    modules=[TestParticles],
    checkdocs = :exports,
    authors="eilifso <eilif.oyre@gmail.com>",
    pages=[
        "Usage" => "index.md",
        "Index" => "index_list.md",
        "Examples" => "examples.md",
        "API" => "api.md",
    ]
)

deploydocs(
    repo = "github.com/ITA-Solar/TestParticles.jl.git",
)