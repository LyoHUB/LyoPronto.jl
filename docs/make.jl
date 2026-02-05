CI = get(ENV, "CI", nothing) == "true" || get(ENV, "GITHUB_TOKEN", nothing) !== nothing
using LyoPronto
using Documenter
using Literate

@info "Using Literate to generate examples"
for file in ["fitting_mannitol.jl", "fitting_rf_mannitol.jl", "all_recipes.jl", "utilities.jl"]
    Literate.markdown((@__DIR__)*"/example/$file", (@__DIR__)*"/src/generated", documenter=true)
end

@info "Building Documentation"
makedocs(;
    sitename = "LyoPronto.jl",
    modules = [LyoPronto],
    # This argument is only so that the sequence of pages in the sidebar is configured
    # By default all markdown files in `docs/src` are expanded and included.
    pages = [
        "Home" => "index.md",
        "Example, conventional lyo" => "generated/fitting_mannitol.md",
        "Example, microwave-assisted lyo" => "generated/fitting_rf_mannitol.md",
        "Other tools" => "generated/utilities.md",
        "Plot recipes" => "generated/all_recipes.md",
        "Reference" => "alldocstrings.md",
    ],
    # Don't worry about what `CI` does in this line.
    format = Documenter.HTML(prettyurls = CI),
)

@info "Deploying Documentation"
if CI
    deploydocs(
        # `repo` MUST be set correctly. Once your GitHub name is set
        # the auto-generated documentation will be hosted at:
        # https://PutYourGitHubNameHere.github.io/LyoPronto.jl/dev/
        # (assuming you have enabled `gh-pages` deployment)
        repo = "github.com/LyoHUB/LyoPronto.jl.git",
        target = "build",
        push_preview = true,
        devbranch = "main",
    )
end

@info "Finished with Documentation"
