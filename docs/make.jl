CI = get(ENV, "CI", nothing) == "true" || get(ENV, "GITHUB_TOKEN", nothing) !== nothing
using LyoPronto
using Documenter
using Literate

@info "Using Literate to generate example"
Literate.markdown((@__DIR__)*"/example/fitting_mannitol.jl", "./src/generated", documenter=true)

@info "Building Documentation"
makedocs(;
    sitename = "LyoPronto.jl",
    modules = [LyoPronto],
    # This argument is only so that the sequence of pages in the sidebar is configured
    # By default all markdown files in `docs/src` are expanded and included.
    pages = [
        "Home" => "index.md",
        "Example, conventional lyo" => "generated/fitting_mannitol.md",
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
