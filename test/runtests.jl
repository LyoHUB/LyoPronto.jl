using Test

# Run test suite
println("Starting tests")
ti = time()

@testset "LyoPronto.jl tests" begin
    @test 1 == 1
end

ti = time() - ti
println("\nTest took total time of:")
println(round(ti/60, digits = 3), " minutes")
