using LyoPronto
using Test

# Run test suite
println("Starting tests")
ti = time()

@testset "RpFormFit" begin
    R0, A1, A2 = 1.0u"cm^2*hr*Torr/g", 14u"cm*hr*Torr/g", 1u"1/cm"
    Rp = RpFormFit(R0, A1, A2)
    @test_throws Unitful.DimensionError Rp(0)
    @test_throws Unitful.DimensionError Rp(5)
    @test Rp(5.0u"cm") == R0 + A1*5.0u"cm"/(1 + A2*5.0u"cm")
    @test Rp(0.0u"cm") == R0 
end

@testset "RampedVariable" begin
    pch = RampedVariable(150u"mTorr")
    @test pch(-Inf) == 150u"mTorr"
    @test pch(-Inf) == pch(0)
    @test pch(-Inf) == pch(15u"hr")
    Tsh = RampedVariable([228.15, 248.15]u"K", 1u"K/minute")
    @test_throws Unitful.DimensionError Tsh(0) 
    @test Tsh(0u"s") == 228.15u"K"
    @test Tsh(Inf*u"s") == 248.15u"K"
    @test Tsh(10u"minute") == 238.15u"K"
    @test Tsh(20u"minute") == 248.15u"K"
    P_per_vial = @test_logs (:warn, "Ramp rate given with probably the wrong sign, changing its sign") RampedVariable([40u"W", 20u"W", 10u"W"], [1u"W/minute", Inf*u"W/minute"], [1u"hr"])
    @test_throws Unitful.DimensionError P_per_vial(0)
    @test P_per_vial(0u"s") == 40u"W"
    @test P_per_vial(10u"minute") == 30u"W"
    @test P_per_vial(20u"minute") == 20u"W"
    @test P_per_vial(50u"minute") == 20u"W"
    @test P_per_vial(79u"minute") == 20u"W"
    @test P_per_vial(81u"minute") == 10u"W"
    @test P_per_vial(Inf*u"minute") == 10u"W"

end


@testset "Simulation test: sucrose conventional" begin
    include("test_sucrose.jl")
end

ti = time() - ti
println("\nTest took total time of:")
println(round(ti/60, digits = 3), " minutes")
