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


@testset "PrimaryDryFit: API for construction" begin
    t1 = collect(range(0.0u"hr", 10.0u"hr", length=5))
    T1a = collect(range(220.0u"K", 230.0u"K", length = length(t1)))
    T1b = collect(range(220.0u"K", 230.0u"K", length = length(t1)-1))
    T2 = collect(range(220.0u"K", 230.0u"K", length = length(t1)-2))
    t_end = 12.0u"hr"
    T1_iend = [length(T1a), length(T1b)]
    T2_iend = [length(T2)]
    # First, check that all constructors work as desired
    master1 = PrimaryDryFit(t1, (T1a, T1b), T1_iend, (T2,), T2_iend, t_end)
    @test PrimaryDryFit(t1, (T1a, T1b), T2, t_end) == master1
    master2 = PrimaryDryFit(t1, (T1a, T1b), T1_iend, (T2,), T2_iend, missing)
    @test PrimaryDryFit(t1, (T1a, T1b), T2,) == master2
    master3 = PrimaryDryFit(t1, (T1a, T1b), T1_iend, T2[end], missing, missing)
    @test PrimaryDryFit(t1, (T1a, T1b), T2[end]) == master3
    master4 = PrimaryDryFit(t1, (T1a, T1b), T1_iend, T2[end], missing, t_end)
    @test PrimaryDryFit(t1, (T1a, T1b), T2[end], t_end) == master4
    master5 = PrimaryDryFit(t1, (T1a, T1b), T1_iend, missing, missing, t_end)
    @test PrimaryDryFit(t1, (T1a, T1b), t_end) == master5
    master6 = PrimaryDryFit(t1, (T1a,), [length(T1a)], missing, missing, missing)
    @test PrimaryDryFit(t1, (T1a,)) == master6
    @test PrimaryDryFit(t1, T1a) == master6

    # Check that wrong constructions error
    @test_throws MethodError PrimaryDryFit(T1a, t_end)
    @test_throws MethodError PrimaryDryFit(t1, t_end)
    @test_throws MethodError PrimaryDryFit(t1, T1a, T2, T1b, t_end)
end


@testset "Simulation test against Python results: sucrose conventional" begin
    include("test_sucrose.jl")
end

@testset "Kv-Rp optimization" begin
    include("test_KRp_opt.jl")
end

ti = time() - ti
println("\nTest took total time of:")
println(round(ti/60, digits = 3), " minutes")
