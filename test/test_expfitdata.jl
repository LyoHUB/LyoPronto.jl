using LyoPronto
using Test


t1 = collect(range(0.0u"hr", 10.0u"hr", length=5))
T1a = collect(range(220.0u"K", 230.0u"K", length = length(t1)))
T1b = collect(range(220.0u"K", 230.0u"K", length = length(t1)-1))
T2 = collect(range(220.0u"K", 230.0u"K", length = length(t1)-2))
t_end = 12.0u"hr"

@testset "Data object constructors" begin
    # TfData: default t_range is 1:length(Tf), t defaults to missing
    tf1 = TfData(T1a)
    @test tf1.Tf == T1a
    @test tf1.t_range == 1:length(T1a)
    @test ismissing(tf1.t)
    # TfData: explicit t_range (sub-window of the time vector)
    tf2 = TfData(T1b, 2:length(T1b)+1)
    @test tf2.t_range == 2:length(T1b)+1
    # TfData: with its own time vector
    tf3 = TfData(T1a; t=t1)
    @test tf3.t == t1

    # TvwSeriesData: default t_range
    tvw1 = TvwSeriesData(T2)
    @test tvw1.Tvw == T2
    @test tvw1.t_range == 1:length(T2)
    @test ismissing(tvw1.t)
    # TvwSeriesData: explicit t_range
    tvw2 = TvwSeriesData(T2, 1:length(T2))
    @test tvw2.t_range == 1:length(T2)

    # TvwEndData: a single endpoint temperature
    tvwe = TvwEndData(T2[end])
    @test tvwe.Tvw_end == T2[end]

    # EndTimeData: a single time
    te1 = EndTimeData(t_end)
    @test te1.t_end == t_end
    # EndTimeData: a tuple of two times (a window)
    te2 = EndTimeData((10.0u"hr", 12.0u"hr"))
    @test te2.t_end == (10.0u"hr", 12.0u"hr")
end

@testset "Data object validation" begin
    # TfData: Tf must be a vector of temperatures
    @test_throws ArgumentError TfData(T1a[1], 1:1)
    @test_throws ArgumentError TfData(collect(range(0.0, 1.0, length=5)), 1:5)
    # TfData: t_range must be a range of integers
    @test_throws ArgumentError TfData(T1a, [1, 2, 3, 4, 5])
    # TfData: t_range length must match Tf length
    @test_throws ArgumentError TfData(T1a, 1:4)
    # TfData: t must be missing or a vector of times
    @test_throws ArgumentError TfData(T1a; t=5.0u"hr")
    @test_throws ArgumentError TfData(T1a; t=collect(1.0:5.0))

    # TvwSeriesData: Tvws must be a vector of temperatures
    @test_throws ArgumentError TvwSeriesData(T2[1], 1:1)
    @test_throws ArgumentError TvwSeriesData(collect(range(0.0, 1.0, length=3)), 1:3)
    # TvwSeriesData: t_range length must match Tvws length
    @test_throws ArgumentError TvwSeriesData(T2, 1:2)

    # TvwEndData: must be a single temperature
    @test_throws ArgumentError TvwEndData(220.0)
    @test_throws ArgumentError TvwEndData(collect(range(220.0u"K", 230.0u"K", length=3)))

    # EndTimeData: must be a time or a tuple of two times
    @test_throws ArgumentError EndTimeData(12.0)
    @test_throws ArgumentError EndTimeData((12.0, 13.0))
end

@testset "time traits" begin
    # Simplest cases
    td = TfData(T1a)
    tvd = TvwSeriesData(T1a)
    @test LyoPronto.time_bound_data(td)
    @test !LyoPronto.has_timevec(td)
    @test !LyoPronto.nontrivial_t_range(td)
    @test LyoPronto.time_bound_data(tvd)
    @test !LyoPronto.has_timevec(tvd)
    @test !LyoPronto.nontrivial_t_range(tvd)
    # No time vectors anyway
    @test !LyoPronto.time_bound_data(TvwEndData(T2[end]))
    @test !LyoPronto.time_bound_data(EndTimeData(t_end))
    # Own time vector
end

@testset "Container constructors" begin
    # Tuple form: two Tf series + one Tvw series + t_end
    master1 = ExpFitData(t1, (TfData(T1a), TfData(T1b), TvwSeriesData(T2), EndTimeData(t_end)))
    @test ExpFitData(t1, (TfData(T1a), TfData(T1b), TvwSeriesData(T2), EndTimeData(t_end))) == master1
    # Varargs form: two Tf series + one Tvw series, no t_end
    master2 = ExpFitData(t1, TfData(T1a), TfData(T1b), TvwSeriesData(T2))
    @test ExpFitData(t1, TfData(T1a), TfData(T1b), TvwSeriesData(T2)) == master2
    # Single object
    master3 = ExpFitData(t1, TfData(T1a))
    @test ExpFitData(t1, (TfData(T1a),)) == master3

    # The data tuple holds the expected object types, in the given order
    @test master1.data isa Tuple
    @test all([o isa T for (o, T) in zip(master1.data, [TfData, TfData, TvwSeriesData, EndTimeData])])
    @test all([o isa T for (o, T) in zip(master2.data, [TfData, TfData, TvwSeriesData])])
    @test all([o isa T for (o, T) in zip(master3.data, [TfData])])

    # Container validation
    @test_throws ArgumentError ExpFitData(collect(1.0:5.0), (TfData(T1a),))
    @test_throws ArgumentError ExpFitData(t1, ())
    @test_throws ArgumentError ExpFitData(t1, (T1a,))
end

@testset "fit_t" begin
    # When obj.t is missing, use the container's t
    tf = TfData(T1a)
    efd = ExpFitData(t1, (tf,))
    @test LyoPronto.fit_t(efd, tf) == t1
    # When obj.t is set, use the object's own t
    t_alt = collect(range(0.0u"hr", 4.0u"hr", length=length(T1a)))
    tf_t = TfData(T1a; t=t_alt)
    efd2 = ExpFitData(t1, (tf_t,))
    @test LyoPronto.fit_t(efd2, tf_t) == t_alt
    # For objects without a time vector, always use the container's t
    te = EndTimeData(t_end)
    efd3 = ExpFitData(t1, (TfData(T1a), te))
    @test LyoPronto.fit_t(efd3, te) == t1
end

@testset "Equality" begin
    @test TfData(T1a) == TfData(T1a)
    @test !(TfData(T1a) == TfData(T1b))
    @test TfData(T1a, 1:5) == TfData(T1a)
    @test !(TfData(T1a, 1:5) == TfData(T1a, 2:6))
    @test TvwSeriesData(T2) == TvwSeriesData(T2)
    @test TvwEndData(T2[end]) == TvwEndData(T2[end])
    @test EndTimeData(t_end) == EndTimeData(t_end)
    @test ExpFitData(t1, (TfData(T1a),)) == ExpFitData(t1, (TfData(T1a),))
    @test !(ExpFitData(t1, (TfData(T1a),)) == ExpFitData(t1, (TfData(T1b),)))
end

@testset "show" begin
    @test string(TfData(T1a)) == "TfData(5 pts, t_range=1:5)"
    @test string(TvwSeriesData(T2)) == "TvwSeriesData(3 pts, t_range=1:3)"
    @test string(TvwEndData(T2[end])) == "TvwEndData(230.0 K)"
    @test string(EndTimeData(t_end)) == "EndTimeData(12.0 hr)"
    @test string(ExpFitData(t1, (TfData(T1a), TfData(T1b)))) == "ExpFitData(5 t pts; TfData, TfData)"
end

# ============================================================================
# Model setup: shared parameters for the Pikal (2-state) and RF (3-state) models
# ============================================================================
# The Pikal model has 2 states: u = [hf (cm), Tf (K)].
# The RF model has 3 states: u = [mf (g), Tf (K), Tvw (K)]. Because the RF
# model's third state is a temperature, it additionally supports TvwSeriesData
# and TvwEndData (which read sol[3, ...]).
#
# Both models are built from the same shared parameter values where possible;
# only the RF-specific fields (vial mass, fill mass, RF power, dielectric and
# coupling constants) differ.

vialsize = "6R"
rad_i, rad_o = get_vial_radii(vialsize)
Ap = π*rad_i^2
Av = π*rad_o^2
csolid = 0.06u"g/mL"
ρsolution = 1u"g/mL"
Rp = RpFormFit(0.8u"cm^2*Torr*hr/g", 14.0u"cm*Torr*hr/g", 1.0u"1/cm")
Vfill = 3u"mL"
pch = RampedVariable(70u"mTorr")
Tsh = RampedVariable([-15u"°C", 10u"°C"].|>u"K", 0.5u"K/minute")
Kshf = RpFormFit(2.75e-4u"cal/s/K/cm^2", 8.93e-4u"cal/s/K/cm^2/Torr", 0.46u"1/Torr")
hf0 = Vfill/Ap

# --- Pikal model (2 states) ---
po = ParamObjPikal((
    (Rp, hf0, csolid, ρsolution),
    (Kshf, Av, Ap),
    (pch, Tsh)
))

# --- RF model (3 states), reusing the shared values above ---
m_v_rf = get_vial_mass(vialsize)
m_f0_rf = Vfill * ρsolution
Bf_rf = 2.0e7u"Ω/m^2"
Bvw_rf = 0.9e7u"Ω/m^2"
Kvwf_rf = 10.0u"W/K/m^2"
f_RF_rf = 8u"GHz"
P_per_vial_rf = RampedVariable(10u"W"/17 * 0.54)
po_rf = ParamObjRF((
    (Rp, hf0, csolid, ρsolution),
    (Kshf, Av, Ap),
    (pch, Tsh, P_per_vial_rf),
    (m_f0_rf, LyoPronto.cp_ice, m_v_rf, LyoPronto.cp_gl),
    (f_RF_rf, LyoPronto.eppf, LyoPronto.epp_gl),
    (Kvwf_rf, Bf_rf, Bvw_rf),
))

# Pre-interpolated solutions (saved at the experimental time points).
sol_conv = solve(ODEProblem(po), LyoPronto.odealg_chunk2; saveat=ustrip.(u"hr", t1))
sol_rf = solve(ODEProblem(po_rf), LyoPronto.odealg_chunk2; saveat=ustrip.(u"hr", t1))
# Non-interpolated solutions (not saved at the experimental time points).
sol_conv_ni = solve(ODEProblem(po), LyoPronto.odealg_chunk2)
sol_rf_ni = solve(ODEProblem(po_rf), LyoPronto.odealg_chunk2)

# Per-solution test data. Each entry is
#   (name, paramobj, sol, efd)
solutions = (
    ("Pikal", po, sol_conv, ExpFitData(t1, TfData(T1a), EndTimeData(t_end))),
    ("Pikal no-saveat", po, sol_conv_ni, ExpFitData(t1, TfData(T1a), EndTimeData(t_end))),
    ("RF", po_rf, sol_rf, ExpFitData(t1, TfData(T1a), TvwSeriesData(T2), EndTimeData(t_end))),
    ("RF no-saveat", po_rf, sol_rf_ni, ExpFitData(t1, TfData(T1a), TvwSeriesData(T2), EndTimeData(t_end))),
)

# --- Tests for obj_exp, err_exp, err_exp! (both models) ---
@testset "obj_exp, err_exp, err_exp! for $name" for (name, po_m, sol, efd) in solutions
    # Pikal: TfData (5) + EndTimeData (1) = 6.
    # RF: TfData (5) + TvwSeriesData (3) + TvwEndData (1) + EndTimeData (1) = 10.
    n = num_errs(efd)
    @test n == (occursin("Pikal", name) ? 6 : 9)

    # obj_exp should compute a finite, non-negative objective
    obj_val = obj_exp(sol, efd)
    @test obj_val isa Float64
    @test isfinite(obj_val)
    @test obj_val >= 0.0

    # obj_exp with Val(NaN) should return Inf
    @test obj_exp(Val(NaN), efd) == Inf

    # err_exp should return finite residuals of the correct length
    errs = err_exp(sol, efd)
    @test length(errs) == n
    @test all(isfinite, errs)

    # err_exp! in-place version should match
    errs2 = zeros(n)
    err_exp!(errs2, sol, efd)
    @test errs2 ≈ errs

    # err_exp! with Val(NaN) fills with Inf
    errs3 = zeros(n)
    err_exp!(errs3, Val(NaN), efd)
    @test all(isinf, errs3)

    # err_exp with Val(NaN)
    errs4 = err_exp(Val(NaN), efd)
    @test all(isinf, errs4)

    # Wrong-length errs should throw
    @test_throws ErrorException err_exp!(zeros(n-1), sol, efd)
    @test_throws ErrorException err_exp!(zeros(n+1), sol, efd)

    # obj_exp with per-dimension weights
    obj_weighted = obj_exp(sol, efd;
        weights=loss_weighting(Tf=2.0u"K^-2", t=2.0u"hr^-2"))
    @test isfinite(obj_weighted)
    @test obj_weighted >= 0.0

    # err_exp with tweight kwarg
    errs_tw = err_exp(sol, efd; weights=residual_weighting(t=1.0u"hr^-1"))
    @test length(errs_tw) == n
    @test all(isfinite, errs_tw)
    @test ~all(iszero, errs_tw)

    # obj_exp with verbose (should not error, but should log some info)
    obj_verbose = @test_logs (:info, r"loss call") match_mode=:any obj_exp(sol, efd; verbose=true)
    @test isfinite(obj_verbose)
end

# Depends on the `solutions` constructed for the above test
@testset "compare non-interp to interp" begin
    obj_interp = obj_exp(solutions[1][3], solutions[1][4])
    obj_ni = obj_exp(solutions[2][3], solutions[2][4]) 
    @test isapprox(obj_interp, obj_ni; rtol=1e-8)

    obj_interp = obj_exp(solutions[3][3], solutions[3][4])
    obj_ni = obj_exp(solutions[4][3], solutions[4][4]) 
    @test isapprox(obj_interp, obj_ni; rtol=1e-8)

end

@testset "Extensible residual weighting" begin
    struct PressureDatum <: LyoPronto.AbstractExpDatum end
    LyoPronto.resid_name(::PressureDatum) = :pe
    LyoPronto.time_bound_data(::PressureDatum) = false
    LyoPronto.obj_exp_datum(sol, dat::PressureDatum; verbose=false) = 1.0u"Pa^2"
    pressure_weight = LyoPronto.loss_weighting(pe=3.0u"Pa^-2")
    @test pressure_weight[LyoPronto.resid_name(PressureDatum())] == 3.0u"Pa^-2"
    @test LyoPronto.obj_exp(sol_conv, ExpFitData((1:5)u"hr", PressureDatum()); weights=pressure_weight) == 3.0
end
