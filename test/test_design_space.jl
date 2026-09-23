
Ap, Av = π .*get_vial_radii("6R").^2  # cross-sectional area inside the vial
csolid = 0.06u"g/mL" # g solute / mL solution
ρsolution = 1u"g/mL" # g/mL total solution density
R0 = 0.8u"cm^2*Torr*hr/g"
A1 = 28.0u"cm*Torr*hr/g"
A2 = 1.0u"1/cm"
Rp = RpFormFit(R0, A1, A2)
# Cycle parameters
Vfill = 3u"mL" # ml
pch = RampedVariable(70u"mTorr")
Tsh = RampedVariable([-15u"°C", 10u"°C"].|>u"K", 0.5u"K/minute")
KC = 2.75e-4u"cal/s/K/cm^2"
KP = 8.93e-4u"cal/s/K/cm^2/Torr"
KD = 0.46u"1/Torr"
Kshf = RpFormFit(KC, KP, KD)
# Computed parameters based on above
hf0 = Vfill / Ap
po = ParamObjPikal((
    (Rp, hf0, csolid, ρsolution),
    (Kshf, Av, Ap),
    (pch, Tsh)
))

@testset "Design Space Tests" begin
    pch_range = 50u"mTorr":50u"mTorr":300u"mTorr"
    tab_sh = @time LyoPronto.ds_calculate_shelf_isotherm(20.0u"°C", pch_range, po)
    @testset "Shelf Isotherm Tests" begin

        @test LyoPronto.new_Tsh_ramp(po.Tsh, 5u"degC") == RampedVariable([po.Tsh.setpts[1], 5u"degC"|>u"K"], po.Tsh.ramprates[1])
        # Test the shelf isotherm produces po_i with all the same Tsh, different pch
        @test all([po.Tsh == tab_sh.po[1].Tsh for po in tab_sh.po])
        @test all(tab_sh.po .|> (x -> x.Tsh(Inf*u"hr")) .== 20.0u"°C"|>u"K")
        @test all(tab_sh.po .|> (x -> x.pch(0u"s")) .== pch_range)

        # Test that each solution has different statistics
        for name in [:max_md, :ave_md, :end_md, :max_Tf, :dry_time]
            vals = getproperty(tab_sh, name)
            @test length(unique(vals)) == length(vals)
        end
        # Test that max_md is greater than ave_md for all solutions
        @test all(tab_sh.max_md .> tab_sh.ave_md)
        # Test that end_md is less than max_md for all solutions
        @test all(tab_sh.end_md .< tab_sh.max_md)
    end

    Tcrit = -15.0u"°C"
    tab_pr = @time LyoPronto.ds_calculate_product_isotherm(Tcrit, pch_range, po)
    @testset "Product Isotherm Tests" begin
        # Test that each solution has different statistics, except max_Tf
        for name in [:max_md, :ave_md, :end_md, :dry_time]
            vals = getproperty(tab_pr, name)
            @test length(unique(vals)) == length(vals)
        end
        # Test that max_md is greater than ave_md for all solutions
        @test all(tab_pr.max_md .> tab_pr.ave_md)
        # Test that end_md is less than max_md for all solutions
        @test all(tab_pr.end_md .< tab_pr.max_md)

        # For Pikal model, flux should be monotonically decreasing with pch
        @test all(diff(tab_pr.max_md) .< 0u"g/hr")
        @test all(diff(tab_pr.ave_md) .< 0u"g/hr")
        @test all(tab_pr.max_Tf .≈ Tcrit)
    end
end
