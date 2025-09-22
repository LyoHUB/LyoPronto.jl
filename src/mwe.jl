using NonlinearSolve
function nl(du, u, p)
    du[1] = u[1] - p[1]
    du[2] = cos(u[1]) - cbrt(p[2])
    du[3] = u[2]^3 - p[1]
    du[4] = sin(u[2]) - sqrt(p[2])
end
function nl_oop(u, p) 
    du = similar(u, 4)
    nl(du, u, p)
    return du
end
nlf = NonlinearFunction{true}(nl, resid_prototype=zeros(4))
nlf = NonlinearFunction{false}(nl_oop)
prob = NonlinearLeastSquaresProblem(nlf, [0.0, 0.0], [1.0, 2.0])

sol1 = solve(prob, LevenbergMarquardt(), reltol=1e-2, abstol=1e-2)
sol2 = solve(prob, LevenbergMarquardt(), reltol=1e-8, abstol=1e-8)

sol1 == sol2 # true
sol1.stats == sol2.stats # false

x = range(0, 10, length=20)
y0 = rand(20) .+ x.^2
using Random
function varying_length(u, p)
    len = rand(19:29)
    pred_y = @. u[1] + u[2]*x + u[3]*x^2
    if len > length(x)
        return pred_y .- y0
    else
        return pred_y[1:len] .- y0[1:len]
    end
end

Random.seed!(1234)
nlvary = NonlinearLeastSquaresProblem{false}(varying_length, [0.0, 0.0, 0.0])
solvary = solve(nlvary, LevenbergMarquardt(), show_trace=Val(true))

