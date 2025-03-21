using CairoMakie
using OrdinaryDiffEq
using LaTeXStrings

# Defining the parameters
α = 8e-10 #[mol L⁻¹ s⁻¹]
Kb = 5e-5 #[mol L⁻¹] (equivalent water concentration)

μₘ = 1e-5 #[s⁻¹]
Ka = 5e-4 #[mol L⁻¹]
Yd = 0.3 #[mol C_biomass/mol CD_substrate]
γa = 1/25 # stoichiometric coefficient of e-acceptor in the anabolic reaction
γc = 4/5 # stoichiometric coefficient of e-acceptor in the catabolic reaction
Ya = 1/(γa+(1/Yd-1)*γc) # [mol C_biomass/mol CA_substrate]
k_dec = 1.15e-6 #[s⁻¹]
Kd = 1e-6 #[mol L⁻¹]
η = 0.1 # [mol CD_substrate/mol C_biomass]
Ks = 0.22*(α/(k_dec*(1/Yd-η)) - Kb)

# Defining the ODE
retentostat!(du, u, p, t) = begin
    Cd, Ca, B, Ec, Eh, P, dPdEh = u
    rmax, rc, re, α, k_dec, Ke, Ka, Kd, Id, Kp, deltag, deltac, deltah, Q, V, c_in, cd_in = p
    rate = rmax * Ca / (Ca + Ka) * Cd / (Cd + Kd)*B * Ec
    ec_rate = rc * Ca / (Ca + Ka) * Cd / (Cd + Kd)*B * P/ (P + Kp)
    eh_rate = re * Ca / (Ca + Ka)* Id / (Cd + Id) * B  * P/ (P + Kp)
    hydrolysis = α * Ec / (Ec + Ke)
    du[1] = Q/V*cd_in -rate + hydrolysis - Q/V*Cd
    du[2] = Q/V*(c_in) - (4//5)*rate - Q/V*Ca
    du[3] = 0
    du[4] = ec_rate - k_dec*Ec
    du[5] = eh_rate - k_dec*Eh 
    du[6] = rate*deltag - ec_rate*deltac - eh_rate*deltah
    du[7] = α*Ke/(Eh + Ke)^2
end

u0 = [1e-8, 0., 1e-3, 0., 0., 1e-4, 0.]
p = [μₘ/Yd, μₘ*2, μₘ, α*5, k_dec*5, Kb, Ka, Kd, 5e-5, Kd, 1800, 2, 2, 6e-7, 1e-2, 1e-3, 0]
tspan = (0., 86400. *10)
prob = ODEProblem(retentostat!, u0, tspan, p)
sol = solve(prob, Rosenbrock23(), reltol=1e-16, abstol=1e-16)

# Plot the results
labels = [L"C_D", L"C_A", L"B", L"E_C", L"E_H", L"P", L"\frac{dP}{dE_H}"]
colors = [:blue, :red, :green, :purple, :orange, :black, :brown]
fig = Figure(size = (500, 700))
ax = Axis(fig[1, 1], xlabel="Time [hrs]", ylabel="Concentration [mol/L]")
ax2 = Axis(fig[2, 1], xlabel="Time [hrs]", ylabel="Concentration [mol/L]")
ax3 = Axis(fig[3, 1], xlabel="Time [hrs]", ylabel=L"\frac{dP}{dE_H}")
ax4 = Axis(fig[4, 1], xlabel="Time [hrs]", ylabel="Power [W]")
u = vcat(sol.u'...)
dPdEh = @. μₘ/Yd * 1800 * u[:, 2]/(u[:, 2] + Ka) * Kd/(u[:, 1] + Kd)^2 * u[:, 4] * u[:, 3] * u[:, 7] - 2 * μₘ * u[:, 2]/(u[:, 2] + Ka) * 5e-5 / (u[:,1] + 5e-5) * 1e-3 * u[:,7]/(u[:,7] + Kd)
for i in 1:2
    lines!(ax, sol.t./3600, u[:,i], label=labels[i], color=colors[i])
end
for i in 4:5
    lines!(ax2, sol.t./3600, u[:,i], label=labels[i], color=colors[i])
end
lines!(ax3, sol.t./3600, dPdEh, label=labels[7], color=colors[7])
lines!(ax4, sol.t./3600, u[:,6], label=labels[6], color=colors[6])
for axx in [ax, ax2, ax3, ax4]
   axislegend(axx, position=:rt)
end
resize_to_layout!(fig)
fig


condition(u, t, integrator) = t > 86400*0.5
affect!(integrator) = integrator.p[end] = 0
p[end] = 1e-3
cb = DiscreteCallback(condition, affect!;
    save_positions = (false, true),)
prob = ODEProblem(retentostat!, u0, tspan, p)
sol = solve(prob, Rosenbrock23(), callback=cb, reltol=1e-8, abstol=1e-8)

# Plot the results
labels = [L"C_D", L"C_A", L"B", L"E_C", L"E_H", L"P", L"\frac{dP}{dE_H}"]
colors = [:blue, :red, :green, :purple, :orange, :black, :brown]
fig = Figure(size = (500, 700))
ax = Axis(fig[1, 1], xlabel="Time [hrs]", ylabel="Concentration [mol/L]")
ax2 = Axis(fig[2, 1], xlabel="Time [hrs]", ylabel="Concentration [mol/L]")
ax3 = Axis(fig[3, 1], xlabel="Time [hrs]", ylabel=L"\frac{dP}{dE_H}")
ax4 = Axis(fig[4, 1], xlabel="Time [hrs]", ylabel="Power [W]")
u = vcat(sol.u'...)
dPdEh = @. μₘ/Yd * 1800 * u[:, 2]/(u[:, 2] + Ka) * Kd/(u[:, 1] + Kd)^2 * u[:, 4] * u[:, 3] * u[:, 7] - 2 * μₘ * u[:, 2]/(u[:, 2] + Ka) * 5e-5 / (u[:,1] + 5e-5) * 1e-3 * u[:,7]/(u[:,7] + Kd)
for i in 1:2
    lines!(ax, sol.t./3600, u[:,i], label=labels[i], color=colors[i])
end
for i in 4:5
    lines!(ax2, sol.t./3600, u[:,i], label=labels[i], color=colors[i])
end
lines!(ax3, sol.t./3600, dPdEh, label=labels[7], color=colors[7])
lines!(ax4, sol.t./3600, u[:,6], label=labels[6], color=colors[6])
for axx in [ax, ax2, ax3, ax4]
   axislegend(axx, position=:rt)
end
resize_to_layout!(fig)
fig