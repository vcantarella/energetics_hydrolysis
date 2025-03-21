using DrWatson
@quickactivate "energetics_hydrolysis"

using Symbolics
using CairoMakie
using DifferentialEquations
using LinearSolve
using DataFrames


@variables α k_dec Yd η Kb Ya rss0 bss0 Er Eg Eh Em m Ca Cd Ka Kd t rm B


D = Differential(t)
D_Cd = α/(Kb + B) - rm * Cd/(Kd + Cd) * Ca/(Ka + Ca)

Symbolics.d