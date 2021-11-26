include("fit_model.jl")
@info "Included model code"

using CSV, DifferentialEquations,
        DataFrames, BlackBoxOptim, Plots, 
		FileIO, JLD2, DataFramesMeta, Statistics,
		DiffEqParamEstim, Optim, StatsPlots,
		Distributions, DiffEqBayes

@info "Packages added"

#Lectura de la base de datos
base_datos    = "base_julia.csv"
datos_todos   = DataFrame(CSV.File(base_datos))

#Fecha inicial
#datos_todos = @subset(datos_todos, :t .< 100)
dia_min     = minimum(datos_todos[!,:t])
dia_max     = maximum(datos_todos[!,:t]) 

#Tiempo a fitear
tspan = (convert(Float64, dia_min), convert(Float64, dia_max))

#Quitamos los missings por si las dudas
candidate_params = repeat([0.0], 21)
#candidate_params[length(candidate_params)] = 350
#candidate_params[length(candidate_params) - 1] = 100
eval_params = 0

#Datos iniciales
datos_init = @subset(datos_todos, :t .== dia_min)

lsim   = convert(Int64, maximum(tspan))
t      = collect(range(1.0, stop=lsim, length=lsim - 1))

#Búsqueda 
S = 1 - 1.3*datos_init[1,:Positivos] - datos_init[1,:Hospitalizados] - 
	datos_todos[7,:Positivos] - datos_init[1,:UCI] - datos_init[1,:Defunciones]
initial_state = vcat([S, datos_todos[7,:Positivos], 
						0.3*datos_init[1,:Positivos], 
						datos_init[1,:Positivos], 
						datos_init[1,:Hospitalizados], 
						datos_init[1,:UCI], 
						datos_init[1,:Defunciones], 
						datos_init[1,:Positivos], 
						datos_init[1,:Hospitalizados], 
						datos_init[1,:UCI]])

#Transformamos a datos
dats          = Matrix(datos_todos[:,[:Defunciones,:Positivos,:Hospitalizados,:UCI]]) 						

ode_algorithm = ImplicitEuler()
prob          = ODEProblem(fit_model, initial_state, tspan)
obj           = build_loss_objective(prob, ode_algorithm, L2Loss(t, dats'), save_idxs = [7,8,9,10], maxiter = 1_000)


#Search objective function
search_obj    = repeat([(0.0,1.0)], length(candidate_params))
#search_obj[length(candidate_params) - 1] = (-180.0,180.0)
#search_obj[length(candidate_params)] = (0.0,365.0)

#Optimize
prob_bb  = bboptimize(obj,
				SearchRange = search_obj,
				NumDimensions = length(search_obj),
				Method = :adaptive_de_rand_1_bin_radiuslimited,
				MaxSteps = 10_000,
				NThreads = Threads.nthreads() - 2,
				TraceMode = :verbose)
parameter_nm = best_candidate(prob_bb)
eval_nm      = best_fitness(prob_bb)


min_vals = fill(0.0, convert(Int64, length(parameter_nm)))
max_vals = fill(1.0, convert(Int64, length(parameter_nm)))
#min_vals[length(candidate_params) - 1] = -180.0
#max_vals[length(candidate_params) - 1] = 180.0
#max_vals[length(candidate_params)] = 365.0

optim_space   = Optim.optimize(obj, min_vals, max_vals, parameter_nm,
				Fminbox(LBFGS()),
				Optim.Options(show_trace = true, iterations = 100, outer_iterations = 100))
parameter_optim  = Optim.minimizer(optim_space)
eval_optim       = Optim.minimum(optim_space)


@save "parameters_paper.jld2" parameter_optim

ode_algorithm = Rosenbrock32()
prob_bfgs     = ODEProblem(fit_model, initial_state, tspan, parameter_optim)
sol_bfgs      = solve(prob_bfgs, ode_algorithm)
Plots.plot(sol_bfgs, vars = (0,7), title = "Muertos", ylabel = "M(t)", xlabel = "t")
@df datos_todos StatsPlots.scatter!(:t, :Defunciones)
png("intento2.png")

Plots.plot(sol_bfgs, vars = (0,9), title = "Hospitalizados", ylabel = "I1(t)", xlabel = "t")
@df datos_todos StatsPlots.scatter!(:t, :Hospitalizados)
png("hosp.png")

Plots.plot(sol_bfgs, vars = (0,8), title = "Pos", ylabel = "I1(t)", xlabel = "t")
@df datos_todos StatsPlots.scatter!(:t, :Positivos)
png("pos.png")

#https://storopoli.io/Bayesian-Julia/pages/12_epi_models/

priors = [truncated(Normal(parameter_optim[1],0.1),0.0,1.0),
	truncated(Normal(parameter_optim[2],0.1),0.0,1.0),
	truncated(Normal(parameter_optim[3],0.1),0.0,1.0),
	truncated(Normal(parameter_optim[4],0.1),0.0,1.0),
	truncated(Normal(parameter_optim[5],0.1),0.0,1.0),
	truncated(Normal(parameter_optim[6],0.1),0.0,1.0),
	truncated(Normal(parameter_optim[7],0.1),0.0,1.0),
	truncated(Normal(parameter_optim[8],0.1),0.0,1.0),
	truncated(Normal(parameter_optim[9],0.1),0.0,1.0),
	truncated(Normal(parameter_optim[10],0.1),0.0,1.0),
	truncated(Normal(parameter_optim[11],0.1),0.0,1.0),
	truncated(Normal(parameter_optim[12],0.1),0.0,1.0),
	truncated(Normal(parameter_optim[13],0.1),0.0,1.0),
	truncated(Normal(parameter_optim[14],0.1),0.0,1.0),
	truncated(Normal(parameter_optim[15],0.1),0.0,1.0),
	truncated(Normal(parameter_optim[16],0.1),0.0,1.0),
	truncated(Normal(parameter_optim[17],0.1),0.0,1.0),
	truncated(Normal(parameter_optim[18],0.1),0.0,1.0),
	truncated(Normal(parameter_optim[19],0.1),0.0,1.0),
	truncated(Normal(parameter_optim[20],0.1),0.0,1.0),
	truncated(Normal(parameter_optim[21],0.1),0.0,1.0),
	truncated(Normal(parameter_optim[22],0.1),0.0,1.0),
	truncated(Normal(parameter_optim[23],0.1),0.0,1.0),
	truncated(Normal(parameter_optim[24],0.1),-180.0,180.0),
	truncated(Normal(parameter_optim[25],0.1),0.0,365.0)]

@info "Estimating parameters via Turing Inference"
num_samples = 1000
bayesian_result = turing_inference(prob, ode_algorithm, t, dats', priors,
                       num_samples = num_samples, epsilon = 0.1,
					   save_idxs = [7,8,9,10], seed = 238746,
                       progress = true, parallel = true)

muestra = DataFrame(bayesian_result)
CSV.write("Turing_result.csv", muestra)

function prob_func(prob_result,i,repeat)
	remake(prob_result, p = muestra[i, 3:27])
end

@info "Simulating posterior"
prob_bfgs     = ODEProblem(fit_model, initial_state, tspan, parameter_optim)
sol  = solve(EnsembleProblem(prob_bfgs, prob_func = prob_func), ode_algorithm, EnsembleSerial(), trajectories = 1000)

summ = EnsembleSummary(sol,quantiles=[0.025,0.975])


Plots.plot(summ, idxs = (3,),label = "M(t)", xlabel = "t")
@df datos_todos StatsPlots.scatter!(:t, :Defunciones)
png("Turing_Mort.png")
