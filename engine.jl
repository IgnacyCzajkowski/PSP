#UWAGA: DO AKTUALIZACJI 
#Wczytywanie parametrów serii symulacji z pliku params.txt
params_array = Vector{}()
file = open("params.txt", "r")
for line in readlines(file)
    data = split(split(line, "#")[1], " ")
    push!(params_array, data)
end  
close(file) 
if length(params_array[1]) == 2
    erdos_reyni::Bool = true
    n = parse(Int, params_array[1][1])
    p = parse(Float64, params_array[1][2])
    network_params = [n, p]
elseif length(params_array[1]) == 3
    erdos_reyni = false
    n = parse(Int, params_array[1][1])
    n0 = parse(Int, params_array[1][2])
    k = parse(Int, params_array[1][3])
    network_params = [n, n0, k]
end

observer_count::Int = parse(Int, params_array[2][1])
if params_array[3][1] == "gamma"
    use_gamma::Bool = true
    beta = parse(Float64, params_array[4][1])
    gamma_start = parse(Float64, params_array[5][1])
    gamma_step = parse(Float64, params_array[5][2])
    i_max = parse(Int, params_array[5][3])
elseif params_array[3][1] == "beta"
    use_gamma = false
    beta_start = parse(Float64, params_array[4][1])
    beta_step = parse(Float64, params_array[4][2])
    i_max = parse(Int, params_array[4][3])
    gamma = parse(Float64, params_array[5][1])
end

j_max::Int = parse(Int, params_array[6][1])

#Wywołanie odpowiedniej serii symulacji 
if use_gamma == true
    main_gamma(beta, network_params, observer_count, gamma_start, gamma_step, i_max, j_max)
elseif use_gamma == false
    main_beta(gamma, network_params, observer_count, beta_start, beta_step, i_max, j_max)    
end