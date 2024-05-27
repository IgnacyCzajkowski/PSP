include("utils.jl")
#using PrettyTables
#=
inf_num = 5
pop = 200
prob_er::Float16 = 0.3
#=
inf_prob_vec::Matrix{Float64} = zeros(5, 1)
inf_prob_vec[1] = 1.0
inf_prob_vec[2] = 0.8
inf_prob_vec[3] = 0.6
inf_prob_vec[4] = 0.3
inf_prob_vec[5] = 0.0
=#

N = generate_network(pop, 3, 2, inf_num)

inf_prob_vec = generate_betas(inf_num, 0.9, 2.0)
=#

function algorithm_test(N::Network, inf_prob_vec::Matrix{Float64}, gamma::Float64, observer_count::Int) 
    adjacency_matrix = Graphs.LinAlg.adjacency_matrix(N.graph)
    obs_indxs, observers_times_matrix = getObservers(N, observer_count)
    time_step::Int = 1
    while !all(x -> x != Int(floatmax(Float16)), observers_times_matrix)
        N_temp_state = copy(N.network_state)
        N = @set N.network_state = interact_witch_closest_fsir(N, inf_prob_vec, gamma, adjacency_matrix)
        actuateObservers(N, N_temp_state, obs_indxs, observers_times_matrix, time_step)
        time_step += 1
        if time_step > 100000   #Warunek brzegowy symulacji
            println("Utknieto w pętli")
            break
        end
    end
    distances_matrix = getDistanceFromObservers(N, obs_indxs)
    score_matrix = getScore(distances_matrix, observers_times_matrix)
    prec, rank = analizeScore(N, score_matrix)
    return prec, rank
end    

function main(betas_vect, network_params, observer_count::Int, gamma_start::Float64, gamma_step::Float64, i_max::Int, j_max::Int)
    file = open("data.txt", "w")
    for i in 1:i_max
        rank_avg_vect_kor = Vector{Float64}()
        prec_avg_vect_kor = Vector{Float16}()
        gamma_vect = Vector{Float64}()
        gamma = gamma_start + gamma_step * (i - 1)  
        push!(gamma_vect, gamma)
        for j in 1:j_max

            if length(network_params) == 2
                nodes::Int = network_params[1]
                prob::Float16 = network_params[2]
                N::Network = generate_network(nodes, prob, inf_num) 
            elseif length(network_params) == 3
                nodes = network_params[1]
                base_nodes::Int = network_params[2]
                edges::Int = network_params[3]
                N = generate_network(nodes, base_nodes, edges, inf_num) 
            end
                prec_kor, rank_kor = algorithm_test(N, betas_vect, gamma, observer_count)
                push!(prec_avg_vect_kor, prec_kor)
                push!(rank_avg_vect_kor, rank_kor)
                #resetExistingNetwork(N)
        end
            prec_avg_kor = sum(prec_avg_vect_kor) / length(prec_avg_vect_kor)
            std_dev_prec_kor = std(prec_avg_vect_kor) / length(prec_avg_vect_kor)
            rank_avg_kor = sum(rank_avg_vect_kor) / length(rank_avg_vect_kor)
            std_dev_rank_kor = std(rank_avg_vect_kor) / length(rank_avg_vect_kor)
            println("srednia Precyzja (korelacyjny):  ", prec_avg_kor, " +/-: ", std_dev_prec_kor)
            println("sredni Ranking (korelacyjny):  ", rank_avg_kor, " +/-: ", std_dev_rank_kor)

            write(file, string(gamma) * " " * string(prec_avg_kor) * " " * string(rank_avg_kor) * " " * string(std_dev_prec_kor) * " " * string(std_dev_rank_kor) * " \n")
    end
    close(file)
end

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
    betas_vect = Vector{Float64}()
    for beta_string in params_array[4]
        push!(betas_vect, parse(Float64, beta_string))
    end     
    inf_num = length(betas_vect)
    betas_vect = reshape(betas_vect, :, 1)
       
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

#algorithm_test(N, inf_prob_vec, 0.0)
main(betas_vect, network_params, observer_count, gamma_start, gamma_step, i_max, j_max)