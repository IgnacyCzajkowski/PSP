using Graphs
using GraphPlot
using Setfield
using Plots
using Statistics
using StatsBase
using LinearAlgebra

using PrettyTables

#Struktura implementująca sieć
struct Network
    graph::Graph
    network_state::Matrix{Int}
    source_idx_matrix::Matrix{Int}
end

#Struktura implementująca obserwatora
mutable struct Observer
    idx::Int
    t::Int 
end

#Funkcja inicjalizująca początkowy stan sieci
function initialize_info_source(network_state::Matrix{Int})
    len = size(network_state)[2]
    inf_num = size(network_state)[1]
    source_indx_matrix::Matrix{Int} = zeros(inf_num, 1)
    for inf in 1:inf_num
        idx = rand(1:len)
        network_state[inf, idx] = 1
        source_indx_matrix[inf, 1] = idx
    end
    return source_indx_matrix 
end

#Funkcja tworząca sieć typu E-R o podanych parametrach i inicjalizująca stan początkowy sieci
function generate_network(nodes::Int, prob::Float16, inf_number::Int)
    G::Graph = erdos_renyi(nodes, prob)
    network_state::Matrix{Int} = zeros(inf_number, nodes)
    for idx in 1:nodes
        #Zapewnienie braku odizolowanych węzłów
        if length(all_neighbors(G, idx)) == 0
            rand_idx = idx
            while   rand_idx == idx
                rand_idx = rand(1:nodes)
            end     
            add_edge!(G, idx, rand_idx)
        end    
    end    
    source_idx_matrix::Matrix{Int} = initialize_info_source(network_state) 
    N = Network(G, network_state, source_idx_matrix)
    return N
end

#Funkcja tworząca sieć typu B-A o podanych parametrach i inicjalizująca stan początkowy sieci
function generate_network(nodes::Int, base_nodes::Int, edges::Int, inf_number::Int)
    G = barabasi_albert(nodes, base_nodes, edges, complete = true)
    network_state::Matrix{Int} = zeros(inf_number, nodes)
    source_idx_matrix::Matrix{Int} = initialize_info_source(network_state) 
    N = Network(G, network_state, source_idx_matrix)
    return N
end

#Funkcja do aktualizacji stanu węzła o indeksie indx w modelu SI
function  interact_witch_closest(N::Network, indx::Int, inf_prob_loc::Float64, inf_number::Int)
    #neighbors = all_neighbors(N.graph, indx)
    adjency_matrix = Matrix(Graphs.LinAlg.adjacency_matrix(N.graph))


    #=
    k_interact = 0
    if N.network_state[indx] == 0
        for i in neighbors
            if N.network_state[i] == 1
                k_interact = k_interact + 1
            end
        end
        if rand() < 1 - (1 - inf_prob_loc) ^ k_interact
            return 1
        end
    end
    =#
    return N.network_state[indx]
end

#Funkcja do aktualizacji stanu węzła o indeksie indx w modelu FSIR
function  interact_witch_closest_fsir(N::Network, inf_prob_vec::Matrix{Float64}, gamma::Float64, adjency_matrix)
    #=
    neighbors = all_neighbors(N.graph, indx)
    k_interact = 0
    if N.network_state[indx] == 0
        for i in neighbors
            if N.network_state[i] == 1
                k_interact = k_interact + 1
            end
        end
        if rand() < 1 - (1 - inf_prob_loc / k_interact^gamma) ^ k_interact
            return 1
        end
    end
    return N.network_state[indx]
    =#
    neighbours_matrix::Matrix{Int} =  N.network_state * adjency_matrix
    inf_overload_matrix = [count(x -> x != 0, col) for col in eachcol(neighbours_matrix .* (1 .- N.network_state))]
    inf_overload_matrix = [x == 0 ? 1 : x for x in inf_overload_matrix]
    inf_overload_matrix = (1 ./ inf_overload_matrix) .^ gamma
    #print(inf_overload_matrix) #
    
    
    #=
    print("state_matrix: ") #
    pretty_table(N.network_state; display_size=(10000, 10000))

    print("adjency_matrix: ") #
    pretty_table(adjency_matrix; display_size=(10000, 10000))

    print("neighbours_matrix: ") #
    pretty_table(neighbours_matrix; display_size=(10000, 10000))
    =#

    #0 to 1 and 1 to 0 flip: 1-x x{0,1}
    transition_matrix = 1 .- (1 .- inf_prob_vec) .^ neighbours_matrix
    transition_matrix = inf_overload_matrix' .* transition_matrix
    
    new_network_state::Matrix{Int} = N.network_state .+ (1 .- N.network_state) .* Int.(rand() .< transition_matrix)
    return new_network_state
end

#Funkcja do inicjalizacji bet z zadanego rozkladu normalnego
function generate_betas(inf_num::Int, mu::Float64, st::Float64)
    betas = randn(inf_num)
    betas = reshape(betas, inf_num, 1)
    betas = mu .+ st .* betas
    return 0.1 .+ 0.8 .* (betas .- minimum(betas)) ./ (maximum(betas) - minimum(betas))
end     

#Funkcja losująca i inicjalizująca wektor obserwatorów UWAGA: dodano tablice czasów obserwatorów
function getObservers(N::Network, l::Int)
    #obs = Vector{Observer}()
    obs_indxs = sample(1:size(N.network_state, 2), l, replace=false)
    #=
    for idx in a
        o = Observer(idx, Int(floatmax(Float16)))  
        if idx == N.source_idx
            o = Observer(idx, 0)   
        end
        push!(obs, o)
    end
    return obs
    =#
    observers_times_matrix = Int(floatmax(Float16)) .* transpose(1 .- N.network_state[:, obs_indxs])
    return obs_indxs, observers_times_matrix
end

#Funkcja służąca aktualizowanie stanu obserwatorów w trakcie symulacji
function actuateObservers(N_new::Network, N_old_state::Matrix{Int}, obs_indxs, observers_times_matrix, time::Int)
    #= 
    for point in obs
        if point.t == Int(floatmax(Float16))  
            if N_new.network_state[point.idx] == 1 && N_old_state[point.idx] == 0
                point.t = time
            end
        end
    end
    =#
    obs_changed_matrix = transpose(N_new.network_state[:, obs_indxs]) .+ transpose(N_old_state[:, obs_indxs]) .== 1
    observers_times_matrix[obs_changed_matrix] .= time
end

#Funkcja zwracająca wektor odległości od obserwatorów węzła o indeksie idx
function getDistanceFromObservers(N::Network, obs_indxs)
    d_all = Vector{Vector{Float64}}()
    for idx in 1:size(N.network_state, 2)
        d = Vector{Float64}()
        ds = desopo_pape_shortest_paths(N.graph, idx)
        for ob_indx in obs_indxs
            if ds.dists[ob_indx] > floatmax(Float16)
                push!(d, floatmax(Float16))
            else    
                push!(d, float(ds.dists[ob_indx]))
            end    
        end
        d = collect(d')
        d = reshape(d, 1, :)
        push!(d_all, vec(d))
    end    
    distances_matrix = vcat(d_all'...)
    return transpose(distances_matrix)
end

#Funkcja zwracająca wyniki algorytmu korelacyjnego dla wszystkich węzłów sieci
function getScore(distances_matrix, observers_times_matrix)
    score_matrix = cor(distances_matrix, observers_times_matrix)
    #=
    score = Vector{Float64}()
    t = Vector{Float64}()
    for point in obs
        push!(t, float(point.t))
    end
    for i in 1:length(N.network_state)
        d = getDistanceFromObservers(N, obs, i)
        sc::Float64 = cor(t, d)
        if isnan(sc)
            sc = -1.0   
        end
        push!(score, sc)
    end
    return score
    =#
    return score_matrix 
end

#Funkcja wyznaczająca precyzję i ranking w symulacji na podstawie wektora wyników
function analizeScore(N::Network, score_matrix)
    prec_vect = Vector{Float16}()
    rank_vect = Vector{Float16}()
    for (i, score_i) in enumerate(eachcol(score_matrix))
        solutions = Vector{Int}()
        for j in 1:length(score_i)
            if abs(score_i[j] - maximum(score_i)) < 0.001
                push!(solutions, j)
            end
        end

        src_score_i = score_i[N.source_idx_matrix[i]]
        if N.source_idx_matrix[i] in solutions
            prec = 1.0 / length(solutions)
        else
            prec = 0.0
        end
        rank = maximum(findall(x -> x == src_score_i, sort(score_i, rev=true)))
        
        push!(prec_vect, prec)
        push!(rank_vect, rank)
    end 
    return mean(prec_vect), mean(rank_vect) 
end

#Funkcja przeprowadzająca pojedyńczą symulacje UWAGA:T0 TRZEBA MOCNO ZMIENIC
function algorithm(N::Network, inf_prob_loc::Float64, gamma::Float64, observer_count::Int) 
    obs = getObservers(N, observer_count) 
    all_obs_infected::Bool = false
    time_step::Int = 1
    
    while all_obs_infected == false
        N_temp_vect = copy(N.network_state)
        N = @set N.network_state = get_next_step(N, inf_prob_loc, gamma)
        actuateObservers(N::Network, N_temp_vect, obs, time_step)
        all_obs_infected = true
        
        for observer in obs
            if observer.t == Int(floatmax(Float16)) 
                all_obs_infected = false
            end
        end    
        time_step += 1
        if time_step > 100000   #Warunek brzegowy symulacji
            println("Utknieto w pętli")
            break
        end
    end
   
    prec_kor, rank_kor = analizeScore(N, getScore(N, obs)) 
    return prec_kor, rank_kor, obs, N
      
end


#Funkcja resetująca stan sieci UWAGA: DO AKTUALIZACJI
function resetExistingNetwork(N::Network)
    new_network_state = Vector{Int}()
    for i in 1:length(N.network_state)
        if i == N.source_idx
            push!(new_network_state, 1)
        else
            push!(new_network_state, 0)
        end
    end
    N = @set N.network_state = new_network_state
end


### 2 OSTATNIE FUNKCJE DOPIERO PO NAPRAWIENIU POPRZEDNICH
#Funkcja przeprowadzająca serię symulacji w trybie beta i zapisująca wyniki symulacji do pliku data.txt