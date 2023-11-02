import os
import networkx as nx
import random
import math
import time


# read edges from file and calculate summary metrics
def read_edges_from_file(file_path):
    edges = []
    vertices = set()
    with open(file_path, 'r', encoding='utf-8') as file:
        lines = file.readlines()
        print(lines[0])
        for line in lines[1:]:
            values = line.strip().split()

            if len(values) == 2:
                u, v = values  
                vertices.add(u)
                vertices.add(v)
            else:
                print(f"Skipping line: {line.strip()}")
            if u == v:
                continue
            else:
                edges.append((u, v))
        E = len(edges)
        N = len(vertices)
        k = 2 * E / N
        delta = 2 * E / (N * (N - 1))
        print("Nr of edges: " + str(E))
        print("Nr of vertices: " + str(N))
        print("k: " + str(k))
        print("delta: " + str(delta))
    return edges


# create the network using networkx
def create_network(file_path):
    edges = read_edges_from_file(file_path)
    G = nx.Graph()
    G.add_edges_from(edges)
    return G


# calculate mean local clustering coefficient for network G
def calculate_mean_local_clustering_coefficient(G):
    return sum(nx.clustering(G).values()) / len(G)


def p_value_null_model(G, T, original_clus_coeff):
    '''
    Computes p-value for null model graph

    Parameters:
        G: input graph
        T: n. iterations for monte carlo estimation

    Return:
        p-value
    '''
    # Counter to count number of null model graphs where CW >= original CW
    counter, counter_org, counter_random, counter_asc, counter_desc = 0, 0, 0, 0, 0
    # Parameters for Erdős-Rényi model
    n = len(G.nodes)
    e = len(G.edges)
    M = math.ceil(0.1 * n)

    # Generate null model graphs and calculate their local clustering coefficients
    for _ in range(T):
        p = 2 * e / (n * (n - 1))
        st = time.time()
        null_model_graph = nx.erdos_renyi_graph(n, p)
        null_clustering_coefficient = calculate_mean_local_clustering_coefficient(null_model_graph)
        et = time.time()
        # get the execution time
        elapsed_time = et - st
        print('Execution time null_clustering_coefficient - considering the graph creation:', elapsed_time, 'seconds')
        print("Null Clustering Coefficient:", null_clustering_coefficient)

        if null_clustering_coefficient >= original_clus_coeff:
            counter += 1

        # original
        st = time.time()
        null_clustering_coefficient_org = mean_local_clustering_coefficient_opt(null_model_graph, M)
        et = time.time()
        # get the execution time
        elapsed_time = et - st
        print('Execution time null_clustering_coefficient_org:', elapsed_time, 'seconds')
        print("Null Clustering Coefficient optim original order:", null_clustering_coefficient_org)

        if null_clustering_coefficient_org >= original_clus_coeff:
            counter_org += 1

        # random
        st = time.time()
        null_clustering_coefficient_random = mean_local_clustering_coefficient_opt(null_model_graph, "rand", M)
        et = time.time()
        # get the execution time
        elapsed_time = et - st
        print('Execution time null_clustering_coefficient_random:', elapsed_time, 'seconds')
        print("Null Clustering Coefficient optim random order:", null_clustering_coefficient_random)

        if null_clustering_coefficient_random >= original_clus_coeff:
            counter_random += 1

        # asc
        st = time.time()
        null_clustering_coefficient_asc = mean_local_clustering_coefficient_opt(null_model_graph, "asc", M)
        et = time.time()
        # get the execution time
        elapsed_time = et - st
        print('Execution time null_clustering_coefficient_asc:', elapsed_time, 'seconds')
        print("Null Clustering Coefficient optim asc order:", null_clustering_coefficient_asc)

        if null_clustering_coefficient_asc >= original_clus_coeff:
            counter_asc += 1

        # desc
        st = time.time()
        null_clustering_coefficient_desc = mean_local_clustering_coefficient_opt(null_model_graph, "desc", M)
        et = time.time()
        # get the execution time
        elapsed_time = et - st
        print('Execution time null_clustering_coefficient_desc:', elapsed_time, 'seconds')
        print("Null Clustering Coefficient optim desc order:", null_clustering_coefficient_desc)

        if null_clustering_coefficient_desc >= original_clus_coeff:
            counter_desc += 1

    p_value = float(counter) / T
    print(
        f"Probability that Null Model Clustering coefficient is greater or equal than the original clustering coefficient is: {p_value: .4f}")
    print(
        f"Probability that Null Model Clustering coefficient with optim original order is greater or equal than the original clustering coefficient is: {float(counter_org) / T: .4f}")
    print(
        f"Probability that Null Model Clustering coefficient with optim random order is greater or equal than the original clustering coefficient is: {float(counter_random) / T: .4f}")
    print(
        f"Probability that Null Model Clustering coefficient with optim asc order is greater or equal than the original clustering coefficient is: {float(counter_asc) / T: .4f}")
    print(
        f"Probability that Null Model Clustering coefficient with optim desc order is greater or equal than the original clustering coefficient is: {float(counter_desc) / T: .4f}")

    return p_value


# switching model
def switch_edges(graph, Q):
    fail, success = 0, 0
    num_edges = len(graph.edges())

    adjacency_matrix = nx.to_numpy_array(graph, dtype=int)

    E = list(graph.edges())

    node_to_index = {node: i for i, node in enumerate(graph.nodes())}
    edgelist = [(node_to_index[u], node_to_index[v]) for u, v in graph.edges()]

    # Calculating all potential switches simultaneously
    edge1_id = random.choices(range(num_edges), k=Q * num_edges)
    edge2_id = random.choices(range(num_edges), k=Q * num_edges)

    # Switching for Q*E times
    for i in range(int(num_edges * Q)):
        edge1 = edgelist[edge1_id[i]]
        edge1_prev = edge1
        edge2 = edgelist[edge2_id[i]]
        edge2_prev = edge2

        # Check if switching is possible
        if (edge1[0] != edge2[0] and
                edge1[1] != edge2[1] and
                edge1[0] != edge2[1] and
                edge2[0] != edge1[1] and
                adjacency_matrix[edge2[0], edge1[1]] == 0 and
                adjacency_matrix[edge1[0], edge2[1]] == 0):

            edge1_end = edge1[1]
            edge1 = (edge1[0], edge2[1])
            edge2 = (edge2[0], edge1_end)

            adjacency_matrix[edge1_prev[0], edge1_prev[1]] = 0
            adjacency_matrix[edge1[0], edge1[1]] = 1
            adjacency_matrix[edge2_prev[0], edge2_prev[1]] = 0
            adjacency_matrix[edge2[0], edge2[1]] = 1
            adjacency_matrix[edge1_prev[1], edge1_prev[0]] = 0
            adjacency_matrix[edge1[1], edge1[0]] = 1
            adjacency_matrix[edge2_prev[1], edge2_prev[0]] = 0
            adjacency_matrix[edge2[1], edge2[0]] = 1

            edgelist[edge1_id[i]] = edge1
            edgelist[edge2_id[i]] = edge2

            success += 1
        else:
            fail += 1

    print("# Fails:", fail)
    print("# Success:", success)

    # Create a new graph from the updated edgelist
    new_graph = nx.Graph()
    new_graph.add_edges_from(edgelist)

    return new_graph


def p_value_rand_model(G, T, original_clus_coeff):
    '''
    Computes p-value for randomized model graph

    Parameters:
        G: input graph
        T: n. iterations for monte carlo estimation
        original_clus_coeff: clustering coefficient of original graph

    Return:
        p-value
    '''
    counter, counter_org, counter_random, counter_asc, counter_desc = 0, 0, 0, 0, 0
    N = len(G)
    E = len(G.edges())
    Q = math.ceil(math.log(E))  # adjusted according to coupon collector's problem Q >=logE
    M = math.ceil(0.1 * N)

    # Generate randomized graphs and calculate their local clustering coefficients
    for _ in range(T):
        st = time.time()
        rand_graph = switch_edges(G, Q)
        randomized_clustering_coefficient = calculate_mean_local_clustering_coefficient(rand_graph)
        if randomized_clustering_coefficient >= original_clus_coeff:
            counter += 1
        print("Random Local Clustering Coefficient:", randomized_clustering_coefficient)
        # get the end time
        et = time.time()
        # get the execution time
        elapsed_time = et - st
        print('Execution time of iteration {_}:', elapsed_time, 'seconds')

        # original
        st = time.time()
        randomized_clustering_coefficient_org = mean_local_clustering_coefficient_opt(rand_graph, "None", M)
        print("Random Clustering Coefficient optim original order:", randomized_clustering_coefficient_org)

        if randomized_clustering_coefficient_org >= original_clus_coeff:
            counter_org += 1
        et = time.time()
        # get the execution time
        elapsed_time = et - st
        print('Execution time of iteration {_} optim org:', elapsed_time, 'seconds')

        # random
        st = time.time()
        randomized_clustering_coefficient_random = mean_local_clustering_coefficient_opt(rand_graph, "rand", M)
        print("Random Clustering Coefficient optim random order:", randomized_clustering_coefficient_random)

        if randomized_clustering_coefficient_random >= original_clus_coeff:
            counter_random += 1
        et = time.time()
        # get the execution time
        elapsed_time = et - st
        print('Execution time of iteration {_} optim random:', elapsed_time, 'seconds')

        # asc
        st = time.time()
        randomized_clustering_coefficient_asc = mean_local_clustering_coefficient_opt(rand_graph, "asc", M)
        print("Random Clustering Coefficient optim asc order:", randomized_clustering_coefficient_asc)

        if randomized_clustering_coefficient_asc >= original_clus_coeff:
            counter_asc += 1
        et = time.time()
        # get the execution time
        elapsed_time = et - st
        print('Execution time of iteration {_} optim asc:', elapsed_time, 'seconds')

        # desc
        st = time.time()
        randomized_clustering_coefficient_desc = mean_local_clustering_coefficient_opt(rand_graph, "desc", M)
        print("Random Clustering Coefficient optim desc order:", randomized_clustering_coefficient_desc)

        if randomized_clustering_coefficient_desc >= original_clus_coeff:
            counter_desc += 1
        et = time.time()
        # get the execution time
        elapsed_time = et - st
        print('Execution time of iteration {_} optim desc:', elapsed_time, 'seconds')

    p_value = float(counter) / T
    print(
        f"Probability that Randomized Model Clustering coefficient is greater or equal than the original clustering coefficient is: {p_value: .4f}")
    print(
        f"Probability that Randomized Model Clustering coefficient with optim original order is greater or equal than the original clustering coefficient is: {float(counter_org) / T: .4f}")
    print(
        f"Probability that Randomized Model Clustering coefficient with optim random order is greater or equal than the original clustering coefficient is: {float(counter_random) / T: .4f}")
    print(
        f"Probability that Randomized Model Clustering coefficient with optim asc order is greater or equal than the original clustering coefficient is: {float(counter_asc) / T: .4f}")
    print(
        f"Probability that Randomized Model Clustering coefficient with optim desc order is greater or equal than the original clustering coefficient is: {float(counter_desc) / T: .4f}")

    return p_value


def mean_local_clustering_coefficient_opt(G, sorting=None, M=None):
    """
    Calculate the optimized version of mean local clustering coefficient

    Args:
        G: network of which compute mean local clustering coeff.
        sorting: type of sorting of the nodes
            (Original ordering,
            Random ordering of vertices (by generating a uniformly random permutation of the vertices),
            Increasing order by degree,
            Decreasing order by degree.)

    Returns:
        float: mean_local_clustering_coefficient
    """
    if M == None:
        M = math.ceil(0.1 * len(G))

    final_graph = None

    if sorting == "rand":
        # Step 1: Obtain the list of vertices
        rand_graph = G.copy()
        vertices = list(rand_graph.nodes())

        # Step 2: Randomly shuffle the list of vertices
        random.shuffle(vertices)

        # Step 3: Create a new graph with vertices in the random order
        # final_graph = rand_graph.subgraph(vertices)
        final_graph = nx.Graph()
        final_graph.add_nodes_from(vertices)
        final_graph.add_edges_from(rand_graph.edges(data=True))

    elif sorting == "asc":
        nodes_by_degree = sorted(G.degree, key=lambda x: x[1])
        sorted_nodes = [node for node, _ in nodes_by_degree]  # Extract node identifiers

        # Create a sorted list of edges based on the sorted nodes
        sorted_edges = [(u, v) for u, v in G.edges() if u in sorted_nodes and v in sorted_nodes]

        final_graph = nx.Graph()
        final_graph.add_nodes_from(sorted_nodes)
        final_graph.add_edges_from(sorted_edges)

    elif sorting == "desc":
        nodes_by_degree = sorted(G.degree, key=lambda x: x[1], reverse=True)
        sorted_nodes = [node for node, _ in nodes_by_degree]  # Extract node identifiers

        # Create a sorted list of edges based on the sorted nodes
        sorted_edges = [(u, v) for u, v in G.edges() if u in sorted_nodes and v in sorted_nodes]

        final_graph = nx.Graph()
        final_graph.add_nodes_from(sorted_nodes)
        final_graph.add_edges_from(sorted_edges)

    else:  # original ordering
        final_graph = G.copy()

    # Get the list of nodes from the original graph
    all_nodes = list(final_graph.nodes())

    # Create a subgraph containing only the first M nodes
    first_M_nodes = all_nodes[:M]
    subgraph = nx.Graph()
    subgraph.add_nodes_from(first_M_nodes)
    subgraph.add_edges_from(final_graph.edges(data=True))

    # Calculate the mean local clustering coefficient (CWS) for the first M nodes:
    mean_local_clustering_coefficient = sum(nx.clustering(subgraph).values()) / len(subgraph)
    return mean_local_clustering_coefficient


# For each network file, load it a create the network and the adj_matrix
networks = []
languages = []

folder_path = "data/"
# read the file and save the networks
for filename in os.listdir(folder_path):
        if filename.endswith('.txt'):
            file_path = os.path.join(folder_path, filename)
            language = filename.replace("_syntactic_dependency_network.txt", "")
            languages.append(language)
            print("")
            print(language)
            print(filename)
            network = create_network(file_path)
            networks.append(network)

for network in networks:
    print(len(network))
    # Calculate local clustering coefficient for the original graph
    original_clustering_coefficient = calculate_mean_local_clustering_coefficient(network)
    print("Original Mean Local Clustering Coefficient:", original_clustering_coefficient)
    print("-----------------------------------")
    # alpha = 0.05 # significance level
    # T>20 since 1/T << alpha if we want Monte Carlo procedures to work
    T = 21

    st = time.time()
    p_value_null = p_value_null_model(network, T, original_clustering_coefficient)
    et = time.time()

    elapsed_time = et - st
    print(f"p-value for null model: {p_value_null}")
    print('Execution time of null model:', elapsed_time, 'seconds')
    print("-----------------------------------")

    st = time.time()
    p_value_rand = p_value_rand_model(network, T, original_clustering_coefficient)
    et = time.time()
    elapsed_time = et - st
    print(f"p-value for rand model: {p_value_rand}")
    print('Execution time of rand model:', elapsed_time, 'seconds')
    print("-----------------------------------")
