
# Virus: no, low, high
# Virus = 2


import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.cm as cm
import seaborn as sns
import copy

def BN(components, x, operations, conv_iter=5, n_iter=250):
    '''Generates a boolean network based on the provided components

    Parameters
    ----------
    components : list
    conv_iter : int, optional
        Number of iterations to run the network for to "simulate" convergence (default is 5)
    n_iter : int, optional
        Number of samples to take for each iteration (default is 250)
        Averaged at the end of the function
    
    Raises
    ----------
        NotImplementedError
            If the component is not implemented
            If an 'if' statement was not executed (check names)
        ValueError
            If the number of orders or iterations is less than 1
            If the number of orders or iterations is not an integer
    '''

    if components is None:
        raise NotImplementedError("No components provided")
    if (conv_iter < 1) or (n_iter < 1):
        raise ValueError("Number of iterations and orders must be positive")
    if (conv_iter % 1 != 0) or (n_iter % 1 != 0):
        raise ValueError("Number of iterations and orders must be integers")
    
    temp_mat = np.zeros((n_iter, len(components), conv_iter))

    for order in np.arange(n_iter):
        x = {component: False for component in components}
        x["Virus"] = True
        x["ACE2"] = True

        for i in np.arange(n_iter):
            indices = np.random.permutation(len(components))
            for idx in indices:
                key = components[idx]
                if key in operations:
                    operations[key]()
                else:
                    raise NotImplementedError(f"Operation for '{key}' not defined in 'operations'")
            indices = np.sort(indices)
            temp_mat[i, :, order] = [x[components[a]] for a in indices]

    return temp_mat

def geneticalgorithm(operations, n_iter=100):
    '''Randomly generates a boolean network based on the provided operations'''
    if operations is None:
        raise NotImplementedError("No operations provided")
    if n_iter < 1:
        raise ValueError("Number of iterations must be positive")
    if n_iter % 1 != 0:
        raise ValueError("Number of iterations must be an integer")
    

def generate_random_function(operations):
    KEYS = list(operations.keys())
    OPERATORS = ["not", "or", "and", "(", ")"]

    function = []
    open_paren_save = []
    
    while len(function) < 6:  # change 6 to modify function length
        if function and function[-1] not in OPERATORS:  # 30% chance to add parenthesis
            function.append(")")
            if open_paren_save:
                function.insert(open_paren_save.pop(), "(")
            
        token_type = np.random.choice(["KEY", "OPERATOR"])
        if token_type == "KEY" and (not function or function[-1] in OPERATORS + ["("]):
            function.append(KEYS[np.random.randint(0, len(KEYS))])
        elif token_type == "OPERATOR" and function and function[-1] not in OPERATORS + ["("]:
            function.append(OPERATORS[np.random.randint(0, len(OPERATORS))])
            if function[-2] == "(":
                open_paren_save.append(len(function) - 2)  # save the position of open parenthesis to match them later
    
    # Make sure all parentheses are closed
    while len(open_paren_save):
        function.pop()  # remove the last opening parenthesis
        function.pop(open_paren_save.pop())  # remove the corresponding closing parenthesis
        
    return function

def create_lambda_function(operations, func):
    KEYS = list(operations.keys())
    OPERATORS = ["not", "or", "and", "(", ")"]

    expr = ""
    for i in range(0, len(func)):
        if func[i] in KEYS:  # key
            expr += 'x["{}"]'.format(func[i])
        elif func[i] in OPERATORS:  # operator
            if func[i] == "not":
                expr = "not " + expr  # prepend not to the expression
            else:
                expr += ' ' + func[i] + ' '
        elif func[i] == "(" or func[i] == ")":  # parentheses
            expr += func[i]
    return lambda x: x.__setitem__(func[0], eval(expr))

def fitness_function(operations, func, target_output):
    KEYS = list(operations.keys())

    x = {}
    for key in KEYS:
        x[key] = True  # initial value
    lambda_func = create_lambda_function(operations, func)
    lambda_func(x)
    return abs(x[func[0]] - target_output)

def mutate(operations, child):
    KEYS = list(operations.keys())
    OPERATORS = ["not", "or", "and", "(", ")"]

    num_mutations = np.random.randint(1, len(child) // 2 + 1)
    for _ in range(num_mutations):
        mutation_point = np.random.randint(1, len(child))
        if mutation_point % 2 == 0:  # Mutate key
            child[mutation_point] = KEYS[np.random.choice(len(KEYS))]
        else:  # Mutate operator
            child[mutation_point] = OPERATORS[np.random.choice(len(OPERATORS))]

def generate_population(operations, size):
    return [generate_random_function(operations) for _ in range(size)]

def crossover(parent1, parent2):
    crossover_point = np.random.randint(1, len(parent1) - 1)
    child1 = parent1[:crossover_point] + parent2[crossover_point:]
    child2 = parent2[:crossover_point] + parent1[crossover_point:]
    return child1, child2

def run_genetic_algorithm(operations, population_size, num_generations, target_output):
    population = generate_population(operations, population_size)
    
    for generation in range(num_generations):
        print(f"Generation {generation}")
        
        # Evaluate fitness of population
        fitnesses = [fitness_function(operations, func, target_output) for func in population]
        
        # Select best individuals
        best_individuals = sorted(zip(population, fitnesses), key=lambda x: x[1])[:population_size//2]
        
        # Create next generation
        next_generation = []
        for individual, fitness in best_individuals:
            next_generation.append(individual)
            print(f"Fitness: {fitness}, Function: {individual}")
        
        # Crossover and mutation
        while len(next_generation) < population_size:
            parent1 = np.random.choice(best_individuals)[0]
            parent2 = np.random.choice(best_individuals)[0]
            child1, child2 = crossover(parent1, parent2)
            mutate(child1)
            mutate(child2)
            next_generation.append(child1)
            next_generation.append(child2)
        
        population = next_generation
    
    return best_individuals[0][0]  # Return the best individual from the last generation

def main():
    np.random.seed(0)
    fig, ax = plt.subplots(1, 1, figsize=(10, 7))
    plt.rcParams['figure.dpi'] = 300

    components = ["Virus", "Viral_Repl", "ACE2", "ANG_2", 
                  "ANG_2_T1R", "ADAM_17", "TLR4", "RIG1", "NFKB", "SIL6R",
                  "IKKB a/b", "TNF", "IRF3", "STAT1", "STAT3", "IL6", "IL6R",
                  "ISG", "C_FLIP", "IFN a/b", "NLRP3", "CASP1", "FOXO3A",
                  "IFNR", "BCL_2", "tBid", "Bax_Bak", "CASP9", "ROS", "TNFR",
                  "FADD", "Pyroptosis", "IL1", "IL1R", "MLKL", "Necroptosis",
                  "RIPK1&3", "CASP8", "Apoptosis", "PI3K", "AKT",
                  "TSC2", "mTORC1", "mTORC2", "CREB_1",
                  "Autophagy", "MAPK_p38", "HIF_1a"]
    components.sort()
    comp_edit = copy.deepcopy(components)

    x = {component: False for component in components}
    operations = {
        "ACE2": lambda: x.__setitem__("ACE2", not x["Virus"] or x["FOXO3A"]),
        "ADAM_17": lambda: x.__setitem__("ADAM_17", x["ANG_2_T1R"]),
        "AKT": lambda: x.__setitem__("AKT", x["mTORC2"] or x["PI3K"] or x["FOXO3A"]),
        "ANG_2": lambda: x.__setitem__("ANG_2", not x["ACE2"]),
        "ANG_2_T1R": lambda: x.__setitem__("ANG_2_T1R", x["ANG_2"]),
        "Apoptosis": lambda: x.__setitem__("Apoptosis", x["CASP8"] or x["CASP9"]),
        "Autophagy": lambda: x.__setitem__("Autophagy", not x["mTORC1"]),
        "BCL_2": lambda: x.__setitem__("BCL_2", (x["NFKB"] or x["CREB_1"] or x["HIF_1a"])),
        "Bax_Bak": lambda: x.__setitem__("Bax_Bak", not x["BCL_2"]),
        "CASP1": lambda: x.__setitem__("CASP1", x["NLRP3"]),
        "CASP8": lambda: x.__setitem__("CASP8", (x["FADD"] and not x["C_FLIP"]) and not x["FOXO3A"]),
        "CASP9": lambda: x.__setitem__("CASP9", x["Bax_Bak"] or x["tBid"]),
        "C_FLIP": lambda: x.__setitem__("C_FLIP", x["NFKB"] and not x["FOXO3A"]),
        "CREB_1": lambda: x.__setitem__("CREB_1", x["AKT"]),
        "FADD": lambda: x.__setitem__("FADD", x["TNFR"]),
        "FOXO3A": lambda: x.__setitem__("FOXO3A", (x["STAT3"] or x["MAPK_p38"]) and not (x["IKKB a/b"] or x["AKT"])),
        "HIF_1a": lambda: x.__setitem__("HIF_1a", (x["NFKB"] or x["mTORC1"] or x['STAT3']) and x["ROS"]), ###added STAT3
        "IFN a/b": lambda: x.__setitem__("IFN a/b", (x["IRF3"]) and not (x["Viral_Repl"])),
        "IFNR": lambda: x.__setitem__("IFNR", x["STAT1"]),
        "IL1": lambda: x.__setitem__("IL1", x["CASP1"] or x["CASP8"] or x["NFKB"]),
        "IL1R": lambda: x.__setitem__("IL1R", x["IL1"]),
        "IL6": lambda: x.__setitem__("IL6", x["NFKB"] or x["MAPK_p38"]),
        "IL6R": lambda: x.__setitem__("IL6R", x["IL6"]),
        "IRF3": lambda: x.__setitem__("IRF3", x["RIG1"]),
        "IKKB a/b": lambda: x.__setitem__("IKKB a/b", (x["TLR4"] or x["IL1R"] or x["TNFR"] or x["AKT"])),
        "ISG": lambda: x.__setitem__("ISG", x["STAT1"]),
        "MAPK_p38": lambda: x.__setitem__("MAPK_p38", (x["ANG_2_T1R"] or x["TLR4"] or x["ROS"])),
        "MLKL": lambda: x.__setitem__("MLKL", x["RIPK1&3"]),
        "mTORC1": lambda: x.__setitem__("mTORC1", not x["TSC2"]),
        "mTORC2": lambda: x.__setitem__("mTORC2", x["PI3K"]),
        "Necroptosis": lambda: x.__setitem__("Necroptosis", x["MLKL"]),
        "NFKB": lambda: x.__setitem__("NFKB", (x["IKKB a/b"] or x["ROS"]) and not x["FOXO3A"]),
        "NLRP3": lambda: x.__setitem__("NLRP3", x["NFKB"] and x["RIG1"]),
        "Pyroptosis": lambda: x.__setitem__("Pyroptosis", x["CASP1"]),
        "PI3K": lambda: x.__setitem__("PI3K", x["TLR4"] or x["ROS"] or x["IL6R"]),
        "RIG1": lambda: x.__setitem__("RIG1", x["Viral_Repl"]),
        "RIPK1&3": lambda: x.__setitem__("RIPK1&3", (x["RIG1"] or x["TLR4"] or x["FADD"])),
        "ROS": lambda: x.__setitem__("ROS", (x["ANG_2_T1R"] or x["MAPK_p38"]) and not x["FOXO3A"]),
        "SIL6R": lambda: x.__setitem__("SIL6R", x["ADAM_17"] or x["IL6"]),
        "STAT1": lambda: x.__setitem__("STAT1", x["IFNR"]),
        "STAT3": lambda: x.__setitem__("STAT3", x["IL6R"]),
        "tBid": lambda: x.__setitem__("tBid", x["CASP8"] or x["ROS"]),
        "TLR4": lambda: x.__setitem__("TLR4", x["Virus"]),
        "TNF": lambda: x.__setitem__("TNF", x["ADAM_17"] or x["NFKB"] or x["MAPK_p38"]),
        "TNFR": lambda: x.__setitem__("TNFR", x["TNF"]),
        "TSC2": lambda: x.__setitem__("TSC2", not x["AKT"]),
        "Viral_Repl": lambda: x.__setitem__("Viral_Repl", (x["Virus"] and x["mTORC2"]) and not x["ISG"]),
        "Virus": lambda: x.__setitem__("Virus", x["Virus"])
    }

    n_iter = 25
    orders = 250
    # mat = np.zeros((n_iter + 1, len(components), orders))

    best_function = run_genetic_algorithm(operations, 100, 50, True)
    print(best_function)
    assert 1 == 0

    mat = np.average(BN(comp_edit, x, operations, n_iter, orders), axis=2)

    yticklabels = [str(x) for x in 1 + np.arange(n_iter)]
    # # yticklabels.append("Average")

    cmap = cm.get_cmap('icefire_r')

    def truncate_colormap(cmap, minval=0.0, maxval=1.0, n=100):
        new_cmap = colors.LinearSegmentedColormap.from_list(
            'truncated_%s' % cmap.name,
            cmap(np.linspace(minval, maxval, n)))
        return new_cmap

    edited_comp = [x.replace("_", " ").replace("a/b", "α/β")
                   for x in components]
    mat_trunc = mat[0:10]
    yticklabels_trunc = yticklabels[0:10]
    ax = sns.heatmap(mat_trunc, cmap=truncate_colormap(cmap, 0.25, 0.75, n=200),
                     linewidths=.05, xticklabels=edited_comp,
                     yticklabels=yticklabels_trunc, vmin=0, vmax=1, alpha=1)
    ax.tick_params(axis='y', which='major', labelsize=10)

    # colorbar = ax.collections[0].colorbar
    # xmin, xmax, delta = np.min(mat), np.max(mat), 0.1
    # colorbar.set_ticks(np.arange(xmin, xmax + delta, delta))

    ax.set_xlabel('Model Component', fontsize=12)
    ax.set_ylabel('Iteration Number', fontsize=12)
    ax.set_title(
        f'Model Component Activation in COVID-19 with {orders} Samples',
        fontsize=14)

    plt.tight_layout()
    fig.savefig('plot.png')


if __name__ == '__main__':
    main()
    print("Complete")
