"""
Created on Fri Jun 21 10:52:11 2013

@author: Giulio Valentino Dalla Riva
@email: gvd16@uclive.ac.nz
"""
#generate a tree list following a birth and death process
#the process is conditioned to N leaves (with success)
#the list is long rang
def generate_trees_bdN(birth=0.57721,death=0.130357,N=49,rang=314):
    from dendropy import TaxonSet, TreeList, treesim
    taxa = TaxonSet()
    trees = TreeList()
    trees=[treesim.birth_death(birth, death, ntax=N, taxon_set=taxa, repeat_until_success=True) for x in range(rang)]
    return trees

#generate a tree list following a birth and death process
#the process goes on for time=T (with success)
#the list is long rang    
def generate_trees_bdT(birth=0.57721,death=0.130357,T=9.1596,rang=314):
    from dendropy import TaxonSet, TreeList, treesim
    taxa = TaxonSet()
    trees = TreeList()
    trees=[treesim.birth_death(birth, death, max_time=T, taxon_set=taxa, repeat_until_success=True) for x in range(rang)]
    return trees

#generate a tree list following a Kingman process with population size = Pop_size
#the list is long rang     
def generate_trees_Kingman(Pop_size=49,rang=314):
    from dendropy import TaxonSet, TreeList, treesim
    taxa = TaxonSet()
    trees = TreeList()
    trees=[treesim.pure_kingman(taxon_set=taxa, pop_size=Pop_size) for x in range(rang)]
    return trees

#given a tree we evolve, starting from the root, a trait {0,1} along the branches
#the trait evolves following a simple 2 state markov chain
#the probability of changing state from 0 to 1 is p_01
#the probability of changing state from 1 to 0 is p_10 
def evolve_markov_trait(tree, p_01,p_10):
    import numpy
    def zero(p_01):
        if numpy.random.rand() < p_01:
            return 1
        else:
            return 0
    def one(p_10):
        if numpy.random.rand() < p_10:
            return 0
        else:
            return 1
    markov_trait = {0 : zero,
                    1 : one
    }        
    for node in tree.preorder_node_iter():
        if node.parent_node is None:
            node.trait = numpy.random.randint(0,2)
        elif node.parent_node.trait == 0:
            node.trait = markov_trait[0](p_01)
        else:
            node.trait = markov_trait[1](p_10)
    return tree

#evolve a Markovian trait on all the trees of tree_list        
def evolve_markov_trait_list(tree_list,p_01=0.46692016,p_10=0.25029078):
    from dendropy import TreeList
    evolved_trees = TreeList()
    evolved_trees = [evolve_markov_trait(tree,p_01,p_10) for tree in tree_list]
    return evolved_trees

#given a tree we evolve, starting from the root, a vectors of binary traits
#each trait in the vector evolves following a simple 2 state markov chain
#the transition probabilities are given as a 2 times tract_length matrix   
def evolve_markov_traits(tree,tract_length=10, trans_prob=[[.1,.01] for i in range(10)]):
    import numpy
    def zero(p_01):
        if numpy.random.rand() < p_01:
            return 1
        else:
            return 0
    def one(p_10):
        if numpy.random.rand() < p_10:
            return 0
        else:
            return 1
    markov_trait = {0 : zero,
                    1 : one
    }        
    for node in tree.preorder_node_iter():
        if node.parent_node is None:
            node.traits = numpy.zeros(tract_length)
        else:
            node.traits = numpy.zeros(tract_length)
            for i in range(tract_length):
                if node.parent_node.traits[i] == 0:
                    p_01 = trans_prob[i][0]
                    node.traits[i] = markov_trait[0](p_01)
                else:
                    p_10 = trans_prob[i][1]
                    node.traits[i] = markov_trait[1](p_10)
    return tree

#evolve a vector of traits for a list of trees        
def evolve_markov_traits_list(tree_list,tract_length=10, trans_prob=[[.1,.01] for i in range(10)]):
    from dendropy import TreeList
    evolved_trees = TreeList()
    evolved_trees = [evolve_markov_traits(tree,tract_length,trans_prob) for tree in tree_list]
    return evolved_trees

#prune the tree following a random field of bullets model with extinction probability ext_prob
def prune_random(tree, ext_prob = .2):
    import numpy
    delete_list= [ t.taxon for t in tree.leaf_nodes() if numpy.random.rand() < ext_prob]
    tree.prune_taxa(delete_list)
    return tree

#FoB prune a list of trees
def prune_random_list(tree_list, ext_prob):
    from dendropy import TreeList, Tree
    trees_pruned = TreeList()
    def try_prune_random(tree):
        try:
            tree_pruned = prune_random(tree,ext_prob)
            return tree_pruned
        except AttributeError:
            t = Tree()
            return t
    trees_pruned = [try_prune_random(tree) for tree in tree_list]
    return trees_pruned    

#prune a singular trait evolved tree: the extinction probability of a leaf depend
#on its trait value and is specified by ext_0 and ext_1
def prune_trait(tree,ext_0=0.618033,ext_1=0.20205):
    import numpy
    delete_list= [ t.taxon for t in tree.leaf_nodes() if
    (t.trait == 0 and numpy.random.rand() < ext_0) or
    (t.trait == 1 and numpy.random.rand() < ext_1)]
    tree.prune_taxa(delete_list)
    return tree

#prune a list of singular trait evolved trees
def prune_trait_list(tree_list,ext_0=0.618033,ext_1=0.20205):
    from dendropy import TreeList, Tree
    trees_pruned = TreeList()
    def try_prune_tree(tree):
        try:
            tree_pruned = prune_trait(tree,ext_0,ext_1)
            return tree_pruned
        except AttributeError:
            t = Tree()
            return t
    trees_pruned = [try_prune_tree(tree) for tree in tree_list]
    return trees_pruned

#define a leaves extinction probability function based on markov traits   
def leaf_extinction_probability(leaf):
    extinction_probability = sum(leaf.traits)/len(leaf.traits)
    return extinction_probability
    
#prune the tree following using a leaf extinction probability based on traits
def prune_traits(tree):
    import numpy
    delete_list= [ t.taxon for t in tree.leaf_nodes() if numpy.random.rand() < leaf_extinction_probability(t)]
    tree.prune_taxa(delete_list)
    return tree

#prune the tree following using a leaf extinction probability based on traits
def prune_traits_threshold(tree):
    import numpy
    delete_list= [ t.taxon for t in tree.leaf_nodes() if numpy.random.rand() < leaf_extinction_probability(t)]
    tree.prune_taxa(delete_list)
    return tree
    
    
#prune a list of trees based on their traits
def prune_traits_list(tree_list):
    from dendropy import TreeList, Tree
    trees_pruned = TreeList()
    def try_prune_traits(tree):
        try:
            tree_pruned = prune_traits(tree)
            return tree_pruned
        except AttributeError:
            t = Tree()
            return t
    trees_pruned = [try_prune_traits(tree) for tree in tree_list]
    return trees_pruned 

#deep copy a list of trees
def copy_tree_list(tree_list_source):
    from dendropy import Tree
    target = [Tree(tree) for tree in tree_list_source]
    return target

#get the PD of a list of trees
def tree_list_lengths(tree_list):
    lengths = [tree.length() for tree in tree_list]   
    return lengths

#get the PD of a couple of lists of trees, i.e. original and pruned trees
def compare_pd_trees(tree_list_1,tree_list_2):
    pd_lists = ([tree_list_lengths(tree_list_1), tree_list_lengths(tree_list_2)])
    return pd_lists

#overall computation of birth death to N leaves singular trait evolved trees   
def do_the_mfob_bdN(birth=0.57721,death=0.130357,N=49,rang=314):
    trees_list = generate_trees_bdN(birth=birth,death=death,N=N,rang=rang)
    trees_evolved = evolve_markov_traits_list(trees_list)
    trees_for_pruning = copy_tree_list(trees_evolved)
    trees_pruned_traits = prune_traits_list(trees_for_pruning)
    trees_for_pruning = copy_tree_list(trees_evolved)
    Ext_Prob = [N - len(trees_pruned_traits[i].leaf_nodes()) for i in range(rang)]
    trees_pruned_random = prune_random_list(trees_for_pruning,ext_prob=Ext_Prob)
    pd_lists = compare_pd_trees(trees_pruned_random,trees_pruned_traits)
    return (trees_pruned_random, trees_pruned_traits, pd_lists)

#overall computation of birth death to time T singular trait evolved trees    
def do_the_mfob_Kingman(Pop_size=49,rang=314,p_01=0.46692016,p_10=0.25029078,ext_0=0.618033,ext_1=0.20205):
    t = generate_trees_Kingman(Pop_size=Pop_size,rang=314)
    t = evolve_markov_trait_list(t,p_01=p_01,p_10=p_10)
    t_original = copy_tree_list(t)
    t_pruned= prune_trait_list(t,ext_0=ext_0,ext_1=ext_1)
    pd_lists = compare_pd_trees(t_original,t_pruned)
    return (t_original, t_pruned, pd_lists)

#overall computation of Kingman singular trait evolved trees    
def do_the_mfob_bdT(birth=0.57721,death=0.130357,T=9.1596,rang=314,p_01=0.46692016,p_10=0.25029078,ext_0=0.618033,ext_1=0.20205):
    t = generate_trees_bdN(birth=birth,death=death,T=T,rang=rang)
    t = evolve_markov_trait_list(t,p_01=p_01,p_10=p_10)
    t_original = copy_tree_list(t)
    t_pruned = prune_trait_list(t,ext_0=ext_0,ext_1=ext_1)
    pd_lists = compare_pd_trees(t_original,t_pruned)
    return (t_original, t_pruned, pd_lists)
