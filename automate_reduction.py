import networkx as nx
from sympy import symbols
from sympy.logic import SOPform

def get_minterm(d):
    """
    Sorts the dict's keys, then returns a sorted array of [0's and 1's] corresponding to the value of each key
    """
    keys = sorted(list(d.keys()))
    return [int(d[k]) for k in keys]

def hash_dict(d):
    """
    Sorts the dict's keys, then returns the dict in a standardized string form
    """
    keys = sorted(list(d.keys()))
    if len(keys)==0: return "{}"
    s = "{"
    for k in keys[:-1]:
        s += "%s: %s, "%(repr(k), repr(d[k]))
    k = keys[-1]
    s += "%s: %s}"%(repr(k), repr(d[k]))
    return s

def label2dict(label):
    """
    Converts labels like "{ph:0, Farnesol:1, Serum:0, Rapamycin:0}" to dict
    """
    terms = label[1:-1].split(',')
    d = dict()
    for term in  terms:
        key, val = term.split(":")
        d[key.strip()] = val.strip()
    return d

def check_label_equality(label1, label2):
    """
    Labels may correspond to the same dict/values, but come in a different order. This turns them into dicts and compares dict equality
    """
    return label2dict(label1)==label2dict(label2)

def check_branch_equality(G, node1, node2):
    """
    Returns True if the graph topology/labels that are a child of node1 is equivalent to that of node2
    """
    out_edges_1 = [e[1] for e in G.edges(node1)]
    out_edges_2 = [e[1] for e in G.edges(node2)]
    # If the nodes are sinks, their branches don't exist, and are equal
    if len(out_edges_1)==0 and len(out_edges_2)==0: return True

    # If the nodes have different numbers of children, the branches can't be equal
    if not len(out_edges_1)==len(out_edges_2): return False

    # Can we come up with a 1:1 mapping of children from node 1 with children from node 2?
    mapping = dict()
    for e1 in out_edges_1:
        l1 = G.nodes[e1]['label']
        for e2 in out_edges_2:
            l2 = G.nodes[e2]['label']
            if check_label_equality(l1, l2):
                mapping[e1]=e2
                break
        # We didn't find a node in branch2 corresponding to e1 in branch1
        if not e1 in mapping:
            return False
    
    # For each paired child, check if their branches are equivalent
    for e1 in out_edges_1:
        e2 = mapping[e1]
        if not check_branch_equality(G, e1, e2): return False

    # Otherwise, the branches are equal
    return True

G = nx.read_graphml("motif_network_no_attractor.graphml")
unique_nodes = dict()
for node in G.nodes():
    key = hash_dict(label2dict(G.nodes[node]['label']))
    if key in unique_nodes:
        unique_nodes[key].append(node)
    else:
        unique_nodes[key] = [node,]

# --------------------------------------------------------
# SM successions start with source nodes. Here I want to aggregate all source node motifs that lead to identical successions
source_nodes = sorted([i for i in G.nodes() if G.in_degree()[i]==0])
equality_graph = nx.Graph()

for _i in range(len(source_nodes)-1):
    n1 = source_nodes[_i]
    for _j in range(_i+1, len(source_nodes)):
        n2 = source_nodes[_j]
        if check_branch_equality(G, n1, n2):
            equality_graph.add_edge(n1, n2)

aggregated_source_nodes = list(nx.connected_components(equality_graph))
asn_sop = []

for asn in aggregated_source_nodes:
    minterms = []
    for node in asn:
        asn_dict = label2dict(G.nodes[node]['label'])
        minterms.append(get_minterm(asn_dict))
    syms = sorted(list(asn_dict.keys()))
    asn_sop.append(SOPform(symbols(syms), minterms))
    print(asn_sop[-1])
    
# --------------------------------------------------------
# Some stable motifs may not depend on input nodes. These should
# have identical branches, regardless of which asn they stem
# from. In these cases, they can be removed and put as top-level
# stable motifs

# I will pick one asn, pick a random source node (they should all be identical within an asn). I will loop over all its child nodes. Then I will loop over each other asn, and pick a random source node (again, they should be identical). 

found_matching_motifs = dict() # Keep track of whether or not a matching asn has been found

asn_reference = list(aggregated_source_nodes[0])
for edge in G.edges(asn_reference[0]): # Loop over 2nd level nodes
    node = edge[1]
    found_matching_motifs[node] = []

    for _i in range(1, len(aggregated_source_nodes)):
        found_matching_motifs[node].append(False)
        asn = list(aggregated_source_nodes[_i])
        for edge2 in G.edges(asn[0]):
            node2 = edge2[1] # This is a 2nd level node in another asn

            # Are the 2nd level nodes the same?
            if check_label_equality(G.nodes[node]['label'], G.nodes[node2]['label']):
                # Are their sub-branches the same?
                if check_branch_equality(G, node, node2):
                    found_matching_motifs[node][_i-1] = True
                    break
        if not found_matching_motifs[node][_i-1]:
            break

unconditional_motifs = []
for key in found_matching_motifs:
    if not False in found_matching_motifs[key]:
        unconditional_motifs.append(key)

###########################################################
# Aggregate it all into a reduced SM diagram
reducedG = nx.DiGraph()
for sop in asn_sop: reducedG.add_node(repr(sop),label=repr(sop))

# Add branches from the unconditional_motifs
subgraphs = []
for node in unconditional_motifs:
    subgraphs.append(nx.subgraph(G, [node,]+list(nx.descendants(G,node))))

# Add branches from the asm's
# Make sure to exclude branches that correspond to the unconditional subgraphs
for asm, sop in zip(aggregated_source_nodes, asn_sop):
    at_least_one_out_edge = False
    for edge in G.edges(list(asm)[0]): # These point to the level 2 motifs
        node = edge[1]
        # Make sure it is not unconditional
        unconditional = False
        for unc in unconditional_motifs:
            if check_label_equality(G.nodes[node]['label'], G.nodes[unc]['label']):
                unconditional=True
                break
        if unconditional: continue
        reducedG.add_edge(repr(sop), node)
        at_least_one_out_edge = True
        subnodes = [node,]
        subnodes += list(nx.descendants(G,node))
        subgraphs.append(nx.subgraph(G, subnodes))
    if not at_least_one_out_edge:
        reducedG.remove_node(repr(sop))

for subgraph in subgraphs:
    reducedG = nx.compose(reducedG, subgraph)

nx.write_graphml(reducedG, "autoreduced.graphml")