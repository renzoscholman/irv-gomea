import networkx as nx
from pgmpy.models import FactorGraph
from pgmpy.factors.discrete import DiscreteFactor
from pgmpy.inference import BeliefPropagation
import matplotlib.pyplot as plt

# Step 1: Define the graph
G = nx.Graph()
G.add_nodes_from(['A', 'B', 'C', 'D', 'E', 'F'])
edges = [
    ('A', 'B'),
    ('B', 'C'),
    ('C', 'A'),
    ('B', 'D'),
    ('C', 'E'),
    ('E', 'F'),
    ('F', 'C'),
    ('D', 'E')
]
G.add_edges_from(edges)

# Draw the original graph
pos = nx.spring_layout(G)
nx.draw(G, pos, with_labels=True, node_color='lightblue', edge_color='gray')
plt.title("Original Graph with Cycles and Cliques")
plt.savefig("original_graph.png")
plt.close()

# Step 2: Create the factor graph
fg = FactorGraph()
variables = ['A', 'B', 'C', 'D', 'E', 'F']
fg.add_nodes_from(variables)

# Define factors based on true dependencies
phi_AB = DiscreteFactor(variables=['A', 'B'], cardinality=[2, 2],
                        values=[[0.9, 0.1],
                                [0.2, 0.8]])
phi_BC = DiscreteFactor(variables=['B', 'C'], cardinality=[2, 2],
                        values=[[0.3, 0.7],
                                [0.6, 0.4]])
phi_BD = DiscreteFactor(variables=['B', 'D'], cardinality=[2, 2],
                        values=[[0.5, 0.5],
                                [0.4, 0.6]])
phi_CE = DiscreteFactor(variables=['C', 'E'], cardinality=[2, 2],
                        values=[[0.7, 0.3],
                                [0.2, 0.8]])
phi_EF = DiscreteFactor(variables=['E', 'F'], cardinality=[2, 2],
                        values=[[0.6, 0.4],
                                [0.1, 0.9]])
phi_CA = DiscreteFactor(variables=['C', 'A'], cardinality=[2, 2],
                        values=[[0.8, 0.2],
                                [0.3, 0.7]])
phi_DE = DiscreteFactor(variables=['D', 'E'], cardinality=[2, 2],
                        values=[[0.4, 0.6],
                                [0.7, 0.3]])

fg.add_factors(phi_AB, phi_BC, phi_BD, phi_CE, phi_EF, phi_CA, phi_DE)

edges = [
    ('A', phi_AB),
    ('B', phi_AB),
    ('B', phi_BC),
    ('C', phi_BC),
    ('B', phi_BD),
    ('D', phi_BD),
    ('C', phi_CE),
    ('E', phi_CE),
    ('E', phi_EF),
    ('F', phi_EF),
    ('C', phi_CA),
    ('A', phi_CA),
    ('D', phi_DE),
    ('E', phi_DE)
]
fg.add_edges_from(edges)
for node, edge in edges:
    fg.add_edge(node, edge, weight=0)

# Draw the factor graph
plt.figure(figsize=(12, 8))
pos = nx.spring_layout(fg)
variable_nodes = [node for node in fg.nodes() if isinstance(node, str)]
factor_nodes = [node for node in fg.nodes() if not isinstance(node, str)]
# nx.draw(fg.to_markov_model())
# plt.draw()
# nx.draw_networkx(fg, pos, with_labels=True, node_color='lightblue', edge_color='gray')
nx.draw_networkx_nodes(fg, pos, nodelist=variable_nodes, node_color='lightblue', node_shape='o', node_size=1000)
nx.draw_networkx_nodes(fg, pos, nodelist=factor_nodes, node_color='lightgreen', node_shape='s', node_size=1500)
nx.draw_networkx_labels(fg, pos, labels={node: node for node in variable_nodes})
nx.draw_networkx_labels(fg, pos, labels={node: 'ϕ' for node in factor_nodes}, font_color='red')
nx.draw_networkx_edges(fg, pos)
plt.title("Factor Graph Representation")
plt.axis('off')
plt.savefig("factor_graph.png")
plt.close()

# Step 3: Perform inference
inference = BeliefPropagation(fg)

# Marginal probability of A
marginal_A = inference.query(variables=['A'], joint=False)
print("Marginal Probability of A:")
print(marginal_A['A'])

# Joint marginal of A and D
marginal_AD = inference.query(variables=['A', 'D'], joint=True)
print("\nJoint Marginal Probability of A and D:")
print(marginal_AD)

# Conditional independence test
cpd_A_given_BC = inference.query(variables=['A'], evidence={'B': 0, 'C': 0})
print("\nConditional Probability of A given B=0 and C=0:")
print(cpd_A_given_BC)

cpd_A_given_BCD = inference.query(variables=['A'], evidence={'B': 0, 'C': 0, 'D': 1})
print("\nConditional Probability of A given B=0, C=0, D=1:")
print(cpd_A_given_BCD)