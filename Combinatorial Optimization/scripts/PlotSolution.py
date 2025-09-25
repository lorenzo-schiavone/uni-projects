import networkx as nx
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import json
import os
import argparse

parser = argparse.ArgumentParser(description="instance 50_2")
parser.add_argument("input", type=str, help=".")

args = parser.parse_args()
instance = args.input  

# File paths
hsol_path = f"../solutions/localSearch/sol{instance}.txt" 
esol_path = f"../solutions/exact/sol{instance}.txt" 
gsol_path = f"../solutions/genetic/sol{instance}.txt"
# Load positions
pos_folder = f'../data/pos_dict'
pos_path = os.path.join(pos_folder, f'pos_dict{instance}.json')

width=40
height=20


def sol_to_edge(file_path):
    try:
        with open(file_path, 'r') as f:
            val = float(f.readline().strip()) 
            seq = list(map(int, f.readline().strip().split()))  
        edges = list(zip(seq[:-1], seq[1:])) 
        return val, edges
    except:
        print(f"Error loading {file_path}")
        return None, []

def draw_graph(ax, pos, edges, title):
    G = nx.Graph()
    pos = {int(k): v for k, v in pos.items()}
    G.add_nodes_from(pos.keys()) 
    G.add_edges_from(edges)
    padding = 1
    rect = patches.Rectangle((-padding, -padding), width + 2*padding, height +2*padding, color='#005000', alpha=0.8, zorder=-1)
    ax.add_patch(rect)
    
    nx.draw_networkx_edges(G, pos, ax=ax, edge_color='#DAA520', width=1.5) 
    nx.draw_networkx_nodes(G, pos, ax=ax, node_size=10, node_color='#C0C0C0', edgecolors='#FFFFFF')  

    ax.set_aspect('equal')
    ax.set_title(title)   
    ax.axis('off')     

def draw_nodes(ax, pos, title):
    G = nx.Graph()
    pos = {int(k): v for k, v in pos.items()}
    G.add_nodes_from(pos.keys()) 
    padding = 1
    rect = patches.Rectangle((-padding, -padding), width + 2*padding, height +2*padding, color='#005000', alpha=0.8, zorder=-1)
    ax.add_patch(rect)
    
    nx.draw_networkx_nodes(G, pos, ax=ax, node_size=10, node_color='#C0C0C0', edgecolors='#FFFFFF')  

    ax.set_aspect('equal')
    ax.set_title(title)   
    ax.axis('off')     

try:
    with open(pos_path, 'r') as f:
        pos = json.load(f)
except:
    print(f'{pos_path} not found')
    pos = {}
fig, axes = plt.subplots(2,2, figsize=(8, 8)) 
draw_nodes(axes[0][0], pos, "Nodes") 

hval, hedges = sol_to_edge(hsol_path)
if hval is not None:
    draw_graph(axes[0][1], pos, hedges, f'Local Search: objv= {hval}')

eval, eedges = sol_to_edge(esol_path)
if eval is not None:
    draw_graph(axes[1][0], pos, eedges, f'Exact: objv= {eval}')

gval, gedges = sol_to_edge(gsol_path)
if gval is not None:
    draw_graph(axes[1][1], pos, gedges, f'Genetic: objv= {gval}')

plt.tight_layout()
plt.show()

