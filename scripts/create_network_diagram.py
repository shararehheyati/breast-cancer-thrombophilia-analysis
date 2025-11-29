# create_network_diagram.py
import networkx as nx
import matplotlib.pyplot as plt
import pandas as pd

print("=== Creating Protein-Protein Interaction Network Diagram ===")

def create_ppi_network():
    """
    Create a professional PPI network diagram for thrombophilia genes
    """
    # 1. Create graph
    G = nx.Graph()
    
    # 2. Add nodes (thrombophilia genes)
    thrombophilia_genes = [
        'F2', 'F3', 'F5',           # Coagulation factors
        'SERPINC1', 'PROC', 'PROS1', 'THBD',  # Anticoagulants  
        'FGA', 'FGB', 'FGG',        # Fibrinogen chains
        'PLAT', 'PLAU', 'SERPINE1'  # Fibrinolytic system
    ]
    
    for gene in thrombophilia_genes:
        G.add_node(gene)
        print(f"Added node: {gene}")
    
    # 3. Add edges based on interaction scores (from your heatmap data)
    # Using the interaction scores from your ppi_network_heatmap.png
    interactions = [
        # Coagulation factors module (strong connections)
        ('F2', 'F3', 0.81), ('F2', 'F5', 0.74), ('F3', 'F5', 0.86),
        
        # Anticoagulant module (strong connections)  
        ('SERPINC1', 'PROC', 0.77), ('SERPINC1', 'PROS1', 0.79),
        ('SERPINC1', 'THBD', 0.73), ('PROC', 'PROS1', 0.75),
        ('PROS1', 'THBD', 0.89),
        
        # Fibrinogen complex (very strong connections)
        ('FGA', 'FGB', 0.78), ('FGA', 'FGG', 0.74), ('FGB', 'FGG', 0.86),
        
        # Fibrinolytic system (moderate connections)
        ('PLAU', 'PLAT', 0.36), ('PLAU', 'SERPINE1', 0.37),
        ('PLAT', 'SERPINE1', 0.41),
        
        # Cross-module connections (weaker)
        ('F3', 'PROS1', 0.56), ('F3', 'THBD', 0.49),
        ('PROS1', 'PLAT', 0.28), ('THBD', 'PLAT', 0.35)
    ]
    
    for gene1, gene2, weight in interactions:
        G.add_edge(gene1, gene2, weight=weight)
        print(f"Added edge: {gene1}-{gene2} (weight: {weight})")
    
    # 4. Create the visualization
    plt.figure(figsize=(14, 10))
    
    # Define node positions for better layout
    pos = {
        # Coagulation factors (top-left)
        'F2': (1, 5), 'F3': (2, 5), 'F5': (3, 5),
        
        # Anticoagulants (top-right)  
        'SERPINC1': (6, 6), 'PROC': (7, 6), 'PROS1': (8, 6), 'THBD': (9, 6),
        
        # Fibrinogen complex (bottom-left)
        'FGA': (1, 1), 'FGB': (2, 1), 'FGG': (3, 1),
        
        # Fibrinolytic system (bottom-right)
        'PLAT': (6, 2), 'PLAU': (7, 2), 'SERPINE1': (8, 2)
    }
    
    # Define node colors by functional module
    node_colors = {
        'Coagulation': '#FF6B6B',      # Red
        'Anticoagulant': '#4ECDC4',    # Blue  
        'Fibrinogen': '#45B7D1',       # Light Blue
        'Fibrinolytic': '#FFA07A'      # Orange
    }
    
    # Assign colors to nodes
    node_color_list = []
    for node in G.nodes():
        if node in ['F2', 'F3', 'F5']:
            node_color_list.append(node_colors['Coagulation'])
        elif node in ['SERPINC1', 'PROC', 'PROS1', 'THBD']:
            node_color_list.append(node_colors['Anticoagulant'])
        elif node in ['FGA', 'FGB', 'FGG']:
            node_color_list.append(node_colors['Fibrinogen'])
        else:
            node_color_list.append(node_colors['Fibrinolytic'])
    
    # Draw the network
    # Draw nodes
    nx.draw_networkx_nodes(G, pos, 
                          node_color=node_color_list,
                          node_size=1500,
                          alpha=0.9,
                          edgecolors='black',
                          linewidths=2)
    
    # Draw edges with varying widths based on interaction strength
    edge_widths = [G[u][v]['weight'] * 5 for u, v in G.edges()]
    nx.draw_networkx_edges(G, pos,
                          width=edge_widths,
                          alpha=0.7,
                          edge_color='gray')
    
    # Draw labels
    nx.draw_networkx_labels(G, pos, 
                           font_size=10, 
                           font_weight='bold',
                           font_family='sans-serif')
    
    # Add title and legend
    plt.title('Protein-Protein Interaction Network\nThrombophilia Genes in Breast Cancer', 
              fontsize=16, fontweight='bold', pad=20)
    
    # Create legend
    legend_elements = [
        plt.Line2D([0], [0], marker='o', color='w', 
                  markerfacecolor=node_colors['Coagulation'], markersize=10, label='Coagulation Factors'),
        plt.Line2D([0], [0], marker='o', color='w', 
                  markerfacecolor=node_colors['Anticoagulant'], markersize=10, label='Anticoagulants'),
        plt.Line2D([0], [0], marker='o', color='w', 
                  markerfacecolor=node_colors['Fibrinogen'], markersize=10, label='Fibrinogen Complex'),
        plt.Line2D([0], [0], marker='o', color='w', 
                  markerfacecolor=node_colors['Fibrinolytic'], markersize=10, label='Fibrinolytic System')
    ]
    
    plt.legend(handles=legend_elements, loc='upper left', bbox_to_anchor=(0, 1))
    
    # Remove axes
    plt.axis('off')
    
    # Adjust layout and save
    plt.tight_layout()
    plt.savefig('results/figures/ppi_network_diagram.png', 
                dpi=300, bbox_inches='tight', facecolor='white')
    plt.show()
    
    print(f"\n=== NETWORK STATISTICS ===")
    print(f"Nodes: {G.number_of_nodes()}")
    print(f"Edges: {G.number_of_edges()}")
    print(f"Average degree: {sum(dict(G.degree()).values()) / G.number_of_nodes():.2f}")
    
    # Calculate modularity
    from networkx.algorithms.community import greedy_modularity_communities
    communities = list(greedy_modularity_communities(G))
    print(f"Detected communities: {len(communities)}")
    
    return G

# Run the function
if __name__ == "__main__":
    import os
    os.makedirs('results/figures', exist_ok=True)
    
    network = create_ppi_network()
    print(f"\n✅ Network diagram saved as: results/figures/ppi_network_diagram.png")