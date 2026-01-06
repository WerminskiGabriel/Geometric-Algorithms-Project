from typing import List, Tuple
from structures import Node, Separator
from geometry import sort_edges_angularly

def build_graph(vertices_coords: List[Tuple[float, float]], edges_indicies: List[Tuple[int, int]]) -> List[Node]:
    nodes = [Node(x, y) for x, y in vertices_coords]
    
    for start_idx, end_idx in edges_indicies:
        u, v = sorted([start_idx, end_idx])
        
        nodes[u].outgoing.append((nodes[v], 1))
        nodes[v].incoming.append((nodes[u], 1))
        
    return nodes

def balance_weights(nodes: List[Node]):
    sort_edges_angularly(nodes)
    
    # forward sweep
    for node in nodes[1:-1]:
        node.weight_in = sum(w for _, w in node.incoming)
        node.weight_out = len(node.outgoing) # on start every weight is set to 1
        
        if node.weight_in > node.weight_out:
            # ecces weight incoming
            target_node, current_w = node.outgoing.pop()
            diff = node.weight_in - node.weight_out
            new_w = current_w + diff
            
            node.outgoing.append((target_node, new_w))
            
            # updating weight in target node
            idx = next(i for i, (n, _) in enumerate(target_node.incoming) if n == node)
            target_node.incoming[idx] = (node, new_w)
            
    # backward sweep
    for node in reversed(nodes[1:-1]):
        node.weight_in = sum(w for _, w in node.incoming)
        node.weight_out = sum(w for _, w in node.outgoing)        
        
        if node.weight_out > node.weight_in:
            # ecces weigth outgoing
            source_node, current_w = node.incoming.pop()
            diff = node.weight_out - node.weight_in
            new_w = current_w + diff
            
            node.incoming.append((source_node, new_w))
            
            # updating weight in source node
            idx = next(i for i, (n, _) in enumerate(source_node.outgoing) if n == node)
            source_node.outgoing[idx] = (node, new_w)
            
def extract_separators(nodes: List[Node]) -> List[Separator]:
    
    def build_single_chain(start_node: Node, separator: Separator):
        current = start_node
        while current.outgoing:
            separator.add_point(current.coords)
            
            # finding leftmost edge with weight > 0
            edges = current.outgoing
            idx = len(edges) - 1
            while idx >= 0 and edges[idx][1] == 0:
                idx -= 1
                
            if idx < 0: break
            
            next_node, weight = edges[idx]
            
            edges[idx] = (next_node, weight - 1)
            current = next_node
        
        separator.add_point(current.coords)
        
    total_chains = sum(w for _, w in nodes[0].outgoing)
    separators = [Separator() for _ in range(total_chains)]
    
    for sep in separators:
        build_single_chain(nodes[0], sep)
        
        points = sep.points
        for i in range(len(points) - 1):
            sep.add_edge(points[i], points[i+1])
            
    return separators
