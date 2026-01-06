from math import isclose
from typing import Tuple, List
from functools import cmp_to_key
from structures import Node

EPSILON = 1e-9

def det(a: Tuple[float, float], b: Tuple[float, float], c: Tuple[float, float]) -> float:
    return (a[0], c[0]) * (b[1] - c[1]) - (a[1] - c[1]) * (b[0] - c[0])

def sort_edges_angularly(nodes: List[Node]):
    def create_comparator(origin: Node, reverse: bool = False):
        origin_coords = origin.coords
        
        def compare(item1, item2):
            node1, _ = item1
            node2, _ = item2
            d = det(origin_coords, node1.coords, node2.coords)
            
            if abs(d) < EPSILON:
                return 0 # coolinear
            
            return -1 if (d > 0) == reverse else 1
        
        return cmp_to_key(compare)
    
    for node in nodes:
        node.outgoing.sort(key=create_comparator(node, reverse=True))
        node.incoming.sort(key=create_comparator(node, reverse=False))
            
        