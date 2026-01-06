from typing import List, Tuple, Optional, Union
from structures import BSTNode, Separator
from geometry import det

def build_bst(separators: List[Separator], parent: Optional[BSTNode] = None) -> Optional[BSTNode]:
    
    if not separators:
        return None
    
    mid = len(separators) // 2
    current_sep = separators[mid]
    
    # creating node
    edges_copy = list(current_sep.edges)
    node = BSTNode(edges_copy, current_sep, parent)
    
    # creating children
    node.left = build_bst(separators[:mid], node)
    node.right = build_bst(separators[mid+1:], node)
    
    return node

def query_tree(point: Tuple[float, float], node: BSTNode, last_visited: BSTNode) -> Separator:
    
    if node is None:
        return last_visited.separator if last_visited else None
    
    px, py = point
    
    # binsearch for finding segment in separator alligned with point
    alligned_segment= None
    
    for p1, p2 in node.segments:
        min_y, max_y = min(p1[1], p2[1]), max(p1[1], p2[1])
        
        if min_y <= py <= max_y:
            d = det(p1, p2, point)
            
            if d < 0: # point to the right
                return query_tree(point, node.right, node)
            elif d > 0: # point to the left
                return query_tree(point, node.left, node)
            else: # collinear
                min_x, max_x = min(p1[0], p2[0]), max(p1[0], p2[0])
                if min_x <= px <= max_x:
                    return node.separator
    
    # default case  
    return query_tree(point, node.right, node)

def find_bounding_separators(found_sep: Separator, all_separators: List[Separator], point: Tuple[float, float]):
    
    try:
        idx = all_separators.index(found_sep)
    except ValueError:
        return None
    
    _, py = point
    
    for p1, p2 in found_sep.edges:
        if min(p1[1], p2[1]) <= py <= max(p1[1], p2[1]):
            d = det(p1, p2, point)
            if d < 0:
                return all_separators[idx + 1] if idx + 1 < len(all_separators) else None
            elif d > 0:
                return all_separators[idx - 1] if idx -1 >= 0 else None
            
    return None

def get_region_polygon(sep1: Separator, sep2: Separator, point: Tuple[float, float]) -> List:
    
    if not sep1 or not sep2: return []
    
    edges1 = sep1.edges
    edges2 = sep2.edges
    unique_edges = []
    
    s1_set = set(edges1)
    s2_set = set(edges2)
    
    for e in edges1:
        if e not in s2_set: unique_edges.append(e)
    for e in edges2:
        if e not in s1_set: unique_edges.append(e)
        
    return unique_edges
            
            
        
        