#!/usr/bin/env python
# coding: utf-8

# # Projekt: Lokalizacja punktu w przestrzeni dwuwymiarowej – Metoda Separatorów
# **Autorzy:** Gabriel Wermiński, Jacek Łoboda
# **Przedmiot:** Algorytmy Geometryczne
# 
# ## 1. Wstęp Teoretyczny
# Problem lokalizacji punktu polega na znalezieniu regionu $R_i$ w podziale płaszczyzny, który zawiera dany punkt zapytania. W tym projekcie prezentujemy **Metodę separatorów**.
# 
# ### Główne założenia metody:
# 1.  **Monotoniczność:** Podział płaszczyzny składa się z obszarów monotoniczny. Jeśli obszary nie są monotoniczne, wymagana jest regularyzacja (np. triangulacja).
# 2.  **Separatory:** Konstruujemy zbiór łańcuchów monotonicznych, które uporządkowują obszary od lewej do prawej.
# 
# Struktura danych to drzewo binarne, gdzie węzły reprezentują łańcuchy, a liście reprezentują regiony. Zapytanie o punkt polega na porównaniu położenia punktu względem kolejnych łańcuchów w drzewie, co pozwala na osiągnięcie logarytmicznego czasu zapytania.
# 

# In[130]:


import math
import random
from typing import List , Tuple , Optional , Union
from functools import cmp_to_key
import matplotlib as plt
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from bitalg.visualizer.main import Visualizer
from time import time
import os


# In[131]:


eps = 1e-24


class Vertex :
    def __init__( self , x_coord: float , y_coord: float ) :
        self.x = x_coord
        self.y = y_coord
        self.adj_out: List[Tuple['Vertex' , int]] = []
        self.adj_in: List[Tuple['Vertex' , int]] = []
        self.accumulated_weight_in = 0
        self.accumulated_weight_out = 0

    @property
    def coords( self ) -> Tuple[float , float] :
        return self.x , self.y

    def __repr__( self ) :
        return f"V({self.x:.2f}, {self.y:.2f})"


class MonotoneChain :
    def __init__( self ) :
        self.path_vertices: List[Tuple[float , float]] = []
        self.path_segments: List[Tuple[Tuple[float , float] , Tuple[float , float]]] = []
        self.id = -1  # NOWE: ID potrzebne do mapowania regionów

    def add_node( self , coords: Tuple[float , float] ) :
        self.path_vertices.append( coords )

    def add_segment( self , start: Tuple[float , float] , end: Tuple[float , float] ) :
        self.path_segments.append( (start , end) )


class SearchTreeNode :
    def __init__( self , segments: List , chain_ref: MonotoneChain , parent: Optional['SearchTreeNode'] = None ) :
        self.left_child: Optional[SearchTreeNode] = None
        self.right_child: Optional[SearchTreeNode] = None
        self.parent = parent
        self.segments = segments
        self.chain_ref = chain_ref


# In[132]:


def ensure_image_directory( directory ) :
    if not os.path.exists( directory ) :
        os.makedirs( directory )


# ## 2. Przetwarzanie Grafu i Wagi Planarne
# Aby zbudować drzewo łańcuchów, musimy najpierw odpowiednio skierować krawędzie i nadać im wagi. Algorytm wykonuje dwa przejścia:
# 1.  Propagacja wag z dołu do góry.
# 2.  Korekta wag z góry na dół.
# 
# Celem jest ustalenie przepływu tak, aby każdy łańcuch mógł zostać wyodrębniony jako ścieżka monotoniczna od źródła do ujścia grafu.

# In[133]:


def cross_product( o: Tuple[float , float] , a: Tuple[float , float] , b: Tuple[float , float] ) -> float :
    return (a[0] - o[0]) * (b[1] - o[1]) - (a[1] - o[1]) * (b[0] - o[0])


def build_graph( raw_vertices: List[Tuple[float , float]] , raw_edges: List[Tuple[int , int]] ) -> List[Vertex] :
    graph_nodes = [Vertex( x , y ) for x , y in raw_vertices]
    for start_idx , end_idx in raw_edges :
        # Krawędzie zawsze w górę (lub w prawo przy równym Y)
        # Sortujemy tak, aby u był "mniejszy" (niżej/lewo)
        p1 = raw_vertices[start_idx]
        p2 = raw_vertices[end_idx]

        if p1[1] < p2[1] or (math.isclose( p1[1] , p2[1] , abs_tol = eps ) and p1[0] < p2[0]) :
            u , v = start_idx , end_idx
        else :
            u , v = end_idx , start_idx

        node_u = graph_nodes[u]
        node_v = graph_nodes[v]
        node_u.adj_out.append( (node_v , 1) )
        node_v.adj_in.append( (node_u , 1) )
    return graph_nodes


def compute_planar_weights( graph: List[Vertex] ) :
    for node in graph :
        center = node.coords

        def angular_comparator( edge1 , edge2 ) :
            p1 = edge1[0].coords
            p2 = edge2[0].coords
            cp = cross_product( center , p1 , p2 )
            if math.isclose( cp , 0 , abs_tol = eps ) : return 0
            return -1 if cp > 0 else 1

        node.adj_out.sort( key = cmp_to_key( angular_comparator ) )
        node.adj_in.sort( key = cmp_to_key( lambda a , b : -1 * angular_comparator( a , b ) ) )

    # "Forward pass" (dół -> góra)
    # Sortujemy topologicznie (po Y,X ) dla poprawności przepływu
    sorted_nodes = sorted( graph , key = lambda v : (v.y , v.x) )

    for v in sorted_nodes :
        # Pomiń jeśli to źródło
        if not v.adj_in and v != sorted_nodes[0] : pass

        v.accumulated_weight_in = sum( w for _ , w in v.adj_in )
        if not v.adj_in : v.accumulated_weight_in = 0  # Źródło ma 0 wejścia z definicji przepływu, ale generuje łańcuchy

        v.accumulated_weight_out = len( v.adj_out )

        if v.accumulated_weight_in > v.accumulated_weight_out :
            if v.adj_out :
                target_node , current_w = v.adj_out.pop( 0 )
                new_weight = current_w + v.accumulated_weight_in - v.accumulated_weight_out
                v.adj_out.insert( 0 , (target_node , new_weight) )

                for idx , (neighbor , w) in enumerate( target_node.adj_in ) :
                    if neighbor == v :
                        target_node.adj_in[idx] = (v , new_weight)
                        break

    # Backward pass (góra -> dół)
    for v in reversed( sorted_nodes ) :
        v.accumulated_weight_in = sum( w for _ , w in v.adj_in )
        v.accumulated_weight_out = sum( w for _ , w in v.adj_out )

        if v.accumulated_weight_out > v.accumulated_weight_in :
            if v.adj_in :
                source_node , current_w = v.adj_in.pop( 0 )
                new_weight = current_w + v.accumulated_weight_out - v.accumulated_weight_in
                v.adj_in.insert( 0 , (source_node , new_weight) )

                for idx , (neighbor , w) in enumerate( source_node.adj_out ) :
                    if neighbor == v :
                        source_node.adj_out[idx] = (v , new_weight)
                        break


# ## 3. Generowanie Łańcuchów Monotonicznych i Struktura Wyszukiwania
# Algorytm zachłannie buduje łańcuchy, "konsumując" wagi krawędzi. Powstałe łańcuchy są następnie organizowane w zbalansowane drzewo poszukiwań binarnego (BST).
# 
# **Liście drzewa:** Odpowiadają regionom.
# **Węzły wewnętrzne:** Odpowiadają łańcuchom (separatorom) .
# Drzewo to pozwala na szybkie określenie, po której stronie separatora znajduje się punkt.

# In[134]:


def generate_monotone_chains( graph: List[Vertex] ) -> List[MonotoneChain] :
    # Źródło to wierzchołek z najmniejszym Y
    source_node = min( graph , key = lambda v : (v.y , v.x) )
    total_chains = sum( w for _ , w in source_node.adj_out )
    chains = [MonotoneChain( ) for _ in range( total_chains )]

    for idx , chain in enumerate( chains ) :
        chain.id = idx
        current_v = source_node

        while True :
            chain.add_node( current_v.coords )
            if not current_v.adj_out :
                break

            chosen_idx = -1
            # Wybierz skrajnie prawą krawędź z wagą > 0
            for i in range( len( current_v.adj_out ) - 1 , -1 , -1 ) :
                neighbor , weight = current_v.adj_out[i]
                if weight > 0 :
                    chosen_idx = i
                    break

            if chosen_idx == -1 : break

            next_v , w = current_v.adj_out[chosen_idx]
            current_v.adj_out[chosen_idx] = (next_v , w - 1)
            current_v = next_v

        # Generowanie separatorów
        verts = chain.path_vertices
        for k in range( len( verts ) - 1 ) :
            chain.add_segment( verts[k] , verts[k + 1] )

    return chains


def create_search_structure( chains: List[MonotoneChain] , parent: Optional[SearchTreeNode] = None ) -> Optional[
    SearchTreeNode] :
    if not chains : return None
    mid_idx = len( chains ) // 2
    median_chain = chains[mid_idx]

    node = SearchTreeNode( median_chain.path_segments , median_chain , parent )
    node.left_child = create_search_structure( chains[:mid_idx] , node )
    node.right_child = create_search_structure( chains[mid_idx + 1 :] , node )
    return node


# In[135]:


def precompute_regions( chains: List[MonotoneChain] ) -> dict :
    """
    Tworzy słownik regionów
    """
    region_map = { }

    for idx in range( len( chains ) - 1 ) :
        chain_left = chains[idx]
        chain_right = chains[idx + 1]

        verts_a = chain_left.path_vertices
        verts_b = chain_right.path_vertices

        bubbles = []
        i , j = 0 , 0
        last_common_a_idx = 0
        last_common_b_idx = 0

        while i < len( verts_a ) and j < len( verts_b ) :
            va = verts_a[i]
            vb = verts_b[j]

            is_same = math.isclose( va[0] , vb[0] , abs_tol = eps ) and \
                      math.isclose( va[1] , vb[1] , abs_tol = eps )

            if is_same :
                if i > last_common_a_idx or j > last_common_b_idx :
                    path_a = verts_a[last_common_a_idx : i + 1]
                    path_b = verts_b[last_common_b_idx : j + 1]

                    all_verts = path_a + path_b
                    min_y = min( v[1] for v in all_verts )
                    max_y = max( v[1] for v in all_verts )

                    region_edges = []
                    for k in range( len( path_a ) - 1 ) :
                        region_edges.append( (path_a[k] , path_a[k + 1]) )
                    for k in range( len( path_b ) - 1 ) :
                        edge = (path_b[k] , path_b[k + 1])
                        if edge not in region_edges :
                            region_edges.append( edge )

                    bubbles.append( {
                        "y_range" : (min_y , max_y) ,
                        "edges" : region_edges
                    } )

                last_common_a_idx = i
                last_common_b_idx = j
                i += 1
                j += 1
            else :
                if va[1] < vb[1] or (math.isclose( va[1] , vb[1] , abs_tol = eps ) and va[0] < vb[0]) :
                    i += 1
                else :
                    j += 1

        bubbles.sort( key = lambda b : b["y_range"][0] )

        searchable_list = []
        for b in bubbles :
            searchable_list.append( (b["y_range"][1] , b) )

        region_map[(chain_left.id , chain_right.id)] = searchable_list

    return region_map


def query_region_fast( chain_left: MonotoneChain , chain_right: MonotoneChain ,
                       point: Tuple[float , float] ,
                       region_map: dict ) -> List[Tuple] :
    """
    Znajduje region w O(log K) używając binary search.
    """
    if chain_left is None or chain_right is None :
        return []

    key = (chain_left.id , chain_right.id)
    if key not in region_map :
        key = (chain_right.id , chain_left.id)
        if key not in region_map :
            return []

    bubbles_list = region_map[key]
    py = point[1]

    # Binary search
    low = 0
    high = len( bubbles_list )

    while low < high :
        mid = (low + high) // 2
        if bubbles_list[mid][0] < py :
            low = mid + 1
        else :
            high = mid

    insert_idx = low

    # Czy punkt mieści się w zakresie Y
    if insert_idx < len( bubbles_list ) :
        region_data = bubbles_list[insert_idx][1]
        min_y , max_y = region_data["y_range"]

        if min_y <= py <= max_y + eps :
            return region_data["edges"]

    return []


# ## 4. Logika Zapytania (Point Location)
# Dla zadanego punktu $p$ algorytm schodzi w dół drzewa. W każdym węźle sprawdzamy relację punktu względem separatora:
# * Jeśli punkt jest na lewo od separatora -> idziemy do lewego dziecka.
#  * Jeśli punkt jest na prawo od separatora -> idziemy do prawego dziecka.
# 
# Test relacji wykorzystuje iloczyn wektorowy (`cross_product`) oraz wyszukiwanie binarne na segmentach separatora.

# In[136]:


def find_position_relative_to_chain( point: Tuple[float , float] , node: SearchTreeNode ) -> int :
    px , py = point
    segments = node.segments
    target_segment = None

    left_idx , right_idx = 0 , len( segments ) - 1

    while left_idx <= right_idx :
        mid_idx = (left_idx + right_idx) // 2
        p1 , p2 = segments[mid_idx]
        y_min , y_max = p1[1] , p2[1]

        if y_min - eps <= py <= y_max + eps :
            target_segment = (p1 , p2)
            break
        elif py < y_min :
            right_idx = mid_idx - 1
        else :
            left_idx = mid_idx + 1

    if target_segment is None :
        return 1

    cp = cross_product( target_segment[0] , target_segment[1] , point )

    if math.isclose( cp , 0 , abs_tol = eps ) :
        min_x = min( target_segment[0][0] , target_segment[1][0] )
        max_x = max( target_segment[0][0] , target_segment[1][0] )
        if min_x <= px <= max_x :
            return 0
        return -1 if px > target_segment[0][0] else 1

    return -1 if cp < 0 else 1


def query_search_tree( point: Tuple[float , float] , node: Optional[SearchTreeNode] ,
                       closest_left: Optional[MonotoneChain] = None ,
                       closest_right: Optional[MonotoneChain] = None ) -> Union[
    MonotoneChain , Tuple[MonotoneChain , MonotoneChain]] :
    if node is None :
        return (closest_left , closest_right)

    position = find_position_relative_to_chain( point , node )

    if position == 0 :
        return node.chain_ref

    if position < 0 :
        return query_search_tree( point , node.right_child , closest_left = node.chain_ref ,
                                  closest_right = closest_right )
    else :
        return query_search_tree( point , node.left_child , closest_left = closest_left ,
                                  closest_right = node.chain_ref )


# In[137]:


def separators_method_point_location_algorithm_visualiser( raw_vertices , raw_edges , point ) :
    vis = Visualizer( )
    vis.add_point( raw_vertices , color = "black" )

    segments = []
    for u , v in raw_edges :
        p1 = raw_vertices[u]
        p2 = raw_vertices[v]
        segments.append( (p1 , p2) )
    vis.add_line_segment( segments , color = "gray" )

    found_edges = run_point_location( raw_vertices , raw_edges , point )

    vis.add_point( point , color = "green" )
    if found_edges :
        vis.add_line_segment( found_edges , color = "red" , linewidth = 3 )

    return vis , found_edges


def animate_point_location( raw_vertices: List[Tuple[float , float]] ,
                            raw_edges: List[Tuple[int , int]] ,
                            query_point: Tuple[float , float] ) :
    """
    Tworzy animację działania algorytmu lokalizacji punktu.
    """

    graph = build_graph( raw_vertices , raw_edges )
    compute_planar_weights( graph )
    separators = generate_monotone_chains( graph )
    bst_root = create_search_structure( separators )

    vis = Visualizer( )

    vis.add_point( raw_vertices , color = "black" , s = 5 )

    all_segments = []
    for u , v in raw_edges :
        p1 = raw_vertices[u]
        p2 = raw_vertices[v]
        all_segments.append( (p1 , p2) )
    vis.add_line_segment( all_segments , color = "lightgray" , linewidth = 1 )

    vis.add_point( [query_point] , color = "green" , s = 20 )

    current_node = bst_root

    while current_node is not None :

        chain_segments = current_node.segments
        chain_fig = vis.add_line_segment( chain_segments , color = "orange" , linewidth = 2 )

        # binsearch na separatorze
        px , py = query_point
        target_segment = None
        left_idx , right_idx = 0 , len( chain_segments ) - 1

        while left_idx <= right_idx :
            mid_idx = (left_idx + right_idx) // 2
            p1 , p2 = chain_segments[mid_idx]
            y_min , y_max = p1[1] , p2[1]

            if y_min - eps <= py <= y_max + eps :
                target_segment = (p1 , p2)
                break
            elif py < y_min :
                right_idx = mid_idx - 1
            else :
                left_idx = mid_idx + 1

        seg_fig = None
        position = 1

        if target_segment :
            seg_fig = vis.add_line_segment( [target_segment] , color = "blue" , linewidth = 4 )

            cp = cross_product( target_segment[0] , target_segment[1] , query_point )
            if math.isclose( cp , 0 , abs_tol = eps ) :
                position = 0
            else :
                position = -1 if cp < 0 else 1

        if seg_fig :
            vis.remove_figure( seg_fig )

        vis.remove_figure( chain_fig )

        if position == 0 :
            break
        elif position < 0 :
            current_node = current_node.right_child
        else :
            current_node = current_node.left_child

    found_edges = run_point_location( raw_vertices , raw_edges , query_point )
    if found_edges :
        vis.add_line_segment( found_edges , color = "red" , linewidth = 3 )

    return vis


def run_point_location( vertices: List[Tuple[float , float]] , edges: List[Tuple[int , int]] ,
                        query_point: Tuple[float , float] ) :
    graph = build_graph( vertices , edges )
    compute_planar_weights( graph )
    separators = generate_monotone_chains( graph )

    # Budowanie
    region_map = precompute_regions( separators )

    bst_root = create_search_structure( separators )

    result = query_search_tree( query_point , bst_root , closest_left = separators[0] , closest_right = separators[-1] )

    if isinstance( result , MonotoneChain ) :
        segments = result.path_segments
        l , r = 0 , len( segments ) - 1
        while l <= r :
            mid = (l + r) // 2
            seg = segments[mid]
            if seg[0][1] <= query_point[1] <= seg[1][1] :
                return [seg]
            elif query_point[1] < seg[0][1] :
                r = mid - 1
            else :
                l = mid + 1
        return result.path_segments

    # Punkt w regionie, binary search na bąbelkach
    sep_left , sep_right = result
    return query_region_fast( sep_left , sep_right , query_point , region_map )


# ## 5. Obsługa Dowolnych Wielokątów (Budowanie)
# Ponieważ metoda łańcuchów wymaga regionów monotonicznych, dowolne wielokąty (w tym wklęsłe) muszą zostać poddane obróbce. Stosujemy tutaj metodę **Ear Clipping** do triangulacji wielokątów.
# 
# Powstały graf trójkątów może być niespójny (zawierać "wyspy"), dlatego funkcja `patch_graph_connectivity` dodaje wirtualne krawędzie łączące, tworząc spójny graf wymagany przez algorytm.

# In[138]:


def cross_product( o , a , b ) :
    return (a[0] - o[0]) * (b[1] - o[1]) - (a[1] - o[1]) * (b[0] - o[0])


def is_point_in_triangle( p , a , b , c ) :
    cp1 = cross_product( a , b , p )
    cp2 = cross_product( b , c , p )
    cp3 = cross_product( c , a , p )
    has_neg = (cp1 < 0) or (cp2 < 0) or (cp3 < 0)
    has_pos = (cp1 > 0) or (cp2 > 0) or (cp3 > 0)
    return not (has_neg and has_pos)


def get_polygon_signed_area( points ) :
    area = 0.0
    for i in range( len( points ) ) :
        j = (i + 1) % len( points )
        area += points[i][0] * points[j][1]
        area -= points[j][0] * points[i][1]
    return area / 2.0


def triangulate_indices( indices , vertices ) :
    local_indices = indices[:]
    triangles = []
    n = len( local_indices )
    max_iter = n * n
    count = 0

    while len( local_indices ) > 2 :
        if count > max_iter :
            break
        n_curr = len( local_indices )
        ear_found = False
        for i in range( n_curr ) :
            prev_idx = local_indices[(i - 1) % n_curr]
            curr_idx = local_indices[i]
            next_idx = local_indices[(i + 1) % n_curr]

            p_prev = vertices[prev_idx]
            p_curr = vertices[curr_idx]
            p_next = vertices[next_idx]

            tolerance = 1e-8
            if cross_product( p_prev , p_curr , p_next ) > tolerance :
                is_ear = True
                for k in range( n_curr ) :
                    test_idx = local_indices[k]
                    if test_idx in (prev_idx , curr_idx , next_idx) : continue
                    if is_point_in_triangle( vertices[test_idx] , p_prev , p_curr , p_next ) :
                        is_ear = False
                        break
                if is_ear :
                    triangles.append( (prev_idx , curr_idx , next_idx) )
                    local_indices.pop( i )
                    ear_found = True
                    break
        if not ear_found : local_indices.pop( 0 )
        count += 1
    return triangles


def process_polygons_to_mesh( polygons_list ) :
    """
    Tworzy siatkę, sortuje wierzchołki rosnąco po Y, X
    """
    temp_vertices = []
    vertex_map = { }
    all_edges = set( )

    def get_vertex_index( pt ) :
        pt_tuple = (round( pt[0] , 6 ) , round( pt[1] , 6 ))
        if pt_tuple not in vertex_map :
            vertex_map[pt_tuple] = len( temp_vertices )
            temp_vertices.append( pt_tuple )
        return vertex_map[pt_tuple]

    for poly in polygons_list :
        if len( poly ) < 3 : continue
        poly_indices = [get_vertex_index( p ) for p in poly]

        original_coords = [temp_vertices[i] for i in poly_indices]
        if get_polygon_signed_area( original_coords ) < 0 :
            poly_indices.reverse( )

        tris = triangulate_indices( poly_indices , temp_vertices )

        for (i1 , i2 , i3) in tris :
            all_edges.add( tuple( sorted( (i1 , i2) ) ) )
            all_edges.add( tuple( sorted( (i2 , i3) ) ) )
            all_edges.add( tuple( sorted( (i3 , i1) ) ) )

    indexed_vertices = list( enumerate( temp_vertices ) )

    indexed_vertices.sort( key = lambda x : (x[1][1] , x[1][0]) )

    sorted_vertices = []
    old_to_new_map = { }

    for new_idx , (old_idx , pt) in enumerate( indexed_vertices ) :
        sorted_vertices.append( pt )
        old_to_new_map[old_idx] = new_idx

    final_edges = []
    for u , v in all_edges :
        new_u = old_to_new_map[u]
        new_v = old_to_new_map[v]
        final_edges.append( tuple( sorted( (new_u , new_v) ) ) )

    final_edges.sort( )

    return sorted_vertices , final_edges


# ## 6. Proceduralne Grafy

# In[139]:


def generate_tri_grid( width , height ) :
    vertices = []
    edges = []
    for y in range( height + 1 ) :
        for x in range( width + 1 ) :
            vertices.append( (float( x ) , float( y )) )

    def get_idx( x , y ) :
        return y * (width + 1) + x

    for y in range( height ) :
        for x in range( width ) :
            u = get_idx( x , y )
            right = get_idx( x + 1 , y )
            top = get_idx( x , y + 1 )
            top_right = get_idx( x + 1 , y + 1 )
            edges.append( (u , right) )
            edges.append( (u , top) )
            edges.append( (right , top_right) )
            edges.append( (top , top_right) )
            edges.append( (u , top_right) )
    unique_edges = set( )
    for u , v in edges :
        if u < v :
            unique_edges.add( (u , v) )
        else :
            unique_edges.add( (v , u) )
    return vertices , list( unique_edges ) , 0


def generate_quad_grid( width , height , skew = 0.0 ) :
    vertices = []
    edges = []
    for y in range( height + 1 ) :
        shift = y * skew
        for x in range( width + 1 ) :
            vertices.append( (float( x + shift ) , float( y )) )

    def get_idx( x , y ) :
        return y * (width + 1) + x

    for y in range( height + 1 ) :
        for x in range( width ) :
            u = get_idx( x , y )
            v = get_idx( x + 1 , y )
            edges.append( (u , v) )
    for x in range( width + 1 ) :
        for y in range( height ) :
            u = get_idx( x , y )
            v = get_idx( x , y + 1 )
            edges.append( (u , v) )
    unique_edges = set( )
    for u , v in edges :
        if u < v :
            unique_edges.add( (u , v) )
        else :
            unique_edges.add( (v , u) )
    return vertices , list( unique_edges ) , skew


def manual_delaunay_generator( width = 10 , height = 10 , num_points = 15 ) :
    vertices = []
    vertices.extend(
        [(0.0 , 0.0) , (float( width ) , 0.0) , (0.0 , float( height )) , (float( width ) , float( height ))] )

    for _ in range( num_points ) :
        x = random.uniform( 0.1 , width - 0.1 )
        y = random.uniform( 0.1 , height - 0.1 )
        vertices.append( (x , y) )

    n = len( vertices )
    edges = set( )

    for i in range( n ) :
        for j in range( i + 1 , n ) :
            for k in range( j + 1 , n ) :
                x1 , y1 = vertices[i]
                x2 , y2 = vertices[j]
                x3 , y3 = vertices[k]

                D = 2 * (x1 * (y2 - y3) + x2 * (y3 - y1) + x3 * (y1 - y2))

                if abs( D ) < 1e-5 : continue

                Ux = ((x1 ** 2 + y1 ** 2) * (y2 - y3) + (x2 ** 2 + y2 ** 2) * (y3 - y1) + (x3 ** 2 + y3 ** 2) * (
                        y1 - y2)) / D
                Uy = ((x1 ** 2 + y1 ** 2) * (x3 - x2) + (x2 ** 2 + y2 ** 2) * (x1 - x3) + (x3 ** 2 + y3 ** 2) * (
                        x2 - x1)) / D

                r_sq = (Ux - x1) ** 2 + (Uy - y1) ** 2
                center = (Ux , Uy)

                is_valid = True
                for m in range( n ) :
                    if m == i or m == j or m == k : continue
                    dist_sq = (vertices[m][0] - center[0]) ** 2 + (vertices[m][1] - center[1]) ** 2

                    if dist_sq < r_sq - 1e-5 :
                        is_valid = False
                        break

                if is_valid :
                    edges.add( tuple( sorted( (i , j) ) ) )
                    edges.add( tuple( sorted( (j , k) ) ) )
                    edges.add( tuple( sorted( (k , i) ) ) )

    edge_list = list( edges )

    # fallback
    if len( edge_list ) == 0 :
        for i in range( len( vertices ) - 1 ) :
            edge_list.append( (i , i + 1) )
        edge_list.append( (len( vertices ) - 1 , 0) )

    return vertices , edge_list , 0


#get_ipython().run_line_magic('matplotlib', 'tk')

def draw_planar_graph( ) :
    """
    Funkcja interaktywna do rysowania grafu planarnego.
    LPM: Dodaj punkt
    PPM: Zamknij region
    SCROLL: Zakończ i zwróć Graph = [[R1], [R2], ...]
    """

    all_regions = []
    current_region = []
    all_vertices = []
    snap_threshold = 5

    fig , ax = plt.subplots( figsize = (8 , 8) )
    ax.set_title( "LPM: Punkt | PPM: Zamknij R | SCROLL: Zakończ i zwróć dane" )
    ax.set_xlim( 0 , 100 )
    ax.set_ylim( 0 , 100 )
    ax.grid( True , alpha = 0.5 )

    snap_cursor , = ax.plot( [] , [] , 'go' , alpha = 0.5 , markersize = 10 , zorder = 10 )
    temp_artists = []

    def get_distance( p1 , p2 ) :
        return math.sqrt( (p1[0] - p2[0]) ** 2 + (p1[1] - p2[1]) ** 2 )

    def find_nearest_vertex( x , y ) :
        best_v = (x , y)
        min_dist = float( 'inf' )
        # Szukaj w już zapisanych regionach i w obecnym
        for v in all_vertices + current_region :
            dist = get_distance( (x , y) , v )
            if dist < min_dist :
                min_dist = dist
                best_v = v

        if min_dist <= snap_threshold :
            return best_v
        return (x , y)

    def on_move( event ) :
        if event.inaxes == ax :
            pos = find_nearest_vertex( event.xdata , event.ydata )
            snap_cursor.set_data( [pos[0]] , [pos[1]] )
            fig.canvas.draw_idle( )

    def on_click( event ) :
        if event.inaxes != ax : return

        if event.button == 2 :
            plt.close( fig )
            return

        if event.button == 3 :
            if len( current_region ) >= 3 :
                p_start = current_region[0]
                p_end = current_region[-1]
                ax.plot( [p_end[0] , p_start[0]] , [p_end[1] , p_start[1]] , 'b-' , linewidth = 2 )

                all_regions.append( list( current_region ) )
                for p in current_region :
                    if p not in all_vertices :
                        all_vertices.append( p )

                current_region.clear( )
                temp_artists.clear( )
                fig.canvas.draw( )
            return

        if event.button == 1 :
            pos = find_nearest_vertex( event.xdata , event.ydata )
            current_region.append( pos )

            pt , = ax.plot( pos[0] , pos[1] , 'ro' , zorder = 5 )
            temp_artists.append( pt )

            if len( current_region ) > 1 :
                prev = current_region[-2]
                ln , = ax.plot( [prev[0] , pos[0]] , [prev[1] , pos[1]] , 'k-' , alpha = 0.7 )
                temp_artists.append( ln )

            fig.canvas.draw( )

    fig.canvas.mpl_connect( 'button_press_event' , on_click )
    fig.canvas.mpl_connect( 'motion_notify_event' , on_move )

    plt.show( )

    return all_regions


def det( a , b , c ) :
    # >0 CCW | <0 CW
    return (a[0] - c[0]) * (b[1] - c[1]) - (a[1] - c[1]) * (b[0] - c[0])


def is_CCW( graph ) :
    n = len( graph )
    for i in range( 1 , n ) :
        prev = graph[i - 1]
        curr = graph[i]
        next = graph[(i + 1) % n]

        if det( prev , curr , next ) > 0 :
            return False
    return True


def unify_regions_CCW( graph ) :
    new_graph = []
    for R in graph :
        if not is_CCW( R ) :
            new_graph.append( R[: :-1] )
        else :
            new_graph.append( R )
    return new_graph


def export_planar_graph( graph , filename ) :
    filePath = os.path.join( "graphs" , filename )
    with open( filePath , 'w' ) as f :
        for region in graph :
            parts = []
            for x , y in region :
                parts.append( f"{x:.8f},{y:.8f}" )
            f.write( " ".join( parts ) + "\n" )


def import_planar_graph( filename ) :
    filePath = os.path.join( "graphs" , filename )
    graph = []
    with open( filePath , 'r' ) as f :
        for line in f :
            clean_line = line.strip( )
            if not clean_line :
                continue

            region = []
            point_strings = clean_line.split( " " )

            for ps in point_strings :
                if "," in ps :
                    x_str , y_str = ps.split( "," )
                    region.append( (float( x_str ) , float( y_str )) )

            graph.append( region )
    return graph


def convert_regions_to_graph( regions ) :
    all_vertices = []
    vertex_map = { }
    edges = []

    for region in regions :
        for point in region :
            if point not in all_vertices :
                all_vertices.append( point )

    all_vertices.sort( key = lambda x : (x[1] , x[0]) )

    for idx , point in enumerate( all_vertices ) :
        vertex_map[point] = idx

    for region in regions :
        for i in range( len( region ) ) :
            p1 = region[i]
            p2 = region[(i + 1) % len( region )]
            idx1 = vertex_map[p1]
            idx2 = vertex_map[p2]

            edge = tuple( sorted( [idx1 , idx2] ) )
            if edge not in edges :
                edges.append( edge )

    return all_vertices , edges


