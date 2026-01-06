import matplotlib.pyplot as plt
import math
import os

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

    # Elementy graficzne
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

        # "scroll" - zakończ
        if event.button == 2 :
            plt.close( fig )
            return

        # prawy - Zamknij region
        if event.button == 3:
            if len( current_region ) >= 3 :
                # Rysuj linię domykającą
                p_start = current_region[0]
                p_end = current_region[-1]
                ax.plot( [p_end[0] , p_start[0]] , [p_end[1] , p_start[1]] , 'b-' , linewidth = 2 )

                # Dodaj do listy regionów
                all_regions.append( list( current_region ) )
                for p in current_region :
                    if p not in all_vertices :
                        all_vertices.append( p )

                current_region.clear( )
                temp_artists.clear( )
                fig.canvas.draw( )
            return

        # Lewy-Dodaj punkt
        if event.button == 1 :
            pos = find_nearest_vertex( event.xdata , event.ydata )
            current_region.append( pos )

            # Rysuj punkt i linię
            pt , = ax.plot( pos[0] , pos[1] , 'ro' , zorder = 5 )
            temp_artists.append( pt )

            if len( current_region ) > 1 :
                prev = current_region[-2]
                ln , = ax.plot( [prev[0] , pos[0]] , [prev[1] , pos[1]] , 'k-' , alpha = 0.7 )
                temp_artists.append( ln )

            fig.canvas.draw( )

    # Podpięcie eventów
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



if __name__ == "__main__" :

    all_regions = draw_planar_graph()
    print( "Pol = ", all_regions)
    print()
    v , e = convert_regions_to_graph( all_regions )
    print( "P = ", v )
    print( "E = ", e )
    """
    filename = "graf5"
    export_planar_graph( unify_regions_CCW( all_regions ) , filename )
    print( import_planar_graph( filename ) )
    """