eps = 1e-24

class Point:
    def __init__(self, x, y):
        self.x = x
        self.y = y
        self.points_in = []
        self.points_out = []
        self.weights_in = 0
        self.weights_out = 0
        
    def __eq__(self, other):
        return abs(self.x - other.x) <= eps and abs(self.y - other.y) <= eps
    
    def __str__(self):
        return f"({self.x}, {self.y})"
    
    def __repr__(self):
        return f"({self.x}, {self.y})"
    
class Edge:
    def __init__(self, point1, point2):
        self.point1 = point1
        self.point2 = point2
        
        def __str__(self):
            return f"({point1.x}, {point1.y}), ({point2.x}, {point2.y})"
        
        def __repr__(self):
            return f"({point1.x}, {point1.y}), ({point2.x}, {point2.y})"
    
class Separator:
    def __init__(self):
        self.points = []
        self.edges = []
        
        def add_point(self, point):
            self.points.append(point)
            
        def add_edge(self, edge):
            self.edge.append(edge)
            
class BSTNode:
    def __init__(self, separator):
        self.saparator = separator
        self.left = None
        self.right = None