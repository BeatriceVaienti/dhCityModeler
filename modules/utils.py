import math
import copy
import cadquery as cq
from queue import PriorityQueue
from math import inf
import time

def remove_nearly_duplicate_points(points, tolerance=0.1):
    """
    Removes nearly duplicate points from a list of points

    Args:  
        - points (list): list of points
        - tolerance (float): tolerance to consider two points as equal
    Returns:
        - new_points (list): list of points without nearly duplicates
    """
    new_points = []
    check_points = set()
    places = round(-math.log10(tolerance))
    for p in points:
        round_p = tuple(round(coord, places) for coord in p)
        if round_p not in check_points:
            check_points.add(round_p)
            new_points.append(p)
    return new_points

def remove_duplicate_points(points):
    """
    Removes duplicate points from a list of points
    Args:
        - points (list): list of points
    Returns:
        - new_points (list): list of points without duplicates
    """
    points_copy=copy.deepcopy(points)
    new_points = []
    for p in points_copy:
        if p not in new_points:
            new_points.append(p)
    return new_points

def remove_points_in_parallel(points):
    """
    Removes the midpoint in three aligned points (angle=pi)
    or the midpoint in two coincident lines (angle=0).

    Args:
        - points (list): list of points

    Returns:
        - new_points (list): list of points with parallel points removed
    """
    new_points = []
    for i in range(len(points)):
        prev_point = points[i-1]
        curr_point = points[i]
        next_point = points[(i+1)%len(points)]
        AB = [curr_point[j] - prev_point[j] for j in range(3)]
        BC = [next_point[j] - curr_point[j] for j in range(3)]
        if math.isclose(AB[0]*BC[1], AB[1]*BC[0]) and math.isclose(AB[1]*BC[2], AB[2]*BC[1]) and math.isclose(AB[0]*BC[2], AB[2]*BC[0]):
            continue
        new_points.append(curr_point)
    return new_points

def remove_nearly_parallel(points, tolerance=0.001):
    """
    given a list of points it removes the mid point in three almost aligned points (angle circa pi)
    or the mid point in two almost coincident lines (angle=0)
    """
    points_copy=copy.deepcopy(points)
    for i in range(0, len(points_copy), 1):
        points3=(points_copy[i:i+3] if i+3<= len(points_copy) else points_copy[i:]+points_copy[:i+3-len(points_copy)])
        An=points3[0]
        Bn=points3[1]
        Cn=points3[2]
        BA=cq.Vector(Bn[0]-An[0], Bn[1]-An[1], Bn[2]-An[2])
        BC=cq.Vector(Bn[0]-Cn[0], Bn[1]-Cn[1], Bn[2]-Cn[2])
        angle=angle_between(BA,BC)
        #take into account the case of the last point
        if i==len(points_copy)-2:
            if (angle > math.pi-tolerance and angle < math.pi+tolerance) or (angle < tolerance and angle > -tolerance):
                points_copy.pop(0)
        elif (i<len(points_copy)-2):
            if (angle > math.pi-tolerance and angle < math.pi+tolerance) or (angle < tolerance and angle > -tolerance):
                points_copy.pop(i+1)
    return points_copy

def fix_orientation(points):
    """
    Checks if the orientation of the points is clockwise or counter clockwise and returns the points in the correct order to perform the roof generation
    """ 
    points_copy=copy.deepcopy(points)
    wire=cq.Workplane("XY").polyline(points_copy).close().wire().val()
    z_dir = cq.Face.makeFromWires(wire).normalAt().z
    if z_dir<0:
        points_copy=points_copy[::-1]
    return points_copy

def fix_points_get_z(points):
    """
    the function takes a list of points and returns a list of points with z=0 and the z value of the first point
    moreover the points are cleaned from duplicates, nearly duplicates, parallel points and the orientation is fixed

    Args:
        - points (list): list of points

    Returns:
        - points_copy (list): list of points with z=0
        - z_base (float): z value of the first point
    """
    points_copy=copy.deepcopy(points)
    try:
        if type(points_copy[0]) is not list:
            points_copy=[list(point) for point in points_copy]
    except:
        print(points_copy)
        
    if len(points_copy[0])<3:
        for point in points_copy:
            point.append(0)   
    z_base=points_copy[0][2]
    for point in points_copy:
        point[2]=0

    #cleaning the list of points
    points_copy=remove_duplicate_points(points_copy)
    points_copy=remove_nearly_duplicate_points(points_copy)
    points_copy=remove_points_in_parallel(points_copy)
    points_copy=fix_orientation(points_copy)
    return points_copy, z_base


def angle_between(v1, v2):
    """
    Returns the angle in radians between vectors 'v1' and 'v2':
    
    Args:
        - v1 (Vector): first vector
        - v2 (Vector): second vector
    
    Returns:
        - angle (float): angle between v1 and v2 in radians
    """
    v1_u = v1.normalized()
    v2_u = v2.normalized()
    v1_v2=  v1_u.x*v2_u.x + v1_u.y*v2_u.y + v1_u.z*v2_u.z 
    if v1_v2>1:
        v1_v2=1
    return math.acos(v1_v2)


def cross_sign(A, B):
    """
    Calculates the sign of the cross product between two vectors A and B:
    
    Args:
        - A (Vector): first vector
        - B (Vector): second vector

    Returns:
        - (bool): True if cross is positive, False if negative or zero
    """
    x1=A.x
    y1=A.y
    x2=B.x
    y2=B.y
    return x1 * y2 > x2 * y1


def perpendicular_vector(B,A):
    """ 
    Given two points B and A and the vector from A to B (AB)
    returns the unit vector perpendicular to it on the same xy plane

    Prodotto vettoriale tra AB e Z:
	#i               j           k        i               j       
	#AB[0]    AB[1]     0       AB[0]    AB[1]   
	#0             0          1         0              0
    """
    #vettore da A a B:
    AB=[B[0]-A[0], B[1]-A[1], 0]
    
    P=[AB[1],-AB[0], 0]
    LengthP = math.sqrt(P[0]*P[0] + P[1]*P[1]+ P[2]*P[2])
    return cq.Vector(P[0]/LengthP  ,P[1]/LengthP  ,P[2]/LengthP)

def get_sphere_radius(base_center, base_radius, p_v):
        """
        given a base circle and a vertical point between the 0 and the radius of the base circle
        returns the radius of the sphere passing through the base circle and through the vertical point

        Args:
            - base_center (tuple): center of the base circle
            - base_radius (float): radius of the base circle
            - p_v (tuple): vertical point positioned on xy of the center and with a z that can be <= to the radius of the base circle
        Returns:
            - sphere_radius (float): radius of the sphere passing through the base circle and through the vertical point
        """

        base_circle = cq.Workplane("XY", base_center).circle(base_radius)
        ## given a base circle we want to find the sphere passing through the base circle and through the vertical point between the 0 and the radius of the base circle
        #the vertical point has the xy coordinates of the center and the z = radius_base*p_v
        vertical_point=[base_center[0],base_center[1],base_radius*p_v]
        vertical_point_tuple=(base_center[0],base_center[1],base_radius*p_v)
        B=vertical_point_tuple
    
        #the first point is the first point of the base circle
        first_point=base_circle.edges().vertices().val().Center()
        first_point_tuple=(first_point.x, first_point.y, first_point.z)
        A=first_point_tuple
    
        last_point=base_circle.edges().val().positionAt(0.5)
        last_point_tuple=(last_point.x,last_point.y,last_point.z)
        C=last_point_tuple
    
        #calculate BC AB and AC (vectors)
        BC=cq.Vector((C[0]-B[0],C[1]-B[1],C[2]-B[2]))
        BA=cq.Vector((A[0]-B[0],A[1]-B[1],A[2]-B[2]))
        #calculate angle between BC and BA
        angle=angle_between(BC,BA)
        sphere_radius=base_radius/math.sin(angle)
        return sphere_radius

"""
the following function are used in order to find the center and radius of the maximum center inscribed in a concave polygon ("poles of inaccessibility" problem)
the scripts were taken from the polylabel python library: https://github.com/Twista/python-polylabel
which is improving the original algorithm by Garcia-Castellano and Umberto Lombardo (2008)
"""

def _point_to_polygon_distance(x, y, polygon):
    inside = False
    inf = float("inf")
    min_dist_sq = inf

    for ring in polygon:
        b = ring[-1]
        for a in ring:

            if ((a[1] > y) != (b[1] > y) and
                    (x < (b[0] - a[0]) * (y - a[1]) / (b[1] - a[1]) + a[0])):
                inside = not inside

            min_dist_sq = min(min_dist_sq, _get_seg_dist_sq(x, y, a, b))
            b = a

    result = math.sqrt(min_dist_sq)
    if not inside:
        return -result
    return result


def _get_seg_dist_sq(px, py, a, b):
    x = a[0]
    y = a[1]
    dx = b[0] - x
    dy = b[1] - y

    if dx != 0 or dy != 0:
        t = ((px - x) * dx + (py - y) * dy) / (dx * dx + dy * dy)

        if t > 1:
            x = b[0]
            y = b[1]

        elif t > 0:
            x += dx * t
            y += dy * t

    dx = px - x
    dy = py - y

    return dx * dx + dy * dy


class Cell(object):
    def __init__(self, x, y, h, polygon):
        self.h = h
        self.y = y
        self.x = x
        self.d = _point_to_polygon_distance(x, y, polygon)
        self.max = self.d + self.h * math.sqrt(2)
    def __lt__(self, other):
            return self.max < other.max

    def __lte__(self, other):
        return self.max <= other.max

    def __gt__(self, other):
        return self.max > other.max

    def __gte__(self, other):
        return self.max >= other.max

    def __eq__(self, other):
        return self.max == other.max


def _get_centroid_cell(polygon):
    area = 0
    x = 0
    y = 0
    points = polygon[0]
    b = points[-1]  # prev
    for a in points:
        f = a[0] * b[1] - b[0] * a[1]
        x += (a[0] + b[0]) * f
        y += (a[1] + b[1]) * f
        area += f * 3
        b = a
    if area == 0:
        return Cell(points[0][0], points[0][1], 0, polygon)
    return Cell(x / area, y / area, 0, polygon)

    pass


def polylabel(polygon, precision=0.1, debug=False, with_distance=True):
    # find bounding box
    first_item = polygon[0][0]
    min_x = first_item[0]
    min_y = first_item[1]
    max_x = first_item[0]
    max_y = first_item[1]
    for p in polygon[0]:
        if p[0] < min_x:
            min_x = p[0]
        if p[1] < min_y:
            min_y = p[1]
        if p[0] > max_x:
            max_x = p[0]
        if p[1] > max_y:
            max_y = p[1]

    width = max_x - min_x
    height = max_y - min_y
    cell_size = min(width, height)
    h = cell_size / 2.0

    cell_queue = PriorityQueue()

    if cell_size == 0:
        if with_distance:
            return [min_x, min_y], None
        else:
            return [min_x, min_y]

    # cover polygon with initial cells
    x = min_x
    while x < max_x:
        y = min_y
        while y < max_y:
            c = Cell(x + h, y + h, h, polygon)
            y += cell_size
            cell_queue.put((-c.max, time.time(), c))
        x += cell_size

    best_cell = _get_centroid_cell(polygon)

    bbox_cell = Cell(min_x + width / 2, min_y + height / 2, 0, polygon)
    if bbox_cell.d > best_cell.d:
        best_cell = bbox_cell

    num_of_probes = cell_queue.qsize()
    while not cell_queue.empty():
        _, __, cell = cell_queue.get()

        if cell.d > best_cell.d:
            best_cell = cell

            if debug:
                print('found best {} after {} probes'.format(
                    round(1e4 * cell.d) / 1e4, num_of_probes))

        if cell.max - best_cell.d <= precision:
            continue

        h = cell.h / 2
        c = Cell(cell.x - h, cell.y - h, h, polygon)
        cell_queue.put((-c.max, time.time(), c))
        c = Cell(cell.x + h, cell.y - h, h, polygon)
        cell_queue.put((-c.max, time.time(), c))
        c = Cell(cell.x - h, cell.y + h, h, polygon)
        cell_queue.put((-c.max, time.time(), c))
        c = Cell(cell.x + h, cell.y + h, h, polygon)
        cell_queue.put((-c.max, time.time(), c))
        num_of_probes += 4

    if debug:
        print('num probes: {}'.format(num_of_probes))
        print('best distance: {}'.format(best_cell.d))
    if with_distance:
        return [best_cell.x, best_cell.y], best_cell.d
    else:
        return [best_cell.x, best_cell.y]

