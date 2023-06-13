"""
This scripts manage the download and processing of the terrain data.
in order for it to work the following dependencies are necessary:
- numpy
- pyproj
- pydelatin # conda install -c conda-forge pydelatin
- scipy # conda install -c conda-forge scipy
- cjio # optional: needed only for cleaning unused vertices
"""
import requests
import json
import numpy as np 
from pyproj import Transformer
from pydelatin import Delatin 
import modules.encoder as encoder
import copy
from pyproj import CRS


def rescale_positions(vertices: np.ndarray, bounds, flip_y: bool = True):
    """
    This function was modified from PyDelatin in order to use float64 instead of float32
    and thus avoid truncation

    Rescale positions to bounding box
    Args:
        - vertices: vertices output from Delatin
        - bounds: linearly rescale position values to this extent, expected to
          be [minx, miny, maxx, maxy].
        - flip_y: (bool) Flip y coordinates. Can be useful since images'
          coordinate origin is in the top left.
    Returns:
        (np.ndarray): ndarray of shape (-1, 3) with positions rescaled. Each row
        represents a single 3D point.
    """
    out = np.zeros(vertices.shape, dtype=np.float64)  ## changed from float32 to float64 to avoid truncation
    
    tile_size = vertices[:, :2].max()
    minx, miny, maxx, maxy = bounds
    x_scale = (maxx - minx) / tile_size
    y_scale = (maxy - miny) / tile_size

    if flip_y:
        scalar = np.array([x_scale, -y_scale])
        offset = np.array([minx, maxy])
    else:
        scalar = np.array([x_scale, y_scale])
        offset = np.array([minx, miny])

    # Rescale x, y positions
    out[:, :2] = vertices[:, :2] * scalar + offset
    out[:, 2] = vertices[:, 2]

    return out

def get_cityjson_bounds(cityjson2D):
    """
    Gets the 2D CityJSON bounding box as xmin, ymin, xmax, ymax 
    all in the original CRS

    Args:
        - cityjson2D: the 2D CityJSON file 
    Returns:
        - xmin_original
        - ymin_original
        - xmax_original
        - ymax_original
        (the two corners of the 2D bounding box of the JSON file)
    """
    GEOG_EXTENT= cityjson2D['metadata']['geographicalExtent']
    xmin_original, ymin_original = GEOG_EXTENT[0], GEOG_EXTENT[1] # xmin and ymin in the original CRS
    xmax_original, ymax_original = GEOG_EXTENT[3], GEOG_EXTENT[4] # xmax and ymax in the original CRS
    return xmin_original, ymin_original, xmax_original, ymax_original

def get_elevations(cityjson2D, service='opentopodata', apikey=None, grid_min_distance = 30): 
    """
    Get the elevation data corresponding to a certain bounding box (xmin and ymin)
    and with a given grid density (regardless of the shape of the bounding rectangle, 
    the grid is divided in an equal number of rows and columns)
    Args:
        - cityjson2D: the 2D CityJSON file 
        - apikey: the GOOGLE API key
        - grid_min_distance: desired minimum distance between two points of the grid in meters
        - service: the elevation service to use. Can be 'google' or 'opentopodata'
    Returns:
        - elevations: a list containing the elevation of the points in the grid
        - n_points: the number of points per side
    """
    xmin_original, ymin_original, xmax_original, ymax_original = get_cityjson_bounds(cityjson2D)
    min_side = min(xmax_original-xmin_original, ymax_original-ymin_original) # we get the min side of the bounding box
    n_points = int(min_side/grid_min_distance) # we compute the number of points we want to have in the grid
    in_epsg = cityjson2D['metadata']['referenceSystem'].split('/')[-1] # Define the EPSG codes for the original CRS and the target CRS
    crs_in = CRS.from_user_input(in_epsg)
    crs_out = CRS.from_user_input("4326")
    transformer = Transformer.from_crs(crs_in, crs_out, always_xy=False)  # pyproj transformer
    lats = np.linspace(ymin_original, ymax_original, n_points)
    lons = np.linspace(xmin_original, xmax_original, n_points)
    xx, yy = np.meshgrid(lons, lats)
    points = np.column_stack((xx.ravel(), yy.ravel()))
    elevations = []
    query = ""
    counter = 0
    N = 100
    if service == 'opentopodata' or apikey is None:
        for point in points:
            #we want to split the query in chunks of 100 points max
            lat, lon = transformer.transform(point[0], point[1])
            query += f"{lat},{lon}|"
            service_url = f'https://api.opentopodata.org/v1/srtm30m?locations='+query
            counter += 1
            if counter == N:
                response = requests.get(service_url)
                data = json.loads(response.text)
                for i in range(counter):
                    elevation = data["results"][i]["elevation"]
                    elevations.append(elevation)
                counter = 0
                query = ""
        if counter > 0:
            response = requests.get(service_url)
            data = json.loads(response.text)
            for i in range(counter): 
                elevation = data["results"][i]["elevation"]
                elevations.append(elevation)
        if True and len(elevations)%100 == 0:
            print("len_elevations", len(elevations))
    elif service == 'google' and apikey is not None:
        for point in points:
            lat, lon = transformer.transform(point[0], point[1])
            service_url = f"https://maps.googleapis.com/maps/api/elevation/json?locations={lat},{lon}&key={apikey}"
            response = requests.get(service_url)
            data = json.loads(response.text)
            elevation = data["results"][0]["elevation"]
            elevations.append(elevation)
            if True and len(elevations)%100 == 0:
                print(len(elevations))
    return elevations, n_points

def create_TIN(elevations, n_points, tin_max_error = 1):
    """
    Creates a TIN (Triangular Irregular Network) corresponding to the elevation data
    using pydelatin (A Python wrapper of hmm (of which Delatin is a port) for fast terrain mesh generation.)
    [https://github.com/kylebarron/pydelatin]

    Args:
        - elevations: a list containing the elevation of the points in the grid
        - n_points: the number of points per side
        - tin_max_error: (float) the maximum error of the mesh.
    Returns:
        - vertices: (ndarray of shape (-1, 3)): the interleaved 3D coordinates of each vertex, 
                    e.g. [[x0, y0, z0], [x1, y1, z1], ...].
        - triangles:  (ndarray of shape (-1, 3)): represents indices within the vertices array. 
                        So [0, 1, 3, ...] would use the first, second, and fourth vertices within 
                        the vertices array as a single triangle. 
    """
    terrain = np.array(elevations).reshape(n_points, n_points).astype(np.float32)
    tin = Delatin(terrain, max_error=tin_max_error)
    vertices, triangles = tin.vertices, tin.triangles
    return vertices, triangles

def get_geographical_extent(vertices):
    """
    Gets the geographical extents of a set of vertices 

    Args: 
        - vertices: list of vertices
    Returns:
        - [minx, miny, minz, maxx, maxy, maxz]: list of the two corners of the bounding box. 
    """
    minx = min([x[0] for x in vertices])
    miny = min([x[1] for x in vertices])
    minz = min([x[2] for x in vertices])
    maxx = max([x[0] for x in vertices])
    maxy = max([x[1] for x in vertices])
    maxz = max([x[2] for x in vertices])
    return [minx, miny, minz, maxx, maxy, maxz]

def find_triangle_in_tin(vertices, triangles, x, y):
    vertices_copy = copy.deepcopy(vertices)
    triangles_copy = copy.deepcopy(triangles)

    for i, triangle in enumerate(triangles_copy):
        v1, v2, v3 = triangle
        v1_x, v1_y, _ = vertices_copy[v1]
        v2_x, v2_y, _ = vertices_copy[v2]
        v3_x, v3_y, _ = vertices_copy[v3]

        # Check if the point (x, y) is inside the triangle using barycentric coordinates
        denom = (v2_y - v3_y) * (v1_x - v3_x) + (v3_x - v2_x) * (v1_y - v3_y)

        alpha = ((v2_y - v3_y) * (x - v3_x) + (v3_x - v2_x) * (y - v3_y)) / denom
        beta = ((v3_y - v1_y) * (x - v3_x) + (v1_x - v3_x) * (y - v3_y)) / denom
        gamma = 1 - alpha - beta

        if 0 <= alpha <= 1 and 0 <= beta <= 1 and 0 <= gamma <= 1:
            return i

        # Check with a different vertex order
        denom = (v3_y - v1_y) * (v2_x - v1_x) + (v1_x - v3_x) * (v2_y - v1_y)

        alpha = ((v3_y - v1_y) * (x - v1_x) + (v1_x - v3_x) * (y - v1_y)) / denom
        beta = ((v1_y - v2_y) * (x - v1_x) + (v2_x - v1_x) * (y - v1_y)) / denom
        gamma = 1 - alpha - beta

        if 0 <= alpha <= 1 and 0 <= beta <= 1 and 0 <= gamma <= 1:
            return i

        # Check with a different vertex order
        denom = (v1_y - v2_y) * (v3_x - v2_x) + (v2_x - v1_x) * (v3_y - v2_y)

        alpha = ((v1_y - v2_y) * (x - v2_x) + (v2_x - v1_x) * (y - v2_y)) / denom
        beta = ((v2_y - v3_y) * (x - v2_x) + (v3_x - v2_x) * (y - v2_y)) / denom
        gamma = 1 - alpha - beta

        if 0 <= alpha <= 1 and 0 <= beta <= 1 and 0 <= gamma <= 1:
            return i
    return None   

def interpolate_z_in_triangle(vertices, triangle, x, y):
    """
    Interpolates the z-coordinate of a point (x, y) within a triangle in a TIN mesh.

    Args:
        - vertices: List of vertices in the TIN mesh.
        - triangle: Triangle in the TIN mesh, specified as a tuple of vertex indices.
        - x: x-coordinate of the point.
        - y: y-coordinate of the point.

    Returns:
        - Interpolated z-coordinate of the point (x, y).
    """

    v1, v2, v3 = triangle
    x1, y1, z1 = vertices[v1]
    x2, y2, z2 = vertices[v2]
    x3, y3, z3 = vertices[v3]

    # Calculate barycentric coordinates of the point (x, y)
    denom = (y2 - y3) * (x1 - x3) + (x3 - x2) * (y1 - y3)
    alpha = ((y2 - y3) * (x - x3) + (x3 - x2) * (y - y3)) / denom
    beta = ((y3 - y1) * (x - x3) + (x1 - x3) * (y - y3)) / denom
    gamma = 1 - alpha - beta

    # Interpolate z-coordinate using barycentric coordinates
    z = alpha * z1 + beta * z2 + gamma * z3

    return z

def get_z_value(x, y, vertices, triangles):
    """
    Given a point (x,y) and the vertices of the terrain, 
    evaluates the height of the point projected on the terrain mesh

    Args:
        - x
        - y
        - vertices
    Returns:
        - z
    """
#get the list of the z values of the vertices
    z_values = []
    for vertex in vertices: #VERTICES OF THE MESH 
        z_values.append(vertex[2])
    # Find the index of the triangle that contains the (x, y) point
    i = find_triangle_in_tin(vertices, triangles, x, y)
    if i is None:
        raise ValueError("Point outside of triangulation")
    else:
        # Get the vertices of the triangle
        #triangle_vertices = vertices_array[triangles[i]]
        triangle = triangles[i]
        z = interpolate_z_in_triangle(vertices, triangle, x, y)
        return z

def check_terrain_presence(cityjson2d):
    """
    Checks whether the Cityjson2D already contains a terrain file 
    Args:
        - cityjson2d
    Returns:
        - terrain_presence: (bool) True if already present
    """
    terrain_presence=False
    for object in cityjson2d['CityObjects']:
        if cityjson2d['CityObjects'][object]['type'] == 'TINRelief':
            terrain_presence=True
            break
    return terrain_presence

def get_z_translation(footprint_vertices, terrain_vertices, triangles):
    """
    Given a building and a terrain, evaluates the z translation needed to put the building on the terrain.
    The z value corresponding to all vertices of the footprint is evaluated and the minimum is returned.

    Args:
        - cityjson2d_terrain: the cityjson2d containing the terrain and the non-translated buildings
        - id_building: the id of the building to evaluate
        - terrain_vertices: the vertices of the terrain 
    Returns:
        - z_final: the z translation needed to put the building on the terrain
    """
    #we check if the LOD0 of the building is already at the z level of the terrain

    z_value = footprint_vertices[0][2]
    z_values = []
    z_final = None
    if z_value == 0:
        at_terrain_level = False
        for vertex in footprint_vertices:  #LOOKING FOR XY VALUES OF THE VERTICES OF THE BUILDING
            x, y = vertex[0], vertex[1]
            z_new = get_z_value(x, y, terrain_vertices, triangles)
            z_values.append(z_new)
        
        min_translation = min(z_values)
        max_translation = max(z_values)
        difference = max_translation - min_translation
        #we translate the building of z_final
        return min_translation, difference

def get_terrain_from_cityjson_bounds(cityjson2D, service, apikey, grid_min_distance = 20, tin_max_error = 1):
    """
    this function could be furtherly divided!
    Args:
        - cityjson2D: the cityjson2D file+
        - service: the service to use to get the elevation data (default: 'opentopodata')
        - apikey: your google api key if you want to use the 'google' service
        - grid_min_distance: minimum distance between two points of the grid (default: 20)
        - tin_max_error: maximum error allowed for the TIN triangulation (default: 1)
    Returns:
        - vertices: the vertices of the terrain
        - triangles: the triangles of the terrain
    """
    elevations, n_points = get_elevations(cityjson2D, service, apikey, grid_min_distance) #we make the request to the GOOGLE api to get the elevation data
    vertices, triangles = create_TIN(elevations, n_points, tin_max_error) #We generate a TIN out of the elevation grid

    return vertices,triangles

def insert_terrain_in_cityjson(cityjson2D, service = 'opentopodata', apikey = None, grid_min_distance = 30, tin_max_error = 1):
    """
    Inserts the terrain in the cityjson2D file and translates the buildings accordingly

    Args:
        - cityjson2D: the cityjson2D file
        - service: the service to use to get the elevation data (default: 'opentopodata')
        - apikey: your google api key if you want to use the 'google' service
        - grid_min_distance: minimum distance between points of the grid (default: 20m)
        - tin_max_error: maximum error allowed for the TIN
    Returns:
        - cityjson2D: the cityjson2D file with the terrain inserted
    """
    terrain_presence = check_terrain_presence(cityjson2D)
    if terrain_presence == True:
        #find the terrain and delete it
        for object in cityjson2D['CityObjects']:
            if cityjson2D['CityObjects'][object]['type'] == 'TINRelief':
                id_terrain = object
                break
        del cityjson2D['CityObjects'][id_terrain]
        #put the buildings back at level 0
        for object in cityjson2D['CityObjects']:
            if cityjson2D['CityObjects'][object]['type'] == 'Building':
                for vertex in cityjson2D['vertices']:
                    vertex[2] = 0

    vertices, triangles = get_terrain_from_cityjson_bounds(cityjson2D, service, apikey, grid_min_distance, tin_max_error)
    old_cityjson2D = copy.deepcopy(cityjson2D)
    GEOG_EXTENT= cityjson2D['metadata']['geographicalExtent']
    SCALE_FACTOR = cityjson2D['transform']['scale']
    TRANSL_FACTOR =  cityjson2D['transform']['translate']
    bounds= [GEOG_EXTENT[0], GEOG_EXTENT[1], GEOG_EXTENT[3], GEOG_EXTENT[4]] #bounds of the cityjson file

    rescaled_vertices = rescale_positions(vertices, bounds, flip_y=True).tolist()
    converted_all_vertices = ((np.array(rescaled_vertices)-np.array(TRANSL_FACTOR))/np.array(SCALE_FACTOR)).astype(int)
    converted_all_vertices = [list([int(x) for x in i]) for i in converted_all_vertices]

    triangles_list = triangles.tolist()
    new_triangles_list = []
    for triangle in triangles_list:
        #flip the triangle (it's a list)
        new_triangles_list.append([triangle]) #we need to add a list inside the list to make it work with the cityjson format
    new_triangles_list #indeces of the vertices encoded right for cityjson but we need to remap them according to the new index numbers

    old_vertices = cityjson2D['vertices']  #we need to add the vertices to the cityjson file when the vertices are not already present
    new_indeces_list = []
    for vertex in converted_all_vertices:
        if vertex not in old_vertices:
            old_vertices.append(vertex)
            new_indeces_list.append(len(old_vertices)-1)
        else:
            new_indeces_list.append(old_vertices.index(vertex))
    
    #now we need to update the triangles list with the new indeces
    for i in range(len(new_triangles_list)):
        new_triangles_list[i] = [[new_indeces_list[x] for x in new_triangles_list[i][0]]]
    # now we need to reverse each triangle in new_triangles_list
    for i in range(len(new_triangles_list)):
        new_triangles_list[i][0].reverse()
    terrain_cityobject = create_terrain_cityobject(rescaled_vertices, new_triangles_list)

    cityjson2D['CityObjects'].update(terrain_cityobject) #add the new object to the cityjson file
    cityjson2D = translate_lod0_geometries(cityjson2D, old_cityjson2D, converted_all_vertices, triangles)
    return cityjson2D

def create_terrain_cityobject(rescaled_vertices, triangles_list):
    """
    Creates the terrain TINRelief CityObject from the rescaled vertices and the list of triangles
    Args:
        - rescaled_vertices
        - triangles_list
    Returns:
        - terrain_cityobject
    """
    geographical_extent = get_geographical_extent(rescaled_vertices) # we get the geographical extent of the new terrain
 
    terrain_cityobject= {
        "terrain": {
            "type": "TINRelief", 
            "geographicalExtent": geographical_extent,
            "geometry": [{
                "type": "CompositeSurface",
                "lod": "1",
                "boundaries": triangles_list
            }]    
        }
    }
    return terrain_cityobject

def translate_lod0_geometries(cityjson2D, old_cityjson2D, vertices,triangles):
    """
    Iterates over all the objects in the cityjson file and translates it of a vertical value
    in order to make the building match the terrain level using the get_z_translation function
    then the vertices are updated.
    when buildings are adjacent, the common vertices may still need to be translated differently
    thus new vertices are added to the cityjson file 
    Args:
    Returns:

    """
    for id in cityjson2D['CityObjects'].keys():
        if cityjson2D['CityObjects'][id]['type'] != 'TINRelief':
            vertices_footprint = encoder.get_footprint_vertices(old_cityjson2D, id)

            z_translation, difference = get_z_translation(vertices_footprint, vertices, triangles)

            z_translation = int(z_translation)
            difference = float(difference)*cityjson2D['transform']['scale'][0]
            geometry = cityjson2D['CityObjects'][id]['geometry']
            cityjson2D['CityObjects'][id]['attributes']['terrainDifference'] = {'value': str(difference)} 
            #check if length is 0 otherwise it means we're not working on a 2D file and we have to translate all vertices (it's not impossible but it means iterating over all vertices and common vertices start being tricky!)
            if len(geometry)==1:
                lod0_boundary = cityjson2D['CityObjects'][id]['geometry'][0]['boundaries'][0][0]
                counter = 0
                for index in lod0_boundary:
                    vertex = cityjson2D['vertices'][index]
                    if vertex[2]==0:
                        vertex[2]=z_translation
                    elif vertex[2]!=0 and vertex[2]!=z_translation:   #attention: if two footprints share a point we still want them to be translated differently so we need to check whether the level is still 0 or not. if not a new vertex has to be created and added.
                        new_vertex = copy.deepcopy(vertex)
                        new_vertex[2] = z_translation
                        cityjson2D['vertices'].append(new_vertex)
                        cityjson2D['CityObjects'][id]['geometry'][0]['boundaries'][0][0][counter] = len(cityjson2D['vertices']) - 1
                        
                    counter += 1
    return cityjson2D