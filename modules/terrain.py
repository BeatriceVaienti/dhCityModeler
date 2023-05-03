# for this module to work you need to install numpy, pyproj and pydelatin
"""
This scripts manage the download and processing of the terrain data.
in order for it to work the following dependencies are necessary:
- numpy
- pyproj
- pydelatin # conda install -c conda-forge pydelatin
- scipy # conda install -c conda-forge scipy
- meshio # optional: needed only for direct export of obj terrain
- cjio # optional: needed only for cleaning unused vertices
"""
import requests
import json
import numpy as np 
from pyproj import Transformer
from pydelatin import Delatin 
from scipy.spatial import Delaunay 
import numpy as np
from scipy.interpolate import LinearNDInterpolator
import modules.encoder as encoder
import meshio 
import copy
import subprocess
import cjio 

def rescale_positions(vertices: np.ndarray, bounds, flip_y: bool = False):
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

def get_elevations(cityjson2D, apikey, grid_density = 20): 
    """
    Get the elevation data corresponding to a certain bounding box (xmin and ymin)
    and with a given grid density (regardless of the shape of the bounding rectangle, 
    the grid is divided in an equal number of rows and columns)
    Args:
        - cityjson2D: the 2D CityJSON file 
        - apikey: the GOOGLE API key
        - grid_density: desired minimum distance between two points of the grid in meters
    Returns:
        - elevations: a list containing the elevation of the points in the grid
        - n_points: the number of points per side
    """
    xmin_original, ymin_original, xmax_original, ymax_original = get_cityjson_bounds(cityjson2D)

    min_side = min(xmax_original-xmin_original, ymax_original-ymin_original) # we get the min side of the bounding box
    n_points = int(min_side/grid_density) # we compute the number of points we want to have in the grid
    
    in_epsg = "EPSG:" + cityjson2D['metadata']['referenceSystem'].split('/')[-1] # Define the EPSG codes for the original CRS and the target CRS
    out_epsg = "EPSG:4326" # EPSG used by the Google Elevation API
    transformer = Transformer.from_crs(in_epsg, out_epsg, always_xy=True)  # pyproj transformer

    lats = np.linspace(ymin_original, ymax_original, n_points)
    lons = np.linspace(xmin_original, xmax_original, n_points)
    xx, yy = np.meshgrid(lons, lats)
    points = np.column_stack((xx.ravel(), yy.ravel()))
    # Query the Google Elevation API for the elevations of the grid points
    elevations = []
    for point in points:
        lat, lon = transformer.transform(point[1], point[0])
        service_url = f"https://maps.googleapis.com/maps/api/elevation/json?locations={lat},{lon}&key={apikey}"
        response = requests.get(service_url)
        data = json.loads(response.text)
        elevation = data["results"][0]["elevation"]
        elevations.append(elevation)
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


def export_TIN(vertices, triangles, folder_path, filename): # set filename without format!
    """
    Exports the mesh corresponding to a TIN triangulation as an .obj file 

    Args:
        - vertices
        - triangles
        - path: path to the folder where to save the file
        - filename: name of the file without the extension

    """
    cells = [("triangle", triangles)]
    mesh = meshio.Mesh(vertices, cells)
    full_path = folder_path + "/" + filename + ".obj"
    mesh.write(full_path)

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

def create_delaunay(vertices):
    """
    Creates a delaunay tassellation of the vertices using Scipy
    Args:
        - vertices 
    Returns:
        - triangulation
    """
    vertices_array = np.array(vertices)
    tri = Delaunay(vertices_array[:, :2])
    return tri

def get_z_value(x, y, vertices):
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
    tri = create_delaunay(vertices)   # Create a Delaunay triangulation of the vertices
    
    # Find the index of the triangle that contains the (x, y) point
    i = tri.find_simplex([x, y])
    if i < 0:
        raise ValueError("Point outside of triangulation")
    else:
        # Get the vertices of the triangle
        #triangle_vertices = vertices_array[triangles[i]]
        interpolator = LinearNDInterpolator(tri, z_values)
        z = interpolator(x, y)
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


def get_z_translation(cityjson2d_terrain, id_building, terrain_vertices):
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
    vertices_footprint = encoder.get_footprint_vertices(cityjson2d_terrain, id_building)
    z_value = vertices_footprint[0][2]
    z_final = None
    if z_value == 0:
        #print('the LOD0 of the building is at level 0')
        at_terrain_level = False
        for vertex in vertices_footprint:  #LOOKING FOR XY VALUES OF THE VERTICES OF THE BUILDING
            x, y = vertex[0], vertex[1]
            z_new = get_z_value(x, y, terrain_vertices)
            if z_final == None:
                z_final = z_new
            elif z_new < z_final:
                z_final = z_new
        #we translate the building of z_final
        return z_final  


def get_terrain_from_cityjson_bounds(cityjson2D, apikey, grid_density = 20, tin_max_error = 1):
    """
    this function could be furtherly divided!
    Args:
    Returns:
    """
    elevations, n_points = get_elevations(cityjson2D, apikey, grid_density) #we make the request to the GOOGLE api to get the elevation data
    vertices, triangles = create_TIN(elevations, n_points, tin_max_error) #We generate a TIN out of the elevation grid

    return vertices,triangles
    
def insert_terrain_in_cityjson(cityjson2D, apikey, grid_density = 20, tin_max_error = 1):
    """
    Inserts the terrain in the cityjson2D file and translates the buildings accordingly

    Args:
        - cityjson2D: the cityjson2D file
        - apikey: your google api key
        - grid_density: number of points per side of the grid
        - tin_max_error: maximum error allowed for the TIN
    Returns:
        - cityjson2D: the cityjson2D file with the terrain inserted
    """
    vertices, triangles = get_terrain_from_cityjson_bounds(cityjson2D, apikey, grid_density, tin_max_error)
    old_cityjson2D = copy.deepcopy(cityjson2D)
    GEOG_EXTENT= cityjson2D['metadata']['geographicalExtent']
    SCALE_FACTOR = cityjson2D['transform']['scale']
    TRANSL_FACTOR =  cityjson2D['transform']['translate']
    bounds= [GEOG_EXTENT[0], GEOG_EXTENT[1], GEOG_EXTENT[3], GEOG_EXTENT[4]] #bounds of the cityjson file

    rescaled_vertices = rescale_positions(vertices, bounds).tolist()
    converted_all_vertices = ((np.array(rescaled_vertices)-np.array(TRANSL_FACTOR))/np.array(SCALE_FACTOR)).astype(int)
    converted_all_vertices = [list([int(x) for x in i]) for i in converted_all_vertices]

    triangles_list = triangles.tolist()
    new_triangles_list = []
    for triangle in triangles_list:
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

    terrain_cityobject = create_terrain_cityobject(rescaled_vertices, new_triangles_list)

    cityjson2D['CityObjects'].update(terrain_cityobject) #add the new object to the cityjson file
    cityjson2D = translate_lod0_geometries(cityjson2D, old_cityjson2D, converted_all_vertices)
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

def translate_lod0_geometries(cityjson2D, old_cityjson2D, vertices):
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
            z_translation = int(get_z_translation(old_cityjson2D, id, vertices))
            geometry = cityjson2D['CityObjects'][id]['geometry']
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

def remove_unused_vertices(cityjson2d, folder_path="./", name="cityjson"): #careful, it saves the file . also it may not be necessary now.
    """
    Uses cjio to remove unused vertices from the cityjson file
    it saves the file in the folder_path with the name specified
    Args:  
        - cityjson2d: the cityjson2d file
        - folder_path: the path where the file will be saved
        - name: the name of the file without extension
    Returns:
        - cityjson2d_new: the cleaned cityjson
    """
    full_path = folder_path + name + ".json"
    with open(full_path, 'w') as outfile:
        json.dump(cityjson2d, outfile)
    subprocess.check_output(["cjio", full_path, "vertices_clean"])
    with open(full_path) as f:
        cityjson2d_new = json.load(f)
    return cityjson2d_new