import modules.parameters as parameters
import modules.roofs as roofs
import modules.walls as walls
import modules.utils as utils
import cadquery as cq
import networkx as nx
import numpy as np
import time
import copy
import logging
import math
import random


def scale_back_vertices(cityjson, vertices):
    """
    The function takes as input the CityJSON file
    which stores information about the initial scaling applied to the model
    and applies that scaling back to the vertices of the footprint
    This way we have the vertices expressed in meters

    Args:
    - cityjson: a CityJSON file (dictionary)
    - vertices: a list of 3D coordinates (list of lists)

    Returns:
    - a new list of 3D coordinates, with the vertices expressed in the original unit of measure
    """
    scale = cityjson['transform']['scale']
    vertices_copy = [[int(v[0] * (1 / scale[0])), int(v[1] * (1 / scale[1])), int(v[2] * (1 / scale[2]))] for v in vertices]

    return vertices_copy

def gdf_to_CityJSON(gdf):
    """ 
    The function uses the functions gdf_to_CityJSON_geometry and gdf_to_CityJSON_attributes 
    to create a 2D CityJSON file from a GeoDataFrame (just the encoding of the 2D geometry and the accepted attributes)
    the gdf is also cleaned before being converted to CityJSON
    remember to change the mapping in gdf_to_CityJSON_geometry if you want to map other attributes from your gdf


    Args:
        gdf (geopandas.GeoDataFrame): The GeoDataFrame to be converted.

    Returns:
        str: The resulting CityJSON file.
    """
    
    gdf_new=clean_gdf(gdf)
    cityjson=gdf_to_CityJSON_geometry(gdf_new)
    cityjson=gdf_to_CityJSON_attributes(gdf_new, cityjson)
    return cityjson

def clean_gdf(gdf):
    """
    Remove rows with empty geometries and convert MultiPolygons to Polygons.

    Args:
        gdf (geopandas.GeoDataFrame): The GeoDataFrame to be cleaned.

    Returns:
        geopandas.GeoDataFrame: The cleaned GeoDataFrame.
    """
    gdf = gdf[gdf.geometry.notnull()] # Drop rows with empty geometries
    gdf['geometry'] = gdf['geometry'].apply(lambda geom: geom.geoms[0] if geom.geom_type == 'MultiPolygon' else geom) # Convert MultiPolygons to Polygons

    return gdf

def gdf_to_CityJSON_attributes(gdf, cityjson, parameters_to_check = ["height", "numberOfFloors", "roof.type", 'roof.parameters.slope', 'floorHeight']):
    """
    Converts the attributes of a GeoDataFrame to CityJSON attributes.
    Args:
        - gdf (geopandas.GeoDataFrame): The GeoDataFrame
        - cityjson (dict): The CityJSON dictionary to be updated
        - parameters_to_check (list): The list of parameters to be extracted
    Returns:
        - cityjson (dict): The updated CityJSON dictionary
    """
    for i in range(len(gdf)):
        
        all_attributes = {}
        for par in parameters_to_check:
            attribute_dict = build_attribute_dict(gdf, i, [par])
            all_attributes = recursive_dict_update(all_attributes, attribute_dict)
        
        city_object_id = "id-" + str(i)
        city_object = cityjson["CityObjects"][city_object_id]
        city_object_attributes = city_object.get("attributes", {})
        city_object_attributes.update(all_attributes)
        city_object_attributes["updatedGeometry"] = True
        city_object["attributes"] = city_object_attributes
    return cityjson

def gdf_to_CityJSON_geometry(gdf):
    """ 
    this function doesn't add any fields (unless its a necessary field for the 3D modelling). 
    It's working as a MAPPER between a geojson geometry to a cityjson geometry
    Args:
        - gdf (geopandas.GeoDataFrame): The GeoDataFrame to be converted.

    Returns:
        - cityjson: The resulting CityJSON file.
    """
    epsg = gdf.crs.to_epsg()
    all_vertices = []
    city_objects = []

    for _,row in gdf.iterrows():

        geometry_coords_new=[]
        geometry_coords = row['geometry'].exterior.coords

        #for every point in geometry_coords from the element check if the element is equal to 
        # the previous one and in that case dont append it to the new list of coordiates
        
        geometry_coords_new = remove_duplicate_points(geometry_coords)
        twod_geometry_indices = coords_to_index_list(geometry_coords_new, all_vertices)[::-1]
        #Preparing attributes to be appended to object dictionary
        all_attributes = {}

        #Creating the CityObject dictionary 
        object = {'type': 'Building',
                'attributes':all_attributes,
                'geometry': [{'type': 'MultiSurface',
                                'lod': '0',
                                'boundaries': [[twod_geometry_indices]]}]
                        }
        city_objects.append(object)

    # Converting the vertices as indicated by CityJSON format (https://www.cityjson.org/specs/1.1.3/#transform-object)
    SCALE_FACTOR = (0.001, 0.001, 0.001)
    GEOG_EXTENT = (
        min([i[0] for i in all_vertices]),
        min([i[1] for i in all_vertices]),
        min([i[2] for i in all_vertices]),
        max([i[0] for i in all_vertices]),
        max([i[1] for i in all_vertices]),
        max([i[2] for i in all_vertices])
    )
    TRANSL_FACTOR = (GEOG_EXTENT[0], GEOG_EXTENT[1], GEOG_EXTENT[2])

    converted_all_vertices = ((np.array(all_vertices)-np.array(TRANSL_FACTOR))/np.array(SCALE_FACTOR)).astype(int)
    converted_all_vertices = [list([int(x) for x in i]) for i in converted_all_vertices]

    cityjson = create_cityjson_dictionary(SCALE_FACTOR, TRANSL_FACTOR, GEOG_EXTENT, epsg, city_objects, converted_all_vertices)

    return cityjson

def extract_attribute(gdf, i, par):
    """
    Extracts a single attribute from the GeoDataframe and creates a single attribute from it.
    Args:
        - gdf (geopandas.GeoDataFrame): The GeoDataFrame 
        - i (int): The index of the row to be extracted
        - par (str): The name of the parameter to be extracted
    Returns:
        - attribute: The resulting CityJSON parameter as a parameter dataclass
    
    """
    attribute = None
    paradata = None
    sources = []
    
    par_name = par.split(".")[-1]
    attribute = parameters.Parameter(name=par_name)
    paradata = parameters.Paradata()
    source = parameters.Source()
    
    for j, element in enumerate(gdf.iloc[i]):
        field_name = gdf.iloc[i].index[j]
        if type(element) == np.int64:
            element = int(element)
        if field_name.startswith(par):
            field_subname = field_name.split(par + ".")[-1]
            attribute.fill(field_subname, element)
            
            if field_name.startswith(par + ".paradata"):
                field_paradata = field_name.split("paradata.")[-1]
                paradata.fill(field_paradata, element)
                attribute.fill("paradata", paradata.create_dictionary())
                
            if field_name.startswith(par + ".sources"):
                source_num = int(field_name.split("sources.")[1].split(".")[0])
                if source_num != len(sources):
                    sources.append(source.create_dictionary())
                    source = parameters.Source(name=None)
                source_info = field_name.split("sources.")[1].split(".")[1]
                source.fill(source_info, element)
    
    sources.append(source.create_dictionary())
    attribute.fill("sources", sources)
    
    return attribute

def remove_duplicate_points(coords):
    """
    Removes consecutive duplicate points in a coordinate list.

    :param coords: list of coordinates
    :type coords: list[list[float]]
    :return: coordinate list with consecutive duplicate points removed
    :rtype: list[list[float]]
    """
    unique_coords = [coords[0]]
    for coord in coords[1:]:
        if coord != unique_coords[-1]:
            unique_coords.append(coord)
    return unique_coords

def coords_to_index_list(coords, vertices):
    """
    Determines which are the vertices of a geometry and returns the indices of those vertices in the vertices list.
    If they are not already in the vertices list, adds them.
    Args:
        - coords: list of coordinates of the geometry
        - vertices: list of vertices of the CityJSON
    Returns:
        - vert_indices: list of indices of the vertices of the geometry
    """
    vert_indices = []
    for point in coords:
        #converting to integer, and adding z=0
        point = [point[0], point[1], 0]
        if point not in vertices:
            vertices.append(point)
        vert_indices.append(vertices.index(point))
    return vert_indices[:-1]

def round_and_scale_coords(multisurface, CityJSON):
    """
    Rounds the coordinates in a multisurface back to the decimal places 
    contained in the CityJSON (scale value) and scales it back to an integer 
    calling the "scale_back_vertices" function.
    
    Args:
        multisurface (list): A multisurface in meters.
        CityJSON (dict): A dictionary that contains the scale factor used to 
        scale the vertices.

    Returns:
        list: The multisurface with rounded and scaled coordinates.
    """
    multisurface_copy = copy.deepcopy(multisurface)
    scale = CityJSON['transform']['scale']
    decimal_places = max(-math.floor(math.log10(scale[i])) for i in range(3))

    for f in range(len(multisurface_copy)):
        for w in range(len(multisurface_copy[f])):
            for v in range(len(multisurface_copy[f][w])):
                multisurface_copy[f][w][v]=list(multisurface_copy[f][w][v])
                for c in range(len(multisurface_copy[f][w][v])):
                    multisurface_copy[f][w][v][c]=round(multisurface_copy[f][w][v][c],decimal_places)
            multisurface_copy[f][w]=(scale_back_vertices(CityJSON, multisurface_copy[f][w]))
    return multisurface_copy

def scale_vertices(cityjson, vertices):
    """
    Scale the given vertices based on the scaling factor specified in the CityJSON file.

    Args:
        cityjson (dict): A CityJSON dictionary containing the scaling information.
        vertices (list): A list of vertex coordinates in the format [(x1, y1, z1), (x2, y2, z2), ...].

    Returns:
        list: A new list of scaled vertex coordinates.
    """
    scaled_vertices = copy.deepcopy(vertices)
    scaling_factor = cityjson['transform']['scale']

    for i, vertex in enumerate(scaled_vertices):
        scaled_vertices[i] = tuple(round(coord * scaling_factor[dim], int(-math.log10(scaling_factor[dim]))) for dim, coord in enumerate(vertex))

    return scaled_vertices

def vertices_to_index(multisurface, cityjson):
    """
    Given a multisurface geometry as a list of lists of coordinate tuples,
    and a CityJSON object with a 'vertices' list, the function converts
    the coordinates of each vertex in the multisurface to its corresponding
    index in the 'vertices' list in CityJSON, and returns the new multisurface
    with indices instead of coordinates.

    Args:
    - multisurface (list): a multisurface geometry as a list of lists of
        coordinate tuples
    - CityJSON (dict): a CityJSON object with a 'vertices' list

    Returns:
    - multisurface_indices (list): the multisurface geometry with each coordinate
        replaced by its corresponding index in the 'vertices' list in CityJSON
    """
    vertices = cityjson['vertices']
    multisurface_indices = []

    for surface in multisurface:
        surface_indices = []
        for ring in surface:
            ring_indices = []
            for point in ring:
                if point in vertices:
                    index = vertices.index(point)
                else:
                    vertices.append(point)
                    index = len(vertices) - 1
                ring_indices.append(index)
            surface_indices.append(ring_indices)
        multisurface_indices.append(surface_indices)

    return multisurface_indices

def create_cityjson_dictionary(SCALE_FACTOR, TRANSL_FACTOR, GEOG_EXTENT, epsg, city_objects, converted_vertices):
    """
    Creates a CityJSON dictionary with the required fields.

    Args:
        SCALE_FACTOR (list): The scale factor.
        TRANSL_FACTOR (list): The translation factor.
        GEOG_EXTENT (list): The geographical extent.
        epsg (int): The EPSG code.
        city_objects (dict): The CityObjects dictionary.
        converted_vertices (list): The vertices list.
    Returns:
        dict: The CityJSON dictionary.
    """
    
    cityjson = {
        "type": "CityJSON",
        "version": "1.1",
        "transform": {
        "scale": [
            SCALE_FACTOR[0],
            SCALE_FACTOR[1],
            SCALE_FACTOR[2]
        ],
        "translate": [
            TRANSL_FACTOR[0],
            TRANSL_FACTOR[1],
            TRANSL_FACTOR[2]
        ]
        },
        "metadata": {
        "geographicalExtent": [
            GEOG_EXTENT[0], GEOG_EXTENT[1], GEOG_EXTENT[2],
            GEOG_EXTENT[3], GEOG_EXTENT[4], GEOG_EXTENT[5]
        ],
        "referenceSystem": "https://www.opengis.net/def/crs/EPSG/0/"+str(epsg)
        },
        
        "CityObjects": {f'id-{i}' : obj for i,obj in enumerate(city_objects)},
        "vertices": converted_vertices
    }

    return cityjson

def build_attribute_dict(gdf, i, parameters_to_check):
    """
    Builds a dictionary of attributes from the GeoDataFrame.
    Args:
        - gdf (geopandas.GeoDataFrame): The GeoDataFrame
        - i (int): The index of the row to be extracted
        - parameters_to_check (list): The list of parameters to be extracted
    Returns:
        - attribute_dict: The resulting dictionary of attributes

    """
    attribute_dict = {}
    for par in parameters_to_check:
        nested_structure = par.split(".")
        if len(nested_structure) == 1:
            attribute_dict[nested_structure[0]] = extract_attribute(gdf, i, par).create_dictionary()
        else:
            cur_dict = attribute_dict
            for k in range(len(nested_structure)):
                if nested_structure[k] not in cur_dict:
                    cur_dict[nested_structure[k]] = {}
                if k == len(nested_structure) - 1:
                    cur_dict[nested_structure[k]] = extract_attribute(gdf, i, par).create_dictionary()
                cur_dict = cur_dict[nested_structure[k]]
    return attribute_dict

def recursive_dict_update(dict1, dict2):
    """
    Recursively update dict1 with dict2. If a key in dict2 is already present in dict1,
    recursively update the value associated with that key in dict1.
    Args:
        - dict1 (dict): The dictionary to be updated
        - dict2 (dict): The dictionary to update dict1 with
    Returns:
        - result_dict (dict): The updated dictionary
    """
    result_dict = copy.deepcopy(dict1)
    for key in dict2:
        if key in result_dict and isinstance(result_dict[key], dict) and isinstance(dict2[key], dict):
            result_dict[key] = recursive_dict_update(result_dict[key], dict2[key])
        else:
            result_dict[key] = dict2[key]
    return result_dict

def get_footprint_vertices(cityjson, id):
    """
    A function that loads a 2D or 3D CityJSON and given a index of a building 
    it returns the LOD0 geometry (i.e. the buildings footprint) 
    replacing the vertices indeces with the actual coordinates:
    given a CityJSON file and the ID of a CityObject,
    the function returns the coordinates of the vertices of the footprint
    NB: the value that is returned is the same one that is stored inside of the cityJSON and does not apply the translation and the scaling of the CityJSON
    therefore the coords are all integers

    Args:
        cityjson (dict): A 2D or 3D CityJSON file.
        id (str): The ID of the CityObject whose footprint coordinates are to be returned.

    Returns:
        coords (list of lists): The coordinates of the vertices of the footprint, in the form
                                 [[x1, y1, z1], [x2, y2, z2], ...].
    """
    #get the vertices of the footprint whose id is given
    vertices = cityjson['vertices']
    foundLOD0=False

    #get the coordinates of the vertices of the footprint
    #coords = [vertices[i] for i in CityJSON['CityObjects'][id]['geometry'][0]['boundaries'][0]]
    for geom in cityjson['CityObjects'][id]['geometry']:
        if geom['lod']=='0':   
            coords= [vertices[i] for i in  geom['boundaries'][0][0]]
            foundLOD0=True
        
        if not foundLOD0:
            coords=[]
    return coords[::-1]

def add_edge_to_graph(e, G, wire_edges_starts):
    """
    Given a wire edge, adds the edge to the graph G 
    only if the edge is not composed by coincident vertices,
    and returns the updated graph and the list of the starts of the wire edges

    Args:
        - e (wire edge): a wire edge (CQ object)
        - G (graph): a graph (networkx object)
        - wire_edges_starts (list): the list of the starts of the wire edges

    Returns:
        - G (graph): the graph with the new edge added
        - wire_edges_starts (list): the list of the starts of the wire edges
    """
    start, end = e.Vertices()
    if start != end:
        start = start.toTuple()
        end = end.toTuple()
        wire_edges_starts.append(start)
        G.add_edge(start, end)
    return G, wire_edges_starts

def check_direction_change(wire_sorted_vertices, face_normal_rounded):
    """
    Checks if the direction of the wire needs to be changed based on the dot product between the sum of cross products
    of the wire edges and the face normal.

    Args:
    - wire_sorted_vertices: a list of vertices that define the wire
    - face_normal_rounded: the normal vector of the face, rounded to a small number of decimal places

    Returns:
    - dot_product: the dot product between the sum of cross products of the wire edges and the face normal
    """
    cross_product_sum = np.zeros(3)
    for i in range(len(wire_sorted_vertices) // 2):
        indeces=sorted(random.sample(range(len(wire_sorted_vertices)), 3))
        vec_edge=wire_sorted_vertices[indeces[1]][0]-wire_sorted_vertices[indeces[0]][0], wire_sorted_vertices[indeces[1]][1]-wire_sorted_vertices[indeces[0]][1], wire_sorted_vertices[indeces[1]][2]-wire_sorted_vertices[indeces[0]][2]
        vec_edge_2=wire_sorted_vertices[indeces[2]][0]-wire_sorted_vertices[indeces[1]][0], wire_sorted_vertices[indeces[2]][1]-wire_sorted_vertices[indeces[1]][1], wire_sorted_vertices[indeces[2]][2]-wire_sorted_vertices[indeces[1]][2]
        cross_product = np.cross(vec_edge, vec_edge_2)
        cross_product_sum+=cross_product
    dot_product=np.dot(cross_product_sum, face_normal_rounded)
    return dot_product

def build_wire_graph(face_edges, index, wire_edges_starts):
    """
    Builds a graph from the edges of a face and the index of the wire that the edges belong to.

    Args:
        - face_edges: a list of edges that define the face
        - index: the index of the wire that the edges belong to

    Returns:
        - G: a graph that represents the face
        - wire_edges_starts: a list of the starting vertices of the edges of the wire
    """
    G = nx.Graph() 
    for e in face_edges[index]:
        G, wire_edges_starts = add_edge_to_graph(e, G, wire_edges_starts)
    return G, wire_edges_starts

def process_plane_face(face, processed_faces):
    """
    Given a face and a list of the faces that were already processed
    the function returns a list of vertices that define the face
    according to CityJSON structure for multisurfaces and sorted along the faces boundary

    Args:
        - face: a face of a cqSolid
        - processed_faces: a list of the faces that were already processed
    Returns:
        - face_vertices_new: a list of vertices that define the face
    """
    wire = cq.Workplane().add(face).wires().vals()
    face_normal = face.normalAt().normalized().toTuple()
    face_normal_rounded = (round(face_normal[0], 3), round(face_normal[1], 3), round(face_normal[2], 3))
    face_vertices_new = []
    face_edges = []
    more_than_two_vertices = 0

    for index, w in enumerate(wire):
        face_edges.append(cq.Workplane().add(w).edges().vals())
        wire_edges_starts = []
        
        len_face_edges = cq.Workplane().add(w).edges().size()
        if len_face_edges <= 2:
            continue

        G, wire_edges_starts = build_wire_graph(face_edges, index, wire_edges_starts)
        wire_sorted_vertices = list(nx.dfs_preorder_nodes(G, wire_edges_starts[0]))

        if len(wire_sorted_vertices) <= 2:
            continue  # Ignore wires with 2 or fewer vertices

        more_than_two_vertices += 1  
        wire_sorted_vertices_new = copy.deepcopy(wire_sorted_vertices)
        if check_direction_change(wire_sorted_vertices, face_normal_rounded) < 0:
            wire_sorted_vertices_new.reverse()
        face_vertices_new.append(wire_sorted_vertices_new)

    if more_than_two_vertices > 0 and face not in processed_faces:
        return face_vertices_new
    return None

def process_triangle(vertices, triangle, curved_surface_center):
    """
    Processes a triangle and returns a list of vertices that define the triangle
    according to CityJSON structure for multisurfaces and sorted along the faces boundary

    Args:
        - vertices: a list of vertices that define the face
        - triangle: a triangle from the face tesselation
        - curved_surface_center: the center of the curved surface
    
    Returns:
        - face_vertices_new: a list of vertices that define the face
    """
    face_vertices = [vertices[triangle[i]].toTuple() for i in range(3)]
    # remove triangles that have two vertices that are the same:
    if face_vertices[0]!=face_vertices[1] and face_vertices[0]!=face_vertices[2] and face_vertices[1]!=face_vertices[2]:
        #create a polygon from the vertices and find its center
        triangle_polyline=cq.Workplane().polyline([vertices[triangle[0]], vertices[triangle[1]], vertices[triangle[2]]])
        triangle_center=triangle_polyline.val().Center().toTuple()
        sub_face=cq.Face.makeFromWires(triangle_polyline.wire().val())
        
        # Calculate the normal vector of the triangle
        vector_side_1=face_vertices[1][0]-face_vertices[0][0], face_vertices[1][1]-face_vertices[0][1], face_vertices[1][2]-face_vertices[0][2]
        vector_side_2=face_vertices[2][0]-face_vertices[1][0], face_vertices[2][1]-face_vertices[1][1], face_vertices[2][2]-face_vertices[1][2]
        normal_np=list(np.cross(vector_side_1, vector_side_2))
        normal = cq.Vector(normal_np).normalized()

        # Determine if the normal vector points towards the curved surface center
        vector_to_center=curved_surface_center[0]-triangle_center[0], curved_surface_center[1]-triangle_center[1], curved_surface_center[2]-triangle_center[2]
        vector_to_center=cq.Vector(vector_to_center).normalized()
        direction_check=np.dot(normal.toTuple(), vector_to_center.toTuple())

        if direction_check>0:
            face_vertices.reverse()
        return face_vertices, sub_face
    else:
        return None, None

def get_multisurface(cqSolid, tessellation_level=8):
    """
    The function takes as input a cqSolid
    and returns a list of the faces of the cqSolid
    each face is composed by a list of vertices
    according to CityJSON structure for multisurfaces and sorted along the faces boundary

    Args:
        - cqSolid: a cqSolid
        - tessellation_level: the level of tessellation of the curved surfaces
    
    Returns:
        - multisurface: a list of the faces of the cqSolid
    """
    multisurface = []           # Initialize an empty list of multisurfaces
    processed_faces = set()     # Initialize a set to keep track of processed faces

    for face in cqSolid.faces().vals():

        if face.geomType() == 'PLANE' or face.geomType() == 'BSPLINE':
            face_vertices_new = process_plane_face(face, processed_faces)

            if face_vertices_new:
                multisurface.append(face_vertices_new)
                processed_faces.add(face)
        else:
            curved_surface_center=face.Center().toTuple() 
            try: 
                vertices, triangles = face.tessellate(tessellation_level)
            except:
                vertices, triangles = face.tessellate(1)

            for triangle in triangles:
                face_vertices_new, sub_face = process_triangle(vertices, triangle, curved_surface_center)
                if face_vertices_new is not None and sub_face not in processed_faces:
                    multisurface.append([face_vertices_new])
                    processed_faces.add(sub_face)

    return multisurface

def fill_geom_attribute(attribute_name, inherited_geometric_attributes, geometric_attributes):
    """
    The function checks if a geometric attribute is already present in the CityJSON file (inherited geometric attributes).
    If not, it adds it to the geometric_attributes dictionary with a None value (or with a 0 in the case of the terrain_difference attribute); 
    This is an internal function that is called by generate_LOD1_geom_attributes.

    Args:
        - attribute_name: name of the attribute
        - inherited_geometric_attributes: dictionary of inherited geometric attributes
        - geometric_attributes: dictionary of geometric attributes
    Returns:
        - geometric_attributes: dictionary of geometric attributes
    """
    if attribute_name not in inherited_geometric_attributes:
        if attribute_name == 'terrainDifference':
            geometric_attributes[attribute_name] = {'value': '0'}
        else:
            geometric_attributes[attribute_name] = {'value': None}
    else:
        geometric_attributes[attribute_name] = inherited_geometric_attributes[attribute_name]
    return geometric_attributes

def generate_LOD1_geom_attributes(id, cityjson):
    """
    The function checks if all the attributes that are necessary to generate a LOD1 geometry are present in the CityJSON.
    If not, it adds them to the geom_attributes dictionary with a None value; # NaN value in dataframe 
    the value will be then generated through the procedural modelling script. # todo: NO it will be filled in here
    This is an internal function that is called by generate_LOD1.

    Args:
        - id: id of the building
        - cityjson: CityJSON file
    Returns:
        - geomAttributes: dictionary of geometric attributes
    """
    # Retrieve the geometric features that are already present for this building in the CityJSON
    inherited_geom_attributes = cityjson['CityObjects'][id].get('attributes', {})
    
    # Initialize the dictionary of geometry attributes with the building ID
    geom_attributes = {"id": id}
    
    
    possible_attributes_lod1 = ['height', 'numberOfFloors', 'floorHeight', 'terrainDifference']
    for attribute in possible_attributes_lod1:
        geom_attributes = fill_geom_attribute(attribute, inherited_geom_attributes, geom_attributes)

    return geom_attributes

def generate_LOD1(id, cityjson):
    """
    The function generates the LOD1 geometry of a building and updates the geom_attributes dictionary.
    The dictionary takes track of which values had to be randomly generated following the Historical CityJSON extension schema.
    Args:
        - id: id of the building
        - cityjson: CityJSON file
    Returns:
        - LOD1: LOD1 geometry of the building

    """
    geom_attributes=generate_LOD1_geom_attributes(id, cityjson)
    
    vertices= get_footprint_vertices(cityjson, id)      #get the vertices of the footprint
    
    new_vertices=scale_vertices(cityjson, vertices)

    logging.info("geom_attributes['height']:")
    logging.info(geom_attributes['height'])
    LOD1, new_attributes= walls.base_walls(new_vertices, height= geom_attributes['height']['value'], numberOfFloors= geom_attributes['numberOfFloors']['value'], terrain_difference= float(geom_attributes['terrainDifference']['value']))
    
    geom_attributes.update(new_attributes)  

    return LOD1, geom_attributes

def generate_LOD2(id, cityjson, LOD1, geom_attributes_LOD1):
    """
    Generate the LOD2 geometry and attributes of a building
    Args:
        id (str): id of the building
        cityjson (CityJSON): CityJSON object of the building
        LOD1 : LOD1 geometry of the building
        geom_attributes_LOD1 (dict): LOD1 geometry attributes of the building
    Returns:
        LOD2 geometry and attributes of the building
    """
    # Extract geometry attributes and generate roof attributes
    geomAttributes=geom_attributes_LOD1
    roof_geomAttributes=generate_LOD2_geom_attributes(id, cityjson)
    roof_geomAttributes.update(geomAttributes)
    total_geom_attributes = roof_geomAttributes
    roof_type = total_geom_attributes['roof']['type']['value']
    roof_parameters = total_geom_attributes['roof']['parameters']

    # Generate the roof geometry depending on its type
    coords_top_face=[]
    top_faces=LOD1.faces(">Z").edges().vertices().vals()
    for i in range(len(top_faces)):
        coords_top_face.append(top_faces[i].toTuple())
    if roof_type =="flat":
        roof, roof_new_attributes=roofs.flat_roof(
            coords_top_face, 
            roof_parameters['baseFloorThickness']['value'], 
            roof_parameters['railingHeight']['value'], 
            roof_parameters['railingWidth']['value'])

        total_geom_attributes['roof']['parameters'].update(roof_new_attributes)
    elif roof_type =="domed":
        roof, roof_new_attributes=roofs.domed_roof(
            coords_top_face, 
            roof_parameters['baseFloorThickness']['value'], 
            roof_parameters['railingHeight']['value'], 
            roof_parameters['railingWidth']['value'], 
            roof_parameters['domeVerticalParameter']['value'], 
            roof_parameters['domeHorizontalParameter']['value'] )

        total_geom_attributes['roof']['parameters'].update(roof_new_attributes)
    elif roof_type=="hip":
        roof_base = utils.remove_nearly_parallel(coords_top_face)
        roof, roof_new_attributes=roofs.hip_roof(
            roof_base, 
            roof_parameters['slope']['value'])

        total_geom_attributes['roof']['parameters'].update(roof_new_attributes)
        try:
            shell, shell_new_attributes= roofs.roof_shell(
                roof, 
                total_geom_attributes['roof'], 
                roof_parameters['upperFloorThickness']['value'], 
                roof_parameters['eavesOverhang']['value'])
            total_geom_attributes['roof']['parameters'].update(shell_new_attributes)
        except:
            shell=None
        if shell != None:
            roof=roof.union(shell)
    elif roof_type=="gable":
        roof_base = utils.remove_nearly_parallel(coords_top_face)
        roof, roof_new_attributes=roofs.gable_roof(
            roof_base, 
            slope = roof_parameters['slope']['value'],
            ids_gable = roof_parameters['gableSides']['value'])
        total_geom_attributes['roof']['parameters'].update(roof_new_attributes)

        try:
            shell, shell_new_attributes= roofs.roof_shell(
                roof, 
                total_geom_attributes['roof'], 
                roof_parameters['upperFloorThickness']['value'], 
                roof_parameters['eavesOverhang']['value'])
            total_geom_attributes['roof']['parameters'].update(shell_new_attributes)
        except:
            shell=None
        if shell != None:
            roof=roof.union(shell)
    try:
        LOD2=LOD1.union(roof)
    except:
        LOD2=LOD1
    
    return LOD2, total_geom_attributes

def update_cityjson_geometry(cityjson_modified, id_building):
    """
    The function takes as input a CityJSON file and the id of a building
    and it recalculates the geometry of that building in LOD1 and LOD2
    returning a new CityJSON file with the updated geometry
    Args:
        - cityjson_modified: CityJSON file with the modified attributes
        - id_building: id of the building to be updated
    Returns:
        - cityjson: CityJSON file with the updated geometry
    """
    #check if the object is a building:
    if cityjson_modified['CityObjects'][id_building]['type'] != 'Building':
        raise Exception('The object is not a building')

    cityjson = copy.deepcopy(cityjson_modified)
    sub_dict = cityjson['CityObjects'][id_building]['geometry']

    def update_lod(lod, lod_value):
        multisurface = get_multisurface(lod)
        multisurface2 = round_and_scale_coords(multisurface, cityjson)
        multisurface3 = vertices_to_index(multisurface2, cityjson)

        lod_dict = {
            "type": "MultiSurface",
            "lod": str(lod_value),
            "boundaries": multisurface3
        }

        lod_exists = any(d['lod'] == lod_value for d in sub_dict)

        if lod_exists:
            index = next((i for i, d in enumerate(sub_dict) if d['lod'] == lod_value), None)
            cityjson['CityObjects'][id_building]['geometry'][index] = lod_dict
        else:
            cityjson['CityObjects'][id_building]['geometry'].append(lod_dict)

    LOD1, geomFeatures_LOD1 = generate_LOD1(id_building, cityjson)
    if True:
        LOD2, geomFeatures_LOD2 = generate_LOD2(id_building, cityjson, LOD1, geomFeatures_LOD1)
    else:
        cityjson['CityObjects'][id_building]['attributes'].update(geomFeatures_LOD1)
        update_lod(LOD1, 1)
        return cityjson
    cityjson['CityObjects'][id_building]['attributes'].update(geomFeatures_LOD2)

    update_lod(LOD1, 1)
    update_lod(LOD2, 2)

    return cityjson

def generate_LOD2_geom_attributes(id, CityJSON, roofTypes=['hip', 'gable', 'flat', 'domed']):
    """
    The function checks if all the necessary parameters to generate a roof are present inside the cityObject 
    and if not, it adds them to the geomAttributes dictionary.
    If the type of roof is not present, a random function determines it here, taking track of the possible values among which the value was chosen.
    These possible values can be personalized in the function roofTypes.
    
    Args:
        id (str): the id of the building
        CityJSON (dict): the CityJSON file
        roofTypes (list, optional): the possible values for the roof type. Defaults to ['hip', 'gable', 'flat', 'domed'].
    
    Returns:
        geomAttributes (dict): the dictionary containing the geomAttributes of the building
    """
    geomAttributes = {"id": id}
    inherited_geomAttributes = CityJSON['CityObjects'][id].get('attributes', {})
    geomAttributes.update(inherited_geomAttributes)
    
    if 'roof' not in geomAttributes:
        roofType = random.choice(roofTypes)
        roof = {
            'type': {
                'value': roofType,
                'metadata': {
                    'source': f"random{roofTypes}"
                }
            }
        }
        geomAttributes['roof'] = {'parameters': {}}
        geomAttributes['roof'].update(roof)
    else:
        roofType = geomAttributes['roof']['type']['value']
    
    if 'parameters' not in geomAttributes['roof']:
        geomAttributes['roof']['parameters'] = get_default_parameters(roofType)
    else:
        parameters = geomAttributes['roof']['parameters']
        default_parameters = get_default_parameters(roofType)
        for key in default_parameters.keys():
            if key not in parameters.keys():
                parameters[key] = default_parameters[key]
        remove_extra_parameters(parameters, roofType)

    return geomAttributes

def get_default_parameters(roofType):
    """
    Returns the default parameters based on the given roof type.

    Args:
        roofType (str): the type of roof

    Returns:
        dict: the default parameters
    """
    default_parameters = {}
    if roofType == 'hip':
        default_parameters = {
            'slope': {'value': None},
            'upperFloorThickness': {'value': None},
            'eavesOverhang': {'value': None}
        }
    elif roofType == 'gable':
        default_parameters = {
            'slope': {'value': None},
            'upperFloorThickness': {'value': None},
            'eavesOverhang': {'value': None},
            'gableSides': {'value': None}
        }
    elif roofType == 'flat':
        default_parameters = {
            'baseFloorThickness': {'value': None},
            'railingHeight': {'value': None},
            'railingWidth': {'value': None}
        }
    elif roofType == 'domed':
        default_parameters = {
            'baseFloorThickness': {'value': None},
            'railingHeight': {'value': None},
            'railingWidth': {'value': None},
            'domeVerticalParameter': {'value': None},
            'domeHorizontalParameter': {'value': None}
        }
    return default_parameters

def remove_extra_parameters(parameters, roofType):
    """
    Removes any extra parameters that are not required based on the roof type.

    Args:
        parameters (dict): the parameters of the roof
        roofType (str): the type of roof

    """
    valid_parameters = {
        'hip': ['slope', 'upperFloorThickness', 'eavesOverhang'],
        'gable': ['slope', 'upperFloorThickness', 'eavesOverhang', 'gableSides'],
        'flat': ['baseFloorThickness', 'railingHeight', 'railingWidth'],
        'domed': ['baseFloorThickness', 'railingHeight', 'railingWidth', 'domeVerticalParameter', 'domeHorizontalParameter']
    }
    for key in list(parameters.keys()):
        if key not in valid_parameters.get(roofType, []):
            parameters.pop(key)
