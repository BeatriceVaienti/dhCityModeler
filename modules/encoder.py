import modules.parameters as parameters
import numpy as np
import copy

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


def create_CityJSON_dictionary(SCALE_FACTOR, TRANSL_FACTOR, GEOG_EXTENT, epsg, city_objects, converted_vertices):
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
            counter=0
            
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

    cityjson = create_CityJSON_dictionary(SCALE_FACTOR, TRANSL_FACTOR, GEOG_EXTENT, epsg, city_objects, converted_all_vertices)

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


def gdf_to_CityJSON_attributes(gdf, cityjson, parameters_to_check = ["height", "numberOfFloors", "roof.type", "roof.parameters.slope"]):
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

