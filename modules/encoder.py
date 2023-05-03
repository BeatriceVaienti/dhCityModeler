def clean_gdf(gdf):
    """
    Remove rows with empty geometries and convert MultiPolygons to Polygons.

    Args:
        gdf (geopandas.GeoDataFrame): The GeoDataFrame to be cleaned.

    Returns:
        geopandas.GeoDataFrame: The cleaned GeoDataFrame.
    """
    gdf = gdf[gdf.geometry.notnull()] # Drop rows with empty geometries
    gdf['geometry'] = gdf['geometry'].apply(lambda geom: geom.geoms[0] if geom.type == 'MultiPolygon' else geom) # Convert MultiPolygons to Polygons

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

