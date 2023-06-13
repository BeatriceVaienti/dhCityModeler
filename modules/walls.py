import cadquery as cq
import random
import modules.utils as utils
import copy
import modules.parameters as parameters

def complete_base_walls_parameters(height=None, numberOfFloors=None, numberOfFloors_random_range=[2,4], floorHeight=None, default_floor_height_range=[2.8, 3.5]):
    """
    Completes the necessary parameters used to create the LOD1 geometry
    calculates the height based on the availability of data on the number of floors and the floor height or directly of the height
    inferring the values of the missing parameters from the available ones

    The best case is that where we know the number of floors and the floor height. in this case the floor height is not used and is merely inferred for later used
    if we know the height but not the number of floors or the floor height, we can infer the number of floors by a default floor height and readjust the floor height to obtain an integer number of floors
    if we know the number of floors but not the height or the floor height, we can infer the height by multiplying the number of floors by a default floor height
    if we know the number of floors and the floor height then the overall height is simply inferrecd by multiplying the number of floors by the floor height
    if we don't know anything, the number of floor is assigned randomly in a given range and the height is inferred by multiplying the number of floors by a default floor height 

    Args:
        height (float, optional): the height of the building. Defaults to None.
        numberOfFloors (int, optional): the number of floors of the building. Defaults to None.
        numberOfFloors_random_range (list, optional): the range of the random number of floors. Defaults to [2,4].
        floorHeight (float, optional): the height of a single floor. Defaults to None.
        default_floor_height_range (list, optional): the range of the random floor height. Defaults to [2.8, 3.5].
    
    Returns:
        new_attributes (dict): a dictionary containing the new attributes to be added to the building
    """
    new_attributes = {}
    if height is not None:
        height= float(height)

    if numberOfFloors is not None:
        numberOfFloors = int(numberOfFloors)
        if numberOfFloors == 0:
            numberOfFloors = 1
    
    if numberOfFloors is not None and height is not None:
        numberOfFloors_source = 'user'
        height_source = 'user'
        floorHeight = height / numberOfFloors
        floorHeight_source = parameters.Source(name="inferred", type="automatic", notes="inferred from height and number of floors")
        floorHeight_dict = {"value": str(floorHeight), "sources": [floorHeight_source.create_dictionary()]}
        new_attributes["floorHeight"] = floorHeight_dict
    elif height is not None and floorHeight is None and numberOfFloors is None:
        height_source = "user"
        numberOfFloors_source = parameters.Source(name="inferred", type="automatic", notes="inferred from height and default height of a single floor")
        floorHeight = round(random.uniform(default_floor_height_range[0], default_floor_height_range[1]),3)
        floorHeight_source =  parameters.Source(name="random", type="automatic", notes="Random value picked in the default range: "+str(default_floor_height_range))
        numberOfFloors = round(height / floorHeight)
        if numberOfFloors == 0:
            numberOfFloors = 1
        floorHeight = height / numberOfFloors
        numberOfFloors_dict={"value": str(numberOfFloors), "sources": [numberOfFloors_source.create_dictionary()]}
        floorHeight_dict={"value": str(floorHeight), "sources": [floorHeight_source.create_dictionary()]}
        new_attributes["floorHeight"]=floorHeight_dict
        new_attributes["numberOfFloors"]=numberOfFloors_dict
    elif numberOfFloors is not None and height is None and floorHeight is None:
        numberOfFloors_source = "user"
        height_source = parameters.Source(name="inferred", type="automatic", notes="inferred from the number of floors and default height of a single floor")
        floorHeight = round(random.uniform(default_floor_height_range[0], default_floor_height_range[1]),3)
        floorHeight_source = parameters.Source(name="random", type="automatic", notes="Random value picked in the default range: "+str(default_floor_height_range))
        height = numberOfFloors * floorHeight
        height_dict = {"value": str(height), "sources": [height_source.create_dictionary()]}
        floorHeight_dict={"value": str(floorHeight), "sources": [floorHeight_source.create_dictionary()]}
        new_attributes["floorHeight"]=floorHeight_dict
        new_attributes["height"]=height_dict
    elif numberOfFloors is not None and floorHeight is not None and height is None:
        numberOfFloors_source = "user"
        floorHeight_source = "user"
        height_source = parameters.Source(name="inferred", type="automatic", notes="inferred from the number of floors and assigned height of a single floor")
        height = numberOfFloors * floorHeight
        height_dict = {"value": str(height), "sources": [height_source.create_dictionary()]}
        new_attributes["height"]=height_dict
    elif numberOfFloors is None and floorHeight is None and height is None:
        numberOfFloors_source = parameters.Source(name="inferred", type="automatic", notes="Random number picked in the default range: "+str(numberOfFloors_random_range))
        height_source = parameters.Source(name="inferred", type="automatic", notes="inferred from the number of floors and assigned height of a single floor")
        floorHeight_source = parameters.Source(name="random", type="automatic", notes="Random value picked in the default range: "+str(default_floor_height_range))
        floorHeight = round(random.uniform(default_floor_height_range[0], default_floor_height_range[1]),3)
        numberOfFloors = random.randint(numberOfFloors_random_range[0], numberOfFloors_random_range[1])
        if numberOfFloors == 0:
            numberOfFloors = 1
        height = numberOfFloors * floorHeight
        height_dict = {"value": str(height), "sources": [height_source.create_dictionary()]}
        new_attributes["height"]=height_dict
        floorHeight_dict={"value": str(floorHeight), "sources": [floorHeight_source.create_dictionary()]}
        new_attributes["floorHeight"]=floorHeight_dict
        numberOfFloors_dict={"value": str(numberOfFloors), "sources": [numberOfFloors_source.create_dictionary()]}
        new_attributes["numberOfFloors"]=numberOfFloors_dict

    return new_attributes

def base_walls(points, height=None, numberOfFloors=None, numberOfFloors_range=[2,4], floorHeight=None, terrain_difference = 0, default_floor_height_range=[2.8, 3.5]):
    """
    Extrudes a 2D shape to create a 3D extrusion of it.
    The resulting shape can then be used to detect the top face and generate a roof.
    The function also adds the attributes "height", "numberOfFloors" and "floorHeight" to the building.

    Args:
        points (list): list of points defining the shape to extrude
        height (float, optional): the height of the extrusion. Defaults to None.
        numberOfFloors (int, optional): the number of floors of the building. Defaults to None.
        numberOfFloors_random_range (list, optional): the range of random values to pick the number of floors from. Defaults to [2,4].
        floorHeight (float, optional): the height of a single floor. Defaults to None.
        default_floor_height_range (list, optional): the range of random values to pick the height of a single floor from. Defaults to [2.8, 3.5].

    Returns:
        tuple: a tuple containing the extruded shape and the new attributes

    """

    points_copy=copy.deepcopy(points)
    points_copy, z_base = utils.fix_points_get_z(points_copy)
    new_attributes = complete_base_walls_parameters(height, numberOfFloors, numberOfFloors_range, floorHeight, default_floor_height_range)
    if 'height' in new_attributes:
        new_height = float(new_attributes["height"]["value"])
    else: 
        new_height = height

    total_height = float(terrain_difference) + new_height
    base_walls=cq.Workplane("XY", (0,0,0)).polyline(points_copy).close().extrude(total_height) 
    
    if z_base!=0:  #if the walls are not starting from z=0 we move them back up:
        base_walls = base_walls.translate((0,0,z_base))
    return base_walls, new_attributes