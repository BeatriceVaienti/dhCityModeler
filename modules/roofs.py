import cadquery as cq
import math
import modules.utils as utils
#import modules.polylabel as polylabel
import modules.parameters as parameters
import copy
import random

def get_points_for_loft(Pn, y, vector, slope, x_min):
    """
    This function returns the points that will be used to create the polyline that then will be used for the lofted surface

    Args:
        Pn (list): point (list of three coordinates) on the base shape.     
        y (float): initial arbitrarily high height of the extrusion used for the subtraction operations.
        vector (Vector): vector that defines the direction of the slope. 
        slope (float): slope expressed in percentage [0,1].
        x_min (float): offset distance used to avoid self intersections.
    Returns:
        P1, P2, P3 (list): list of three points that will be used to create the polyline that then will be used for the lofted surface
    """
    P1 = [Pn[0],Pn[1], Pn[2]+y]
    P2 = [Pn[0]+vector.x, Pn[1]+vector.y, Pn[2]+y]
    P3 = [Pn[0]+vector.x, P1[1]+vector.y, Pn[2]+slope*x_min]
    return P1, P2, P3

def create_hip_attributes_dictionary(slope, default_slope_range=[0.3,0.6]):
    """
    This function creates the dictionary that contains the metadata information about the parameters that were used to create the hip roof.
    If the slope is not assigned, the default slope is randomly picked from a range of values (default_slope_range) with a uniform distribution.
    Args:
        slope (float): slope expressed in percentage [0,1].
        default_slope_range (list): range of values from which the slope is randomly picked if not assigned.
    Returns:
        roof_new_attributes (dict): dictionary that contains the metadata information about the parameters that were used to create the roof.
        slope (float): slope expressed in percentage [0,1].
    """

    roof_new_attributes={}
    if slope is None:
        slope = round(random.uniform(default_slope_range[0], default_slope_range[1]),3)
        slope_source = parameters.Source(name="random", type="automatic", notes="Random value picked in the default range: "+str(default_slope_range))
        slope_dict = {"value": str(slope), "sources": [slope_source.create_dictionary()]}
        roof_new_attributes["slope"] = slope_dict
    else:
        slope=float(slope)
    return roof_new_attributes, slope
    
def create_subtraction_loft(An, Bn, Cn, v_a, v_b, v_c, y, slope, x_min):
    """
    The function creates the subtraction loft for the hip roof 
    around a given point Bn and using the two segments AB and BC on its sides

    """
    A1, A2, A3 = get_points_for_loft(An, y, v_a, slope, x_min)
    B1, B2, B3 = get_points_for_loft(Bn, y, v_b, slope, x_min)
    C1, C2, C3 = get_points_for_loft(Cn, y, v_c, slope, x_min)


    loft1=(cq.Workplane("XY",(0,0,0)).polyline([An,A1,A2,A3]).close().polyline([Bn,B1,B2,B3]).close()).loft(ruled=True) 
    loft2=(cq.Workplane("XY",(0,0,0)).polyline([Bn,B1,B2,B3]).close().polyline([Cn,C1,C2,C3]).close()).loft(ruled=True)
    loft1 = loft1.union(loft2, clean=True, glue=True)
    return loft1

def hip_roof(points, slope=None, default_slope_range=[0.3,0.6]):
    """
    The function creates a hip roof given a list of boundary points and the slope expressed in percentage [0,1].
    If not assigned, the default slope is randomly picked from a range of values (default_slope_range)  
    Returns a CadQuery solid and a dictionary with the parameters used to generate the shape and their origin

    Args:
        points (list): list of points that define the boundary of the roof.
        slope (float): slope expressed in percentage [0,1].
        default_slope_range (list): range of values from which the slope is randomly picked if not assigned.
    Returns:
        roof (CadQuery solid): CadQuery solid representing the roof.
        roof_new_attributes (dict): dictionary that contains the metadata information about the parameters that were used to create the roof.
    """
    points_copy=copy.deepcopy(points)
    points_copy, z_base = utils.fix_points_get_z(points_copy)
    #compile the metadata information about the parameters that were used to create the roof
    roof_new_attributes, slope = create_hip_attributes_dictionary(slope, default_slope_range)
    
    y=70 #y is the height of extrusion used for the subtraction operations, it should be higher than the highest ridge. since we are considering roofs we will suppose them to be lower than 100 m
    roof=cq.Workplane("XY",(0,0,0)).polyline(points_copy).close().extrude(y)  #extruded volume, used to subtract the sloped parts and obtain the final roof

    for i in range(0, len(points_copy), 1):
        points4=(points_copy[i:i+4] if i+4<= len(points_copy) else points_copy[i:]+points_copy[:i+4-len(points_copy)]) # starting from the current point take the three subsequent points
        # the three points on the curve are called An Bn and Cn
        An, Bn, Cn, Dn = points4      #we consider the fourth point as well to check if the angles on the sides of BC are both convex 

        # we determine the offset distance in order to avoid self intersections    (x=L*tg(alpha/2) , alpha=angle between BA and BC)
        V_BA = cq.Vector(An[0]-Bn[0], An[1]-Bn[1], An[2]-Bn[2]) # vector from B to A
        V_BC = cq.Vector(Cn[0]-Bn[0], Cn[1]-Bn[1], Cn[2]-Bn[2]) # vector from B to C
        V_CB = cq.Vector(Bn[0]-Cn[0], Bn[1]-Cn[1], Bn[2]-Cn[2]) # vector from C to B
        V_CD = cq.Vector(Dn[0]-Cn[0], Dn[1]-Cn[1], Dn[2]-Cn[2]) # vector from C to D
        alpha =utils.angle_between(V_BA,V_BC) # angle between BA and BC
        L_BA = math.sqrt(V_BA.x**2+ V_BA.y**2 + V_BA.z**2) 
        L_BC = math.sqrt(V_BC.x**2 + V_BC.y**2 + V_BC.z**2)
        B_is_concave =utils.cross_sign(V_BA,V_BC) #checks whether the point Bn corresponds to a concavity
        C_is_concave =utils.cross_sign(V_CB,V_CD) #checks whether the point Cn corresponds to a concavity
        # based on the concavity and on the angle formed by the two vectors we determine the factor that will be multiplied to the lenght of the sides to determine x_A and x_C
        if B_is_concave:
            fac=(math.tan(alpha/2)*10) #if the angle is concave the distance can be arbitrary big, so we multiply the obtained value with a factor 5
        elif round(math.tan(alpha/2),2)>=1 and not B_is_concave:
            fac=(math.tan(alpha/2)*0.99)
        else:
            fac=(math.tan(alpha/2)*0.99)

        x_min = min(L_BA*fac, L_BC*fac, 100) #Offset distance # in case the two adjacent sides are nearly parallel, the distance would tend to infinity, so we set a maximum value of 100

        v_a=utils.perpendicular_vector(An,Bn).multiply(x_min) #offset vector for side AB
        v_c=utils.perpendicular_vector(Bn,Cn).multiply(x_min) #offset vector for side BC
        v1v2=v_a.add(v_c) #sum of v1 and v2, it gives me the correct direction (bisecting vector)
        Lv1v2=math.sqrt(v1v2.x**2+ v1v2.y**2+ v1v2.z**2) #length of the vector sum of v1 and v2
        v_b=v1v2.normalized().multiply(2*x_min**2/Lv1v2) # the lenght of v3 is determined

        loft = create_subtraction_loft(An, Bn, Cn, v_a, v_b, v_c, y, slope, x_min)
        roof = roof.cut(loft,clean=True)

        if not C_is_concave and not B_is_concave:
            lenght_vect=100
            perp_vect=utils.perpendicular_vector(Bn,Cn).multiply(lenght_vect)
            B1, B2, B3 = get_points_for_loft(Bn, y, perp_vect, slope, lenght_vect)
            C1, C2, C3 = get_points_for_loft(Cn, y, perp_vect, slope, lenght_vect)
            extra_subtraction=(cq.Workplane("XY",(0,0,0)).polyline([Bn,B1,B2,B3]).close().polyline([Cn,C1,C2,C3]).close()).loft(ruled=True)
            roof = roof.cut(extra_subtraction,clean=True)
        
    roof = roof.translate((0,0,z_base))  #move the roof back to the height z_base
    return roof, roof_new_attributes

def create_gable_attributes_dictionary(slope, ids_gable, slope_range):
    roof_new_attributes = {}
    if slope is None:
        slope = round(random.uniform(slope_range[0], slope_range[1]),3)
        slope_source = parameters.Source(name="random", type="automatic", notes="Random value picked in the default range: "+str(slope_range))
        slope_dict = {"value": str(slope), "sources": [slope_source.create_dictionary()]}
        roof_new_attributes["slope"] = slope_dict
    else:
        slope=float(slope)

    if ids_gable is None:
        gable_source = "automatic"
        gable_source = parameters.Source(name="default", type="automatic", notes="The sides with the vertical side are automatically identified as the ones with only three vertices")
        gable_dict = {"value": None, "sources": [gable_source.create_dictionary()]}
        roof_new_attributes["idsGable"] = gable_dict

    return roof_new_attributes, slope

def make_triangular_face(A, B, C):
                            wire = cq.Workplane("XY",(0,0,0)).polyline([A,B,C]).close().wire().val()
                            face = cq.Face.makeFromWires(wire)
                            return face 

def gable_roof(points, slope=None, default_slope_range=[0.3,0.6], ids_gable=None):
    """
    Creates a gable roof given a list of boundary points, the slope expressed in percentage [0,1] and optionally the index of the sides that have a vertical side.
    The index is expressed as the first point that is enclosing the segment, following the order of apparition of the points.
    If no index is provided, then the gable sides are automatically identified as those for which the gable roof presents a sloped side that has only three vertices. 
    Returns a CadQuery solid and a dictionary with the parameters used to generate the shape and their origin

    Args:
        points (list): list of points that define the boundary of the roof
        slope (float, optional): slope of the roof expressed in percentage [0,1]. Defaults to None.
        default_slope_range (list, optional): range of values used to randomly pick the slope if no slope is provided. Defaults to [0.3,0.6].
        ids_gable (list, optional): list of indices of the sides that have a vertical side. Defaults to None.

    Returns:
        roof: the CadQuery solid representing the roof and a dictionary with the parameters used to generate the shape and their origin
        roof_new_attributes: dictionary with the parameters used to generate the shape and their origin
    """

    roof_new_attributes, slope = create_gable_attributes_dictionary(slope, ids_gable, slope_range = default_slope_range)
    points_copy=copy.deepcopy(points)
    points_copy, z_base = utils.fix_points_get_z(points_copy)

    y=100 #y is the height of extrusion used for the subtraction operations, it should be higher than the highest ridge. since we are considering roofs we will suppose them to be lower than 100 m

    if ids_gable == None or ids_gable == []:
        roof=hip_roof(points_copy, slope)[0]
        sloped_faces= roof.faces("not|Z").vals()
        for sloped_face in sloped_faces:
            n_sides=len(sloped_face.Edges())
            if n_sides==3:
                sloped_face_cq = cq.Workplane().add(sloped_face)
                lowest_point=sloped_face_cq.vertices(">>Z[0]").val().Center() 
                if lowest_point.z==0.0:  #the lowest point of the sloped face needs to be on the boundary of the roof       
                    B=sloped_face_cq.vertices("<<Z[0]").val().Center() #B is the top vertex of the triangular side
                    AC=sloped_face_cq.vertices("<<Z[1]").vals()

                    if len(AC)==2: # if the sloped face is the correct one, it has only two vertices on the base
                        A=AC[0].Center()
                        C=AC[1].Center()

                        v_CA = cq.Vector(A.x-C.x, A.y-C.y, 0)
                        v_CB = cq.Vector(B.x-C.x, B.y-C.y, 0)
                        gamma = utils.angle_between(v_CA,v_CB)

                        v_AC = cq.Vector(C.x-A.x, C.y-A.y, 0)
                        v_AB = cq.Vector(B.x-A.x, B.y-A.y, 0)
                        alpha = utils.angle_between(v_AC,v_AB)

                        rho = gamma - alpha + (math.pi/2)
                        theta = (math.pi/2) - alpha
                        xi = (math.pi/2) - gamma
                        
                        pos_B1 = (math.sin(xi)*math.sin(gamma))/(math.sin(rho)*math.sin(theta+xi)) # pos_B1 is the reparametrized position of B' on AC (|AB'|/|AC|)
                        x_B1 = A.x + pos_B1*(C.x-A.x)
                        y_B1 = A.y + pos_B1*(C.y-A.y)
                        #if the sides before and after AC are not parallel, then the z of B1 is different from that of B and we need to determine it
                        if 2*alpha+2*gamma != math.pi:
                            L_AB1 = math.sqrt((x_B1-A.x)**2 + (y_B1-A.y)**2)
                            L_AB = math.sqrt((B.x-A.x)**2 + (B.y-A.y)**2)
                            pos_B1z = (L_AB1/L_AB)*(math.sin(math.pi-2*alpha)/math.sin(math.pi-alpha))  # |B1 B2|/|B2 B| projected on plane xy (distance calculated considering the xy values of the points)
                            z_B1 = pos_B1z*(B.z)  # [...]*(B2.z-B.z) but B2.z=0 
                            B1=[x_B1,y_B1,z_B1] 
                        else:
                            B1 = [x_B1, y_B1, B.z]
                
                        face2 = make_triangular_face(A,B1,C)
                        face3 = make_triangular_face(C,B1,B)
                        face4 = make_triangular_face(A,B1,B)
                            
                        shell = cq.Shell.makeShell([sloped_face,face2,face3,face4]).fix()
                        solid = cq.Solid.makeSolid(shell)
                        roof = roof.union(solid, clean=True, glue= True)
    else: #in case we manually defined the gable sides
        roof=cq.Workplane("XY", (0,0,0)).polyline(points_copy).close().extrude(y)  #extruded volume, used to subtract the sloped parts and obtain the final roof
        for i in range(0, len(points_copy), 1):
            points4=(points_copy[i:i+4] if i+4<= len(points_copy) else points_copy[i:]+points_copy[:i+4-len(points_copy)]) # starting from the current point take the three subsequent points
            An, Bn, Cn, Dn = points4 
            # we determine the offset distance in order to avoid self intersections
            # (x=L*tg(alpha/2) , alpha=angle between BA and BC)
            # we determine the offset distance in order to avoid self intersections    (x=L*tg(alpha/2) , alpha=angle between BA and BC)
            V_BA = cq.Vector(An[0]-Bn[0], An[1]-Bn[1], An[2]-Bn[2]) # vector from B to A
            V_BC = cq.Vector(Cn[0]-Bn[0], Cn[1]-Bn[1], Cn[2]-Bn[2]) # vector from B to C
            V_CB = cq.Vector(Bn[0]-Cn[0], Bn[1]-Cn[1], Bn[2]-Cn[2]) # vector from C to B
            V_CD = cq.Vector(Dn[0]-Cn[0], Dn[1]-Cn[1], Dn[2]-Cn[2]) # vector from C to D
            alpha =utils.angle_between(V_BA,V_BC) # angle between BA and BC
            L_BA = math.sqrt(V_BA.x**2+ V_BA.y**2 + V_BA.z**2) 
            L_BC = math.sqrt(V_BC.x**2 + V_BC.y**2 + V_BC.z**2)
            B_is_concave =utils.cross_sign(V_BA,V_BC) #checks whether the point Bn corresponds to a concavity
            C_is_concave =utils.cross_sign(V_CB,V_CD) #checks whether the point Cn corresponds to a concavity
            # based on the concavity and on the angle formed by the two vectors we determine the factor that will be multiplied to the lenght of the sides to determine x_A and x_C
            if B_is_concave:
                fac=(math.tan(alpha/2)*10) #if the angle is concave the distance can be arbitrary big, so we multiply the obtained value with a factor 5
            elif round(math.tan(alpha/2),2)>=1 and not B_is_concave:
                fac=(math.tan(alpha/2)*0.99)
            else:
                fac=(math.tan(alpha/2)*0.99)
       
            x_min = min(L_BA*fac, L_BC*fac, 100) #Offset distance

            v1=utils.perpendicular_vector(An,Bn).multiply(x_min) #offset vector for side AB
            v2=utils.perpendicular_vector(Bn,Cn).multiply(x_min) #offset vector for side BC
            v1v2=v1.add(v2) #sum of v1 and v2, it gives me the correct direction (bisecting vector)           
            Lv1v2=math.sqrt(v1v2.x**2+ v1v2.y**2+ v1v2.z**2) #length of the vector sum of v1 and v2
            v3=v1v2.normalized().multiply(2*x_min**2/Lv1v2) # the lenght of v3 is determined
        
            if i+1 in ids_gable  or (i==len(points_copy)-1 and 0 in ids_gable):
                #V_BC_new=V_BC.normalized().multiply(x_min/math.cos(angle_between(V_BC,v1)))
                v1_new = v1.normalized().multiply(L_BC*math.cos(utils.angle_between(V_BC,v1)))
                x_new = L_BC*math.cos(utils.angle_between(V_BC,v1))
                B1, B2, B3 = get_points_for_loft(Bn, y, V_BC, slope, x_new)
                A1, A2, A3 = get_points_for_loft(An, y, v1_new, slope, x_new)
                #if the current point is the one before the gable point we create a face in An and then extrude it until Bn (from An to Bn)
                #the face is created by the four points An, A1, A2 and A3
                loft = (cq.Workplane("XY",(0,0,0)).polyline([An,A1,A2,A3]).close().polyline([Bn,B1,B2,B3]).close()).loft(ruled=True) 
            elif i in ids_gable:
                v2_new = v2.normalized().multiply(L_BA*math.cos(utils.angle_between(V_BA,v2)))
                x_new = L_BA*math.cos(utils.angle_between(V_BA,v2)) 
                B1, B2, B3 = get_points_for_loft(Bn, y, V_BA, slope, x_new)
                C1, C2, C3 = get_points_for_loft(Cn, y, v2_new, slope, x_new)

                loft = (cq.Workplane("XY",(0,0,0)).polyline([Bn,B1,B2,B3]).close().polyline([Cn,C1,C2,C3]).close()).loft(ruled=True) 
            else:
                A1, A2, A3 = get_points_for_loft(An, y, v1, slope, x_min)
                C1, C2, C3 = get_points_for_loft(Cn, y, v2, slope, x_min)
                B1, B2, B3 = get_points_for_loft(Bn, y, v3, slope, x_min)
	
                loft1=(cq.Workplane("XY",(0,0,0)).polyline([An,A1,A2,A3]).close().polyline([Bn,B1,B2,B3]).close()).loft(ruled=True) 
                loft2=(cq.Workplane("XY",(0,0,0)).polyline([Bn,B1,B2,B3]).close().polyline([Cn,C1,C2,C3]).close()).loft(ruled=True)
                loft=loft1+loft2
            roof=roof-loft 
    #move the roof back to the height z_base
    roof = roof.translate((0,0,z_base))
    return roof, roof_new_attributes

def create_shell_attributes_dictionary(roof_vert_thickness, roof_horiz_offset, default_roof_vert_thickness_range, default_roof_horiz_offset_range):
    """
    Creates a dictionary with the attributes of the shell of the roof

    Args:
        roof_vert_thickness (float): the thickness of the vertical part of the roof
        roof_horiz_offset (float): the offset of the horizontal part of the roof
        default_roof_vert_thickness_range (list): the default range of the roof vertical thickness
        default_roof_horiz_offset_range (list): the default range of the roof horizontal offset

    Returns:
        roof_new_attributes (dict): the dictionary with the attributes of the shell of the roof
    """
    roof_new_attributes = {}
    if roof_vert_thickness is None:
        roof_vert_thickness = round(random.uniform(default_roof_vert_thickness_range[0],default_roof_vert_thickness_range[1]),3)
        roof_vert_thickness_source = parameters.Source(name="random", type="automatic", notes="Random value picked in the default range: "+str(default_roof_vert_thickness_range))
        roof_vert_thickness_dict = {"value": str(roof_vert_thickness), "sources": [roof_vert_thickness_source.create_dictionary()]}
        roof_new_attributes["upperFloorThickness"] = roof_vert_thickness_dict
    else:
        roof_vert_thickness = float(roof_vert_thickness)
    
    if roof_horiz_offset is None:
        roof_horiz_offset = round(random.uniform(default_roof_horiz_offset_range[0],default_roof_horiz_offset_range[1]),3)
        roof_horiz_offset_source = parameters.Source(name="random", type="automatic", notes="Random value picked in the default range: "+str(default_roof_horiz_offset_range))
        roof_horiz_offset_dict = {"value": str(roof_horiz_offset), "sources": [roof_horiz_offset_source.create_dictionary()]}
        roof_new_attributes["eavesOverhang"] = roof_horiz_offset_dict
    else:
        roof_horiz_offset = float(roof_horiz_offset)
    return roof_new_attributes, roof_vert_thickness, roof_horiz_offset

def roof_shell(simple_roof, roof_metadata, roof_vert_thickness = None, roof_horiz_offset=None, default_roof_vert_thickness_range=[0.2,0.3], default_roof_horiz_offset_range=[0.2,0.4], ids_gable='automatic'):
    """
    Creates a more detailed representation for hip and gable roofs
    constructing a shell with a vertical thickness equal to roof_vert_thickness
    And an horizontal offset equal to roof_horiz_offset

    Args:
        simple_roof (cq object): the simple roof
        roof_metadata (dict): the metadata of the roof
        roof_vert_thickness (float, optional): the thickness of the vertical part of the roof. Defaults to None.
        roof_horiz_offset (float, optional): the offset of the horizontal part of the roof. Defaults to None.
        default_roof_vert_thickness_range (list, optional): the default range of the roof vertical thickness. Defaults to [0.2,0.3].
        default_roof_horiz_offset_range (list, optional): the default range of the roof horizontal offset. Defaults to [0.2,0.4].   
        ids_gable (str, optional): the id of the gable roof. Defaults to 'automatic'.
    
    Returns:
        roof (cq object): the shell of the roof
        roof_new_attributes (dict): the dictionary with the attributes of the shell of the roof
    """
    roof_new_attributes, roof_vert_thickness, roof_horiz_offset = create_shell_attributes_dictionary(roof_vert_thickness, roof_horiz_offset, default_roof_vert_thickness_range, default_roof_horiz_offset_range)
    offset=simple_roof.faces("-Z").edges().toPending().offset2D(roof_horiz_offset, "intersection").vertices().vals()
    big_base = [coord.toTuple() for coord in offset]


    if roof_metadata['type']['value']=='hip':
        slope=roof_metadata["parameters"]["slope"]['value']
        big_roof=hip_roof(big_base, slope)[0]
    elif roof_metadata['type']['value']=="gable":
        slope=roof_metadata["parameters"]["slope"]['value']
        if ids_gable == 'automatic' or roof_metadata["parameters"]["ids_gable"]["sources"][0]["name"]=="automatic":
            big_roof=gable_roof(big_base, slope)[0]     
        elif roof_metadata["parameters"]["ids_gable"]["source"][0]["name"]=="user":
            big_roof=gable_roof(big_base, slope, ids_gable=roof_metadata["parameters"]["ids_gable"]["value"])[0]

    outer_surface=big_roof.faces("(not |Z)  and (not #Z)").each(lambda f: f if f.normalAt().z>0 else None # 
        ,False, False, True).faces() 
                                            # select the outer set of sloped faces
    max_ridge_line=outer_surface.edges(">Z").val().locationAt(0.5).toTuple()[0][2]
    min_ridge_line=simple_roof.edges(">Z").val().locationAt(0.5).toTuple()[0][2]

    vert_dist=max_ridge_line-min_ridge_line                                          # get the vertical distance between the two sides of the shell
    inner_surfaces=outer_surface.translate((0,0,-vert_dist))                         # move the outer surface down to get the inner surface
    inner_shell=cq.Shell.makeShell(inner_surfaces.vals())                            # create a shell from the inner surface (join the faces)

    #extrude the inner shell vertically with a height equal to roof_vert_width to get the final roof shell:
    roof_shell=cq.Workplane().add(inner_shell).faces().each(lambda f: cq.Solid.extrudeLinear(f,(0,0,roof_vert_thickness)),True,False,True).solids().combine() 
    return roof_shell, roof_new_attributes

def create_flat_roof_dictionary(floor_width, railing_height, railing_width, default_floor_width_range, default_railing_height_range, default_railing_width_range):
    
    roof_new_attributes = {}
    if floor_width is None:
        floor_width = round(random.uniform(default_floor_width_range[0],default_floor_width_range[1]),3)
        floor_width_source = parameters.Source(name="random", type="automatic", notes="Random value picked in the default range: "+str(default_floor_width_range))
        floor_width_dict = {"value": str(floor_width), "sources": [floor_width_source.create_dictionary()]}
        roof_new_attributes["baseFloorThickness"] = floor_width_dict
    else:
        floor_width=float(floor_width)

    if railing_height is None:
        railing_height = round(random.uniform(default_railing_height_range[0],default_railing_height_range[1]),3)
        railing_height_source = parameters.Source(name="random", type="automatic", notes="Random value picked in the default range: "+str(default_railing_height_range))
        railing_height_dict = {"value": str(railing_height), "sources": [railing_height_source.create_dictionary()]}
        roof_new_attributes["railingHeight"] = railing_height_dict
    else:
        railing_height = float(railing_height)
    
    if railing_width is None:
        railing_width = round(random.uniform(default_railing_width_range[0],default_railing_width_range[1]),3)
        railing_width_source = parameters.Source(name="random", type="automatic", notes="Random value picked in the default range: "+str(default_railing_width_range))
        railing_width_dict = {"value": str(railing_width), "sources": [railing_width_source.create_dictionary()]}
        roof_new_attributes["railingWidth"] = railing_width_dict
    else:
        railing_width = float(railing_width)
    
    return roof_new_attributes, floor_width, railing_height, railing_width  

def flat_roof(points, floor_width=None, railing_height=None, railing_width=None, default_floor_width_range=[0.25,0.35], default_railing_height_range=[0.9,1.1], default_railing_width_range=[0.1,0.2]):
    """
    the function creates a flat roof given a list of boundary points 
    optionally the vertical floor width, the railing height and width can be assigned
    otherwise their values will be randomly picked in the default ranges using a uniform distribution
    Returns a CadQuery solid and a dictionary with the parameters used to generate the shape and their origin
    
    Args:
        points: list of points that define the base of the roof
        floor_width: vertical width of the floor
        railing_height: height of the railing
        railing_width: width of the railing
        default_floor_width_range: default range for the floor width
        default_railing_height_range: default range for the railing height
        default_railing_width_range: default range for the railing width
    Returns:
        roof: CadQuery solid
        roof_new_attributes: dictionary with the parameters used to generate the shape and their origin
    """

    points_copy=copy.deepcopy(points)
    
    points_copy, z_base = utils.fix_points_get_z(points_copy)
    roof_new_attributes, floor_width, railing_height, railing_width = create_flat_roof_dictionary(floor_width, railing_height, railing_width, default_floor_width_range, default_railing_height_range, default_railing_width_range)
    
    roof_base = cq.Workplane("XY",(0,0,0)).polyline(points_copy).close().extrude(floor_width)  
    # take the outer edges of the roof and create an inner offset of railing_width and create a wall from the offset with height railing_height
    try:
        railing = roof_base.faces(">Z").wires().toPending().extrude(railing_height).faces("|Z").shell(-railing_width, "intersection")
        roof = roof_base + railing
    except:
        roof = roof_base
    #move the roof back to the height z_base
    roof = roof.translate((0,0,z_base))
    return roof, roof_new_attributes

def create_domed_roof_dictionary(floor_width, railing_height, railing_width, dome_vert_p, dome_horiz_p, default_floor_width_range, default_railing_height_range, default_railing_width_range, default_dome_vert_p_range, default_dome_horiz_p_range):
    """
    Creates a dictionary with the parameters used to generate the shape

    Args:
        floor_width: vertical width of the floor
        railing_height: height of the railing
        railing_width: width of the railing
        dome_vert_p: vertical radius of the dome (percentage)
        dome_horiz_p: horizontal radius of the dome (percentage)
        default_floor_width_range: default range for the floor width
        default_railing_height_range: default range for the railing height
        default_railing_width_range: default range for the railing width
        default_dome_vert_p_range: default range for the vertical radius of the dome
        default_dome_horiz_p_range: default range for the horizontal radius of the dome
    Returns:
        roof_new_attributes: dictionary with the parameters used to generate the shape
        floor_width: vertical width of the floor
        railing_height: height of the railing
        railing_width: width of the railing
        dome_vert_p: vertical radius of the dome (percentage)
        dome_horiz_p: horizontal radius of the dome (percentage)

    """
    roof_new_attributes = {}
    if floor_width is None:
        floor_width = round(random.uniform(default_floor_width_range[0],default_floor_width_range[1]),3)
        floor_width_source = parameters.Source(name="random", type="automatic", notes="Random value picked in the default range: "+str(default_floor_width_range))
        floor_width_dict = {"value": str(floor_width), "sources": [floor_width_source.create_dictionary()]}
        roof_new_attributes["baseFloorThickness"] = floor_width_dict
    else:
        floor_width=float(floor_width)

    if railing_height is None:
        railing_height = round(random.uniform(default_railing_height_range[0],default_railing_height_range[1]),3)
        railing_height_source = parameters.Source(name="random", type="automatic", notes="Random value picked in the default range: "+str(default_railing_height_range))
        railing_height_dict = {"value": str(railing_height), "sources": [railing_height_source.create_dictionary()]}
        roof_new_attributes["railingHeight"] = railing_height_dict
    else:
        railing_height = float(railing_height)
    
    if railing_width is None:
        railing_width = round(random.uniform(default_railing_width_range[0],default_railing_width_range[1]),3)
        railing_width_source = parameters.Source(name="random", type="automatic", notes="Random value picked in the default range: "+str(default_railing_width_range))
        railing_width_dict = {"value": str(railing_width), "sources": [railing_width_source.create_dictionary()]}
        roof_new_attributes["railingWidth"] = railing_width_dict
    else:
        railing_width = float(railing_width)
    
    if dome_vert_p is None:
        dome_vert_p = round(random.uniform(default_dome_vert_p_range[0],default_dome_vert_p_range[1]),3)
        dome_vert_p_source = parameters.Source(name="random", type="automatic", notes="Random value picked in the default range: "+str(default_dome_vert_p_range))
        dome_vert_p_dict = {"value": str(dome_vert_p), "sources": [dome_vert_p_source.create_dictionary()]}
        roof_new_attributes["domePercentVertRadius"] = dome_vert_p_dict
    else:
        dome_vert_p = float(dome_vert_p)
    
    if dome_horiz_p is None:
        dome_horiz_p = round(random.uniform(default_dome_horiz_p_range[0],default_dome_horiz_p_range[1]),3)
        dome_horiz_p_source = parameters.Source(name="random", type="automatic", notes="Random value picked in the default range: "+str(default_dome_horiz_p_range))
        dome_horiz_p_dict = {"value": str(dome_vert_p), "sources": [dome_horiz_p_source.create_dictionary()]}
        roof_new_attributes["domePercentBaseRadius"] = dome_horiz_p_dict
    else:
        dome_horiz_p = float(dome_horiz_p)
    
    return roof_new_attributes, floor_width, railing_height, railing_width, dome_vert_p, dome_horiz_p

def domed_roof(points, floor_width=None, railing_height=None, railing_width=None, dome_vert_p = None, dome_horiz_p = None, default_floor_width_range=[0.25,0.35], default_railing_height_range=[0,0.5], default_railing_width_range=[0.1,0.2], default_dome_vert_p_range=[0.4,0.7], default_dome_horiz_p_range=[0.7,0.9]):
    """
    the function creates a flat roof given a list of boundary points 
    optionally the vertical floor width, the railing height and width can be assigned
    otherwise their values will be randomly picked in the default ranges using a uniform distribution
    Returns a CadQuery solid and a dictionary with the parameters used to generate the shape and their origin
    
    Args:
        points (list): list of points that define the boundary of the roof
        floor_width (float): the vertical floor width of the roof
        railing_height (float): the height of the railing
        railing_width (float): the width of the railing
        dome_vert_p (float): the vertical radius of the dome (percentage of the height)
        dome_horiz_p (float): the horizontal radius of the dome (percentage of the base radius)

    Returns:
        CadQuery solid: domed roof
        dictionary: dictionary with the parameters used to generate the shape
    """
    points_copy=copy.deepcopy(points)
    
    points_copy, z_base = utils.fix_points_get_z(points_copy)
    roof_new_attributes, floor_width, railing_height, railing_width, dome_vert_p, dome_horiz_p = create_domed_roof_dictionary(floor_width, railing_height, railing_width, dome_vert_p, dome_horiz_p, default_floor_width_range, default_railing_height_range, default_railing_width_range, default_dome_vert_p_range, default_dome_horiz_p_range)
    #parameters:
    p_v = dome_vert_p
    p_h = dome_horiz_p

    base_center, max_radius= utils.polylabel([points], with_distance=True)
    base_center.append(0)
    
    max_radius_old=max_radius
    if railing_height>0:
        max_radius=max_radius -railing_width

    #find the actual radius
    base_radius=max_radius*p_h
    sphere_radius=utils.get_sphere_radius(base_center, base_radius, p_v)
    roof_base = cq.Workplane("XY",(0,0,0)).polyline(points_copy).close().extrude(floor_width)  
    sphere_center=(base_center[0], base_center[1],  - (sphere_radius-base_radius*p_v) + floor_width)
    sphere=cq.Workplane("XY", sphere_center).sphere(sphere_radius).workplane(sphere_radius- (base_radius*p_v)).split(keepTop=True)
    
    if railing_height>0:
        # take the outer edges of the roof and create an inner offset of railing_width and create a wall from the offset with height railing_height
        try:
            railing = roof_base.faces(">Z").wires().toPending().extrude(railing_height).faces("|Z").shell(-railing_width, "intersection")
            roof = roof_base + railing
        except:
            roof = roof_base
    else:
        roof = roof_base
    #move the roof back to the height z_base
    roof= roof.add(sphere)
    roof = roof.translate((0,0,z_base))
    return roof, roof_new_attributes
