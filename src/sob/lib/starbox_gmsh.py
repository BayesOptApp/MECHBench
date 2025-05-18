from src.sob.lib.gmsh_base_meshes import Template_GMSH_Mesh_Constructor
import gmsh
import json
from pathlib import Path
from typing import Union
import os,sys


class Starbox_GMSH(Template_GMSH_Mesh_Constructor):
    def __init__(self, model_name = "Default"):
        super().__init__(model_name)

    @property
    def required_parameters(self)->tuple:
        return ("extrusion_length",
                "elform",
                "nip", 
                "shrf",
                "elsize",
                "id_min",
                "h_level",
                "mid",
                "trigger_depth",
                "trigger_rows",
                "grid",
                "num_layers",
                "cell",
                'n_elements_side'
                )
    

    def __call__(self, var_file:Union[Path,str], save:bool = True, *args, **kwds):
        #super().__call__(*args, **kwds)

        try:
            # Initialize GMSH
            gmsh.initialize(sys.argv,
                            run=False,
                            interruptible=True)
        except:
            print("GMSH couldn't be initialised")

        import numpy as np
        from math import sqrt
        

        # Add the model name
        gmsh.model.add(self.model_name)
        # Load the file given by parameter
        params_dict = self.load_json_file(var_file)

        if not self.check_parameters(params_dict):
            raise ValueError("The parameters to generate the mesh does not" \
            "fulfill the minimal criteria")
        
        # Check from the parameters if the verbosity is set
        verbous = params_dict.pop('gmsh_verbosity', 0)
        
        if verbous == 0:
            gmsh.option.setNumber("General.Terminal", 0)
            gmsh.option.setNumber("General.Verbosity", 0)
        

        # Mesh refinement level
        lv:int = params_dict['h_level']


        # Thickness map
        THICKNESS_MAP = {i + 1: entry["t"] for i, entry in enumerate(params_dict['cell'])}

        # Height of the crashbox
        EXTRUSION_HEIGHT:float = params_dict['extrusion_length']

        # TRIGGER HEIGHT
        TRIGGER_HEIGHT:float = params_dict['elsize']*params_dict['trigger_rows']

        # TRIGGER DEPTH
        depth = params_dict['trigger_depth']


        # POINT SIZE
        lc:float = params_dict['elsize']


        # LATERAL DIVISIONS
        L_D:int = params_dict['n_elements_side']


        grid = [[entry["x"], entry["y"]] for entry in params_dict['grid']]


        #### WRITE THE POINTS ON GMSH

        try:
            gmsh.initialize(sys.argv)
        except:
            print("GMSH couldn't be initialised")


        gmsh.model.add("trial1")
        gmsh.option.setNumber("Mesh.MeshSizeMin", lc)
        gmsh.option.setNumber("Mesh.MeshSizeMax", lc)



        ### ==========================================
        ### Assign the points from 1001 to 1008
        ### ========================================== 

        init_count = 1001
        for ii in range(len(grid)):
            gmsh.model.geo.addPoint(grid[ii][0], grid[ii][1], EXTRUSION_HEIGHT, lc, init_count)
            init_count +=1


        # Append the points to be part of some Physical Group
        gmsh.model.addPhysicalGroup(0,
                                        [a for a in range(1001,1009)],
                                        1,
                                        "Top_Nodes")

        ### ==========================================
        ### Assign the points from 1000001 to 10000008
        ### ========================================== 

        init_count = 1000001
        for ii in range(len(grid)):
            gmsh.model.geo.addPoint(grid[ii][0], grid[ii][1], 0.0, lc, init_count)
            init_count +=1


        # Append the points to be part of some Physical Group
        gmsh.model.addPhysicalGroup(0,
                                        [a for a in range(1000001,1000009)],
                                        2,
                                        "Bottom_Nodes")

        ### ==========================================
        ### Assign the points from 1000009 to 1000016
        ### ========================================== 

        init_count = 1000009
        for ii in range(len(grid)):
            gmsh.model.geo.addPoint(grid[ii][0], grid[ii][1], 84.0, lc, init_count)
            init_count +=1


        ### ==========================================
        ### Assign the points from 1000017 to 1000024
        ### ========================================== 

        init_count = 1000017
        for ii in range(len(grid)):
            gmsh.model.geo.addPoint(grid[ii][0], grid[ii][1], EXTRUSION_HEIGHT-TRIGGER_HEIGHT, lc, init_count)
            init_count +=1

        ### ==========================================
        ### Assign the points from 1000025 to 1000032
        ### ========================================== 

        init_count = 1000025
        for ii in range(len(grid)):
            gmsh.model.geo.addPoint(grid[ii][0], grid[ii][1], 88.0, lc, init_count)
            init_count +=1

        ### ==========================================
        ### Assign the points from 1000033 to 1000040
        ### ========================================== 

        init_count = 1000033
        for ii in range(len(grid)):
            gmsh.model.geo.addPoint(grid[ii][0], grid[ii][1], EXTRUSION_HEIGHT-TRIGGER_HEIGHT-lc, lc, init_count)
            init_count +=1

        ### ==========================================
        ### Assign the points from 1000041 to 1000048
        ### ========================================== 

        init_count = 1000041
        for ii in range(len(grid)):
            gmsh.model.geo.addPoint(grid[ii][0], grid[ii][1], 92.0, lc, init_count)
            init_count +=1
        ### ==========================================
        ### Assign the points from 1000049 to 1000056
        ### ========================================== 

        init_count = 1000049
        for ii in range(len(grid)):
            gmsh.model.geo.addPoint(grid[ii][0], grid[ii][1], 100.0, lc, init_count)
            init_count +=1

        ### ==========================================
        ### Add some intermediate points
        ### ==========================================
        gmsh.model.geo.synchronize()
        init_counter = 1000057
        init_point = 1000025

        for ii in range(4):
            for jj in range(8):

                if jj !=7:
                    p1 = gmsh.model.get_value(0,init_point,[0])
                    p2 = gmsh.model.get_value(0,init_point+1,[0])
                else:
                    p1 = gmsh.model.get_value(0,init_point,[0])
                    p2 = gmsh.model.get_value(0,init_point-7,[0])
                
                # Compute the displacement vector
                vec = 1/L_D * (p2-p1)
                
                if ii in [0,1]:
        
                    

                    # Add the points
                    gmsh.model.geo.add_point(p1[0]+vec[0],
                                            p1[1]+vec[1],
                                            p1[2]+vec[2],
                                            lc,
                                            init_counter)
                    
                    # Add the points
                    gmsh.model.geo.add_point(p2[0]-vec[0],
                                            p2[1]-vec[1],
                                            p2[2]-vec[2],
                                            lc,
                                            init_counter+1)
                
                else:
                    # Compute normal
                    if np.remainder(jj,2)==0:
                        norm_vec = np.linalg.cross(vec,np.asarray([0,0,1]))
                    else:
                        norm_vec = np.linalg.cross(vec,np.asarray([0,0,-1]))
                    norm_vec = norm_vec/np.linalg.norm(norm_vec,2)

                    # Add the points
                    gmsh.model.geo.add_point(p1[0]+vec[0]+norm_vec[0]*depth,
                                            p1[1]+vec[1]+norm_vec[1]*depth,
                                            p1[2]+vec[2]+norm_vec[2]*depth,
                                            lc,
                                            init_counter)
                    
                    # Add the points
                    gmsh.model.geo.add_point(p2[0]-vec[0]+norm_vec[0]*depth,
                                            p2[1]-vec[1]+norm_vec[1]*depth,
                                            p2[2]-vec[2]+norm_vec[2]*depth,
                                            lc,
                                            init_counter+1)
                
                init_point +=1
                init_counter +=2





        del init_counter, init_point

        gmsh.model.geo.synchronize()
        # Append the points to be part of some Physical Group
        gmsh.model.addPhysicalGroup(0,
                                        [a for a in range(1000010,gmsh.model.geo.getMaxTag(0)+1)],
                                        3,
                                        "Infill_Nodes")

        ### ==========================================
        ### Adding lines
        ### ========================================== 

        # Add lines between points at the top
        gmsh.model.geo.addLine(1001, 1002, 1)
        gmsh.model.geo.addLine(1002, 1003, 2)
        gmsh.model.geo.addLine(1003, 1004, 3)
        gmsh.model.geo.addLine(1004, 1005, 4)
        gmsh.model.geo.addLine(1005, 1006, 5)
        gmsh.model.geo.addLine(1006, 1007, 6)
        gmsh.model.geo.addLine(1007, 1008, 7)
        gmsh.model.geo.addLine(1008, 1001, 8)

        # Append the points to be part of some Physical Group
        gmsh.model.addPhysicalGroup(1,
                                        [a for a in range(1,9)],
                                        4,
                                        "Top_Lines")

        # Add lines to the bottom
        gmsh.model.geo.addLine(1000001, 1000002, 9)
        gmsh.model.geo.addLine(1000002, 1000003, 10)
        gmsh.model.geo.addLine(1000003, 1000004, 11)
        gmsh.model.geo.addLine(1000004, 1000005, 12)
        gmsh.model.geo.addLine(1000005, 1000006, 13)
        gmsh.model.geo.addLine(1000006, 1000007, 14)
        gmsh.model.geo.addLine(1000007, 1000008, 15)
        gmsh.model.geo.addLine(1000008, 1000001, 16)

        # Append the points to be part of some Physical Group
        gmsh.model.addPhysicalGroup(1,
                                        [a for a in range(9,17)],
                                        5,
                                        "Bottom_Lines")

        # Add lines at 84.0
        gmsh.model.geo.addLine(1000009, 1000010, 17)
        gmsh.model.geo.addLine(1000010, 1000011, 18)
        gmsh.model.geo.addLine(1000011, 1000012, 19)
        gmsh.model.geo.addLine(1000012, 1000013, 20)
        gmsh.model.geo.addLine(1000013, 1000014, 21)
        gmsh.model.geo.addLine(1000014, 1000015, 22)
        gmsh.model.geo.addLine(1000015, 1000016, 23)
        gmsh.model.geo.addLine(1000016, 1000009, 24)

        # Add lines at 108.0
        gmsh.model.geo.addLine(1000017, 1000018, 25)
        gmsh.model.geo.addLine(1000018, 1000019, 26)
        gmsh.model.geo.addLine(1000019, 1000020, 27)
        gmsh.model.geo.addLine(1000020, 1000021, 28)
        gmsh.model.geo.addLine(1000021, 1000022, 29)
        gmsh.model.geo.addLine(1000022, 1000023, 30)
        gmsh.model.geo.addLine(1000023, 1000024, 31)
        gmsh.model.geo.addLine(1000024, 1000017, 32)

        # Set number of divisions (as transfinite curves) to 13
        for it in range(1,33):
            gmsh.model.geo.mesh.set_transfinite_curve(it,
                                                    L_D+1,
                                                    meshType="Progression",
                                                    coef=1.0)


        # Add the lines from bottom up to 84.0
        gmsh.model.geo.addLine(1000001, 1000009, 33)
        gmsh.model.geo.addLine(1000002, 1000010, 34)
        gmsh.model.geo.addLine(1000003, 1000011, 35)
        gmsh.model.geo.addLine(1000004, 1000012, 36)
        gmsh.model.geo.addLine(1000005, 1000013, 37)
        gmsh.model.geo.addLine(1000006, 1000014, 38)
        gmsh.model.geo.addLine(1000007, 1000015, 39)
        gmsh.model.geo.addLine(1000008, 1000016, 40)

        # Set number of divisions (as transfinite curves) to 22
        for it in range(33,41):
            gmsh.model.geo.mesh.set_transfinite_curve(it,
                                                    22,
                                                    meshType="Progression",
                                                    coef=1.0)

        # Add the lines from top to 108.0
        gmsh.model.geo.addLine(1001, 1000017, 41)
        gmsh.model.geo.addLine(1002, 1000018, 42)
        gmsh.model.geo.addLine(1003, 1000019, 43)
        gmsh.model.geo.addLine(1004, 1000020, 44)
        gmsh.model.geo.addLine(1005, 1000021, 45)
        gmsh.model.geo.addLine(1006, 1000022, 46)
        gmsh.model.geo.addLine(1007, 1000023, 47)
        gmsh.model.geo.addLine(1008, 1000024, 48)

        # Set number of divisions (as transfinite curves) to 4
        for it in range(41,49):
            gmsh.model.geo.mesh.set_transfinite_curve(it,
                                                    4,
                                                    meshType="Progression",
                                                    coef=1.0)
            
        # Add the lines from 84 to 88
        gmsh.model.geo.addLine(1000009, 1000025, 49)
        gmsh.model.geo.addLine(1000010, 1000026, 50)
        gmsh.model.geo.addLine(1000011, 1000027, 51)
        gmsh.model.geo.addLine(1000012, 1000028, 52)
        gmsh.model.geo.addLine(1000013, 1000029, 53)
        gmsh.model.geo.addLine(1000014, 1000030, 54)
        gmsh.model.geo.addLine(1000015, 1000031, 55)
        gmsh.model.geo.addLine(1000016, 1000032, 56)

        # Set number of divisions (as transfinite curves) to 2
        for it in range(49,57):
            gmsh.model.geo.mesh.set_transfinite_curve(it,
                                                    2,
                                                    meshType="Progression",
                                                    coef=1.0)
            

        # Add the lines from 104 to 108
        gmsh.model.geo.addLine(1000017, 1000033, 57)
        gmsh.model.geo.addLine(1000018, 1000034, 58)
        gmsh.model.geo.addLine(1000019, 1000035, 59)
        gmsh.model.geo.addLine(1000020, 1000036, 60)
        gmsh.model.geo.addLine(1000021, 1000037, 61)
        gmsh.model.geo.addLine(1000022, 1000038, 62)
        gmsh.model.geo.addLine(1000023, 1000039, 63)
        gmsh.model.geo.addLine(1000024, 1000040, 64)

        # Set number of divisions (as transfinite curves) to 2
        for it in range(57,65):
            gmsh.model.geo.mesh.set_transfinite_curve(it,
                                                    2,
                                                    meshType="Progression",
                                                    coef=1.0)
            

        # Add the lines from 88 to 92
        gmsh.model.geo.addLine(1000025, 1000041, 65)
        gmsh.model.geo.addLine(1000026, 1000042, 66)
        gmsh.model.geo.addLine(1000027, 1000043, 67)
        gmsh.model.geo.addLine(1000028, 1000044, 68)
        gmsh.model.geo.addLine(1000029, 1000045, 69)
        gmsh.model.geo.addLine(1000030, 1000046, 70)
        gmsh.model.geo.addLine(1000031, 1000047, 71)
        gmsh.model.geo.addLine(1000032, 1000048, 72)

        # Set number of divisions (as transfinite curves) to 2
        for it in range(65,73):
            gmsh.model.geo.mesh.set_transfinite_curve(it,
                                                    2,
                                                    meshType="Progression",
                                                    coef=1.0)
            

        # Add the lines from 100 to 104
        gmsh.model.geo.addLine(1000033, 1000049, 73)
        gmsh.model.geo.addLine(1000034, 1000050, 74)
        gmsh.model.geo.addLine(1000035, 1000051, 75)
        gmsh.model.geo.addLine(1000036, 1000052, 76)
        gmsh.model.geo.addLine(1000037, 1000053, 77)
        gmsh.model.geo.addLine(1000038, 1000054, 78)
        gmsh.model.geo.addLine(1000039, 1000055, 79)
        gmsh.model.geo.addLine(1000040, 1000056, 80)

        # Set number of divisions (as transfinite curves) to 2
        for it in range(73,81):
            gmsh.model.geo.mesh.set_transfinite_curve(it,
                                                    2,
                                                    meshType="Progression",
                                                    coef=1.0)

        # Add the lines from 92 to 100
        # Add the lines from 100 to 104
        gmsh.model.geo.addLine(1000041, 1000049, 81)
        gmsh.model.geo.addLine(1000042, 1000050, 82)
        gmsh.model.geo.addLine(1000043, 1000051, 83)
        gmsh.model.geo.addLine(1000044, 1000052, 84)
        gmsh.model.geo.addLine(1000045, 1000053, 85)
        gmsh.model.geo.addLine(1000046, 1000054, 86)
        gmsh.model.geo.addLine(1000047, 1000055, 87)
        gmsh.model.geo.addLine(1000048, 1000056, 88)

        # Set number of divisions (as transfinite curves) to 3
        for it in range(81,89):
            gmsh.model.geo.mesh.set_transfinite_curve(it,
                                                    3,
                                                    meshType="Progression",
                                                    coef=1.0)

        ### ========================
        # Add the lines at levels 88.0. 92.0, 100,0, 104,0
        ### ========================

        init_counter = 1000057
        init_point = 1000025
        init_curve = 89

        for ii in range(4):
            for jj in range(8):
                if jj !=7:
                    ref1 = init_point
                    ref2 = init_point + 1
                else:
                    ref1 = init_point
                    ref2 = init_point - 7
                
                ref3 = init_counter
                ref4 = init_counter + 1

                # Add curve 1i
                gmsh.model.geo.add_line(ref1,ref3,tag=init_curve)
                gmsh.model.geo.mesh.set_transfinite_curve(init_curve,
                                                2,
                                                meshType="Progression",
                                                coef=1.0)
                
                # Add curve 2i
                gmsh.model.geo.add_line(ref3,ref4,tag=init_curve+1)
                gmsh.model.geo.mesh.set_transfinite_curve(init_curve+1,
                                                L_D-1,
                                                meshType="Progression",
                                                coef=1.0)
                
                # Add curve 3i
                gmsh.model.geo.add_line(ref4,ref2,tag=init_curve+2)
                gmsh.model.geo.mesh.set_transfinite_curve(init_curve+2,
                                                2,
                                                meshType="Progression",
                                                coef=1.0)
                
                init_curve +=3
                init_counter += 2
                init_point +=1


        del init_counter, init_point, init_curve, ref1, ref2, ref3, ref4


        ### ===============================
        # Connect points from 88.0 to 92.0
        ### ===============================
        gmsh.model.geo.synchronize()
        init_curve = gmsh.model.geo.getMaxTag(1) + 1

        # Perform a loop to attach the vertical lines
        init_counter = 1000057

        for jj in range(8):

            ref1 = init_counter
            ref2 = init_counter + 32

            # Connect the points
            # Add curve i
            gmsh.model.geo.add_line(ref1,ref2,tag=init_curve)
            gmsh.model.geo.mesh.set_transfinite_curve(init_curve,
                                            2,
                                            meshType="Progression",
                                            coef=1.0)
            
            gmsh.model.geo.add_line(ref1+1,ref2+1,tag=init_curve+1)
            gmsh.model.geo.mesh.set_transfinite_curve(init_curve+1,
                                            2,
                                            meshType="Progression",
                                            coef=1.0)

                

            init_counter += 2
            init_curve +=2

        ### ===============================
        # Connect points from 92.0 to 100.0
        ### ===============================

        gmsh.model.geo.synchronize()
        init_curve = gmsh.model.geo.getMaxTag(1) + 1
        init_counter = 1000089

        for jj in range(8):

            ref1 = init_counter
            ref2 = init_counter + 16

            # Connect the points
            # Add curve i
            gmsh.model.geo.add_line(ref1,ref2,tag=init_curve)
            gmsh.model.geo.mesh.set_transfinite_curve(init_curve,
                                            3,
                                            meshType="Progression",
                                            coef=1.0)
            
            gmsh.model.geo.add_line(ref1+1,ref2+1,tag=init_curve+1)
            gmsh.model.geo.mesh.set_transfinite_curve(init_curve+1,
                                            3,
                                            meshType="Progression",
                                            coef=1.0)

                

            init_counter += 2
            init_curve +=2

        ### ===============================
        # Connect points from 100.0 to 104.0
        ### ===============================
        gmsh.model.geo.synchronize()
        init_curve = gmsh.model.geo.getMaxTag(1) + 1
        init_counter = 1000105

        for jj in range(8):

            ref1 = init_counter
            ref2 = init_counter -32

            # Connect the points
            # Add curve i
            gmsh.model.geo.add_line(ref1,ref2,tag=init_curve)
            gmsh.model.geo.mesh.set_transfinite_curve(init_curve,
                                            2,
                                            meshType="Progression",
                                            coef=1.0)
            
            gmsh.model.geo.add_line(ref1+1,ref2+1,tag=init_curve+1)
            gmsh.model.geo.mesh.set_transfinite_curve(init_curve+1,
                                            2,
                                            meshType="Progression",
                                            coef=1.0)

                

            init_counter += 2
            init_curve +=2

        del init_curve, init_counter, jj

        gmsh.model.geo.synchronize()
        # Append the points to be part of some Physical Group
        gmsh.model.addPhysicalGroup(1,
                                        [a for a in range(17,gmsh.model.geo.getMaxTag(1)+1)],
                                        6,
                                        "Infill_Lines")

        ### =======================================================
        ### Generate the surfaces
        ### =======================================================

        # Attach bottom to 84
        gmsh.model.geo.addCurveLoop([9, 17, 33, 34], 1, True)
        gmsh.model.geo.addPlaneSurface([1], 1)

        gmsh.model.geo.addCurveLoop([10, 18, 34, 35], 2, True)
        gmsh.model.geo.addPlaneSurface([2], 2)

        gmsh.model.geo.addCurveLoop([11, 19, 35, 36], 3, True)
        gmsh.model.geo.addPlaneSurface([3], 3)

        gmsh.model.geo.addCurveLoop([12, 20, 36, 37], 4, True)
        gmsh.model.geo.addPlaneSurface([4], 4)

        gmsh.model.geo.addCurveLoop([13, 21, 37, 38], 5, True)
        gmsh.model.geo.addPlaneSurface([5], 5)

        gmsh.model.geo.addCurveLoop([14, 22, 38, 39], 6, True)
        gmsh.model.geo.addPlaneSurface([6], 6)

        gmsh.model.geo.addCurveLoop([15, 23, 39, 40], 7, True)
        gmsh.model.geo.addPlaneSurface([7], 7)

        gmsh.model.geo.addCurveLoop([16, 24, 40, 33], 8, True)
        gmsh.model.geo.addPlaneSurface([8], 8)

        # Attach top to 108.0
        gmsh.model.geo.addCurveLoop([1,25,41,42], 9, True)
        gmsh.model.geo.addPlaneSurface([9], 9)

        gmsh.model.geo.addCurveLoop([2,26,42,43], 10, True)
        gmsh.model.geo.addPlaneSurface([10], 10)

        gmsh.model.geo.addCurveLoop([3,27,43,44], 11, True)
        gmsh.model.geo.addPlaneSurface([11], 11)

        gmsh.model.geo.addCurveLoop([4,28,44,45], 12, True)
        gmsh.model.geo.addPlaneSurface([12], 12)

        gmsh.model.geo.addCurveLoop([5,29,45,46], 13, True)
        gmsh.model.geo.addPlaneSurface([13], 13)

        gmsh.model.geo.addCurveLoop([6,30,46,47], 14, True)
        gmsh.model.geo.addPlaneSurface([14], 14)

        gmsh.model.geo.addCurveLoop([7,31,47,48], 15, True)
        gmsh.model.geo.addPlaneSurface([15], 15)

        gmsh.model.geo.addCurveLoop([8,32,48,41], 16, True)
        gmsh.model.geo.addPlaneSurface([16], 16)





        # Attach bottom 84.0 to 88.0
        gmsh.model.geo.addCurveLoop([17,49,50,89,90,91], 17, True)
        gmsh.model.geo.addPlaneSurface([17], 17)

        gmsh.model.geo.addCurveLoop([18,50,51,92,93,94], 18, True)
        gmsh.model.geo.addPlaneSurface([18], 18)

        gmsh.model.geo.addCurveLoop([19,51,52,95,96,97], 19, True)
        gmsh.model.geo.addPlaneSurface([19], 19)

        gmsh.model.geo.addCurveLoop([20,52,53,98,99,100], 20, True)
        gmsh.model.geo.addPlaneSurface([20], 20)

        gmsh.model.geo.addCurveLoop([21,53,54,101,102,103], 21, True)
        gmsh.model.geo.addPlaneSurface([21], 21)

        gmsh.model.geo.addCurveLoop([22,54,55,104,105,106], 22, True)
        gmsh.model.geo.addPlaneSurface([22], 22)

        gmsh.model.geo.addCurveLoop([23,55,56,107,108,109], 23, True)
        gmsh.model.geo.addPlaneSurface([23], 23)

        gmsh.model.geo.addCurveLoop([24,56,49,110,111,112], 24, True)
        gmsh.model.geo.addPlaneSurface([24], 24)

        # Attach 104.0 to 108.0
        gmsh.model.geo.synchronize()
        gmsh.model.geo.addCurveLoop([25,57,58,113,114,115], 25, True)
        gmsh.model.geo.addPlaneSurface([25], 25)

        gmsh.model.geo.addCurveLoop([26,58,59,116,117,118], 26, True)
        gmsh.model.geo.addPlaneSurface([26], 26)

        gmsh.model.geo.addCurveLoop([27,59,60,119,120,121], 27, True)
        gmsh.model.geo.addPlaneSurface([27], 27)

        gmsh.model.geo.addCurveLoop([28,60,61,122,123,124], 28, True)
        gmsh.model.geo.addPlaneSurface([28], 28)

        gmsh.model.geo.addCurveLoop([29,61,62,125,126,127], 29, True)
        gmsh.model.geo.addPlaneSurface([29], 29)

        gmsh.model.geo.addCurveLoop([30,62,63,128,129,130], 30, True)
        gmsh.model.geo.addPlaneSurface([30], 30)

        gmsh.model.geo.addCurveLoop([31,63,64,131,132,133], 31, True)
        gmsh.model.geo.addPlaneSurface([31], 31)

        gmsh.model.geo.addCurveLoop([32,64,57,134,135,136], 32, True)
        gmsh.model.geo.addPlaneSurface([32], 32)


        # Attach bottom 88.0 to 92.0
        init_low = 89
        init_left = 65
        init_right = 185
        init_top = 137
        init_surface = 33

        for ii in range(4):
            
            gmsh.model.geo.addCurveLoop([init_low,init_left,init_right,init_top], init_surface, True)
            gmsh.model.geo.addSurfaceFilling([init_surface],init_surface)

            gmsh.model.geo.addCurveLoop([init_low+1,init_right,init_right+1,init_top+1], init_surface+1, True)
            gmsh.model.geo.addPlaneSurface([init_surface+1],init_surface+1)

            gmsh.model.geo.addCurveLoop([init_low+2,init_right+1,init_left+1,init_top+2], init_surface+2, True)
            gmsh.model.geo.addSurfaceFilling([init_surface+2],init_surface+2)

            gmsh.model.geo.addCurveLoop([init_low+3,init_left+1,init_right+2,init_top+3], init_surface+3, True)
            gmsh.model.geo.addSurfaceFilling([init_surface+3],init_surface+3)

            gmsh.model.geo.addCurveLoop([init_low+4,init_right+2,init_right+3,init_top+4], init_surface+4, True)
            gmsh.model.geo.addPlaneSurface([init_surface+4],init_surface+4)


            if ii !=3: 
                gmsh.model.geo.addCurveLoop([init_low+5,init_right+3,init_left+2,init_top+5], init_surface+5, True)
                gmsh.model.geo.addSurfaceFilling([init_surface+5],init_surface+5)
            
            else:
                gmsh.model.geo.addCurveLoop([init_low+5,init_right+3,init_left-6,init_top+5], init_surface+5, True)
                gmsh.model.geo.addSurfaceFilling([init_surface+5],init_surface+5)


            init_surface += 6
            init_low +=6
            init_top +=6
            init_left+=2
            init_right+=4


        # Attach bottom 92.0 t0 100.0
        init_low = 137
        init_left = 81
        init_right = 201
        init_top = 161
        init_surface = 58

        for ii in range(4):
            
            gmsh.model.geo.addCurveLoop([init_low,init_left,init_right,init_top], init_surface, True)
            gmsh.model.geo.addPlaneSurface([init_surface],init_surface)

            gmsh.model.geo.addCurveLoop([init_low+1,init_right,init_right+1,init_top+1], init_surface+1, True)
            gmsh.model.geo.addPlaneSurface([init_surface+1],init_surface+1)

            gmsh.model.geo.addCurveLoop([init_low+2,init_right+1,init_left+1,init_top+2], init_surface+2, True)
            gmsh.model.geo.addPlaneSurface([init_surface+2],init_surface+2)

            gmsh.model.geo.addCurveLoop([init_low+3,init_left+1,init_right+2,init_top+3], init_surface+3, True)
            gmsh.model.geo.addPlaneSurface([init_surface+3],init_surface+3)

            gmsh.model.geo.addCurveLoop([init_low+4,init_right+2,init_right+3,init_top+4], init_surface+4, True)
            gmsh.model.geo.addPlaneSurface([init_surface+4],init_surface+4)


            if ii !=3: 
                gmsh.model.geo.addCurveLoop([init_low+5,init_right+3,init_left+2,init_top+5], init_surface+5, True)
                gmsh.model.geo.addPlaneSurface([init_surface+5],init_surface+5)
            
            else:
                gmsh.model.geo.addCurveLoop([init_low+5,init_right+3,init_left-6,init_top+5], init_surface+5, True)
                gmsh.model.geo.addPlaneSurface([init_surface+5],init_surface+5)


            init_surface += 6
            init_low +=6
            init_top +=6
            init_left+=2
            init_right+=4


        # Attach bottom 100.0 to 104.0
        init_low = 161
        init_left = 73
        init_right = 217
        init_top = 113
        init_surface = 82

        for ii in range(4):
            
            gmsh.model.geo.addCurveLoop([init_low,init_left,init_right,init_top], init_surface, True)
            gmsh.model.geo.addSurfaceFilling([init_surface],init_surface)

            gmsh.model.geo.addCurveLoop([init_low+1,init_right,init_right+1,init_top+1], init_surface+1, True)
            gmsh.model.geo.addPlaneSurface([init_surface+1],init_surface+1)

            gmsh.model.geo.addCurveLoop([init_low+2,init_right+1,init_left+1,init_top+2], init_surface+2, True)
            gmsh.model.geo.addSurfaceFilling([init_surface+2],init_surface+2)

            gmsh.model.geo.addCurveLoop([init_low+3,init_left+1,init_right+2,init_top+3], init_surface+3, True)
            gmsh.model.geo.addSurfaceFilling([init_surface+3],init_surface+3)

            gmsh.model.geo.addCurveLoop([init_low+4,init_right+2,init_right+3,init_top+4], init_surface+4, True)
            gmsh.model.geo.addPlaneSurface([init_surface+4],init_surface+4)


            if ii !=3: 
                gmsh.model.geo.addCurveLoop([init_low+5,init_right+3,init_left+2,init_top+5], init_surface+5, True)
                gmsh.model.geo.addSurfaceFilling([init_surface+5],init_surface+5)
            
            else:
                gmsh.model.geo.addCurveLoop([init_low+5,init_right+3,init_left-6,init_top+5], init_surface+5, True)
                gmsh.model.geo.addSurfaceFilling([init_surface+5],init_surface+5)


            init_surface += 6
            init_low +=6
            init_top +=6
            init_left+=2
            init_right+=4


        # Set surfaces (so-far) as transfinite
        for ii in range(1,107):

            if (ii < 17 or ii >32): 
                gmsh.model.geo.mesh.setTransfiniteSurface(ii)
                gmsh.model.geo.mesh.setRecombine(2, ii)
            
            elif ii==17:
                gmsh.model.geo.mesh.setTransfiniteSurface(ii,
                                                        "Left",
                                                        [1000009, 1000010, 1000026, 1000025])
                gmsh.model.geo.mesh.setRecombine(2, ii)

            elif ii==18:
                gmsh.model.geo.mesh.setTransfiniteSurface(ii,
                                                        "Left",
                                                        [1000010, 1000011, 1000027, 1000026])
                gmsh.model.geo.mesh.setRecombine(2, ii)
            
            elif ii==19:
                gmsh.model.geo.mesh.setTransfiniteSurface(ii,
                                                        "Left",
                                                        [1000011, 1000012, 1000028, 1000027])
                gmsh.model.geo.mesh.setRecombine(2, ii)

            elif ii==20:
                gmsh.model.geo.mesh.setTransfiniteSurface(ii,
                                                        "Left",
                                                        [1000012, 1000013, 1000029, 1000028])
                gmsh.model.geo.mesh.setRecombine(2, ii)

            elif ii==21:
                gmsh.model.geo.mesh.setTransfiniteSurface(ii,
                                                        "Left",
                                                        [1000013, 1000014, 1000030, 1000029])
                gmsh.model.geo.mesh.setRecombine(2, ii)

            elif ii==22:
                gmsh.model.geo.mesh.setTransfiniteSurface(ii,
                                                        "Left",
                                                        [1000014, 1000015, 1000031, 1000030])
                gmsh.model.geo.mesh.setRecombine(2, ii)
            
            elif ii==23:
                gmsh.model.geo.mesh.setTransfiniteSurface(ii,
                                                        "Left",
                                                        [1000015, 1000016, 1000032, 1000031])
                gmsh.model.geo.mesh.setRecombine(2, ii)
            elif ii==24:
                gmsh.model.geo.mesh.setTransfiniteSurface(ii,
                                                        "Left",
                                                        [1000016, 1000009, 1000025, 1000032])
                gmsh.model.geo.mesh.setRecombine(2, ii)
            elif ii==25:
                gmsh.model.geo.mesh.setTransfiniteSurface(ii,
                                                        "Left",
                                                        [1000033, 1000034, 1000018, 1000017])
                gmsh.model.geo.mesh.setRecombine(2, ii)
            elif ii==26:
                gmsh.model.geo.mesh.setTransfiniteSurface(ii,
                                                        "Left",
                                                        [1000034, 1000035, 1000019, 1000018])
                gmsh.model.geo.mesh.setRecombine(2, ii)
            elif ii==27:
                gmsh.model.geo.mesh.setTransfiniteSurface(ii,
                                                        "Left",
                                                        [1000035, 1000036, 1000020, 1000019])
                gmsh.model.geo.mesh.setRecombine(2, ii)
            elif ii==28:
                gmsh.model.geo.mesh.setTransfiniteSurface(ii,
                                                        "Left",
                                                        [1000036, 1000037, 1000021, 1000020])
                gmsh.model.geo.mesh.setRecombine(2, ii)
            elif ii==29:
                gmsh.model.geo.mesh.setTransfiniteSurface(ii,
                                                        "Left",
                                                        [1000037, 1000038, 1000022, 1000021])
                gmsh.model.geo.mesh.setRecombine(2, ii)
            elif ii==30:
                gmsh.model.geo.mesh.setTransfiniteSurface(ii,
                                                        "Left",
                                                        [1000038, 1000039, 1000023, 1000022])
                gmsh.model.geo.mesh.setRecombine(2, ii)
            elif ii==31:
                gmsh.model.geo.mesh.setTransfiniteSurface(ii,
                                                        "Left",
                                                        [1000039, 1000040, 1000024, 1000023])
                gmsh.model.geo.mesh.setRecombine(2, ii)

            elif ii==32:
                gmsh.model.geo.mesh.setTransfiniteSurface(ii,
                                                        "Left",
                                                        [1000040, 1000024, 1000017, 1000033])
                gmsh.model.geo.mesh.setRecombine(2, ii)
                


        # Remove duplicate elements
        gmsh.model.geo.removeAllDuplicates()

        gmsh.model.geo.synchronize()
        # Append the points to be part of some Physical Group
        gmsh.model.addPhysicalGroup(2,
                                        [a for a in range(1,gmsh.model.geo.getMaxTag(2)+1)],
                                        7,
                                        "Infill_Surfs")

        gmsh.model.geo.synchronize()


        # # Append the points to be part of some Physical Group
        # gmsh.model.addPhysicalGroup(2,
        #                                 [a for a in range(1,gmsh.model.geo.getMaxTag(2)+1)],
        #                                 8,
        #                                 "Infill_Surfs_2")
        #gmsh.model.addPhysicalGroup(1, [1, 2, 4], 5)
        #gmsh.model.addPhysicalGroup(2, [1], name="My surface")

        # We can then generate a 2D mesh...
        gmsh.model.mesh.setAlgorithm(2,1,4)

        gmsh.option.setNumber("Mesh.Smoothing", 100)


        ### =========================================
        ### Generate the mesh
        ### =========================================
        gmsh.model.mesh.generate(dim = 2)
        gmsh.model.mesh.removeDuplicateNodes()
        gmsh.model.mesh.removeDuplicateElements()

        gmsh.model.geo.synchronize()



        ### -------------------------------------- 
        ### Refine depending on the number of dimensions
        ### --------------------------------------
        for _ in range(lv-1):
            gmsh.model.mesh.refine()


        gmsh.plugin.setNumber("SimplePartition", "NumSlicesX", 1)
        gmsh.plugin.setNumber("SimplePartition", "NumSlicesY", 1)
        gmsh.plugin.setNumber("SimplePartition", "NumSlicesZ", params_dict['num_layers'])
        gmsh.plugin.run("SimplePartition")

        # Get the bounding box and then split
        xx1 =  gmsh.model.getNumberOfPartitions()

        ### =========================================
        ### Get the elements within a bounding box
        ### =========================================
        element_list:list = []

        #global_bounding_box = gmsh.model.getBoundingBox()
        for aa in range(xx1):
            # Initialise a sublist
            sub_element_list = []
            element_list.append(sub_element_list)



        # Loop all over the partitioned entities
        entities = gmsh.model.getEntities(dim=2)
        for e in entities:
            # Extract the partition ID
            partitions = gmsh.model.getPartitions(e[0], e[1])

            partition_element_list= gmsh.model.mesh.getElements(e[0], e[1])[1]

            if len(partition_element_list):
                for new_elem in partition_element_list:
                    element_list[partitions[0]-1].extend(new_elem)

            # if e[0]==2:
            #     if len(partitions):
            #         print("Entity " + str(e) + " of type " +
            #                 gmsh.model.getType(e[0], e[1]))
            #         print(" - Partition(s): " + str(partitions))
            #         print(" - Parent: " + str(gmsh.model.getParent(e[0], e[1])))
            #         print(" - Boundary: " + str(gmsh.model.getBoundary([e])))

        # UNPARTITION
        gmsh.model.mesh.unpartition()
        gmsh.model.geo.synchronize()

        # Get the Ids of the Nodes linked to the points at the top
        top_nodesIDs = gmsh.model.mesh.getNodesForPhysicalGroup(0,1)[0]

        # Get the rest of nodes
        rest_nodes = np.setdiff1d(gmsh.model.mesh.getNodes(dim=-1,returnParametricCoord=False)[0],top_nodesIDs)

        new_IDs_rest = np.arange(start=1000001,
                            stop=1000001 + rest_nodes.size,
                            step=1)
        new_IDs_top = np.arange(start=1001,
                            stop=1001+top_nodesIDs.size,
                            step=1)

        # Replace the IDs of these nodes
        all_nodes = gmsh.model.mesh.getNodes(dim=-1,returnParametricCoord=False)[0]
        new_set_IDs = np.hstack((new_IDs_top,new_IDs_rest))
        gmsh.model.mesh.renumberNodes(all_nodes, new_set_IDs)


        ### Delete some variables
        del all_nodes, new_IDs_top, new_IDs_rest, rest_nodes, top_nodesIDs


        # Synchronize
        gmsh.model.geo.synchronize()

        bottom_nodes_ID = gmsh.model.mesh.getNodesForPhysicalGroup(1,5)[0]
        top_nodes_ID = gmsh.model.mesh.getNodesForPhysicalGroup(1,4)[0]

        # Set the part IDS as a vector
        set_pids = np.arange(101,101+xx1,1,dtype=int).tolist()


        with open("mesh.k", "w") as f:
            
            ### NOTE: Write the header of keyword
            f.write("*KEYWORD\n")

            

            # Node sets
            l = '*SET_NODE_LIST\n'
            l += '${0:->9}\n'.format('sid')
            l += '{0:>10}\n'.format(101)
            l += '${0:->9}{1:->10}{2:->10}{3:->10}'.format('nid1',
                                                            'nid2',
                                                            'nid3',
                                                            'nid4')
            
            l += '{0:->10}{1:->10}{2:->10}{3:->10}\n'.format('nid5',
                                                                'nid6',
                                                                'nid7',
                                                                'nid8')
            bottom_nodes_ID.sort()
            for i in range(0, len(bottom_nodes_ID), 8):
                chunk = bottom_nodes_ID[i:i+8]
                for n in chunk:
                    l += '{0:>10}'.format(n)
                l += '\n'
            
            f.write(l)

            ### ++++++++++++++++++++++++++++++++++++
            ### Write the top nodes as part of a group
            ### ++++++++++++++++++++++++++++++++++++

            


            # Node sets
            l = '*SET_NODE_LIST\n'
            l += '${0:->9}\n'.format('sid')
            l += '{0:>10}\n'.format(102)
            l += '${0:->9}{1:->10}{2:->10}{3:->10}'.format('nid1',
                                                            'nid2',
                                                            'nid3',
                                                            'nid4')
            l += '{0:->10}{1:->10}{2:->10}{3:->10}\n'.format('nid5',
                                                                'nid6',
                                                                'nid7',
                                                                'nid8')
            top_nodes_ID.sort()
            for i in range(0, len(top_nodes_ID), 8):
                chunk = top_nodes_ID[i:i+8]
                for n in chunk:
                    l += '{0:>10}'.format(n)
                l += '\n'
            
            f.write(l)
            
            # Nodes
            node_tags, node_coords, _ = gmsh.model.mesh.getNodes()
            f.write("*NODE\n")
            f.write('${0:->7}{1:->19}{2:->19}{3:->19}\n'.format('nid', 'x', 'y', 'z'))
            for i in range(len(node_tags)):
                f.write('{0:>8},{1:18.10f},{2:18.10f},{3:18.10f}\n'.format(
                                                                    int(node_tags[i]),  #nid
                                                                        node_coords[3*i], #x
                                                                        node_coords[3*i+1], #y 
                                                                        node_coords[3*i+2] #z
                                                                    )
                                                                        
                        )
            
            ### ++++++++++++++++++++++++++++++++++
            ### Write the Parts
            ### ++++++++++++++++++++++++++++++++++
            # NOTE: REPEATING THE PROCEDURE FOR EACH REPETITION
            for ii,itt in enumerate(params_dict['cell']):
                f.write('*PART\n')
                f.write(f'Part_{ii+101}\n')
                f.write('')

                f.write('${0:->9}{1:->10}{2:->10}\n'.format('pid', 'secid', 'mid'))
                f.write('{0:>10}{1:>10}{2:>10}\n'.format(itt['cid'], itt['cid'], params_dict['mid']))

                f.write('*SECTION_SHELL\n')
                f.write('${0:->9}{1:->10}{2:->10}{3:->10}\n'.format('secid',
                                                                'elform',
                                                                'shrf',
                                                                'nip'))
                f.write('{0:>10}{1:>10}{2:>10.5f}{3:>10}\n'.format(itt['cid'], 
                                                                params_dict['elform'], #Element formulation (default 2 -> Tsai-Belytschko)
                                                                params_dict['shrf'], #shrf #Shearing factor as in Reissner-Mindlin plate
                                                                params_dict['nip']))  #nip (Integration Points)
                
                f.write('${0:->14}{1:->15}{2:->15}{3:->15}\n'.format('t1', 't2', 't3', 't4'))
                f.write('{0:>15.10f}{1:>15.10f}{2:>15.10f}{3:>15.10f}\n'.format(
                                                                                THICKNESS_MAP[ii+1],
                                                                                THICKNESS_MAP[ii+1],
                                                                                THICKNESS_MAP[ii+1],
                                                                                THICKNESS_MAP[ii+1]
                                                                                ))


                
            # Elements grouped by partition
            f.write("*ELEMENT_SHELL\n")
            f.write('${0:->7}{1:->8}'.format('eid', 'pid'))
            f.write('{0:->8}{1:->8}{2:->8}{3:->8}\n'.format('n1', 'n2', 'n3', 'n4'))

            for ii in range(xx1):
                for eid in element_list[ii]:
                    f.write('{0:>8}{1:>8}'.format(eid, ii+101))
                    # Get the nodes of each element
                    n_eid = gmsh.model.mesh.getElement(eid)[1]
                    f.write('{0:>8}{1:>8}{2:>8}{3:>8}\n'.format(n_eid[0], 
                                                                n_eid[1], 
                                                                n_eid[2], 
                                                                n_eid[3]))

            f.write("*END\n")


        # UNPARTITION
        #gmsh.model.mesh.unpartition()

        # ... and save it to disk
        #gmsh.write("t1.geo_unrolled")
        gmsh.write("mesh1_star_box.vtk")




        self._end_gmsh()