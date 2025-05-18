from src.sob.lib.gmsh_base_meshes import Template_GMSH_Mesh_Constructor
import gmsh
import json
from pathlib import Path
from typing import Union, List
import os,sys


class Crashtube_GMSH(Template_GMSH_Mesh_Constructor):
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
                "grid",
                "cell",
                'n_elements_side',
                'trigger_dict_1',
                'trigger_dict_2'
                )
    

    def __call__(self, var_file:Union[Path,str], save:bool = True, *args, **kwds):

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
            raise ValueError("The parameters to generate the mesh does not " \
            "fulfill the minimal criteria")
        
        # Check from the parameters if the verbosity is set
        verbous = params_dict.pop('gmsh_verbosity', 0)
        
        if verbous == 0:
            gmsh.option.setNumber("General.Terminal", 0)
            gmsh.option.setNumber("General.Verbosity", 0)

        # MESH REFINEMENT LEVEL
        lv:int = params_dict['h_level']

        # POINT SIZE
        lc:float = params_dict['elsize']

        THICKNESS:float = params_dict['cell'][0]["t"]

        # Number of elements divisions on the side
        NUM_ELEM_SIDE:int = params_dict['n_elements_side']


        HEIGHT = params_dict['extrusion_length']



        ### ==========================================
        ### TRIGGER CHARACTERISTICS
        ### ==========================================

        TRIGGERS_FACE_1:dict = params_dict['trigger_dict_1']
    

        TRIGGERS_FACE_2:dict = params_dict['trigger_dict_2']
        

        # INITIAL VERTICES
        #init_grid = [[entry["x"], entry["y"]] for entry in params_dict['grid']]
        # BOUNDING BOX 
        WIDTH = max([entry["y"] for entry in params_dict['grid']]) - min([entry["y"] for entry in params_dict['grid']])
        LENGTH = max([entry["x"] for entry in params_dict['grid']]) - min([entry["x"] for entry in params_dict['grid']])

        elem_width_side_1 = LENGTH/NUM_ELEM_SIDE
        elem_width_side_2 = WIDTH/NUM_ELEM_SIDE

        base_verts:list = [ [LENGTH/2.0,WIDTH/2.0],
                            [LENGTH/2.0-elem_width_side_1,WIDTH/2.0],
                            [-LENGTH/2.0+elem_width_side_1,WIDTH/2.0],
                            [-LENGTH/2.0,WIDTH/2.0],
                            [-LENGTH/2.0,WIDTH/2.0-elem_width_side_2],
                            [-LENGTH/2.0,-WIDTH/2.0+elem_width_side_2],
                            [-LENGTH/2.0,-WIDTH/2.0],
                            [-LENGTH/2.0+elem_width_side_1,-WIDTH/2.0],
                            [LENGTH/2.0-elem_width_side_1,-WIDTH/2.0],
                            [LENGTH/2.0,-WIDTH/2.0],
                            [LENGTH/2.0,-WIDTH/2.0+elem_width_side_2],
                            [LENGTH/2.0,WIDTH/2.0-elem_width_side_2]]
        



        def get_height_points_profile(trigger_dict_1:dict,
                                    trigger_dict_2:dict)->List[List[float]]:
            r"""
            This is a function to generate profiles of points given the triggers positions
            """

            def determine_trigger_profile(yy:float,
                                        hh:float)->List[float]:
                r"""
                Determine the profile of the trigger

                Args:
                --------------------
                yy: (`float`): The centroid of the trigger
                hh: (`float`): The height of the trigger 
                """


                # Compute the offset from centroid to the starting point
                # The middle of the "gap" between a+lc and a+lc+hh should be yy
                offset = lc + hh / 2

                a = yy - offset

                # Snap 'a' to nearest multiple of lc
                a = np.round(a / lc) * lc

                profile = [a,
                        a + lc,
                        a + lc + hh,
                        a + 2 * lc + hh]

                return [int(x) for x in profile]
            
            
            list_face_1 = []
            list_face_2 = []
            # Perform the function from Face 1 (and 3... for now)
            for key,trigg in trigger_dict_1.items():
                y_i = trigg['Y']
                h_i = trigg['HEIGHT']

                list_face_1.append(determine_trigger_profile(y_i,h_i))
            

            for key,trigg in trigger_dict_2.items():
                y_i = trigg['Y']
                h_i = trigg['HEIGHT']

                list_face_2.append(determine_trigger_profile(y_i,h_i))
            
            return list_face_1,list_face_2
        



        #### WRITE THE POINTS ON GMSH

        try:
            gmsh.initialize(sys.argv)
        except:
            print("GMSH couldn't be initialised")


        gmsh.model.add("trial_crash_tube")
        gmsh.option.setNumber("Mesh.MeshSizeMin", lc)
        gmsh.option.setNumber("Mesh.MeshSizeMax", lc)

        ### ========================================
        ### ADD POINTS AND CURVES
        ### ========================================

        # Top Points Layer
        topPointsID = []
        bottomPointsID = []

        # A structure with all the points per layer
        layersPointsID = []

        # Horizontal lines of bottom and top layers
        topLinesID = []
        bottomLinesID = []

        # Generate the points from the top and bottom

        for ii,point in enumerate(base_verts):
            topPointsID.append(gmsh.model.geo.addPoint(point[0],point[1],HEIGHT,lc))


        for ii,point in enumerate(base_verts):
            bottomPointsID.append(gmsh.model.geo.addPoint(point[0],point[1],0.0,lc))

        # Append the points to respective physical groups
        topPoints_phys_tag = gmsh.model.geo.addPhysicalGroup(dim=0,
                                        tags=topPointsID,
                                        name="Top Points"
                                        )

        topPoints_corners_phys_tag = gmsh.model.geo.addPhysicalGroup(dim=0,
                                        tags=[topPointsID[0],
                                            topPointsID[3],
                                            topPointsID[6],
                                            topPointsID[9],
                                            ],
                                        name="Top Points Corners"
                                        )

        bottomPoints_phys_tag = gmsh.model.geo.addPhysicalGroup(dim=0,
                                        tags=bottomPointsID,
                                        name="Bottom Points"
                                        )

        # Set the points profile

        list_1, list_2 = get_height_points_profile(TRIGGERS_FACE_1,TRIGGERS_FACE_2)

        # Get a list of unique points to generate at the profile
        unique_points = np.unique(np.vstack((np.asarray(list_1).ravel(),
                                            np.asarray(list_2).ravel()))).tolist()


        # Add list of unique points in the geometry
        for iPoint in unique_points:
            # Initialize a list to hold the Points IDs of each layer
            iList = []
            for vert in base_verts:
                iList.append(gmsh.model.geo.addPoint(vert[0],vert[1],iPoint,lc))
            
            # Append the iList to the main `layersPointsID`
            layersPointsID.append(iList)


        layersPointsID_phys_tag = gmsh.model.geo.addPhysicalGroup(dim=0,
                                        tags=np.asarray(layersPointsID).ravel().tolist(), # flatten the array
                                        name="Layer Points"
                                        )


        ### ============================================================
        ### WRITE THE LINES
        ### ============================================================

        # Generate lines at bottom and top 
        # top
        for iid in range(len(topPointsID)):
            if iid != len(topPointsID)-1:
                topLinesID.append(gmsh.model.geo.addLine(topPointsID[iid],topPointsID[iid+1]))
            else:
                topLinesID.append(gmsh.model.geo.addLine(topPointsID[iid],topPointsID[0]))

        # bottom
        for iid in range(len(bottomPointsID)):
            if iid != len(bottomPointsID)-1:
                bottomLinesID.append(gmsh.model.geo.addLine(bottomPointsID[iid],bottomPointsID[iid+1]))
            else:
                bottomLinesID.append(gmsh.model.geo.addLine(bottomPointsID[iid],bottomPointsID[0]))




        topLines_phys_tag = gmsh.model.geo.addPhysicalGroup(dim=1,
                                        tags=topLinesID,
                                        name="Top Lines"
                                        )


        bottomLines_phys_tag = gmsh.model.geo.addPhysicalGroup(dim=1,
                                        tags=bottomLinesID,
                                        name="Bottom Lines"
                                        )


        # Append the lines among all the layers
        vertical_lines_layers_ids = []
        horizontal_lines_layers_ids = []

        for ii in range(len(layersPointsID)+1):
            iList = []
            if ii==0:
                for jj in range(len(bottomPointsID)):
                    iList.append(gmsh.model.geo.addLine(bottomPointsID[jj],layersPointsID[ii][jj]))
            elif ii == len(layersPointsID):
                for jj in range(len(topPointsID)):
                    iList.append(gmsh.model.geo.addLine(layersPointsID[ii-1][jj],topPointsID[jj]))
            else:
                for jj in range(len(layersPointsID[ii-1])):
                    iList.append(gmsh.model.geo.addLine(layersPointsID[ii-1][jj],layersPointsID[ii][jj]))
            
            vertical_lines_layers_ids.append(iList)


        # Horizontal lines among the layers
        for ii in range(len(layersPointsID)):
            iList = []
            for jj in range(len(layersPointsID[ii])):
                if jj != len(layersPointsID[ii])-1:
                    iList.append(gmsh.model.geo.addLine(layersPointsID[ii][jj],layersPointsID[ii][jj+1]))
                else:
                    iList.append(gmsh.model.geo.addLine(layersPointsID[ii][jj],layersPointsID[ii][0]))
            
            horizontal_lines_layers_ids.append(iList)

        # Add lines (horizontal and vertical as Physical Groups)

        horizontalLines_phys_tag = gmsh.model.geo.addPhysicalGroup(dim=1,
                                        tags=np.asarray(horizontal_lines_layers_ids).ravel().tolist(),
                                        name="Horizontal Lines"
                                        )

        verticalLines_phys_tag = gmsh.model.geo.addPhysicalGroup(dim=1,
                                        tags=np.asarray(vertical_lines_layers_ids).ravel().tolist(),
                                        name="Vertical Lines"
                                        )
        ### ================================================================
        ### MAKE THE TRIGGERS
        ### ================================================================

        gmsh.model.geo.synchronize()

        # Make the triggers of faces 1 (and 3)
        for ii,ll_ in enumerate(list_1):
            mid1, mid2 = ll_[1] , ll_[2]

            # Extract the depth of the trigger
            iDepth = TRIGGERS_FACE_1[str(ii+1)]["DEPTH"]

            # Match the heights to the vector of internal points
            idxs = np.where(np.isclose(mid1, unique_points))[0], np.where(np.isclose(mid2, unique_points))[0]

            # Get the exact points to change
            mod_points_f1 = []
            mod_points_f3 = []

            for jj in np.arange(int(idxs[0]),int(idxs[1]+1),1,dtype=int):
                #FACE 1
                mod_points_f1.extend(layersPointsID[jj][1:3])
                mod_points_f3.extend(layersPointsID[jj][7:9])


            mod_points_list_f1 = [(0,pp) for pp in mod_points_f1 ]
            mod_points_list_f3 = [(0,pp) for pp in mod_points_f3 ]


            # Generate the displacement on Face 1
            gmsh.model.geo.translate(mod_points_list_f1,
                                    0,
                                    iDepth,
                                    0)
            
            # Generate the displacement on Face 3
            gmsh.model.geo.translate(mod_points_list_f3,
                                    0,
                                    -iDepth,
                                    0)

        # Make the triggers of faces 2 (and 4)
        for ii,ll_ in enumerate(list_2):
            mid1, mid2 = ll_[1] , ll_[2]

            # Extract the depth of the trigger
            iDepth = TRIGGERS_FACE_2[str(ii+1)]["DEPTH"]

            # Match the heights to the vector of internal points
            idxs = np.where(np.isclose(mid1, unique_points))[0], np.where(np.isclose(mid2, unique_points))[0]

            # Get the exact points to change
            mod_points_f2 = []
            mod_points_f4 = []

            for jj in np.arange(int(idxs[0]),int(idxs[1]+1),1,dtype=int):
                #FACE 1
                mod_points_f2.extend(layersPointsID[jj][4:6])
                mod_points_f4.extend(layersPointsID[jj][10:12])


            mod_points_list_f2 = [(0,pp) for pp in mod_points_f2 ]
            mod_points_list_f4 = [(0,pp) for pp in mod_points_f4 ]


            # Generate the displacement on Face 1
            gmsh.model.geo.translate(mod_points_list_f2,
                                    -iDepth,
                                    0,
                                    0)
            
            # Generate the displacement on Face 3
            gmsh.model.geo.translate(mod_points_list_f4,
                                    iDepth,
                                    0,
                                    0)


        ### ================================================================
        ### GENERATE THE SURFACES 
        ### ================================================================

        # Create the lateral surfaces
        lateral_surfaces = []
        for ii in range(len(vertical_lines_layers_ids)):
            iSurface_layer = []
            for jj in range(len(vertical_lines_layers_ids[ii])):

                if jj == len(bottomLinesID)-1:
                    limit_vert = vertical_lines_layers_ids[ii][0]
                else:
                    limit_vert = vertical_lines_layers_ids[ii][jj+1]

                if ii == 0:
                    curve_loop_lateral = gmsh.model.geo.addCurveLoop(
                        [
                            bottomLinesID[jj],
                            vertical_lines_layers_ids[ii][jj],
                            limit_vert,
                            horizontal_lines_layers_ids[ii][jj],
                        ],
                        reorient=True
                    )
                elif ii == len(vertical_lines_layers_ids)-1:
                    curve_loop_lateral = gmsh.model.geo.addCurveLoop(
                        [
                            horizontal_lines_layers_ids[ii-1][jj],
                            vertical_lines_layers_ids[ii][jj],
                            limit_vert,
                            topLinesID[jj],
                        ],
                        reorient=True
                    )
                else:
                    curve_loop_lateral = gmsh.model.geo.addCurveLoop(
                        [
                            horizontal_lines_layers_ids[ii-1][jj],
                            vertical_lines_layers_ids[ii][jj],
                            limit_vert,
                            horizontal_lines_layers_ids[ii][jj],
                        ],
                        reorient=True
                    )

                iSurface_layer.append(gmsh.model.geo.addSurfaceFilling([curve_loop_lateral]))
            
            lateral_surfaces.append(iSurface_layer)

        gmsh.model.geo.synchronize()


        # Set the surfaces as physical groups
        lateral_surfaces_phys_tag = gmsh.model.geo.addPhysicalGroup(dim=2,
                                        tags=np.asarray(lateral_surfaces).ravel().tolist(),
                                        name="Lateral Surfaces"
                                        )


        gmsh.model.geo.synchronize()
        # Set the lines as transfinite curves
        for iid in range(0,len(bottomLinesID),3):
            gmsh.model.geo.mesh.setTransfiniteCurve(bottomLinesID[iid],
                                                nPoints=2,
                                                meshType="Progression",
                                                coef=1.0)
            
            gmsh.model.geo.mesh.setTransfiniteCurve(bottomLinesID[iid+1],
                                                nPoints=11,
                                                meshType="Progression",
                                                coef=1.0)
            
            gmsh.model.geo.mesh.setTransfiniteCurve(bottomLinesID[iid+2],
                                                nPoints=2,
                                                meshType="Progression",
                                                coef=1.0)

        gmsh.model.geo.synchronize()

        # Set the lines as transfinite curves
        for iid in range(0,len(topLinesID),3):
            gmsh.model.geo.mesh.setTransfiniteCurve(topLinesID[iid],
                                                nPoints=2,
                                                meshType="Progression",
                                                coef=1.0)
            
            gmsh.model.geo.mesh.setTransfiniteCurve(topLinesID[iid+1],
                                                nPoints=11,
                                                meshType="Progression",
                                                coef=1.0)
            
            gmsh.model.geo.mesh.setTransfiniteCurve(topLinesID[iid+2],
                                                nPoints=2,
                                                meshType="Progression",
                                                coef=1.0)


        gmsh.model.geo.synchronize()

        # Set the horizontal lines as transfinite curves
        for jj in range(len(horizontal_lines_layers_ids)):

            for iid in range(0,len(horizontal_lines_layers_ids[jj]),3):
                gmsh.model.geo.mesh.setTransfiniteCurve(horizontal_lines_layers_ids[jj][iid],
                                                    nPoints=2,
                                                    meshType="Progression",
                                                    coef=1.0)
                
                gmsh.model.geo.mesh.setTransfiniteCurve(horizontal_lines_layers_ids[jj][iid+1],
                                                    nPoints=11,
                                                    meshType="Progression",
                                                    coef=1.0)
                
                gmsh.model.geo.mesh.setTransfiniteCurve(horizontal_lines_layers_ids[jj][iid+2],
                                                    nPoints=2,
                                                    meshType="Progression",
                                                    coef=1.0)
                
        gmsh.model.geo.synchronize()

        # Set the vertical lines as transfinite curves
        for jj in range(len(vertical_lines_layers_ids)):
            # Get the distance
            if jj==0:
                dist = unique_points[0]
            elif jj == len(vertical_lines_layers_ids)-1:
                dist = HEIGHT - unique_points[-1]
            else:
                dist = unique_points[jj]-unique_points[jj-1]
            
            nn = int(dist/4.0 +1)
            for iid in range(len(vertical_lines_layers_ids[jj])):
                gmsh.model.geo.mesh.setTransfiniteCurve(tag=vertical_lines_layers_ids[jj][iid],
                                                        nPoints=nn,
                                                        meshType="Progression",
                                                        coef=1.00)

        gmsh.model.geo.synchronize()
        # Set these surfaces as transfinite
        for ii in range(len(lateral_surfaces)):
            for jj in range(len(lateral_surfaces[ii])):
                actualSurfId = lateral_surfaces[ii][jj]

                if jj == len(lateral_surfaces[ii])-1:
                    limit_point_idx = 0
                else:
                    limit_point_idx = jj+1

                # Extract the points

                if ii == 0:
                    pointList_ = [bottomPointsID[jj],bottomPointsID[limit_point_idx],
                                layersPointsID[ii][jj],layersPointsID[ii][limit_point_idx]]
                    
                elif ii == len(lateral_surfaces)-1:
                    pointList_ = [topPointsID[jj],topPointsID[limit_point_idx],
                                layersPointsID[ii-1][jj],layersPointsID[ii-1][limit_point_idx]]
                else:
                    pointList_ = [layersPointsID[ii][jj],layersPointsID[ii][limit_point_idx],
                                layersPointsID[ii-1][jj],layersPointsID[ii-1][limit_point_idx]]

                gmsh.model.geo.mesh.setTransfiniteSurface(actualSurfId,
                                                    cornerTags=pointList_,
                                                    arrangement="Left")
                gmsh.model.geo.mesh.setRecombine(2, actualSurfId)


        # Synchronize
        gmsh.model.geo.synchronize()


        # Remove duplicate elements
        gmsh.model.geo.removeAllDuplicates()
        gmsh.model.geo.synchronize()
        gmsh.option.setNumber("Mesh.Smoothing", 100)


        ### =========================================
        ### Generate the mesh
        ### =========================================
        #gmsh.model.mesh.generate(dim = 1)
        gmsh.model.mesh.generate(dim = 2)
        gmsh.model.mesh.removeDuplicateNodes()
        gmsh.model.mesh.removeDuplicateElements()

        gmsh.model.geo.synchronize()



        ### -------------------------------------- 
        ### Refine depending on the number of dimensions
        ### --------------------------------------
        for _ in range(lv-1):
            gmsh.model.mesh.refine()

        # Get the Ids of the Nodes linked to the points at the top
        top_nodesIDs = gmsh.model.mesh.getNodesForPhysicalGroup(0,topPoints_corners_phys_tag)[0]

        # Get the rest of nodes
        rest_nodes = np.setdiff1d(gmsh.model.mesh.getNodes(dim=-1,
                                                        returnParametricCoord=False)[0],
                                                        top_nodesIDs)

        new_IDs_rest = np.arange(
                            start=1000001,
                            stop=1000001 + rest_nodes.size,
                            step=1
                            )

        new_IDs_top = np.arange(start=1001,
                            stop=1001+top_nodesIDs.size,
                            step=1)

        # Replace the IDs of these nodes
        all_nodes = np.hstack((top_nodesIDs.ravel(),rest_nodes.ravel()))
        new_set_IDs = np.hstack(
                                (
                                new_IDs_top,
                                new_IDs_rest
                                )
                                )

        gmsh.model.mesh.renumberNodes(all_nodes, new_set_IDs)


        # ### Delete some variables
        # del all_nodes, new_IDs_top, new_IDs_rest, rest_nodes, top_nodesIDs


        # Synchronize
        gmsh.model.geo.synchronize()

        bottom_nodes_ID = gmsh.model.mesh.getNodesForPhysicalGroup(1,bottomLines_phys_tag)[0]
        top_nodes_ID = gmsh.model.mesh.getNodesForPhysicalGroup(1,topLines_phys_tag)[0]

        # Get the list of all elements
        element_list = np.asarray(
                                    gmsh.model.mesh.getElements(dim=2, tag=-1)[1]
                                ,dtype=int).tolist()

        # Set the part IDS as a vector
        set_pids = [params_dict['cell'][0]['cid']]
        init_sid = params_dict['cell'][0]['cid']


        with open("mesh.k", "w") as f:
            
            ### NOTE: Write the header of keyword
            f.write("*KEYWORD\n")

            ### ++++++++++++++++++++++++++++++++++++++
            ### Write the Part list (1 per partition)
            ### ++++++++++++++++++++++++++++++++++++++
            f.write("*SET_PART_LIST\n")
            f.write('${0:->9}\n'.format('sid'))
            f.write('{0:>10}\n'.format(101))
            f.write('${0:->9}{1:->10}{2:->10}{3:->10}'.format('pid1',
                                                        'pid2',
                                                        'pid3',
                                                        'pid4'))
            f.write('{0:->10}{1:->10}{2:->10}{3:->10}\n'.format('pid5',
                                                            'pid6',
                                                            'pid7',
                                                            'pid8'))

            for i in range(0, len(set_pids), 8):
                chunk = set_pids[i:i+8]
                for p in chunk:
                    f.write('{0:>10}'.format(p))
                f.write('\n')

            # #     r"""
            # #         l += '${0:->9}{1:->10}{2:->10}{3:->10}'.format('pid1',
            # #                                                'pid2',
            # #                                                'pid3',
            # #                                                'pid4')
            # # l += '{0:->10}{1:->10}{2:->10}{3:->10}\n'.format('pid5',
            # #                                                  'pid6',
            # #                                                  'pid7',
            # #                                                  'pid8')
            # #     """

            # #     """
            # #     '{0:>10}'.format(p)
            # #     """


            # ### ++++++++++++++++++++++++++++++++++++
            # ### Write the bottom nodes as part of a group
            # ### ++++++++++++++++++++++++++++++++++++

            

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
            for ii in range(len(set_pids)):
                f.write('*PART\n')
                f.write(f'Part_{ii+init_sid}\n')
                f.write('')

                f.write('${0:->9}{1:->10}{2:->10}\n'.format('pid', 'secid', 'mid'))
                f.write('{0:>10}{1:>10}{2:>10}\n'.format(ii+init_sid, ii+init_sid, 999))

                f.write('*SECTION_SHELL\n')
                f.write('${0:->9}{1:->10}{2:->10}{3:->10}\n'.format('secid',
                                                                'elform',
                                                                'shrf',
                                                                'nip'))
                f.write('{0:>10}{1:>10}{2:>10.5f}{3:>10}\n'.format(ii+init_sid, 
                                                                params_dict['elform'], #Element formulation (default 2 -> Tsai-Belytschko)
                                                                params_dict['shrf'], #shrf #Shearing factor as in Reissner-Mindlin plate
                                                                params_dict['nip']))  #nip (Integration Points)
                
                f.write('${0:->14}{1:->15}{2:->15}{3:->15}\n'.format('t1', 't2', 't3', 't4'))
                f.write('{0:>15.10f}{0:>15.10f}{0:>15.10f}{0:>15.10f}\n'.format(THICKNESS))


                
            # Elements grouped by partition
            f.write("*ELEMENT_SHELL\n")
            f.write('${0:->7}{1:->8}'.format('eid', 'pid'))
            f.write('{0:->8}{1:->8}{2:->8}{3:->8}\n'.format('n1', 'n2', 'n3', 'n4'))

            for ii in range(len(set_pids)):
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
        gmsh.write("mesh3_crash_tube.vtk")



        self._end_gmsh()