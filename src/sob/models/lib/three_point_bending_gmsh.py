from src.sob.models.lib.gmsh_base_meshes import Template_GMSH_Mesh_Constructor
import gmsh
import json
from pathlib import Path
from typing import Union, List
import os,sys


class ThreePointBending_GMSH(Template_GMSH_Mesh_Constructor):
    def __init__(self, model_name = "Default"):
        super().__init__(model_name)

    @property
    def required_parameters(self)->tuple:
        return (
                "elform",
                "nip", 
                "elsize",
                "h_level",
                "thickness_maps",
                 "grid",
                 'nelx',
                 'nely_div',
                 'nelz_div'
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
        

        ### ==========================================
        ### NON GMSH FUNCTIONS
        ### ==========================================

        def segment_elements_by_z(barycenters:list, elements_ids:list, delta=1e-6):
            """
            Segments elements into separate lists by z-position within a delta threshold.

            Returns:
            - segmented_elements: list of arrays of element indices, grouped by z-layer
            """
            z_vals = np.asarray(barycenters)[:, 2]
            sorted_indices = np.argsort(z_vals)
            sorted_z = z_vals[sorted_indices]

            segmented_elements = []
            current_group = [elements_ids[sorted_indices[0]]]
            current_z_ref = sorted_z[0]

            for idx, z in zip(sorted_indices[1:], sorted_z[1:]):
                if abs(z - current_z_ref) <= delta:
                    current_group.append(elements_ids[idx])
                else:
                    segmented_elements.append(np.array(current_group))
                    current_group = [elements_ids[idx]]
                    current_z_ref = z
            segmented_elements.append(np.array(current_group))  # Final group

            return segmented_elements

        ### ==========================================
        ### PREAMBLE
        ### ==========================================

        # Mesh refinement level
        lv:int = params_dict["h_level"]


        THICKNESS_MAPS = params_dict["thickness_maps"]


        # POINT SIZE
        lc:float = params_dict["elsize"]

        left_side_grid = []

        for aa in np.arange(-40.0,41.0,20.0,dtype=int):
            for bb in np.arange(-60.0,61.0,20.0, dtype=int):
                left_side_grid.append((-400,int(bb),int(aa)))

        ### left points ID
        leftPointsID:list = [*range(1,len(left_side_grid)+1)]
        leftCurvesID:list = []

        #### WRITE THE POINTS ON GMSH

        try:
            gmsh.initialize(sys.argv)
        except:
            print("GMSH couldn't be initialised")

        # Try to not print messages to terminal from GMSH
        # gmsh.option.setNumber("General.Terminal", 0)
        # gmsh.option.setNumber("General.Verbosity", 0)



        gmsh.model.add("trial_three_point")
        gmsh.option.setNumber("Mesh.MeshSizeMin", lc)
        gmsh.option.setNumber("Mesh.MeshSizeMax", lc)

        ### ========================================
        ### ADD POINTS AND CURVES
        ### ========================================

        # Generate the points from the left and add to the mesh
        for ii,point in enumerate(left_side_grid):
            gmsh.model.geo.addPoint(point[0],point[1],point[2],lc,ii+1)

        # Generate the lines on the left
        # Loop to get the horizontal lines first
        for ii in range(1,len(leftPointsID)+1,7):
            curList = [*range(ii,ii+7)]
            for jj in range(len(curList)-1):
                # Generate the point
                leftCurvesID.append(gmsh.model.geo.addLine(curList[jj],curList[jj+1]))

        # Loop to get the vertical lines
        for ii in range(1,len(leftPointsID)-6,7):
            curList = [*range(ii,ii+7)]
            for jj in curList:
                # Generate the point
                leftCurvesID.append(gmsh.model.geo.addLine(jj,jj+7))

        # Add the nodes and lines to corresponding physical groups
        gmsh.model.geo.addPhysicalGroup(0,leftPointsID,1,"leftNodes")
        gmsh.model.geo.addPhysicalGroup(1,leftCurvesID,2,"leftCurves")
        gmsh.model.geo.synchronize()

        # Loop all over the left curves and set them as transfinite
        for iCurve in leftCurvesID:
            gmsh.model.geo.mesh.setTransfiniteCurve(iCurve,params_dict['nely_div']+1,"Progression",1)

        # Synchronize
        gmsh.model.geo.synchronize()

        del ii, jj, curList, aa, bb

        ### =============================
        ### EXTRUSION
        ### =============================

        # Generate some dimTags for the extrusion
        dimTagsExtrusion = []

        for iNode in leftPointsID:
            dimTagsExtrusion.append((0,iNode))

        for iCurve in leftCurvesID:
            dimTagsExtrusion.append((1,iCurve))

        del iNode, iCurve
        # Extrude the points and lines from the right side
        outDimTags = gmsh.model.geo.extrude(dimTags=dimTagsExtrusion,
                                            dx=params_dict['extrusion_length'],
                                            dy=0,
                                            dz=0,
                                            recombine=True
                                            )

        # Classify the new DimTags
        newPoints:list = []
        newLines:list =  []
        newSurfs:list =  []

        # Loop all over the new outDimTags
        for iDim,iTag in outDimTags:
            if iTag >=0:
                if iDim ==0:
                    if iTag not in newPoints:
                        newPoints.append(iTag)
                elif iDim ==1:
                    if iTag not in newLines:
                        newLines.append(iTag)
                elif iDim==2:
                    if iTag not in newSurfs:
                        newSurfs.append(iTag)
        # Synchronize
        gmsh.model.geo.synchronize()

        # Classify the elements corresponding the numbers
        rightPointsID:list = newPoints.copy()
        gmsh.model.geo.addPhysicalGroup(0,rightPointsID,3,"rightNodes")

        # Synchronize
        gmsh.model.geo.synchronize()

        # Loop all over the lines and select the longitudinal and the ones from the right side
        longitudinalLinesID:list = []
        rightCurvesID:list = []

        for iLine in newLines:
            if iLine <= 93: 
                longitudinalLinesID.append(iLine)
            else: 
                rightCurvesID.append(iLine)


        # Generate the groups
        gmsh.model.geo.addPhysicalGroup(1,rightCurvesID,4,"rightCurves")
        gmsh.model.geo.addPhysicalGroup(1,longitudinalLinesID,5,"longitudinalCurves")

        # Loop all over the right curves and set them as transfinite
        for iCurve in rightCurvesID:
            gmsh.model.geo.mesh.setTransfiniteCurve(iCurve,params_dict['nely_div']+1,"Progression",1)

        # Loop all over the longitudinal curves and set them as transfinite
        for iCurve in longitudinalLinesID:
            gmsh.model.geo.mesh.setTransfiniteCurve(iCurve,params_dict['nelx']+1,"Progression",1)

        # delete variables
        del newPoints, newLines, iLine, iCurve

        # Synchronize
        gmsh.model.geo.synchronize()

        # Classify the surfaces depending on the positioning

        outer_core_surfaces:list = [97,101,105,109,113,117, # bottom
                                    193,197,201,205,209,213, # top
                                    217,245,273,301, # Right side
                                    241,269,297,325] # left Side

        surface_3:list = [221,249,277,305] # Right surface
        surface_3_prime:list = [237,265,293,321] # Left surface
        surface_2:list = [225,253,281,309] # Right-middle surface
        surface_2_prime:list = [233,261,289,317] # Left-middle surface
        surface_1:list = [229,257,285,313] # Middle surface


        level_2_surfaces:list = [121,125,129,133,137,141]
        level_3_surfaces:list = [145,149,153,157,161,165]
        level_4_surfaces:list = [169,173,177,181,185,189]

        gmsh.model.geo.synchronize()

        # Loop all over the surfaces
        for iSurf in newSurfs:

            # Set the surfaces as transfinite
            gmsh.model.geo.mesh.setTransfiniteSurface(iSurf,
                                                    arrangement="Left")
            
            # Recombine the elements of the surface
            gmsh.model.geo.mesh.setRecombine(2,iSurf)
            # if iSurf >=121 and iSurf <=141:
            #     level_2_surfaces.append(iSurf)
            # elif iSurf >=145 and iSurf <=165:
            #     level_3_surfaces.append(iSurf)
            # elif iSurf >=169 and iSurf <=189:
            #     level_4_surfaces.append(iSurf)
            
            # Set another algorithm for meshing the surf
            gmsh.model.geo.mesh.setAlgorithm(2,iSurf,3)

        # Delete some variables
        del newSurfs, iSurf

        # Synchronize
        gmsh.model.geo.synchronize()

        # Generate new groups based on the surfaces
        gmsh.model.geo.addPhysicalGroup(2,outer_core_surfaces,6,"outerCoreSurfaces")
        gmsh.model.geo.addPhysicalGroup(2,level_2_surfaces,7,"levelTwoSurfaces")
        gmsh.model.geo.addPhysicalGroup(2,level_3_surfaces,8,"levelThreeSurfaces")
        gmsh.model.geo.addPhysicalGroup(2,level_4_surfaces,9,"levelFourSurfaces")
        gmsh.model.geo.addPhysicalGroup(2,surface_1,10,"surface1")
        gmsh.model.geo.addPhysicalGroup(2,surface_2,11,"surface2")
        gmsh.model.geo.addPhysicalGroup(2,surface_2_prime,12,"surface2Prime")
        gmsh.model.geo.addPhysicalGroup(2,surface_3,13,"surface3")
        gmsh.model.geo.addPhysicalGroup(2,surface_3_prime,14,"surface3Prime")


        # Synchronize
        gmsh.model.geo.synchronize()

        # Remove duplicate elements
        gmsh.model.geo.removeAllDuplicates()
        gmsh.model.geo.synchronize()

        #gmsh.option.setNumber("Mesh.Smoothing", 100)


        ### =========================================
        ### Generate the mesh
        ### =========================================
        gmsh.model.mesh.generate(dim = 2)
        #gmsh.model.mesh.setSmoothing()
        gmsh.model.mesh.removeDuplicateNodes()
        gmsh.model.mesh.removeDuplicateElements()


        ### -------------------------------------- 
        ### Refine depending on the number of dimensions
        ### -----------------------------------gmsh.plugin.setNumber("SimplePartition", "NumSlicesX", 1)
        # gmsh.plugin.setNumber("SimplePartition", "NumSlicesY", 1)
        # gmsh.plugin.setNumber("SimplePartition", "NumSlicesZ", 16*2**(lv-1))
        # gmsh.plugin.run("SimplePartition")---
        for _ in range(lv-1):
            gmsh.model.mesh.refine()


        # gmsh.plugin.setNumber("SimplePartition", "NumSlicesX", 1)
        # gmsh.plugin.setNumber("SimplePartition", "NumSlicesY", 1)
        # gmsh.plugin.setNumber("SimplePartition", "NumSlicesZ", 16*2**(lv-1))
        # gmsh.plugin.run("SimplePartition")

        # # Get the bounding box and then split
        #xx1 =  gmsh.model.getNumberOfPartitions()
        #xx1 =  30*(lv)

        ### =========================================
        ### Get the elements inside a list classified by surface/ physical group ID
        ### =========================================
        element_list:list = []

        barycenters_surf_1 = []
        barycenters_surf_2 = []
        barycenters_surf_3 = []
        barycenters_surf_2_prime = []
        barycenters_surf_3_prime = []


        gmsh.model.geo.synchronize()


        for ii in range(6,15):
            # Get the elements belonging to each partition based on surface ID
            iEntities = gmsh.model.getEntitiesForPhysicalGroup(2,ii)
            elem_list = []
            # Get the elements of the surface
            for iEntity in iEntities:
                elems = gmsh.model.mesh.getElements(dim=2, tag=iEntity)
                elem_list.extend(elems[1])
                # print(elems[1])
                if ii >= 10:
                    # Get the barycenter of the elements
                    partial_bary = gmsh.model.mesh.getBarycenters(elementType=3,
                                                                tag=iEntity,
                                                                fast=True,
                                                                primary=False).reshape(-1,3)
                    if ii == 10:
                        # Append the barycenter to the list
                        barycenters_surf_1.extend(partial_bary.copy())
                    elif ii == 11:
                        # Append the barycenter to the list
                        barycenters_surf_2.extend(partial_bary.copy())
                    elif ii == 12:
                        # Append the barycenter to the list
                        barycenters_surf_2_prime.extend(partial_bary.copy())
                    elif ii == 13:
                        # Append the barycenter to the list
                        barycenters_surf_3.extend(partial_bary.copy())
                    elif ii == 14:
                        # Append the barycenter to the list
                        barycenters_surf_3_prime.extend(partial_bary.copy())
            
            element_list.append(np.asarray(elem_list).ravel().tolist())


        surf_1_segmented_elem_list = segment_elements_by_z(barycenters_surf_1,element_list[4])
        surf_2_segmented_elem_list = segment_elements_by_z(barycenters_surf_2,element_list[5])
        surf_2_prime_segmented_elem_list = segment_elements_by_z(barycenters_surf_2_prime,element_list[6])
        surf_3_segmented_elem_list = segment_elements_by_z(barycenters_surf_3,element_list[7])
        surf_3_prime_segmented_elem_list = segment_elements_by_z(barycenters_surf_3_prime,element_list[8])


        # for ii in range(6,15):
        #     # Get the elements belonging to each partition based on surface ID
        #     iEntities = gmsh.model.getEntitiesForPhysicalGroup(2,ii)

        #     # Get the elements of the surface
        #     for iEntity in iEntities:
        #         elems = gmsh.model.mesh.getElements(dim=2, tag=iEntity)


        #     # Get the name of the entity
        #     name = gmsh.model.getPhysicalName(2,ii)
            


        # Loop all over the partitioned entities
        # entities = gmsh.model.getEntities(dim=2)
        # for e in entities:
        #     # Extract the partition ID
        #     partitions = gmsh.model.getPartitions(e[0], e[1])

        #     partition_element_list= gmsh.model.mesh.getElements(e[0], e[1])[1]

        #     if len(partition_element_list):
        #         for new_elem in partition_element_list:
        #             element_list[partitions[0]-1].extend(new_elem)

        #     if e[0]==2:
        #         if len(partitions):
        #             print("Entity " + str(e) + " of type " +
        #                     gmsh.model.getType(e[0], e[1]))
        #             print(" - Partition(s): " + str(partitions))
        #             print(" - Parent: " + str(gmsh.model.getParent(e[0], e[1])))
        #             print(" - Boundary: " + str(gmsh.model.getBoundary([e])))

        # UNPARTITION
        #gmsh.model.mesh.unpartition()
        gmsh.model.geo.synchronize()

        # Replace the IDs of these nodes
        all_nodes = gmsh.model.mesh.getNodes(dim=-1,returnParametricCoord=False)[0]
        new_set_IDs = np.arange(1000001, 1000001 + len(all_nodes), dtype=int)
        gmsh.model.mesh.renumberNodes(all_nodes, new_set_IDs)



        gmsh.model.geo.synchronize()
        leftNodesMeshID = gmsh.model.mesh.getNodesForPhysicalGroup(1,2)[0]
        rightNodesMeshID = gmsh.model.mesh.getNodesForPhysicalGroup(1,4)[0]


        with open("mesh.txt", "w") as f:
            
            ### NOTE: Write the header of keyword
            f.write("#--------------------------------------------------------------------------------------------------|\n")
            f.write("#---1----|----2----|----3----|----4----|----5----|----6----|----7----|----8----|----9----|---10----|\n")
            #f.write("/BEGIN\n")
            #f.write("TUBE_MESH\n")
            #f.write('{0:>10}{1:>10}\n'.format(2023,0))
            #f.write('{0:>20}{1:>20}{2:>20}\n'.format("Mg","mm","s"))

            ### ++++++++++++++++++++++++++++++++++++++
            ### Write the Nodes List
            ### ++++++++++++++++++++++++++++++++++++++
            f.write("#---1----|----2----|----3----|----4----|----5----|----6----|----7----|----8----|----9----|---10----|\n")
            f.write("#--------------------------------------------------------------------------------------------------|\n")
            f.write("#---1----|----2----|----3----|----4----|----5----|----6----|----7----|----8----|----9----|---10----|\n")
            f.write("/NODE\n")

            node_tags, node_coords, _ = gmsh.model.mesh.getNodes()

            f.write('#{0:->10}{1:->20}{2:->20}{3:->20}\n'.format('nid', 'x', 'y', 'z'))
            for i in range(len(node_tags)):
                f.write('{0:>10}{1:20.5f}{2:20.5f}{3:20.5f}\n'.format(
                                                                    int(node_tags[i]),  #nid
                                                                        node_coords[3*i], #x
                                                                        node_coords[3*i+1], #y 
                                                                        node_coords[3*i+2] #z
                                                                    ))
            

            ### ++++++++++++++++++++++++++++++++++
            ### Write the Nodes from the left side
            ### ++++++++++++++++++++++++++++++++++

            f.write("#---1----|----2----|----3----|----4----|----5----|----6----|----7----|----8----|----9----|---10----|\n")
            f.write("#--------------------------------------------------------------------------------------------------|\n")
            f.write("#---1----|----2----|----3----|----4----|----5----|----6----|----7----|----8----|----9----|---10----|\n")


            f.write("/GRNOD/NODE/101\n")
            f.write("LEFT NODES\n")

            l =""

            l += '#{0:->10}{1:->10}{2:->10}{3:->10}{4:->10}'.format('nid1',
                                                            'nid2',
                                                            'nid3',
                                                            'nid4',
                                                            "nid5")

            l += '{0:->10}{1:->10}{2:->10}{3:->10}{4:->10}\n'.format('nid6',
                                                                'nid7',
                                                                'nid8',
                                                                'nid9',
                                                                'nid10')

            leftNodesMeshID.sort()
            for i in range(0, len(leftNodesMeshID), 10):
                chunk = leftNodesMeshID[i:i+10]
                for n in chunk:
                    l += '{0:>10}'.format(n)
                l += '\n'
            
            f.write(l)


            f.write("/GRNOD/NODE/102\n")
            f.write("RIGHT NODES\n")

            l =""

            l += '#{0:->10}{1:->10}{2:->10}{3:->10}{4:->10}'.format('nid1',
                                                            'nid2',
                                                            'nid3',
                                                            'nid4',
                                                            "nid5")

            l += '{0:->10}{1:->10}{2:->10}{3:->10}{4:->10}\n'.format('nid6',
                                                                'nid7',
                                                                'nid8',
                                                                'nid9',
                                                                'nid10')

            rightNodesMeshID.sort()
            for i in range(0, len(rightNodesMeshID), 10):
                chunk = rightNodesMeshID[i:i+10]
                for n in chunk:
                    l += '{0:>10}'.format(n)
                l += '\n'
            
            f.write(l)



            ### ++++++++++++++++++++++++++++++++++
            ### Write the Parts
            ### ++++++++++++++++++++++++++++++++++


            # Write the outer core surfaces
            f.write("#---1----|----2----|----3----|----4----|----5----|----6----|----7----|----8----|----9----|---10----|\n")
            f.write("#--------------------------------------------------------------------------------------------------|\n")
            f.write("#---1----|----2----|----3----|----4----|----5----|----6----|----7----|----8----|----9----|---10----|\n")


            f.write("/PART/101\n")
            f.write("OUTER_CORE_SURFACES\n")
            f.write('#{0:->10}{1:->10}{2:>10}\n'.format('prop_ID', 'mat_ID','subset_id'))
            f.write('{0:>10}{1:>10}{2:>10}\n'.format(1, 1, 0))

            f.write("/SHELL/101\n")
            # Write the element lists
            f.write('#{0:->10}'.format('eid'))
            f.write('{0:->10}{1:->10}{2:->10}{3:->10}\n'.format('n1', 'n2', 'n3', 'n4'))
            for ii in range(len(element_list[0])):

                n_eid = gmsh.model.mesh.getElement(element_list[0][ii])[1]

                f.write('{0:>10}{1:>10}{2:>10}{3:>10}{4:>10}\n'.format(element_list[0][ii],
                                                                    n_eid[0], 
                                                                    n_eid[1], 
                                                                    n_eid[2], 
                                                                    n_eid[3]))
            # Write the level surfaces
            f.write("#---1----|----2----|----3----|----4----|----5----|----6----|----7----|----8----|----9----|---10----|\n")
            f.write("#--------------------------------------------------------------------------------------------------|\n")
            f.write("#---1----|----2----|----3----|----4----|----5----|----6----|----7----|----8----|----9----|---10----|\n")


            f.write("/PART/102\n")
            f.write("LEVEL_SURFACES\n")
            f.write('#{0:->10}{1:->10}{2:>10}\n'.format('prop_ID', 'mat_ID','subset_id'))
            f.write('{0:>10}{1:>10}{2:>10}\n'.format(2, 1, 0))

            f.write("/SHELL/102\n")
            # Write the element lists
            f.write('#{0:->10}'.format('eid'))
            f.write('{0:->10}{1:->10}{2:->10}{3:->10}\n'.format('n1', 'n2', 'n3', 'n4'))
            
            elem_list_level_surfs = []
            for ii in range(1, 4):
                elem_list_level_surfs.extend(element_list[ii])

            for ii in range(len(elem_list_level_surfs)):

                n_eid = gmsh.model.mesh.getElement(elem_list_level_surfs[ii])[1]

                f.write('{0:>10}{1:>10}{2:>10}{3:>10}{4:>10}\n'.format(elem_list_level_surfs[ii],
                                                                    n_eid[0], 
                                                                    n_eid[1], 
                                                                    n_eid[2], 
                                                                    n_eid[3]))
            
            # Write the surface 1 (with sub-references)
            f.write("#---1----|----2----|----3----|----4----|----5----|----6----|----7----|----8----|----9----|---10----|\n")
            f.write("#--------------------------------------------------------------------------------------------------|\n")
            f.write("#---1----|----2----|----3----|----4----|----5----|----6----|----7----|----8----|----9----|---10----|\n")



            # WRITE THE ELEMENTS OF THE SURFACE 1
            f.write("/PART/1000\n")
            f.write("SURFACE_1\n")
            f.write('#{0:->10}{1:->10}{2:>10}\n'.format('prop_ID', 'mat_ID','subset_id'))
            f.write('{0:>10}{1:>10}{2:>10}\n'.format(3, 1, 0))
            f.write("/SHELL/1000\n")

            # Write the element lists
            f.write('#{0:->10}'.format('eid'))
            f.write('{0:->10}{1:->10}{2:->10}{3:->10}{4:->30}{5:->20}\n'.format('n1', 'n2', 'n3', 'n4','phi','thick'))

            for ii in range(len(surf_1_segmented_elem_list)):
                # Get the elements of the surface
                for jj in range(len(surf_1_segmented_elem_list[ii])):
                    n_eid = np.asarray(gmsh.model.mesh.getElement(surf_1_segmented_elem_list[ii][jj])[1])

                    f.write('{0:>10}{1:>10}{2:>10}{3:>10}{4:>10}{5:>30}{6:>20.5f}\n'.format(surf_1_segmented_elem_list[ii][jj],
                                                                        n_eid[0], 
                                                                        n_eid[1], 
                                                                        n_eid[2], 
                                                                        n_eid[3],
                                                                        0.0,
                                                                        THICKNESS_MAPS[0][str(ii+1)]))

            # WRITE THE ELEMENTS OF THE SURFACE 2
            f.write("/PART/2000\n")
            f.write("SURFACE_2\n")
            f.write('#{0:->10}{1:->10}{2:>10}\n'.format('prop_ID', 'mat_ID','subset_id'))
            f.write('{0:>10}{1:>10}{2:>10}\n'.format(3, 1, 0))
            f.write("/SHELL/2000\n")

            # Write the element lists
            f.write('#{0:->10}'.format('eid'))
            f.write('{0:->10}{1:->10}{2:->10}{3:->10}{4:->30}{5:->20}\n'.format('n1', 'n2', 'n3', 'n4','phi','thick'))

            for ii in range(len(surf_2_segmented_elem_list)):
                # Get the elements of the surface
                for jj in range(len(surf_2_segmented_elem_list[ii])):
                    n_eid = np.asarray(gmsh.model.mesh.getElement(surf_2_segmented_elem_list[ii][jj])[1])

                    f.write('{0:>10}{1:>10}{2:>10}{3:>10}{4:>10}{5:>30}{6:>20.5f}\n'.format(surf_2_segmented_elem_list[ii][jj],
                                                                        n_eid[0], 
                                                                        n_eid[1], 
                                                                        n_eid[2], 
                                                                        n_eid[3],
                                                                        0.0,
                                                                        THICKNESS_MAPS[1][str(ii+1)]))



            # WRITE THE ELEMENTS OF THE SURFACE 2 PRIME
            f.write("/PART/2001\n")
            f.write("SURFACE_2_PRIME\n")
            f.write('#{0:->10}{1:->10}{2:>10}\n'.format('prop_ID', 'mat_ID','subset_id'))
            f.write('{0:>10}{1:>10}{2:>10}\n'.format(3, 1, 0))
            f.write("/SHELL/2001\n")

            # Write the element lists
            f.write('#{0:->10}'.format('eid'))
            f.write('{0:->10}{1:->10}{2:->10}{3:->10}{4:->30}{5:->20}\n'.format('n1', 'n2', 'n3', 'n4','phi','thick'))

            for ii in range(len(surf_2_prime_segmented_elem_list)):
                # Get the elements of the surface
                for jj in range(len(surf_2_prime_segmented_elem_list[ii])):
                    n_eid = np.asarray(gmsh.model.mesh.getElement(surf_2_prime_segmented_elem_list[ii][jj])[1])

                    f.write('{0:>10}{1:>10}{2:>10}{3:>10}{4:>10}{5:>30}{6:>20.5f}\n'.format(surf_2_prime_segmented_elem_list[ii][jj],
                                                                        n_eid[0], 
                                                                        n_eid[1], 
                                                                        n_eid[2], 
                                                                        n_eid[3],
                                                                        0.0,
                                                                        THICKNESS_MAPS[2][str(ii+1)]))

            # WRITE THE ELEMENTS OF THE SURFACE 3
            f.write("/PART/3000\n")
            f.write("SURFACE_3\n")
            f.write('#{0:->10}{1:->10}{2:>10}\n'.format('prop_ID', 'mat_ID','subset_id'))
            f.write('{0:>10}{1:>10}{2:>10}\n'.format(3, 1, 0))
            f.write("/SHELL/3000\n")

            # Write the element lists
            f.write('#{0:->10}'.format('eid'))
            f.write('{0:->10}{1:->10}{2:->10}{3:->10}{4:->30}{5:->20}\n'.format('n1', 'n2', 'n3', 'n4','phi','thick'))

            for ii in range(len(surf_3_segmented_elem_list)):
                # Get the elements of the surface
                for jj in range(len(surf_3_segmented_elem_list[ii])):
                    n_eid = np.asarray(gmsh.model.mesh.getElement(surf_3_segmented_elem_list[ii][jj])[1])

                    f.write('{0:>10}{1:>10}{2:>10}{3:>10}{4:>10}{5:>30}{6:>20.5f}\n'.format(surf_3_segmented_elem_list[ii][jj],
                                                                        n_eid[0], 
                                                                        n_eid[1], 
                                                                        n_eid[2], 
                                                                        n_eid[3],
                                                                        0.0,
                                                                        THICKNESS_MAPS[3][str(ii+1)]))
            
            # WRITE THE ELEMENTS OF THE SURFACE 3 PRIME
            f.write("/PART/3001\n")
            f.write("SURFACE_3_PRIME\n")
            f.write('#{0:->10}{1:->10}{2:>10}\n'.format('prop_ID', 'mat_ID','subset_id'))
            f.write('{0:>10}{1:>10}{2:>10}\n'.format(3, 1, 0))
            f.write("/SHELL/3001\n")

            # Write the element lists
            f.write('#{0:->10}'.format('eid'))
            f.write('{0:->10}{1:->10}{2:->10}{3:->10}{4:->30}{5:->20}\n'.format('n1', 'n2', 'n3', 'n4','phi','thick'))

            for ii in range(len(surf_3_prime_segmented_elem_list)):
                # Get the elements of the surface
                for jj in range(len(surf_3_prime_segmented_elem_list[ii])):
                    n_eid = np.asarray(gmsh.model.mesh.getElement(surf_3_prime_segmented_elem_list[ii][jj])[1])

                    f.write('{0:>10}{1:>10}{2:>10}{3:>10}{4:>10}{5:>30}{6:>20.5f}\n'.format(surf_3_prime_segmented_elem_list[ii][jj],
                                                                        n_eid[0], 
                                                                        n_eid[1], 
                                                                        n_eid[2], 
                                                                        n_eid[3],
                                                                        0.0,
                                                                        THICKNESS_MAPS[4][str(ii+1)]))


            # Write a surface group in Radioss
            f.write("#---1----|----2----|----3----|----4----|----5----|----6----|----7----|----8----|----9----|---10----|\n")
            f.write("#--------------------------------------------------------------------------------------------------|\n")
            f.write("#---1----|----2----|----3----|----4----|----5----|----6----|----7----|----8----|----9----|---10----|\n")

            f.write("/SURF/PART/EXT/4\n")
            f.write("SURFACES_INTERFACE\n")
            f.write('{0:>10}{1:>10}{2:>10}{3:>10}{4:>10}{5:>10}{6:>10}\n'.format(101, 102,1000,2000,2001,3000,3001))


            f.write("#enddata\n")
            f.write("/END\n")


        # UNPARTITION
        #gmsh.model.mesh.unpartition()

        # ... and save it to disk
        #gmsh.write("t1.msh")
        gmsh.write("t2.geo_unrolled")
        gmsh.write("mesh2_three_point.vtk")
        gmsh.write("mesh2_three_point.rad")
        #gmsh.write("mesh2.key")



        self._end_gmsh()