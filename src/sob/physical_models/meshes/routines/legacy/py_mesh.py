# -*- coding: utf-8 -*-
"""

Created on October 2016

@author: Pablo Lozano, pablo.lozano@tum.de

Modified by K.Komeilizadeh
"""
from __future__ import print_function #  K.Komeilizadeh
import sys
import numpy as np

def py_mesh(var_file, save=True):

    # Reading variables
    # reads var_file line by line and processes its contents. 
    variables = {} # variables dictionary
    grids = {}
    cells = {}

    make_triggers = False
    multiple_triggers = False
    with open(var_file) as fi:
        for line in fi:

            if len(line) == 1: 
                continue # it is skipped, as it is likely an empty line
            elif line.startswith('#'):
                continue # indicating a comment, also be skipped
            
            # Split the line into a list called tmp using a comma as the delimiter
            tmp = line.strip().split(',') 
            
            # first string is the keyword
            if 'elsize' == tmp[0]:
                variables['elsize'] = float(tmp[1])

            if 'extrusion_length' == tmp[0]:
                variables['extrusion_length'] = float(tmp[1])

            if 'id_min' == tmp[0]:
                variables['id_min'] = int(tmp[1])

            if 'trigger_depth' == tmp[0]:
                variables['trigger_depth'] = float(tmp[1])
                make_triggers = True # ??
            

            if 'trigger_rows' == tmp[0]:
                variables['trigger_rows'] = int(tmp[1])

            # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! newer py_mesh, K.Komeilizadeh                
            if 'trigger' == tmp[0]:
                if 'trigger' not in variables.keys(): # if trigger doesn't exist, create a list
                    variables['trigger'] = []
                variables['trigger'].append(int(tmp[1])) # otherwise add tmp[1] into list                

            if 'mid' == tmp[0]:
                variables['mid'] = int(tmp[1])

            if 'elform' == tmp[0]:
                variables['elform'] = int(tmp[1])

            if 'shrf' == tmp[0]:
                variables['shrf'] = float(tmp[1])

            if 'nip' == tmp[0]:
                variables['nip'] = int(tmp[1])

            if 'database' == tmp[0]:
                if 'database' not in variables.keys():
                    variables['database'] = []
                variables['database'].append(int(tmp[1]))

            # Geometric object is constructed from its bounds: vertices, edges, and face
            # Vertices:
            if 'grid' == tmp[0]:
                gid = int(tmp[1]) # Convert tmp[1] to an integer
                gx = float(tmp[2]) # Convert to an float, x coordinate
                gy = float(tmp[3]) # Convert to an float, y coordinate
                if gid not in grids:
                    grids[gid] = [gx, gy]
                else:
                    print ('ERROR: duplicated grid:', gid)
                    sys.exit(1)
            # Edges:
            if 'cell' == tmp[0]:
                cid = int(tmp[1])
                g0 = int(tmp[2])
                g1 = int(tmp[3])
                t = float(tmp[4])
                if cid not in cells:
                    cells[cid] = {'line': [g0, g1], 'thickness': t}
                # cells dictionary will look like this: 1,2,3 are cell IDs
                # cells = {1: {'line': [0, 1], 'thickness': 0.1},
                #          2: {'line': [2, 3], 'thickness': 0.2},
                #          3: {'line': [4, 5], 'thickness': 0.3}}

                else:
                    print ('ERROR: duplicated cell:', cid)
                    sys.exit(1)
    #print
    # ---- read file finished, all data is stored in "variables", "grids", and "cells"

    # Meshing parameters
    elsize = variables['elsize']
    extrusion_length = variables['extrusion_length']
    id_min = variables['id_min']

    if make_triggers:
        trigger_depth = variables['trigger_depth']
        if 'trigger_rows' in variables:
            trigger_rows = variables['trigger_rows']
        else:
            trigger_rows = 2

    # Material parameters
    mid = variables['mid']

    # Section parameters
    if 'elform' in variables:
        elform = variables['elform']
    else:
        elform = 2

    if 'shrf' in variables:
        shrf = variables['shrf']
    else:
        shrf = 1.0

    if 'nip' in variables:
        nip = variables['nip']
    else:
        nip = 5

    # Mesh generation

    nodes = {}
    shells = {}

    common_nodes = {}

    set_nodes = {}
    set_nodes[101] = []
    set_nodes[102] = []
    
    # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! added K.Komeilizadeh
    nodes_of_parts_at_top = {}  # set of nodes at top that belongs to each part
    nodes_of_parts_at_bottom = {}
    nodes_of_parts_at_top_coordinate = {}
    nodes_of_parts_at_bottom_coordinate = {}
    for pid in sorted(cells.keys()):
        nodes_of_parts_at_top[pid]=[]  
        nodes_of_parts_at_bottom[pid]=[]   
    # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 
    # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! from newer py_mesh, K.Komeilizadeh
    database_grids = variables['database']
    database_nodes = []

    if 'trigger' in variables.keys():
        trigger_cells = variables['trigger']
    else:
        trigger_cells = []
		
    nid = id_min
    eid = id_min

    # "grids"(1001 to 1008) store the ids and coordinates of nodes which define the geometry
    # "cells"(101 to 108): each cell related to two nodes(from 1001 to 1008), and thickness

    # Loop over the cell ids (Edges)
    for pid in sorted(cells.keys()):
        # get the id of nodes which are related to current cell
        grid0 = cells[pid]['line'][0] # g0 (from 1001 to 1008)
        grid1 = cells[pid]['line'][1] # g1
        # get the coordinates of these two nodes
        p0 = np.array(grids[grid0]) # corresponding [x, y] in grids
        p1 = np.array(grids[grid1])
        # compute distance of these two nodes
        length = np.linalg.norm(p1 - p0)
        # compute unit vector pointing in the direction of the line segment.
        line_dir = (p1 - p0)/length
        # calculates a vector perpendicular to the line segment
        line_perp = np.array([-line_dir[1], line_dir[0]])

        '''
        Define the number of elements along the line (num_elements_line) 
        Calculate the element size (elsize_line) along the line based on the specified value of 12.
        '''
        num_elements_line = 12
        elsize_line = length/num_elements_line
        '''
        Modified, the element size along the line adapts to the line length
        
        # num_elements_line = int(round(length/elsize))
        # elsize_line = length/num_elements_line
        '''
        '''
        Define the number of elements along the extrusion length
        ''' 
        num_elements_ext = int(round(extrusion_length/elsize))
        # element size in the extrusion direction
        elsize_ext = extrusion_length/num_elements_ext

        # the node ids are stored in a matrix, rows are for the nodes in extrusion length direction, columns are for the nodes in edge direction
        nodes_ids = np.zeros((num_elements_ext+1, num_elements_line+1),
                             dtype=int)

        if grid0 in common_nodes and not grid1 in common_nodes:
            nodes_ids[:,0] = np.array(common_nodes[grid0])
            for j in range(num_elements_ext+1):
                z = j*elsize_ext
                for i in range(1, num_elements_line+1):
                    pi = p0 + i*elsize_line*line_dir
                    grid = np.append(pi, z)
                    nodes[nid] = grid
                    nodes_ids[j, i] = nid
                    nid += 1
            edge1 = nodes_ids[:,-1]
            common_nodes[grid1] = edge1
        elif grid1 in common_nodes and not grid0 in common_nodes:
            nodes_ids[:,-1] = np.array(common_nodes[grid1])
            for j in range(num_elements_ext+1):
                z = j*elsize_ext
                for i in range(num_elements_line):
                    pi = p0 + i*elsize_line*line_dir
                    grid = np.append(pi, z)
                    nodes[nid] = grid
                    nodes_ids[j, i] = nid
                    nid += 1
            edge0 = nodes_ids[:,0]
            common_nodes[grid0] = edge0
        elif grid0 in common_nodes and grid1 in common_nodes:
            nodes_ids[:,0] = np.array(common_nodes[grid0])
            nodes_ids[:,-1] = np.array(common_nodes[grid1])
            for j in range(num_elements_ext+1):
                z = j*elsize_ext
                for i in range(1, num_elements_line):
                    pi = p0 + i*elsize_line*line_dir
                    grid = np.append(pi, z)
                    nodes[nid] = grid
                    nodes_ids[j, i] = nid
                    nid += 1
        else:
            for j in range(num_elements_ext+1):
                z = j*elsize_ext
                for i in range(num_elements_line+1):
                    pi = p0 + i*elsize_line*line_dir
                    grid = np.append(pi, z)
                    nodes[nid] = grid
                    nodes_ids[j, i] = nid
                    nid += 1
            edge0 = nodes_ids[:,0]
            edge1 = nodes_ids[:,-1]
            common_nodes[grid0] = edge0
            common_nodes[grid1] = edge1

        nid_replace0 = nodes_ids[-1,  0]
        nid_replace1 = nodes_ids[-1, -1]
        nodes[grid0] = nodes.pop(nid_replace0)
        nodes[grid1] = nodes.pop(nid_replace1)
        nodes_ids[-1,  0] = grid0
        nodes_ids[-1, -1] = grid1

        # if grid0 in database_grids:
        #     n = nodes_ids[-1, 0]
        #     if n not in database_nodes:
        #         database_nodes.append(n)
        # if grid1 in database_grids:
        #     n = nodes_ids[-1, -1]
        #     if n not in database_nodes:
        #         database_nodes.append(n)

        for n in nodes_ids[0,:]:
            if n not in set_nodes[101]:
                set_nodes[101].append(n)
                # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! added by K.Komeilizadeh 
                nodes_of_parts_at_bottom[pid].append(n)
                nodes_of_parts_at_bottom_coordinate[n] = nodes[n]
                # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   [pid]=[] 
                
        for n in nodes_ids[-1,:]:
            if n not in set_nodes[102]:
                set_nodes[102].append(n)
                # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! added by K.Komeilizadeh
                nodes_of_parts_at_top[pid].append(n) 
                nodes_of_parts_at_top_coordinate[n] = nodes[n]
                # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
#        # !!!!!!!!!!!!!!!!!!!!!!!!!!!!! original trigger        
#        if make_triggers:
#            trigger_nodes = nodes_ids[-3,1:-1]
#            for i in range(1, trigger_rows):
#                trigger_nodes = np.append(trigger_nodes, nodes_ids[-3-i,1:-1])
#            for gid in trigger_nodes:
#                xy = nodes[gid][:2]
#                z = nodes[gid][2]
#                xy_t = xy + trigger_depth*line_perp
#                grid = np.append(xy_t, z)
#                nodes[gid] = grid
		# !!!!!!!!!!!!!!!!!!!!!!!!!!!!! new pymesh trigger, K.Komeilizadeh
        
        if make_triggers:
            if pid in trigger_cells:
                trigger_nodes = nodes_ids[-6,1:-1]  #!!!!!!!!! changed from -3 to -6
                for i in range(1, trigger_rows):
                    trigger_nodes = np.append(trigger_nodes, nodes_ids[-6-i,1:-1]) #!!!!!!!!! changed from -3 to -6
                for gid in trigger_nodes:
                    xy = nodes[gid][:2]
                    z = nodes[gid][2]
                    xy_t = xy + trigger_depth*line_perp
                    grid = np.append(xy_t, z)
                    nodes[gid] = grid
        
                    
        for i in range(len(nodes_ids)-1):
            for j in range(len(nodes_ids[i])-1):
                n1 = nodes_ids[i  , j  ]
                n2 = nodes_ids[i  , j+1]
                n3 = nodes_ids[i+1, j+1]
                n4 = nodes_ids[i+1, j  ]
                shells[eid] = [pid, n1, n2, n3, n4]
                eid += 1

    # Write mesh
    l = ''

    l += '*KEYWORD\n'

    # # Database nodes
    # l += '*DATABASE_HISTORY_NODE\n'
    # l += '${0:->9}{1:->10}{2:->10}{3:->10}'.format('nid1',
    #                                                'nid2',
    #                                                'nid3',
    #                                                'nid4')
    # l += '{0:->10}{1:->10}{2:->10}{3:->10}\n'.format('nid5',
    #                                                  'nid6',
    #                                                  'nid7',
    #                                                  'nid8')
    # database_nodes.sort()
    # for i in range(0, len(database_nodes), 8):
    #     chunk = database_nodes[i:i+8]
    #     for n in chunk:
    #         l += '{0:>10}'.format(n)
    #     l += '\n'

    # Part sets
    l += '*SET_PART_LIST\n'
    l += '${0:->9}\n'.format('sid')
    l += '{0:>10}\n'.format(101)
    l += '${0:->9}{1:->10}{2:->10}{3:->10}'.format('pid1',
                                                   'pid2',
                                                   'pid3',
                                                   'pid4')
    l += '{0:->10}{1:->10}{2:->10}{3:->10}\n'.format('pid5',
                                                     'pid6',
                                                     'pid7',
                                                     'pid8')
    set_pids = sorted(cells.keys())
    for i in range(0, len(set_pids), 8):
        chunk = set_pids[i:i+8]
        for p in chunk:
            l += '{0:>10}'.format(p)
        l += '\n'

    # Node sets
    for sid in set_nodes:
        l += '*SET_NODE_LIST\n'
        l += '${0:->9}\n'.format('sid')
        l += '{0:>10}\n'.format(sid)
        l += '${0:->9}{1:->10}{2:->10}{3:->10}'.format('nid1',
                                                       'nid2',
                                                       'nid3',
                                                       'nid4')
        l += '{0:->10}{1:->10}{2:->10}{3:->10}\n'.format('nid5',
                                                         'nid6',
                                                         'nid7',
                                                         'nid8')
        set_nodes[sid].sort()
        for i in range(0, len(set_nodes[sid]), 8):
            chunk = set_nodes[sid][i:i+8]
            for n in chunk:
                l += '{0:>10}'.format(n)
            l += '\n'

    # Parts
    for pid in sorted(cells.keys()):
        secid = pid
        thickness = cells[pid]['thickness']

        l += '*PART\n'
        l += 'Part_{0:}\n'.format(pid)
        l += '${0:->9}{1:->10}{2:->10}\n'.format('pid', 'secid', 'mid')
        l += '{0:>10}{1:>10}{2:>10}\n'.format(pid, secid, mid)

        l += '*SECTION_SHELL\n'
        l += '${0:->9}{1:->10}{2:->10}{3:->10}\n'.format('secid',
                                                         'elform',
                                                         'shrf',
                                                         'nip')
        l += '{0:>10}{1:>10}{2:>10}{3:>10}\n'.format(secid, elform, shrf, nip)

        l += '${0:->9}{1:->10}{2:->10}{3:->10}\n'.format('t1', 't2', 't3', 't4')
        l += '{0:>10}{0:>10}{0:>10}{0:>10}\n'.format(thickness)

    # Elements
    l += '*ELEMENT_SHELL\n'
    l += '${0:->7}{1:->8}'.format('eid', 'pid')
    l += '{0:->8}{1:->8}{2:->8}{3:->8}\n'.format('n1', 'n2', 'n3', 'n4')

    for s in sorted(shells.keys()):
        eid = s
        pid = shells[s][0]
        n1 = shells[s][1]
        n2 = shells[s][2]
        n3 = shells[s][3]
        n4 = shells[s][4]
        l += '{0:>8}{1:>8}'.format(eid, pid)
        l += '{0:>8}{1:>8}{2:>8}{3:>8}\n'.format(n1, n2, n3, n4)

    # Nodes
    l += '*NODE\n'
    l += '${0:->7}{1:->16}{2:->16}{3:->16}\n'.format('nid', 'x', 'y', 'z')
    for n in sorted(nodes.keys()):
        nid = n
        x = nodes[n][0]
        y = nodes[n][1]
        z = nodes[n][2]
        l += '{0:>8}{1:16.7f}{2:16.7f}{3:16.7f}\n'.format(nid, x, y, z)
		
	# added !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    #       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!	
    
    l +='*END\n'

    if save:
        with open('mesh.k', 'w') as fo:
            fo.write(l)

    #return l
    # !!!!!!!!!!!!!!!!!!!!!!!!!!!!! added by K.Komeilizadeh
    return nodes_of_parts_at_top, nodes_of_parts_at_bottom,\
           nodes_of_parts_at_top_coordinate, nodes_of_parts_at_bottom_coordinate
    # !!!!!!!!!!!!!!!!!!!!!!!

def main():

    try:
        var_file = sys.argv[1]
    except IndexError:
        var_file = 'py_mesh.input'

    mesh = py_mesh(var_file, save=True)
    print ()
    print ('N o r m a l')

    #return 0
    #!!!!!!!!!!!!!!!!!!!!!!! added by K.Komeilizadeh
    return mesh
    #!!!!!!!!!!!!!!!!!!!!!!!

if __name__ == '__main__':
    main()
