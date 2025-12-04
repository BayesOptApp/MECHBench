import numpy as np
import os

from src.sob.physical_models.fem_settings.abstractFEMSettings import AbstractFEMSettings
from src.sob.physical_models.meshes import ThreePointBendingMesh


class ThreePointBendingModel(AbstractFEMSettings):
    def __init__(self, mesh:ThreePointBendingMesh) -> None:

        self.mesh = mesh
        self.mesh.write_mesh_file()
    
    @property
    def material_card_type(self):
        return "OpenRadioss"
    
    @property
    def wall_n_id(self)->int:
        r"""
        Returns the wall node id
        """
        return 99999


    
    @property
    def wall_loc(self)->float:
        r"""
        Returns the wall location
        """
        return 77.5
    
    @property
    def rigid_mass(self)->float:
        r"""
        Returns the wall mass (in metric tons)
        """
        return 0.086
    
    @property
    def wall_vel(self)->float:
        r"""
        Returns the wall velocity
        """
        return -10000.0
        #return -11111.11111
        #return -10277.777777778
    
    @property
    def impactor_diameter(self)->float:
        r"""
        Returns the impactor diameter
        """
        return 70.0
    
    @property
    def material_density(self)->float:
        r"""
        Returns the material density
        """
        return 2.7E-9
    
    @property
    def mat_young_mod(self)->float:
        r"""
        Returns the material young modulus (in MPa)
        """
        #return 210000
        return 70000
    
    @property
    def mat_poisson_r(self)->float:
        r"""
        Returns the material poisson ratio
        """
        return 0.33
    
    @property
    def impactor_offset(self)->float:
        return 2.50
    
    def absorbed_energy(self):
        r'''
        Returns the initial kinetic energy
        '''
        return (self.rigid_mass*1000)*(self.wall_vel/1000)**2/2

    def merge_files(self, output_file, input_files):
        # Function to merge multiple files into one
        with open(output_file, 'w') as output:
            for input_file in input_files:
                with open(input_file, 'r') as input:
                    output.write(input.read())

    def write_shell_property(self, thickness_list, property_ids):
        adr = os.path.join(os.getcwd(),'shell.rad') 
        inf = open(adr, 'w')
        inf.write('\n#---1----|----2----|----3----|----4----|----5----|----6----|----7----|----8----|----9----|---10----|\n')
        inf.write('#-  6. GEOMETRICAL SETS:\n')
        inf.write('#---1----|----2----|----3----|----4----|----5----|----6----|----7----|----8----|----9----|---10----|\n')
        ########################################################################################
        inf.write("/PROP/SHELL/1\n")
        inf.write("PID_tube\n")
        inf.write("#   Ishell    Ismstr     Ish3n    Idrill                            P_thick_fail\n")
        inf.writelines("{:>10}{:>10}{:>10}{:>10}{:>10}{:>10}{:>10}{:>10}{:>10}{:>10}\n".format('24','1','2', '2','', '', '', '0', '', ''))
        inf.write("#                 hm                  hf                  hr                  dm                  dn\n")
        inf.writelines("{:>10}{:>10}{:>10}{:>10}{:>10}{:>10}{:>10}{:>10}{:>10}{:>10}\n".format('','0','', '0','','0','','0','','0'))
        # inf.write("                   0                   0                   0                   0                   0\n")
        inf.write("#        N   Istrain               Thick              Ashear              Ithick     Iplas\n")
        inf.writelines("{:>10}{:>10}{:>10}{:>10}{:>10}{:>10}{:>10}{:>10}{:>10}{:>10}\n".format('5','0','', '1.8','', '0', '', '1', '1',''))
        #inf.writelines("{:>10}{:>10}{:>10}{:>10}{:>10}{:>10}{:>10}{:>10}{:>10}{:>10}\n".format('5','0','', '1.8','', str(5/6), '', '1', '1',''))
        # inf.write("         5         0                 1.8                   0                   1         1\n")
        inf.write("#---1----|----2----|----3----|----4----|----5----|----6----|----7----|----8----|----9----|---10----|\n")
        ########################################################################################
        for i in range(0, len(property_ids)):
            inf.write('/PROP/SHELL/'+str(property_ids[i])+'\n')
            inf.write('PID_v'+str(i+1)+'\n')
            inf.write("#   Ishell    Ismstr     Ish3n    Idrill                            P_thick_fail\n")
            inf.write("        24         1         2         2                                       0\n")
            inf.write("#                 hm                  hf                  hr                  dm                  dn\n")
            inf.writelines("{:>10}{:>10}{:>10}{:>10}{:>10}{:>10}{:>10}{:>10}{:>10}{:>10}\n".format('','0','', '0','','0','','0','','0'))
            inf.write("#        N   Istrain               Thick              Ashear              Ithick     Iplas\n")
            #inf.writelines("{:>10}{:>10}{:>10}{:>10}{:>10}{:>10}{:>10}{:>10}{:>10}{:>10}\n".format('5','0','', str(thickness_list[i]),'', '1', '', '1', '1',''))
            inf.writelines("{:>10}{:>10}{:>10}{:>10}{:>10}{:>10}{:>10}{:>10}{:>10}{:>10}\n".format('5','0','', str(thickness_list[i]),'', str(5/6), '', '1', '1',''))
            # inf.write("{:>40}\n".format(str(thickness_list[i])))
        ######################################################################################## 
        inf.write("#---1----|----2----|----3----|----4----|----5----|----6----|----7----|----8----|----9----|---10----|\n")   
        inf.write("/PROP/SHELL/7\n")
        inf.write("PID_h\n")
        inf.write("#   Ishell    Ismstr     Ish3n    Idrill                            P_thick_fail\n")
        inf.write("        24         1         2         2                                       0\n")
        inf.write("#                 hm                  hf                  hr                  dm                  dn\n")
        inf.write("                   0                   0                   0                   0                   0\n")
        inf.write("#        N   Istrain               Thick              Ashear              Ithick     Iplas\n")
        #inf.write("         5         0                  .7                   1                   1         1\n")
        inf.write(f"         5         0                  .7                   {str(5/6)}                   1         1\n")
        inf.write("#---1----|----2----|----3----|----4----|----5----|----6----|----7----|----8----|----9----|---10----|\n")
        inf.write("/END\n")
        inf.write("#---1----|----2----|----3----|----4----|----5----|----6----|----7----|----8----|----9----|---10----|\n")
    

    def write_input_file(self, thickness_list, property_ids=[2,3,4,5,6]):
        # 1 -> all 5 shell thickness vary with same value. 
        # 2 -> only first and the last shell thickness vary, other fixed with middle value. 
        # 3 -> the first, middle, and last shell thickness vary. 
        # 4 -> expect for middle shell, the other 4 shell thickness vary.
        # 5 -> all three shell thickness vary.
        
        #self.write_shell_property(thickness_list, property_ids)
        self._write_combined_file()
        self._write_bc_wall()
        self._write_material()
        self._write_dcc()
        self._write_property()
        #self.merge_files("ThreePointBending_0000.rad", ["ThreePointBending_base.rad", "shell.rad"])
    
    def _write_combined_file(self):
        """
        Combines file together. combine is ready to be run via Radioss/OpenRadioss
        """
        adr = os.path.join(os.getcwd(),'ThreePointBending_0000.rad') 
        inf = open(adr, 'w')
        inf.write("#RADIOSS STARTER\n")
        inf.write("#--------------------------------------------------------------------------------------------------|\n")
        inf.write("#---1----|----2----|----3----|----4----|----5----|----6----|----7----|----8----|----9----|---10----|\n")
        inf.write("#--------------------------------------------------------------------------------------------------|\n")
        inf.write("/BEGIN\n")

        inf.write("COMBINE\n")
        #inf.write('{0:>10}{1:>10}\n'.format(2023,0))
        inf.write('      2023         0\n')
        inf.write("                  Mg                  mm                   s\n")
        inf.write("                  Mg                  mm                   s\n")

        inf.write('#------------------------------------------------------------------------------------|\n')
        inf.write('#---1----|----2----|----3----|----4----|----5----|----6----|----7----|----8----|----9----|---10----|\n')
        inf.write('#------------------------------------------------------------------------------------|\n')

        inf.write('/ANALY\n')
        inf.write('#    N2D3D              IPARITH      ISUB\n')
        inf.write('         0                   1         0\n')
        inf.write('#---1----|----2----|----3----|----4----|----5----|----6----|----7----|----8----|----9----|---10----|\n')
        inf.write('/DEF_SOLID\n')
        inf.write('#  I_SOLID    ISMSTR     ICPRE             ITETRA4  ITETRA10      IMAS    IFRAME\n')
        inf.write('         0         0         0                   0         0         0         0\n')
        inf.write('#---1----|----2----|----3----|----4----|----5----|----6----|----7----|----8----|----9----|---10----|\n')
        inf.write('/DEF_SHELL\n')
        inf.write('#  I_SHELL    ISMSTR   ICthick     Iplas   Istrain         -         -     Ish3n     Idril\n')
        inf.write('        24         2         1         1         1                             2         0\n')
        inf.write('#---1----|----2----|----3----|----4----|----5----|----6----|----7----|----8----|----9----|---10----|\n')
        inf.write('/IOFLAG\n')
        inf.write('#     IPRI                         IOUTP    IOUTYY   IROOTYY     IDROT\n')
        inf.write('         0                             0         0         0         0\n')
        inf.write('#---1----|----2----|----3----|----4----|----5----|----6----|----7----|----8----|----9----|---10----|\n')
        inf.write('/SPMD\n')
        inf.write('#   DOMDEC     Nproc              Dkword             Nthread\n')
        inf.write('         0         1                   0                   1\n')
        
    

        inf.write("#--------------------------------------------------------------------------------------------------|\n")
        inf.write("#---1----|----2----|----3----|----4----|----5----|----6----|----7----|----8----|----9----|---10----|\n")
        inf.write("#--------------------------------------------------------------------------------------------------|\n")
        inf.writelines('#include material.txt\n')
        inf.writelines('#include property.txt\n')
        inf.writelines('#include mesh.txt\n')
        inf.writelines('#include bc_wall.txt\n') 
        inf.writelines('#include dcc.txt\n')
        inf.write("#--------------------------------------------------------------------------------------------------|\n")
        inf.write("#---1----|----2----|----3----|----4----|----5----|----6----|----7----|----8----|----9----|---10----|\n")
        inf.write("#--------------------------------------------------------------------------------------------------|\n")
        inf.write('/END')   
        inf.close()

    def _write_bc_wall(self):
        # Writes the boundary condition for the cylinder impactor and constraints movement of
        # the clamped DoF

        adr = os.path.join(os.getcwd(),'bc_wall.txt') 
        inf = open(adr, 'w')
        inf.write("#--------------------------------------------------------------------------------------------------|\n")
        inf.write("#---1----|----2----|----3----|----4----|----5----|----6----|----7----|----8----|----9----|---10----|\n")
        inf.write("#--------------------------------------------------------------------------------------------------|\n")
        #inf.write("/BEGIN\n")
        #inf.write("BOUNDARY_CONDITIONS\n")
        #inf.write('{0:>10}{1:>10}\n'.format(2023,0))
        #inf.write('{0:>20}{1:>20}{2:>20}\n'.format("Mg","mm","s"))
        #inf.write('*KEYWORD\n')
        inf.write('#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$\n')
        inf.write('#                                 Rigid Wall                                  #\n')
        inf.write('#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$\n')
        inf.write('\n\n')
        inf.write('/NODE\n')
        inf.write('#    nid               x               y               z\n')
        inf.write("{:>10}{:>20}{:>20}{:>20}\n".format(str(self.wall_n_id),\
                '0.0', '0.0', str(self.wall_loc)))
        inf.write('/RWALL/CYL/1\n')
        inf.write('IMPACTOR\n')
        inf.write('#{:>9}{:>10}{:>10}{:>10}\n'.format("node_ID","Slide", "grnod_ID1", "grnod_ID2"))
        inf.writelines("{:>10}{:>10}{:>10}{:>10}\n".format(self.wall_n_id,'0','0', '0'))
        inf.write('#           D_search                fric            Diameter                ffac       ifq\n')
        inf.writelines("{:>20}{:>20}{:>20}{:>20}{:>10}\n".format(2*self.impactor_diameter,'1', self.impactor_diameter, '0','0'))
        inf.write('#               Mass                VX_0                VY_0                VZ_0\n')
        inf.writelines("{:>20}{:>20}{:>20}{:>20}\n".format(str(self.rigid_mass), '0.0', '0.0', self.wall_vel))
        inf.write('#               X_M1                Y_M1                Z_M1\n')
        inf.writelines("{:>20}{:>20}{:>20}\n".format('0.0', '100.0', str(self.wall_loc)))
        # inf.write('*RIGIDWALL_PLANAR_ID\n') 
        # inf.write('$#      id\n')          
        # inf.write('         2\n') 
        # inf.write('$#    nsid    nsidex     boxid    offset     birth     death     rwksf\n') 
        # inf.write('         0         0         0     0.000     0.0001.0000E+201.00000000\n')      
        # inf.write('$#      xt        yt        zt        xh        yh        zh      fric      wvel\n')    
        # inf.write('     0.000     0.000     0.000     0.000     0.00010.0000000     0.000     0.000\n')  
        inf.write('\n#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$\n')
        inf.write('#                                Boundary SPC                                 $\n')
        inf.write('#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$\n')
        inf.write('\n\n')  
        inf.write('/BCS/1\n')
        inf.write("LEFT_BC\n")
        inf.write('#  Tra rot   skew_ID  grnod_ID\n')
        inf.write("{:>6}{:>4}{:>10}{:>10}\n".format('111','101','0',101))
        inf.write('/BCS/2\n')
        inf.write("RIGHT_BC\n")
        inf.write('#  Tra rot   skew_ID  grnod_ID\n')
        inf.write("{:>6}{:>4}{:>10}{:>10}\n".format('111','101','0',102))
        
        # ---------- lowest node set id is 101
        # inf.writelines("{:>10}{:>10}{:>10}{:>10}{:>10}{:>10}{:>10}{:>10}\n".format(str(101),'0', '1'
        #             , '1', '1', '1', '1', '1'))
        # inf.write('$\n$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$')
        # inf.write('\n$                                For Output                                     $')
        # inf.write('$\n$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$')
        # inf.write('\n$\n')  
        # inf.write('*DATABASE_HISTORY_NODE\n')
        # inf.write('$#    nid1     nid2     nid3     nid4     nid5     nid6     nid7     nid8\n')
        # ---------- lowest node set id is 101
        # inf.writelines("{:>10}{:>10}{:>10}{:>10}{:>10}{:>10}{:>10}{:>10}\n".format(str(self.wall_n_id), 
        #                                     str(self.mesh.node_starting_id), '0', '0', '0', '0', '0', '0'))
        inf.write("#enddata")
        inf.write('\n/END')
        inf.close()
    
    def _write_material(self):
        adr = os.path.join(os.getcwd(),'material.txt')
        inf = open(adr, 'w')
        inf.write("#--------------------------------------------------------------------------------------------------|\n")
        inf.write("#---1----|----2----|----3----|----4----|----5----|----6----|----7----|----8----|----9----|---10----|\n")
        inf.write("#--------------------------------------------------------------------------------------------------|\n")
        #inf.write("/BEGIN\n")
        #inf.write("MATERIAL\n")
        #inf.write('{0:>10}{1:>10}\n'.format(2023,0))
        #inf.write('{0:>20}{1:>20}{2:>20}\n'.format("Mg","mm","s"))
        # inf.write("#\n/MAT/PLAS_JOHNS/1\n")
        # inf.write("Aluminum\n")
        # inf.write("#              RHO_I\n")
        # inf.write("{:>20}\n".format(str(self.material_density)))
        # inf.write("#                  E                  Nu     Iflag\n")
        # inf.write("{:>20}{:>20}{:>10}\n".format(str(self.mat_young_mod), str(self.mat_poisson_r), '1'))
        # inf.write("#              SIG_Y                 UTS                EUTS           EPS_p_max            SIG_max0\n")
        # #inf.write("{:>20}{:>20}{:>20}{:>20}{:>20}\n".format(300, 600, 0.2, '0.0', '0.0'))
        # inf.write("{:>20}{:>20}{:>20}{:>20}{:>20}\n".format(180, 248.5, 0.4, '0.0', '0.0'))
        # inf.write("#                  c           EPS_DOT_0       ICC   Fsmooth               F_cut               Chard\n")
        # inf.write("{:>20}{:>20}{:>10}{:>10}{:>20}{:>20}\n".format(0, 0, 0, 0, 0, 0))
        # inf.write("#                  m              T_melt              rhoC_p                 T_r\n")
        # inf.write("{:>20}{:>20}{:>20}{:>20}\n".format(0, 0, 0, 0))
        inf.write("/FUNCT/1\n")
        inf.write("Plasticity\n")
        inf.write("#{:>19}{:>20}\n".format('X', 'Y'))
        inf.write("{:>20}{:>20}\n".format(0, 180))
        inf.write("{:>20}{:>20}\n".format(.01, 190))
        inf.write("{:>20}{:>20}\n".format(.02, 197))
        inf.write("{:>20}{:>20}\n".format(.05, 211.5))
        inf.write("{:>20}{:>20}\n".format(.1, 225.8))
        inf.write("{:>20}{:>20}\n".format(.15, 233.6))
        inf.write("{:>20}{:>20}\n".format(.2, 238.5))
        inf.write("{:>20}{:>20}\n".format(.4, 248.5))
        inf.write("#--------------------------------------------------------------------------------------------------|\n")
        inf.write("#---1----|----2----|----3----|----4----|----5----|----6----|----7----|----8----|----9----|---10----|\n")
        inf.write("#--------------------------------------------------------------------------------------------------|\n")
        inf.write("#\n/MAT/COWPER/1\n")
        inf.write("Aluminum\n")
        inf.write("#\n")
        inf.write("#              RHO_I\n")
        inf.write("{:>20}\n".format(str(self.material_density)))
        inf.write("#                  E                  Nu  \n")
        inf.write("{:>20}{:>20}\n".format(str(self.mat_young_mod), str(self.mat_poisson_r)))
        inf.write("#{:>10}{:>20}{:>20}{:>20}{:>20}\n".format("a", "b","n", "C_hard", "sigma_max_0"))
        #inf.write("{:>20}{:>20}{:>20}{:>20}{:>20}\n".format(0, 0, 1.0, 1, str(1e+19)))
        inf.write("{:>20}{:>20}{:>20}{:>20}{:>20}\n".format(0, 0, 1.0, 1, 0))
        inf.write("#{:>10}{:>20}{:>10}{:>10}{:>20}{:>20}\n".format("c", "p","ICC", "F_smooth", "F_cut", "VP"))
        #inf.write("{:>20}{:>20}{:>10}{:>10}{:>20}{:>20}\n".format(0, 1.0, 1, 0, str(1e20), 2))
        inf.write("{:>20}{:>20}{:>10}{:>10}{:>20}{:>20}\n".format(0, 1.0, 1, 0, 0, 2))
        inf.write("#{:>10}{:>20}{:>20}\n".format("e_p^max", "e_t1","e_t2"))
        #inf.write("{:>20}{:>20}{:>20}\n".format(str(1e+20), str(1e+20),str(2*1e+20)))
        inf.write("{:>20}{:>20}{:>20}\n".format(str(0), str(0),str(0)))
        inf.write("#{:>9}{:>30}\n".format("fct_Idy","F_scale_y"))
        inf.write("{:>10}{:>30}\n".format(1,1.0))
        inf.write("#--------------------------------------------------------------------------------------------------|\n")
        inf.write("#---1----|----2----|----3----|----4----|----5----|----6----|----7----|----8----|----9----|---10----|\n")
        inf.write("#--------------------------------------------------------------------------------------------------|\n")
        inf.write("#\n")
        inf.write("#\n")
        inf.write("#--------------------------------------------------------------------------------------------------|\n")
        inf.write("#---1----|----2----|----3----|----4----|----5----|----6----|----7----|----8----|----9----|---10----|\n")
        inf.write("#--------------------------------------------------------------------------------------------------|\n")

        inf.write("#\n")
        inf.write("#--------------------------------------------------------------------------------------------------|\n")
        inf.write("#---1----|----2----|----3----|----4----|----5----|----6----|----7----|----8----|----9----|---10----|\n")
        inf.write("#--------------------------------------------------------------------------------------------------|\n")


        # inf.writelines("""
        # /MAT/PLAS_JOHNS/1
        # Steel_DP600
        # #              RHO_I
        #             7.85E-9                   0
        # #                  E                  Nu     Iflag
        #             210000                  .3         1
        # #              SIG_Y                 UTS                EUTS           EPS_p_max            SIG_max0
        #                 300                 600                  .2                   0                   0
        # #                  c           EPS_DOT_0       ICC   Fsmooth               F_cut               Chard
        #                 0                   0         0         0                   0                   0
        # #                  m              T_melt              rhoC_p                 T_r
        #                 0                   0                   0                   0
        #             """)
        inf.write("#enddata")
        inf.write('\n/END')

        inf.close()
    
    def _write_dcc(self):
        """
        writes Database, Control and Contact keywords 

        Output
            Creates adr.k includes control, contact and database 
        """   
        adr = os.path.join(os.getcwd(),'dcc.txt')

        inf = open(adr, 'w')
        #inf.write("/BEGIN\n")
        #inf.write("SIMULATION_CONTROL\n")
        #inf.write('{0:>10}{1:>10}\n'.format(2023,0))
        #inf.write('{0:>20}{1:>20}{2:>20}\n'.format("Mg","mm","s"))
        # ----------- Control
        inf.write('#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$\n')
        inf.write('#                                   Control                                   $\n')
        inf.write('#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$\n')
        inf.write("/TH/MODE/1\n")
        inf.write("intrusionTrackModes\n")
        inf.write("#     var1      var2      var3      var4      var5      var6      var7      var8      var9     var10\n")
        inf.write("       DEF\n")
        inf.write("#     Obj1      Obj2      Obj3      Obj4      Obj5      Obj6      Obj7      Obj8      Obj9     Obj10\n")
        inf.write("         1        \n")
        inf.write("#---1----|----2----|----3----|----4----|----5----|----6----|----7----|----8----|----9----|---10----|\n")

        inf.write("/TH/NODE/1\n")
        inf.write("intrusionTrack\n")
        inf.write("#     var1      var2      var3      var4      var5      var6      var7      var8      var9     var10\n")
        inf.write("         A         D         V         \n")
        inf.write("#{:>9}{:>10}{:>50}\n".format("NODid", 'Iskew', "NODname"))
        inf.write("{:>10}{:>10}{:>50}\n".format(self.wall_n_id, '0', "IntrusionNode"))
        inf.write("#---1----|----2----|----3----|----4----|----5----|----6----|----7----|----8----|----9----|---10----|\n")

        inf.write("/TH/RWALL/2\n")
        inf.write("TH_RWALL\n")
        inf.write("#     var1      var2      var3      var4      var5      var6      var7      var8      var9     var10\n")
        inf.write("       DEF       \n")
        inf.write("#     Obj1      Obj2      Obj3      Obj4      Obj5      Obj6      Obj7      Obj8      Obj9     Obj10\n")
        inf.write("         1\n")
        inf.write("#---1----|----2----|----3----|----4----|----5----|----6----|----7----|----8----|----9----|---10----|\n")

#/TH/RBODY/3
#TH_RBODY
#     var1      var2      var3      var4      var5      var6      var7      var8      var9     var10
#      DEF       
#     Obj1      Obj2      Obj3      Obj4      Obj5      Obj6      Obj7      Obj8      Obj9     Obj10
#         1         2
        # Write the contact details
        inf.write("#\n")
        inf.write("#---1----|----2----|----3----|----4----|----5----|----6----|----7----|----8----|----9----|---10----|\n")
        inf.write("#--------------------------------------------------------------------------------------------------|\n")
        inf.write("#\n")
        inf.write("/INTER/TYPE25/1\n")
        inf.write("self_contact\n")
        inf.write("# Surf_ID1  Surf_ID2      Istf      Ithe      Igap   Irem_i2                Idel     Iedge\n")
        inf.write("         4         0         0         0         0         0                   0         0\n")
        inf.write("# grnd_IDS                     Gap_scale          %mesh_size           Gap_max_s           Gap_max_m\n")
        inf.write("         0                             0                   0                   0                   0\n")
        inf.write("#              Stmin               Stmax     Igap0    Ishape          Edge_angle\n")
        inf.write("                   0                   0         0         0                   0\n")
        inf.write("#              Stfac                Fric           Tpressfit              Tstart               Tstop\n")
        inf.write("                   0                  .9                   0                   0                   0\n")
        inf.write("#      IBC               IVIS2    Inacti               ViscS    Ithick                          Pmax\n")
        inf.write("       000                   0         0                   0         0                             0\n")
        inf.write("#    Ifric    Ifiltr               Xfreq             sens_ID                                 fric_ID\n")
        inf.write("         0         0                   0                   0                                       0\n")
        inf.write("#enddata")
        inf.write('\n/END')
        inf.close()
    
    def _write_property(self):
        """
        writes the property of the material
        """
        adr = os.path.join(os.getcwd(),'property.txt') 
        inf = open(adr, 'w')
        inf.write("#--------------------------------------------------------------------------------------------------|\n")
        inf.write("#---1----|----2----|----3----|----4----|----5----|----6----|----7----|----8----|----9----|---10----|\n")
        inf.write("#--------------------------------------------------------------------------------------------------|\n")
        #inf.write("/BEGIN\n")
        #inf.write("PROPERTY\n")
        #inf.write('{0:>10}{1:>10}\n'.format(2023,0))
        #inf.write('{0:>20}{1:>20}{2:>20}\n'.format("Mg","mm","s"))
        inf.write("#--------------------------------------------------------------------------------------------------|\n")
        inf.write("#---1----|----2----|----3----|----4----|----5----|----6----|----7----|----8----|----9----|---10----|\n")
        inf.write("#--------------------------------------------------------------------------------------------------|\n")

        inf.write("/PROP/SHELL/1\n")
        inf.write("PID_tube\n")
        inf.write("#   Ishell    Ismstr     Ish3n    Idrill                            P_thick_fail\n")
        inf.write("        24         1         2         2                                       0\n")
        inf.write("#                 hm                  hf                  hr                  dm                  dn\n")
        inf.write("                   0                   0                   0                   0                   0\n")
        inf.write("#        N   Istrain               Thick              Ashear              Ithick     Iplas\n")
        inf.write("         5         1                 1.8                   0                   1         1\n")

        inf.write("#---1----|----2----|----3----|----4----|----5----|----6----|----7----|----8----|----9----|---10----|\n")
        inf.write("/PROP/SHELL/2\n")
        inf.write("PID_h\n")
        inf.write("#   Ishell    Ismstr     Ish3n    Idrill                            P_thick_fail\n")
        inf.write("        24         1         2         2                                       0\n")
        inf.write("#                 hm                  hf                  hr                  dm                  dn\n")
        inf.write("                   0                   0                   0                   0                   0\n")
        inf.write("#        N   Istrain               Thick              Ashear              Ithick     Iplas\n")
        inf.write("         5         1                 0.7                   0                   1         1\n")
        inf.write("#---1----|----2----|----3----|----4----|----5----|----6----|----7----|----8----|----9----|---10----|\n")

        inf.write("/PROP/SHELL/3\n")
        inf.write("PID_V\n")
        inf.write("#   Ishell    Ismstr     Ish3n    Idrill                            P_thick_fail\n")
        inf.write("        24         1         2         2                                       0\n")
        inf.write("#                 hm                  hf                  hr                  dm                  dn\n")
        inf.write("                   0                   0                   0                   0                   0\n")
        inf.write("#        N    Istrain              Thick              Ashear              Ithick     Iplas\n")
        inf.write("{:>10}{:>10}{:>20}{:>20.5f}{:>20}{:>10}\n".format(5, 1, 1.8, 5/6, 1,1))
        #inf.write(f"        5                           1.8             {0.8336}                  1         1\n")
        inf.write("#---1----|----2----|----3----|----4----|----5----|----6----|----7----|----8----|----9----|---10----|\n")


        inf.write("#enddata")
        inf.write('\n/END')
        inf.close()