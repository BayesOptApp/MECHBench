import os
import numpy as np

from src.sob.physical_models.fem_settings.abstractFEMSettings import AbstractFEMSettings
from src.sob.physical_models.meshes import StarBoxMesh




class StarBoxModel(AbstractFEMSettings):
    def __init__(self, mesh:StarBoxMesh, **kwargs) -> None:
        self.mesh = mesh
        self.units = mesh.units
        mesh.write_mesh_file()
        # mesh.units =  '  kg  mm  ms  kN  GPa  kN-mm'
        
        # Define default parameters for load_impactor
        self.impactor_defaults = {
            'wall_n_id': 999999,
            'wall_loc': self.mesh.extrusion_length+1,
            'wall_mass': 250.0,
            'wall_vel': 7.0,
        }
        self.rigid_mass = self.impactor_defaults['wall_mass']
        # Define default parameters for load_database
        self.database_defaults = {
            'end_time': 45.0,
            'database_dtime': 0.05,
            'consider_d3plot': True, # cons_d3plot
            'd3plot_dtime': 0.5, # d3plot_dtime
            'cons_d3thdt' : False,
            'd3thdt_dtime': 10000, # d3thdt_dt
            'shell_warping': 1,
            'binary_ascii': 2,
            'write_cshell':True,
            'tdel': 0.0
            # Add other default parameters here
        }
        
        self.shell_warping = 1    # BWC, lsdyna default is 2. if there is warping set it to 1
        self.binary_ascii = 2  # 1: only ascii   2: only binary   3: both ascii and binary
        self.write_cshell = True   
        self.tdel = 0.0

        self.intout = 'STRESS'
        self.nodout = 'STRESS' 

        self.intp_db = 3
        self.sigflg = 1
        self.epsflg = 1
        self.rltflg = 1

        # Define default parameters for load_material
        self.material_defaults = {
            'mat_id': 999,
            'mat_density': 7.83E-6,
            'mat_young_mod': 200.0,
            'mat_poisson_r': 0.3,
            'mat_yield_initial': 0.366,
            'mat_tang_mod': 0.0,
            'mat_failure_pstrain': 1.0E+21,
            'mat_cowper_symond_c': 40.0,
            'mat_cowper_symond_p': 5.0,
            'mat_vp_rate_efffect': 1,
            'mat_load_curve_id': 1,
            'mat_effective_plastic_strain_stress': np.array([[0., 0.366], [2.5e-2, 0.4240],
                                                             [4.9e-2, 0.476], [7.2e-2, 0.507],
                                                             [9.5e-2, 0.529], [0.118, 0.546],
                                                             [0.140, 0.559], [0.182, 0.584]]),
            # Add other default parameters here
        }

        self._load_impactor(**kwargs)
        self._load_database(**kwargs)
        self._load_material(**kwargs)

        # ------- nodal forces
        self.write_nod_force_top = True
        self.write_nod_force_bottom = True

    def mass(self):
        return self.mesh.volume()*self.mat_density
    
    def absorbed_energy(self):
        # initial kinetic energy
        return self.wall_mass*self.wall_vel**2/2
    
    @property
    def material_card_type(self):
        return "LS_Dyna"
    
    @property
    def impactor_offset(self)->float:
        return 1.00
    
    def _load_impactor(self, **kwargs):
        for key, value in kwargs.items():
            setattr(self, key, value)
        for key, value in self.impactor_defaults.items():
            setattr(self, key, kwargs.get(key, value))

    def _load_database(self, **kwargs):
        for key, value in kwargs.items():
            setattr(self, key, value)
        for key, value in self.database_defaults.items():
            setattr(self, key, kwargs.get(key, value))

    def _load_material(self, **kwargs):
        for key, value in kwargs.items():
            setattr(self, key, value)
        for key, value in self.material_defaults.items():
            setattr(self, key, kwargs.get(key, value))

    def _write_bc_wall(self):
        # ----------------------------------------------------------- rigid walls
        adr = os.path.join(os.getcwd(),'bc_wall.k') 
        inf = open(adr, 'w')
        inf.write('*KEYWORD\n')
        inf.write('$\n$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$')
        inf.write('\n$                                 Rigid Wall                                  $')
        inf.write('$\n$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$')
        inf.write('\n$\n')
        inf.write('*NODE\n')
        inf.write('$    nid               x               y               z      tc      rc\n')
        inf.write("{:>8}{:>16}{:>16}{:>16}{:>8}{:>8}\n".format(str(self.wall_n_id),\
                '0.0', '0.0', str(self.wall_loc), '0', '0'))
        # inf.write('*RIGIDWALL_PLANAR_MOVING_FORCES_ID\n')
        # #inf.write('*RIGIDWALL_PLANAR_MOVING_ID\n')
        # inf.write('$#      id\n')
        # inf.write("{:>10}\n".format('1'))
        # inf.write('$#    nsid    nsidex     boxid    offset     birth     death     rwksf  \n')
        # # inf.writelines("{:>10}{:>10}{:>10}{:>10}{:>10}{:>10}{:>10}\n".format('0','0', \
        # #     '0', '0.0', '0.0', '1.00E20', '1.0'))
        # inf.writelines("{:>10}{:>10}{:>10}{:>10}{:>10}{:>10}{:>10}\n".format(self.wall_n_id,'0', \
        #      '0', '1.00E20', '0.0', '1.00E20', '1.0'))
        # inf.write('$#      xt        yt        zt        xh        yh        zh      fric      wvel\n')
        # inf.writelines("{:>10}{:>10}{:>10}{:>10}{:>10}{:>10}{:>10}{:>10}\n".format('0.0','0.0', \
        #     str(self.wall_loc), '0.0', '0.0', '0.0', '1.0', '0.0'))
        # inf.write('$#    mass        v0\n')
        # inf.write("{:>10}{:>10}\n".format(str(self.wall_mass),str(self.wall_vel)))    
        # inf.write('$#    soft      ssid        n1        n2        n3        n4\n')
        # inf.writelines("{:>10}{:>10}{:>10}{:>10}{:>10}{:>10}\n".format('0','0', \
        #     str(self.wall_n_id), '0', '0', '0'))
        # inf.write('*RIGIDWALL_PLANAR_ID\n') 
        # #inf.write('*RIGIDWALL_PLANAR_FORCES_ID\n') 
        # inf.write('$#      id\n')          
        # inf.write("{:>10}\n".format('2')) 
        # inf.write('$#    nsid    nsidex     boxid    offset     birth     death     rwksf\n') 
        # inf.writelines("{:>10}{:>10}{:>10}{:>10}{:>10}{:>10}{:>10}\n".format('0','0', \
        #     '0', '0.0', '0.0', '1.00E20', '1.0'))      
        # inf.write('$#      xt        yt        zt        xh        yh        zh      fric      wvel\n')    
        # inf.writelines("{:>10}{:>10}{:>10}{:>10}{:>10}{:>10}{:>10}{:>10}\n".format('0.0','0.0', \
        #     0.0, '0.0', '0.0', '1.0', '1.0', '0.0'))  
        inf.write('$\n$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$')
        inf.write('\n$                                Boundary SPC                                 $')
        inf.write('$\n$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$')
        inf.write('\n$\n')  
        inf.write('*BOUNDARY_SPC_SET\n')
        inf.write('$#    nsid       cid      dofx      dofy      dofz     dofrx     dofry     dofrz\n')
        # ---------- lowest node set id is 101
        # inf.writelines("{:>10}{:>10}{:>10}{:>10}{:>10}{:>10}{:>10}{:>10}\n".format(str(101),'0', '1'
        #             , '1', '1', '1', '1', '1'))
        inf.writelines("{:>10}{:>10}{:>10}{:>10}{:>10}{:>10}{:>10}{:>10}\n".format(str(101),'0', '1'
                     , '1', '1', '1', '1', '1'))
        inf.write('$\n$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$')
        inf.write('\n$                                For Output                                     $')
        inf.write('$\n$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$')
        inf.write('\n$\n')  
        inf.write('*DATABASE_HISTORY_NODE\n')
        inf.write('$#    nid1     nid2     nid3     nid4     nid5     nid6     nid7     nid8\n')
        # ---------- lowest node set id is 101
        inf.writelines("{:>10}{:>10}{:>10}{:>10}{:>10}{:>10}{:>10}{:>10}\n".format(str(self.wall_n_id), 
                                            str(self.mesh.node_starting_id), '0', '0', '0', '0', '0', '0'))
        inf.write('$\n*END')
        inf.close()

    def _write_dcc(self):
        """
        writes Database, Control and Contact keywords 

        Inputs
            #adr          # addresss of the to be created output file
            endtim        # endtime
            database_dtime   # dt for all databases except than d3polt and d3thdt  
            d3plot_dtime     # dt foe d3plot
            d3thdt_dt     # dt for d3thdt
            binary_ascii  # output format which can be binary, ascii or both
            write_cshell  # writing control_shell or not
            # --- extendend binary database inputs
            intp_db      # maxint in extended binary
            sigflg       # flag for including (1) stress tensor in shell database
            epsflg       # flag for including (1) effective plastic strain in shell database
            rltflg       # flag for including (1) resultant stress in shell database
            intout       # output stress/strain at all integration points in ELOUT
            nodout       # output stress/strain at connectivity point in ELOUT
            CONTROL_ACCURACY is not set as material is not direction dependent. otherwise
            add CONTROL_ACCURACY and set INN to appropriate value
        Output
            creted adr.k includes control, contact and database 
        """   
        adr = os.path.join(os.getcwd(),'dcc.k')
        inf = open(adr, 'w')
        inf.write('*KEYWORD\n')
        # ----------- Control
        inf.write('$\n$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$')
        inf.write('\n$                                   Control                                   $')
        inf.write('$\n$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$')
        inf.write('\n$\n')
        if  self.write_cshell == True:
            inf.writelines('*CONTROL_SHELL\n')
            inf.writelines('$#  wrpang     esort     irnxx    istupd    theory       bwc     miter      proj\n')
            inf.writelines("{:>10}{:>10}{:>10}{:>10}{:>10}{:>10}{:>10}{:>10}\n".format('20.0', 
                        '1', '-1', '0', '2', str(int(self.shell_warping)), '1', '0'))
            inf.writelines('$# rotascl    intgrd    lamsht    cstyp6    tshell\n')
            inf.writelines("{:>10}{:>10}{:>10}{:>10}{:>10}\n".format('1.0', '0', '0', '1', '0'))
            inf.writelines('$# psstupd   sidt4tu     cntco    itsflg    irquad \n')
            inf.writelines("{:>10}{:>10}{:>10}{:>10}{:>10}\n".format('0', '0', '0', '0', '2'))
            inf.writelines('$#  nfail1    nfail4   psnfail    keepcs     delfr   drcpsid    drcprm \n')
            inf.writelines("{:>10}{:>10}{:>10}{:>10}{:>10}{:>10}{:>10}\n".format('1', 
                        '0', '0', '0', '0', '0', '1.0'))    
        else:
            print ('\nCONTROL_SHELL was not written')
        
        inf.writelines('*CONTROL_ENERGY\n')
        inf.writelines('$     hgen      rwen    slnten     rylen\n')
        inf.writelines('         2         2         1         1\n')

        ### ++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        ### CONTROL FOR PARALLEL ARITHMETIC
        ### ++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        inf.writelines('*CONTROL_PARALLEL\n')
        inf.writelines('$    -----     -----     CONST     -----\n')
        inf.writelines('         0         0         1         1\n')
        
        ### ++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        ###
        ### ++++++++++++++++++++++++++++++++++++++++++++++++++++++++

        inf.writelines('*CONTROL_TERMINATION\n')  
        inf.writelines('$   endtim    endcyc     dtmin    endeng    endmas\n')
        inf.writelines("{:>10}{:>10}{:>10}{:>10}{:>10}\n".format(str(self.end_time), '0', '0.0', '0.0', '1.0E8'))
        
        # ----------- Database
        inf.write('$\n$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$')
        inf.write('\n$                                  Database                                   $')
        inf.write('$\n$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$')
        inf.write('\n$\n')
        bin_asc = str(int(self.binary_ascii))

        # inf.writelines('*DATABASE_FORMAT\n')            
        # inf.writelines('$#      FORMAT N\n') 
        # inf.writelines("{:>10}\n".format("", "1"))


        # inf.writelines('*DATABASE_GLSTAT\n')            
        # inf.writelines('$#      dt    binary      lcur     ioopt\n')
        # inf.writelines("{:>10}{:>10}{:>10}{:>10}\n".format(str(self.database_dtime), bin_asc, '0', '1'))

        # inf.writelines('*DATABASE_MATSUM\n')            
        # inf.writelines('$#      dt    binary      lcur     ioopt\n')
        # inf.writelines("{:>10}{:>10}{:>10}{:>10}\n".format(str(self.database_dtime), bin_asc, '0', '1')) 

        # !! Cannot recogenized by OpenRadioss
        # !! nodfor gives time histories of contact forces at nodes.
        # inf.writelines('*DATABASE_NODFOR\n')            
        # inf.writelines('$#      dt    binary      lcur     ioopt\n')
        # inf.writelines("{:>10}{:>10}{:>10}{:>10}\n".format(str(self.database_dtime), bin_asc, '0', '1')) 
        
        ### ++++++++++++++++++++++++++++++++++++++
        ### DATABASE CROSS SECTION_SET
        ### ++++++++++++++++++++++++++++++++++++++
        # inf.writelines('*DATABASE_CROSS_SECTION_SET_ID\n')            
        # inf.writelines('$#      csid    title \n')
        # inf.writelines("{:>10}{:>20}\n".format('1', "TOP PLANE"))
        # inf.writelines('$#      nsid    hsid    bsid    ssid   tsid     dsid \n')
        # inf.writelines("{:>10}\n".format('101'))

        # inf.writelines('*DATABASE_BINARY_INTFOR\n')            
        # inf.writelines('$#      dt\n')
        # inf.writelines("{:>10}\n".format(str(self.database_dtime)))

        ### ++++++++++++++++++++++++++++++++++++++
        ### DATABASE CROSS SECTION_SET
        ### ++++++++++++++++++++++++++++++++++++++

        inf.writelines('*DATABASE_NODOUT\n')            
        inf.writelines('$#      dt    binary      lcur     ioopt   option1   option2 \n')
        inf.writelines("{:>10}{:>10}{:>10}{:>10}{:>10}{:>10}\n".format(str(self.database_dtime), bin_asc, '0', '1','0.0','0'))

        inf.writelines('*DATABASE_RWFORC\n')            
        inf.writelines('$#      dt    binary      lcur     ioopt\n')
        inf.writelines("{:>10}{:>10}{:>10}{:>10}\n".format(str(self.database_dtime), bin_asc, '0', '1'))

        # inf.writelines('*DATABASE_RCFORC\n')            
        # inf.writelines('$#      dt    binary      lcur     ioopt\n')
        # inf.writelines("{:>10}{:>10}{:>10}{:>10}\n".format(str(self.database_dtime), bin_asc, '0', '1'))  
        # 
        # inf.writelines('*DATABASE_SECFORC\n')            
        # inf.writelines('$#      dt    binary      lcur     ioopt\n')
        # inf.writelines("{:>10}{:>10}{:>10}{:>10}\n".format(str(self.database_dtime), bin_asc, '0', '1'))

        inf.writelines('*DATABASE_RCFORC\n')            
        inf.writelines('$#      dt    binary      lcur     ioopt\n')
        inf.writelines("{:>10}{:>10}{:>10}{:>10}\n".format(str(self.database_dtime), bin_asc, '0', '1'))

        inf.writelines('*DATABASE_ELOUT\n')            
        inf.writelines('$#      dt    binary      lcur     ioopt\n')
        inf.writelines("{:>10}{:>10}{:>10}{:>10}\n".format(str(self.database_dtime), bin_asc, '0', '1'))

        # inf.writelines('*DATABASE_SPCFORC\n')            
        # inf.writelines('$#      dt    binary      lcur     ioopt\n')
        # inf.writelines("{:>10}{:>10}{:>10}{:>10}\n".format(str(self.database_dtime), bin_asc, '0', '1'))

        if self.consider_d3plot == True:
            inf.writelines('*DATABASE_BINARY_D3PLOT\n')    
            inf.writelines('$#      dt      lcdt      beam     npltc    psetid\n')    
            inf.writelines("{:>10}{:>10}{:>10}{:>10}{:>10}\n".format(str(self.d3plot_dtime), '0', '0', '0', '0'))    
            inf.writelines('$#   ioopt\n');inf.writelines('         0\n')              

        if self.cons_d3thdt == True:
            inf.writelines('*DATABASE_BINARY_D3THDT\n')    
            inf.writelines('$#      dt      lcdt      beam     npltc    psetid\n')    
            inf.writelines("{:>10}{:>10}{:>10}{:>10}{:>10}\n".format(str(self.d3thdt_dtime), '0', '0', '0', '0'))    

        inf.writelines('*DATABASE_EXTENT_BINARY\n') 
        inf.writelines('$#   neiph     neips    maxint    strflg    sigflg    epsflg    rltflg    engflg\n')
        inf.writelines("{:>10}{:>10}{:>10}{:>10}{:>10}{:>10}{:>10}{:>10}\n".format('0','0', \
        str(self.intp_db), '0', str(self.sigflg), str(self.epsflg), str(self.rltflg),'1'))    
        inf.writelines('$#  cmpflg    ieverp    beamip     dcomp      shge     stssz    n3thdt   ialemat\n') 
        inf.writelines("{:>10}{:>10}{:>10}{:>10}{:>10}{:>10}{:>10}{:>10}\n".format('0','0', \
        '0', '1', '1', '1', '2','1')) 
        inf.writelines('$# nintsld   pkp_sen      sclp     hydro     msscl     therm    intout    nodout\n')
        inf.writelines("{:>10}{:>10}{:>10}{:>10}{:>10}{:>10}{:>10}{:>10}\n".format('0','0', \
        '1.0', '0', '0', '0', self.intout, self.nodout))
        inf.writelines('$#    dtdt    resplt     neipb\n')    
        inf.writelines("{:>10}{:>10}{:>10}\n".format('0', '0', '0'))
        # !! Cannot recogenized by OpenRadioss
        # !! nodfor gives time histories of contact forces at nodes.
        # inf.writelines('*DATABASE_MASSOUT\n') 
        # inf.writelines('$#   setid     ndflg     rbflg\n')     
        # inf.writelines("{:>10}{:>10}{:>10}\n".format('0', '1', '0')) 
        
    # ----------- Contact
        inf.write('$\n$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$')
        inf.write('\n$                                   Contact                                   $')
        inf.write('$\n$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$')
        inf.write('\n$\n')
        inf.writelines('*CONTACT_AUTOMATIC_SINGLE_SURFACE\n')     
        inf.writelines('$#    ssid      msid     sstyp     mstyp    sboxid    mboxid       spr       mpr\n')       
        inf.writelines("{:>10}{:>10}{:>10}{:>10}{:>10}{:>10}{:>10}{:>10}\n".format('0','0', \
        '0', '0', '0', '0', '0','0')) 
        inf.writelines('$#      fs        fd        dc        vc       vdc    penchk        bt        dt\n')       
        inf.writelines("{:>10}{:>10}{:>10}{:>10}{:>10}{:>10}{:>10}{:>10}\n".format('0.08','0.8', \
        '0.0', '0.0', '0.0', '0', '0.0','1.0E20'))    
        inf.writelines('$#     sfs       sfm       sst       mst      sfst      sfmt       fsf       vsf\n')
        inf.writelines("{:>10}{:>10}{:>10}{:>10}{:>10}{:>10}{:>10}{:>10}\n".format('0.0','1.0', \
        '0.0', '0.0', '1.0', '1.0', '1.0','1.0')) 
        inf.write('$\n')
        inf.write('*END')   
        inf.close()

    def _write_material(self):
        adr = os.path.join(os.getcwd(),'material.k')
        inf = open(adr, 'w')
        inf.write('$#  units:' + self.units + '\n')
        inf.write('*KEYWORD\n')    
        inf.write('*MAT_PIECEWISE_LINEAR_PLASTICITY\n') 
        inf.write('$#     mid        ro         e        pr      sigy      etan      fail      tdel\n')
        inf.write("{:>10}{:>10}{:>10}{:>10}{:>10}{:>10}{:>10}{:>10}\n".format(str(int(self.mat_id)),\
                    str(self.mat_density), str(self.mat_young_mod), str(self.mat_poisson_r), str(self.mat_yield_initial), 
                    str(self.mat_tang_mod), str(self.mat_failure_pstrain), str(self.tdel)))
        inf.write('$#       c         p      lcss      lcsr        vp\n') 
        inf.write("{:>10}{:>10}{:>10}{:>10}{:>10}\n".format(str(int(self.mat_cowper_symond_c)),\
                str(self.mat_cowper_symond_p), str(int(self.mat_load_curve_id)), '0', 
                str(int(self.mat_vp_rate_efffect))))
        inf.write('$#    eps1      eps2      eps3      eps4      eps5      eps6      eps7      eps8\n')
        inf.write("{:>10}{:>10}{:>10}{:>10}{:>10}{:>10}{:>10}{:>10}\n".format('0.0','0.0', \
                '0.0', '0.0', '0.0', '0.0', '0.0','0.0'))  
        inf.write('$#     es1       es2       es3       es4       es5       es6       es7       es8\n')
        inf.write("{:>10}{:>10}{:>10}{:>10}{:>10}{:>10}{:>10}{:>10}\n".format('0.0','0.0', \
                '0.0', '0.0', '0.0', '0.0', '0.0','0.0'))     
        inf.write('*DEFINE_CURVE\n') 
        inf.write('$#    lcid      sidr       sfa       sfo      offa      offo    dattyp\n') 
        inf.write("{:>10}{:>10}{:>10}{:>10}{:>10}{:>10}{:>10}\n".format(str(int(self.mat_load_curve_id)),
                '0', '1.0', '1.0', '0.0', '0.0', '0'))   
        inf.write('$#                a1                  o1\n')
        for i in range(np.size(self.mat_effective_plastic_strain_stress,0)):
            inf.write("{:>20}{:>20}\n".format(str(self.mat_effective_plastic_strain_stress[i,0]),
                    str(self.mat_effective_plastic_strain_stress[i,1])))     
        inf.write('*END')   
        inf.close()

    def _write_nodal_force_top(self):
        """
        writes set of node list to calculate the external nodal forces on the top
        of the structure. In py_mesh node set at the bottom have id = 102, this will be
        the id of the node set for which the nodal force group will be calculated.
        Inputs                          
        """
        adr = os.path.join(os.getcwd(),'nodal_force_top.k') 
        inf = open(adr, 'w')
        inf.write('*KEYWORD\n')
        if self.write_nod_force_top == True:
            inf.write('$\n$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$')
            inf.write('\n$                           Nodal forces at top                                $')
            inf.write('\n$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$\n')   
            inf.write('$#\n')
            added_id = 101
            inf.write('**DATABASE_CROSS_SECTION_SET_ID\n')
            inf.write('$#      csid    title \n')
            inf.write("{:>10}{:>70}\n".format(str(added_id), "TOP PLANE"))
            inf.write('$#      nsid    hsid    bsid    ssid   tsid     dsid \n')
            inf.write("{0:>10}{1:>10}{1:>10}{1:>10}{1:>10}{1:>10}\n".format(str(added_id),0))
            #inf.write('*DATABASE_NODAL_FORCE_GROUP\n')
            #inf.write('$#    nsid       cid\n')
            #inf.write("{:>10}{:>10}\n".format(str(added_id),'0'))
            #inf.write('$\n')
            # 
            # inf.write('*DATABASE_HISTORY_NODE_SET\n')  
            # inf.write('$#    nsid       cid\n')
            # inf.write("{:>10}{:>10}\n".format(str(added_id),'0'))
            # inf.write('$\n')

        inf.write('*END')   
        inf.close()       

    def _write_nodal_force_bottom(self):
        """
        writes set of node list to calculate the external nodal forces on the bottom
        of the structure. In py_mesh node set at the bottom have id = 102, this will be
        the id of the node set for which the nodal force group will be calculated.
        Inputs                          
        """
        adr = os.path.join(os.getcwd(),'nodal_force_bottom.k') 
        inf = open(adr, 'w')
        inf.write('*KEYWORD\n')
        if self.write_nod_force_bottom == True:
            inf.write('$\n$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$')
            inf.write('\n$                          Nodal forces at bottom                              $')
            inf.write('\n$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$\n')   
            inf.write('$#\n')
            added_id = 102
            inf.write('*DATABASE_NODAL_FORCE_GROUP\n')
            inf.write('$#    nsid       cid\n')
            inf.write("{:>10}{:>10}\n".format(str(added_id),'0'))        
            inf.write('$\n')                   
        inf.write('*END')   
        inf.close()
    
    def _write_radioss_output_control(self)->None:
        r""" This is a trial function to write a
        Radioss themed output"""

        adr = os.path.join(os.getcwd(),'radioss_control_output.rad') 
        inf = open(adr, 'w')

        inf.write('###################################################################################')
        inf.write('\n#                          TRIAL RWALL DEFINITION                              $')
        inf.write('\n####################################################################################\n')   
        inf.write('##\n')

        inf.write('/RWALL/PLANE/1\n')
        inf.write('IMPACTOR\n')
        inf.write('#{:>10}{:>10}{:>10}{:>10}\n'.format("node_ID","Slide", "grnod_ID1", "grnod_ID2"))
        inf.writelines("{:>10}{:>10}{:>10}{:>10}\n".format(self.wall_n_id,'0','0', '0'))
        inf.write('#           D_search                fric            Diameter                ffac       ifq\n')
        inf.writelines("{:>20}{:>20}{:>20}{:>20}{:>10}\n".format('200','1', '150', '0','0'))
        inf.write('#               Mass                VX_0                VY_0                VZ_0\n')
        inf.writelines("{:>20}{:>20}{:>20}{:>20}\n".format(str(self.wall_mass), '0.0', '0.0', -self.wall_vel))
        inf.write('#               X_M1                Y_M1                Z_M1\n')
        inf.writelines("{:>20}{:>20}{:>20}\n".format('0.0', 0.0, 0.0))

        inf.write('/RWALL/PLANE/2\n')
        inf.write('GROUND\n')
        inf.write('#{:>10}{:>10}{:>10}{:>10}\n'.format("node_ID","Slide", "grnod_ID1", "grnod_ID2"))
        inf.writelines("{:>10}{:>10}{:>10}{:>10}\n".format('0','0',101, '0'))
        inf.write('#           D_search                fric            Diameter                ffac       ifq\n')
        inf.writelines("{:>20}{:>20}{:>20}{:>20}{:>10}\n".format('200','1', '150', '0','0'))
        inf.write('#               X_M                Y_M                Z_M\n')
        inf.writelines("{:>20}{:>20}{:>20}\n".format('0.0', 0.0, 0.0))
        inf.write('#               X_M1                Y_M1                Z_M1\n')
        inf.writelines("{:>20}{:>20}{:>20}\n".format('0.0', 0.0, 1.0))
        #### SOME WRITINGS ABOUT THE SECTION
        # inf.write(f'/TH/SECTIO/{101}\n')
        # inf.write('##   thgroup_name\n')
        # inf.write('BOTTOM_NODES_FORCES\n')
        # inf.write("{:>10}{:>10}{:>10}{:>10}{:>10}{:>10}\n".format(
        #     'FNX', 'FNY', 'FNZ', 'FTX', 'FTY', 'FTZ'))  
        # inf.write("{:>10}{:>10}{:>10}{:>10}{:>10}{:>10}\n".format(
        #     10001, 10002, 10003, 10004, 10005, 10006))
        # # GET INFORMATION ABOUT RIGID WALL
        # inf.write(f'/TH/RWALL/{100}\n')
        # inf.write('##   thgroup_name\n')
        # inf.write('DA_RIGID_WALL\n')
        # inf.write("{:>10}{:>10}{:>10}{:>10}{:>10}{:>10}\n".format(
        #     'FNX', 'FNY', 'FNZ', 'FTX', 'FTY', 'FTZ'))  
        # inf.write('#\n')
        # inf.write(f'/TH/RWALL/{101}\n')
        # inf.write('##   thgroup_name\n')
        # inf.write('DA_RIGID_WALL_2\n')
        # inf.write("{:>10}{:>10}{:>10}{:>10}{:>10}{:>10}\n".format(
        #     'FNX', 'FNY', 'FNZ', 'FTX', 'FTY', 'FTZ'))  
        inf.write('#\n')
        inf.write("#enddata\n")              
        inf.write('/END\n')   
        inf.close()


    def write_input_files(self):
        self._write_bc_wall()
        self._write_dcc()
        self._write_nodal_force_top()
        self._write_nodal_force_bottom() # Cannot recogenized by OpenRadioss
        self._write_material()
        self._write_radioss_output_control()
        """
        Combines file together. combine is ready to be run via LS_Dyna
        """
        adr = os.path.join(os.getcwd(),'combine.k') 
        inf = open(adr, 'w')
        inf.write('*KEYWORD\n')
        inf.write('$ UNITS\n')
        inf.write("*CONTROL_UNITS\n")
        inf.write('{:>10}{:>10}{:>10}\n'.format('mm','ms','kg'))
        inf.writelines('*INCLUDE\n')
        inf.writelines('mesh.k\n')
        inf.writelines('material.k\n')
        inf.writelines('bc_wall.k\n')
        inf.writelines('*INCLUDE_RADIOSS\n')
        inf.writelines('radioss_control_output.rad\n')
        inf.writelines('*INCLUDE\n')
        inf.writelines('dcc.k\n')
        inf.writelines('nodal_force_top.k\n')
        
        

        
        
        # inf.writelines('nodal_force_bottom.k\n')
        inf.write('*END')   
        inf.close()