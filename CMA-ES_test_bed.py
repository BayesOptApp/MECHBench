"""
This script works by using the CMA-ES library to optimize the StarBox Problem 

"""

from src import sob
# Import the basic libraries
import os
import numpy as np
import csv
import pandas as pd
import multiprocessing as mp
import shutil
import time
from typing import List, Callable, Iterable

# Import the platform library
import platform

# Import IOH extended library
import ioh

# Import the CMA-ES Library
import cma

# Import the Modular CMA-ES library (Jacob de Nobel, Diederick Vermetten)
from modcma import c_maes
from modcma import AskTellCMAES

##### -----------------------------------------------------------
##### ------------------CONSTANTS--------------------------------
##### -----------------------------------------------------------


# THIS IS TO SET THE LAMBDA TO ADJUST THE RESTRICTION ON THE PROBLEM
LAMDA:float = 3.981071705534969283e+02
INTRUSION_PRIME:float = 60.00


# These are the default parameters for CMA-ES
BUDGET:int = 5 # Manage a budget of 
DIMENSIONS:list = [1,3,5]
CURDIM:int = DIMENSIONS[0]
SIGMA_0:float = (0.25/1)*10
N_RUNS:int =  1 # Number of runs
MAX_RESTARTS_DEFAULT:int = 5

# DEFAULT LOGGER CONSTANTS
DEFAULT_TRIGGERS:list = [ioh.logger.trigger.ALWAYS]

# LIBRARY TO USE
# If set to '1' then use Niko Hansen's Library
# If '2', then use Modular CMA-ES library (Jacob de Nobel, Diederick Vermetten)
CMA_ES_HAND:int = 2

# Change this to define the case to evaluate the population size
POP_CASE = 1


class Starbox_problem(ioh.problem.RealSingleObjective):

    # This is a handle to be modified
    linux_system = not (platform.system() == 'Windows')
    if linux_system:
        #batch_file_path = "/media/feifan/TSHIBA/12_GitHub/OpenRadioss/linux_scripts_mk3/openradioss_run_script_ps.sh"
        batch_file_path = "/home/ivanolar/Documents/OpenRadioss2/OpenRadioss_linux64/OpenRadioss/Tools/openradioss_gui/runopenradioss.py"
        
    else:
        #batch_file_path = "D:/OpenRadioss/win_scripts_mk3/openradioss_run_script_ps.bat"
        batch_file_path = "D:/OpenRadioss/win_scripts_mk3/openradioss_run_script_ps.bat"

    def __init__(self, n_variables: int = 5, instance: int = 1,is_minimisation:bool = True,
                 opt_dimensionality:int=5, lamda:float = LAMDA, reference_intrusion:float = INTRUSION_PRIME):
        
        bounds = ioh.iohcpp.RealBounds(n_variables, -5, 5)
        optimum = ioh.iohcpp.RealSolution([0]* n_variables, 0.0)
        super().__init__('Star_Box_Problem_Sim', n_variables, instance, is_minimisation, bounds, [],  optimum)
  
        self.__cur_intrusion:float = -np.inf
        self.__cur_sea:float = -np.inf

        self.__opt_dimensionality:int = opt_dimensionality

        self.__lamda:float = lamda

        self.__reference_intrusion:float = reference_intrusion

        # Generate a problem instance
        self.__prob_inst:sob.StarBox = sob.get_problem(model_type=1,
                                                       dimension=n_variables,
                                                       output_data=["specific_energy",'intrusion'],
                                                       batch_file_path=Starbox_problem.batch_file_path,
                                                       sequential_id_numbering=False)
 
    def evaluate(self, X:np.ndarray):
        
        # Reshape the array to be read by the sob objects
        X_ = X.tolist()
        
        responses:list = self.__prob_inst(X_,problem_id=self.state.evaluations+1)
        self.__cur_sea = responses[0]
        self.__cur_intrusion = responses[1]

        funct_eval:float = self.__cur_sea - self.__lamda*np.abs(self.__cur_intrusion-self.__reference_intrusion)

        if self.meta_data.optimization_type == "MAX":
            return funct_eval
        else:
            return -1* funct_eval
        
    
    @property
    def cur_intrusion(self)->float:
        return self.__cur_intrusion
    
    @property
    def cur_sea(self)->float:
        return self.__cur_sea
    
    @property
    def opt_dimensionality(self)->int:
        return self.__opt_dimensionality
    
    @opt_dimensionality.setter
    def opt_dimensionality(self,new_opt_dimensionality:int)-> None:
        self.__opt_dimensionality = new_opt_dimensionality
    
    @property
    def lamda(self)->float:
        return self.__lamda
    
    @lamda.setter
    def lamda(self,new_lamda:float)->None:
        self.__lamda = new_lamda
    
    @property
    def reference_instrusion(self)->float:
        return self.__reference_intrusion
    
    @reference_instrusion.setter
    def reference_instrusion(self,new_reference_instrusion:float)->None:
        self.__reference_intrusion = new_reference_instrusion

    @property
    def prob_inst(self)->sob.StarBox:
        return self.__prob_inst
    
    @property
    def blocking_algorithm(self)->bool:
        return self.__blocking_algorithm
    
    @blocking_algorithm.setter
    def blocking_algorithm(self,new_blocking_set:bool)->None:
        self.__blocking_algorithm = new_blocking_set


def modify_array(x_0:np.ndarray):
    # This is the decorator factory function
    def decorator(func):
        def wrapper(x_inp:np.ndarray):
            # Apply modifications based on `modification_type`
            arr = np.array(x_inp).ravel()
            x_01 = np.array(x_0).ravel()
            if x_inp.size == 1:
                modified_array = np.array([x_01[0],x_01[1],x_01[2],x_01[3],arr[0]])
            elif x_inp.size == 3:
                modified_array = np.array([arr[0],arr[1],x_01[2],x_01[3],arr[2]])
            else:
                modified_array = np.array(x_inp).ravel()  # No modification if type is unknown

            # Call the original function with the modified array
            return func(modified_array)
        return wrapper
    return decorator


def return_simulation_setup(initial_val:np.ndarray,changing_dimensions:int=1,
                            initial_budget:int=BUDGET,full_sampler:bool=True)->cma.CMAOptions:
    
    if not initial_val.size  == 5:
        raise ValueError("The initial value should be of size 5")

    # Fill the fixed variables property
    if changing_dimensions == 1:
        fixed_vars = {
            0:initial_val[0,0],
            1:initial_val[0,1],
            2:initial_val[0,2],
            3:initial_val[0,3],
        }
    elif changing_dimensions ==3:
        fixed_vars = {
            2:initial_val[0,2],
            3:initial_val[0,3],
        }
    else:
        fixed_vars = {}

    # Initialize the CMA-ES Object
    # Options

    if full_sampler:
        samp:cma.sampler = cma.sampler.GaussFullSampler
    else:
        samp:cma.sampler = cma.sampler.GaussStandardConstant

    opts:cma.CMAOptions = cma.CMAOptions()
    opts.set({'bounds':[-5.0,5.0],
            'tolfun': 1e-12,
            'maxfevals':initial_budget,
            'CMA_sampler':samp,
            'fixed_variables':fixed_vars,
    })

    # Return the options of the CMA-ES algorithm

    return opts

def adjust_initial_input(x0:np.ndarray,dim:int):

    x_01:np.ndarray = x0.ravel()

    if not x0.size==5:
        raise ValueError("The size of the array must be 5")
    if dim ==1:
        return np.array([x_01[-1]]).ravel()
    elif dim ==3:
        return np.array([x_01[0], x_01[1],x_01[-1]]).ravel()
    elif dim ==5:
        return x_01
    else:
        raise ValueError("Dimension is badly set")


def return_initial_mu_lamda(dim:int, case:int)->List[int]:

    if dim not in (1,3,5):
        raise ValueError("The dimension should be an integer equal to 1,3 or 5")
    
    # Some lambda functions to compute the parent and offspring sizes
    default_lamda = lambda d: int(np.floor(4+3*np.log(d)))
    default_mu = lambda lamdda: int(np.ceil(lamdda/2))

    # Now evaluate the cases
    if case == 0:
        # The computation is let free depending on the dimension
        lamda = default_lamda(dim)
        return lamda,default_mu(lamda)
    
    elif case == 1:
        lamda = default_lamda(1)

    elif case == 2:
        lamda = default_lamda(3)

    elif case == 3:
        lamda = default_lamda(5)
    

    return default_mu(lamda), lamda


def main():

    dim = CURDIM

    mu_0, lamda_0 = return_initial_mu_lamda(dim=dim,case=POP_CASE)

    if CMA_ES_HAND ==1:
        extra_str = "Hansen"
    elif CMA_ES_HAND ==2:
        extra_str = "Vermetten-De_Nobel"

    logger_ioh = ioh.logger.Analyzer(
        root=os.getcwd(),                  # Store data in the current working directory
        folder_name=f"Star_Box_CMA_ES_{dim}D_dual",       # in a folder named: 'my-experiment'
        algorithm_name=f"CMA-ES_{extra_str}",    # meta-data for the algorithm used to generate these results
        store_positions=True,               # store x-variables in the logged files
        triggers= DEFAULT_TRIGGERS,
    )

    # Initialize a seed to instantiate the simulator
    x0 = np.random.uniform(-5,5,(1,5))

   
    # Initialize a new IOH instance
    prob:ioh.problem.RealSingleObjective = Starbox_problem(n_variables=5,
                                                           instance=1,
                                                           is_minimisation=True,
                                                           opt_dimensionality=dim)


    # Return the evolutionary strategy setup
    opts:cma.CMAOptions = return_simulation_setup(initial_val=x0,
                                                          changing_dimensions=dim,
                                                          initial_budget=BUDGET)

    # Get the logger to read the properties from the object
    # This has to be called before attaching the logger to the ioh.problem instance
    # since the properties will not be watched by the logger
    logger_ioh.watch(prob,['cur_intrusion','cur_sea'])
    logger_ioh.add_run_attributes(prob,['opt_dimensionality','lamda'])
    
 

    prob.attach_logger(logger_ioh)
    

    if CMA_ES_HAND == 1:
        #Evaluate the test
        # The handle 'eval_initial_x' is to force the algorithm to 
        # evaluate the seed solution
        result:cma.evolution_strategy.CMAEvolutionStrategyResult = cma.fmin(prob,x0,SIGMA_0,
                                                                            opts,restarts=MAX_RESTARTS_DEFAULT,
                                                                            eval_initial_x=True)
    
    elif CMA_ES_HAND ==2:


        @modify_array(x0)
        def obj_func_1(x_, func:Starbox_problem=prob) -> float:
            return func(x_)
        
        
        obj_func_1(x0)


        module_1 = c_maes.parameters.Modules()

        if MAX_RESTARTS_DEFAULT == 0:
            module_1.restart_strategy = c_maes.options.RestartStrategy.NONE
        
        else:
            module_1.restart_strategy = c_maes.options.RestartStrategy.IPOP
            
        #module_1.restart_strategy = c_maes.options.RESTART
        #module_2.restart_strategy = c_maes.options.RESTART

        module_1.bound_correction = c_maes.options.CorrectionMethod.SATURATE
        

        module_1.matrix_adaptation = c_maes.options.MatrixAdaptationType.COVARIANCE

        x01 = adjust_initial_input(x0=x0,dim=dim)

        opts = c_maes.parameters.Settings(dim=dim, modules=module_1,x0=x01, 
                                            sigma0 =SIGMA_0, budget=BUDGET, verbose=True,
                                            mu0 = mu_0, lambda0 = lamda_0,
                                            lb=np.array([-5.0]*dim), 
                                            ub=np.array([5.0]*dim))
        

        parameters = c_maes.Parameters(opts)

        cma1 = c_maes.ModularCMAES(parameters)


        cma1.run(obj_func_1)
        #cma1.run(ff)
        #while not cma1.break_conditions():
        #    cma1.step(ff)


    
    # Reset the IOH problem holder
    prob.reset()

    # Detach the logger from the problem
    prob.detach_logger()


    


if __name__ == '__main__':
    os.chdir("/home/ivanolar/Documents/Test_Output")
    main()