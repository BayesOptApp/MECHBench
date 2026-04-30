from src.sob import get_problem
import unittest

# =================================================
# CONSTANTS
# =================================================
DEFAULT_RUNNER_OPTIONS = {
    'np': 1,
    'nt': 1
}

class TestMECHBench(unittest.TestCase):

    def test_generate_multiple_instances(self):
        # Generate the first StarBox instance and run simulation
        
        # Get the problem
        a = get_problem(model_type=1,
                        dimension=1,
                        runner_options=DEFAULT_RUNNER_OPTIONS,
                        output_data='mass')
        
        # Generate an input deck for the problem
        a(variable_array=[1],deck_id=1)

        self.assertEqual(a.dimension, 1)
        

        # Generate the second StarBox instance
        b = get_problem(model_type=1,
                        dimension=5,
                        runner_options=DEFAULT_RUNNER_OPTIONS,
                        output_data='mass')
        
        b(variable_array=[1,2,3,4,5],deck_id=2)
        
        self.assertEqual(b.dimension, 5)
    
    def test_check_mass(self):
        a = get_problem(model_type=1,
                        dimension=3,
                        runner_options=DEFAULT_RUNNER_OPTIONS,
                        output_data='mass')
        
        self.assertIsNotNone(a(variable_array=[1,2,3],deck_id=1))

        b = get_problem(model_type=1,
                        dimension=5,
                        runner_options=DEFAULT_RUNNER_OPTIONS,
                        output_data='mass')
        
        self.assertIsNotNone(b(variable_array=[1,2,3,4,5],deck_id=2))

        c = get_problem(model_type=1,
                        dimension=5,
                        runner_options=DEFAULT_RUNNER_OPTIONS,
                        output_data='mass')
        
        self.assertIsNotNone(c(variable_array=[1,2,3,4,2],deck_id=3))

    def test_check_multi_input(self):

        a = get_problem(model_type=1,
                        dimension=3,
                        runner_options=DEFAULT_RUNNER_OPTIONS,
                        output_data='mass')
        
        self.assertIsNotNone(a(variable_array=[1,2,3],deck_id=1))
        self.assertIsNotNone(a(variable_array=[1,2,5],deck_id=2))

    def test_check_absorbed_energy(self):
        a = get_problem(model_type=1,
                        dimension=3,
                        runner_options=DEFAULT_RUNNER_OPTIONS,
                        output_data='absorbed_energy')
        
        self.assertIsNotNone(a(variable_array=[1,2,3],deck_id=1))
        self.assertIsNotNone(a(variable_array=[1,2,5],deck_id=2))

    def test_check_intrusion(self):
        a = get_problem(model_type=1,
                        dimension=3,
                        runner_options=DEFAULT_RUNNER_OPTIONS,
                        output_data='intrusion')
        
        self.assertAlmostEqual(a(variable_array=[1,2,3],deck_id=1), 46.69148, places=5)
        self.assertAlmostEqual(a(variable_array=[1,2,5],deck_id=2), 37.2778, places=5)

    def test_check_crashtube_mesh_generation(self):
        a = get_problem(model_type=3,
                        dimension=2,
                        runner_options=DEFAULT_RUNNER_OPTIONS,
                        output_data='mass')
        
        self.assertIsNotNone(a(variable_array=[1,2],deck_id=1))

        b = get_problem(model_type=3,
                        dimension=3,
                        runner_options=DEFAULT_RUNNER_OPTIONS,
                        output_data='mass')
        
        self.assertIsNotNone(b(variable_array=[1,2,3],deck_id=2))

        c = get_problem(model_type=3,
                        dimension=4,
                        runner_options=DEFAULT_RUNNER_OPTIONS,
                        output_data='mass')
        
        self.assertIsNotNone(c(variable_array=[1,2,3,4],deck_id=3))
        
        d = get_problem(model_type=3,
                        dimension=5,
                        runner_options=DEFAULT_RUNNER_OPTIONS,
                        output_data='mass')
        
        self.assertIsNotNone(d(variable_array=[1,2,3,4,5],deck_id=4))

    def test_check_crashtube_problem(self):
        a = get_problem(model_type=3,
                        dimension=4,
                        runner_options=DEFAULT_RUNNER_OPTIONS,
                        output_data='intrusion')
        
        self.assertIsNotNone(a(variable_array=[1,2,3,4],deck_id=1))
    
    def test_check_threepointbending_problem(self):
        a = get_problem(model_type=2,
                        dimension=4,
                        runner_options=DEFAULT_RUNNER_OPTIONS,
                        output_data='mass')
        
        self.assertIsNotNone(a(variable_array=[1,2,3,4],deck_id=1))

        b = get_problem(model_type=2,
                        dimension=4,
                        runner_options=DEFAULT_RUNNER_OPTIONS,
                        output_data='absorbed_energy')
        
        self.assertIsNotNone(b(variable_array=[-1,2,3,-2],deck_id=2))

        c = get_problem(model_type=2,
                        dimension=4,
                        runner_options=DEFAULT_RUNNER_OPTIONS,
                        output_data='intrusion')
        
        self.assertIsNotNone(c(variable_array=[0,-4,2,0],deck_id=3))

    def test_check_mass_output(self):
        instance_ = get_problem(model_type=3,
                                dimension=6,
                                runner_options=DEFAULT_RUNNER_OPTIONS,
                                output_data='mass')
        
        self.assertIsNotNone(instance_(variable_array=[-1,-3,4,-5,2,3],deck_id=1))
        self.assertIsNotNone(instance_(variable_array=[-1,-3,2,-5,3,1],deck_id=2))
        self.assertIsNotNone(instance_(variable_array=[-1,-3,2,-5,-5,5],deck_id=3))


if __name__ == '__main__':
    unittest.main()

