import unittest
from Omdrun.mdp import mdp_parser
from openmm import unit

class MyTestCase(unittest.TestCase):
    def test_mdp_nvt(self):
        mdp_inputs = mdp_parser().read("01-mdps/nvt_charmm.mdp")
        self.assertEqual(mdp_inputs.pcoupltype, None)

    def test_mdp_npt(self):
        mdp_inputs = mdp_parser().read("01-mdps/npt_charmm.mdp")
        self.assertEqual(mdp_inputs.pcoupltype, "membrane")
        self.assertEqual(mdp_inputs.tau_p, 25)
        self.assertEqual(mdp_inputs.ref_p, [1.0] * unit.bar)


if __name__ == '__main__':
    unittest.main()
