
from openmm import openmm, unit


class mdp_parser:
    def __init__(self):
        self.integrator = openmm.LangevinIntegrator
        self.dt = 2 * unit.femtosecond
        self.nstmaxh = 1000                # time up check interval
        self.nsteps = 5000                 # number of step
        self.ncycle = 0                    # number of cycle
        self.nstxout_compressed = 5000     # save xtc trajectory every X step, 0 means no saving
        self.nstdcd = 0                    # save xtc trajectory every X step, 0 means no saving
        self.nstenergy = 1000              # save csv file every X step
        self.tau_t = 2.0 * unit.picosecond # 1/gama, inverse friction constant
        self.ref_t = 298 * unit.kelvin     # reference temperature
        self.gen_vel = False               #
        self.gen_temp = 298 * unit.kelvin  #
        self.restraint = False             #
        self.res_fc = 1000                 # restraint force constant, in kJ/mol/nm^2
        self.pcoupltype = None             # can be "None", "isotropic", "semiisotropic/membrane", "anisotropic"
        self.ref_p = 1.0 * unit.bar        #
        self.nstpcouple = 25               # in steps
        self.surface_tension = 0.0         # in kJ/mol/nm^2

    def read(self, input_mdp):
        with open(input_mdp) as f:
            lines = f.readlines()
        for line in lines:
            if line.find(';') >= 0: line = line.split(';')[0]
            line = line.strip()
            if "=" in line:
                segments = line.split('=')
                input_param = segments[0].lower().strip().replace("-","_")
                inp_val = segments[1].strip()
                if input_param == "integrator":
                    if   inp_val == "LangevinIntegrator":       self.integrator = openmm.LangevinIntegrator
                    elif inp_val == "LangevinMiddleIntegrator": self.integrator = openmm.LangevinMiddleIntegrator
                    else: raise ValueError(f"{inp_val} is not support in mdp_parser")
                if input_param == "dt":                 self.dt = float(inp_val) * unit.picosecond
                if input_param == "nstmaxh":            self.nstmaxh = int(inp_val)
                if input_param == "nsteps":             self.nsteps = int(inp_val)
                if input_param == "ncycle":             self.ncycle = int(inp_val)
                if input_param == "nstxout_compressed": self.nstxout_compressed = int(inp_val)
                if input_param == "nstdcd":             self.nstdcd = int(inp_val)
                if input_param == "nstenergy":          self.nstenergy = int(inp_val)
                if input_param == "tau_t":              self.tau_t = float(inp_val) * unit.picosecond
                if input_param == "ref_t":              self.ref_t = float(inp_val) * unit.kelvin
                if input_param == "gen_vel":
                    if   inp_val.lower() in ["yes", "on" ]: self.gen_vel = True
                    elif inp_val.lower() in ["no",  "off"]:  self.gen_vel = False
                    else : raise ValueError(f"{inp_val} is not a valid input for gen_vel")
                if input_param == "gen_temp":   self.gen_temp = float(inp_val) * unit.kelvin
                if input_param == "restraint":
                    if   inp_val.lower() in ["yes", "on" ]:  self.restraint = True
                    elif inp_val.lower() in ["no",  "off"]:  self.restraint = False
                    else : raise ValueError(f"{inp_val} is not a valid input for restraint")
                if input_param == "res_fc":     self.res_fc = float(inp_val)
                if input_param == "pcoupltype":
                    if inp_val.lower() in ["isotropic", "membrane", "anisotropic"]:
                        self.pcoupltype = inp_val.lower()
                    elif inp_val.lower()== "semiisotropic":
                        self.pcoupltype = "membrane"
                    elif inp_val.lower() == "none":
                        self.pcoupltype = None
                    else:
                        raise ValueError(f"{inp_val} is not a valid input for pcoupltype")
                if input_param == "ref_p":      self.ref_p = [float(i) for i in inp_val.split()] * unit.bar
                if input_param == "nstpcouple": self.nstpcouple = int(inp_val)
                if input_param == "surface_tension": self.surface_tension = float(inp_val)
        return self


