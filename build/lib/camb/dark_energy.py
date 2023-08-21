from .baseconfig import F2003Class, fortran_class, numpy_1d, CAMBError, np, \
    AllocatableArrayDouble, f_pointer
from ctypes import c_int, c_double, byref, POINTER, c_bool
import numpy as np

max_num_of_fluids = 4
max_num_of_params = 5

class DarkEnergyModel(F2003Class):
    """
    Abstract base class for dark energy model implementations.
    """
    _fields_ = [
        ("__is_cosmological_constant", c_bool),
        ("__num_perturb_equations", c_int),
        ("num_of_components", c_int)]

    def validate_params(self):
        return True


class DarkEnergyEqnOfState(DarkEnergyModel):
    """
    Abstract base class for models using w and wa parameterization with use w(a) = w + (1-a)*wa parameterization,
    or call set_w_a_table to set another tabulated w(a). If tabulated w(a) is used, w and wa are set
    to approximate values at z=0.

    See :meth:`.model.CAMBparams.set_initial_power_function` for a convenience constructor function to
    set a general interpolated P(k) model from a python function.

    """
    _fortran_class_module_ = 'DarkEnergyInterface'
    _fortran_class_name_ = 'TDarkEnergyEqnOfState'

    _fields_ = [
        ("w", c_double, "w(0)"),
        ("wa", c_double, "-dw/da(0)"),
        ("cs2", c_double, "fluid rest-frame sound speed squared"),
        ("use_tabulated_w", c_bool, "using an interpolated tabulated w(a) rather than w, wa above"),
        ("__no_perturbations", c_bool, "turn off perturbations (unphysical, so hidden in Python)")
    ]

    _methods_ = [('SetWTable', [numpy_1d, numpy_1d, POINTER(c_int)])]

    def set_params(self, w=-1.0, wa=0, cs2=1.0):
        """
         Set the parameters so that P(a)/rho(a) = w(a) = w + (1-a)*wa

        :param w: w(0)
        :param wa: -dw/da(0)
        :param cs2: fluid rest-frame sound speed squared
        """
        self.w = w
        self.wa = wa
        self.cs2 = cs2
        self.validate_params()

    def validate_params(self):
        if not self.use_tabulated_w and self.wa + self.w > 0:
            raise CAMBError('dark energy model has w + wa > 0, giving w>0 at high redshift')

    def set_w_a_table(self, a, w):
        """
        Set w(a) from numerical values (used as cublic spline). Note this is quite slow.

        :param a: array of scale factors
        :param w: array of w(a)
        :return: self
        """
        if len(a) != len(w):
            raise ValueError('Dark energy w(a) table non-equal sized arrays')
        if not np.isclose(a[-1], 1):
            raise ValueError('Dark energy w(a) arrays must end at a=1')

        self.f_SetWTable(a, w, byref(c_int(len(a))))
        return self


@fortran_class
class DarkEnergyFluid(DarkEnergyEqnOfState):
    """
    Class implementing the w, wa or splined w(a) parameterization using the constant sound-speed single fluid model
    (as for single-field quintessense).

    """

    _fortran_class_module_ = 'DarkEnergyFluid'
    _fortran_class_name_ = 'TDarkEnergyFluid'

    def validate_params(self):
        super().validate_params()
        if not self.use_tabulated_w:
            if self.wa and (self.w < -1 - 1e-6 or 1 + self.w + self.wa < - 1e-6):
                raise CAMBError('fluid dark energy model does not support w crossing -1')


@fortran_class
class DarkEnergyPPF(DarkEnergyEqnOfState):
    """
    Class implementating the w, wa or splined w(a) parameterization in the PPF perturbation approximation
    (`arXiv:0808.3125 <https://arxiv.org/abs/0808.3125>`_)
    Use inherited methods to set parameters or interpolation table.

    """
    # cannot declare c_Gamma_ppf directly here as have not defined all fields in DarkEnergyEqnOfState (TCubicSpline)
    _fortran_class_module_ = 'DarkEnergyPPF'
    _fortran_class_name_ = 'TDarkEnergyPPF'


@fortran_class
class MultiFluidDE(DarkEnergyModel):
    """
    Class for implementing multiple fluid components in dark energy.
    """
    _fortran_class_module_ = 'MultiFluidDE'
    _fortran_class_name_ = 'TMultiFluidDE'

    _fields_ = [
        ("DebugLevel", c_int, "Verbosity of the module"),
        ("models", c_int*max_num_of_fluids, "Which models to use for each fluid."),
        ("w0", c_double, "CPL parameter"),
        ("wa", c_double, "CPL parameter"),
        ("w1", c_double, "Second binned w"),
        ("w2", c_double, "Third binned w"),
        ("w3", c_double, "Fourth binned w"),
        ("z1", c_double, "Bin limit"),
        ("z2", c_double, "Bin limit"),
        ("z3", c_double, "Bin limit"),
        ("zc", c_double, "EDE critical redshift."),
        ("fde_zc", c_double, "EDE energy contribution at critical redshift."),
        ("theta_i", c_double, "Initial field value divided by the potential frequency."),
        ("wn", c_double, "EDE transition EoS parameter."),
        ("n", c_double, "Power law of the axion/monomial potential."),
        ("grhonode_zc", c_double, "Energy density of everything except DE at zc."),
        ("freq", c_double, "Oscillation frequency of the axion field."),
        ("fac1", c_double, "Power law of the axion/monomial potential."),
        ("fac2", c_double, "Energy density of everything except DE at zc."),
        ("fac3", c_double, "Oscillation frequency of the axion field.")
    ]

    def set_params(self, num_of_components, models, w0 = -1, wa = 0,
                   zc = 3000, fde_zc = 0, wn = 1, theta_i = 1,
                   w1 = -1, w2 = -1, w3 = -1, z1 = 0.7, z2 = 1.4, z3 = 2.1):
        """
         Set dark energy fluid parameters.
        """
        # JVR - C array types
        INTARRAYFLUIDS = c_int * max_num_of_fluids

        if len(models)!=max_num_of_fluids:
            raise ValueError("Models array must have length %d" % max_num_of_fluids)
        models_in_ctypes = INTARRAYFLUIDS(models[0], models[1], models[2], models[3])

        self.num_of_components = num_of_components
        self.models = models_in_ctypes
        self.w0 = w0
        self.wa = wa
        self.w1 = w1
        self.w2 = w2
        self.w3 = w3
        self.z1 = z1
        self.z2 = z2
        self.z3 = z3
        self.zc = zc
        self.fde_zc = fde_zc
        self.theta_i = theta_i
        self.wn = wn

@fortran_class
class AxionEffectiveFluid(DarkEnergyModel):
    """
    Example implementation of a specifc (early) dark energy fluid model
    (`arXiv:1806.10608 <https://arxiv.org/abs/1806.10608>`_).
    Not well tested, but should serve to demonstrate how to make your own custom classes.
    """
    _fields_ = [
        ("w_n", c_double, "effective equation of state parameter"),
        ("fde_zc", c_double, "energy density fraction at z=zc"),
        ("zc", c_double, "decay transition redshift (not same as peak of energy density fraction)"),
        ("theta_i", c_double, "initial condition field value")]

    _fortran_class_name_ = 'TAxionEffectiveFluid'
    _fortran_class_module_ = 'DarkEnergyFluid'

    def set_params(self, w_n, fde_zc, zc, theta_i=None):
        self.w_n = w_n
        self.fde_zc = fde_zc
        self.zc = zc
        if theta_i is not None:
            self.theta_i = theta_i


# base class for scalar field quintessence models
class Quintessence(DarkEnergyModel):
    r"""
    Abstract base class for single scalar field quintessence models.

    For each model the field value and derivative are stored and splined at sampled scale factor values.

    To implement a new model, need to define a new derived class in Fortran,
    defining Vofphi and setting up initial conditions and interpolation tables (see TEarlyQuintessence as example).

    """
    _fields_ = [
        ("DebugLevel", c_int),
        ("astart", c_double),
        ("integrate_tol", c_double),
        ("sampled_a", AllocatableArrayDouble),
        ("phi_a", AllocatableArrayDouble),
        ("phidot_a", AllocatableArrayDouble),
        ("__npoints_linear", c_int),
        ("__npoints_log", c_int),
        ("__dloga", c_double),
        ("__da", c_double),
        ("__log_astart", c_double),
        ("__max_a_log", c_double),
        ("__ddphi_a", AllocatableArrayDouble),
        ("__ddphidot_a", AllocatableArrayDouble),
        ("__state", f_pointer)
    ]
    _fortran_class_module_ = 'Quintessence'


@fortran_class
class EarlyQuintessence(Quintessence):
    r"""
    Example early quintessence (axion-like, as arXiv:1908.06995) with potential

     V(\phi) = m^2f^2 (1 - cos(\phi/f))^2 + \Lambda_{cosmological constant}

    """

    _fields_ = [
        ("n", c_double, "power index for potential"),
        ("f", c_double, r"f/Mpl (sqrt(8\piG)f); only used for initial search value if use_zc is True"),
        ("m", c_double, "mass parameter in reduced Planck mass units; "
                        "only used for initial search value of use_zc is True"),
        ("theta_i", c_double, "phi/f initial field value"),
        ("frac_lambda0", c_double, "fraction of dark energy in cosmologicla constant today (approximated as 1)"),
        ("use_zc", c_bool, "solve for f, m to get specific critical reshidt zc and fde_zc"),
        ("zc", c_double, "reshift of peak fractional early dark energy density"),
        ("fde_zc", c_double, "fraction of early dark energy density to total at peak"),
        ("npoints", c_int, "number of points for background integration spacing"),
        ("min_steps_per_osc", c_int, "minimumum number of steps per background oscillation scale"),
        ("fde", AllocatableArrayDouble, "after initialized, the calculated backgroundearly dark energy "
                                        "fractions at sampled_a"),
        ("__ddfde", AllocatableArrayDouble)
    ]
    _fortran_class_name_ = 'TEarlyQuintessence'

    def set_params(self, n, f=0.05, m=5e-54, theta_i=0.0, use_zc=True, zc=None, fde_zc=None):
        self.n = n
        self.f = f
        self.m = m
        self.theta_i = theta_i
        self.use_zc = use_zc
        if use_zc:
            if zc is None or fde_zc is None:
                raise ValueError("must set zc and fde_zc if using use_zc")
            self.zc = zc
            self.fde_zc = fde_zc


# short names for models that support w/wa
F2003Class._class_names.update({'fluid': DarkEnergyFluid, 'ppf': DarkEnergyPPF})
