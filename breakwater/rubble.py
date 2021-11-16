import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from tabulate import tabulate


from breakwater.equipment.equipment import Truck, Vessel, Excavator, Crane, Barge, PlateFeeder
from .core import substructure
from .core.bishop import Bishop
from .core.overtopping import rubble_mound
from .core.scour import scour_protection
from .core.stability import hudson, vandermeer
from .core.toe import toe_stability
from .utils.breakwater_divide import GaussianA, divide_cross_section, maxDiff
from .utils.exceptions import (
    InputError,
    NotSupportedError,
    RockGradingError,
    user_warning,
    EquipmentError,
)

class RubbleMound:
    """General Rubble Mound breakwater class

    Makes a conceptual design for the substructure of a rubble mound
    breakwater. The class computes the necessary nominal diameter and
    rock class of the underlayer, and a filter layer if one is needed.
    Depending on the rock class it is possible that a new variant is
    generated, these are identified by an a, b, c or d. In the attribute
    :py:attr:`variantIDs` a list of generated variants is stored.
    Furthermore, it computes the required crest freeboard of the
    breakwater by using the maximum allowed overtopping discharge from
    the LimitStates with the formulas from EurOtop (2018). Additionally,
    the toe of the breakwater is designed with toe stability formula of
    Van der Meer (1998).

    Parameters
    ----------
    Dn50 : float
        Dn50 of the armour [m]
    core_material : dict
        {class: str, Dn50: float}
    rho : float
        density of the material of the armour layer [kg/m³]
    rho_w : float
        density of water [kg/m³]
    armour_layer : str
        type of armour layer, fully supported armour layers are Rock,
        Xbloc and XblocPlus. Full support means that the rules for the
        underlayer have only been included for these types of armour
        layer. In case another armour layer is used the rule for the
        underlayer, the :py:obj:`filter_rule` must manually be set.
    layers : int
        number of layers in the armour layer
    LimitStates : list
        list of :py:class:`LimitState` ULS, SLS or other defined limit
        states defined with :py:class:`LimitState`. When designing with
        one limit state it must still be entered as a list.
    Grading : :py:class:`RockGrading`
        standard rock grading defined in the NEN-EN 13383-1 or a user
        defined rock grading
    safety : float, optional, default: 1
        safety factor of design (number of standard deviations from the
        mean). Positive values increase the safety, negative values
        decrease the safety of the breakwater
    layers_underlayer : int, optional, default: 2
        number of layers in the underlayer
    slope_toe : tuple, optional, default: (2, 3)
        slope of the toe
    B_toe : float, optional, default: None
        width of the top of the toe in meters. By default the width of
        toe is taken as 3 * Dn50_toe.
    slope : tuple, optional, default: None
        Slope of the armour layer (V, H). For example a slope of 3V:4H
        is defined as (3, 4)
    B : float, optional, default: None
        Crest width [m]
    beta : float, optional, default: 0
        angle between direction of wave approach and a line normal to
        the breakwater (degrees).
    filter_rule : {'Rock', 'Xbloc', 'XblocPlus'}, optional, default: None
        filter rule to use for the substructure of the breakwater, for
        Rock, Xbloc and XblocPlus the correct filter rule is
        automatically selected. In case another type of armour layer is
        used one of these filter rules must be chosen.
    Soil : :py:class:`Soil`, optional, default: None
        by default Soil is None, which means that the geotechnical checks
        are not performed. By specifying a Soil object, the geotechnical
        checks are automatically performed.
    phi : float, optional, default: 40
        internal friction angle of rock [degrees]
    id : int, optional, default: None
        add a unique id to the breakwater

    Attributes
    ----------
    logger : dict
        dict of warnings and messages
    structure : dict
        dictionary with the computed Dn50, rock class, and average Dn50
        of the rock class for each layer and the toe. This dictionary
        includes all variants, use :py:meth:`get_variant` to get the
        parameters of one specific variant. Alternatively,
        :py:meth:`print_variants` can be used to print the details of
        one, multiple or all variants.
    alpha : float
        slope of the structure in radians
    id : int
        unique id of the breakwater
    variantIDs : list
        list with the IDs of the variants generated for this rubble
        mound breakwater.
    Rc : float
        the crest freeboard of the structure [m]
    width_scour : float
        the required length of the scour protection [m]
    bishop : :py:class:`Bishop`
        the :py:class:`Bishop` that was normative in the computation
    F_norm : float
        normative factor of safety computed with :py:class:`Bishop`
    """

    def __init__(
        self,
        Dn50,
        core_material,
        rho,
        rho_w,
        armour_layer,
        layers,
        LimitStates,
        Grading,
        safety=1,
        layers_underlayer=2,
        slope_toe=(2, 3),
        structure_type=None,
        B_toe=None,
        slope=None,
        B=None,
        beta=0,
        filter_rule=None,
        Soil=None,
        phi=40,
        id=None,
        **kwargs,
    ):
        """See help(RubbleMound) for more info"""
        # if not from RockRubbleMound or ConcreteRubbleMound
        if (
            not isinstance(self, RockRubbleMound)
            and not isinstance(self, ConcreteRubbleMound)
            and not isinstance(self, ConcreteRubbleMoundRevetment)
        ):
            # not called from a child class
            # therefore some attributes must be set
            self.logger = {"INFO": [], "WARNING": []}
            self.structure = {
                "armour": {
                    "computed Dn50": Dn50,
                    "class": None,
                    "class Dn50": Dn50,
                    "state": None,
                }
            }
            self.alpha = np.arctan(slope[0] / slope[1])
            self.id = id
            self.variantIDs = ["a", "b", "c", "d"]

        # set attribute of bishop and normative F
        self.bishop = None
        self.F_norm = None

        Dn50_core = core_material['Dn50']
        self.Grading = Grading
        self.depth_area = None

        # set input as private attribute
        self._input_arguments = {
            "structure_type": structure_type,
            "slope": slope,
            "slope_toe": slope_toe,
            "B": B,
            "armour": armour_layer,
            "Grading": Grading,
            "core_class": core_material['class'],
            "Dn50_core": core_material['Dn50'],
        }

        # check for supported armour layers, and if filter_rule is set
        supported = substructure._supported_armour_layers()
        if armour_layer in supported:
            filter_rule = armour_layer
        elif filter_rule is None:
            supported_rules = ", ".join(supported)
            raise NotSupportedError(
                (
                    f"Filter rule for {armour_layer} is not implemented, set "
                    f"filter rule to use with filter_rule to {supported_rules}"
                )
            )

        # set the LimitStates as private attribute to use in plot
        self._LimitStates = LimitStates

        # design the first underlayer of the breakwater
        # rho is the density of the material of the armour layer
        computed_dn_u = substructure.underlayer(
            Dn_armour=Dn50, armour_layer=filter_rule, rho=rho, rho_rock=Grading.rho
        )

        # set empty lists to store design values
        class_underlayer, class_dn_u = [], []
        class_filter, class_dn_f, computed_dn_f = [], [], []

        for i, dn in enumerate(computed_dn_u):
            rock_class = Grading.get_class(dn)

            if rock_class in class_underlayer:
                # if rock class is already in the underlayer there is
                # no need to generate an additional variant
                continue
            else:
                # new design and values must thus be saved
                new_underlayer = True

                class_dn = Grading.get_class_dn50(rock_class)
                class_underlayer.append(rock_class)
                class_dn_u.append(class_dn)

                # design a second underlayer/filter layer
                dn_filter_range = substructure.filter_layers(Dn=class_dn, rho=rho)
                for dn_filter in dn_filter_range:
                    if dn_filter > Dn50_core:
                        # computed dn of the filter is larger than the
                        # core, a filter layer is thus needed
                        rock_class = Grading.get_class(dn_filter)
                        if rock_class in class_filter and not new_underlayer:
                            # if rock class is already in the filter layer
                            # there is no need to generate an additional
                            # variant
                            continue
                        else:
                            # new design and values must thus be saved
                            class_dn = Grading.get_class_dn50(rock_class)
                            computed_dn_f.append(dn_filter)
                            class_filter.append(rock_class)
                            class_dn_f.append(class_dn)
                            new_underlayer = False
                    else:
                        # a filter layer is not needed
                        if new_underlayer:
                            # set all values for the filter layer to None
                            computed_dn_f.append(None)
                            class_filter.append(None)
                            class_dn_f.append(None)
                            # set new_underlayer to False so that this
                            # action is not repeated
                            new_underlayer = False
                        else:
                            # still designing for the same underlayer
                            # so no need to set new values
                            continue

        # add underlayer and filter to the structure
        self.structure["underlayer"] = {
            "computed Dn50": computed_dn_u,
            "class": class_underlayer,
            "class Dn50": class_dn_u,
            "state": "see armour",
            "layers": layers_underlayer,
        }
        if any(class_filter):
            self.structure["filter layer"] = {
                "computed Dn50": computed_dn_f,
                "class": class_filter,
                "class Dn50": class_dn_f,
                "state": "see armour",
                "layers": layers_underlayer,
            }

        # determine number of variants
        self.variantIDs = self.variantIDs[: len(class_filter)]

        # add generated variants to the logger
        if len(class_underlayer) == 2:
            self.logger["INFO"].append(
                "two rock classes possible for the underlayer, generated new "
                "variant b"
            )
            if len(class_filter) == 4:
                self.logger["INFO"].append(
                    "two rock classes possible for the filter layer, "
                    "generated new variant c"
                )
                self.logger["INFO"].append(
                    "two rock classes possible for the filter layer, "
                    "generated new variant d"
                )
        if len(class_filter) == 3:
            self.logger["INFO"].append(
                "two rock classes possible for the filter layer, generated "
                "new variant c"
            )

        # design toe of the structure
        delta_rock = (Grading.rho - rho_w) / rho_w

        # make a first estimate for the height of the toe
        h_toe = self._estimate_htoe()

        # set variables to store results from normative LimitState
        Dn50_toe, state_toe = 0, 0

        for i, LimitState in enumerate(LimitStates):
            # get values from the LimitState
            Hs = LimitState.get_Hs(definition="H13")
            h = LimitState.h
            Nod = LimitState["Nod"]

            # make first estimate for the water level above the toe
            ht_estimate = h - h_toe

            # use while loop since ht depends on dn50 of the toe
            # first set temporary values for the while loop
            Dn50_toe_temp = 0
            compute_toe = True
            counter = 0

            while compute_toe:
                Dn50_toe_computed = toe_stability(
                    Hs=Hs, h=LimitState.h, ht=ht_estimate, Delta=delta_rock, Nod=Nod
                )

                # check for convergence
                if abs(Dn50_toe_computed - Dn50_toe_temp) < 0.05 or counter > 50:
                    # value has converged, so break loop
                    compute_toe = False

                # replace old value with the new one
                Dn50_toe_temp = Dn50_toe_computed

                # make new estimate for the water level above the toe
                ht_estimate = h - self._estimate_htoe(Dn50=Dn50_toe_computed)

                counter += 1
            # check if computed Dn50 of current LimitState is larger
            # than current normative Dn50
            if Dn50_toe_temp > Dn50_toe:
                # if larger the normative Dn50 must be changed
                Dn50_toe = Dn50_toe_temp
                state_toe = i

        class_toe = Grading.get_class(Dn50_toe)
        class_Dn50 = Grading.get_class_dn50(class_toe)

        self.structure["toe"] = {
            "computed Dn50": Dn50_toe,
            "class": class_toe,
            "class Dn50": class_Dn50,
            "state": state_toe,
        }

        # determine crest height
        self.Rc, self._state_overtopping = 0, 0

        # check if crest width is more than 3*Dn50 of armour
        if B >= 3 * Dn50:
            # no action required
            pass
        else:
            # get armour_layer
            if armour_layer == "Rock":
                material = "armourstones"
            else:
                material = "units"

            user_warning(
                (
                    f"Given crest width is smaller than three {material}, it is "
                    "advised to increase the width to at least "
                    f"{np.round(3*Dn50, 2)} m"
                )
            )

        # convert beta from deg to rad for overtopping computation
        beta = beta * np.pi / 180.0

        for i, LimitState in enumerate(LimitStates):
            Hm0 = LimitState.get_Hs(definition="Hm0")
            xi = LimitState.surf_similarity(alpha=self.alpha, number="spectral")

            Rc_temp = rubble_mound(
                Hm0=Hm0,
                q=LimitState["q"],
                xi_m_min_1=xi,
                alpha=self.alpha,
                beta=beta,
                gamma_b=1,
                gamma_v=1,
                gam_star=1,
                Gc=B,
                Dn50=Dn50,
                armour_layer=armour_layer,
                layers=layers,
                permeability="permeable",
                safety=safety,
            )

            # check if computed Rc of current LimitState is larger
            # than current normative Rc
            if Rc_temp > self.Rc:
                # if larger the normative Rc must be changed
                self.Rc = Rc_temp
                self._state_overtopping = i

        # now that the toe is designed the width of the toe can be computed
        # check if a width has been given
        if B_toe is None:
            # set toe width to 3 times Dn50
            self.B_toe = 3 * self.structure["toe"]["class Dn50"]
        else:
            # user specified width
            self.B_toe = B_toe

            # check if given width is larger than 3 stones
            if B_toe < 3 * self.structure["toe"]["class Dn50"]:
                user_warning(
                    (
                        "given width of the toe is smaller than three times the "
                        "Dn50 of the toe"
                    )
                )

        # Compute required scour protection
        self.width_scour = 0
        for LimitState in LimitStates:
            w = scour_protection(L=LimitState.L(period="Tm"), slope=slope)

            # check if larger than previous value
            if w >= self.width_scour:
                # set w as new width scour
                self.width_scour = w

        # check if a soil has been given
        if Soil is not None:
            # show warning as implementation is not yet verified
            user_warning(
                (
                    f"The implementation of Bishop into {type(self).__name__} "
                    "has not yet been verified with an example"
                )
            )

            # check if called from ConcreteRubbleMound
            if isinstance(self, ConcreteRubbleMound):
                # raise warning that armour is modelled as rock
                user_warning("The armour layer is modelled as Rock in Bishop")

            # compute the height of the breakwater
            h = LimitStates[self._state_overtopping].h + self.Rc

            # compute horizontal length of the slope
            x = slope[1] * h / slope[0]

            # determine number of slices, with slice width of 1 m
            num_slices = np.round(h, 0)

            n = 0.4

            # compute volumetric weights
            gamma_w = rho_w * 9.81 / 1000
            gamma_r = (1 - n) * Grading.rho * 9.81 / 1000
            gamma_r_sat = gamma_r + n * gamma_w

            # set variable to store normative factor of safety
            self.F_norm = 10 * 10

            # iterate over the LimitStates
            for LimitState in LimitStates:
                # create bishop object
                slip = Bishop(point2=(x, h), wlev=LimitState.h)

                # add soil and rock layer
                slip.add_layer(
                    gamma=Soil.gamma,
                    gamma_sat=Soil.gamma_sat,
                    c=Soil.c,
                    phi=Soil.phi * 180 / np.pi,
                    name="Subsoil",
                    ymin=-50,
                    ymax=0,
                )
                slip.add_layer(
                    gamma=gamma_r,
                    gamma_sat=gamma_r_sat,
                    c=0,
                    phi=phi,
                    name="Rock",
                    ymin=0,
                    ymax=h,
                )

                # compute factor of safety
                slip.compute(num_slices=int(num_slices), gamma_w=gamma_w)

                # check if normative factor of safety
                if slip.circles[slip.normative].F < self.F_norm:
                    # update normative F and bishop object as attribute
                    self.F_norm = slip.circles[slip.normative].F
                    self.bishop = slip

        self.rho = rho

    def _estimate_htoe(self, Dn50=0):
        """Method to estimate the height of the toe"""
        # set ht variable
        ht = 0

        # get normative layer thickness
        for i, id in enumerate(self.variantIDs):
            # get structure
            structure = self.get_variant(id)

            # get thickness of the layer
            t_armour = self._layer_thickness(
                "armour", self._input_arguments["armour"], structure
            )
            t_underlayer = self._layer_thickness("underlayer", "Rock", structure)
            t_filter = self._layer_thickness("filter layer", "Rock", structure)

            ht_est = t_underlayer + t_filter

            if Dn50 != 0:
                ht_est += np.ceil(t_armour / Dn50) * Dn50

            # check if larger than previous estimate
            if ht_est > ht:
                ht = ht_est

        return ht

    def _validate_variant(self, variants):
        """Validate the input of the variant

        Parameters
        ----------
        variants : tuple
            variantIDs given as args input
        """
        # check if a variant is specified
        if not variants:
            # no input is given so get the valid args for variant
            valid_args = self.variantIDs
            valid_args.append("all")

            # raise error
            valid = ", ".join(valid_args)
            raise InputError(
                "did not specify which variants to plot, possible arguments "
                f"are {valid}"
            )

        # check if input all is in variants
        if "all" in variants:
            # set specified variants to all variantIDs
            variants = tuple(self.variantIDs)

        # return the variants
        return variants

    def _layer_thickness(self, layer, material, structure):
        """Compute the thickness of the layer

        Parameters
        ----------
        layer : {armour, underlayer, filter layer}
            name of the layer for which the layer thickness must be
            computed
        material : str
            material of which the layer is made
        structure : dict
            parameters and values (Dn50, Rock class) for one variant

        Returns
        -------
        float
            thickness of the layer [m]
        """
        # check if layer is in structure
        if layer in structure:
            # get number of layers and Dn50 from the structure
            layers = structure[layer]["layers"]
            Dn50 = structure[layer]["class Dn50"]

            # get the layer coefficient
            kt = substructure.layer_coefficient(
                material, layers=layers, placement="standard"
            )

            # compute and return the thickness of the layer
            return kt * layers * Dn50

        else:
            # set layer thickness to 0, as layer is not in the structure
            return 0

    def _layers(self, variantID):
        """compute the coordinates of all layers

        Parameters
        ----------
        variantID : str
            identifier of the variant, see :py:attr:`variantIDs` for a
            list of all generated variants.

        Returns
        -------
        dict
            coordinates of all layers
        """
        # get structure of the current variant
        structure = self.get_variant(variantID)

        # set empty dict to store coordinates in
        coordinates = {}

        # get the slope
        V, H = self._input_arguments["slope"]
        V_toe, H_toe = self._input_arguments["slope_toe"]

        # compute constant to transformthickness of layer to x and
        # y coordinates, switched V and H because orthogonality
        transform_x = V / np.sqrt(V ** 2 + H ** 2)
        transform_y = H / np.sqrt(V ** 2 + H ** 2)

        # get the height and width of the structure
        height = self._LimitStates[self._state_overtopping].h + self.Rc
        B = self._input_arguments["B"]

        # determine thickness of the layers
        t_armour = self._layer_thickness(
            "armour", self._input_arguments["armour"], structure
        )
        t_underlayer = self._layer_thickness("underlayer", "Rock", structure)
        t_filter = self._layer_thickness("filter layer", "Rock", structure)

        if t_filter == 0:
            # add scour protection below underlayer on sea side
            t_scour = 2 * self._input_arguments["Dn50_core"]
        else:
            # filter layer will be used as scour protection
            t_scour = 0

        # check structure type: breakwater or revetment.
        # breakwater is covered with armour on both sides
        if self._input_arguments["structure_type"] == "breakwater":
            # compute armour layer
            armour_y1 = t_filter + t_underlayer + t_scour
            armour_y2 = height
            armour_y3 = height

            armour_x1 = -0.5 * B - H * (armour_y3 - armour_y1) / V
            armour_x2 = -0.5 * B
            armour_x3 = 0.5 * B
            armour_x4 = 0.5 * B + H * armour_y3 / V

            # compute line between armour and underlayer
            arm_under_y1 = t_filter + t_underlayer + t_scour
            arm_under_y2 = t_filter + t_underlayer + t_scour
            arm_under_y3 = height - t_armour
            arm_under_y4 = height - t_armour

            arm_under_x5 = (
                armour_x4 - H * (t_armour * transform_y) / V - t_armour * transform_x
            )
            arm_under_x4 = arm_under_x5 - H * arm_under_y4 / V
            arm_under_x3 = -arm_under_x4
            arm_under_x2 = -H * (arm_under_y3 - arm_under_y2) / V + arm_under_x3

            # compute lower line of the underlayer
            under_y1 = t_filter + t_scour
            under_y2 = t_filter + t_scour
            under_y3 = height - t_armour - t_underlayer
            under_y4 = height - t_armour - t_underlayer

            under_x5 = (
                arm_under_x5
                - H * (t_underlayer * transform_y) / V
                - t_underlayer * transform_x
            )
            under_x4 = under_x5 - H * under_y4 / V
            under_x3 = -under_x4
            under_x2 = -H * (under_y3 - under_y2) / V + under_x3

            # compute the toe
            Dn50 = structure["toe"]["class Dn50"]

            # check if armour layer is made out of rock
            # as there is a different toe structure for rock and armour units
            if self._input_arguments["armour"] == "Rock":
                # compute point where armour layer intersect with the toe
                x = (abs(arm_under_x2 - armour_x1) * V_toe / H_toe) / (
                    V / H + V_toe / H_toe
                )
                y = abs(arm_under_x2 - armour_x1) * V * V_toe / (H * V_toe + V * H_toe)

                # determine height of the toe
                htoe = np.ceil(y / Dn50) * Dn50

                toe_low = t_underlayer + t_filter + t_scour
                toe_top = toe_low + htoe

                toe_x4 = arm_under_x2
                toe_x3 = toe_x4 - H_toe * htoe / V_toe
                toe_x2 = toe_x3 - self.B_toe
                toe_x1 = toe_x2 - H_toe * htoe / V_toe

                # change first points of the armour layer
                # so that the armour layer starts at the intersection
                armour_x1 = armour_x1 + x
                armour_y1 = armour_y1 + y

            else:
                # armour units
                htoe = np.ceil(t_armour / Dn50) * Dn50

                toe_low = t_underlayer + t_filter + t_scour
                toe_top = toe_low + htoe

                toe_x4 = armour_x1
                toe_x3 = armour_x1 + H * htoe / V
                toe_x2 = toe_x3 - self.B_toe
                toe_x1 = toe_x2 - H_toe * htoe / V_toe

            # check if there is a filter layer
            if t_filter != 0:
                # set points of end of the layers
                arm_under_x1 = toe_x1 - 3 * structure["underlayer"]["class Dn50"]
                arm_under_x0 = arm_under_x1 - H * t_underlayer / V
                arm_under_y0 = t_filter

                # compute the lower line of the filter layer
                filter_y3 = height - t_armour - t_underlayer - t_filter
                filter_y4 = height - t_armour - t_underlayer - t_filter

                filter_x5 = (
                    under_x5 - H * (t_filter * transform_y) / V - t_filter * transform_x
                )
                filter_x4 = filter_x5 - H * filter_y4 / V
                filter_x3 = -filter_x4
                filter_x2 = -H * filter_y3 / V + filter_x3
                filter_x1 = arm_under_x0 - self.width_scour
                filter_x0 = filter_x1 - t_filter * H / V

            else:
                # no filter layer
                arm_under_x1 = toe_x1 - 3 * structure["underlayer"]["class Dn50"]
                arm_under_x0 = arm_under_x1 - H * t_underlayer / V
                arm_under_y0 = t_scour

                # add scour protection
                scour_x2 = arm_under_x0 - self.width_scour
                scour_x1 = scour_x2 - t_scour * H / V

            # add lines to the coordinates
            # are added from left to right and back
            coordinates["armour"] = {
                "x": [
                    armour_x1,
                    armour_x2,
                    armour_x3,
                    armour_x4,
                    arm_under_x5,
                    arm_under_x4,
                    arm_under_x3,
                    arm_under_x2,
                    armour_x1,
                ],
                "y": [
                    armour_y1,
                    armour_y2,
                    armour_y3,
                    0,
                    0,
                    arm_under_y4,
                    arm_under_y3,
                    arm_under_y2,
                    armour_y1,
                ],
            }

            coordinates["underlayer"] = {
                "x": [
                    arm_under_x0,
                    arm_under_x1,
                    arm_under_x2,
                    arm_under_x3,
                    arm_under_x4,
                    arm_under_x5,
                    under_x5,
                    under_x4,
                    under_x3,
                    under_x2,
                    arm_under_x0,
                ],
                "y": [
                    arm_under_y0,
                    arm_under_y1,
                    arm_under_y2,
                    arm_under_y3,
                    arm_under_y4,
                    0,
                    0,
                    under_y4,
                    under_y3,
                    under_y2,
                    arm_under_y0,
                ],
            }

            # check if there is a filter to add to coordinates
            if t_filter != 0:
                # add filter to coordinates

                coordinates["filter layer"] = {
                    "x": [
                        filter_x0,
                        filter_x1,
                        under_x2,
                        under_x3,
                        under_x4,
                        under_x5,
                        filter_x5,
                        filter_x4,
                        filter_x3,
                        filter_x2,
                        filter_x0,
                    ],
                    "y": [
                        0,
                        t_filter,
                        under_y2,
                        under_y3,
                        under_y4,
                        0,
                        0,
                        filter_y4,
                        filter_y3,
                        0,
                        0,
                    ],
                }

                # add core to coordinates
                coordinates["core"] = {
                    "x": [filter_x2, filter_x3, filter_x4, filter_x5, filter_x2],
                    "y": [0, filter_y3, filter_y4, 0, 0],
                }
            else:
                coordinates["core"] = {
                    "x": [
                        scour_x1,
                        scour_x2,
                        under_x2,
                        under_x3,
                        under_x4,
                        under_x5,
                        scour_x1,
                    ],
                    "y": [0, t_scour, t_scour, under_y3, under_y4, 0, 0],
                }

            # add the coordinates of the toe
            coordinates["toe"] = {
                "x": [toe_x1, toe_x2, toe_x3, toe_x4, toe_x1],
                "y": [toe_low, toe_top, toe_top, toe_low, toe_low],
            }

        # check structure type: breakwater or revetment. default is breakwater
        # revetment is covered with armour only on sea side
        if self._input_arguments["structure_type"] == "revetment":
            # compute armour layer
            armour_y1 = t_filter + t_underlayer + t_scour
            armour_y2 = height
            armour_y3 = height

            armour_x1 = -0.5 * B - H * (armour_y3 - armour_y1) / V
            armour_x2 = -0.5 * B
            armour_x3 = 0.5 * B
            armour_x4 = 0.5 * B + H * armour_y3 / V

            # compute line between armour and underlayer
            arm_under_y1 = t_filter + t_underlayer + t_scour
            arm_under_y2 = t_filter + t_underlayer + t_scour
            arm_under_y3 = height - t_armour
            arm_under_y4 = height - t_armour

            arm_under_x5 = (
                armour_x4 - H * (t_armour * transform_y) / V - t_armour * transform_x
            )
            arm_under_x4 = arm_under_x5 - H * arm_under_y4 / V
            arm_under_x3 = -arm_under_x4
            arm_under_x2 = -H * (arm_under_y3 - arm_under_y2) / V + arm_under_x3

            # compute lower line of the underlayer
            under_y1 = t_filter + t_scour
            under_y2 = t_filter + t_scour
            under_y3 = height - t_armour - t_underlayer
            under_y4 = height - t_armour - t_underlayer

            under_x5 = (
                arm_under_x5
                - H * (t_underlayer * transform_y) / V
                - t_underlayer * transform_x
            )
            under_x4 = under_x5 - H * under_y4 / V
            under_x3 = -under_x4
            under_x2 = -H * (under_y3 - under_y2) / V + under_x3

            # compute the toe
            Dn50 = structure["toe"]["class Dn50"]

            # check if armour layer is made out of rock
            # as there is a different toe structure for rock and armour units
            if self._input_arguments["armour"] == "Rock":
                # compute point where armour layer intersect with the toe
                x = (abs(arm_under_x2 - armour_x1) * V_toe / H_toe) / (
                    V / H + V_toe / H_toe
                )
                y = abs(arm_under_x2 - armour_x1) * V * V_toe / (H * V_toe + V * H_toe)

                # determine height of the toe
                htoe = np.ceil(y / Dn50) * Dn50

                toe_low = t_underlayer + t_filter + t_scour
                toe_top = toe_low + htoe

                toe_x4 = arm_under_x2
                toe_x3 = toe_x4 - H_toe * htoe / V_toe
                toe_x2 = toe_x3 - self.B_toe
                toe_x1 = toe_x2 - H_toe * htoe / V_toe

                # change first points of the armour layer
                # so that the armour layer starts at the intersection
                armour_x1 = armour_x1 + x
                armour_y1 = armour_y1 + y

            else:
                # armour units
                htoe = np.ceil(t_armour / Dn50) * Dn50

                toe_low = t_underlayer + t_filter + t_scour
                toe_top = toe_low + htoe

                toe_x4 = armour_x1
                toe_x3 = armour_x1 + H * htoe / V
                toe_x2 = toe_x3 - self.B_toe
                toe_x1 = toe_x2 - H_toe * htoe / V_toe

            # check if there is a filter layer
            if t_filter != 0:
                # set points of end of the layers
                arm_under_x1 = toe_x1 - 3 * structure["underlayer"]["class Dn50"]
                arm_under_x0 = arm_under_x1 - H * t_underlayer / V
                arm_under_y0 = t_filter

                # compute the lower line of the filter layer
                filter_y3 = height - t_armour - t_underlayer - t_filter
                filter_y4 = height - t_armour - t_underlayer - t_filter

                filter_x5 = (
                    under_x5 - H * (t_filter * transform_y) / V - t_filter * transform_x
                )
                filter_x4 = filter_x5 - H * filter_y4 / V
                filter_x3 = -filter_x4
                filter_x2 = -H * filter_y3 / V + filter_x3
                filter_x1 = arm_under_x0 - self.width_scour
                filter_x0 = filter_x1 - t_filter * H / V

            else:
                # no filter layer
                arm_under_x1 = toe_x1 - 3 * structure["underlayer"]["class Dn50"]
                arm_under_x0 = arm_under_x1 - H * t_underlayer / V
                arm_under_y0 = t_scour

                # add scour protection
                scour_x2 = arm_under_x0 - self.width_scour
                scour_x1 = scour_x2 - t_scour * H / V

            # add lines to the coordinates
            # are added from left to right and back
            coordinates["armour"] = {
                "x": [
                    armour_x1,
                    armour_x2,
                    under_x4,
                    under_x4,
                    arm_under_x3,
                    arm_under_x2,
                    armour_x1,
                ],
                "y": [
                    armour_y1,
                    armour_y2,
                    armour_y3,
                    arm_under_y3,
                    arm_under_y3,
                    arm_under_y2,
                    armour_y1,
                ],
            }

            coordinates["underlayer"] = {
                "x": [
                    arm_under_x0,
                    arm_under_x1,
                    arm_under_x2,
                    arm_under_x3,
                    under_x4,
                    under_x4,
                    under_x3,
                    under_x2,
                    arm_under_x0,
                ],
                "y": [
                    arm_under_y0,
                    arm_under_y1,
                    arm_under_y2,
                    arm_under_y3,
                    arm_under_y4,
                    under_y4,
                    under_y3,
                    under_y2,
                    arm_under_y0,
                ],
            }

            # check if there is a filter to add to coordinates
            if t_filter != 0:
                # add filter to coordinates

                coordinates["filter layer"] = {
                    "x": [
                        filter_x0,
                        filter_x1,
                        under_x2,
                        under_x3,
                        under_x4,
                        under_x5,
                        filter_x5,
                        filter_x4,
                        filter_x3,
                        filter_x2,
                        filter_x0,
                    ],
                    "y": [
                        0,
                        t_filter,
                        under_y2,
                        under_y3,
                        under_y4,
                        0,
                        0,
                        filter_y4,
                        filter_y3,
                        0,
                        0,
                    ],
                }

                # add core to coordinates
                coordinates["core"] = {
                    "x": [filter_x2, filter_x3, filter_x4, filter_x5, filter_x2],
                    "y": [0, filter_y3, filter_y4, 0, 0],
                }
            else:
                coordinates["core"] = {
                    "x": [
                        scour_x1,
                        scour_x2,
                        under_x2,
                        under_x3,
                        under_x4,
                        under_x5,
                        scour_x1,
                    ],
                    "y": [0, t_scour, t_scour, under_y3, under_y4, 0, 0],
                }

            # add the coordinates of the toe
            coordinates["toe"] = {
                "x": [toe_x1, toe_x2, toe_x3, toe_x4, toe_x1],
                "y": [toe_low, toe_top, toe_top, toe_low, toe_low],
            }
        # return the coordinates of the specified variants
        return coordinates


    def install_all_equipments(
        self,
        coords,
        equip,
        layer,
        grading,
        mass,
        end_lay,
        area,
        key,
        max_height,
        max_height_xcorner,
        length_top,
        cost_key
    ):
        """
        Function which updates the depth_area dict based on if an equipment can install a section
        Parameters
        ----------
        coords: list
            coordinates of the section
        equip: object
            An equipment object
        layer: str
            Which layer are we installing (e.g. core)
        grading: str
            What is the grading of the layer
        mass: float
            The mass of the material of the layer
        end_lay: float
            Upper y-coordinate of the section
        area: float
            Area of the section
        key: str
            the range of the section (e.g. '0.0-1.0')
        max_height: float
            Maximum height of already build layer
        max_height_xcorner: float
            Corner of heightest build layer
        length_top: float
            Length of highest build layer at max_height
        cost_key: str
            Either CO2 or cost

        Returns
        -------
        tuple
        """
        x, y = list(zip(*coords))
        x, y = np.array(x), np.array(y)
        slope = self._input_arguments["slope"]
        xcorner_layer = max(coords)[0]
        length_top_sec = maxDiff(x[np.where(y == max(y))])


        # If the equipment is not yet evaluated create a dict, else continue with the already present dict created earlier
        if equip.name not in self.depth_area[layer]["Area_yrange"][key]["equipment"].keys():
            self.depth_area[layer]["Area_yrange"][key]["equipment"][equip.name] = {}


        # Check to which class the equipment belongs and whether to install
        if isinstance(equip, Truck):
            if equip.install(
                layer=layer,
                grading_layer=grading,
                ymax=max_height,
                section_coords=coords,
            ):

                time = (
                    area / equip.get_production_rate(layer=layer, grading_layer=grading)

                )

                price = area *  equip.get_price(layer=layer, grading_layer=grading, key= cost_key)


                self.depth_area[layer]["Area_yrange"][key]["equipment"][equip.name][cost_key] = price
                self.depth_area[layer]["Area_yrange"][key]["equipment"][equip.name]['time'] = time

                self.depth_area[layer]["Area_yrange"][key]["color"] = "g"
                # Is the height of the layer higher than the current max layer? replace top and min_x values
                if end_lay > max_height:
                    max_height = end_lay
                    max_height_xcorner = xcorner_layer
                    length_top = length_top_sec

        elif isinstance(equip, PlateFeeder):

            if equip.install(
                layer=layer,
                grading_layer=grading,
                ymax=max_height,
                section_coords=coords,
            ):

                time = (
                    area / equip.get_production_rate(layer=layer, grading_layer=grading)

                )

                price = area *  equip.get_price(layer=layer, grading_layer=grading, key= cost_key)


                self.depth_area[layer]["Area_yrange"][key]["equipment"][equip.name][cost_key] = price
                self.depth_area[layer]["Area_yrange"][key]["equipment"][equip.name]['time'] = time

                self.depth_area[layer]["Area_yrange"][key]["color"] = "g"
                # Is the height of the layer higher than the current max layer? replace top and min_x values
                if end_lay > max_height:
                    max_height = end_lay
                    max_height_xcorner = xcorner_layer
                    length_top = length_top_sec

        elif isinstance(equip, Vessel):

            if equip.install(section_coords=coords, layer=layer, grading_layer=grading):
                time = (
                    area / equip.get_production_rate(layer=layer, grading_layer=grading)

                )
                price = equip.get_price(layer=layer, grading_layer=grading, key= cost_key) * area


                self.depth_area[layer]["Area_yrange"][key]["equipment"][equip.name][cost_key] = price
                self.depth_area[layer]["Area_yrange"][key]["equipment"][equip.name]['time'] = time

                self.depth_area[layer]["Area_yrange"][key]["color"] = "g"

                # Is the height of the layer higher than the current max layer? replace top and min_x values
                if end_lay > max_height:
                    max_height = end_lay
                    max_height_xcorner = xcorner_layer
                    length_top = length_top_sec

        elif isinstance(equip, Excavator):

            if equip.install(
                layer=layer,
                grading_layer=grading,
                ymax=max_height,
                section_coords=coords,
                xmax_top=max_height_xcorner,
                mass=mass,
                length_top=length_top,
            ):
                time = (
                    area / equip.get_production_rate(layer=layer, grading_layer=grading)

                )
                price = equip.get_price(layer=layer, grading_layer=grading, key= cost_key) * area



                self.depth_area[layer]["Area_yrange"][key]["equipment"][equip.name][cost_key] = price
                self.depth_area[layer]["Area_yrange"][key]["equipment"][equip.name]['time'] = time

                self.depth_area[layer]["Area_yrange"][key]["color"] = "g"

                if end_lay > max_height:
                    max_height = end_lay
                    max_height_xcorner = xcorner_layer
                    length_top = length_top_sec

        elif isinstance(equip, Crane):
            if equip.install(
                layer=layer,
                grading_layer=grading,
                ymax=max_height,
                xmax_top=max_height_xcorner,
                section_coords=coords,
                mass=mass,
            ):

                time = (
                    area / equip.get_production_rate(layer=layer, grading_layer=grading)

                )
                price = equip.get_price(layer=layer, grading_layer=grading, key= cost_key) * area



                self.depth_area[layer]["Area_yrange"][key]["equipment"][equip.name][cost_key] = price
                self.depth_area[layer]["Area_yrange"][key]["equipment"][equip.name]['time'] = time

                self.depth_area[layer]["Area_yrange"][key]["color"] = "g"

                if end_lay > max_height:
                    max_height = end_lay
                    max_height_xcorner = xcorner_layer
                    length_top = length_top_sec

        elif isinstance(equip, Barge):

            if equip.install(
                layer=layer,
                grading_layer=grading,
                slope=slope,
                section_coords=coords,
                xmax_top=max_height_xcorner,
                mass=mass,
            ):

                price = (
                    equip.get_price(layer=layer, grading_layer=grading, key= cost_key) * area
                )

                time = area / (
                    equip.get_production_rate(
                        layer=layer, grading_layer=grading
                    )

                )

                self.depth_area[layer]["Area_yrange"][key]["equipment"][equip.name][cost_key] = price
                self.depth_area[layer]["Area_yrange"][key]["equipment"][equip.name]['time'] = time

                self.depth_area[layer]["Area_yrange"][key]["color"] = "g"

                if end_lay > max_height:
                    max_height = end_lay
                    max_height_xcorner = xcorner_layer
                    length_top = length_top_sec

        return max_height, max_height_xcorner, length_top

    def equipment_dataframe(self):

        d = {
            (i, j): self.depth_area[i]["Area_yrange"][j]
            for i in self.depth_area.keys()
            for j in self.depth_area[i]["Area_yrange"].keys()
        }

        d2 = {}
        area = []
        for layer, section in d.keys():
            n1, n2 = round(float(section.split("-")[0]), 2), round(
                float(section.split("-")[1]), 2
            )
            d2[(layer, f"{n1}-{n2}")] = d[(layer, section)]['equipment']
            area.append(sum(d[(layer, section)]['area']))

        mux = pd.MultiIndex.from_tuples(d2.keys())
        df = pd.DataFrame(list(d2.values()), index=mux)
        df['area'] = area
        cols = [c for c in df.columns if c != 'area']
        cols.insert(0, 'area')
        df = df[cols]
        df = df.applymap(lambda x: None if x == {} else x)


        return df

    def optimal_equipment_set(self, dataframe, optimize_on ="cost"):

        """

        Parameters
        ----------
        dataframe
        optimize_on: str or list
            one or more of 'CO2', 'cost' or 'time'

        Returns
        -------

        """

        df = dataframe

        df.index = np.arange(0, len(df))

        # All the sections
        all_sections = df.index
        # All the equipment
        equipment = [c for c in df.columns if c != 'area']
        equipment_layers = {}
        # Which equipment can install which sections
        for i in range(len(equipment)):
            equipment_layers[equipment[i]] = []
            for index, row in df.iterrows():
                if bool(row[equipment[i]]):
                        equipment_layers[equipment[i]].append(index)

        combinations_dict = []
        combinations_lst = []

        all_equips_dict = {}

        for equip in equipment:

            i = 0
            installed = equipment_layers[equip].copy()
            equips_dict = {}
            equips_dict[equip] = installed.copy()
            equips_lst = [equip]

            # can't install a single layer
            if len(installed) != 0:
                while list(installed) != list(all_sections) and i < len(equipment):

                    rest_equip = sorted(equipment_layers.keys(),
                                        key=lambda k: len(set(equipment_layers[k]).difference(set(installed))),
                                        reverse= True)

                    most_new = rest_equip[0]

                    newlayers = [l for l in equipment_layers[most_new] if l not in installed]
                    installed.extend(newlayers)
                    equips_dict[rest_equip[0]] = newlayers
                    equips_lst.append(rest_equip[0])
                    installed = sorted(installed)
                    i += 1

                if list(installed) == list(all_sections):
                    equips_lst = sorted(equips_lst)
                    if equips_lst not in combinations_dict:
                        combinations_dict.append(equips_dict)
                        combinations_lst.append(equips_lst)

                all_equips_dict[f'{equips_lst}'] = equips_dict

        cost_combinations = {}

        for combi in combinations_dict:
            s = ''
            cost = 0
            CO2 = 0
            time = 0
            for key, value in combi.items():
                s += (key + ' + ')
                for section in value:
                    d = df[key].iloc[section]
                    if 'cost' in d:
                        cost += d['cost']
                    if 'CO2' in d:
                        CO2 += d['CO2']
                    if 'time' in d:
                        time += d['time']

            cost_combinations[s[:-2]] = {'cost': cost, 'CO2': CO2, 'time': time}

        df_optimal = pd.DataFrame.from_dict(cost_combinations, orient='index')
        df_optimal = df_optimal.sort_values(by= optimize_on)
        df_optimal['cost'] = df_optimal['cost'].round(2)
        df_optimal['CO2'] = df_optimal['CO2'].round(5)
        df_optimal['time'] = df_optimal['time'].round(5)



        return df_optimal

    def equipment_cost(self, *variants, equipment, cost, CO2, optimize_on, plot_error = True):
        """
        Compute the cost of the breakwater sections using equipment
        Parameters
        ----------
        id: str
            variant of breakwater
        areas:
            areas of the sections
        equipment: list
            list of equipments from the equipment class
        structure: dict
            structure of the breakwater
        cost: bool
            calculate the cost
        CO2: bool
            calculate the CO2 emission
        optimize_on: str or list
            one or more of 'CO2', 'cost' or 'time'
        Returns
        -------
            dictionary with the prices per layer
        Raises
        ------
        EquipmentError
            if the given equipment is not able to fill all the sections
        """

        if type(optimize_on) == str:
            if optimize_on not in ['CO2', 'cost', 'time']:
                raise KeyError('The optimize_on argument should be one of CO2, cost or time')
        elif type(optimize_on) == list:
            if not set(optimize_on).issubset(['CO2', 'cost', 'time']):
                raise KeyError('The optimize_on argument should be one of CO2, cost or time')
        else:
            raise KeyError('The argument optimize_on should be either list or string')

        cost_variables = []

        if cost:
            cost_variables.append('cost')
        if CO2:
            cost_variables.append('CO2')

        df_variants = {}

        variants = self._validate_variant(variants)

        optimal_equipment = {}

        for id in variants:
            depth_area = self._layers(id)
            slope = self._input_arguments["slope"]
            areas = self.area(id)
            structure = self.get_variant(id)
            self.depth_area = divide_cross_section(depth_area= depth_area, slope= slope)
            for cost_var in cost_variables:
                # Max height of already built sections and the corner coordinate
                max_height = 0
                max_height_xcorner = 0
                length_top = 0

                build_order = ['core', 'filter layer', 'underlayer', 'armour', 'toe']
                for layer in build_order:
                    if layer in areas.keys():
                        # What is the grading of the layer?
                        grading = None
                        if layer == "core":
                            grading = self._input_arguments['core_class']
                            mass = self.rho * (self._input_arguments['Dn50_core'])**3
                        # use armour units
                        elif layer == "armour" and self._input_arguments["armour"] != "Rock":
                            grading = self.structure['armour']['Mass']
                            mass = int(self.structure["armour"]["class"] * 2400 / 1000)

                        else:
                            grading = structure[layer]["class"]
                            mass = self.Grading.rosin_rammler(grading, 0.90) / 1000

                        for key, value in self.depth_area[layer]["Area_yrange"].items():
                            # Find the start and the end y-coordinate of the layer and the area of this section
                            start_lay, end_lay = float(key.split("-")[0]), float(key.split("-")[1])
                            area = self.depth_area[layer]["Area_yrange"][key]["area"]
                            coords = self.depth_area[layer]["Area_yrange"][key]["coordinates"]

                            for equip in equipment:
                                for p in range(0, len(coords), 2):
                                    # Left and right of the breakwater are not entirely symmetric so install seperate
                                    x = coords[p]
                                    y = coords[p + 1]
                                    if p == 0 and area[0] != 0:
                                        area1 = area[0]
                                    else:
                                        area1 = area[1]

                                    section_coords = list(zip(x, y))

                                    (
                                        max_height,
                                        max_height_xcorner,
                                        length_top,
                                    ) = self.install_all_equipments(
                                        coords=section_coords,
                                        equip=equip,
                                        layer=layer,
                                        grading=grading,
                                        mass=mass,
                                        end_lay=end_lay,
                                        area=area1,
                                        key=key,
                                        max_height=max_height,
                                        max_height_xcorner=max_height_xcorner,
                                        length_top=length_top,
                                        cost_key= cost_var
                                    )

            df_variants[id] = self.equipment_dataframe()
            cols = [c for c in df_variants[id].columns if c != 'area']
            df2 = df_variants[id][cols].dropna(how= 'all')

            # #if one row contains (exception of area column) only None the section can't be installed
            if len(df2) != len(df_variants[id]):
                raise EquipmentError( f"variant '{id}' can't be installed completely",
                                      other = self,
                                      id = id,
                                      areas = areas,
                                      depth_area = self.depth_area,
                                      plot_error = plot_error)
            else:
                optimal_equipment[id] = self.optimal_equipment_set(dataframe= df_variants[id],
                                                                   optimize_on = optimize_on)

        return optimal_equipment


    def cost_CO2_or_money(self, *variants, core_price, transport_cost, unit_price, output, key):
        """

        Parameters
        ----------
            *variants : str
            IDs of the variants to plot, see :py:attr:`variantIDs` for
            a list of all generated variants. If 'all' is in the
            arguments, all variants will be plotted.
        core_price : dict
            cost of the core material per m³
            {'cost': ... [EUR/m3], 'CO2': ... [kg/m3]}
        unit_price : dict
            the cost of an armour unit per m³ for each unit volume
            { V1: {'cost': ... [EUR/m3], 'CO2': ... [kg/m3]}}
        transport_cost : dict
            the cost to transport a m³ of rock from the quarry to the
            project location
            {'cost': ... [EUR/m3], 'CO2': ... [kg/m3]}
        output : {variant, layer, average}
            format of the output dict, variant returns the total cost
            of each variant, layer the cost of each layer for each
            variant and average returns the average cost.
        key: str
            Either CO2 or cost

        Returns
        -------

        """

        Grading = self._input_arguments['Grading']

        # set empty dict to store the output in
        cost = {}
        variant_price = {}
        # iterate over the generated variants
        for id in variants:
            # get the areas and structure of the variants
            areas = self.area(id)
            structure = self.get_variant(id)
            # iterate over the layers to price each layer
            for layer, area in areas.items():
                if layer is 'core':
                    # core is not included in the structure dict
                    price = ((core_price[key] + transport_cost[key])
                              * area)

                elif (self._input_arguments['armour'] is not 'Rock'
                        and layer is 'armour'):
                    # concrete armour units
                    V_armour = structure['armour']['class']
                    price = area * unit_price[V_armour][key]

                else:
                    # layer of the breakwater
                    rock_class = structure[layer]['class']

                    # get the price per meter
                    price = ((Grading[rock_class][key] + transport_cost[key])
                              * area)

                # add to dict
                variant_price[layer] = price

            # add to cost dict
            if output is 'variant' or output is 'average':
                # add total cost of all layers
                if key == 'cost':
                    cost[id] = np.round(np.sum(list(variant_price.values())), 2)
                else:
                    cost[id] = np.round(np.sum(list(variant_price.values())), 4)
            elif output is 'layer':
                # add the cost of each layer
                cost[id] = variant_price
            else:
                # invalid input
                raise NotSupportedError(
                    (f'Cost can\'t be exported as {output}, must be variant, '
                      'layer or average'))

        return cost

    def material_cost(
            self, *variants, core_price, unit_price, transport_cost,
            output='variant'):

        """ Compute the cost per meter for each variant
        Method to compute the cost of each generated variant, the cost
        is computed per meter. Cost can be either in terms of CO2 emission
        or monetary values

        Parameters
        ----------
        *variants : str
            IDs of the variants to plot, see :py:attr:`variantIDs` for
            a list of all generated variants. If 'all' is in the
            arguments, all variants will be plotted.
        core_price : dict
            cost of the core material per m³
            {'cost': ... [EUR/m3], 'CO2': ... [kg/m3]}
        unit_price : dict
            the cost of an armour unit per m³
            {'cost': ... [EUR/m3], 'CO2': ... [kg/m3]}
        transport_cost : dict
            the cost to transport a m³ of rock from the quarry to the
            project location
            {'cost': ... [EUR/m3], 'CO2': ... [kg/m3]}
        output : {variant, layer, average}
            format of the output dict, variant returns the total cost
            of each variant, layer the cost of each layer for each
            variant and average returns the average cost.
        Returns
        -------
        dict
            the cost
        Raises
        ------
        RockGradingError
            if no pricing is included in the given RockGrading
        """

        # check if transport cost have been given
        if transport_cost is None:
            # no cost thus set cost to zero
            transport_cost = {'cost': 0, 'CO2': 0}

        # validate variants
        variants = self._validate_variant(variants)

        # get the grading, and check if the cost has been added
        Grading = self._input_arguments['Grading']

        if 'cost' in Grading[list(Grading.grading.keys())[0]]:
            # pricing has been added
            pass
        else:
            # pricing has not been added, raise error
            raise RockGradingError('There is no pricing in the RockGrading')

        dicts = [core_price, transport_cost]

        for d in dicts:
            if type(d) != dict:
                raise KeyError('The arguments core_price, unit_price and transport_cost should be of the type dict and formatted as: formatted as:  {"cost: val[EUR/m3], "CO2": val[kg/m3]}')
            elif 'cost' not in d.keys() or 'CO2' not in d.keys():
                raise KeyError('The keys for core_price, unit_price and transport_cost should be "cost" and "CO2"')


        cost_material, CO2_material = {}, {}

        for id in variants:
            cost_material[id] = 0
            CO2_material[id] = 0
        
        if core_price['cost'] != None:
            
            cost_material = self.cost_CO2_or_money(*variants, core_price = core_price, transport_cost = transport_cost,
                                                   unit_price = unit_price, output = output, key = 'cost')
            # check if average must be computed
            if output is 'average':
                # compute average cost
                cost_material = {'average': np.round(np.average(list(cost_material.values())), 2)}
    
                # check if the cost have only been computed for 1 variant
                if len(variants) == 1:
                    # print user_warning and change key into variant
                    cost_material[variants[0]] = cost_material.pop('average')
                    user_warning(
                        ('Computing the average for one variantID, changed key '
                         'average in dict with the specified variantID'))
                    
        if core_price['CO2'] != None:
            
            CO2_material = self.cost_CO2_or_money(*variants, core_price = core_price, transport_cost = transport_cost,
                                                   unit_price = unit_price, output = output, key = 'CO2')
            # check if average must be computed
            if output is 'average':
                # compute average cost
                CO2_material = {'average': np.round(np.average(list(CO2_material.values())), 4)}
    
                # check if the cost have only been computed for 1 variant
                if len(variants) == 1:
                    # print user_warning and change key into variant
                    CO2_material[variants[0]] = CO2_material.pop('average')
                    user_warning(
                        ('Computing the average for one variantID, changed key '
                         'average in dict with the specified variantID'))

        return cost_material, CO2_material

    def total_cost(
        self,
        *variants,
        core_price,
        unit_price,
        equipment=None,
        transport_cost= None,
        optimize_on = ['cost', 'time'],
        output="variant",
        plot_error = True,
    ):
        """Compute the cost for either the material or CO2 footprint per meter for each variant

        Method to compute the cost of each generated variant, the cost
        is computed per meter

        Parameters
        ----------
        *variants : str
            IDs of the variants to plot, see :py:attr:`variantIDs` for
            a list of all generated variants. If 'all' is in the
            arguments, all variants will be plotted.
        equipment: lst
            list of equipment out of Equipment class
        core_price : dict
            cost of the core material per m³
            {'cost': ... [EUR/m3], 'CO2': ... [kg/m3]}
        unit_price : dict
            the cost of an armour unit per m³
            {'cost': ... [EUR/m3], 'CO2': ... [kg/m3]}
        transport_cost : dict
            the cost to transport a m³ of rock from the quarry to the
            project location
            {'cost': ... [EUR/m3], 'CO2': ... [kg/m3]}
        output : {variant, layer, average}
            format of the output dict, variant returns the total cost
            of each variant, layer the cost of each layer for each
            variant and average returns the average cost.

        Returns
        -------
        tuple
            the cost

        Raises
        ------
        RockGradingError
            if no pricing is included in the given RockGrading
        EquipmentError
            if the given equipment is not able to fill all the sections
        """
        # validate variants
        if transport_cost == None:
            transport_cost = {'cost': 0, 'CO2': 0}
            
        if type(optimize_on) == str:
            if optimize_on not in ['CO2', 'cost', 'time']:
                raise KeyError('The optimize_on argument should be one of CO2, cost or time')
        elif type(optimize_on) == list:
            if not set(optimize_on).issubset(['CO2', 'cost', 'time']):
                raise KeyError('The optimize_on argument should be one of CO2, cost or time')
        else:
            raise KeyError('The argument optimize_on should be either list or string')

        variants = self._validate_variant(variants)
        total_cost = {}
        total_CO2 = {}
        total_duration = {}
        opt_equipment = {}

        cost_material, CO2_material = self.material_cost(*variants, core_price= core_price, unit_price= unit_price,
                                                         transport_cost= transport_cost, output= output)


        for id in variants:
            total_cost[id] = cost_material[id]
            total_CO2[id] = CO2_material[id]
            total_duration[id] = 0
            opt_equipment[id] = None

        if equipment != None:

            calc_cost_equip, calc_CO2_equip = True, True

            for equip in equipment:
                if not isinstance(equip, Barge):
                    for key, item in equip.design_type.items():
                        if item['cost'] == None:
                            calc_cost_equip = False
                        if item['CO2'] == None:
                            calc_CO2_equip = False

            # if the variant can be installed pop from the total_cost dict
            for id in variants:

                try:
                    optimal_set = self.equipment_cost(id,
                                              equipment= equipment,
                                              cost = calc_cost_equip,
                                              CO2 = calc_CO2_equip,
                                              optimize_on = optimize_on,
                                              plot_error = plot_error)

                    cost_equip = optimal_set[id]['cost'].iloc[0]
                    CO2_equip = optimal_set[id]['CO2'].iloc[0]
                    duration = optimal_set[id]['time'].iloc[0]
                    op_equipment = optimal_set[id].index[0]

                except EquipmentError:
                    total_cost.pop(id)
                    total_CO2.pop(id)
                    total_duration.pop(id)
                    opt_equipment.pop(id)

                    user_warning(f'variant {id} can not be installed with the given equipment')

                    continue

                else:

                    total_cost[id] += cost_equip
                    total_CO2[id] += CO2_equip
                    total_duration[id] = duration
                    opt_equipment[id] = op_equipment


        return total_cost, total_CO2, total_duration, opt_equipment


    def get_variant(self, variantID):
        """Get the dimensions for the specified variant

        Parameters
        ----------
        variantID : str
            identifier of the variant, see :py:attr:`variantIDs` for a
            list of all generated variants.

        Returns
        -------
        dict
            Parameters and values (Dn50, Rock class) for one variant

        Raises
        ------
        KeyError
            If there is no variant with the given identifier
        """
        variant = {"armour": {}, "underlayer": {}}
        VariantCount = len(self.variantIDs)

        if variantID in self.variantIDs:
            if variantID == "a":
                key_u = 0
                key_f = 0
            elif variantID == "b":
                if VariantCount == 4:
                    key_u = 0
                    key_f = 1
                elif "filter layer" in self.structure:
                    key_f = 1
                    if VariantCount == 2:
                        key_u = 0
                    else:
                        key_u = 1
                elif VariantCount == 2:
                    key_u = 1
            elif variantID == "c":
                key_u = 1
                if VariantCount >= 3:
                    key_f = 2
            elif variantID == "d":
                key_u = 1
                key_f = 3
        else:
            raise KeyError(
                f"Variant with ID = {variantID} is not a variant, "
                f"generated variants are: {self.variantIDs}"
            )

        # add armour layer
        variant["armour"] = self.structure["armour"]

        # add underlayer
        underlayer = self.structure["underlayer"]

        for param in underlayer:
            if param == "state" or param == "layers":
                variant["underlayer"][param] = underlayer[param]
            else:
                variant["underlayer"][param] = underlayer[param][key_u]

        # add filter layer
        if "filter layer" in self.structure:
            filter_layer = self.structure["filter layer"]
            variant["filter layer"] = {}
            for i, param in enumerate(filter_layer):
                # check because values of these are not list
                if param == "state" or param == "layers":
                    # no need to use key_f, because constant for all
                    variant["filter layer"][param] = filter_layer[param]
                else:
                    # param is list and thus use key_f to get correct value
                    variant["filter layer"][param] = filter_layer[param][key_f]

            # check if valid filter layer or all None
            # slice because last two are state and layers
            sliced_values = list(variant["filter layer"].values())[:-2]
            if not any(sliced_values):
                # all values are None, thus delete filter layer
                del variant["filter layer"]

        # add toe
        if "toe" in self.structure:
            variant["toe"] = self.structure["toe"]

        return variant

    def plot(self, *variants, wlev=None, save_name=None, equipment= None):
        """Plot the cross section of the specified breakwater(s)

        Parameters
        ----------
        *variants : str
            IDs of the variants to plot, see :py:attr:`variantIDs` for
            a list of all generated variants. If 'all' is in the
            arguments, all variants will be plotted.
        wlev : str, optional, default: None
            label of the :py:class:`LimitState` from which the water
            level will be plotted. If no value is specified the water
            level from the normative limit state is used, which is the
            normative LimitState from the crest freeboard computation.
        save_name : str, optional, default: None
            if given the cross section is not shown but saved with the
            given name

        Raises
        ------
        InputError
            If no variants are specified or if the label of wlev is not
            a valid label of a :py:class:`LimitState`
        KeyError
            If there is no variant with the given identifier
        """
        # validate variants
        variants = self._validate_variant(variants)

        if wlev is None:
            wlev = self._state_overtopping
        else:
            for i, LimitState in enumerate(self._LimitStates):
                if LimitState.label == wlev:
                    wlev = i
                    break

        # check if wlev has been changed from LimitState label into index
        if isinstance(wlev, str):
            # wlev is still a string so not changed, which means that
            # the specified wlev is not a specified LimitState
            raise InputError("There is no LimitState with the given label")

        V, H = self._input_arguments["slope"]

        plt.figure(figsize=(10, 5))

        for i, id in enumerate(variants):
            # set subplot
            if len(variants) == 2:
                plt.subplot(1, 2, i + 1)
            elif len(variants) >= 3:
                plt.subplot(2, 2, i + 1)

            # get the coordinates
            coordinates = self._layers(id)

            # set xlim_max and xlim_min variable
            xlim_max, xlim_min = 0, 0

            # plot lines

            for layer, lines in coordinates.items():
                plt.plot(lines["x"], lines["y"], color="k")

                # check largest value for xlim
                if np.max(lines["x"]) >= xlim_max:
                    # set max as xlim_max
                    xlim_max = np.max(lines["x"])

                # check smallest value for xlim
                if np.min(lines["x"]) <= xlim_min:
                    # set min as xlim_min
                    xlim_min = np.min(lines["x"])

            # bottom + wlev
            x_wlev_max = 0.5 * self._input_arguments["B"] + H * self.Rc / V

            plt.axhline(y=0, color="k", linewidth=2)
            plt.hlines(
                y=self._LimitStates[wlev].h,
                xmin=xlim_min * 1.2,
                xmax=-x_wlev_max,
                color="b",
            )
            if not self._input_arguments["structure_type"] == "revetment":
                plt.hlines(
                    y=self._LimitStates[wlev].h,
                    xmin=x_wlev_max,
                    xmax=xlim_max * 1.2,
                    color="b",
                )

            # set xlim and ylim
            ymax = (self._LimitStates[self._state_overtopping].h + self.Rc) * 1.2
            plt.xlim(xlim_min * 1.2, xlim_max * 1.2)
            plt.ylim(-0.5, ymax)

            # add title to the plot
            if save_name is None:
                if isinstance(self.id, int):
                    title = "Cross section of rubble mound breakwater " f"{self.id}{id}"
                else:
                    title = f"Cross section of rubble mound breakwater {id}"
            else:
                name = save_name.split("/")[-1]
                title = f"Cross section of {name}"

            # get equipment to show in design explorer
            if equipment != None:
                plt.figtext(0.5, 0.01, f"Used equipment: {equipment}", ha="center", fontsize=10)



            plt.title(title)

            plt.gca().set_aspect("equal", adjustable="box")
            plt.grid()

        plt.tight_layout()

        if save_name is not None:
            plt.savefig(f"{save_name}.png")
            plt.close()
        else:
            plt.show()

    def print_variant(self, *variants, decimals=3):
        """Print the details for the specified variant(s)

        This method will print the computed Dn50, rock class, average
        Dn50 of the class, normative LimitState for all layers and
        specified variants.

        Parameters
        ----------
        *variants : str
            IDs of the variants to plot, see :py:attr:`variantIDs` for
            a list of all generated variants. If 'all' is in the
            arguments, all variants will be plotted.
        decimals : int, optional, default: 3
            number of decimal places to round to

        Raises
        ------
        InputError
            If no arguments are specified
        KeyError
            If there is no variant with the given identifier
        """
        # validate variants
        variants = self._validate_variant(variants)

        # loop over the variants
        for id in variants:
            variant = self.get_variant(id)

            # set the name of the table
            if isinstance(self.id, int):
                table_name = f"  Variant {self.id}{id}"
            else:
                table_name = f"  Variant {id}"
            print(table_name)
            table = []

            # iterate over the layers
            for i, (layer, dimensions) in enumerate(variant.items()):
                if i == 0:
                    # set headers of the table
                    headers = list(dimensions.keys())
                    headers.insert(0, "Layer")
                # add empty list for each row (layer)
                table.append([])
                table[i].append(layer)
                table[i].extend(dimensions.values())
                # print the label of the LimitState instead of the index
                state = dimensions["state"]
                index_state = headers.index("state")
                if isinstance(table[i][index_state], int):
                    table[i][index_state] = self._LimitStates[state].label

            print(
                tabulate(table, headers, tablefmt="github", floatfmt=(f".{decimals}f"))
            )
            print("")
            print(
                f"Rc = {np.round(self.Rc, decimals=decimals)} m, designed "
                f"with {self._LimitStates[self._state_overtopping].label} "
                "limit state"
            )
            print("\n")

    def print_logger(self, level="warnings"):
        """Print messages and warnings from the logger

        Parameters
        ----------
        msg_level : {'info', 'warnings'}, optional, default: 'warnings'
            specify print level, highest level is warnings and lowest
            level is info. Note that the info level will also print all
            warnings
        """
        # check if correct input has been given
        if level.lower() not in ["warnings", "info"]:
            raise NotSupportedError(
                f"{level} not implemented, must be info or warnings"
            )

        # print logger
        for type, messages in self.logger.items():
            if level.lower() == "warnings":
                if type == "INFO":
                    continue
            print(f"{type}:")
            if messages:
                for message in messages:
                    print(message)
            else:
                print(f"no {type.lower()} messages in log")
            print()

    def area(self, variantID):
        """Compute the area of all layers

        Method computes the area of each layer using Gauss's area
        formula. Which is given by the following formula:

        .. math::
           \\mathbf{A}=\\frac{1}{2} | \\sum_{i=1}^{n-1} x_{i} y_{i+1}
           +x_{n} y_{1}-\\sum_{i=1}^{n-1} x_{i+1} y_{i}-x_{1} y_{n} |

        Parameters
        ----------
        variantID : str
            identifier of the variant, see :py:attr:`variantIDs` for a
            list of all generated variants.

        Returns
        -------
        dict
            dict with the area of each layer
        """
        # get the coordinates of the layers
        coordinates = self._layers(variantID)

        # iterate over the layers to compute the area
        area = {}
        for layer, coord in coordinates.items():
            # get the x and y coordinates
            x = coord["x"]
            y = coord["y"]

            # use Gauss's area formula
            A = GaussianA(x, y)

            # add to area dict
            area[layer] = A

        return area


class RockRubbleMound(RubbleMound):
    """Design a breakwater with Rock as armour layer

    Makes a conceptual design for a conventional rubble mound breakwater
    with rock as the armour layer, for one or several limit states. The
    following computations are performed:

    - The armour layer is designed with the Van der Meer formulas for
      deep and shallow water (van der Meer, 1988; van Gent et al., 2003).
    - The underlayer is designed by using the rules for the underlayer
    - A filter layer is designed if one is needed, depends on
      :py:obj:`Dn50_core`
    - The toe is designed with the toe stability formula of
      Van der Meer (1998).
    - The crest freeboard is computed with the formula from EurOtop
      (2018)
    - The required width of the scour protection with Sumer and Fredsoe
      (2000)
    - If a :py:class:`Soil` is specified, a slip circle analysis is
      performed with :py:class:`Bishop`

    .. note::
       Depending on the input it might be that more rock classes are
       possible for the underlayer (and filter layer). In case the
       upper bound of the underlayer rule results in a different rock
       class as the lower bound, a new variant is generated. See
       :py:attr:`variantIDs` for a list of generated variants.

    .. note::
       The notional permeability, P, in the van der Meer formula is set
       to a constant value of 0.4. This is due to the fact that the
       substructure of the breakwater will always have an underlayer.

    Parameters
    ----------
    slope : tuple
        Slope of the armour layer (V, H). For example a slope of 3V:4H
        is defined as (3, 4)
    slope_foreshore : tuple
        slope of the foreshore (V, H). For example a slope of 1:100 is
        defined as (1, 100)
    rho_w : float
        density of water [kg/m³]
    B : float
        Crest width [m]
    N : int
        Number of incident waves at the toe of the structure [-]
    LimitState : :py:class:`LimitState` or list of :py:class:`LimitState`
        ULS, SLS or another limit state defined with
        :py:class:`LimitState`
    Grading : :py:class:`RockGrading`
        standard rock grading defined in the NEN-EN 13383-1 or a user
        defined rock grading
    core_material : dict
        core_class and Dn50
    safety : float, optional, default: 1
        safety factor of design (number of standard deviations from the
        mean)
    slope_toe : tuple, optional, default: (2,3)
        slope of the toe
    B_toe : float, optional, default: None
        width of the top of the toe in meters. By default the width of
        toe is taken as 3 * Dn50_toe.
    beta : float, optional, default: 0
        angle between direction of wave approach and a line normal to
        the breakwater (degrees).
    layers : int, optional, default: 2
        number of layers in the armour layer
    layers_underlayer : int, optional, default: 2
        number of layers in the underlayer
    vdm : {min, max, avg}, optional, default: max
        value to return in case both the deep and shallow water formula
        are valid. min for the lowest value, max for the highest value
        and avg for the average value, default is max.
    Soil : :py:class:`Soil`, optional, default: None
        by default Soil is None, which means that the geotechnical checks
        are not performed. By specifying a Soil object, the geotechnical
        checks are automatically performed.
    phi : float, optional, default: 40
        internal friction angle of rock [degrees]
    id : int, optional, default: None
        add a unique id to the breakwater

    Attributes
    ----------
    logger : dict
        dict of warnings and messages
    structure : dict
        dictionary with the computed Dn50, rock class, and average Dn50
        of the rock class for each layer and the toe. This dictionary
        includes all variants, use :py:meth:`get_variant` to get the
        parameters of one specific variant. Alternatively,
        :py:meth:`print_variant` can be used to print the details of
        one, multiple or all variants.
    alpha : float
        slope of the structure in radians
    id : int
        unique id of the breakwater
    variantIDs : list
        list with the IDs of the variants generated for this rubble
        mound breakwater.
    Rc : float
        the crest freeboard of the structure [m]
    width_scour : float
        the required length of the scour protection [m]
    """

    def __init__(
        self,
        slope,
        slope_foreshore,
        rho_w,
        B,
        N,
        LimitState,
        Grading,
        core_material,
        safety=1,
        slope_toe=(2, 3),
        structure_type="breakwater",
        B_toe=None,
        beta=0,
        layers=2,
        layers_underlayer=2,
        vdm="max",
        Soil=None,
        phi=40,
        id=None,
        **kwargs,
    ):
        """See help(RockRubbleMound) for more info"""
        # set logger and structure
        self.logger = {"INFO": [], "WARNING": []}
        self.structure = {}

        # compute angles
        self.alpha = np.arctan(slope[0] / slope[1])
        slope_foreshore = np.arctan(slope_foreshore[0] / slope_foreshore[1])

        # set id and variantIDs
        self.id = id
        self.variantIDs = ["a", "b", "c", "d"]

        # compute relative buoyant density
        delta_rock = (Grading.rho - rho_w) / rho_w

        # set attribute for vandermeer to check validity
        # notional permeability is fixed due to substructure
        self._vandermeer = {"N": N, "P": 0.4, "Delta": delta_rock}

        # convert single LimitState to list if needed
        if isinstance(LimitState, list):
            LimitStates = LimitState
        else:
            LimitStates = [LimitState]

        # set temporary values to check changes
        Dn50, state = 0, 0

        for i, LimitState in enumerate(LimitStates):
            # design armour layer
            Dn50_temp = vandermeer(
                LimitState=LimitState,
                Delta=delta_rock,
                P=0.4,
                N=N,
                alpha=self.alpha,
                slope_foreshore=slope_foreshore,
                val=vdm,
                safety=safety,
                logger=self.logger,
            )

            # check if computed Dn50 of current LimitState is larger
            # than current normative Dn50
            if Dn50_temp > Dn50:
                # if larger the normative Dn50 must be changed
                Dn50 = Dn50_temp
                state = i

        class_armour = Grading.get_class(Dn50)

        # check dn85/dn15 range with rosin_rammler
        M15 = Grading.rosin_rammler(class_=class_armour, y=0.15)
        dn15 = (M15 / Grading.rho) ** (1 / 3)
        M85 = Grading.rosin_rammler(class_=class_armour, y=0.85)
        dn85 = (M85 / Grading.rho) ** (1 / 3)
        self._vandermeer["gradation"] = dn85 / dn15

        class_Dn50 = Grading.get_class_dn50(class_armour)
        self.structure["armour"] = {
            "computed Dn50": Dn50,
            "computed Dn50": Dn50,
            "class": class_armour,
            "class Dn50": class_Dn50,
            "state": state,
            "layers": layers,
        }

        # design underlayer, filter layer and crest height
        super().__init__(
            Dn50=class_Dn50,
            core_material=core_material,
            rho=Grading.rho,
            rho_w=rho_w,
            armour_layer="Rock",
            layers=layers,
            LimitStates=LimitStates,
            Grading=Grading,
            safety=safety,
            layers_underlayer=layers_underlayer,
            slope_toe=slope_toe,
            structure_type=structure_type,
            B_toe=B_toe,
            slope=slope,
            B=B,
            beta=beta,
            id=id,
            Soil=Soil,
            phi=phi,
            **kwargs,
        )

    def __str__(self):
        return (
            f"id.{self.id}: breakwater with rock as armour layer, "
            f"and variants: {self.variantIDs}"
        )

    def cost(
        self,
        *variants,
        equipment=None,
        core_price=None,
        transport_cost,
        output="variant",
    ):
        """Compute the cost for the material or the CO2 footprint per meter for each variant

        Method to compute the cost of each generated variant, the cost
        is computed per meter. The cost of the rocks must be specified
        in the RockGrading. If transport cost are not included in the
        price of rocks or core_price it can be given with the argument
        transport_cost.

        Parameters
        ----------
        *variants : str
            IDs of the variants to plot, see :py:attr:`variantIDs` for
            a list of all generated variants. If 'all' is in the
            arguments, all variants will be plotted.
        type: {'Material, 'CO2'}
            Indicate whether the costs are calculated for the material or the CO2
        transport_cost : float, optional, default: None
            the cost to transport a m³ of rock from the quarry to the
            project location
        output : {variant, layer, average}
            format of the output dict, variant returns the total cost
            of each variant, layer the cost of each layer for each
            variant and average returns the average cost.

        Returns
        -------
        dict
            the cost

        Raises
        ------
        RockGradingError
            if no pricing is included in the given RockGrading
        """
        # compute the cost of the concept
        cost = self.total_cost(
            *variants,
            equipment=equipment,
            core_price=core_price,
            unit_price=0,
            transport_cost=transport_cost,
            output=output,
        )

        return cost

    def check_validity(self, decimals=3):
        """Check if the used parameters are within the validity range

        The Van der Meer equations for deep and shallow water are
        empirical equations, meaning that they are based on experiments.
        Therefore, the Van der Meer equations are, strictly speaking,
        only valid if the parameters are within the range of the
        parameters from the experiments. This method prints a table from
        which the validity range, used value and if the parameter is
        within validity range can be read.

        Parameters
        ----------
        decimals : int, optional, default: 3
            number of decimals
        """
        # set table
        headers = ["Parameter    ", "Validity range", "Used value", "in range?"]
        table_vdm, table_toe = [], []

        # get the normative LimitState
        state = self.structure["armour"]["state"]
        state_toe = self.structure["toe"]["state"]
        LimitState = self._LimitStates[state]

        # check the validity ranges of the van der meer formula
        for info in self.logger["INFO"]:
            if LimitState.label in info:
                # msg from the logger with the normative LimitState
                if "vandermeer_deep" in info:
                    # vandermeer_deep was used
                    print(
                        (
                            " Range of validity of parameters in deep water "
                            "formulae by Van der Meer"
                        )
                    )

                    # check the slope of the structure
                    slope = self._input_arguments["slope"]
                    table_vdm.append(
                        ["tan(alpha)", "1:6 - 1:1.5", f"{slope[0]}:{slope[1]}"]
                    )
                    if slope[0] / slope[1] >= 1 / 6 and slope[0] / slope[1] <= 1 / 1.5:
                        table_vdm[0].append("yes")
                    else:
                        table_vdm[0].append("NO")

                    # check the number of waves
                    N = self._vandermeer["N"]
                    table_vdm.append(["N", "< 7500", np.round(N, decimals)])
                    if N < 7500:
                        table_vdm[1].append("yes")
                    else:
                        table_vdm[1].append("NO")

                    # check wave steepness s_om (based on Tm)
                    s_om = LimitState.s(number="mean")
                    table_vdm.append(["s_om", "0.01 - 0.06", np.round(s_om, decimals)])
                    if s_om >= 0.01 and s_om <= 0.06:
                        table_vdm[2].append("yes")
                    else:
                        table_vdm[2].append("NO")

                    # surf_similarity parameter based on Tm
                    xi_m = LimitState.surf_similarity(alpha=self.alpha, number="mean")
                    table_vdm.append(["xi_m", "0.7 - 7", np.round(xi_m, decimals)])
                    if xi_m >= 0.7 and xi_m <= 7:
                        table_vdm[3].append("yes")
                    else:
                        table_vdm[3].append("NO")

                    # delta
                    delta = self._vandermeer["Delta"]
                    table_vdm.append(["Delta", "1 - 2.1", np.round(delta, decimals)])
                    if delta >= 1 and delta <= 2.1:
                        table_vdm[4].append("yes")
                    else:
                        table_vdm[4].append("NO")

                    # relative water depth at the toe
                    Hs = LimitState.get_Hs(definition="H13")
                    h = LimitState.h
                    table_vdm.append(["h/Hs", "> 3", np.round(h / Hs, decimals)])
                    if h / Hs > 3:
                        table_vdm[5].append("yes")
                    else:
                        table_vdm[5].append("NO")

                    # notional permeability
                    P = self._vandermeer["P"]
                    table_vdm.append(["P", "0.1 - 0.6", np.round(P, decimals)])
                    if P >= 0.1 and P <= 0.6:
                        table_vdm[6].append("yes")
                    else:
                        table_vdm[6].append("NO")

                    # armourstone gradation
                    gradation = self._vandermeer["gradation"]
                    table_vdm.append(
                        ["Dn85/Dn15", "< 2.5", np.round(gradation, decimals)]
                    )
                    if gradation < 2.5:
                        table_vdm[7].append("yes")
                    else:
                        table_vdm[7].append("NO")

                    # damage-storm duration ratio
                    Sd = LimitState["Sd"]
                    damage_ratio = Sd / np.sqrt(N)
                    table_vdm.append(
                        ["Sd/sqrt(N)", "< 0.9", np.round(damage_ratio, decimals)]
                    )
                    if damage_ratio < 0.9:
                        table_vdm[8].append("yes")
                    else:
                        table_vdm[8].append("NO")

                    # stability number
                    Dn50 = self.structure["armour"]["computed Dn50"]
                    stability = Hs / (delta * Dn50)
                    table_vdm.append(
                        ["Hs/(Delta*Dn50)", "1 - 4", np.round(stability, decimals)]
                    )
                    if stability >= 1 and stability <= 4:
                        table_vdm[9].append("yes")
                    else:
                        table_vdm[9].append("NO")

                    # damage level parameter
                    table_vdm.append(["Sd", "1 - 20", np.round(Sd, decimals)])
                    if Sd > 1 and Sd < 30:
                        table_vdm[10].append("yes")
                    else:
                        table_vdm[10].append("NO")

                elif "vandermeer_shallow" in info:
                    # vandermeer_shallow was used
                    print(
                        (
                            " Range of validity of parameters in shallow water "
                            "formulae by Van der Meer"
                        )
                    )

                    # check the slope of the structure
                    slope = self._input_arguments["slope"]
                    table_vdm.append(
                        ["tan(alpha)", "1:4 - 1:2", f"{slope[0]}:{slope[1]}"]
                    )
                    if slope[0] / slope[1] >= 1 / 4 and slope[0] / slope[1] <= 1 / 2:
                        table_vdm[0].append("yes")
                    else:
                        table_vdm[0].append("NO")

                    # check the number of waves
                    N = self._vandermeer["N"]
                    table_vdm.append(["N", "< 3000", np.round(N, decimals)])
                    if N < 3000:
                        table_vdm[1].append("yes")
                    else:
                        table_vdm[1].append("NO")

                    # check wave steepness s_om (based on Tm)
                    s_om = LimitState.s(number="mean")
                    table_vdm.append(["s_om", "0.01 - 0.06", np.round(s_om, decimals)])
                    if s_om >= 0.01 and s_om <= 0.06:
                        table_vdm[2].append("yes")
                    else:
                        table_vdm[2].append("NO")

                    # surf_similarity parameter based on Tm
                    xi_m = LimitState.surf_similarity(alpha=self.alpha, number="mean")
                    table_vdm.append(["xi_m", "1 - 5", np.round(xi_m, decimals)])
                    if xi_m >= 1 and xi_m <= 5:
                        table_vdm[3].append("yes")
                    else:
                        table_vdm[3].append("NO")

                    # surf_similarity parameter based on Tm-1.0
                    xi_m_min_1 = LimitState.surf_similarity(
                        alpha=self.alpha, number="spectral"
                    )
                    table_vdm.append(
                        ["xi_m_min_1", "1.3 - 6.5", np.round(xi_m_min_1, decimals)]
                    )
                    if xi_m_min_1 >= 1.3 and xi_m_min_1 <= 6.5:
                        table_vdm[4].append("yes")
                    else:
                        table_vdm[4].append("NO")

                    # wave height ratio
                    Hs = LimitState.get_Hs(definition="H13")
                    wave_ratio = LimitState["H2_per"] / Hs
                    table_vdm.append(
                        ["H2%/Hs", "1.2 - 1.4", np.round(wave_ratio, decimals)]
                    )
                    if wave_ratio >= 1.2 and wave_ratio <= 1.4:
                        table_vdm[5].append("yes")
                    else:
                        table_vdm[5].append("NO")

                    # Deep-water wave height ratio
                    h = LimitState.h
                    table_vdm.append(["Ho/h", "0.25 - 1.5"])
                    try:
                        deep_water_ratio = LimitState["Ho"] / h
                        table_vdm[6].append(np.round(deep_water_ratio, decimals))
                        if deep_water_ratio >= 0.25 and deep_water_ratio <= 1.5:
                            table_vdm[6].append("yes")
                        else:
                            table_vdm[6].append("NO")
                    except KeyError:
                        table_vdm[6].append(f"Ho not in {LimitState.label}")

                    # armourstone gradation
                    gradation = self._vandermeer["gradation"]
                    table_vdm.append(
                        ["Dn85/Dn15", "1.4 - 2.0", np.round(gradation, decimals)]
                    )
                    if gradation >= 1.4 and gradation <= 2:
                        table_vdm[7].append("yes")
                    else:
                        table_vdm[7].append("NO")

                    # core material - armour ratio
                    Dn50 = self.structure["armour"]["computed Dn50"]
                    Dn50_core = self._input_arguments["Dn50_core"]
                    ratio_core = Dn50_core / Dn50
                    table_vdm.append(
                        ["Dn50-core/Dn50", "0 - 0.3", np.round(ratio_core, decimals)]
                    )
                    if ratio_core <= 0.3:
                        table_vdm[8].append("yes")
                    else:
                        table_vdm[8].append("NO")

                    # stability number
                    delta = self._vandermeer["Delta"]
                    stability = Hs / (delta * Dn50)
                    table_vdm.append(
                        ["Hs/(Delta*Dn50)", "0.5 - 4.5", np.round(stability, decimals)]
                    )
                    if stability >= 0.5 and stability <= 4.5:
                        table_vdm[9].append("yes")
                    else:
                        table_vdm[9].append("NO")

                    # damage level parameter
                    Sd = LimitState["Sd"]
                    table_vdm.append(["Sd", "< 30", Sd])
                    if Sd < 30:
                        table_vdm[10].append("yes")
                    else:
                        table_vdm[10].append("NO")

        # check validity ranges of toe formula
        # check ratio water depth above toe over water depth
        Dn50_toe = self.structure["toe"]["computed Dn50"]
        htoe = np.ceil(self._estimate_htoe() / Dn50_toe) * Dn50_toe
        ht = h - htoe

        ratio_depth_toe = ht / h
        table_toe.append(["ht/h", "0.4 - 0.9", np.round(ratio_depth_toe, decimals)])
        if ratio_depth_toe >= 0.4 and ratio_depth_toe <= 0.9:
            table_toe[0].append("yes")
        else:
            table_toe[0].append("NO")

        # ratio ht over Dn50 of the toe
        ratio_dn = ht / Dn50_toe
        table_toe.append(["ht/Dn50", "3 - 25", np.round(ratio_dn, decimals)])
        if ratio_dn >= 3 and ratio_dn <= 25:
            table_toe[1].append("yes")
        else:
            table_toe[1].append("NO")

        # print the table
        print(tabulate(table_vdm, headers, tablefmt="github"))
        print()
        print(" Range of validity of parameters in toe formula")
        print(tabulate(table_toe, headers, tablefmt="github"))


class ConcreteRubbleMound(RubbleMound):
    """Design a breakwater with concrete armour units as armour layer

    Makes a conceptual design for a conventional rubble mound breakwater
    with armour units as the armour layer, for one or several limit
    states. The following computations are performed:

    - The armour layer is designed with the Hudson formula (Hudson, 1959),
      with :math:`H_{1/3}` as wave height, in accordance with the design
      guidelines for Xbloc and XblocPlus (Delta Marine Consultants, 2018).
    - The underlayer is designed by using the rules for the underlayer
    - A filter layer is designed if one is needed, depends on
      :py:obj:`Dn50_core`
    - The toe is designed with the toe stability formula of
      Van der Meer (1998).
    - The crest freeboard is computed with the formula from EurOtop
      (2018)
    - The required width of the scour protection with Sumer and Fredsoe
      (2000)
    - If a :py:class:`Soil` is specified, a slip circle analysis is
      performed with :py:class:`Bishop`

    .. note::
       Depending on the input it might be that more rock classes are
       possible for the underlayer (and filter layer). In case the
       upper bound of the underlayer rule results in a different rock
       class as the lower bound, a new variant is generated. See
       :py:attr:`variantIDs` for a list of generated variants.

    .. warning::
       Currently only the rules for the underlayer of Rock, Xbloc and
       XblocPlus have been implemented. However, it is possible to
       design with another type of armour units. In this case the filter
       rule must manually be set to Rock, Xbloc or XblocPlus.

    Parameters
    ----------
    slope : tuple
        Slope of the armour layer (V, H). For example a slope of 3V:4H
        is defined as (3, 4)
    slope_foreshore : tuple
        slope of the foreshore (V, H). For example a slope of 1:100 is
        defined as (1, 100)
    B : float
        Crest width [m]
    rho_w : float
        density of water [kg/m³]
    LimitState : :py:class:`LimitState` or list of :py:class:`LimitState`
        ULS, SLS or another limit state defined with
        :py:class:`LimitState`
    ArmourUnit : obj
        armour unit class which inherits from :py:class:`ConcreteArmour`,
        for instance :py:class:`Xbloc` or :py:class:`XblocPlus`
    Grading : :py:class:`RockGrading`
        standard rock grading defined in the NEN-EN 13383-1 or a user
        defined rock grading
    core_material : dict
        core_class and Dn50
    safety : float, optional, default: 1
        safety factor of design (number of standard deviations from the
        mean)
    slope_toe : tuple, optional, default: (2,3)
        slope of the toe
    B_toe : float, optional, default: None
        width of the top of the toe in meters. By default the width of
        toe is taken as 3 * Dn50_toe.
    beta : float, optional, default: 0
        angle between direction of wave approach and a line normal to
        the breakwater (degrees).
    layers : int, optional, default: 1
        number of layers in the armour layer
    layers_underlayer : int, optional, default: 2
        number of layers in the underlayer
    filter_rule : {'Rock', 'Xbloc', 'XblocPlus'}, optional, default: None
        filter rule to use for the substructure of the breakwater, for
        Rock, Xbloc and XblocPlus the correct filter rule is
        automatically selected. In case another type of armour layer is
        used one of these filter rules must be chosen.
    Soil : :py:class:`Soil`, optional, default: None
        by default Soil is None, which means that the geotechnical checks
        are not performed. By specifying a Soil object, the geotechnical
        checks are automatically performed.
    phi : float, optional, default: 40
        internal friction angle of rock [degrees]
    id : int, optional, default: None
        add a unique id to the breakwater

    Attributes
    ----------
    logger : dict
        dict of warnings and messages
    structure : dict
        dictionary with the computed Dn50, rock class, and average Dn50
        of the rock class for each layer and the toe. This dictionary
        includes all variants, use :py:meth:`get_variant` to get the
        parameters of one specific variant. Alternatively,
        :py:meth:`print_variant` can be used to print the details of
        one, multiple or all variants.
    alpha : float
        slope of the structure in radians
    id : int
        unique id of the breakwater
    variantIDs : list
        list with the IDs of the variants generated for this rubble
        mound breakwater.
    Rc : float
        the crest freeboard of the structure [m]
    width_scour : float
        the required length of the scour protection [m]
    """

    def __init__(
        self,
        slope,
        slope_foreshore,
        B,
        rho_w,
        LimitState,
        ArmourUnit,
        Grading,
        core_material,
        safety=1,
        slope_toe=(2, 3),
        structure_type="breakwater",
        B_toe=None,
        beta=0,
        layers=1,
        layers_underlayer=2,
        filter_rule=None,
        Soil=None,
        phi=40,
        id=None,
        **kwargs,
    ):
        """See help(ConcreteRubbleMound) for more info"""
        # set logger and structure
        self.logger = {"INFO": [], "WARNING": []}
        self.structure = {}

        # Density
        self.rho = ArmourUnit.rho

        # compute angles
        self.alpha = np.arctan(slope[0] / slope[1])
        slope_foreshore = np.arctan(slope_foreshore[0] / slope_foreshore[1])

        # compute relative buoyant density
        self.id = id
        self.variantIDs = ["a", "b", "c", "d"]

        # compute relative buoyant density
        delta_armour = (ArmourUnit.rho - rho_w) / rho_w

        # convert single LimitState to list if needed
        if isinstance(LimitState, list):
            LimitStates = LimitState
        else:
            LimitStates = [LimitState]

        # determine the type of armour layer
        armour_layer = ArmourUnit.name
        if armour_layer in ["Xbloc", "XblocPlus"]:
            # the formula for Xbloc and XblocPlus is based on the
            # Hudson formula with a fixed slope of 3:4
            angle = np.arctan(3 / 4)
            filter_rule = armour_layer
        else:
            angle = self.alpha

        # raise error if no filter rule is specified
        if filter_rule is None:
            supported_rules = ", ".join(substructure._supported_armour_layers())
            raise NotSupportedError(
                (
                    f"Filter rule for {armour_layer} is not implemented, set "
                    f"filter rule to use with filter_rule to {supported_rules}"
                )
            )

        # set temporary values to check changes
        dn, state = 0, 0

        for i, LimitState in enumerate(LimitStates):
            # unpack Hydraulic Conditions
            Hs = LimitState.get_Hs(definition="H13")

            # design armour layer
            dn_temp = hudson(H=Hs, Kd=ArmourUnit.kd, Delta=delta_armour, alpha=angle)

            # check if computed Dn50 of current LimitState is larger
            # than current normative Dn50
            if dn_temp > dn:
                # if larger the normative Dn50 must be changed
                dn = dn_temp
                state = i

        V_required = ArmourUnit.get_class(dn)

        self.structure["armour"] = {
            "computed Dn50": dn,
            "class": V_required,
            "class Dn50": V_required ** (1 / 3),
            "state": state,
            "layers": layers,
            'Mass': ArmourUnit.units[V_required]['M']
        }

        # design underlayer, filter layer and crest height
        super().__init__(
            Dn50=V_required ** (1 / 3),
            core_material=core_material,
            rho=ArmourUnit.rho,
            rho_w=rho_w,
            armour_layer=armour_layer,
            layers=layers,
            LimitStates=LimitStates,
            Grading=Grading,
            safety=safety,
            layers_underlayer=layers_underlayer,
            slope_toe=slope_toe,
            structure_type=structure_type,
            B_toe=B_toe,
            slope=slope,
            B=B,
            beta=beta,
            id=id,
            filter_rule=filter_rule,
            Soil=Soil,
            phi=phi,
            **kwargs,
        )

        # check if a correction factor must be applied for Xbloc or
        # XblocPlus (DMC design guidelines for Xbloc and XblocPlus)
        correction = getattr(ArmourUnit, "correction_factor", None)
        if callable(correction):
            # make dict to pass to method correction_factor
            params = {
                "Hs": LimitStates[state].get_Hs("H13"),
                "h": LimitStates[state].h,
                "Rc": self.Rc,
                "occurrence_hs": False,
                "slope": self.alpha,
                "slope_foreshore": slope_foreshore,
                "permeability": "permeable",
                "Dn": dn,
                "layers": layers,
                "B": B,
                "beta": beta,
            }

            correction_factor = ArmourUnit.correction_factor(
                logger=self.logger, **params
            )

        else:
            self.logger["INFO"].append(
                "Given ArmourUnit did not have a method to compute a "
                "correction factor. See documentation how to define a method "
                "to compute a correction factor for a custom armour unit"
            )

            correction_factor = 1.0

        # check if a correction factor must be applied
        if correction_factor != 1:
            new_dn = correction_factor ** (1 / 3) * dn

            # replave values in structure with new ones for armour
            new_class = ArmourUnit.get_class(new_dn)
            self.structure["armour"] = {
                "computed Dn50": new_dn,
                "class": new_class,
                "class Dn50": new_class ** (1 / 3),
                "state": state,
                "layers": layers,
                'Mass': ArmourUnit.units[V_required]['M']
            }

            # check if the class of the units has changed
            if new_class != V_required:
                self.logger["INFO"].append(
                    "correction factor resulted in new class, so design will "
                    "be changed"
                )

                # because of the new class the design for the underlayer
                # and filter layer must be changed as well
                super().__init__(
                    Dn50=new_class ** (1 / 3),
                    core_material=core_material,
                    rho=ArmourUnit.rho,
                    rho_w=rho_w,
                    armour_layer=armour_layer,
                    layers=layers,
                    LimitStates=LimitStates,
                    Grading=Grading,
                    safety=safety,
                    layers_underlayer=layers_underlayer,
                    slope_toe=slope_toe,
                    structure_type=structure_type,
                    B_toe=B_toe,
                    slope=slope,
                    B=B,
                    beta=beta,
                    id=id,
                    filter_rule=filter_rule,
                    Soil=Soil,
                    phi=phi,
                    **kwargs,
                )

            else:
                # same class so design does not have to be changed
                self.logger["INFO"].append(
                    "correction factor did not result in a new class, so "
                    "design will not be changed"
                )

    def __str__(self):
        return (
            f"id.{self.id}: breakwater with armour units as armour layer, "
            f"and variants: {self.variantIDs}"
        )

    def cost(
        self,
        *variants,
        equipment=None,
        core_price=None,
        unit_price=None,
        transport_cost,
        output="variant",
    ):
        """Compute the cost per meter for each variant for the materials or the CO2 footprint

        Method to compute the cost of each generated variant, the cost
        is computed per meter. The cost of the rocks in the substructure
        must be specified in the RockGrading. If transport cost are not
        included in the price of rocks or core_price it can be given
        with the argument transport_cost.

        .. note::
           The transport_cost are not added to the price of the armour
           layer. The assumption has been made that the cost of
           producing and transporting the armour units is included in
           the unit_price.

        Parameters
        ----------
        *variants : str
            IDs of the variants to plot, see :py:attr:`variantIDs` for
            a list of all generated variants. If 'all' is in the
            arguments, all variants will be plotted.
        type: {'Material, 'CO2'}
            Indicate whether the costs are calculated for the material or the CO2
        unit_price : float
            the cost of an armour unit per m³
        transport_cost : float, optional, default: None
            the cost to transport a m³ of rock from the quarry to the
            project location
        output : {variant, layer, average}
            format of the output dict, variant returns the total cost
            of each variant, layer the cost of each layer for each
            variant and average returns the average cost.

        Returns
        -------
        dict
            the cost

        Raises
        ------
        RockGradingError
            if no pricing is included in the given RockGrading
        """
        # compute the cost of the concept
        cost = self.total_cost(
            *variants,
            equipment=equipment,
            core_price=core_price,
            unit_price=unit_price,
            transport_cost=transport_cost,
            output=output,
        )

        return cost


class ConcreteRubbleMoundRevetment(RubbleMound):

    """Design a revetment with concrete armour units as armour layer

    Makes a conceptual design for a conventional rubble mound revetment
    with armour units as the armour layer, for one or several limit
    states. The following computations are performed:

    - The armour layer is designed with the Hudson formula (Hudson, 1959),
      with :math:`H_{1/3}` as wave height, in accordance with the design
      guidelines for Xbloc and XblocPlus (Delta Marine Consultants, 2018).
    - The underlayer is designed by using the rules for the underlayer
    - A filter layer is designed if one is needed, depends on
      :py:obj:`Dn50_core`
    - The toe is designed with the toe stability formula of
      Van der Meer (1998).
    - The crest freeboard is computed with the formula from EurOtop
      (2018)
    - The required width of the scour protection with Sumer and Fredsoe
      (2000)
    - If a :py:class:`Soil` is specified, a slip circle analysis is
      performed with :py:class:`Bishop`

    .. note::
       Depending on the input it might be that more rock classes are
       possible for the underlayer (and filter layer). In case the
       upper bound of the underlayer rule results in a different rock
       class as the lower bound, a new variant is generated. See
       :py:attr:`variantIDs` for a list of generated variants.

    .. warning::
       Currently only the rules for the underlayer of Rock, Xbloc and
       XblocPlus have been implemented. However, it is possible to
       design with another type of armour units. In this case the filter
       rule must manually be set to Rock, Xbloc or XblocPlus.

    Parameters
    ----------
    slope : tuple
        Slope of the armour layer (V, H). For example a slope of 3V:4H
        is defined as (3, 4)
    slope_foreshore : tuple
        slope of the foreshore (V, H). For example a slope of 1:100 is
        defined as (1, 100)
    B : float
        Crest width [m]
    rho_w : float
        density of water [kg/m³]
    LimitState : :py:class:`LimitState` or list of :py:class:`LimitState`
        ULS, SLS or another limit state defined with
        :py:class:`LimitState`
    ArmourUnit : obj
        armour unit class which inherits from :py:class:`ConcreteArmour`,
        for instance :py:class:`Xbloc` or :py:class:`XblocPlus`
    Grading : :py:class:`RockGrading`
        standard rock grading defined in the NEN-EN 13383-1 or a user
        defined rock grading
    core_material : dict
        core_class and Dn50
    safety : float, optional, default: 1
        safety factor of design (number of standard deviations from the
        mean)
    slope_toe : tuple, optional, default: (2,3)
        slope of the toe
    B_toe : float, optional, default: None
        width of the top of the toe in meters. By default the width of
        toe is taken as 3 * Dn50_toe.
    beta : float, optional, default: 0
        angle between direction of wave approach and a line normal to
        the breakwater (degrees).
    layers : int, optional, default: 1
        number of layers in the armour layer
    layers_underlayer : int, optional, default: 2
        number of layers in the underlayer
    filter_rule : {'Rock', 'Xbloc', 'XblocPlus'}, optional, default: None
        filter rule to use for the substructure of the breakwater, for
        Rock, Xbloc and XblocPlus the correct filter rule is
        automatically selected. In case another type of armour layer is
        used one of these filter rules must be chosen.
    Soil : :py:class:`Soil`, optional, default: None
        by default Soil is None, which means that the geotechnical checks
        are not performed. By specifying a Soil object, the geotechnical
        checks are automatically performed.
    phi : float, optional, default: 40
        internal friction angle of rock [degrees]
    id : int, optional, default: None
        add a unique id to the breakwater

    Attributes
    ----------
    logger : dict
        dict of warnings and messages
    structure : dict
        dictionary with the computed Dn50, rock class, and average Dn50
        of the rock class for each layer and the toe. This dictionary
        includes all variants, use :py:meth:`get_variant` to get the
        parameters of one specific variant. Alternatively,
        :py:meth:`print_variant` can be used to print the details of
        one, multiple or all variants.
    alpha : float
        slope of the structure in radians
    id : int
        unique id of the breakwater
    variantIDs : list
        list with the IDs of the variants generated for this rubble
        mound breakwater.
    Rc : float
        the crest freeboard of the structure [m]
    width_scour : float
        the required length of the scour protection [m]
    """

    def __init__(
        self,
        slope,
        slope_foreshore,
        B,
        rho_w,
        LimitState,
        ArmourUnit,
        Grading,
        core_material,
        safety=1,
        slope_toe=(2, 3),
        structure_type="revetment",
        B_toe=None,
        beta=0,
        layers=1,
        layers_underlayer=2,
        filter_rule=None,
        Soil=None,
        phi=40,
        id=None,
        **kwargs,
    ):
        """See help(ConcreteRubbleMound) for more info"""
        # set logger and structure
        self.logger = {"INFO": [], "WARNING": []}
        self.structure = {}

        # compute angles
        self.alpha = np.arctan(slope[0] / slope[1])
        slope_foreshore = np.arctan(slope_foreshore[0] / slope_foreshore[1])

        # compute relative buoyant density
        self.id = id
        self.variantIDs = ["a", "b", "c", "d"]

        self.rho = ArmourUnit.rho

        # compute relative buoyant density
        delta_armour = (ArmourUnit.rho - rho_w) / rho_w

        # convert single LimitState to list if needed
        if isinstance(LimitState, list):
            LimitStates = LimitState
        else:
            LimitStates = [LimitState]

        # determine the type of armour layer
        armour_layer = ArmourUnit.name
        if armour_layer in ["Xbloc", "XblocPlus"]:
            # the formula for Xbloc and XblocPlus is based on the
            # Hudson formula with a fixed slope of 3:4
            angle = np.arctan(3 / 4)
            filter_rule = armour_layer
        else:
            angle = self.alpha

        # raise error if no filter rule is specified
        if filter_rule is None:
            supported_rules = ", ".join(substructure._supported_armour_layers())
            raise NotSupportedError(
                (
                    f"Filter rule for {armour_layer} is not implemented, set "
                    f"filter rule to use with filter_rule to {supported_rules}"
                )
            )

        # set temporary values to check changes
        dn, state = 0, 0

        for i, LimitState in enumerate(LimitStates):
            # unpack Hydraulic Conditions
            Hs = LimitState.get_Hs(definition="H13")

            # design armour layer
            dn_temp = hudson(H=Hs, Kd=ArmourUnit.kd, Delta=delta_armour, alpha=angle)

            # check if computed Dn50 of current LimitState is larger
            # than current normative Dn50
            if dn_temp > dn:
                # if larger the normative Dn50 must be changed
                dn = dn_temp
                state = i

        V_required = ArmourUnit.get_class(dn)

        self.structure["armour"] = {
            "computed Dn50": dn,
            "class": V_required,
            "class Dn50": V_required ** (1 / 3),
            "state": state,
            "layers": layers,
            'Mass': ArmourUnit.units[V_required]['M']
        }

        # design underlayer, filter layer and crest height
        super().__init__(
            Dn50=V_required ** (1 / 3),
            core_material=core_material,
            rho=ArmourUnit.rho,
            rho_w=rho_w,
            armour_layer=armour_layer,
            layers=layers,
            LimitStates=LimitStates,
            Grading=Grading,
            safety=safety,
            layers_underlayer=layers_underlayer,
            slope_toe=slope_toe,
            structure_type=structure_type,
            B_toe=B_toe,
            slope=slope,
            B=B,
            beta=beta,
            id=id,
            filter_rule=filter_rule,
            Soil=Soil,
            phi=phi,
            **kwargs,
        )

        # check if a correction factor must be applied for Xbloc or
        # XblocPlus (DMC design guidelines for Xbloc and XblocPlus)
        correction = getattr(ArmourUnit, "correction_factor", None)
        if callable(correction):
            # make dict to pass to method correction_factor
            params = {
                "Hs": LimitStates[state].get_Hs("H13"),
                "h": LimitStates[state].h,
                "Rc": self.Rc,
                "occurrence_hs": False,
                "slope": self.alpha,
                "slope_foreshore": slope_foreshore,
                "permeability": "permeable",
                "Dn": dn,
                "layers": layers,
                "B": B,
                "beta": beta,
            }

            correction_factor = ArmourUnit.correction_factor(
                logger=self.logger, **params
            )

        else:
            self.logger["INFO"].append(
                "Given ArmourUnit did not have a method to compute a "
                "correction factor. See documentation how to define a method "
                "to compute a correction factor for a custom armour unit"
            )

            correction_factor = 1.0

        # check if a correction factor must be applied
        if correction_factor != 1:
            new_dn = correction_factor ** (1 / 3) * dn

            # replave values in structure with new ones for armour
            new_class = ArmourUnit.get_class(new_dn)
            self.structure["armour"] = {
                "computed Dn50": new_dn,
                "class": new_class,
                "class Dn50": new_class ** (1 / 3),
                "state": state,
                "layers": layers,
                'Mass': ArmourUnit.units[V_required]['M']
            }

            # check if the class of the units has changed
            if new_class != V_required:
                self.logger["INFO"].append(
                    "correction factor resulted in new class, so design will "
                    "be changed"
                )

                # because of the new class the design for the underlayer
                # and filter layer must be changed as well
                super().__init__(
                    Dn50=new_class ** (1 / 3),
                    core_material= core_material,
                    rho=ArmourUnit.rho,
                    rho_w=rho_w,
                    armour_layer=armour_layer,
                    layers=layers,
                    LimitStates=LimitStates,
                    Grading=Grading,
                    safety=safety,
                    layers_underlayer=layers_underlayer,
                    slope_toe=slope_toe,
                    structure_type=structure_type,
                    B_toe=B_toe,
                    slope=slope,
                    B=B,
                    beta=beta,
                    id=id,
                    filter_rule=filter_rule,
                    Soil=Soil,
                    phi=phi,
                    **kwargs,
                )

            else:
                # same class so design does not have to be changed
                self.logger["INFO"].append(
                    "correction factor did not result in a new class, so "
                    "design will not be changed"
                )

    def __str__(self):
        return (
            f"id.{self.id}: breakwater with armour units as armour layer, "
            f"and variants: {self.variantIDs}"
        )

    def cost(
        self,
        *variants,
        equipment=None,
        core_price=None,
        unit_price=None,
        transport_cost= None,
        output="variant",
    ):
        """Compute the cost per meter for each variant for the materials or the CO2 footprint

        Method to compute the cost of each generated variant, the cost
        is computed per meter. The cost of the rocks in the substructure
        must be specified in the RockGrading. If transport cost are not
        included in the price of rocks or core_price it can be given
        with the argument transport_cost.

        .. note::
           The transport_cost are not added to the price of the armour
           layer. The assumption has been made that the cost of
           producing and transporting the armour units is included in
           the unit_price.

        Parameters
        ----------
        *variants : str
            IDs of the variants to plot, see :py:attr:`variantIDs` for
            a list of all generated variants. If 'all' is in the
            arguments, all variants will be plotted.
        type: {'Material, 'CO2'}
            Indicate whether the costs are calculated for the material or the CO2
        unit_price : float
            the cost of an armour unit per m³
        transport_cost : float, optional, default: None
            the cost to transport a m³ of rock from the quarry to the
            project location
        output : {variant, layer, average}
            format of the output dict, variant returns the total cost
            of each variant, layer the cost of each layer for each
            variant and average returns the average cost.

        Returns
        -------
        dict
            the cost

        Raises
        ------
        RockGradingError
            if no pricing is included in the given RockGrading
        """
        # compute the cost of the concept

        cost = self.total_cost(
            *variants,
            equipment=equipment,
            core_price=core_price,
            unit_price=unit_price,
            transport_cost=transport_cost,
            output=output,
        )

        return cost
