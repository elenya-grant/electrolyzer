"""
This module defines the Hydrogen Electrolyzer control code.
"""
import numpy as np

from electrolyzer import Electrolyzer


class ElectrolyzerSupervisor:
    def __init__(self, electrolyzer_dict, control_type, dt=1):

        # H2 farm parameters
        self.n_stacks = electrolyzer_dict["n_stacks"]
        self.n_cells = electrolyzer_dict["n_cells"]
        self.cell_area = electrolyzer_dict["cell_area"]

        # TODO query electrolyzer for this
        self.stack_rating_kW = electrolyzer_dict["stack_rating_kW"]
        self.stack_rating = self.stack_rating_kW * 1e3

        # TODO query electrolyzer for this
        self.stack_min_power = (10 / 100) * self.stack_rating

        # TODO remove
        self.stack_input_voltage = electrolyzer_dict["stack_input_voltage"]
        self.temperature = electrolyzer_dict["temperature"]
        self.dt = dt

        # Controller storage variables

        # array of stack activation status 0 for inactive, 1 for active
        self.active = np.zeros(self.n_stacks)

        # array of stack waiting status 0 for active or inactive, 1 for waiting
        self.waiting = np.zeros(self.n_stacks)

        # only for sequential controller TODO: find sneakier place to initialize this
        self.active_constant = np.zeros(self.n_stacks)
        self.variable_stack = 0  # again, only for sequential controller
        self.stack_rotation = []
        self.stacks_on = 0
        self.stacks_waiting = 0
        self.stacks_off = []
        self.stacks_waiting_vec = np.zeros((self.n_stacks))
        self.deg_state = np.zeros(self.n_stacks)

        self.stacks = self.create_electrolyzer_stacks()  # initialize stack objects

        self.control_type = control_type
        """
        --- Current control_type Options ---

        Rotation-based electrolyzer action schemes:
            'power sharing rotation': power sharing, rotation
            'sequential rotation': sequentially turn on electrolzyers, rotate
                electrolyzer roles based on set schedule (i.e. variable electrolyzer,
                etc.)

        Degredation-based electrolyzer action schemes:
            'even split eager deg': power sharing, eager to turn on electrolyzers
            'even split hesitant deg': power sharing
            'sequential even wear deg': sequentially turn on electrolzyers, distribute
                wear evenly
            'sequential single wear deg': sequentially turn on electrolyzers, put all
                degradation on single electrolyzer
            'baseline deg': sequtntially turn on and off electrolyzers but only when you
                have to
        """
        if "sequential" in self.control_type:
            # TODO: current filter width hardcoded at 5 min, make an input
            self.filter_width = round(1200 / self.dt)

            # TODO: decide how to initialize past_power
            self.past_power = [0]

        # delete all these eventually they are just for troubleshooting and plotting
        self.P_indv_store = []
        self.active_store = []
        self.deg_state_store = []
        self.waiting_store = []
        self.active_actual_store = []
        self.H2_store = []
        self.unused_power = []

    def create_electrolyzer_stacks(self):
        # initialize electrolyzer objects
        stacks = []
        for i in range(self.n_stacks):
            stacks.append(
                Electrolyzer(self.n_cells, self.cell_area, self.temperature, dt=self.dt)
            )
            self.stack_rotation.append(i)
            print(
                "electrolyzer stack ",
                i + 1,
                "out of ",
                self.n_stacks,
                "has been initialized",
            )
        return stacks

    def control(self, power_in):
        """
        Inputs:
            power_in: power (W) to be consumed by the H2 farm every time step
        Returns:
            H2_mass_out: mass of h2 (kg) produced during each time step
            H2_mass_flow_rate: mfr of h2 (kg/s) during that time step
            power_left: power error (W) between what the stacks were
                supposed to consume and what they actually consumed
            curtailed_wind: power error (W) between available power_in and
                what the stacks are commanded to consume
        """

        # calculate stack power distribution
        if self.control_type == "power sharing rotation":
            stack_power, curtailed_wind = self.power_sharing_rotation(power_in)
        elif self.control_type == "sequential rotation":
            stack_power, curtailed_wind = self.sequential_rotation(power_in)
        elif self.control_type == "even split eager deg":
            stack_power = self.distribute_power_equal_eager(power_in)
            curtailed_wind = 0
        elif self.control_type == "even split hesitant deg":
            stack_power = self.distribute_power_equal_hesitant(power_in)
            curtailed_wind = 0
        elif self.control_type == "sequential even wear deg":
            stack_power = self.distribute_power_sequential_even_wear(power_in)
            curtailed_wind = 0
        elif self.control_type == "sequential single wear deg":
            stack_power = self.distribute_power_sequential_single_wear(power_in)
            curtailed_wind = 0
        elif self.control_type == "baseline deg":
            stack_power = self.baseline_controller(power_in)
            curtailed_wind = 0

        # Query stacks for their status and turn them on or off as needed
        on_or_waiting = np.zeros(self.n_stacks)

        if "deg" in self.control_type:
            for i in range(self.n_stacks):
                if self.stacks[i].stack_on or self.stacks[i].stack_waiting:
                    on_or_waiting[i] = 1

            # which stacks the controller thinks are on and which are actually on
            mismatch = self.active - on_or_waiting

            for i in range(len(mismatch)):
                # this means the controller wants an electrolyzer on and that
                # electrolyzer isnt on or waiting
                if mismatch[i] == 1:
                    self.stacks[i].turn_stack_on()
                elif mismatch[i] == 0:
                    pass
                # this means the controller wants an electrolyzer off and
                # the electrolyzer is on
                elif mismatch[i] == -1:
                    self.stacks[i].turn_stack_off()

            for i in range(self.n_stacks):
                if self.stacks[i].stack_waiting:
                    self.waiting[i] = 1
                else:
                    self.waiting[i] = 0

        active_actual = np.zeros(self.n_stacks)
        for i in range(self.n_stacks):
            if self.stacks[i].stack_on:
                active_actual[i] = 1

        power_left = 0
        H2_mass_out = 0
        self.stacks_waiting = 0
        self.stacks_on = 0
        H2_mass_flow_rate = np.zeros((self.n_stacks))

        # simulate 1 time step for each stack
        for i in range(self.n_stacks):
            H2_mfr, H2_mass_i, power_left_i = self.stacks[i].run(stack_power[i])

            self.deg_state[i] = self.stacks[i].V_degradation

            # Update stack status
            if self.stacks[i].stack_on:
                self.stacks_on += 1
                self.active[i] = 1
                on_or_waiting[i] = 1
            if self.stacks[i].stack_waiting:
                self.waiting[i] = 1
                self.stacks_waiting_vec[i] = 1
                on_or_waiting[i] = 1
            else:
                self.waiting[i] = 0
                self.stacks_waiting_vec[i] = 0

            H2_mass_flow_rate[i] = H2_mfr
            H2_mass_out += H2_mass_i
            power_left += power_left_i

        curtailed_wind = max(0, power_in - (np.dot(on_or_waiting, stack_power)))

        self.P_indv_store.append(stack_power)
        self.active_store.append(np.copy(self.active))
        self.waiting_store.append(np.copy(self.waiting))
        self.deg_state_store.append(np.copy(self.deg_state))
        self.active_actual_store.append(np.copy(active_actual))
        self.H2_store.append(np.copy(H2_mass_flow_rate))
        self.unused_power.append(np.copy(power_left))

        return H2_mass_out, H2_mass_flow_rate, power_left, curtailed_wind

    def power_sharing_rotation(self, power_in):
        # Control strategy that shares power between all electrolyzers equally
        if sum(self.active + self.waiting) == 0:
            P_indv = np.ones(1) * power_in / self.n_stacks
        else:
            P_indv = (
                np.ones(1) * power_in / sum(self.active + self.waiting)
            )  # divide the power evenely amongst electrolyzers
        P_indv_kW = P_indv / 1000

        stacks_supported = min(power_in // (self.stack_rating / 2), self.n_stacks)
        diff = int(stacks_supported - sum(self.active + self.waiting))

        # Power sharing control #
        #########################
        if diff > 0:
            # elif P_indv_kW > (0.8 * self.stack_rating_kW):
            if diff > 1 or P_indv_kW > (0.8 * self.stack_rating_kW):
                if sum(self.waiting) == 0 and sum(self.active) != self.n_stacks:
                    for i in range(0, diff):
                        ij = 0 + i
                        while self.active[self.stack_rotation[ij]] > 0:
                            ij += 1
                        self.turn_on_stack(self.stack_rotation[ij])

        if diff < 0:
            if P_indv_kW < (0.2 * self.stack_rating_kW):
                if sum(self.active) > 0:
                    self.turn_off_stack(self.stack_rotation[0])
                    self.stack_rotation = self.stack_rotation[1:] + [
                        self.stack_rotation[0]
                    ]

        new_stack_power = (
            np.ones((self.n_stacks)) * power_in / sum(self.active + self.waiting)
        )
        if (new_stack_power[0] / 1000) > (self.stack_rating_kW):
            curtailed_wind = (
                new_stack_power[0] - (self.stack_rating_kW * 1000)
            ) * self.stacks_on
            new_stack_power = np.ones((self.n_stacks)) * (self.stack_rating_kW * 1000)
        else:
            curtailed_wind = 0

        return new_stack_power, curtailed_wind

    def sequential_rotation(self, power_in):
        # Control strategy that fills up the electrolyzers sequentially

        P_indv = np.ones((self.n_stacks))

        p_in_kw = power_in / 1000

        n_full = p_in_kw // self.stack_rating_kW
        left_over_power = p_in_kw % self.stack_rating_kW
        stack_difference = int(n_full - sum(self.active + self.waiting))
        elec_var = self.stack_rotation[0]
        # other_elecs = [x for x in self.active if x not in self.stacks_off]

        # calculate the slope of power_in
        # (1) update past_power with current input
        if len(self.past_power) < self.filter_width:
            temp = np.zeros((len(self.past_power) + 1))
            temp[0:-1] = self.past_power[:]
            temp[-1] = np.copy(power_in)
            self.past_power = np.copy(temp)
        else:
            temp = np.zeros((len(self.past_power)))
            temp[0:-1] = self.past_power[1:]
            temp[-1] = np.copy(power_in)
            self.past_power = np.copy(temp)
        # (2) apply filter to find slope
        slope = power_in - np.mean(self.past_power)
        slope = np.mean(
            (np.mean(self.past_power[1:]) - np.mean(self.past_power[0:-1])) / self.dt
        )
        print(slope)

        if stack_difference >= 0:
            # P_indv = P_indv * self.stack_rating_kW * 1000
            P_indv = P_indv * 0
            # for i in self.active:
            #     if i > 0:
            P_indv[self.active > 0] = self.stack_rating_kW * 1000
            curtailed_wind = (stack_difference * self.stack_rating_kW) + left_over_power
            if (
                sum(self.waiting) == 0
                and sum(self.active) != self.n_stacks
                and left_over_power > (0.15 * self.stack_rating_kW)
            ):
                for i in range(0, stack_difference):
                    ij = 0 + i
                    while self.active[self.stack_rotation[ij]] > 0:
                        ij += 1
                    self.turn_on_stack(self.stack_rotation[ij])
        if stack_difference < 0:
            curtailed_wind = 0
            P_indv = P_indv * 0

            # if n_full > 0
            # for i in self.active:
            #     if i > 0:
            #         P_indv[i] = self.stack_rating_kW * 1000
            P_indv[self.active > 0] = self.stack_rating_kW * 1000

            # for i in self.stacks_off:
            #     P_indv[i] = 0
            if stack_difference > -2:
                if sum(self.waiting) == 0 and sum(self.active) != self.n_stacks:
                    ij = 0
                    while self.active[self.stack_rotation[ij]] > 0:
                        ij += 1
                    on_stack = self.stack_rotation[ij]
                    self.turn_on_stack(on_stack)
                    P_indv[on_stack] = self.stack_rating_kW * 1000

            elif (
                (
                    left_over_power < (0.1 * self.stack_rating_kW)
                    and sum(self.waiting) == 0
                )
            ) or (stack_difference < -1 and sum(self.waiting) == 0):
                if sum(self.active) > 0 and slope < 0:
                    self.turn_off_stack(self.stack_rotation[0])
                    self.stack_rotation = self.stack_rotation[1:] + [
                        self.stack_rotation[0]
                    ]
                    P_indv[elec_var] = 0
                    curtailed_wind = left_over_power
            elif (
                left_over_power < (0.1 * self.stack_rating_kW) and sum(self.waiting) > 0
            ):
                P_indv[self.waiting > 0] = (
                    (left_over_power + self.stack_rating_kW) * 1000 / 2
                )
                P_indv[elec_var] = (left_over_power + self.stack_rating_kW) * 1000 / 2
            # TODO : Find a way to turn on electrolyzers ahead of time based on wind
            # power signal slope
            elif left_over_power > (0.8 * self.stack_rating_kW) and slope > 0:
                # if or stack_difference
                if sum(self.waiting) == 0 and sum(self.active) != self.n_stacks:
                    ij = 0
                    while self.active[self.stack_rotation[ij]] > 0:
                        ij += 1
                    self.turn_on_stack(self.stack_rotation[ij])
                    # self.turn_on_stack()
                if sum(self.waiting) > 0:
                    P_indv[self.waiting > 0] = left_over_power * 1000 / 2
                    P_indv[elec_var] = left_over_power * 1000 / 2
                else:
                    P_indv[elec_var] = left_over_power * 1000
            else:
                if sum(self.waiting) > 0:
                    # if
                    P_indv[self.waiting > 0] = self.stack_rating_kW * 1000
                    P_indv[elec_var] = left_over_power * 1000
                else:
                    P_indv[elec_var] = left_over_power * 1000
        return P_indv, curtailed_wind * 1000

    def distribute_power_equal_eager(self, power_in):
        n_active = min([self.n_stacks, int(np.floor(power_in / self.stack_min_power))])
        if n_active > 0:
            # calculate curtailed wind here
            P_i = min([power_in / n_active, self.stack_rating])
            # curtailed_wind = max(0, (power_in/n_active) - self.stack_rating)
        else:
            P_i = 0

        if n_active == sum(self.active):
            pass  # do not need to turn on or off
        elif n_active > sum(self.active):
            diff = int(n_active - sum(self.active))
            self.active += self.get_healthiest_inactive(
                self.active, self.deg_state, diff
            )
        elif n_active < sum(self.active):
            diff = int(sum(self.active) - n_active)
            self.active *= self.get_illest_active(self.active, self.deg_state, diff)

        P_indv = P_i * self.active

        return P_indv

    def distribute_power_equal_hesitant(self, power_in):
        # dont turn on another electrolyzer until the other electrolyzers are all at
        # rated. turn off when all electrolyzers are below min power from previous step
        n_active = int(sum(self.active))

        # need to turn on another one
        if power_in > (self.stack_rating * n_active + self.stack_min_power):
            # number of stacks that need to be turned on
            diff = np.ceil(
                (power_in - n_active * self.stack_rating) / self.stack_rating
            )
            diff = min([diff, self.n_stacks - n_active])
            n_active += diff
        # not enough power to run all electrolyzers at minimum so we have to turn
        # some off
        elif power_in < self.stack_min_power * n_active:
            # want the number of elecs where np.ceil(P*min_rated)
            diff = n_active - np.floor(power_in / self.stack_min_power)
            n_active -= diff

        if n_active > 0:
            P_i = min([power_in / n_active, self.stack_rating])
        else:
            P_i = 0

        if n_active == sum(self.active):
            pass  # do not need to turn on or off
        elif n_active > sum(self.active):
            # need to turn on this many electrolzyers pick the healthiest elecs as
            # the next ones to turn on but only pick the healthiest from the ones
            # that are turned off
            diff = int(n_active - sum(self.active))
            self.active += self.get_healthiest_inactive(
                self.active, self.deg_state, diff
            )

        elif n_active < sum(self.active):
            diff = int(sum(self.active) - n_active)  # need to turn off this many
            self.active *= self.get_illest_active(self.active, self.deg_state, diff)

        P_indv = P_i * self.active

        return P_indv

    def distribute_power_sequential_even_wear(self, power_in):
        P_indv = np.zeros(self.n_stacks)

        n_active = np.min(
            [
                self.n_stacks,
                np.ceil((power_in - self.stack_min_power) / self.stack_rating),
            ]
        )

        diff = int(n_active - sum(self.active))

        if n_active > sum(self.active):
            stacks_to_turn_on = self.get_healthiest_inactive(
                self.active, self.deg_state, diff
            )
            # for robustitude, this index should be [0][0] but it throws error so
            # come back to this
            self.variable_stack = np.nonzero(stacks_to_turn_on)[0][0]
            self.active += stacks_to_turn_on
            self.active_constant = np.copy(self.active)
            self.active_constant[self.variable_stack] = 0
        elif n_active < sum(self.active):
            self.active[self.variable_stack] = 0

            # need this to be get illest constant
            stacks_to_turn_off = self.get_illest_active(
                self.active_constant, self.deg_state, np.abs(diff)
            )

            # these are the illest stacks from constant - the illest of these should
            # turn into variable
            self.variable_stack = np.nonzero(stacks_to_turn_off - 1)[0]

            # this only works if there is just one stack being turned off at at time
            self.active_constant[self.variable_stack] = 0

        variable_stack_P = np.min(
            [
                self.stack_rating,
                power_in - (sum(self.active_constant)) * self.stack_rating,
            ]
        )
        if variable_stack_P < self.stack_min_power:
            variable_stack_P = 0

        P_indv = self.stack_rating * self.active_constant
        P_indv[self.variable_stack] = variable_stack_P

        return P_indv

    def distribute_power_sequential_single_wear(self, power_in):
        # if we are trying a lot of different control strategies, maybe it would be
        # better to keep them in a seperate file and call them as a module
        P_indv = np.zeros(self.n_stacks)

        n_active = np.min(
            [
                self.n_stacks,
                np.ceil((power_in - self.stack_min_power) / self.stack_rating),
            ]
        )

        diff = int(n_active - sum(self.active))

        if n_active > sum(self.active):
            stacks_to_turn_on = self.get_healthiest_inactive(
                self.active, self.deg_state, diff
            )

            # for robustitude, this index should be [0][0] but it throws error so
            # come back to this
            # self.variable_stack = np.nonzero(stacks_to_turn_on)[0][0]
            self.variable_stack = 0
            self.active += stacks_to_turn_on
            self.active_constant = np.copy(self.active)
            self.active_constant[self.variable_stack] = 0
        elif n_active < sum(self.active):
            # self.active[self.variable_stack] = 0

            # need this to be get illest constant
            stacks_to_turn_off = self.get_illest_active(
                self.active_constant, self.deg_state, np.abs(diff)
            )

            # these are the illest stacks from constant - the illest of
            # these should turn into variabl
            # self.variable_stack = np.nonzero(stacks_to_turn_off-1)[0]

            # this only works if there is just one stack being turned off at at time
            # self.active_constant[self.variable_stack] = 0
            self.active_constant *= stacks_to_turn_off
            self.active *= stacks_to_turn_off

        variable_stack_P = np.min(
            [
                self.stack_rating,
                power_in - (sum(self.active_constant)) * self.stack_rating,
            ]
        )
        if variable_stack_P < self.stack_min_power:
            variable_stack_P = 0

        P_indv = self.stack_rating * self.active_constant
        P_indv[self.variable_stack] = variable_stack_P

        return P_indv

    def baseline_controller(self, power_in):
        """
        Hesitant to turn on, hesitant to turn off
        """

        p_avail = power_in

        # turn some on
        if (power_in > np.sum(self.active) * self.stack_rating) & (
            power_in > self.stack_min_power
        ):
            turn_on_ind = int(min([sum(self.active), self.n_stacks - 1]))
            self.stacks[turn_on_ind].turn_stack_on()
            self.active[turn_on_ind] = 1

        # turn some off
        elif power_in < np.sum(self.active) * self.stack_min_power:
            turn_off_ind = int(sum(self.active) - 1)
            self.stacks[turn_off_ind].turn_stack_off()
            self.active[turn_off_ind] = 0

        P_indv = self.stack_min_power * self.active
        p_avail -= sum(P_indv)

        for i in range(self.n_stacks):

            if p_avail >= (self.stack_rating - self.stack_min_power):
                P_indv[i] += self.stack_rating - self.stack_min_power
                p_avail -= self.stack_rating - self.stack_min_power
            elif p_avail < (self.stack_rating - self.stack_min_power):
                P_indv[i] += p_avail
                p_avail = 0

        return P_indv

    def turn_off_stack(self, off_stack):
        # print('turn off stack')

        # off_stack = self.stack_rotation[0]
        # self.stacks_on.remove(off_stack)
        self.active[off_stack] = 0
        self.waiting[off_stack] = 0
        # self.stacks_off.append(off_stack)
        self.stacks[off_stack].turn_stack_off()
        # self.stacks_on = self.stacks_on - 1

    def turn_on_stack(self, on_stack):
        # print('turn on stack')

        # on_stack = self.stacks_off[0]
        # self.stacks_off.remove(on_stack)
        self.stacks[on_stack].turn_stack_on()
        self.waiting[on_stack] = 1
        # print(on_stack)

    def get_healthiest_inactive(self, active, deg_state, n_activate):
        # TODO remove the arrays passed to this method since they are stored in
        # self already
        inact = np.nonzero(active - 1)[0]
        ds = deg_state[inact]
        deg_inds = np.argsort(ds)
        temp = np.zeros_like(active)
        temp[inact[deg_inds[0:n_activate]]] = 1
        return temp  # active + temp will turn on the stacks at 1s

    def get_illest_active(self, active, deg_state, n_deactivate):
        act = np.nonzero(active)[0]
        ds = deg_state[act]
        deg_inds = np.flip(np.argsort(ds))
        temp = np.ones_like(active)
        temp[act[deg_inds[0:n_deactivate]]] = 0

        # active * temp will turn off the stacks at 0s
        return temp

    # want to be able to call this from time to time TODO
    # def calculater_LCOH_from_current_state()
    # roughly, find the slope of degradation over hydrogen then predict the lifetime
