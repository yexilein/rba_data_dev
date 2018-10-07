"""Module generating subtilis model."""

# python 2 compatibility
from __future__ import absolute_import, division, print_function

# imports
import sys
from os import path
import math

sys.path = [path.join(sys.path[0], '../../rba')] + sys.path
import rba  # noqa
from rba.utils.inject_efficiencies import EfficiencyInjecter  # noqa


def main():
    builder = rba.ModelBuilder('params.in')
    subtilis = builder.build_model()
    subtilis.set_medium('data/curated_medium.tsv')
    EfficiencyInjecter(subtilis, 'data/catalytic_activity_medium_2.csv')
    apply_old_stoichiometries(subtilis.enzymes)
    add_flagella_constraint(subtilis)
    subtilis.write()


def add_flagella_constraint(subtilis):
    subtilis.targets.target_groups.append(flagella_activation())
    for fn in flagella_activation_functions():
        subtilis.parameters.functions.append(fn)
    subtilis.parameters.aggregates.append(flagella_activation_aggregate())


def flagella_activation():
    target_group = rba.xml.TargetGroup('flagella_activation')
    target = rba.xml.TargetReaction('Th')
    target.value = 'flagella_proton_flux'
    target_group.reaction_fluxes.append(target)
    return target_group


def flagella_activation_functions():
    return [
        rba.xml.Function('flagella_speed', 'constant', {'CONSTANT': 5.81}),
        rba.xml.Function('flagella_h_consumption', 'constant',
                         {'CONSTANT': 0.9415}),
        rba.xml.Function('number_flagella', 'linear',
                         {'LINEAR_COEF': 4.5197, 'LINEAR_CONSTANT': 3.7991,
                          'X_MIN': 0.25, 'X_MAX': 1.6,
                          'Y_MIN': float('-Inf'), 'Y_MAX': float('Inf')})
        ]


def flagella_activation_aggregate():
    aggregate = rba.xml.Aggregate('flagella_proton_flux', 'multiplication')
    aggregate.function_references.append(
        rba.xml.FunctionReference('number_flagella')
    )
    aggregate.function_references.append(
        rba.xml.FunctionReference('flagella_speed')
    )
    aggregate.function_references.append(
        rba.xml.FunctionReference('flagella_h_consumption')
    )
    return aggregate


def add_zero_cost_flags(enzymes):
    with open('data/zero_cost.csv', 'r') as f:
        zero_cost = []
        for line in f:
            zero_cost += line.rstrip('\n').split('\t')
    for enzyme in enzymes.enzymes:
        id_ = old_name(enzyme.id)
        if id_ in zero_cost:
            enzyme.zero_cost = True
            zero_cost.remove(id_)


def read_activities(f, medium):
    fn_type = None
    result = {}
    for line in f:
        token = line.rstrip('\n').split('\t')
        if token[0] == 'function':
            # format 'function id type'
            if token[1] == medium:
                fn_type = token[2]
        else:
            # format 'reaction_id fn_id param1_id param1_value ... paramn_value
            if token[1] == medium:
                reaction = token[0]
                params = iter(token[2:])
                result[reaction] = {
                    id_: value for id_, value in zip(params, params)
                    }
    if not fn_type:
        raise UserWarning('Could not retrieve medium {}.'.format(medium))
    return result


def old_name(new_name):
    if new_name == 'R_maintenance_atp':
        return 'Eatpm'
    else:
        return new_name.rsplit('_', 1)[0]


def apply_old_stoichiometries(enzymes):
    with open('data/stoichiometry.csv', 'r') as f:
        data = {}
        for line in f:
            [species, sto] = line.rstrip('\n').split('\t')
            data[species] = float(sto)
    for enzyme in enzymes.enzymes:
        for sr in enzyme.machinery_composition.reactants:
            try:
                sr.stoichiometry = data[sr.species]
            except KeyError:
                pass

if __name__ == "__main__":
    main()
