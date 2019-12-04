# In this file we will implement FBA using cobrapy - 
# https://cobrapy.readthedocs.io/en/stable/building_model.html

from __future__ import print_function

from cobra import Model, Reaction, Metabolite




def create_metabolite(met_id, formula, name, compartment):
    new_met = Metabolite(met_id, formula, name, compartment)
    return new_met

# metabolites contains both the metabolite information and it's coefficient
def create_reaction(rxn_id, name, subsystem, lower_bound, upper_bound, metabolites_d3):
    reaction = Reaction(rxn_id)
    reaction.name = name
    reaction.subsystem = subsystem
    reaction.lower_bound = lower_bound
    reaction.upper_bound = upper_bound
    for metabolite_d2 in metabolites_d3:
        metabolite_data = metabolite_d2[0]
        metabolite_coeff = metabolite_d2[1]
        new_met = create_metabolite(metabolite_data)
        reaction.add_metabolites({ new_met: metabolite_coeff })

    return reaction


def create_model_from_rxns(name, rxns_d5, objective_rxn_id):
    model = Model(name)
    for rxn_d4 in rxns_d5:
        new_rxn = create_rxn(rxn_d4)
        model.add_reactions([new_rxn])
    model.objective = objective_rxn_id
    return model

def run_fba(model_name, rxns_d5, objective_rxn_id, media):
    model = create_model_from_rxns(model_name, rxns_d5, objective_rxn_id)
    solution = model.optimize()

    print(solution.objective_value) # The objective value
    print(solution.status) # The status from the linear programming solver
    print(solution.fluxes) # A pandas series with flux indexed by reaction identifier.
    print(solution.shadow_prices) # A pandas series with shadow price indexed by the metabolite id
   

    



def online_example():
    model = Model('example_model')
    reaction = Reaction('30AS140')
    reaction.name = '3 oxoacyl acyl carrier protein synthase n C140 '
    reaction.subsystem = 'Cell Envelope Biosynthesis'
    reaction.lower_bound = 0.  # This is the default
    reaction.upper_bound = 1000.  # This is the default


    #creating the metabolites
    ACP_c = Metabolite(
    'ACP_c',
    formula='C11H21N2O7PRS',
    name='acyl-carrier-protein',
    compartment='c')

    omrsACP_c = Metabolite(
    '3omrsACP_c',
    formula='C25H45N2O9PRS',
    name='3-Oxotetradecanoyl-acyl-carrier-protein',
    compartment='c')
    co2_c = Metabolite('co2_c', formula='CO2', name='CO2', compartment='c')

    malACP_c = Metabolite(
    'malACP_c',
    formula='C14H22N2O10PRS',
    name='Malonyl-acyl-carrier-protein',
    compartment='c')
    h_c = Metabolite('h_c', formula='H', name='H', compartment='c')

    ddcaACP_c = Metabolite(
    'ddcaACP_c',
    formula='C23H43N2O8PRS',
    name='Dodecanoyl-ACP-n-C120ACP',
    compartment='c')


    #Adding stoichiometric coefficients to reaction.
    reaction.add_metabolites({
    malACP_c: -1.0,
    h_c: -1.0,
    ddcaACP_c: -1.0,
    co2_c: 1.0,
    ACP_c: 1.0,
    omrsACP_c: 1.0
    })

    #Gene Reaction Rule: A Boolean representation of the gene
    # Requirements for this reaction to be active.
    reaction.gene_reaction_rule = '( STM2378 or STM1197 )'
    

    #Adding reaction to the model:
    model.add_reactions([reaction])

    # Creating an objective for the model:
    model.objective = '30AS140'


    #Solution object has several attributes
    solution = model.optimize()

    print(solution.objective_value) # The objective value
    print(solution.status) # The status from the linear programming solver
    print(solution.fluxes) # A pandas series with flux indexed by reaction identifier.
    print(solution.shadow_prices) # A pandas series with shadow price indexed by the metabolite id



online_example()
