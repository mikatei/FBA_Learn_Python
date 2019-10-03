# In this file we make use of optlang
# More info here: https://optlang.readthedocs.io/en/latest/

from optlang import Model, Variable, Constraint, Objective

# You can declare the symbolic variables here with upper and lower bounds:

'''

x1 = Variable('x1', lb=0, ub = 100)

'''

#S is the stoichiomatrix as passed in by numpy
# objective function is the last variable (v)
# upperbound needs to be an int
# Constraints is as long as the amount of compounds
# flux_bounds is as long as the amount of reactions. It is a d2_list
# flux_bounds = [[lower bounds],[upper bounds]]
def stoichiomatrix_solution(S, flux_bounds, objective_index, objective_direction):

    #We make a variable 'v-(index)' for each reaction (column) in the matrix:
    variables = make_variables(S, flux_bounds)
    constraints = make_constraints(S, variables)
    obj = make_objective(objective_index, objective_direction, variables)

    model= Model(name='Stoichiomatrix')
    model.objective = obj
    model.add(constraints)
    status = model.optimize()

    return [status, model]


# This function makes the variables
def make_variables(S, flux_bounds):
    variables = []
    row_1 = S[0]
    for i in range(len(row_1)):
        v = Variable('v-' + str(i+1), lb = flux_bounds[0][i], ub = flux_bounds[1][i])
        variables.append(v)
    print(variables)
    return variables


def make_constraints(S, variables):
    #Creating the constraints, one per compound:
    constraints = []
    for row in S:
        constraint_sum = 0
        for i in range(len(row)):
            constraint_sum += row[i]*variables[i]
        c = Constraint(constraint_sum, lb=0, ub =0)
        constraints.append(c)
    return constraints


def make_objective(objective_index, objective_direction, variables):

    #The objective is just to either Maximize or Minimize a Variable.
    obj_var = variables[objective_index]
    print("Objective variable name: " + obj_var.name)
    obj = Objective(variables[objective_index], direction = objective_direction)

    return obj


def model_print(model):
    print("status:", model.status)
    #print("objective variable name: " + model.objective.name)
    print("objective value:", model.objective.value)
    print("----------")
    print(model.variables.items())
    for var_name, var in model.variables.items():
        print(var_name, "=", var.primal)


def make_fluxes(model):
    #fluxes holds the names and their values, then we sort by that and make the fluxes array
    fluxes = []
    for var_name, var in model.variables.items():
        fluxes.append([int(var_name[2:]),var.primal])
    fluxes.sort(key = lambda fluxes: fluxes[0])
    flux_array = []
    for flux in fluxes:
        flux_array.append(flux[1])
    return flux_array
