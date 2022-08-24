#!/usr/bin/env python

# Jack Dewey
# CSC 445 Programming Project

import sys
import argparse
import numpy
import re
min_rounding_val = 0.0000000069


class Iteration():
    '''
    Class is used to track basic and non-basic variables
    '''
    def __init__(self, seq):
        self.index = {var:index for index, var in enumerate(seq)}
        self.lo_vars = list(seq)

    def replacing(self, variable1, variable2):
        i = self.index.pop(variable1)
        self.lo_vars[i] = variable2
        self.index[variable2] = i

    def adding(self, to_add):
        i = len(self.lo_vars)
        self.lo_vars.append(to_add)
        self.index[to_add] = i


def filereader():
    '''
    :return: Returns a list of 3 lists.
                Index 0 is a list of objective variable
                Index 1 is list of constraints
                Index 2 is list of all right hand side constant
    This function reads from a text file and creates appropriate output for simplex solving
    '''
    with open(sys.argv[1], 'r') as f:
        obj = f.readline()
        obj = obj.strip('\n')
        obj = obj.strip()
        obj = obj.split()
        obj = [float(i) for i in obj]
        constraint = []
        constant_set = []
        for line in f:
            new_constraint = line
            input = new_constraint.strip('\n')
            input = input.split()
            input = [float(i) for i in input]
            constant = input.pop()
            constant_set.append(constant)
            constraint.append(input)
    f.close()
    return(obj, constraint, constant_set)


def simplex_start(constraint, constant_set, obj):
    '''
    :param constraint: List of constraints in form of [1,2,3] for x1,x2,x3
    :param constant_set: List of constants in form of [4,5,6], these constants refer to the <= inequalities
    :param obj: List of objective functions, each constant [7,8,9] refer to 5x1, 6x2, 7x3 in the objective function
    :return: Returns a list of optimal values for each constraint, in form of [10, 11, 12] as a solution
                for x1=10, x2=11, x3=12

    This is the core of the program. The simplex function is primarily handled here, with components being handed off
    to subfunctions, like any pivoting and the initial creation of tableau.
    '''

    # Create an auxiliary problem
    auxiliary = do_auxiliary(constraint, constant_set, obj)
    if auxiliary:
        col, row, constraint, constant_set, obj, pivoted = auxiliary

    # If we could not build an auxiliary problem, then it must be infeasible.
    else:
        print("infeasible", end='')
        return None

    # This while loop checks each variable in the objective function. If it sees there are any improvements to be
    # made, then it continues while looping. Otherwise it breaks out and returns optimal solution.
    # Our program applies bland's rule, as it always selects the lowest eligible index
    while any([objective > min_rounding_val for objective in obj]):
        individual_obj, go_next = next((index_of_objfunc, col.lo_vars[index_of_objfunc])
                                        for index_of_objfunc in range(len(obj))
                                        if obj[index_of_objfunc] > min_rounding_val)

        fkmanplswork = [constant_set[i] / constraint[i][individual_obj]
                        if constraint[i][individual_obj] > min_rounding_val else float("inf")
                        for i in range(len(constant_set))]

        min_enum, check_unbound = min(enumerate(fkmanplswork), key = lambda x: x[1])

        # We can determine if it's unbounded if we ever find an infinity
        if check_unbound == float("inf"):
            print("unbounded", end='')
            return None

        # Since it is not unbounded, we can pivot
        else:
            k = row.lo_vars[min_enum]
            col, row, constraint, constant_set, obj, \
            pivoted = pivot(col, row, constraint, constant_set, obj, pivoted, k, go_next)

    # If we reach this point, we must have broken out of the while loop. We can only break out of while loop if
    # each objective variable is greater than zero, thus we have reached optimal solution
    print("optimal")
    return [(constant_set[row.index[index]] if (index in row.lo_vars) else 0) for index in range(len(obj))]

def pivot(rcol, rrow, constraint, constant_set, obj_func, pivoted, k, go_next):
    '''
    :param rcol: Columns
    :param rrow: Rows
    :param constraint: All constraints
    :param constant_set: All constants
    :param obj_func: All objective variables from objfunc
    :param pivoted: Result from previous pivot or from initial build_aux_tableau
    :param k: Index for row of pivot
    :param go_next: Index for col of pivot
    :return: The result of pivoting by k, go_next in the form of rcol, rrow, constraint, constant_set, obj_func, pivoted
                after being pivoted

    This function does the actual pivoting, then returns the result of the pivot
    Is called by both do_auxiliary and simplex_start
    '''
    # The actual pivoting
    row = rrow.index[k]
    col = rcol.index[go_next]
    len_constants = len(constant_set)
    len_obj = len(obj_func)

    # Updating values after dividing
    divide_by = constraint[row][col]
    constant_set[row] /= divide_by
    constraint[row][col] = 1.0
    constraint[row] = [x / divide_by for x in constraint[row]]

    # Determine which coefficients are not in the pivoted row or column
    for i in range(len_constants):
        if i == row: continue
        constant_set[i] -= constraint[i][col] * constant_set[row]
        for j in range(len_obj):
            if j == col: continue
            constraint[i][j] -= constraint[i][col] * constraint[row][j]
        constraint[i][col] /= -divide_by
    pivoted += obj_func[col] * constant_set[row]
    for j in range(len_obj):
        if j == col: continue
        obj_func[j] -= obj_func[col] * constraint[row][j]

    # Update the variables in basis/not basis
    obj_func[col] /= -divide_by
    rcol.replacing(go_next, k)
    rrow.replacing(k, go_next)
    return rcol, rrow, constraint, constant_set, obj_func, pivoted

def do_auxiliary(constraint, constant_set, obj_func):
    '''
    :param constraint: All constraints
    :param constant_set: All constants
    :param obj_func: All objective variables from objfunc
    :return: Returns information for a tableau after creating any needed auxiliary variables
    This function creates and solves auxiliary program, then returning a valid lp to the main solver program.
    This function also calls pivot() to solve the aux problem.
    '''

    len_constants = len(constant_set)
    len_obj = len(obj_func)

    rcol = Iteration(range(len_obj))
    rrow = Iteration(range(len_obj, len_obj + len_constants))
    k = min(range(len(constant_set)), key = lambda x: constant_set[x])
    if constant_set[k] > min_rounding_val:
        return rcol, rrow, constraint, constant_set, obj_func, 0

    # Begin creating the auxiliary
    obj = [x for x in obj_func]
    obj_func = [0] * (len_obj) + [-1]
    constraint = [row + [-1] for row in constraint]
    rcol.adding(len_obj + len_constants)
    rcol, rrow, constraint, constant_set, obj_func, pivoted = \
        pivot(rcol, rrow, constraint, constant_set, obj_func, 0, len_obj + k, len_obj + len_constants)

    # Solve the auxiliary problem using pivots
    # Our program applies bland's rule, as it always selects the lowest eligible index
    while any([x > min_rounding_val for x in obj_func]):
        index_2, e = next((j, rcol.lo_vars[j]) for j in range(len(obj_func)) if obj_func[j] > min_rounding_val)
        temp = [constant_set[index_1] / constraint[index_1][index_2] if constraint[index_1][index_2] > min_rounding_val else float("inf")
                for index_1 in range(len(constant_set))]
        index_4, val = min(enumerate(temp), key = lambda x: x[1])
        if val == float("inf"):
            return list()
        else:
            indiv_row_index = rrow.lo_vars[index_4]
            rcol, rrow, constraint, constant_set, obj_func, pivoted = \
                pivot(rcol, rrow, constraint, constant_set, obj_func, pivoted, indiv_row_index, e)

    # Then need to remove all the auxiliaries
    temp_list = [(constant_set[rrow.index[index_3]] if (index_3 in rrow.index) else 0) for index_3 in range(len_obj + len_constants + 1)]
    if abs(temp_list[len_obj+len_constants]) < min_rounding_val:

        # We move auxiliary variable row out of basis
        if len_obj+len_constants in rrow.index:
            index_4, indiv_row_index = rrow.index[len_obj + len_constants], len_obj + len_constants
            e = next(e for e, j in rcol.index.items() if abs(constraint[index_4][j]) > min_rounding_val)
            rcol, rrow, constraint, constant_set, obj_func, pivoted =\
                pivot(rcol, rrow, constraint, constant_set, obj_func, pivoted, indiv_row_index, e)

        # Then afterwards remove the corresponding column for aux
        col = rcol.index[len_obj + len_constants]
        constraint = [row[:col] + row[col + 1:] for row in constraint]
        obj_func = [0] * (len_obj)
        pivoted = 0
        rcol = Iteration(rcol.lo_vars[:col] + rcol.lo_vars[col + 1:])

        # Finally we then reinstate the basic variables into obj function
        for index_4, val in enumerate(obj):
            if abs(val) < min_rounding_val: continue
            elif index_4 in rrow.lo_vars:
                row_index = rrow.index[index_4]
                temp = [-val * x for x in constraint[row_index]]
                for index_2 in range(len(obj_func)):
                    obj_func[index_2] += temp[index_2]
                pivoted += val * constant_set[row_index]
            else:
                row_index = rcol.index[index_4]
                obj_func[row_index] += val
        return rcol, rrow, constraint, constant_set, obj_func, pivoted
    else:
        # There is no solution, thus must be infeasible.
        # Commented out print here as line 72 will (should) always activate instead
        #print("infeasible")
        return None




def main():
    '''
    This file, main.py, is called via command line of "py main.py <INSERT TEXT FILE>". The text files are the provided
    .txt files in the input folder of CSC 445 Programming Project
    :return: Nothing. Prints either, infeasible, unbounded,
                or optimal \n <OPTIMAL VALUE> \n <OPTIMAL Xn VALUES>, depending on the solution
    '''
    result = filereader()
    objective_func, constraints, constants = result
    return_val = simplex_start(constraints, constants, objective_func)
    result2 = filereader()
    objective_vals = result2[0]
    if return_val != None:
        do_sum = zip(return_val, objective_vals)
        print(sum(x*y for x,y in do_sum))

        # Specifications asked for minimum 7 sigfig, so I'm just gonna do 9 to be safe
        rounded = ([round(num, 9) for num in return_val])
        for val in rounded:
            print(val, end=' ')






if __name__ == '__main__':
    main()