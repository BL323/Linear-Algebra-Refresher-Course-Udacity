from decimal import Decimal, getcontext
from copy import deepcopy
from vector import Vector
from plane import Plane

getcontext().prec = 30


class LinearSystem(object):

    ALL_PLANES_MUST_BE_IN_SAME_DIM_MSG = 'All planes in the system should live in the same dimension'
    NO_SOLUTIONS_MSG = 'No solutions'
    INF_SOLUTIONS_MSG = 'Infinitely many solutions'

    def __init__(self, planes):
        try:
            d = planes[0].dimension
            for p in planes:
                assert p.dimension == d

            self.planes = planes
            self.dimension = d

        except AssertionError:
            raise Exception(self.ALL_PLANES_MUST_BE_IN_SAME_DIM_MSG)


    def swap_rows(self, row1, row2):
        # should add check here to ensure row 1 & 2 are not out of range
        tmp = self[row1]
        self[row1] = self[row2]
        self[row2] = tmp       

    def multiply_coefficient_and_row(self, coefficient, row):
        plane = self[row] 
        scaledConstTerm = plane.constant_term * coefficient
        scaledNormalVec = plane.normal_vector.multiply(coefficient)
        self[row] = Plane(scaledNormalVec, scaledConstTerm)


    def add_multiple_times_row_to_row(self, coefficient, row_to_add, row_to_be_added_to):
        planeToMutate = self[row_to_add]
        targetPlane = self[row_to_be_added_to]       
        constTerm = (planeToMutate.constant_term * coefficient) + targetPlane.constant_term
        normalVec = planeToMutate.normal_vector.multiply(coefficient).plus(targetPlane.normal_vector)
        self[row_to_be_added_to] = Plane(normalVec, constTerm)
        
    def swap_with_row_below_for_nonZero_coeff_if_able(self, row, col):
        num_equations = len(self)

        for k in range(row+1, num_equations):
            coeff = MyDecimal(self[k].normal_vector[col])
            if(not coeff.is_near_zero()):
                self.swap_rows(row, k)
                return True
        
        return False
        

    def clear_coefficients_below(self, row, col):
        num_equations = len(self)
        beta = MyDecimal(self[row].normal_vector[col])    

        for k in range(row+1, num_equations):
            n = self[k].normal_vector
            gamma = n[col]
            alpha = -gamma/beta
            self.add_multiple_times_row_to_row(alpha, row, k)

    def clear_coefficients_above(self, row, col):
        for k in range(row)[::-1]:
            n = self[k].normal_vector
            alpha = -(n[col])
            self.add_multiple_times_row_to_row(alpha, row, k)
    
    def compute_triangular_form(self):
        system = deepcopy(self)

        num_equations = len(system)
        num_variables = system.dimension
        j = 0
        for i in range(num_equations):
            while(j < num_variables):
                # coeff in row i, j-th term
                c = system[i].normal_vector[j]
                if(MyDecimal(c).is_near_zero()):
                    swap_succeeded = system.swap_with_row_below_for_nonZero_coeff_if_able(i, j)  
                    if not swap_succeeded:
                        j += 1
                        continue

                system.clear_coefficients_below(i, j)
                j += 1
                break        
        return system

    def scale_row_to_make_coefficient_equal_one(self, row, col):
        n = self[row].normal_vector
        beta = Decimal('1.0') / n[col]
        self.multiply_coefficient_and_row(beta, row)
         
    def compute_rref(self):
        tf = self.compute_triangular_form()
        num_equations = len(tf)

        indexes = tf.indices_of_first_nonzero_terms_in_each_row()

        for i in range(num_equations)[::-1]:
            j = indexes[i]
            if j < 0:
                continue
            
            tf.scale_row_to_make_coeff_sinuglar(i, j)
            tf.clear_coefficients_above(i, j)

        return tf

    def scale_row_to_make_coeff_sinuglar(self, row, col):
        n = self[row].normal_vector
        beta = Decimal('1.0') / n[col]
        self.multiply_coefficient_and_row(beta, row)

    def indices_of_first_nonzero_terms_in_each_row(self):
        num_equations = len(self)
        num_variables = self.dimension

        indices = [-1] * num_equations

        for i,p in enumerate(self.planes):
            try:
                indices[i] = p.first_nonzero_index(p.normal_vector)
            except Exception as e:
                if str(e) == Plane.NO_NONZERO_ELTS_FOUND_MSG:
                    continue
                else:
                    raise e

        return indices

    def perform_gaussian_elimination(self):
        rref = self.compute_rref()
        
        # check rref for contradictions or too few pivots
        rref.raise_if_contradictions_exist()
        rref.raise_if_too_few_pivots()    
        
        num_variables = rref.dimension
        solution_coordinates = [rref.planes[i].constant_term for i in range(num_variables)]
    
        return Vector(solution_coordinates);

    def compute_solution(self):
        try:
            return self.perform_gaussian_elimination();
        except Exception as ex:
            if str(ex) == self.NO_SOLUTIONS_MSG:
                return self.NO_SOLUTIONS_MSG
            elif str(ex) == self.INF_SOLUTIONS_MSG:
                return self.INF_SOLUTIONS_MSG
            else:
                raise ex
            
        
    def raise_if_contradictions_exist(self):
        for p in self.planes:
            try:
                p.first_nonzero_index(p.normal_vector)
                
            except Exception as ex:
                if str(ex) == "No nonzero elements found":
                    constant_term = MyDecimal(p.constant_term)
                    if not constant_term.is_near_zero():
                        raise Exception(self.NO_SOLUTIONS_MSG)
                else:
                    raise ex
                    
    def raise_if_too_few_pivots(self):
        pivot_indicies = self.indices_of_first_nonzero_terms_in_each_row()
        num_pivots = sum([1 if index >= 0 else 0 for index in pivot_indicies])
        num_variables = self.dimension
        
        if(num_pivots < num_variables):
            raise Exception(self.INF_SOLUTIONS_MSG)

    def __len__(self):
        return len(self.planes)


    def __getitem__(self, i):
        return self.planes[i]


    def __setitem__(self, i, x):
        try:
            assert x.dimension == self.dimension
            self.planes[i] = x

        except AssertionError:
            raise Exception(self.ALL_PLANES_MUST_BE_IN_SAME_DIM_MSG)


    def __str__(self):
        ret = 'Linear System:\n'
        temp = ['Equation {}: {}'.format(i+1,p) for i,p in enumerate(self.planes)]
        ret += '\n'.join(temp)
        return ret


class MyDecimal(Decimal):
    def is_near_zero(self, eps=1e-10):
        return abs(self) < eps


