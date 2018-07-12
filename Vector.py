class Vector(object):
    def __init__(self, coordinates):
        try:
            if not coordinates:
                raise ValueError
            self.coordinates = tuple(coordinates)
            self.dimension = len(coordinates)

        except ValueError:
            raise ValueError('The coordinates must be nonempty')

        except TypeError:
            raise TypeError('The coordinates must be an iterable')

    def __str__(self):
        return 'Vector: {}'.format(self.coordinates)

    def __eq__(self, vec):
        return self.coordinates == vec.coordinates
    
    def Plus(self, vec):
        new_coordinates = [x + y for x,y in zip(self.coordinates, vec.coordinates)]
        return Vector(new_coordinates)
            
    def Minus(self, vec):
        new_coordinates = [x - y for x,y in zip(self.coordinates, vec.coordinates)]
        return Vector(new_coordinates)
            
    def Multiply(self, scalar):
        new_coordinates = [x * scalar for x in self.coordinates]                
        return Vector(new_coordinates)

    def Magnitude(self):
        return sum([x ** 2 for x in self.coordinates])**(1.0/2)

    def Normalisation(self):
        try:
            normalisationFactor = 1 / self.Magnitude()
            return Vector([x * normalisationFactor for x in self.coordinates])
        
        except ZeroDivisionError:
            raise Exception("Cannot normalise a vector with zero values")