from abc import ABC, abstractmethod
import re


class Equation(ABC):
    """
    Abstract base class for all types of equations.

    Each equation has a degree and a type. The degree is the highest power of the
    variable(s) in the equation, and the type is a string that describes the
    equation (e.g. "Linear Equation", "Quadratic Equation", etc.).
    """
    degree: int
    type: str

    def __init__(self, *args):
        """
        Initializes an Equation object.

        Args:
            *args: coefficients of the equation, in descending order of degree

        Raises:
            TypeError: if the number of coefficients is not equal to the degree
                of the equation plus one
            TypeError: if any of the coefficients are not of type 'int' or 'float'
            ValueError: if the highest degree coefficient is zero
        """
        if (self.degree + 1) != len(args):
            raise TypeError(
                f"'Equation' object takes {self.degree + 1} positional arguments but {len(args)} were given"
            )
        if any(not isinstance(arg, (int, float)) for arg in args):
            raise TypeError("Coefficients must be of type 'int' or 'float'")
        if args[0] == 0:
            raise ValueError("Highest degree coefficient must be different from zero")
        self.coefficients = {(len(args) - n - 1): arg for n, arg in enumerate(args)}

    def __init_subclass__(cls):
        """
        Raises errors if the subclass does not define the required attributes.
        """
        if not hasattr(cls, "degree"):
            raise AttributeError(
                f"Cannot create '{cls.__name__}' class: missing required attribute 'degree'"
            )
        if not hasattr(cls, "type"):
            raise AttributeError(
                f"Cannot create '{cls.__name__}' class: missing required attribute 'type'"
            )

    def __str__(self):
        """
        Converts the equation to a string.

        The string is formatted as a mathematical equation, with the coefficients
        and variables separated by spaces. The equation is then followed by '= 0'.
        """
        terms = []
        for n, coefficient in self.coefficients.items():
            if not coefficient:
                continue
            if n == 0:
                terms.append(f'{coefficient:+}')
            elif n == 1:
                terms.append(f'{coefficient:+}x')
            else:
                terms.append(f"{coefficient:+}x**{n}")
        equation_string = ' '.join(terms) + ' = 0'
        return re.sub(r"(?<!\d)1(?=x)", "", equation_string.strip("+"))

    @abstractmethod
    def solve(self):
        """
        Solves the equation.

        Returns:
            list: A list of solutions to the equation. If the equation has no
                solutions, an empty list is returned.
        """
        pass

    @abstractmethod
    def analyze(self):
        """
        Analyzes the equation and returns information about it.

        Returns:
            dict: A dictionary containing information about the equation.
        """
        pass


class LinearEquation(Equation):
    """
    Represents a linear equation in the form of
    ax + b = 0.
    """
    degree = 1
    type = 'Linear Equation'

    def solve(self):
        """
        Solves the equation.

        Returns:
            float: The solution to the equation
        """
        a, b = self.coefficients.values()
        x = -b / a
        return [x]

    def analyze(self):
        """
        Analyzes the equation and returns information about it.

        Returns:
            dict: A dictionary containing the slope and intercept of the
            equation, as well as other information.
        """
        slope, intercept = self.coefficients.values()
        return {'slope': slope, 'intercept': intercept}


class QuadraticEquation(Equation):
    """
    Represents a quadratic equation in the form of
    ax^2 + bx + c = 0.
    """
    degree = 2
    type = 'Quadratic Equation'

    def __init__(self, *args):
        """
        Initializes a QuadraticEquation object.

        Args:
            a (int or float): Coefficient of the x^2 term
            b (int or float): Coefficient of the x term
            c (int or float): Constant term
        """
        super().__init__(*args)
        a, b, c = self.coefficients.values()
        self.delta = b ** 2 - 4 * a * c

    def solve(self):
        """
        Solves the equation.

        Returns:
            list: A list of solutions to the equation. If the equation has no
            solutions, an empty list is returned.
        """
        if self.delta < 0:
            return []
        a, b, _ = self.coefficients.values()
        x1 = (-b + (self.delta) ** 0.5) / (2 * a)
        x2 = (-b - (self.delta) ** 0.5) / (2 * a)
        if self.delta == 0:
            return [x1]

        return [x1, x2]

    def analyze(self):
        """
        Analyzes the equation and returns information about it.

        Returns:
            dict: A dictionary containing the following information about the
            equation:
                x (float): The x-coordinate of the vertex of the parabola
                y (float): The y-coordinate of the vertex of the parabola
                min_max (str): Whether the parabola has a minimum or maximum
                concavity (str): Whether the parabola opens upwards or downwards
        """
        a, b, c = self.coefficients.values()
        x = -b / (2 * a)
        y = a * x ** 2 + b * x + c
        if a > 0:
            concavity = 'upwards'
            min_max = 'min'
        else:
            concavity = 'downwards'
            min_max = 'max'
        return {'x': x, 'y': y, 'min_max': min_max, 'concavity': concavity}


def solver(equation):
    """
    This function takes an Equation object and returns a string which includes
    the equation, its solutions, and additional details about the equation.
    """
    if not isinstance(equation, Equation):
        raise TypeError("Argument must be an Equation object")

    output_string = f'\n{equation.type:-^24}\n'
    output_string += f'\n\n{equation!s:^24}\n\n'
    output_string += f'{"Solutions":-^24}\n\n'

    results = equation.solve()
    result_list = []
    match results:
        # If there are no real roots, create a list with a single string
        case []:
            result_list = ['No real roots']
        # If there is only one real root, create a list with a single string
        case [x]:
            result_list = [f'x = {x:+.3f}']
        # If there are two real roots, create a list with two strings
        case [x1, x2]:
            result_list = [f'x1 = {x1:+.3f}', f'x2 = {x2:+.3f}']

    for result in result_list:
        output_string += f'{result:^24}\n'

    output_string += f'\n{"Details":-^24}\n\n'

    details = equation.analyze()
    details_list = []

    match details:
        case {'slope': slope, 'intercept': intercept}:
            # For linear equations, we have the slope and intercept
            details_list = [f'slope = {slope:>16.3f}', f'y-intercept = {intercept:>10.3f}']
        case {'x': x, 'y': y, 'min_max': min_max, 'concavity': concavity}:
            # For quadratic equations, we have the x and y coordinates of the vertex
            # and the concavity and whether the vertex is a minimum or a maximum
            coord = f'({x:.3f}, {y:.3f})'
            details_list = [f'concavity = {concavity:>12}', f'{min_max} = {coord:>18}']

    for detail in details_list:
        output_string += f'{detail}\n'

    return output_string
