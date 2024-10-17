from equation import (LinearEquation, QuadraticEquation, solver)

lin_eq = LinearEquation(2, 3)
quad_eq = QuadraticEquation(1, 2, 1)

print(solver(lin_eq))
print(solver(quad_eq))
