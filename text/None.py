from cvxopt.base import matrix

Vb = ['10.5', '11.4']
if isinstance(Vb[0], str):
    for item in Vb:
        Vb = float(Vb[item])

print(Vb)