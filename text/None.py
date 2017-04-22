from cvxopt.base import matrix, sparse

# Vb = ['10.5', '11.4']
# if isinstance(Vb[0], str):
#     for item in Vb:
#         Vb = float(Vb[item])
#
# print(Vb)
A = [[1,2], [4,5]]

for i in range(2):
    A[i][i] += 1
    print(A)



# A[0][1] += -complex(1, 2)
# print(matrix(A))
#
# print(5/4)
