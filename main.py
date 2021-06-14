"""https://github.com/yardenho/part4_q35.git"""

def neville(pointsList, m, n, X):
    """
    :param pointsList: a list of points
    :param m: the start point index
    :param n: the end point index
    :param X: an x value for it we want to find the y value in the function that is created from the given points
    :return: a proximity of f(x)
    """
    if m is n:  # the start index and end index are the same (the same point)
        print("P" + str(m) + " = " + str(pointsList[m][1]))
        return pointsList[m][1]   # return the y value of the point
    else:
        Xm = pointsList[m][0]  # the start point x value
        Xn = pointsList[n][0]  # the end point x value
        res = (X - Xm) * neville(pointsList, m+1, n, X) - (X - Xn) * neville(pointsList, m, n-1, X)
        res /= (Xn - Xm)
        print("P" + str(m) +"," + str(n) + " = " + str(res))
        return res   # the result (f(X))


# polynomial method

def newMat(numRow, numCol):
    """
    :param numRow: the number of rows in the mat
    :param numCol: the number of columns in the mat
    :return: a zero matrix in the required size
    """
    mat = []   # the zero matrix the function returns
    for i in range(numRow):
        mat.append([])    # create a new row
        for j in range(numCol):
            mat[i].append(0)    # fill the row with
    return mat


def printMatrix(a, name=None):
    """
    :param a: a matrix to print
    :return: prints in matrix format
    """
    print("-----------------------------")
    if name is not None:
        print(name + " = ")
    for i in range(len(a)):
        if i is len(a) - 1:
            print(" " + str(a[i]) + "]")
        elif i is 0:
            print("[" + str(a[i]))
        else:
            print(" " + str(a[i]))
    print("")

#start LUD composition

def LUdecomposition(mat, b):
    """
    :param mat: the coefficients matrix
    :param b: the solution vector
    :return: none
    """
    inversL, L, counter = toUpperTriangularMat(mat)  # calculate the L matrix and the inverse L
    print("-----------------------------")
    string = "L = "
    for i in range(1, counter):
        if i < counter - 1:
            string += "invJ" + str(i) + " * "
        else:
            string += "invJ" + str(i)
    print(string)
    printMatrix(L, "L")  # print the L matrix
    print("-----------------------------")
    string = "inverse L = "
    for i in range(counter - 1, 0, -1):
        if i > 1:
            string += "J" + str(i) + " * "
        else:
            string += "J" + str(i)
    print(string)
    printMatrix(inversL, "inverse L")
    print("-----------------------------")
    print("U = invL * A")
    U = multMat(inversL, mat)  # calculate the U matrix
    inversU, counter = FromUpperToInvers(U, createIdentityMatrix(len(U)))  # calculate thr inverse of U
    print("-----------------------------")
    string = "inverse U = "
    for i in range(counter - 1, 0, -1):
        string += "E" + str(i) + " * "
    string += "U"
    print(string)
    printMatrix(U, "U")  # print the U matrix
    print("-----------------------------")
    print("x = invU * invL * b")
    x = multMat(inversL, b)  # finding the result vector
    x = multMat(inversU, x)  # finding the result vector
    printMatrix(x, "x")  # print the X matrix
    return x


def FromUpperToInvers(A, inverseMat, counter=1):
    """
    :param A: upper matrix
    :param inverseMat: the matrix that will become the inverse
    :return: Inverse matrix
    """
    elemntarMat = createIdentityMatrix(len(A))  # identity matrix
    for i in range(len(A) - 1, -1, -1):  # run over the columns
        for j in range(i):  # run over the lines above the pivot
            elemntarMat[j][i] = -(A[j][i] / A[i][i])
            A = multMat(elemntarMat, A)
            inverseMat = multMat(elemntarMat, inverseMat)
            printMatrix(elemntarMat, "E" + str(counter))
            counter += 1
            elemntarMat[j][i] = 0
        if A[i][i] != 1:  # convert the pivots to one
            elemntarMat[i][i] = 1 / A[i][i]
            A = multMat(elemntarMat, A)
            inverseMat = multMat(elemntarMat, inverseMat)
            printMatrix(elemntarMat, "E" + str(counter))
            counter += 1
            elemntarMat[i][i] = 1
    return inverseMat, counter


def toUpperTriangularMat(A):
    """
    :param A: the matrix we want to turn into a upper triangle matrics
    :return: the multiplication of the elementary matrics that create U
    """
    L = createIdentityMatrix(len(A))  # create indetity matrix
    InL = createIdentityMatrix(len(A))  # create indetity matrix
    counter = 1
    for i in range(len(A)):  # run over the lines
        if A[i][i] == 0:  # if the pivot is 0
            for j in range(i + 1, len(A)):  # run over the columns
                if A[j][i] != 0:  # if the element under the pivot is not 0
                    elementary = linesExchange(A, i, j)
                    L = multMat(L, elementary)  # make lines exchange and multiply
                    printMatrix(elementary, "J" + str(counter))
                    InL = multMat((linesExchange(A, i, j)), InL)  # make lines exchange and multiply
                    printMatrix(identity, "inverse J" + str(counter))
                    A = multMat((linesExchange(A, i, j)), A)
                    counter += 1
                    break
        if A[i][i] != 0:  # check if B is regular
            for j in range(i + 1, len(A)):  # run over the columns
                identity = createIdentityMatrix(len(A))
                identity[j][i] = -(A[j][i] / A[i][i])  # elementary matrix
                printMatrix(identity, "J" + str(counter))
                InL = multMat(identity, InL)  # L^(-1)
                A = multMat(identity, A)
                identity[j][i] *= -1  # changing the element in order to find L
                printMatrix(identity, "inverse J" + str(counter))
                L = multMat(L, identity)
                counter += 1
    return InL, L, counter


def linesExchange(A, line1, line2):
    """
    :param A: A matrix
    :param line1: A line
    :param line2: The line we want to exchange with
    :return: elementry matrix
    """
    idendityMax = createIdentityMatrix(len(A))  # create identity matrix

    # exchange the members in line1
    temp = idendityMax[line1][line1]
    idendityMax[line1][line1] = idendityMax[line2][line1]
    idendityMax[line2][line1] = temp

    # exchange the members in line2
    temp = idendityMax[line2][line2]
    idendityMax[line2][line2] = idendityMax[line1][line2]
    idendityMax[line1][line2] = temp
    return idendityMax


def createIdentityMatrix(size):
    """
    :param size: the size of the square matrix
    :return:
    """
    identityMat = newMat(size, size)  # create a zero matrix in the required size
    for index in range(size):  # go over the main diagonal
        identityMat[index][index] = 1  # change the elements in the main diagonal to 1
    return identityMat


def multMat(A, B):
    """
    :param A: a matrix in sise n*m
    :param B: a mtrix in size m*k
    :return: A*B  (in size n*k)
    """
    if len(A[1]) == len(B):  # check if A and B have the same number of rows and columns
        C = newMat(len(A), len(B[0]))  # the matrix the function returns
        for i in range(len(C)):
            for j in range(len(C[1])):
                for k in range(len(B)):
                    C[i][j] += A[i][k] * B[k][j]
        return C
    else:
        return None  # the multiplication  is impossible


def polynomial(pointsList, X):
    """
    :param pointsList: the list of the points
    :param X: the point that we want to find her approximate value
    :return: the approximate value of X
    """
    mat = newMat(len(pointsList), len(pointsList))
    for i in range(len(mat)):
        for j in range(len(mat[0])):
            mat[i][j] = pow(pointsList[i][0], j)  # The coefficient matrix
    print("A = ")
    printMatrix(mat)
    b = newMat(len(pointsList), 1)
    for i in range(len(b)):
        b[i][0] = pointsList[i][1]  # The column vector of constant terms
    print("b = ")
    printMatrix(b)
    matRes = LUdecomposition(mat, b)  # returns the solution matrix
    print("p" + str(len(pointsList)) + "(x) = ", end="")
    string = ""
    for i in range(len(matRes)):
        string += str(matRes[i][0])
        if i is not 0:
            string += " * x^" + str(i)
        if i is not len(matRes) - 1:
            string += " + "
    print(string + "\n")
    # calc mat
    res = 0
    for i in range(len(matRes)):
        res += matRes[i][0] * pow(X, i)  # calc the y value for the requested x
    return res


def calcFinalResult(result, epsilon, day, hour, minutes):
    """
    :param result: the result
    :param epsilon: the epsilon we decided on for the question
    :param day: the day of the calculation
    :param hour: the hour of the calculation
    :param minutes: the minutes of the calculation
    :return: the result by the requested format
    """
    stringRes = str(result)  # cast the result to string
    i = 0
    while stringRes[i] is not ".":  # run over the string while we get to the point
        i += 1  # count how many digits there is before the point
    i += 1
    count = 1
    while epsilon < 1:  # checking how digit needs after the point
        epsilon *= 10
        count += 1
    differ = len(stringRes)
    while differ < i + count:   # fill difference with zero if the number is to short after the point
        stringRes += "0"
        differ += 1
    stringRes = stringRes[:i + count] + "00000" + day + hour + minutes
    return stringRes


def Driver():
    point_list = [[1.2, 3.5095], [1.3, 3.6984], [1.4, 3.9043], [1.5, 4.1293], [1.6, 4.3756]]
    X = 1.37
    epsilon = 10 ** -4
    print("==== Neville Method ====")
    n = neville(point_list, 0, len(point_list) - 1, X)
    print("Final result:\nf(" + str(X) + ") = " + str(calcFinalResult(n, epsilon, '13', '18', '59')))
    print("\n==== Polynomial Method ====\n")
    p = polynomial(point_list, X)
    print("Final result:\nf(" + str(X) + ") = " + str(calcFinalResult(p,epsilon ,'13', '18', '59')))

    if abs(n - p) < epsilon:
        print("\n* The difference between the two methods is smaller than the epsilon")
    else:
        print("\n* The difference between the two methods is bigger than the epsilon - needs another method")


Driver()