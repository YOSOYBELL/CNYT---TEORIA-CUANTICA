import math
import numpy

"""1. El sistema debe calcular la probabilidad de encontrarlo en una posición en particular."""


def sumaCuadrados(c):
    suma = abs(c) ** 2
    return suma


def normaVector(v):
    suma = 0
    for i in range(len(v)):
        suma += sumaCuadrados(v[i])
    return math.sqrt(suma)


def superposition(vec, pos):
    """
        Retorna la probabilidad de posicion dado un vector
        (list) -> float
    """
    c = sumaCuadrados(vec[pos])
    norma = normaVector(vec)
    super = c / norma ** 2
    return round(super, 4)


# print("Probabilidad: ", superposition(vec=[(-3 - 1j), (-2j), 1j, 2], pos=2))

"""2. El sistema si se le da otro vector Ket debe buscar la probabilidad de transitar del primer vector al segundo. """


def transition(vec1, vec2):
    """ Retorna la amplitud de una transición
        (list), (list) -> complex"""
    productoInterno = numpy.inner(vec1, vec2)
    normaVec1 = normaVector(vec1)
    normaVec2 = normaVector(vec2)
    ket = productoInterno / (normaVec1 * normaVec2)
    return numpy.round(ket)


# print("Amplitud de transición: ", transition(vec1=[(-1j), 1], vec2=[1, (-1j)]))

########################################################################################

"""2. Ahora con una matriz que describa un observable y un vector ket, el sistema revisa que la matriz sea
hermitiana, y si lo es, calcula la media y la varianza del observable en el estado dado. """


def tranpuestaMatriz(mat):
    # 7. Transpuesta de una matriz/vector
    for i in range(len(mat)):
        for j in range(len(mat)):
            if j < i:
                nueva = mat[i][j]
                mat[i][j] = mat[j][i]
                mat[j][i] = nueva
    return mat


def conjugadaMatriz(mat):
    matConj = numpy.conjugate(mat)
    return matConj


def adjuntaMatriz(mat):
    # 9. Adjunta (daga) de una matriz/vector
    matriz = conjugadaMatriz(mat)
    matriz = tranpuestaMatriz(mat)
    return matriz


def probarHermitiana(mat):
    adjunta = adjuntaMatriz(mat)
    if mat == adjunta:
        return True
    else:
        return False


def unitaria(mat):
    return numpy.allclose(numpy.eye(mat.shape[0]), mat.H * mat)


def proobarHermitiana(mat):
    conj = numpy.conjugate(mat)
    trans = numpy.transpose(conj)
    boo = numpy.array_equal(mat, trans)
    return boo


def accionMatrizVector(v, mat):
    v = numpy.array(v)
    mat = numpy.array(mat)
    return mat.dot(v)


def innerProduct(vec1, vec2):
    pd = 0
    for i in range(len(vec1)):
        pd += vec2[i] * vec1[i].conjugate()
    return numpy.round(pd, 2)


def media(observable, vecEstado):
    media = 0
    ishermitiana = proobarHermitiana(observable)
    if ishermitiana == True:
        prod = accionMatrizVector(observable, vecEstado)
        media = innerProduct(prod, vecEstado)
    return media


obs = [[1, -1j], [1j, 2]]
vecEstado = [math.sqrt(2) / 2 + 0j, (math.sqrt(2) / 2) * 1j]


# print("La media es: ", media(obs, vecEstado))


def varianza(observable, vecEstado):
    mediaC = media(observable, vecEstado)
    newMat = numpy.identity(len(observable))
    for i in range(len(newMat)):
        for j in range(len(newMat)):
            newMat[i][j] = newMat[i][j] * mediaC
    resta = numpy.matrix(observable) - numpy.matrix(newMat)
    resta = numpy.array(resta)
    producto = numpy.dot(resta, resta)
    var = media(producto, vecEstado)
    return var


# print("La varianza es: ", varianza(obs, vecEstado))

"""El sistema calcula los valores propios del observable y la probabilidad de que el sistema
    transite a alguno de los vectores propios después de la observación. """


def valPropio(mat):
    mat = numpy.matrix(mat)
    return numpy.linalg.eig(mat)


def probabilidadTransicion(observable, posicion):
    """
    (list) -> float
    """
    valoresPropios = valPropio(observable)
    probabilidad = superposition(valoresPropios, posicion)
    return probabilidad


# print("Los valores propios son: ", valPropio(obs), "y la probabilidad es ", probabilidadTransicion(obs, 2) )

""" COMPLETAR LIBRERIA """


def ejercicio_431(observable):
    """ Find all the possible states the system described in Exercise 4.2.2
    can transition into after a measurement has been carried out.

    (list) -> float """
    spinX = [1, 0]
    spinX = numpy.array(spinX)
    prod = accionMatrizVector(observable, spinX)
    probabilidad = superposition(prod, 0)
    return probabilidad


observable = [[0, 1/2], [1/2, 0]]
# print(ejercicio431(observable))


def ejercicio_432(observable):
    """Perform the same calculations as in the last example, using Exercise 4.3.1.
    Then draw the probability distribution of the eigenvalues as in the previous example.
    (list) -> float"""

    spinX = [1, 0]
    spinX = numpy.array(spinX)
    valores_propios = valPropio(observable)
    vectores = normaVector(valores_propios)
    producto = accionMatrizVector(observable, spinX)
    p1 = normaVector(innerProduct(producto, valores_propios[0]))
    p2 = normaVector(innerProduct(producto, valores_propios[1]))

    sum = p1 * valores_propios[0] + p2 * valores_propios[1]
    return sum

# print(ejercicio_432(observable))

def ejercicio_441(u1, u2):
    """Verify that are unitary matrices. Multiply them and verify that their product is also unitary.
    (list), (list) -> boolean, boolean, boolean"""
    u1 = numpy.matrix(u1)
    u2 = numpy.matrix(u2)
    isU1 = unitaria(u1)
    isU2 = unitaria(u2)
    producto = numpy.dot(u1, u2)
    productoUni = unitaria(producto)
    return isU1, isU2, productoUni


u1 = [[0, 1], [1, 0]]
u2 = [[math.sqrt(2) / 2, math.sqrt(2) / 2], [math.sqrt(2) / 2, (math.sqrt(2) / 2) * -1]]


# print(ejercicio_441(u1, u2))


def ejercicio_442(estadoInicial, mat):
    """Go back to Example 3.3.2 (quantum billiard ball), keep the same
    initial state vector [1, 0, 0, 0]T, but change the unitary map to """
    matriz = numpy.matrix(mat)
    estado = accionMatrizVector(estadoInicial, matriz)
    probabilidad = probabilidadTransicion(estado, 3)
    return probabilidad

unitaryMap = [[0, 1 / math.sqrt(2), 1 / math.sqrt(2), 0], [(1 / math.sqrt(2)) * 1j, 0, 0, 1 / math.sqrt(2)],
              [1 / math.sqrt(2), 0, 0, (1 / math.sqrt(2)) * 1j], [0, 1 / math.sqrt(2), (1 / math.sqrt(2)) * -1, 0]]
initialVector = [1, 0, 0, 0]

# print(ejercicio_442(initialVector, unitaryMap))