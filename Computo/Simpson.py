# Daianna Gonzalez Padilla.
# Computo cientifico, 2do Semestre 2021, LCG: Proyecto Final
# Programa 7: Metodo de integracion por la regla de Simpson (1/3)

# Diccionario de variables:

# empezar permite reiniciar el programa
# n es el numero de bandas para calcular el area bajo la curva
# fx guarda las funcion dada como input, f toma tal input y lo ve como una expresion en terminos de x
# a es el limite inferior del intervalo de busqueda y b el superior
# h es el tamano de los subintervalos o bandas
# d es el numero de decimales de exactitud
# k es la constante que multiplica a cada yi
# y es cada f(a+hi) con h el incremento
# A es la aproximacion del area bajo la curva de la funcion f en el intervalo a--b
# i es un contador
# s es un acumulador
# xi y y permiten graficar la funcion tabulando en un intervalo dado
# X y Y son arreglos para mostrar el area bajo la curva en el intervalo a--b


import sympy as sp
from sympy import *
import numpy as np
import matplotlib.pyplot as plt


print("Programa para aproximar la integral una funcion en un intervalo\n")

# Cada que se oprima 1 se empieza el programa
empezar = int(input("Para comenzar oprima 1\n"))
while empezar == 1:

    # Funcion para pedir el numero par de subintervalos, si no es par se pide otro
    def n_bandas():
        n = int(input("Ingrese el numero de bandas, debe ser par:\n"))
        if not n % 2 == 0:
            print("El numero de bandas debe ser par, intente nuevamente\n")
            n = n_bandas()
        return (n)


    # Se pide la funcion f(x) que se quiere integrar y el intervalo a-b para ello, asi como el numero de bandas
    fx = input("Inserte la funcion que quiere integrar, asegurese de usar las expresiones correctas y permitidas: \n")
    x = sp.symbols('x')
    f = simplify(fx)
    a = float(input("Ingrese el limite inferior del intervalo:\n"))
    b = float(input("Ingrese el limite superior  del intervalo:\n"))
    n = n_bandas()
    d = int(input("Numero de decimales:"))

    # Se calcula h: el tamano de los subintervalos
    h = (b - a) / n
    # Se efectua la suma yo + 4y1 + 2y2 + 4y3 + ... + 4y(n-1) +yn, en total n+1 terminos
    s = 0
    i = 0
    while i <= n:
        # Cada yi es f(a+hi): se evalua la funcion para cada x empezando en a e incrementando en h
        y = round(f.evalf(subs={'x': (a + (h * i))}),d+1)

        # Los coeficientes k de yo y yn son 1, de los yi con i par son 2 e impar 4
        if i == 0 or i == n:
            k = 1
        elif i % 2 == 0:
            k = 2
        elif not i % 2 == 0:
            k = 4
        s = s + (y * k)
        i = i + 1
    # Se calcula el area total sumando las areas bajo las curvas para cada pareja de subintervalos y al final
    # A=h/3* suma de los k(yi)
    A = round((h / 3) * s,d)
    print("La integral de la funcion en el intervalo dado es aprox: ", A)

    #Se imprime la grafica de la integral de la funcion en tal intervalo
    #Se evaluan valores de x en f para graficar la funcion
    try:
        xi = np.linspace(a-1, b+1, 1000)
        y = [f.evalf(subs={'x': i}) for i in xi]
        plt.plot(xi, y)
        plt.axhline(0, color="black")
        plt.axvline(0, color="black")
        #Arreglos con los puntos en la curva que delimitan el area bajo la curva
        X = np.arange(int(a), int(b)+0.1, 0.1)
        X = np.array(X, dtype=float)
        Y = [f.evalf(subs={'x': i}) for i in X]
        Y = np.array(Y, dtype=float)
        # Se definen los limites de la grafica mostrada
        plt.xlim(int(a) - 1, int(b) + 1)
        plt.ylim(float(min(Y))-1, float(max(Y))+1)
        #Se rellena el area bajo la curva en el intervalo de a -> b
        plt.fill_between(X, Y, color='r')
        plt.title('Integral definida de la funcion', fontsize=10)
        plt.show()
    # Si la funcion no es continua en el a-1, b+1, no se imprime la grafica
    except TypeError:
            pass



    empezar = int(input("Para ingresar otra funcion oprima 1, para salir 2\n"))
print("Gracias por usar este programa")