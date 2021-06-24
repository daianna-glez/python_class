# Daianna Gonzalez Padilla.
# Computo cientifico, 2do Semestre 2021, LCG: Proyecto Final
# Programa 10: Metodo Runge-Kutta para resolver ecuaciones diferenciales

# Diccionario de variables:

# empezar permite reiniciar el programa
# fx guarda las funcion dada como input, f toma tal input y lo ve como una expresion en terminos de x y de y
# xo es el el valor en x del valor inicial
# yo es el valor en y del valor inicial
# h es el valor del ancho de cada banda
# p es el valor en x para el cual se quiere resolver la ecuacion diferencial
# d es el numero de decimales de exactitud
# xi son los x's dados por xi+h en cada iteracion
# yi,yi2 son valores de y's consecutivos dados por y(n+1)=yn + 1/6(k1 + 2k2 + 2k3 + k4)
# k1,k2,k3,k4 son valores dados para cada yi y permiten obtener el siguiente y



import sympy as sp
from sympy import *

print("Programa para aproximar la solucion de una E.D. de 1er orden\n")

# Cada que se oprima 1 se empieza el programa
empezar = int(input("Para comenzar oprima 1\n"))
while empezar == 1:

    # Se pide la funcion del lado dereccho de la igualdad de la E.D; es decir f(x,y)
    fx = input("Inserte la derivada de la funcion que quiere encontrar, asegurese de usar las expresiones correctas "
           "y permitidas: \n")
    x = sp.symbols('x')
    y = sp.symbols('y')
    f = simplify(fx)
    # Se pide el valor inicial y(xo)=y1, xn-x(n+1)=h que determina la precision del resultado, asi como el x para el
    # cual se quiere resolver (p)
    xo = float(input("Ingrese el valor en x de su punto inicial:\n"))
    yo = float(input("Ingrese el valor en y de su punto inicial:\n"))
    h = float(input("Ingrese el valor del ancho de cada banda, debe ser suficientemente pequeno:\n"))
    p = float(input("Ingrese el valor de x con el cual quiere resolver:\n"))
    d = int(input("Numero de decimales:"))

    # Calculo de los yi hasta que se llegue al valor de x dado (p):
    # Inicialmente se toma xo y yo para obtener y1
    xi = xo
    yi = yo
    i = 1
    while not xi == p:
        # Se obtienen los y(n+1)=yn + 1/6(k1 + 2k2 + 2k3 + k4) hasta que xi valga p
        # k1=h*f(xn,yn)
        k1 = round(h * f.evalf(subs={'x': xi, 'y': yi}), d + 1)
        # k2=h*f(xn+h/2 , yn+k1/2)
        k2 = round(h * f.evalf(subs={'x': xi + (h / 2), 'y': yi + (k1 / 2)}), d + 1)
        # k3=h*f(xn+h/2 , yn+k2/2)
        k3 = round(h * f.evalf(subs={'x': xi + (h / 2), 'y': yi + (k2 / 2)}), d + 1)
        # k4=h*f(xn+h, yn+k3)
        k4 = round(h * f.evalf(subs={'x': xi + h, 'y': yi + k3}), d + 1)
        yi2 = round(yi + ((1 / 6) * (k1 + (2 * k2) + (2 * k3) + k4)), d + 1)
        yi = yi2
        xi = round(xi + h, 1)
        i = i + 1


    print("Resultado de la resolucion con el x dado: ", round(yi2,d))
    empezar = int(input("Para ingresar otra derivada oprima 1, para salir 2\n"))
print("Gracias por usar este programa")