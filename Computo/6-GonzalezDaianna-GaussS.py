# Daianna Gonzalez Padilla.
# Computo cientifico, 2do Semestre 2021, LCG: Proyecto Final
# Programa 6: Metodo de Gauss-Seidel para resolucion de sistemas de ecuaciones lineales

# Diccionario de variables:

# empezar permite reiniciar el programa
# n es el numero de incognitas de las ecuaciones e igualmente el numero de ecuaciones
# k es el numero de decimales de exactitud
# coef es un diccionario que contiene los indices de cada entrada con su respectivo valor (coeficiente)
# incog es un diccionario con las incognitas Xi del sistema
# result es un diccionario con el numero de cada ecuacion y su valor del lado derecho de la igualdad
# paro es una bandera para determinar cuando se deja de iterar
# i, j y e son contadores e indican indices
# suma es un acumulador
# Xi1 guarda el valor de la incognita  Xi obtenida una iteracion antes de la actual
# converg es una bandera para determinar si el sistema es diagonalmente dominante y por tanto convergente




print("Programa para resolver sistemas de ecuaciones lineales con solucion unica por Gauss Seidel\n")
# Cada que se oprima 1 se empieza el programa
empezar = int(input("Para comenzar oprima 1\n"))
while empezar == 1:

    coef = {}
    incog = {}
    result = {}

    #Funcion para determinar si el sistema es diagonalmente dominante
    def convergencia(coef):
        i=1
        s=0
        converg=True
        # El sistema converge si |aii|> suma(|aij|) con aii el coeficiente de la diagonal y aij los restantes en la
        # ecuacion i, si lo es se regresa 1, si no 0
        while i<=n:
            j=1
            while j<=n:
                if not i==j:
                    s=s+abs(coef[(i,j)])
                j=j+1
            if coef[(i,i)]<= s:
                converg=False
                break
            i=i+1
        if converg==True:
            return (1)
        else:
            return (0)



    # Se pide la cantidad de ecuaciones del sistema (y de incognitas)
    print("Ingrese el sistema, este debe tener solucion unica\n")
    n = int(input("Ingrese la cantidad de ecuaciones/incognitas del sistema: "))
    # Se piden los coeficientes de las ecuaciones y se guardan en el diccionario coef junto a sus indices i,j
    # correspondientes
    i = 1
    while i <= n:
        j = 1
        while j <= n:
            print("Ingrese el coeficiente de la incognita ", j, "de la ecuacion ", i)
            coef[(i, j)] = float(input())
            j = j + 1
        # Se piden los resultados de las i ecuaciones y se guardan en result
        print("Ingrese el resultado de la ecuacion", i)
        result[i] = float(input())
        # Se inicializan las incognitas en cero para la iteracion 0, se guardan en incog donde cada incog[i] es un Xi
        incog[i] = 0
        i = i + 1
    # Se pide el numero de decimales exactos
    k = int(input("Ingrese el numero de decimales exactos: "))

    #Funcion para obtener las aproximaciones de las incognitas
    def aproximaciones(incog,coef,result):
        # Se hacen los despejes de las incognitas Xi para cada iteracion:
        # Obtenemos los valores de Xi= ((resultado de la ecuacion i)-(suma de los productos de los coeficientes i,j con
        # las incognitas Xj)) / coeficiente i,i de la diagonal
        paro = False
        # Las iteraciones paran cuando se han alcanzado los k decimales iguales en los valores de las incognitas para
        # dos iteraciones sucesivas
        while paro == False:
            i = 1
            e = 0
            while i <= n:
                j = 1
                suma = 0
                # Se guarda el valor de Xi obtenido en la iteracion m antes de pasar a la m+1
                Xi1 = incog[i]
                while j <= n:
                    # Las incognitas Xi se obtienen en la iteracion m con los Xj de la iteracion m-1
                    if not i == j:
                        suma = suma + (coef[i, j] * incog[j])
                    j = j + 1
                # Valor de la incognita i en la iteracion m:
                incog[i] = (result[i] - suma) / coef[i, i]

                # Si Xi(en m+1)*(10**k)= Xi(en m)*(10**k), es decir, que tienen los mismos k decimales, e incrementa en 1
                if int(incog[i] * (10 ** k)) == int(Xi1 * (10 ** k)):
                    e = e + 1
                i = i + 1

            # Si para todas las incognitas Xi en una iteracion m pasa que e incrementa 1, o sea que todas poseen los
            # mismos k decimales que las obtenidas en la iteracion anterior, entonces e==n y las iteraciones se detienen
            # cuando paro==True
            if e == n:
                paro = True
            else:
                paro = False


        # Se imprimen los valores obtenidos para las n incognitas
        i = 1
        while i <= n:
            print("Incognita X", i, ":", round(incog[i], k))
            i = i + 1

    #Si el sistema es diagonalmente dominante se obtienen las aproximaciones
    if convergencia(coef) == 1:
        aproximaciones(incog,coef,result)

    # Si el sistema no es diagonalmente dominante, se evalua si converge con base en la reduccion o amplificacion del
    # error
    else:
        try:
            aproximaciones(incog,coef,result)
        # Si el sistema no converge, se itera infinitamente y hay un error OverflowError
        except OverflowError as ex:
            print("El sistema no coverge\n")
            empezar = int(input("Para reacomodar las ecuaciones/salir oprima 1\n"))
        #Si alguna entrada en la diagonal se deben reacomodar las ecuaciones o descartar al sistema
        except ZeroDivisionError as ex:
            print("La diagonal no puede contener ceros\n")
            empezar = int(input("Para reacomodar las ecuaciones/salir oprima 1\n"))


    empezar = int(input("Para ingresar otro sistema oprima 1, para salir 2\n"))
print("Gracias por usar este programa")