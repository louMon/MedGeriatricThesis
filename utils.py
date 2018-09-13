import pandas as pd
import numpy as np

#Generar un arreglo con los tamanhos de las diferentes bases
def generateArrayBase(campos, arrMedicamentos):
    arrBase = []
    cantEstados = 0     
    for j in range(0, campos):
        if(j == 0):
            cantEstados = arrMedicamentos[j] ** (j+1)
        else:
            cantEstados = (arrBase[j-1] * arrMedicamentos[j])
        arrBase.append(cantEstados)
    return arrBase

#Generar una lista de verdad por variable dada, segun sea su base
def generateArrayState(i,filas, maxEstado, cantVecesPorEstado,arrBase,arrMedicamentos):
    medlist = []
    contEstado = 1
    contadorInterno = 1
    for j in range (0, filas):                 
        medlist.append(maxEstado)
        #Lo hago para verificar si ya se puede reiniciar con el mayor
        if(contadorInterno == arrBase[i]):
            maxEstado = arrMedicamentos[i] -1
            contadorInterno = 1
            contEstado = 1
        else:
            #Lo hago para verificar si ya se puede cambiar de estado a un menor
            if(contEstado == cantVecesPorEstado):
                maxEstado -=1
                contEstado =1
            else:
                contEstado +=1
            contadorInterno +=1
    return medlist

#Devuelve la cantidad de filas de la tabla de verdad
def lenLista(arrMedicamentos):    
    filas = 1 
    for i in arrMedicamentos:
        filas  *=i 
    return filas

#Genera dataframe de las posibles combinaciones de medicamentos
def createCombinationDataframe(arrEstadosMedicamentos,arrHeader):
    campos = len(arrEstadosMedicamentos)
    dicc = {}
    i = 0   
    #Obteniendo el tamanho de la base por cada posicion del medicamento
    arrBase = generateArrayBase(campos, arrEstadosMedicamentos)   
    #Obteniendo la cantidad de filas
    filas   = lenLista(arrEstadosMedicamentos)         
    #Obteniendo la gran tabla de verdad
    for key in arrHeader:  
        maxEstado = arrEstadosMedicamentos[i] -1 
        cantVecesPorEstado = arrBase[i]/arrEstadosMedicamentos[i]
        medlist = generateArrayState(i,filas, maxEstado,cantVecesPorEstado,arrBase,arrEstadosMedicamentos)
        dicc[key] = medlist
        i+=1   
    df = pd.DataFrame(dicc)
    return df

#Genera dataframe de las probabilidades de presencia o ausencia del sintoma
def initProbabilisticDataframe(arrEstadosMedicamentos, listEstados, listProbabilidades):
    df= pd.DataFrame(columns=listEstados)
    filas   = lenLista(arrEstadosMedicamentos)
    for i in range(filas):
        df2 = pd.DataFrame([listProbabilidades], columns=listEstados)
        df =df.append(df2)
    df=df.reset_index(drop=True)    
    return df

#Devuelve el arreglo de los indices correspondientes a cada combinacion de medicamentos
def getIndexToSet(arrayIndex):
    arr = []
    for i in arrayIndex:
        arr.append(i)
    return arr

#Colocar la probabilidad dado la efectividad y evidencia para un medicamento
def setProbabilisticValue(dfSintoma, diccEstados, indexMedicamentos):
    new_df = pd.DataFrame(diccEstados, index=indexMedicamentos)
    dfSintoma.update(new_df)
    
#Obtener la lista de valores de la evidencia    
def getListaEvidencia(df):
    lista =[]
    # Iteración por filas del DataFrame:
    for indice_fila, fila in df.transpose().iterrows():
        lista.append(fila.values.tolist())
    return lista

