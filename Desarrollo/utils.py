import pandas as pd
import numpy as np
import random

#########################################################################################################Utils to genetic algorithm

#Se escoge al padre
def chooseParent(population,ts):
    popsize = len(population)    
    rivalParent = getArrayRivals(popsize,ts)   
    maxFitness = 0
    iParent=0
    for index in rivalParent:
        if (population[index].fitness > maxFitness):
            maxFitness = population[index].fitness 
            iParent = index      
    return iParent
    

#Obtiene el arreglo de rivales
def getArrayRivals(popsize,ts):
    rival = []
    for i in range(ts): 
        rival_index = random.randint(0,popsize-1)
        if rival_index not in rival: 
            rival.append(rival_index)
    return rival


#Setea el estado del cromosoma
def setStatesChromosome(diccMedicines, individuo):
    chromosome = individuo.chromosome
    n = len(chromosome)
    i = 0
    for key in diccMedicines:
        diccMedicines[key] = chromosome[i]
        i+=1
    return diccMedicines
        

#Agrega la evidencia en el diccionario de antecedentes
def createEvidence(dicc1, dicc2):
    evidence = dicc1.copy()
    evidence.update(dicc2)
    return evidence


#Cambia el alelo del gen
def changeState(val, numStates):
    arrStates = []
    for i in range(0,numStates):
        arrStates.append(i)
    arrStates.remove(val)
    new_value = random.choice(arrStates)
    return new_value
    
#######################################################################################################
##Utils to bayesian network


#Forma el single query para el dataframe (solo las opciones senhaladas como prendidas, las demas apagadas)
def formQuery(arrAll, arrChoose):
    query=""
    arrBit = []
    print(arrAll)
    for idx,val in enumerate(arrAll):
        print(arrChoose)
        if(arrAll[idx] in arrChoose):
            arrBit.append(1)
        else:
            arrBit.append(0)
        print(arrBit)
    for idx,val in enumerate(arrAll):
        if(idx == (len(arrAll) -1)):
            query+= str(arrAll[idx]) + '==' + str(arrBit[idx])
        else:
            query+= str(arrAll[idx]) + '==' + str(arrBit[idx])+ ' & '
    return query
    
#Forma el listado de probabilidades entre solo una variable 
def formSingleProbability(arrState, arrProbability):
    probability_string=""
    dicc = {}
    for idx,val in enumerate(arrState):
        dicc[val] = [arrProbability[idx]]
    return dicc

    

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

def formDoubleEvidence(arrVariableEvidence):
    arrDoubleEvidence = []
    for i,val_i in enumerate(arrVariableEvidence):
        for j, val_j in enumerate(arrVariableEvidence):
            if(val_i != val_j):
                if(j>=i):
                    aux = []
                    aux.append(val_i)
                    aux.append(val_j)
                    arrDoubleEvidence.append(aux)     
    return arrDoubleEvidence
        

def generateProbabilisticList(arrStateEvidence, arrVariableEvidence,arrProbabilisticEvidence, 
                               arrStateVarible, arrProbabilistic):
    
    df_ANT = createCombinationDataframe(arrStateEvidence, arrVariableEvidence) 
    df     = initProbabilisticDataframe(arrStateEvidence, arrStateVarible, arrProbabilistic)
    
    #single evidence (solo una variable prendida)
    for i,val in enumerate(arrVariableEvidence):
        query_string = formQuery(arrVariableEvidence, val)
        ef_evidence  = df_ANT.query(query_string)
        arr_index    = getIndexToSet(ef_evidence.index)
        dicc = formSingleProbability(arrStateVarible,arrProbabilisticEvidence[i])
        setProbabilisticValue(df, dicc, arr_index)
    
    arrDoubleVaribleEvidence = formDoubleEvidence(arrVariableEvidence)
    
    #double evidence (combinacion de dos variables prendidas)
    for i, val in enumerate(arrDoubleVaribleEvidence):
        query_string = formQuery(arrDoubleVaribleEvidence, val)
        print(query_string)
        ef_evidence  = df_ANT.query(query_string)
        arr_index    = getIndexToSet(ef_evidence.index)
        dicc = formSingleProbability(arrStateVarible,arrProbabilisticEvidence[i])
        setProbabilisticValue(df, dicc, arr_index)
    
    values=getListaEvidencia(df)
    return values
