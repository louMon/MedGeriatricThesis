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

def createHeader(campos):
    arrH=[]
    for i in range(campos):
        arrH.append('med'+str(i))
    return arrH

def createDataframe(arrMedicamentos,arrHeader):
    campos = len(arrMedicamentos)
    filas = 1 
    dicc = {}
        
    #Obteniendo la base por cada posicion del medicamento
    arrBase   = generateArrayBase(campos, arrMedicamentos)
    
    #Obteniendo la cantidad de filas
    for i in arrMedicamentos:
        filas  *=i 
        
    #Obteniendo la gran tabla de verdad
    i = 0
    for key in arrHeader:  
        maxEstado = arrMedicamentos[i] -1 
        cantVecesPorEstado = arrBase[i]/arrMedicamentos[i]
        medlist = generateArrayState(i,filas, maxEstado,cantVecesPorEstado,arrBase,arrMedicamentos)
        dicc[key] = medlist
        i+=1   
    df = pd.DataFrame(dicc)
    return df
