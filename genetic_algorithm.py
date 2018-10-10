from pgmpy.models import BayesianModel
from pgmpy.factors.discrete import TabularCPD
from pgmpy.inference import BeliefPropagation
from pgmpy.factors.discrete import DiscreteFactor

from utils import *
import pandas as pd
import numpy as np
import multiprocessing
import random

class IndividualPathology:

    """ Clase destinada al sistema patologico (cardiologia, endocrinologia, etc) """

    def __init__(self, systemPathology, diccBackground, diccMedicines,arrStates,arrSymtopm, arrWeigth, bayesianNetwork):
        self.systemPathology = systemPathology
        self.diccBackground ={}
        self.diccMedicines  ={}
        
        self.arrStates       = arrStates[:]
        self.arrSymtopm      = arrSymtopm[:]
        self.arrWeigth       = arrWeigth[:]
        self.bayesianNetwork = bayesianNetwork
        
        for index in diccBackground:  
            self.diccBackground[index]  = diccBackground[index]
        
        for index in diccMedicines:           
            self.diccMedicines[index]  = diccMedicines[index]
            
            
class Individual:
    """ Clase abstracta para individuos de un algoritmo evolutivo. """

    def __init__(self, chromosome):
        self.chromosome = chromosome

    def crossover(self, other):
        "Retorna un nuevo individuo cruzando self y other."
        raise NotImplementedError
        
    def mutate(self):
        "Cambia los valores de algunos genes."
        raise NotImplementedError
        
        
        
class IndividualMedicines(Individual):
    """ Clase que implementa el individuo en el problema de combinacion de medicamentos """

    def __init__(self, chromosome):
        self.chromosome = chromosome[:]
        self.fitness    = -1

    def crossover_onepoint(self, other):
        """ Retorna dos nuevos individuos del cruzamiento de un punto entre self y other """
        c = random.randrange(len(self.chromosome))
        ind1 = Individual_medicines(self.chromosome[:c] + other.chromosome[c:])
        ind2 = Individual_medicines(other.chromosome[:c] + self.chromosome[c:])
        return [ind1, ind2]   
    
    def crossover_uniform(self, other):
        chromosome1 = []
        chromosome2 = []
        """ Retorna dos nuevos individuos del cruzamiento uniforme entre self y other """
        for i in range(len(self.chromosome)):
            if random.uniform(0, 1) < 0.5:
                chromosome1.append(self.chromosome[i])
                chromosome2.append(other.chromosome[i])
            else:
                chromosome1.append(other.chromosome[i])
                chromosome2.append(self.chromosome[i])
        ind1 = Individual_medicines(chromosome1)
        ind2 = Individual_medicines(chromosome2)
        return [ind1, ind2]     

    def mutate_position(self, arrStates):
        """ Cambia aleatoriamente el estado del medicamento (prescribir, prescribir bajas condiciones, no prescribir) """
        for i,val in enumerate(self.chromosome):
            if (random.uniform(0, 1) < 0.25):
                new_value = changeState(val, arrStates[i])
                self.chromosome[i] = new_value
        return self.chromosome        

def fitnessFunction(params):
    """Retorna el fitness de un cromosoma en el problema de combinacion de medicamentos """
    index, pathologyInd, individuo = params

    bp = BeliefPropagation(pathologyInd.bayesianNetwork)
    
    mergeSymEff      = pathologyInd.arrSymtopm + pathologyInd.arrEffects    
    newDiccMedicines = setStatesChromosome(pathologyInd.diccMedicines, individuo) 
    
    evidence = createEvidence(pathologyInd.diccBackground,pathologyInd.diccMedicines)   
    phi      = bp.query(variables=mergeSymEff, evidence = evidence)

    sumAbsen = 0
    sumTotal = 0
    
    for i,val in enumerate(mergeSymEff):
        arrProb = phi[val].values
        large = len(arrProb)
        sumAbsen += arrProb[large-1]*pathologyInd.arrWeigth[i]
        sumTotal += 1*pathologyInd.arrWeigth[i]  
    fitness = sumAbsen / sumTotal
    individuo.fitness = fitness

    return [index, fitness]


def evaluate_population(pathologyInd, population, fitness_fn, pool_size):
    """ Evalua una poblacion de individuos con la funcion de fitness pasada """
    popsize    = len(population)

    param_arr = []
    for i in range(popsize):
        if population[i].fitness == -1:  
            param_arr.append([i, pathologyInd, population[i]])

    pool = multiprocessing.Pool(pool_size)
    pool_results = pool.map(fitness_fn, param_arr)
    
    for fitness_arr in pool_results:
        population[fitness_arr[0]].fitness = fitness_arr[1]
    pool.close()
    pool.join()


def select_parents_tourney(population,ts):
    """ Seleccion de los mejores individuos o conjunto de medicamentos """
    
    population_select = population.copy()
    iParent1 = chooseParent(population_select,ts)
    population_select.pop(iParent1)
    iParent2 = chooseParent(population_select,ts)    
      
    return (population[iParent1], population[iParent2])

def select_survivors(population, offspring_population, numsurvivors):
    
    """ Selecciona sobrevivientes para la sgte generacion, dada la poblacion actual y poblacion de hijos """
    next_population = []
    population.extend(offspring_population)
    isurvivors = sorted(range(len(population)), key=lambda i: population[i].fitness, reverse=True)[:numsurvivors]
    for i in range(numsurvivors): next_population.append(population[isurvivors[i]])
    return next_population


def genetic_algorithm(pathologyInd, population, fitness_fn, ngen=100, pmut=0.1,ts=2):
    """ Algoritmo Genetico 
    Evalua la poblacion inicial, selecciona las parejas para el cruzamiento, 
    realiza el cruzamiento, luego la mutacion, evalua la poblacion descendencia,
    selecciona la siguiente poblacion de individuos """
    
    popsize    = len(population)

    evaluate_population(pathologyInd,population, fitness_fn)
    ibest       = sorted(range(len(population)), key=lambda i: population[i].fitness, reverse=True)[:1]
    bestfitness = [population[ibest[0]].fitness]
    print("Poblacion inicial, best_fitness = {}".format(population[ibest[0]].fitness))
    for g in range(ngen):
        
        mating_pool = []
        for i in range(int(popsize/2)): mating_pool.append(select_parents_tourney(population,ts)) 

        offspring_population = []
        for i in range(len(mating_pool)): 
            offspring_population+= mating_pool[i][0].crossover_onepoint(mating_pool[i][1])
        
        for i in range(len(offspring_population)):
            if random.uniform(0, 1) < pmut: 
                offspring_population[i].chromosome = offspring_population[i].mutate_position(pathologyInd.arrStates)

        evaluate_population(pathologyInd, offspring_population, fitness_fn) 

        population = select_survivors(population, offspring_population, popsize)

        ibest = sorted(range(len(population)), key=lambda i: population[i].fitness, reverse=True)[:1]
        bestfitness.append(population[ibest[0]].fitness)
        print("generacion {}, best_fitness = {}".format(g, population[ibest[0]].fitness))
        print("Combinacion de medicamentos= {}".format(population[ibest[0]].chromosome))
    
    return population[ibest[0]], bestfitness  

def genetic_search_medicines(ts,pathologyName, background, medicines, states, symtopms, 
                             weights, effects, bayesianNetwork, fitness_fn, num_medicines=11, 
                             popsize=20, ngen=30, pmut=0.25):
    population = []

    for i in range(popsize):
        chromosome = [random.randint(0,1) for j in range(1,(num_medicines+1))]
        random.shuffle(chromosome)
        population.append( IndividualMedicines(chromosome) )
    
    pathologyInd = IndividualPathology(pathologyName, background, medicines,states, 
                                       symtopms, weights, effects, bayesianNetwork)

    return genetic_algorithm(pathologyInd, population, fitness_fn, ngen, pmut,ts)

    