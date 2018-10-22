from flask import Flask, request, jsonify
from pgmpy.readwrite import BIFReader
from genetic_algorithm import *

app = Flask(__name__)


@app.route('/processmedication',methods=['POST'])
def processmedication():
    data = request.get_json()
    sistema = data['sistema_medico']
    antecedentes = data['antecedentes']
    medicamentos = data['medicamentos']
    estados_medicamentos = data['estados_medicamentos']
    diagnostico  = data['diagnostico']
    pesos_diagnostico = data['pesos_diagnostico']
    
    reader = BIFReader('cardiologia.bif')
    bn0 = reader.get_model()
    
    best_ind, bestfitness = genetic_search_medicines(3, sistema,antecedentes, medicamentos,estados_medicamentos, diagnostico, pesos_diagnostico, bn0, fitnessFunction)
    arr = print_medicamentos(best_ind.chromosome,medicamentos)

    return jsonify({'med_recomendados':arr})

if __name__ == '__main__':
    app.run(debug=True)

    
    