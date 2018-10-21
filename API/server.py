from flask import Flask, jsonify
import psycopg2
# Import smtplib for the actual sending function
import smtplib
import random
import sys
import json
import requests

from genetic_algorithm import *

app = Flask(__name__)

post = [{
        'sistema_medico': 'Cardiologia',
        'antecedentes': { 'AN_BRADICARDIA':0,
                          'AN_DISFUNCION_VENTRI':1,
                          'AN_HIPONATREMIA':1,
                          'AN_HIPOPOTASEMIA':1,
                          'AN_EXTRENIMIENTO':0,
                          'AN_MEDIDA_PRESION':0
                          },
        'medicamentos': {'MR_VERAPAMILO':0,
                          'MR_AMLODIPINO':0,
                          'MR_CLORTALIDONA':0,
                          'MR_ENALAPRIL':0,
                          'MR_IRBESARTAN':0,
                          'MR_LOSARTAN':0,
                          'MR_BISOPROLOL':0,
                          'MR_CAPTOPRIL':0
                        },
        'diagnostico' : {'SA_INSUFICIENCIA_CARDIACA':0,
                         'SA_HIPERTENSION_ARTERIAL':1,
                         'SA_ANGINA_PECHO':0
                        }
        }]


@app.route("/")
def hello():
    return "<h1> "+ post +"</h1>"


if __name__ == '__main__':
    app.run(debug=True)

    
    