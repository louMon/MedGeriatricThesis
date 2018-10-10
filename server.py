from flask import Flask, jsonify
import psycopg2

# Import smtplib for the actual sending function
import smtplib
import random
import sys
import json
import requests

app = Flask(__name__)
conn = psycopg2.connect(host="localhost",database="medicalSecure_development", user="pguser", password="1795lula")
cur = conn.cursor()