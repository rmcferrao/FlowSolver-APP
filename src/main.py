from flask import Flask, redirect, url_for, render_template, request, session, flash, send_from_directory, after_this_request
from datetime import timedelta, datetime 
import tempfile
import os
from algorithms.simulationWrapper import simu_wrapper
from flask_sqlalchemy import SQLAlchemy
from flask_caching import Cache

__author__ = "RafaelFerrao"

config = {
    "DEBUG": True,          # some Flask specific configs
    "CACHE_TYPE": "simple", # Flask-Caching related configs
    "CACHE_DEFAULT_TIMEOUT": 20
}

APP_ROOT = os.path.dirname(os.path.abspath(__file__))
app = Flask(__name__)

app.config['SEND_FILE_MAX_AGE_DEFAULT'] = 0
app.config.from_mapping(config)
cache = Cache(app)

app.config['FENICS_TEMP_FOLDER'] =  tempfile.mkdtemp()
app.config['FENICS_TEMP_FOLDER_IMAGES'] = app.config['FENICS_TEMP_FOLDER'] + '/images'

app.permanent_session_lifetime = timedelta(seconds=10)

@app.route("/")
def home():
    return render_template("flowSolver.html")

@app.route("/index")
def index():
    return render_template("index.html")

@app.route("/flowSolver", methods=["POST", "GET"])
def flow_solver():
    if request.method == 'GET':
        return render_template('flowSolver.html')

    elif request.method == 'POST':
        
        img_file = request.files['image-file']
        binImagePath = os.path.join(app.config['FENICS_TEMP_FOLDER'] , 'binImage.png') 
        img_file.save(binImagePath)

        dict_inp = request.form
        for key, val in zip(dict_inp.keys(), dict_inp.values()):
            print("Key: {} -> Value: {} and Type {} ----------".format(key,val, type(val)))
  
        basic_inputs = [dict_inp['meshRefinement'], 
                        dict_inp['pixelFreeFlowValue'], 
                        dict_inp['flowEquation'], 
                        dict_inp['inletPressure'], 
                        dict_inp['outletPressure'], 
                        dict_inp['fluidViscosity']
                        ]

        basic_inputs.insert(0, binImagePath)

        try:
            if dict_inp['flowEquation'] == 'Brinkman':
                k = dict_inp['rockMatrixPerm']
                simu_wrapper(*basic_inputs, k = k)

            elif dict_inp['flowEquation'] == 'Stokes':
                refine_mesh = True
                simu_wrapper(*basic_inputs, del_cells = refine_mesh)

            return redirect(url_for('flowResults'))

        except EnvironmentError as err:
            print('OS error: {}'.format(err))
            return redirect(url_for('SimulationError', error=err))

@app.route("/SimulationError/<string:error>", methods=["GET"])
def SimulationError(error):
    return render_template('SimulationError.html', error=error)    
        
@app.route("/flowResults", methods=["POST", "GET"])
def flowResults():
    if request.method == "GET":
        imgs_name = os.listdir(app.config['FENICS_TEMP_FOLDER_IMAGES'])
        return render_template('flowResults.html', image_names=imgs_name)

@app.route("/flowResults/<path:filename>", methods=["POST", "GET"])
def send_image(filename):
    return send_from_directory(app.config['FENICS_TEMP_FOLDER_IMAGES'], filename)


@app.route("/flowResults/download", methods=["POST", "GET"])
def download_zip_results():
    zip_file = 'simulation_results.zip'
    return send_from_directory(app.config['FENICS_TEMP_FOLDER'], zip_file, as_attachment=True)

@app.after_request
def add_header(response):
    # response.cache_control.no_store = True
    response.headers['Cache-Control'] = 'no-store, no-cache, must-revalidate, post-check=0, pre-check=0, max-age=0'
    response.headers['Pragma'] = 'no-cache'
    response.headers['Expires'] = '-1'
    return response

@app.errorhandler(404)
def page_not_found(e):
    return render_template("Error404.html")


if __name__ == "__main__":
    app.run(debug=True, port=8000)