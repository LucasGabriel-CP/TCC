import gurobipy as gp
from pathlib import Path
import yaml

config = None
with open(Path("gurobi_conf/conf.yml"), "r", encoding="utf-8") as file:
                config = yaml.safe_load(file)

options = {
    "WLSACCESSID": config['WLSACCESSID'],
    "WLSSECRET":  config['WLSSECRET'],
    "LICENSEID": config['LICENSEID'],
}

with gp.Env(params=options) as env, gp.Model(env=env) as model:
    
    model.optimize()