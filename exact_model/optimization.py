import gurobipy as gp
from pathlib import Path
import yaml
from instance_vars import InstanceVars

class Optmizer:
    def __init__(self) -> None:
        config = None
        with open(Path("gurobi_conf/conf.yml"), "r", encoding="utf-8") as file:
                        config = yaml.safe_load(file)

        self.options = {
            "WLSACCESSID": config['WLSACCESSID'],
            "WLSSECRET":  config['WLSSECRET'],
            "LICENSEID": config['LICENSEID'],
        }
        self.env = gp.Env(params=self.options)
        self.model = gp.Model(env=self.env)

        self.data = InstanceVars()
        self.X = None
        self.Y = None
        self.Z = None
        self.z = None
        self.gap = None
        self.runtime = None

    def add_data(self, filename):
        self.data.read_instance(filename)

    def define_vars(self):
        self.X = self.model.addVars(self.data.V, self.data.V, self.data.T, vtype=gp.GRB.BINARY, name='X')
        self.Y = self.model.addVars(self.data.V, self.data.L, self.data.T, vtype=gp.GRB.BINARY, name='Y')

    def constraint_1(self):
        """ Cliente deve ser atendido por soh uma facilidade. """
        for t in range(self.data.T):
            for j in self.data.Dt[t]:
                exp = gp.LinExpr()
                for i in range(self.data.V):
                    exp += self.X[i, j, t]
                self.model.addLConstr(exp == 1)

    def constraint_2(self):
        """ Numeros de facilidades simultaneas nao podem passar K. """
        for t in range(self.data.T):
            exp = gp.LinExpr()
            for i in range(self.data.V):
                for l in range(self.data.L):
                    p = max(t-self.data.centers[l]+1, 0)
                    for tt in range(p, t+1):
                        exp += self.Y[i, l, tt]
            self.model.addLConstr(exp <= self.data.K)
    
    def constraint_3(self):
        """ Todos clientes precisam ser atendidos. """
        for t in range(self.data.T):
            for i in range(self.data.V):
                for j in self.data.Dt[t]:
                    exp = gp.LinExpr()
                    for l in range(self.data.L):
                        p = max(t-self.data.centers[l]+1, 0)
                        for tt in range(p, t+1):
                            exp += self.Y[i, l, tt]
                    self.model.addLConstr(exp >= self.X[i, j, t])

    def constraint_4(self): # (NOVA)
        """ Nao pode ter mais de uma facilidade no mesmo lugar. """
        for t in range(self.data.T):
            for i in range(self.data.V):
                exp = gp.LinExpr()
                for l in range(self.data.L):
                    p = max(t-self.data.centers[l]+1, 0)
                    for tt in range(p, t+1):
                        exp += self.Y[i, l, tt]
                self.model.addLConstr(exp <= 1)

    def LKM(self):
        exp = gp.LinExpr()
        for t in range(self.data.T):
            for i in range(self.data.V):
                for j in self.data.Dt[t]:
                    exp += self.data.graph[i][j] * self.X[i, j, t]

        self.model.setObjective(
            exp,
            gp.GRB.MINIMIZE
        )

    def LKC(self):
        self.z = self.model.addVar(lb=0.0, vtype=gp.GRB.CONTINUOUS, name="obj")

        # Faazer Z ser a maior distancia
        for t in range(self.data.T):
            for i in range(self.data.V):
                for j in self.data.Dt[t]:
                    self.model.addLConstr(self.X[i,j,t]*self.data.graph[i][j] <= self.z)

        self.model.setObjective(self.z, gp.GRB.MINIMIZE)


    def lesgo(self, time_limit: float = 3600, problem: str = "LKC"):
        self.define_vars()
        self.constraint_1()
        self.constraint_2()
        self.constraint_3()
        self.constraint_4()
        if problem == "LKC":
            self.LKC()
        elif problem == "LKM":
            self.LKM()


        self.model.Params.TimeLimit = time_limit
        self.model.Params.LogToConsole = 0
        self.model.optimize()

        self.Z = self.model.ObjVal
        self.gap = self.model.MIPGap * 100
        self.runtime = self.model.Runtime
