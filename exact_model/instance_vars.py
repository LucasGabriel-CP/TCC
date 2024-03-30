from pathlib import Path

class InstanceVars:
    def __init__(self) -> None:
        self.T = None
        self.V = None
        self.L = None
        self.K = None
        self.graph = None
        self.centers = None
        self.Dt = None
        self.filename = Path()
        self.path = Path("C:/Users/lucas/Documents/TCC/data/lck_instances")

    def read_instance(self, _filename):
        self.filename = _filename
        self.instantiate()

    def instantiate(self):
        with open(self.path / self.filename, 'r') as fp:
            line = fp.readline()
            line = line.split()
            #read V,T,L,k
            self.V = int(line[0])
            self.T = int(line[1])
            self.L = int(line[2])
            self.K = int(line[3])
            fp.readline()

            #read d(i,j)
            self.graph = [[] for _ in range(self.V)]
            for i in range(self.V):
                line = fp.readline()
                line = line.split()
                for j in line:
                    self.graph[i].append(int(j))
            fp.readline()

            #read delta_l
            line = fp.readline()
            line = line.split()
            self.centers = []
            for i in line:
                self.centers.append(int(i))
            fp.readline()

            #read D_t
            self.Dt = [[] for _ in range(self.T)]
            for i in range(self.T):
                line = fp.readline()
                line = line.split()
                for j in line[1:]:
                    self.Dt[i].append(int(j))
