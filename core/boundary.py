class MyException(Exception):
    pass

class Domain:
    def __init__(self,list):
        '''Domain([b,c]). BEM in [-oo,b]; FEM in [b,c]; BEM in [c,oo]'''
        if len(list)!=2:
            raise MyException("Wrong number of boundaries of domain. Expecting 4.")
        self.lBEM=list[0]
        self.FEM=[list[0],list[1]]
        self.FEM_width=list[1]-list[0]
        self.uBEM=list[1]

    def __call__(self,i):
        if i<0:
            return self[0]
        elif i<=self.Ninput:
            if i not in self.store:
                self.store[i]=self[0]+i*self.width
            return self.store[i]
        else:
            return self[1]

    def __getitem__(self,i):
        return self.FEM[i]

    def extend(self,Ninput):
        self.store={}
        self.Ninput=Ninput
        self.width=float(self.FEM_width)/(Ninput+1)
