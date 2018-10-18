#! /usr/bin/env python3

import os

varset = []

class Data:
    def __init__(self):
        self.R_0    = []
        self.A      = []
        self.u      = []
        self.v      = []
        self.w      = []
        self.sigma  = 0
        self.n      = 0

    def __gt__(self, other):
        return self.n > other.n
        


def dump_res(data):
    with open('r_0.txt', 'w') as f:
        for i in data:
            print(i.n, *i.R_0, file=f)

    with open('A.txt', 'w') as f:
        for i in data:
            print(i.n, *i.A, file=f)

    with open('u.txt', 'w') as f:
        for i in data:
            print(i.n, *i.u, file=f)

    with open('v.txt', 'w') as f:
        for i in data:
            print(i.n, *i.v, file=f)

    with open('w.txt', 'w') as f:
        for i in data:
            print(i.n, *i.w, file=f)

    with open('sigma.txt', 'w') as f:
        for i in data:
            print(i.n, i.sigma, file=f)



def parse_result(name):
    res = Data()
    try: 
        f = open(name, 'r')
    except IOError as e:
        print("No file ", name)
    else:
        r_0line = f.readline()
        for i in range(3):
            res.R_0.append(float(r_0line.split()[i]))

        r_0line = f.readline()
        res.sigma = float(r_0line.split()[0])

        r_0line = f.readline()
        res.n = int(r_0line.split()[0])
        
        r_0line = f.readline()
        for i in range(3):
            res.u.append(float(r_0line.split()[i]))
        
        r_0line = f.readline()
        for i in range(3):
            res.v.append(float(r_0line.split()[i]))
        
        r_0line = f.readline()
        for i in range(3):
            res.w.append(float(r_0line.split()[i]))
        
        r_0line = f.readline()
        for i in range(3):
            res.A.append(float(r_0line.split()[i]))

    return res

def walk():
    res = []
    dirs = os.listdir('.')
    for d in dirs:
        name = d + "/unfresult.txt"
        data = parse_result(name)
        if data.sigma != 0:
            res.append(data)
    return sorted(res)
        

def main():
    res = walk()
    dump_res(res)


main()
        
