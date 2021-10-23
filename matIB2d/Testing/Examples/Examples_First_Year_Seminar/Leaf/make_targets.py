#! /usr/bin/python
# quick script to write the target point file for my experiment

import sys

def output_string(num, strength):
    s = "{}\n".format(num)
    for i in range(1,num+1):
        s = s + "{} {}\n".format(i, strength)
    return s

def main():
    if len(sys.argv)==1:
        print("Need file")
    elif len(sys.argv)==2:
        fname = sys.argv[1]
        expname = fname.split(".")[0]
        targetname = expname + ".target"
        print("Output name: {}.".format(targetname))
        with open(fname, 'r') as f:
            vertexnum = int(f.readline())
        print("{} vertices in file.".format(vertexnum))
        out = output_string(vertexnum, 1e7)
        with open(targetname, "w") as f:
            f.write(out)

if __name__ == '__main__': main()
