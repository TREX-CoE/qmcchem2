#!/usr/bin/env python3

from mpi4py import MPI
import sys
import gzip
import random
import math
import subprocess

admc_exec = "/home/scemama/qmcchem/src/MAIN/admc"
n_walk_per_proc = 10

def start():
    return subprocess.Popen(
        [ admc_exec, sys.argv[1] ],
        stdin=subprocess.PIPE,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE)


def read(process,len_walk):
    line = process.stdout.readline().decode("utf-8").strip()
    walk_num = int(line)
    walkers = []
    print(walk_num)
    for k in range(walk_num):
        w = []
        for i in range(len_walk):
            line = process.stdout.readline().decode("utf-8").strip()
            w.append( line )
        w = '\n'.join(w)
        walkers.append(w)

    _, E, W = process.stdout.readline().decode("utf-8").split()
    return walkers, float(E), float(W)


def write(process, message):
    process.stdin.write(f"{message}\n".encode("utf-8"))
    process.stdin.flush()


def terminate(process):
    process.stdin.close()
    process.terminate()
    process.wait(timeout=0.2)

def print_energy(EnergyWeight, Energy2Weight, Weight, N):
    e  = EnergyWeight / Weight
    e2 = Energy2Weight / Weight
    err = math.sqrt(abs(e*e - e2) / max(1,(N-1)) )
    print("%f +/- %f"%(e, err))
    return err

def main():
    try:
      input_dir = sys.argv[1]
    except:
        print("syntax: argv[0] [FILE]")
        sys.exit(-1)

    # Pool of electron coordinates
    with gzip.open(input_dir+"/electrons/elec_coord_pool.gz","r") as f:
        data = f.read().decode("utf-8").split()

    len_walk = int(data[1])*int(data[2])
    icount = 0
    buffer = []
    walkers = []
    for d in data[4:]:
        buffer.append(d)
        icount += 1
        if (icount == len_walk):
            walkers.append(buffer)
            buffer = []
            icount = 0

    walkers = [ '\n'.join(x) for x in walkers ]
    do_loop = True

    EnergyWeight  = 0.
    Energy2Weight = 0.
    Weight        = 0.
    NSamples      = 0.

    # Start processes
    proc = start()
    while do_loop:

        # Once every 1000, shuffle the list of walkers
        if random.random() < 0.01:
            print("SHUFFLE")
            random.shuffle(walkers)

        # Pick new walkers
        new_coords = walkers[:n_walk_per_proc]
        walkers = walkers[n_walk_per_proc:]

        # Send new walkers to the process
        write(proc, '\n'.join(new_coords))

        # Fetch new walkers from the process
        new_coords, e_new, w_new = read(proc, len_walk)
        walkers += new_coords

        # Print energy
        ew = e_new * w_new
        EnergyWeight  += ew
        Energy2Weight += e_new * ew
        Weight        += w_new
        NSamples      += 1.
        print (len(walkers))
        err = print_energy(EnergyWeight, Energy2Weight, Weight, NSamples)

        if err < 1.e-3:
                do_loop = False

    terminate(proc)
    return




if __name__ == "__main__":
    main()
