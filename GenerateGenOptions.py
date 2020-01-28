#!/usr/bin/env python2
#encoding: UTF-8

import getopt
import sys

def usage():
    print("The usage:");
    print("python GenerateGenOptions.py -n Nsim ");

def main():

# Default values of generator parameters
    Nsim = 10000
    Eb = 10.6
    tLim = -1.2
    EgMin = 4
    EgMax = 10.6
    Q2Cut = 0.02
    LUND = 0


    try:
        opts, args = getopt.getopt(sys.argv[1:], "n:e:t:lho", ["help", "output=", "Egmin=", "Egmax=", "q2Cut=", "LUND"])
    except getopt.GetoptError as err:
        # print help information and exit:
        print str(err)  # will print something like "option -a not recognized"
        usage()
        sys.exit(2)
    output = None
    verbose = False
    for o, a in opts:
        if o == "-n":
            Nsim = a
        elif o == "-e":
            Eb = a
        elif o == "-t":
            if(float(a) > 0):
                print("You provided tMax as " + a)
                print("tMax should be a negative number.")
                print("Exiting")
                sys.exit()               
            tLim = a
        elif o in ("--Egmin"):
            EgMin = a;
        elif o in ("--Egmax"):
            EgMax = a;
        elif o in ("--q2Cut"):
            Q2Cut = a;
        elif o in ("-l", "--LUND"):
            LUND = 1
        elif o in ("-h", "--help"):
            usage()
            sys.exit()
        else:
            assert False, "unhandled option"


    optFile = open("GenOptions.dat","w")
    
    optFile.write("Nsim     " + str(Nsim) + "\n");
    optFile.write("Eb       " + str(Eb) + "\n");
    optFile.write("tLim     " + str(tLim) + "\n");
    optFile.write("EgMin    " + str(EgMin) + "\n");
    optFile.write("EgMax    " + str(EgMax) + "\n");
    optFile.write("Q2Cut    " + str(Q2Cut) + "\n");
    optFile.write("LUND     " + str(LUND));



if __name__ == "__main__":
    main();
