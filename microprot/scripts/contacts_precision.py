#!/usr/bin/env python
# -*- coding: utf-8 -*-

from optparse import OptionParser
import numpy as np
from find_contacts import calc_distance

__author__ = "Tomasz Kosciolek"
__version__ = "1.01b"
__last_update__ = "07/04/2016"

r'''
WORK NOTES:
add verbose flag option
add option to calculate HA contacts (need to re-assess the method to calculate
distances) when printing output make sure only filename is included (no path)

!!!!!!!!!!!!!!!
something is wrong with sorting when PPV values > 10.0 !!!!!!!!!
!!!!!!!!!!!!!!!
'''


def calc_prec(tr, fa):
    return "{0:.2f}".format(float(tr)/(float(tr)+float(fa)))


def ReadContacts(inp_fp, topx, contype):
    infile = open(inp_fp, 'r')
    lines = infile.readlines()[0:]
    infile.close()

    # make sure it's a PFRMAT RR file
    if not lines[0].startswith('PFRMAT RR'):
        raise ValueError('Contacts file is not in PFRMAT (CASP) format!')

    if contype == "all":
        contype_coff = 4
        contype_max = 10e8
    elif contype == "lr":
        contype_coff = 23
        contype_max = 10e8
    elif contype == "sr":
        contype_coff = 4
        contype_max = 23

    contacts = []

    for i in lines:
        line = i.split()
        if i.startswith("PFRMAT"):
            continue
        if (abs(int(line[1]) - int(line[0])) > contype_coff) and \
           (abs(int(line[0]) - int(line[1])) < contype_max):
            contacts.append([line[0], line[1], line[4]])

    contacts = np.array(contacts)

    # sort contacts list by predicted probability
    contacts_sorted = contacts[np.argsort(contacts[:, 2])][::-1]

    topX_contacts = contacts_sorted[0:int(topx)]

    aacid1 = topX_contacts[:, 0:2].tolist()
    aacid2 = topX_contacts[:, 2].tolist()

    # return AA1 and AA2 values
    return aacid1, aacid2


def FindPDBContacts(con_fp, pdb_fp, out_fp, aacid, coff, mind, topl, contype):

    finput = open(pdb_fp, 'r')
    if out_fp is not None:
        foutput = open(out_fp, 'w')
    else:
        foutput = None

    lines = finput.readlines()
    finput.close()

    # starting true and false values
    true = 0
    false = 0

    # extract coords
    coords = []
    for i in lines:
        line = i.split()
        if line[0] == 'ATOM' and aacid == line[2] and line[3] != 'GLY':
            if i[21] != ' ':
                coords.append([int(line[5]),
                               float(line[6]),
                               float(line[7]),
                               float(line[8])])
            else:
                coords.append([int(line[4]),
                               float(line[5]),
                               float(line[6]),
                               float(line[7])])
        elif line[0] == 'ATOM' and line[2] == 'CA' and line[3] == 'GLY':
            if i[21] != ' ':
                coords.append([int(line[5]),
                               float(line[6]),
                               float(line[7]),
                               float(line[8])])
            else:
                coords.append([int(line[4]),
                               float(line[5]),
                               float(line[6]),
                               float(line[7])])

    # read contacts
    PredCons = []
    PredPPV = []
    ReadTopX = int(len(coords)/float(topl))

    # make sure it is a 2xN array
    PredCons, PredPPV = ReadContacts(con_fp, ReadTopX, contype)

    # calc distances
    for i in range(len(coords)):
        actual = coords[i][0]
        for j in coords[i:]:
            if j[0] >= actual + int(mind):
                distance = calc_distance(coords[i][1], j[1],
                                         coords[i][2], j[2],
                                         coords[i][3], j[3])
                if (distance <= float(coff)) and \
                   ([str(actual), str(j[0])] in PredCons):

                    # write output only if provided
                    if foutput:
                        foutput.write(str(actual) + "\t" + str(j[0]) + "\t" +
                                      str(PredPPV[PredCons.index([str(actual),
                                          str(j[0])])]) + ' TRUE\n')
                    true += 1

                elif (distance > float(coff)) and \
                     ([str(actual), str(j[0])] in PredCons):
                    if foutput:
                        foutput.write(str(actual) + "\t" + str(j[0]) + "\t" +
                                      str(PredPPV[PredCons.index([str(actual),
                                          str(j[0])])]) + ' FALSE\n')
                    false += 1
        actual = []
    if foutput:
        foutput.close()
    print(str(pdb_fp) + " Top-L/" + str(topl) + " precision: " +
          str(calc_prec(true, false)))


def main():
    parser = OptionParser(usage=("usage: %prog [OPTIONS] [psicov_file] "
                                 "[pdb] OPTIONAL:[output]"))
    parser.add_option("-a", dest='aminoacid', default='CB',
                      help="CA, CB etc. (Default: CB)")
    parser.add_option("-c", dest='cutoff', default='8',
                      help="Cutoff (Default: 8)")
    parser.add_option("-d", dest='mindist', default='5',
                      help="Minimum sequence separation (Default: 5)")
    parser.add_option("-l", dest='topL', default='1',
                      help="Top-L/x contacts (Default: 1)")
    parser.add_option("-t", dest='conType', default='all',
                      help="Type of contacts [all | sr | lr] (Default: all)")
    (options, args) = parser.parse_args()
    aminoacid = options.aminoacid
    cutoff = options.cutoff
    mindist = options.mindist
    topL = options.topL
    conType = options.conType

    if len(args) != 2:
        if len(args) != 3:
            parser.error("Incorrect number of arguments!... -h for help")
        else:
            outfile = args[2]
    else:
        outfile = None

    confile = args[0]
    infile = args[1]
    FindPDBContacts(confile, infile, outfile, aminoacid,
                    cutoff, mindist, topL, conType)

if __name__ == "__main__":
    main()
