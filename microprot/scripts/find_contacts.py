#!/usr/bin/env python

from optparse import OptionParser
from contacts import calc_distance, get_PDB_coordinates

__author__ = "Tomasz Kosciolek"
__version__ = "1.01b"
__last_update__ = "07/04/2016"


def find_contacts(out_fh, coords, seq_sep, contact_coff):
    # calc distances
    for i in range(len(coords)):
        actual = coords[i][0]
        for j in coords[i:]:
            if j[0] >= actual + int(seq_sep):
                distance = calc_distance(coords[i][1], j[1],
                                         coords[i][2], j[2],
                                         coords[i][3], j[3])
                if distance <= float(contact_coff):
                    out_fh.write('%s\t%s\t%s\n' % (str(actual),
                                 str(j[0]), str(distance)))
        actual = []
    out_fh.close()


def main():
    parser = OptionParser(usage="usage: %prog [OPTIONS] [pdb] [output]")
    parser.add_option("-a", dest='aminoacid', default='CB',
                      help="CA, CB etc. (Default: CB)")
    parser.add_option("-c", dest='cutoff', default='8',
                      help="Cutoff (Default: 8)")
    parser.add_option("-d", dest='seq_separation', default='5',
                      help="Minimum sequence separation (Default: 5)")
    (options, args) = parser.parse_args()

    if len(args) != 2:
        parser.error("Incorrect number of arguments!... -h for help")

    finp = open(args[0], 'r')
    foutput = open(args[1], 'w')

    # extract coords
    coords = get_PDB_coordinates(finp, options.aminoacid)

    find_contacts(foutput, coords, options.seq_separation, options.cutoff)

    print("Done.")

if __name__ == "__main__":
    main()
