#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
import click


r'''
Functions:
==========

read_PDB_coordinates
--------------------
    input:
        inp_fh  :   file handle
        aacid   :   str
    output:
        coords  :   list (4d)

read_contact_predictions
------------------------
    input:
        inp_fh  :   file handle
        topX    :   int
        contype :   str
        min_sep :   int
    output:
        aacids  :   list (2d)
        ppv     :   list

find_PDB_contacts
-----------------
TODO: add option to calculate HA (heavy atom) contacts
      (need to re-assess the method to calculate distances)

    input:
        coords      :   list (4d)
        out_fh      :   file handle
        contact_coff:   int
        seq_sep     :   int
        con_type    :   str
    output:
        _contact    :   list (3d - aa1, aa2, dist)

contact_precision
-----------------
    input:
        coords      :   list (4d)
        pred_cons   :   list (2d)
        pred_ppv    :   list
        out_fh      :   file handle
        con_coff    :   int
        minsep      :   int
        verbose     :   bool
    output:
        contact_precision   :   float

_contacts
---------
command line interface
'''


def _calc_prec(tr, fa):
    return (float(tr)/(float(tr)+float(fa)))


def _calc_distance(xA, xB, yA, yB, zA, zB):
    return np.sqrt(np.power(xA-xB, 2) +
                   np.power(yA-yB, 2) +
                   np.power(zA-zB, 2))


def _topN_contacts(coords, topL):
    return int(len(coords)/float(topL))


def _contype(contact_type, sequence_sep):
    if contact_type == "all":
        con_min = sequence_sep
        con_max = 10e9
    elif contact_type == "lr":
        con_min = 24
        con_max = 10e9
    elif contact_type == "sr":
        con_min = sequence_sep
        con_max = 23
    else:
        raise ValueError('contact type should be either `all`, `sr` or `lr`')
    return con_min, con_max


def read_PDB_coordinates(inp_fh, aacid='CB'):
    '''
    Read PDB coordinates from ATOM recuds in a `inp_fh` file.
    Returns a 4D list of coordinates with residue number
    and (x, y, z) coordinates.
    '''
    coords = []
    lines = inp_fh.readlines()
    for i in lines:
        line = i.split()
        if line[0] == 'ATOM' and aacid == line[2] \
           and line[3] != 'GLY':
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
    return coords


def read_contact_predictions(inp_fh, topX=1e10, contype='all', min_sep=5):
    '''
    Read `topx` predicted contacts of type `contype` from `inp_fp`
    and return lists of interacting residues `aacids` and predicted PPVs `ppv`
    '''

    lines = inp_fh.readlines()[0:]

    # make sure it's a PFRMAT RR file
    if not lines[0].startswith('PFRMAT RR'):
        raise ValueError('Contacts file is not in PFRMAT (CASP) format!')

    contype_coff, contype_max = _contype(contype, min_sep)

    contacts = []

    for i in lines:
        line = i.split()
        if i.startswith("PFRMAT"):
            continue
        if (abs(int(line[1]) - int(line[0])) >= contype_coff) and \
           (abs(int(line[0]) - int(line[1])) <= contype_max):
            contacts.append([line[0], line[1], line[4]])

    contacts = np.array(contacts)

    # sort contacts list by predicted probability
    contacts_sorted = contacts[np.argsort(contacts[:, 2].astype(float))][::-1]

    topX_contacts = contacts_sorted[0:int(topX)]

    aacids = topX_contacts[:, 0:2].astype(int).tolist()
    ppv = topX_contacts[:, 2].astype(float).tolist()

    # return aminoacids (AA1 and AA2 values) and predicted PPV values
    return aacids, ppv


def find_PDB_contacts(coords, out_fh=None,
                      contact_coff=8, seq_sep=5, con_type='all'):
    '''
    for `coords` find contacts that are within `contact_coff`
    at least `seq_sep` residues apart.
    If `out_fh` is provided, write the resuls to file.
    If not, return a generator.
    '''
    con_min, con_max = _contype(con_type, seq_sep)

    # calc distances
    for i in range(len(coords)):
        actual = coords[i][0]
        for j in coords[i:]:
            if (j[0] >= actual + int(con_min)) and \
               (j[0] <= actual + int(con_max)):
                distance = _calc_distance(coords[i][1], j[1],
                                          coords[i][2], j[2],
                                          coords[i][3], j[3])
                if distance <= float(contact_coff):
                    _contact = (str(actual), str(j[0]), float(distance))
                    _con_fmt = '%s\t%s\t%.4f\n' % _contact
                    if out_fh:
                        out_fh.write(_con_fmt)
                    else:
                        yield _contact
        actual = []
    if out_fh:
        out_fh.close()


def contact_precision(coords, pred_cons, pred_ppv, out_fh=None, con_coff=8,
                      minsep=5, verbose=False):
    '''
    Calculates contact precision (float) for `coords` using `pred cons` list
    and `pred_ppv` predicted probability values.
    Optionally, the function can annotate the contact predictions with TRUE or
    FALSE values and write to file `out_fh` and/or print to StdOut
    (verbose flag)
    '''
    # starting true and false values
    true_con = 0
    false_con = 0

    # calc distances
    for i, AAs in enumerate(pred_cons):
        aa1, aa2 = AAs
        if abs(aa1-aa2) >= minsep:
                distance = _calc_distance(coords[aa1-1][1], coords[aa2-1][1],
                                          coords[aa1-1][2], coords[aa2-1][2],
                                          coords[aa1-1][3], coords[aa2-1][3])
                if distance <= float(con_coff):
                    # write output only if provided
                    output = ('%s\t%s\t%.4f\t%s' %
                              (aa1, aa2, pred_ppv[i], 'TRUE'))
                    if out_fh:
                        out_fh.write(output+'\n')
                    if verbose:
                        print(output)
                    true_con += 1
                elif distance > float(con_coff):
                    output = ('%s\t%s\t%.4f\t%s' %
                              (aa1, aa2, pred_ppv[i], 'FALSE'))
                    if out_fh:
                        out_fh.write(output+'\n')
                    if verbose:
                        print(output)
                    false_con += 1
    if out_fh:
        out_fh.close()
    return _calc_prec(true_con, false_con)


# RUN FROM COMMAND LINE
@click.command()
@click.option('--mode', default='precision',
              type=click.Choice(['precision', 'find']),
              help='Run mode: \n \
                    precision - calculate contact precision \n \
                    find - find contacts in a PDB file')
@click.option('--aminoacid', '-a', default='CB',
              help='Contact atom (CA, CB, etc.)')
@click.option('--cutoff', '-c', default=8,
              help='Contact distance cutoff')
@click.option('--minsep', '-s',  default=5,
              help='Minimum sequence separation')
@click.option('--topl', '-l', default=1,
              help='Top-L/x contacts')
@click.option('--contype', '-t', default='all',
              type=click.Choice(['all', 'sr', 'lr']),
              help='Type of contacts (all OR short-range OR long-range)')
@click.option('--verbose', '-v', is_flag=True, help='Verbose mode')
@click.argument('confile', nargs=1, type=click.Path(exists=True))
@click.argument('infile', nargs=1, type=click.Path(exists=True))
@click.argument('outfile', nargs=1, default=None, type=click.Path(),
                required=False)
def _contacts(mode, confile, infile, outfile, aminoacid,
              cutoff, minsep, topl, contype, verbose):

    if contype == 'all':
        contacts_type = 'all'
    elif contype == 'sr':
        contacts_type = 'short-range (<= 23 residues)'
    elif contype == 'lr':
        contacts_type = 'long-range (> 23 residues)'

    if outfile:
        out_cons_fh = open(outfile, 'w')
    else:
        out_cons_fh = None

    coords = read_PDB_coordinates(open(infile, 'r'), aminoacid)

    if mode == 'precision':
        # read contacts
        ReadTopN = _topN_contacts(coords, topl)

        # make sure it is a 2xN array
        pred_contacts, pred_PPVs = read_contact_predictions(open(confile, 'r'),
                                                            ReadTopN, contype,
                                                            minsep)

        pred_precision = contact_precision(coords, pred_contacts, pred_PPVs,
                                           out_cons_fh, cutoff, minsep,
                                           verbose)
        print('precision: %.4f for %s top-L/%i %s contacts' %
              (pred_precision, infile.split('/')[-1],
               int(topl), contacts_type))

    elif mode == 'find':
        for cons in find_PDB_contacts(coords, out_cons_fh,
                                      cutoff, minsep, contype):
            print('%s\t%s\t%.2f' % (str(cons[0]), str(cons[1]),
                                    float(cons[2])))


if __name__ == "__main__":
    _contacts()
