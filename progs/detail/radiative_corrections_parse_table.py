#!/usr/bin/env python3

import re
import numpy as np

'''
The intention of this script is to parse the ascii files provided in the source of the paper arxiv:1711.11001
and dump them in a formatted way into a header file which can be included in the PDalitzCorrections plugin to
create TGraph2D objects from it. The correction factors are applied as a weight to eta and eta' Dalitz decays
'''

def read_file(path):
    decimal = re.compile(r"\d?\.\d*")
    with open(path) as f:
        corrections, corrs = [], {}
        for line in f:
            if line.startswith('#') and 'y = ' in line:
                y = decimal.search(line).group(0)
                y = float(y)
                corrs[y] = []
            elif not line.startswith('#') and line.strip():
                vals = line.split()
                x = float(vals[0])
                vals = [vals[i] for i in [3, 4, 8]]
                corr = sum(map(float, vals))
                corrs[y].append(x)
                corrections.append((x, y, corr))

    print('parsed y values:', corrs.keys())
    print('individual data points per y value:')
    for y in corrs.keys():
        print('{:.2f}: {:3d}'.format(y, len(corrs[y])))
        # only keep first and last element for x min/max
        corrs[y] = corrs[y][::len(corrs[y])-1]

    return corrections, corrs

def write_table(corrections, path):
    with open(path, 'w+') as f:
        for val in corrections:
            f.write('%g %g %g\n' % val)

def write_header(corrections, limits, path, name):
    with open(path, 'w+') as f:
        f.write('#ifndef __%s_CORRECTIONS__\n' % name.upper())
        f.write('#define __%s_CORRECTIONS__\n\n' % name.upper())
        f.write('#include <vector>\n\n')
        x, y, corr = zip(*corrections)
        print(limits.values())
        x_tuples = [val for tuple in limits.values() for val in tuple]
        f.write('namespace %s {\n' % name)
        f.write('    static std::vector<double> x = {%s};\n'
                % ', '.join(map(str, x)))
        f.write('    static std::vector<double> y = {%s};\n'
                % ', '.join(map(str, y)))
        f.write('    static std::vector<double> corr = {%s};\n\n'
                % ', '.join(map(str, corr)))
        f.write('    static std::vector<double> y_vals = {%s};\n'
                % ', '.join(map(str, limits.keys())))
        f.write('    static std::vector<double> x_tuples = {%s};\n}\n\n'
                % ', '.join(map(str, x_tuples)))
        f.write('#endif')
    print('C++ header written to file:', path)

def main():
    # eta -> e+ e- g
    corr, limits = read_file('/home/sascha/1711.11001/anc/eta_e')
    write_header(corr, limits, 'eta_dilepton_radiative_corrections.h', 'eta_ee')
    # eta -> mu+ mu- g
    corr, limits = read_file('/home/sascha/1711.11001/anc/eta_mu')
    write_header(corr, limits, 'eta_dimuon_radiative_corrections.h', 'eta_mumu')
    # eta' -> e+ e- g
    corr, limits = read_file('/home/sascha/1711.11001/anc/etap_e')
    write_header(corr, limits, 'etap_dilepton_radiative_corrections.h', 'etap_ee')
    # eta' -> mu+ mu- g
    corr, limits = read_file('/home/sascha/1711.11001/anc/etap_mu')
    write_header(corr, limits, 'etap_dimuon_radiative_corrections.h', 'etap_mumu')

if __name__ == '__main__':
    main()
