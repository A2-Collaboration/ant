#!/usr/bin/env python3
# -*- coding: latin-1 -*-
# vim: set ai ts=4 sw=4 sts=4 noet fileencoding=utf-8 ft=python

'''
This Python script is intended to run the simulation chain using
a specified MC generator like Ant-pluto, Ant-cocktail, or Ant-mcgun,
and the corresponding a2geant package on blaster.
Therefore a configuration file or command arguments are used
in order to create the jobs and submit them to the specified queue.
As a fallback an arbitrary generator together with AddFlags can be used.
'''

import os, sys
import re
import errno
import argparse
import datetime
import subprocess
import tempfile
from contextlib import contextmanager
from shutil import rmtree
from os.path import abspath, dirname, join as pjoin
from distutils.spawn import find_executable
from math import ceil
# import own modules (parse Pluto string, config file, provide colored output)
from color import print_color, print_error
import parse_pluto_string
from parse_config import Settings, find_config, read_config


@contextmanager
def cd(newdir, cleanup=lambda: True):
    """cd method to change directory as a decorator.
    due to context manager usage, directory is reverted
    even after an Exception is thrown"""
    prevdir = os.getcwd()
    os.chdir(os.path.expanduser(newdir))
    try:
        yield
    finally:
        os.chdir(prevdir)
        cleanup()

@contextmanager
def tempdir():
    """use decorator and context manager to create a temporary directory
    cd to temp directory and yield path, clean up is Exception safe"""
    dirpath = tempfile.mkdtemp()
    def cleanup():
        """delete created temp directory"""
        rmtree(dirpath)
    with cd(dirpath, cleanup):
        yield dirpath

def unquote(string):
    """
    If a string has single or double quotes around it, remove them.
    Make sure the pair of quotes match.
    If a matching pair of quotes is not found, return the string unchanged.
    """
    if (string[0] == string[-1]) and string.startswith(("'", '"')):
        return string[1:-1]
    return string

def check_path(path, create=False, write=True):
    """Check if given path exists and is readable as well as writable if specified;
    if create is true and the path doesn't exist, the directories will be created if possible"""
    path = os.path.expanduser(path)
    exist = os.path.isdir(path)
    if not exist and create:
        print("Directory '%s' does not exist, it will be created now" % path)
        # try to create the directory; if it should exist for whatever reason,
        # ignore it, otherwise report the error
        try:
            os.makedirs(path)
        except OSError as exception:
            if exception.errno == errno.EACCES:
                print_error("[ERROR] You don't have the permission to create directories in '%s'" % dirname(path))
                return False
            elif exception.errno != errno.EEXIST:
                raise
        return True
    elif not exist:
        print_error("[ERROR] Directory '%s' does not exist" % path)
        return False
    else:
        if not is_readable(path):
            print_error("[ERROR] Directory '%s' is not readable" % path)
            return False
        if write and not is_writable(path):
            print_error("[ERROR] Directory '%s' is not writable" % path)
            return False
    return True

def check_file(path, file):
    """Check if a file in a certain path exists and is readable;
    if the file argument is omitted only path is checked, it's recommended to use
    check_path instead for this case"""
    path = os.path.expanduser(path)
    if file is None:
        if not os.path.isfile(path):
            print_error("[ERROR] The file '%s' does not exist!" % (path))
            return False
        else:
            if not is_readable(path):
                print_error("[ERROR] The file '%s' is not readable!" % (path))
                return False
            return True
    path = get_path(path, file)
    if not os.path.isfile(path):
        print_error("[ERROR] The file '%s' does not exist!" % path)
        return False
    else:
        if not is_readable(path):
            print_error("[ERROR] The file '%s' is not readable!" % (path))
            return False
        return True

def check_permission(path, permission):
    """Check if the given permission is allowed for path"""
    if os.path.exists(path):
        return os.access(path, permission)
    else:
        return False

def is_readable(path):
    """Check if path is readable"""
    return check_permission(path, os.R_OK)

def is_writable(path):
    """Check if path is writable"""
    return check_permission(path, os.W_OK)

def is_executable(path):
    """Check if path is executable"""
    return check_permission(path, os.X_OK)

def get_path(path, file=None):
    """Get the absolute path of a given path and an optional file"""
    if not file:
        return abspath(os.path.expanduser(path))
    return abspath(os.path.expanduser(pjoin(path, file)))

def format_channel(channel, spaces=True):
    """Format a decay string using unicode symbols,
    if spaces is True only the initial and final states are shown"""
    replace = [
        ('dilepton', 'γ*'),
        ('dimoun', 'γ*'),
        ('pi', 'π'),
        ('etap', 'eta\''),
        ('eta', 'η'),
        ('mu', 'µ'),
        ('omega', 'ω'),
        ('rho', 'ρ'),
        ('phi', 'φ'),
        ('Sigma', 'Σ'),
        ('Delta', 'Δ'),
        ('Lambda', 'Λ'),
        ('Omega', 'Ω'),
        ('Xi', 'Ξ'),
        ('g', 'γ'),
        ('0', '⁰'),
        ('+', '⁺'),
        ('-', '⁻')
    ]

    for i, j in replace:
        channel = channel.replace(i, j)

    if spaces:
        chan = re.split(r'_', channel)
        chan = [chan[0], chan[-1]]  # remove intermediate states

        try:
            channel = "  {0:<4s} -->  {1}".format(*chan)
        except:
            channel = "  " + channel
    else:
        channel = channel.replace('_', ' --> ')

    return channel

def unit_prefix(number):
    """Format big number by using unit prefixes"""
    if number >= 1000000000:
        if str(number).count('0') >= 9:
            return re.sub(r"000000000$", "G", str(number))
        else:
            return str(number/1E9) + 'G'
    elif number >= 1000000:
        if str(number).count('0') >= 6:
            return re.sub(r"000000$", "M", str(number))
        else:
            return str(number/1E6) + 'M'
    elif number >= 1000:
        if str(number).count('0') >= 3:
            return re.sub(r"000$", "k", str(number))
        else:
            return str(number/1E3) + 'k'
    else:
        return str(number)

def timestamp():
    """Return a string with the current time stamp"""
    return str(datetime.datetime.now()).split('.')[0]

def max_file_number(lst):
    """Get the maximum file number from a list of files with the format *_[0-9]+.root"""
    if not lst:
        return 0
    num = re.compile(r'^.+_(\d+)\.root$')
    nrs = [int(num.search(l).group(1)) for l in lst if num.search(l) is not None]
    if not nrs:
        return 0
    else:
        nrs.sort()
        return nrs[-1]

def check_simulation_files(settings, channel):
    """Do some sanity checks if the existing simulation files seem to be okay
    and return the maximum file number"""
    mcgen_files = [f for f in os.listdir(settings.get('MCGEN_DATA'))
                   if f.startswith(settings.get('MCGEN_PREFIX')) and channel in f]
    geant_files = [f for f in os.listdir(settings.get('GEANT_DATA'))
                   if f.startswith(settings.get('GEANT_PREFIX')) and channel in f]
    max_mcgen = max_file_number(mcgen_files)
    max_geant = max_file_number(geant_files)
    if max_geant > max_mcgen:
        print_color("   [Warning]", 'YELLOW')
        print("There are more Geant4 simulation files than MC generated\nfiles for channel %s"
              % format_channel(channel, False))
        input("Will continue by pressing any key ")
    elif max_geant < max_mcgen:
        print_color("   [Warning]", 'YELLOW')
        print("There are more MC generated files than Geant4 simulated\nfiles for channel %s"
              % format_channel(channel, False))
        input("Will continue by pressing any key ")

    return max(max_mcgen, max_geant)

def list_file_amount(settings, events=False):
    """List the simulated amount of files per channel;
    if events is True and PyROOT can be used, the total amount of events will be determined, too"""
    mcgen_data = settings.get('MCGEN_DATA')
    geant_data = settings.get('GEANT_DATA')
    if not os.path.exists(mcgen_data):
        sys.exit('MC generated data path %s not found' % mcgen_data)
    if not os.path.exists(geant_data):
        sys.exit('Geant path %s not found' % geant_data)
    mcgen_files = [f for f in os.listdir(mcgen_data) if f.startswith(settings.get('MCGEN_PREFIX'))]
    geant_files = [f for f in os.listdir(geant_data) if f.startswith(settings.get('GEANT_PREFIX'))]

    # try to determine amount of Cocktail and Particle gun files
    cocktail_files = [f for f in mcgen_files if 'cocktail' in f.lower()]
    if cocktail_files:
        print('Cocktail: %d files' % len(cocktail_files))

    gun_files = [f for f in mcgen_files if 'gun' in f.lower()]
    if gun_files:
        print('Particle gun: %d files' % len(gun_files))

    # check if Pluto files are left
    pluto_files = set(mcgen_files) - set(cocktail_files) - set(gun_files)
    if not pluto_files:
        return

    # if Pluto files left, gather channels and determine the amount per channel
    channels = []
    decay = re.compile(r'^%s_(\S+)_\d+\.root$' % settings.get('MCGEN_PREFIX'))
    for name in pluto_files:
        channel = decay.search(name).group(1) if decay.search(name) else ''
        if channel and channel not in channels:
            channels.append(channel)

    print('Amount of simulated %s per channel:' % ('events' if events else 'files'))
    for channel in channels:
        pluto_channel = [f for f in pluto_files if channel in f]
        geant_channel = [f for f in geant_files if channel in f]
        max_pluto = max_file_number(pluto_channel)
        max_geant = max_file_number(geant_channel)
        maximum = max(max_pluto, max_geant)
        if maximum > 0:
            if not events:
                print(' {0:<20s} -- {1:>4d} files'.format(format_channel(channel), maximum))
            else:
                try:
                    from ROOT import TFile, TTree
                except ImportError:
                    print_error('[ERROR] PyROOT is not available, the total amount of '
                                'simulated events cannot be determined.')
                    return
                sum_events = 0
                #from ROOT import TFile, TTree, TH1
                for name in pluto_channel:
                    filename = get_path(mcgen_data, name)
                    current = TFile(filename)
                    if not current.IsOpen():
                        print_error("The file '%s' could not be opened" % filename)
                        continue
                    elif not current.GetListOfKeys().GetSize():
                        print_error("Found no directory in file '%s'" % current.GetName())
                        continue
                    name = current.GetListOfKeys().First().GetName()
                    tree = current.Get(name)
                    sum_events += tree.GetEntriesFast()
                print(' {0:<20s} -- {1:>4d} files,  total {2:>8s} events'.format(format_channel(channel), maximum, unit_prefix(sum_events)))

def check_directory(settings, value, force, verbose, relative=None, write=True):
    """Check if a given (relative) directory exists and is writable, update settings accordingly"""
    # check if the given output directory exists
    path = settings.get(value)
    if verbose:
        print('Check the specified "%s" directory %s' % (value, path))
    if relative is not None:
        path = get_path(relative, path)
    else:
        path = get_path(path)
    if not check_path(path, force, write):
        if verbose and not force:
            print('        Please make sure the specified directory exists.')
        return False
    if verbose:
        print('Path found:', path)

    settings.set(value, path)
    return True

def check_directories(settings, force=False, verbose=False):
    """Check if all the needed directories exist and are writable"""
    # check if the given output path exists
    if not check_directory(settings, 'OUTPUT_PATH', force, verbose):
        return False

    # check if the given log output path exists
    if not check_directory(settings, 'LOG_DATA', force, verbose, settings.get('OUTPUT_PATH')):
        return False

    # check if the given output path for MC generated events exists
    if (not settings.get('GEANT_ONLY') and
            not check_directory(settings, 'MCGEN_DATA', force, verbose, settings.get('OUTPUT_PATH'))):
        return False

    # check if the given geant output path exists
    if (not settings.get('MCGEN_ONLY') and
            not check_directory(settings, 'GEANT_DATA', force, verbose, settings.get('OUTPUT_PATH'))):
        return False

    # check if the given a2geant path exists
    if (not settings.get('MCGEN_ONLY') and
            settings.get('A2_GEANT_PATH') and
            not check_directory(settings, 'A2_GEANT_PATH', False, verbose, write=False)):
        return False

    # check if the given generator path exists
    if (not settings.get('GEANT_ONLY') and
            settings.get('GENERATOR_PATH') and
            not check_directory(settings, 'GENERATOR_PATH', False, verbose, write=False)):
        return False

    return True

def check_bin(path, file):
    """Check if a binary exists and is executable"""
    path = get_path(path, file)
    if not os.path.exists(path):
        print_error('[ERROR] The file "%s" does not exist!' % path)
        return False
    if not is_executable(path):
        print_error('[ERROR] The file "%s" is not executable!' % path)
        return False
    return path

def check_mcgen_bin(settings, generator_path='', verbose=False):
    """Check if the MC generator binary exists, return the absolute path"""
    generator = settings.get('GENERATOR')
    generator_path = settings.get('GENERATOR_PATH')
    if generator_path:
        if verbose:
            print('Searching for MC generator in %s' % generator_path)
        generator = check_bin(generator_path, generator)
        if generator and verbose:
            print('Found %s' % settings.get('GENERATOR'))
        if not generator:
            print_error('[ERROR] %s not found in %s!' % (settings.get('GENERATOR'), generator_path))
            sys.exit(1)
    else:
        generator = find_executable(generator)
        if not generator:
            print_error('[ERROR] %s not found!' % settings.get('GENERATOR'))
            if verbose:
                print("%s couldn't be found within your $PATH variable" % settings.get('GENERATOR'))
            sys.exit(1)
        else:
            generator = abspath(generator)
            if verbose:
                print('%s found: %s', (settings.get('GENERATOR'), generator))

    return generator

def check_geant_bin(settings, verbose=False):
    """Check if the Geant binary exists,
    do some simple checks related to the a2geant package,
    return the absolute path"""
    geant = None

    geant_path = settings.get('A2_GEANT_PATH')

    # if a custom Geant binary is given, skip all the checks related to the
    # conveniences provided by the Ant branch of the a2geant4 repository
    if settings.get('GEANT_BINARY'):
        if geant_path:
            if verbose:
                print('Searching for the Geant binary in %s' % geant_path)
            geant = check_bin(geant_path, settings.get('GEANT_BINARY'))
            if not geant:
                print_error('[ERROR] Geant executable "%s" not found!' % settings.get('GEANT_BINARY'))
                sys.exit(1)
            elif verbose:
                print('Geant executable "%s" found in %s' % (settings.get('GEANT_BINARY'), geant_path))
        else:
            geant = find_executable(settings.get('GEANT_BINARY'))
            if not geant:
                print_error("[ERROR] The Geant binary '%s' couldn't be found within your $PATH variable")
                if verbose:
                    print("No A2_GEANT_PATH defined in the config file, please provide a proper path or\n"
                            "or add the directory containing your Geant binary to your $PATH variable")
                sys.exit(1)
            else:
                if verbose:
                    print('Geant executable found:', geant)

        return geant

    # default way performing more checks regarding helper scripts and features/settings
    if geant_path:
        if verbose:
            print('Searching for the A2 binary in %s' % geant_path)
        if not check_bin(geant_path, 'A2'):
            print_error('[ERROR] A2 Geant executable not found!')
            if check_bin(geant_path, 'A2Geant4'):
                print_color("Please set 'GEANT_BINARY' in your config to 'A2Geant4' since you're using the new Geant version", 'YELLOW')
            sys.exit(1)
        elif verbose:
            print('A2 executable found in %s' % geant_path)
        if check_file(geant_path, 'runA2Geant'):
            geant = check_bin(geant_path, 'runA2Geant')
        else:
            geant = check_bin(geant_path, 'runGeant.sh')
        if not geant:
            print_error('[ERROR] The runA2Geant or runGeant.sh script could not be found or used!')
            sys.exit(1)
    else:
        geant = find_executable('A2')
        # fallback solution if not in $PATH
        if not geant:
            print_color("[WARNING] The A2 Geant binary couldn't be found within your $PATH variable", 'YELLOW')
            print('          try to use fallback solution')
            fallback = '/home/wagners/opt/a2geant/build/bin'
            if check_file(fallback, 'A2'):
                geant = check_bin(fallback, 'A2')
        if not geant:
            print_error('[ERROR] A2 executable not found!')
            if verbose:
                print("The A2 Geant binary couldn't be found neither in your $PATH variable nor in the fallback solution")
            sys.exit(1)
        else:
            if verbose:
                print('A2 executable found:', geant)
            geant_path = dirname(abspath(geant))
            if check_file(geant_path, 'runA2Geant'):
                geant = check_bin(geant_path, 'runA2Geant')
            else:
                geant = check_bin(geant_path, 'runGeant.sh')
            if not geant:
                print_error('[ERROR] The runA2Geant or runGeant.sh script could not be found or used!')
                sys.exit(1)

    # check if Geant version is used which can read in Pluto files (without pluto2mkin converter)
    if os.path.exists(get_path(geant_path, 'pluto2mkin')):
        print_error('[ERROR] pluto2mkin converter found in %s' % geant_path)
        print("        It's highly recommended to use the Ant branch of the a2geant repository.")
        sys.exit(1)

    # check target length in A2 Geant4 DetectorSetup.mac
    geant_base = geant_path.replace('/build/bin', '')
    geant_macros = get_path(geant_base, 'macros')
    if not check_file(geant_macros, 'DetectorSetup.mac'):
        print("        No 'DetectorSetup.mac' macro found in the Geant macros directory.")
    target_length = ''
    with open(get_path(geant_macros, 'DetectorSetup.mac'), 'r') as mac:
        for line in mac:
            if '/A2/det/setTargetLength' in line:
                target_length = line.split()[1]
    if float(target_length) < 10.:
        print_color("[WARNING] The target length specified in the 'DetectorSetup.mac' macro", 'YELLOW')
        print_color('          in the Geant macros directory is smaller than the usual lH2 target', 'YELLOW')
        print_color('          size: %s cm. If you consider to use a smeared z vertex, make sure' % target_length, 'YELLOW')
        print_color('          the specified target length is correctly set.', 'YELLOW')
        print()

    return geant

def check_binaries(settings, generator_path='', verbose=False):
    """Check if the needed binaries exist, return the absolute paths to them"""
    generator, geant = None, None

    # first of all check if the specified qsub binary exists
    if not find_executable(settings.get('QSUB_BIN')):
        print_error('[ERROR] The binary %s could not be found!' % settings.get('QSUB_BIN'))
        sys.exit(1)

    # skip the MC generator related checks if run in Geant-only mode
    # and set the generator to True (pass later checks for valid generator path)
    if not settings.get('GEANT_ONLY'):
        generator = check_mcgen_bin(settings, generator_path, verbose)
    else:
        generator = True

    # skip the Geant-related checks in case of MCgen-only mode
    if settings.get('MCGEN_ONLY'):
        return generator, True

    geant = check_geant_bin(settings, verbose)

    return generator, geant

def sanity_check_cocktail(settings):
    """Check if the given settings for Ant-cocktail seem okay"""
    if settings.get('GEANT_ONLY'):
        return True

    setup = settings.get('COCKTAIL_SETUP')
    binning = int(settings.get('COCKTAIL_BINNING'))

    if setup and not setup.startswith('Setup_'):
        print_color("[WARNING] The specified detector setup doesn't start with 'Setup_'", 'YELLOW')

    if not setup and not binning:
        print_error('[ERROR] No Setup and no binning for the beam energy specified!')
        print_error('        Please make sure to provide one option')
        return False
    if setup and binning:
        print_error('[ERROR] You provided both a setup and energy binning for the Cocktail')
        print_error('        Please make sure to provide only one of these options')
        return False

    return True

def sanity_check_mcgun(settings):
    """Check if the given settings for Ant-mcgun seem okay"""
    if settings.get('GEANT_ONLY'):
        return True

    theta_min, theta_max = settings.get('GUN_THETA').split()
    try:
        theta_min = float(theta_min)
    except ValueError:
        print_error('[ERROR] theta_min is not a float value!')
        return False
    try:
        theta_max = float(theta_max)
    except ValueError:
        print_error('[ERROR] theta_max is not a float value!')
        return False

    if theta_min > theta_max:
        print_error('[ERROR] theta_max is smaller than theta_min')
        return False

    if settings.get('GUN_OPENING'):
        try:
            float(settings.get('GUN_OPENING'))
        except ValueError:
            print_error('[ERROR] provided opening angle is not a float value!')
            return False

    return True

def sanity_check_channels(generator, channels):
    """check if the provided channel string match the specified generator"""
    if 'Ant-cocktail' in generator:
        wrong = [c[0] for c in channels if not c[0].lower().startswith('"cocktail"')]
        if wrong:
            print_error("[ERROR] The following channels don't match the requirements for Ant-cocktail")
            print(*wrong, sep='\n')
            print('The correct syntax is: "Cocktail" <#files> <#events>')
            return False
    if 'Ant-mcgun' in generator:
        wrong = [c[0] for c in channels if not c[0].lower().startswith('"gun:')]
        if wrong:
            print_error("[ERROR] The following channels don't match the requirements for Ant-mcgun")
            print(*wrong, sep='\n')
            print('The correct syntax is: "Gun: <space-separated particles>" <#files> <#events>')
            return False
    # for Ant-pluto case just check that no definitions for Cocktail and Gun are present
    if 'Ant-pluto' in generator:
        wrong = [c[0] for c in channels if (c[0].lower().startswith('"gun') or
                                            c[0].lower().startswith('"cocktail'))]
        if wrong:
            print_error("[ERROR] The following channels don't match the requirements for Ant-pluto")
            print(*wrong, sep='\n')
            print('The correct syntax for Ant-pluto is: "Pluto string" <#files> <#events>')
            if any('Gun' in l or 'gun' in l for l in wrong):
                print('Or you may want to use Ant-mcgun as your MC generator')
            if any('Cocktail' in l or 'cocktail' in l for l in wrong):
                print('Or you may want to use Ant-cocktail as your MC generator')
            return False

    return True

def input_digit(input_msg, max_retries=4, fail_msg='Invalid input, this channel will be skipped'):
    """Show the user an input dialogue to enter a digit, return 0 if it failed after max_retries"""
    num, count = 0, 0
    while True:
        count += 1
        try:
            num = int(input(input_msg + ' '))
            break
        except ValueError:
            if count < max_retries:
                print("Your input wasn't a number, please try again:")
            else:
                print(fail_msg)
                break
    return num

def simulation_dialogue():
    """Show a simulation dialogue in which the user can enter the channels to be simulated
    if they're not specified in the config file"""
    channels = []

    positive_responses = ['y', 'Y', 'j', 'J', 'yes', 'Yes']
    negative_responses = ['n', 'N', 'no', 'No']
    allowed_responses = positive_responses + negative_responses
    opt = input("\nDo you want to enter channels which should be simulated? [y/n]: ")
    while opt not in allowed_responses:
        opt = input("You've entered an invalid response! Please try again: ")

    if opt in negative_responses:
        return channels

    while opt in positive_responses:
        if not channels:
            print('Please enter a channel which should be simulated:')
            print('Note: The syntax has to be exactly the Pluto syntax for the reaction,\n'
                  '      e.g. "pi0 [g g]" for the decay of a pi0 into two photons.\n'
                  '      The recoil proton is taken into account, type only the desired reaction')
        else:
            opt = input("Do you want to enter another channel? [y/n]: ")
            while opt not in allowed_responses:
                opt = input("You've entered an invalid response! Please try again: ")
            if opt in negative_responses:
                continue

        chan = 'p '
        chan += input(chan)
        n_files = input_digit("How much files should be generated for this channel?", max_retries=3)
        if not n_files:
            print("This channel will be skipped")
            continue
        n_events = input_digit("How mich events should be generated per file?")
        if not n_events:
            print("This channel will be skipped")
            continue
        channels.append([chan, n_files, n_events])

    print_color("You've entered %d channels for the simulation process" % len(channels), 'BLUE')
    return channels

def get_simulation_files(source):
    """Obtain a file list from a given source
    the source could be a directory containing ROOT files or a file with a list of file paths"""
    files = []
    source = os.path.expanduser(source)
    if os.path.isdir(source):
        files = [get_path(abspath(source), f) for f in os.listdir(source) if
                 (os.path.isfile(pjoin(source, f)) and f.lower().endswith('.root'))]
    elif os.path.isfile(source):
        with open(source, 'r') as lst:
            files = [l.strip() for l in lst.readlines() if not l.startswith('#') and l.split()]
    else:
        print_error('[ERROR] Couldn\'t handle given input "%s"', source)
        sys.exit('No files for simulation found')

    return files

def get_decay_string(channel, level=1):
    """Get a decay string for a certain channel, print warning if proton is missing"""
    channel = channel.strip('"')

    # Ant-cocktail
    if channel.lower().startswith('cocktail'):
        return 'cocktail'

    # Ant-mcgun
    if channel.lower().startswith('gun:'):
        channel = channel.split(':')[-1].strip()
        channel = parse_pluto_string.get_initial_state(channel)
        channel = parse_pluto_string.particle_list_to_string(channel)
        return channel + '-gun'

    # default: assume Ant-pluto reaction channel syntax
    if channel.startswith('p '):
        channel = channel[2:]
    else:
        print_color('[WARNING] recoil proton missing in decay string: %s' % channel, 'WARNING')
    channel = parse_pluto_string.get_decay_string(channel, level)

    return channel

def get_file_name(prefix, channel, number, ext='root'):
    """Return the file string"""
    return '%s_%s_%04d.%s' % (prefix, channel, number, ext)

def create_sub(log_file, job_tag, job_number, settings):
    """Prepare the qsub command for the job submission"""
    qsub_cmd = settings.get('QSUB_BIN')
    qsub_cmd += ' -m ' + settings.get('QSUB_MAIL')
    qsub_cmd += " -M %s@kph.uni-mainz.de" % os.getenv('USER')
    qsub_cmd += " -N %s/%d" % (job_tag, job_number)
    qsub_cmd += " -j oe -o %s" % log_file
    qsub_cmd += " -z -q %s -V -p %d" % (settings.get('QUEUE'), settings.get('PRIORITY'))
    qsub_cmd += " -l ncpus=1,walltime=%s" % settings.get('WALLTIME')
    qsub_cmd += " %s" % settings.get('QSUB_EXTRA')

    return qsub_cmd

def create_mcgen_cmd(settings, generator, reaction, mcgen_file, events=1):
    """Prepare the MC generator command depending on the given generator"""
    emin = settings.get('Emin')
    emax = settings.get('Emax')
    setup = settings.get('COCKTAIL_SETUP')
    binning = settings.get('COCKTAIL_BINNING')
    addflags = settings.get('AddFlags')

    mcgen_cmd = generator
    if 'Ant-cocktail' in generator:
        flags = ''
        if setup:
            flags = '-s %s' % setup
        elif binning:
            flags = '--Emin %f --Emax %f -N %d' % (emin, emax, binning)
        mcgen_cmd += ' -o %s -n %d %s' % (mcgen_file, events, flags)
    elif 'Ant-mcgun' in generator:
        particles = reaction.replace('Gun:', '').strip('"')
        for particle in particles.split():
            mcgen_cmd += ' -p %s' % particle
        mcgen_cmd += ' -o %s -n %d --Emin %f --Emax %f' % (mcgen_file, events, emin, emax)
        if settings.get('GUN_THETA'):
            mcgen_cmd += ' --thetaMin {} --thetaMax {}'.format(*settings.get('GUN_THETA').split())
        if settings.get('GUN_OPENING'):
            mcgen_cmd += ' --openingAngle {}'.format(settings.get('GUN_OPENING'))
    elif 'Ant-pluto' in generator:
        mcgen_cmd += ' --reaction %s -o %s -n %d --Emin %f --Emax %f --no-bulk' \
                    % (reaction, mcgen_file, events, emin, emax)
    mcgen_cmd += ' ' + addflags

    return mcgen_cmd

def create_geant_cmd(settings, geant, mcgen_file, geant_file, test_job=False):
    """Create the runGeant command used to start a2geant"""
    # in case of MCgen-only mode, just return true to do nothing
    if settings.get('MCGEN_ONLY'):
        return 'true'

    geant_cmd = '%s %s %s' % (geant, mcgen_file, geant_file)
    if settings.get('GEANT_BINARY'):
        geant_cmd = '%s --if="%s" --of="%s"' % (geant, mcgen_file, geant_file)

    flags = settings.get('GeantFlags')
    # in case of a test job simulate only 1 event
    if test_job:
        flags += ' -n 1'
    if not flags:
        return geant_cmd

    # split geant flags, preserving spaces in between quotes
    delimiter = ' '
    pattern = re.compile(r'''((?:[^%s"']|"[^"]*"|'[^']*')+)''' % delimiter)
    flags = pattern.split(flags)[1::2]

    # check if geant flags contain regex
    regex = [unquote(flag) for flag in flags if unquote(flag).startswith('s')]
    for reg in regex:
        geant_cmd += " '%s' " % reg

    # remove regex from flags list
    flags = [flag for flag in flags if not any(flag for reg in regex if reg in flag)]
    geant_cmd += ' ' + ' '.join(flags)

    return geant_cmd

def test_process(cmd, time=None):
    """Try to run a process and check its return code
    if something went wrong, print command output
    if time is given, kill process after time expired"""
    # use shell=True, otherwise the command passed to Popen to execute,
    # including cmd.split(' ', 1), produces errors (probably due to reaction string)
    proc = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE,
                            stderr=subprocess.PIPE)

# Python 3.3 needed for timeout, most recent version on blaster is 3.1
#    if time:
#        if time < 0:
#            print_error('given time is negative!')
#            return False
#        try:
#            outs, errs = proc.communicate(timeout=time)
#            if proc.returncode:
#                print_error('non-zero returncode for command: %s' % cmd)
#                print_color('error output:', 'YELLOW')
#                print(errs.decode('utf-8'))
#                print_color('standard output:', 'YELLOW')
#                print(outs.decode('utf-8'))
#                return False
#            else:
#                return True
#        except TimeoutExpired:
#            proc.kill()
#            outs, errs = proc.communicate()
#            print_error('[ERROR] given timeout expired, process killed')
#            print_color('error output:', 'YELLOW')
#            print(errs.decode('utf-8'))
#            print_color('standard output:', 'YELLOW')
#            print(outs.decode('utf-8'))
#            return False

    outs, errs = proc.communicate()
    if proc.returncode:
        print_error('[ERROR] non-zero returncode for command: %s' % cmd)
        print_color('error output:', 'YELLOW')
        print(errs.decode('utf-8'))
        print_color('standard output:', 'YELLOW')
        print(outs.decode('utf-8'))
        return False

    return True

def run_test_job(settings, simulation, generator, geant):
    """Run a job with one event locally and check if it works"""
    first_job = next(iter(simulation or []), None)
    if not first_job:
        print_error('[ERROR] Unable to retrieve information of first job to be submitted')
        return False

    if settings.get('GEANT_ONLY'):
        filename = first_job
        out = os.path.basename(filename).rsplit('.', 1)[0]

        with tempdir() as tmp_path:
            geant_file = get_path(tmp_path, get_file_name(settings.get('GEANT_PREFIX'), out, 0))
            geant_cmd = create_geant_cmd(settings, geant, filename, geant_file, test_job=True)

            if not test_process(geant_cmd):
                return False

        return True

    decay_string, reaction, _, _, _ = first_job

    with tempdir() as tmp_path:
        mcgen_file = get_path(tmp_path, get_file_name(settings.get('MCGEN_PREFIX'), decay_string, 0))
        geant_file = get_path(tmp_path, get_file_name(settings.get('GEANT_PREFIX'), decay_string, 0))

        mcgen_cmd = create_mcgen_cmd(settings, generator, reaction, mcgen_file)
        geant_cmd = create_geant_cmd(settings, geant, mcgen_file, geant_file, test_job=True)

        if not test_process(mcgen_cmd):
            return False
        if not settings.get('MCGEN_ONLY') and not test_process(geant_cmd):
            return False

    return True

def submit_job(cmd, log_file, job_tag, job_number, settings):
    """Submit a job command"""
    qsub = create_sub(log_file, job_tag, job_number, settings)
    echo_args = ['echo', '-e', '%s' % cmd]
    proc = subprocess.Popen(echo_args, stdout=subprocess.PIPE)
    qsub_proc = subprocess.check_output(qsub.split(), stdin=proc.stdout)
#    if qsub_proc:
#        print('non-empty output of job-submission, maybe something went wrong...')
#        print(qsub_proc)
#        print()

def submit_jobs(settings, simulation, generator, geant, tag, total, verbose, length=20):
    """Create the MC generator and geant commands, submit the jobs, show progress bar"""
    # for progress bar
    bar = '='
    empty = ' '
    point = total/100
    increment = total/length
    job = 0

    mcgen_data = settings.get('MCGEN_DATA')
    geant_data = settings.get('GEANT_DATA')
    log_data = settings.get('LOG_DATA')
    time = timestamp()
    submit_log = 'submit_%s.log' % time.replace(' ', '_').replace(':', '.')[:-3]  # remove seconds
    submitted = ['Submitting %d jobs on %s\n\n' % (total, time)]

    if settings.get('GEANT_ONLY'):
        submitted.append('Gathered %d files to submit\n\n' % total)
        submitted.append('Used qsub command: %s\n\n' % create_sub('"%s_logfile"' % tag, tag, 42, settings))

        for i, filename in enumerate(simulation):
            job += 1
            out = os.path.basename(filename).rsplit('.', 1)[0]
            geant_file = get_path(geant_data, get_file_name(settings.get('GEANT_PREFIX'), out, i))
            log = get_path(log_data, get_file_name(tag, out, i, 'log'))
            geant_cmd = create_geant_cmd(settings, geant, filename, geant_file)
            job_cmd = ''
            if verbose:
                geant_echo = 'echo; echo %s; echo %s; echo' \
                        % ('Running Geant simulation with command:', geant_cmd.replace('"', '\"'))
                job_cmd = '%s; %s' % (geant_echo, geant_cmd)
            else:
                job_cmd = geant_cmd
            submit_job(job_cmd, log, tag, job, settings)
            submitted.append('%s\n' % geant_cmd)

            # progress bar
            sys.stdout.write('\r')
            fill = int(job/increment)
            # use ceil to ensure that 100% is reached, just in case of low precision
            sys.stdout.write("[%s%s] %3d%%" % (bar*fill, empty*(length-fill), ceil(job/point)))
            sys.stdout.flush()
        print()

        with open(get_path(settings.get('OUTPUT_PATH'), submit_log), 'w') as log:
            log.writelines(submitted)

        return

    total_events = 0
    for channel, _, files, events, _ in simulation:
        amount = files*events
        chnl = "  {0:<30s} {1:>4d} files per {2:>4s} events (total {3:>4s} events)\n" \
               .format(channel.replace('_', ' --> '), files, unit_prefix(events), unit_prefix(amount))
        submitted.append(chnl)
        total_events += amount
    submitted.append(" Total %s events in %d files\n\n" % (unit_prefix(total_events), total))
    submitted.append('\nUsed qsub command: %s\n\n' % create_sub('"%s_logfile"' % tag, tag, 42, settings))

    for decay_string, reaction, files, events, number in simulation:
        for i in range(1, files+1):
            job += 1
            mcgen_file = get_path(mcgen_data, get_file_name(settings.get('MCGEN_PREFIX'), decay_string, number+i))
            geant_file = get_path(geant_data, get_file_name(settings.get('GEANT_PREFIX'), decay_string, number+i))
            log = get_path(log_data, get_file_name(tag, decay_string, number+i, 'log'))
            mcgen_cmd = create_mcgen_cmd(settings, generator, reaction, mcgen_file, events)
            geant_cmd = create_geant_cmd(settings, geant, mcgen_file, geant_file)
            job_cmd = ''
            if verbose:
                mcgen_echo = 'echo; echo %s; echo %s; echo' \
                        % ('Running MC generation with command:', mcgen_cmd.replace('"', '\"'))
                geant_echo = 'echo; echo; echo %s; echo; echo %s; echo %s; echo' \
                        % ('-+'*50, 'Running Geant simulation with command:', geant_cmd.replace('"', '\"'))
                job_cmd = '%s; %s; %s; %s' \
                        % (mcgen_echo, mcgen_cmd, \
                        geant_echo, geant_cmd)
            else:
                job_cmd = '%s; %s' % (mcgen_cmd, geant_cmd)
            submit_job(job_cmd, log, tag, job, settings)
            submitted.append('%s; %s\n' % (mcgen_cmd, geant_cmd))

            # progress bar
            sys.stdout.write('\r')
            fill = int(job/increment)
            # use ceil to ensure that 100% is reached, just in case of low precision
            sys.stdout.write("[%s%s] %3d%%" % (bar*fill, empty*(length-fill), ceil(job/point)))
            sys.stdout.flush()
    print()

    with open(get_path(settings.get('OUTPUT_PATH'), submit_log), 'w') as log:
        log.writelines(submitted)

def is_valid_file(parser, arg):
    """Helper function for argparse to check if a file exists"""
    if not os.path.isfile(arg):
        parser.error('The file %s does not exist!' % arg)
    else:
        #return open(arg, 'r')
        return arg

def is_valid_dir(parser, arg):
    """Helper function for argparse to check if a directory exists"""
    if not os.path.isdir(os.path.expanduser(arg)):
        parser.error('The directory %s does not exist!' % arg)
    else:
        return os.path.expanduser(arg)


def main():
    """Main function: process all information, prepare simulation process, submit jobs"""
    parser = argparse.ArgumentParser(description='Submit simulations on blaster, configuration '
                                     'is done via a config file "sim_settings" in ~ or . '
                                     'and/or a list of command line arguments')
    parser.add_argument('-c', '--config', nargs=1, metavar='config_file',
                        dest='config',# required=True,
                        type=lambda x: is_valid_file(parser, x),
                        help='Optional: Specify a custom config file')
    parser.add_argument('-o', '--output', nargs=1, metavar='output_directory',
                        type=lambda x: is_valid_dir(parser, x),
                        help='Optional: Custom output directory')
    parser.add_argument('-g', '--generator', nargs=1, type=str, metavar='MC generator binary',
                        help='Optional: If the generator binary is not specified in the config '
                        'or should be overwritten, use this option to specify the binary manually')
    parser.add_argument('-p', '--path-generator', nargs=1, type=str, metavar='MC generator binary path',
                        help='Optional: If the generator binary (e.g. Ant-pluto) can\'t be found '
                        'in your $PATH variable, use this option to specify the path manually')
    parser.add_argument('--mcgen-only', dest='only_mcgen', action='store_true',
                        help='Only run the generation of MC files with the specified MC generator')
    parser.set_defaults(only_mcgen=False)
    parser.add_argument('--geant-only', dest='only_geant', action='store_true',
                        help='Only run the Geant simulation; Geant-readable input files required. '
                        'The files which should be simulated can be given either with the option -i '
                        'or a list of files can be piped to this submit script.')
    parser.set_defaults(only_geant=False)
    parser.add_argument('-i', '--input', nargs=1, metavar='Geant Input',
                        help='Optional: This option is only needed if the Geant-only mode is used. '
                        'The input can be either a list of files for Geant or a directory containing '
                        'ROOT files for the Geant simulation. Alternatively the list of files '
                        'can be piped to this script while omitting this option.')
    parser.add_argument('-l', '--list', nargs='?', const=True,
                        help='List the amount of existing files per channel '
                        'in the output directory and exit; if "all" is specified '
                        'as an optional argument the amount of events will be listed')
    parser.add_argument('-e', '--example-config', nargs='?', const='example_settings',
                        help='Export an example of default settings and exit. '
                        'If no file is specified, the output will be written to "example_settings"')
    parser.add_argument('-a', '--add-flags', nargs=1, type=str, metavar='"Additional flags"',
                        help='Optional: Define additional flags which will be passed to the '
                        'MC generator. Flags defined in the settings file will be overwritten. '
                        'All additional flags must be given as one quoted string!')
    parser.add_argument('-t', '--tag', nargs=1, type=str, metavar='Job Tag',
                        help='Optional: Specify a job tag, default is Sim')
    parser.add_argument('-w', '--walltime', type=int, nargs=1, metavar='walltime',
                        help='Walltime for jobs, time in hours')
    parser.add_argument('-q', '--queue', type=str, nargs=1, metavar='queue',
                        help='Queue which should be used for the jobs')
    parser.add_argument('-f', '--force', action='store_true',
                        help='Force creation of directories if they do not exist '
                        'or overwrite existing files')
    parser.add_argument('-v', '--verbose', action='store_true',
                        help='Print additional output')

    args = parser.parse_args()
    verbose = args.verbose
    force = args.force

    settings = None
    if args.example_config:
        print_color('[INFO] Write example settings to file "%s"' % args.example_config, 'BLUE')
        settings = Settings()
        if not settings.export(args.example_config, force):
            print_error('[ERROR] Creating example settings "%s" failed.' % args.example_config)
        sys.exit(0)

    config = None
    channels = []
    if args.config:
        config = args.config[0]
        if verbose:
            print('Use config file %s' % config)
    else:
        if verbose:
            print('Try to find a default config file')
        config = find_config(verbose=verbose)
        if not config:
            print_color('No config file found, use default values.', 'YELLOW')
            settings = Settings()
    if config:
        print_color('The config file "%s" will be used.' % config, 'BLUE')
        if verbose:
            print('Read config file %s' % config)
        with open(config, 'r') as conf:
            settings, channels = read_config(conf, args.only_geant)

    if args.list:
        if verbose and args.list == 'all':
            print('Trying to determine the total amount of simulated events')
        if args.output:
            settings.set('OUTPUT_PATH', get_path(args.output[0]))
        list_file_amount(settings, args.list == 'all')
        sys.exit(0)

    if verbose:
        print('The following default settings were determined:')
        settings.print()
        if channels:
            print('The following channels have been found:')
            print(channels)

    if args.output:
        if not check_path(args.output[0], force):
            sys.exit('The output directory %s cannot be used' % args.output[0])
        output = get_path(args.output[0])
        print_color('Setting custom output directory: %s' % output, 'GREEN')
        settings.set('OUTPUT_PATH', output)

    # check if the script is run in Geant-only or MCgen-only mode
    if args.only_mcgen and args.only_geant:
        print_error('[ERROR] Both Geant-only and MCgen-only mode are activated. Terminate...')
        sys.exit(1)
    if args.only_mcgen:
        print_color('Run in MCgen-only mode', 'BLUE')
    if args.only_geant:
        print_color('Run in Geant-only mode', 'BLUE')
    settings.set('MCGEN_ONLY', args.only_mcgen)
    settings.set('GEANT_ONLY', args.only_geant)

    if not args.only_geant and args.path_generator:
        path_generator = get_path(args.path_generator[0])
        print_color('Setting custom path for MC generator %s' % path_generator, 'GREEN')
        settings.set('GENERATOR_PATH', path_generator)

    if not check_directories(settings, force, verbose):
        sys.exit(1)
    settings.set('MCGEN_DATA', get_path(settings.get('OUTPUT_PATH'), settings.get('MCGEN_DATA')))
    settings.set('GEANT_DATA', get_path(settings.get('OUTPUT_PATH'), settings.get('GEANT_DATA')))
    settings.set('LOG_DATA', get_path(settings.get('OUTPUT_PATH'), settings.get('LOG_DATA')))

    if not args.only_geant and args.generator:
        generator = args.generator[0]
        print_color('Use custom MC generator %s' % generator, 'GREEN')
        settings.set('GENERATOR', generator)

    if not args.only_geant and args.add_flags:
        print_color('Set custom flags to pass to the MC generator: %s' % args.add_flags[0], 'GREEN')
        settings.set('AddFlags', args.add_flags[0])

    if not args.only_mcgen and settings.get('GEANT_BINARY'):
        print_color('Custom Geant binary specified in config file: %s' % settings.get('GEANT_BINARY'), 'GREEN')
        print('Please make sure to specify all needed flags for Geant in "GeantFlags"')
        print('Currently set Geant flags from config file:', settings.get('GeantFlags'))
        print('Cancel to modify the config with Ctrl+C')

    mc_generator, geant = check_binaries(settings, verbose)
    if not mc_generator or not geant:
        sys.exit(1)

    if not args.only_geant and not channels:
        print_color('[Warning] No channels specified in the config file', 'YELLOW')
        print_color('          Use -e to export example settings or -h for help', 'YELLOW')
        channels = simulation_dialogue()
        if not channels:
            sys.exit('No channels entered for simulation, exiting.')

    if args.walltime:
        print_color('Setting custom walltime to %d hours' % args.walltime[0], 'GREEN')
        settings.set('WALLTIME', '%d:00:00' % args.walltime[0])
    if args.queue:
        print_color('Setting custom queue to %s' % args.queue[0], 'GREEN')
        settings.set('QUEUE', args.queue[0])

    tag = args.tag[0] if args.tag else 'Sim'
    if verbose:
        print('Use the tag "%s" for the submitted jobs' % tag)

    if verbose:
        print('These are the updated settings:')
        settings.print()

    if not args.only_geant and 'Ant-cocktail' in mc_generator:
        if not sanity_check_cocktail(settings):
            sys.exit(1)
    if not args.only_geant and 'Ant-mcgun' in mc_generator:
        if not sanity_check_mcgun(settings):
            sys.exit(1)

    if not args.only_geant and not sanity_check_channels(mc_generator, channels):
        sys.exit(1)

    simulation = []
    if args.only_geant:
        files = []
        if verbose:
            print('Gather the files for the Geant simulation...')
        if args.input:
            if verbose:
                print('Scan %s for simulation files' % args.input[0])
            files = get_simulation_files(args.input[0])
        else:
            if verbose:
                print('Try to read input from stdin')
            if not sys.stdin.isatty():
                files = [l.strip() for l in sys.stdin.readlines() if not l.startswith('#')]
            else:
                print('No piped data found...')
                if verbose:
                    print('You may want to specify input with -i')
        if not files:
            print_error('[ERROR] No input files given for simulation')
            sys.exit(1)
        simulation = files
        if verbose:
            print('The following files have been found for simulation:')
            print(*simulation, sep='\n')
    else:
        if verbose:
            print('Determining the existing amount of files...')
        for channel in channels:
            decay_string = get_decay_string(channel[0])
            max_number = check_simulation_files(settings, decay_string)
            simulation.append([decay_string, channel[0], channel[1], channel[2], max_number])

    if args.only_geant:
        print_color(str(len(simulation)) + ' files found to run through Geant', 'BLUE')
        print(' Files will be stored in ' + settings.get('GEANT_DATA'))

        # run a test job for the first file to be submitted and check the output
        print_color('Test provided commands', 'BLUE')
        print(' Running first test job locally . . .')
        if not run_test_job(settings, simulation, mc_generator, geant):
            print_error('[ERROR] Test job failed, aborting job submission')
            if verbose:
                print('Please check your config file and make sure to provide absolute paths')
            sys.exit(1)
        print_color('Test job successful', 'BLUE')

        # start the job submission
        print('Start submitting jobs, total', len(simulation))
        submit_jobs(settings, simulation, mc_generator, geant, tag, len(simulation), verbose)
        print_color('Done!', 'GREEN')

        sys.exit(0)

    print_color(str(len(simulation)) + ' channels configured. '
                'The following simulation will take place:', 'BLUE')
    total_files, total_events = 0, 0
    for channel, _, files, events, _ in simulation:
        total = files*events
        chnl = channel
        if 'pluto' in mc_generator.lower():
            chnl = format_channel(channel)
        print("{0:<20s} {1:>4d} files per {2:>4s} events (total {3:>4s} events)"
              .format(chnl, files, unit_prefix(events), unit_prefix(total)))
        total_files += files
        total_events += total
    print(" Total %s events in %d files" % (unit_prefix(total_events), total_files))
    print(" Files will be stored in " + settings.get('OUTPUT_PATH'))

    # run a test job for the first command to be submitted and check the output
    print_color('Test provided commands', 'BLUE')
    print(' Running first test job locally . . .')
    if not run_test_job(settings, simulation, mc_generator, geant):
        print_error('[ERROR] Test job failed, aborting job submission')
        sys.exit(1)
    print_color('Test job successful', 'BLUE')

    # start the job submission
    print('Start submitting jobs, total', total_files)
    submit_jobs(settings, simulation, mc_generator, geant, tag, total_files, verbose)
    print_color('Done!', 'GREEN')


if __name__ == '__main__':
    try:
        main()
    except KeyboardInterrupt:
        print('\nCtrl+C detected, will abort simulation process')
        sys.exit(0)
    except Exception as e:
        print_error('An error occured during execution:')
        print(e)
        if '-v' or '--verbose' in sys.argv:
            raise e
        sys.exit(1)
