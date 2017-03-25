#!/usr/bin/env python
# vim: set ai ts=4 sw=4 sts=4 noet fileencoding=utf-8 ft=python

import os
import re
import errno
import configparser
from color import print_error, print_color

CONFIG = ['sim_settings', '~/sim_settings']

def isfloat(string):
    """Checks if a given string is a float"""
    try:
        float(string)
        return True
    except ValueError:
        return False

class Settings():
    """Class to store and process settings"""
    def __init__(self, settings=None):
        self.__settings = {
            'OUTPUT_PATH': '.',
            'MCGEN_DATA': 'mcgen',
            'GEANT_DATA': 'geant',
            'LOG_DATA': 'log',
            'A2_GEANT_PATH': '',
            'QSUB_BIN': 'qsub',
            'QSUB_MAIL': 'a',
            'QUEUE': 'dflt',
            'WALLTIME': '12:00:00',
            'PRIORITY': 0,
            'GENERATOR': 'Ant-pluto',
            'GENERATOR_PATH': '',
            'Emin': 1420,
            'Emax': 1580,
            'COCKTAIL_SETUP': '',
            'COCKTAIL_BINNING': 0,
            'GUN_THETA': '0 180',
            'GUN_OPENING': '',
            'GeantFlags': '',
            'AddFlags': '',
            }

        if settings and len(settings[0]) == 2:
            for setting in settings:
                if setting[0] not in self.__settings:
                    print_error('[ERROR] Setting "%s" not valid! Will be skipped.' % setting[0])
                    continue
                if not setting[1] and self.__settings[setting[0]]:
                    print_color('[INFO] No value given for "{0}". Default value "{1}" will '
                                'be used'.format(setting[0], self.__settings[setting[0]]), 'BLUE')
                    continue
                if setting[0] in ('PRIORITY', 'COCKTAIL_BINNING'):
                    if not setting[1].isdigit():
                        print_error('[ERROR] The given value for "%s" is not an integer, use '
                                    'default value %d instead.' % (setting[1], self.__settings[setting[1]]))
                        continue
                    self.__settings[setting[0]] = int(setting[1])
                    continue
                if setting[0] in ('Emin', 'Emax'):
                    if not isfloat(setting[1]):
                        print_error('[ERROR] The given energy "%s" is not a float! Use default '
                                    'value %.1f instead.' % (setting[1], self.__settings[setting[0]]))
                        continue
                    self.__settings[setting[0]] = float(setting[1])
                    continue
                self.__settings[setting[0]] = setting[1]
        elif settings and len(settings[0]) != 2:
            print_error('[ERROR] Given settings in the wrong format, '
                        'list with tuples (setting, value) expected.')

    def print(self):
        """Print all current settings"""
        for key, val in self.__settings.items():
            print("%s: %s" % (key, val))

    def export(self, file_name, force=False):
        """Export the current settings to the given file"""
        file_name = os.path.expanduser(file_name)
        path = os.path.dirname(file_name)
        if not path:
            path = os.getcwd()
            file_name = os.path.join(path, file_name)
        if not os.path.exists(path) and not force:
            print_error('[ERROR] The path "%s" for exporting the settings does not exist!' % path)
            return False
        elif not os.path.exists(path):
            print_color('[INFO] The path "%s" does not exist, it will be created now.' % path, 'BLUE')
            try:
                os.makedirs(path)
            except OSError as exception:
                if exception.errno == errno.EACCES:
                    print_error('[ERROR] You do not have the permission to create directories in "%s"'
                                % os.path.dirname(path))
                    return False
                elif exception.errno != errno.EEXIST:
                    raise
        if not os.access(path, os.W_OK):
            print_error('[ERROR] The given path "%s" os not writable!' % path)
            return False
        if os.path.isfile(file_name):
            if not force:
                print_error('[ERROR] The specified file "%s" already exists, aborting export.' % file_name)
                return False
            else:
                print_color('[INFO] The file "%s" already exists, it will be overwritten.' % file_name, 'BLUE')
        with open(file_name, 'w') as file:
            file.write('%s\n' % '[settings]')
            file.write('%s\n' % '# output path, current directory will be used if missing')
            file.write('%s = %s\n' % ('OUTPUT_PATH', self.__settings['OUTPUT_PATH']))
            file.write('%s\n' % '# directory relative to output path above to store MC generated files')
            file.write('%s = %s\n' % ('MCGEN_DATA', self.__settings['MCGEN_DATA']))
            file.write('%s\n' % '# relative path to output path above '
                       'where the Geant4 files should be stored')
            file.write('%s = %s\n' % ('GEANT_DATA', self.__settings['GEANT_DATA']))
            file.write('%s\n' % '# log directory relative to output path above, '
                       'log will be used if empty')
            file.write('%s = %s\n' % ('LOG_DATA', self.__settings['LOG_DATA']))
            file.write('%s\n' % '# path to the a2geant binaries, $PATH is used if empty')
            file.write('%s = %s\n\n' % ('A2_GEANT_PATH', self.__settings['A2_GEANT_PATH']))
            file.write('%s\n' % '# some default settings for jobs')
            file.write('%s = %s\n' % ('QSUB_BIN', self.__settings['QSUB_BIN']))
            file.write('%s\n' % '# mail when job aborts (a), begins (b), ends (e), no mails (n)')
            file.write('%s = %s\n' % ('QSUB_MAIL', self.__settings['QSUB_MAIL']))
            file.write('%s = %s\n' % ('QUEUE', self.__settings['QUEUE']))
            file.write('%s = %s\n' % ('WALLTIME', self.__settings['WALLTIME']))
            file.write('%s = %s\n\n' % ('PRIORITY', self.__settings['PRIORITY']))
            file.write('%s\n' % '# simulation specific settings')
            file.write('%s\n' % '# the MC generator which should be used, i.e.')
            file.write('%s\n' % '# Ant-pluto, Ant-cocktail, or Ant-mcgun')
            file.write('%s\n' % '# Note: If you specify a different MC generator, please')
            file.write('%s\n' % '#       make sure that its output is readable by Geant!')
            file.write('%s\n' % '#       Flags can be passed via AddFlags option.')
            file.write('%s = %s\n' % ('GENERATOR', self.__settings['GENERATOR']))
            file.write('%s\n' % '# leave the generator path blank to use the $PATH variable')
            file.write('%s = %s\n' % ('GENERATOR_PATH', self.__settings['GENERATOR_PATH']))
            file.write('%s\n' % '# minimum energy of the photon beam')
            file.write('%s: %s\n' % ('Emin', self.__settings['Emin']))
            file.write('%s\n' % '# maximum energy of the photon beam')
            file.write('%s: %s\n' % ('Emax', self.__settings['Emax']))
            file.write('%s\n' % '# Ant-cocktail specific settings')
            file.write('%s\n' % '# Define only ONE of the following two settings!')
            file.write('%s\n' % '# Setup which should be used, e.g. "Setup_2014_07_EPT_Prod"')
            file.write('%s: %s\n' % ('COCKTAIL_SETUP', self.__settings['COCKTAIL_SETUP']))
            file.write('%s\n' % '# Binning for the beam energy, min and max energy as defined above')
            file.write('%s: %s\n' % ('COCKTAIL_BINNING', self.__settings['COCKTAIL_BINNING']))
            file.write('%s\n' % '# Ant-mcgun specific settings, min and max energy used as defined above')
            file.write('%s\n' % '# covered theta range in degree in the format min_theta max_theta')
            file.write('%s: %s\n' % ('GUN_THETA', self.__settings['GUN_THETA']))
            file.write('%s\n' % '# opening angle between particles in degree')
            file.write('%s: %s\n' % ('GUN_OPENING', self.__settings['GUN_OPENING']))
            file.write('%s\n' % '# additional flags passed to runGeant (which calls a2geant), '
                       'for example regex to replace information in detector macro setup')
            file.write('%s\n' % "# like 's~^(/A2/det/setTargetLength).*~$1 5 cm~'")
            file.write('%s: %s\n\n' % ('GeantFlags', self.__settings['GeantFlags']))
            file.write('%s\n' % '# additional flags passed to the generator, for example --flatEbeam')
            file.write('%s: %s\n\n' % ('AddFlags', self.__settings['AddFlags']))
            file.write('%s\n' % '[channels]')
            file.write('%s\n' % '# channels which should be simulated with Ant-pluto, line has to start with ";", '
                       'given in the syntax used in Pluto (do not forget the recoil proton!), '
                       'the amount of files and the number of events per file')
            file.write('%s\n' % '#;"p pi0 [g g]" 10 100000')
            file.write('%s\n' % '# in case of Ant-cocktail, start the string with "Cocktail", '
                       'followed by the amount of files and the number of events per file')
            file.write('%s\n' % '#;"Cocktail" 100 10000')
            file.write('%s\n' % '# for the particle gun Ant-mcgun, start the string with "Gun: <particle_list>", '
                       'particle list separated by spaces, followed by the amount of files '
                       'and the number of events per file')
            file.write('%s\n' % '#;"Gun: g g" 100 10000')

        return True

    def get(self, key):
        """Return a specific setting if it exists"""
        if key in self.__settings:
            return self.__settings[key]
        else:
            raise KeyError('The requested key "%s" does not exist!' % key)

    def set(self, key, value):
        """Set a specific setting if key exists"""
        if key not in self.__settings:
            raise KeyError('The key "%s" you want to update does not exist!' % key)
        else:
            self.__settings[key] = value

def check_file(path, verbose=False):
    """Check if a file exists"""
    path = os.path.expanduser(path)
    if not os.path.isfile(path):
        if verbose:
            print_error("[ERROR] The file '%s' does not exist!" % path)
        return None
    else:
        return path

def read_config(config_file):
    """Read the given config file and return the parsed settings and channels"""
    config = configparser.ConfigParser()
    config.optionxform = str  # preserve case of the options
    config.read(config_file.name)
    settings = list(config.items('settings'))
    settings = Settings(settings)
    if not config.has_section('channels'):
        print_color('[WARNING] No channels specified in %s' % config_file.name, 'YELLOW')
        return settings, []
    lines = [line.rstrip('\n') for line in config_file.readlines() if not line.startswith('#') and line.split()]  # last part excludes empty lines
    channels = lines[lines.index('[channels]')+1:]
    if not channels:
        print_color('[WARNING] No channels specified in %s' % config_file.name, 'YELLOW')
        return settings, []
    delimiter = ' '
    pattern = re.compile(r'''((?:[^%s"']|"[^"]*"|'[^']*')+)''' % delimiter)  # split using delimiter outside of quotations [http://stackoverflow.com/a/2787064]
    chnl = []
    for channel in channels:
        option = pattern.split(channel.lstrip(';'))[1::2]
        if len(option) != 3:
            print_error('[ERROR] Wrong number of arguments for channel %s' % option[0])
            print('        This channel will be skipped')
            continue
        if option[1] is '0' or option[2] is '0':
            #print('Skip channel', format_channel(option[0], False))
            print('Skip channel', option[0])
            continue
        if not option[1].isdigit() or not option[2].isdigit():
            print_error('[ERROR] The amount of files and events to simulate have to be integers!')
            print('        Check your settings for channel %s' % option[0])
            continue
        chnl.append((option[0], int(option[1]), int(option[2])))
    return settings, chnl

def find_config(config_places=CONFIG, verbose=False):
    """Checks given locations for a sim_settings config file,
    return path of the config if found or None otherwise"""
    if isinstance(config_places, str):
        config_places = [config_places]
    elif not isinstance(config_places, list):
        if verbose:
            print_error('[ERROR] Expected argument to be a list or string!')
        return None
    config = None
    for conf in config_places:
        if verbose:
            print('Checking config file %s' % conf)
        config = check_file(conf, verbose)
        if config:
            if verbose:
                print('File found: %s' % config)
            break

    return config


def main():
    """Main method to test reading and exporting config"""
    config = find_config()
    if not config:
        print('No config found!')
        return

    channels = []
    settings = Settings()
    settings_export = 'default_settings'
    if not settings.export(settings_export):
        print('Creating example settings "%s" failed.' % settings_export)
    with open(config, 'r') as conf:
        settings, channels = read_config(conf)

    settings.print()
    if not channels:
        print('No channels found!')
        return
    print(channels)


if __name__ == "__main__":
    main()
