# vim: set ai ts=4 sw=4 sts=4 noet fileencoding=utf-8 ft=python

'''
This module provides functions to color string in the terminal output
and provides a logging class with colored strings
'''

__version__ = '1.0'

# Colored output
#The background is set with 40 plus the number of the color, and the foreground with 30
BLACK, RED, GREEN, YELLOW, BLUE, MAGENTA, CYAN, WHITE = range(8)

#These are the sequences need to get colored ouput
RESET_SEQ = "\033[0m"
COLOR_SEQ = "\033[1;%dm"
BOLD_SEQ = "\033[1m"

COLORS = {
    'WARNING': YELLOW,
    'INFO': WHITE,
    'DEBUG': BLUE,
    'CRITICAL': YELLOW,
    'ERROR': RED,
    'RED': RED,
    'GREEN': GREEN,
    'YELLOW': YELLOW,
    'BLUE': BLUE,
    'MAGENTA': MAGENTA,
    'CYAN': CYAN,
    'WHITE': WHITE
}

def color_string(string, color):
    fmt = ''
    if color in COLORS:
        fmt = COLOR_SEQ % (30 + COLORS[color]) + string + RESET_SEQ
    elif color in range(8):
        fmt = COLOR_SEQ % (30 + color) + string + RESET_SEQ
    else:
        fmt = string
    return fmt

def bold_string(string):
    return BOLD_SEQ + string + RESET_SEQ

def print_color(string, color):
    print(color_string(string, color))

import sys

def print_error(string):
    print(color_string(string, RED), file=sys.stderr)

import logging

class ColoredFormatter(logging.Formatter):
    def __init__(self, *args, **kwargs):
        # can't do super(...) here because Formatter is an old school class
        logging.Formatter.__init__(self, *args, **kwargs)

    def format(self, record):
        levelname = record.levelname
        color     = COLOR_SEQ % (30 + COLORS[levelname])
        record.levelname = levelname.center(8)
        message   = logging.Formatter.format(self, record)
        message   = message.replace("$RESET", RESET_SEQ)\
                           .replace("$BOLD",  BOLD_SEQ)\
                           .replace("$COLOR", color)
        for k, v in COLORS.items():
            message = message.replace("$" + k,    COLOR_SEQ % (v+30))\
                             .replace("$BG" + k,  COLOR_SEQ % (v+40))\
                             .replace("$BG-" + k, COLOR_SEQ % (v+40))
        return message + RESET_SEQ

# Custom logger class with multiple destinations
class ColoredLogger(logging.Logger):
    FORMAT = "[%(asctime)s] [$COLOR%(levelname)s$RESET]  $BOLD%(message)s$RESET"
    DATEFORMAT = '%Y-%m-%d %H:%M:%S'
    def __init__(self, name):
        logging.Logger.__init__(self, name, logging.DEBUG)

        color_formatter = ColoredFormatter(fmt=self.FORMAT, datefmt=self.DATEFORMAT)

        console = logging.StreamHandler()
        console.setFormatter(color_formatter)

        self.addHandler(console)
        return

