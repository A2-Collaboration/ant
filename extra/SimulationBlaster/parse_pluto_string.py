#!/usr/bin/env python3
# vim: set ai ts=4 sw=4 sts=4 noet fileencoding=utf-8 ft=python

import re
from itertools import groupby

def parenthetic_contents(string):
    """Generate parenthesized contents in string as pairs (level, contents)"""
    stack = []
    for i, char in enumerate(string):
        if char == '[':
            stack.append(i)
        elif char == ']' and stack:
            start = stack.pop()
            yield (len(stack), string[start + 1: i])

def sub_spaces(string):
    """Substitute multiple whitespaces with a single one plus remove trailing spaces
    and spaces after/before opening/closing brackets"""
    string = ' '.join(string.split())
    string = re.sub(r'(\S*)\s*(\])', r'\1\2', string)
    return re.sub(r'(\[)\s*(\S*)', r'\1\2', string)

def particle_list_to_string(plist, group=2):
    """Converts a given list with particles into a string
    same particles are grouped if greater than group=2, the second argument's default value"""
    grouped_list = [[i, len(list(j))] for i, j in groupby(sorted(plist))]
    string = ''
    for part, count in grouped_list:
        if count > group:
            string += '%d%s' % (count, part)
        else:
            string += part*count
    return string.replace("eta'", 'etap')

def get_initial_state(channel):
    """Returns a list containing the initial state particles"""
    lst = list(parenthetic_contents(channel))
    for l in lst:
        if l[0] == 0:
            channel = channel.replace(l[1], '')
    channel = channel.replace('[]', '')
    return re.split(r'\s+', channel.strip())

def get_final_intermediate_states(channel):
    """Parses the Pluto decay string and sorts particles into a list of final state particles
    as well as a dictionary containing intermediate states for each level
        (containing the FS of the corresponding level)"""
    channel = sub_spaces(channel)
    lst = list(parenthetic_contents(channel))
    max_idx = sorted(lst, key=lambda l: l[0])[-1][0]

    intermediate_states = {}
    final_state = []

    brackets = re.compile(r'\[(\s*)\]')
    empty_brackets = re.compile(r'(\[\s*\])')
    intermediate = re.compile(r'\[?(\S+)\s*\[\s*\]')
    intermediate_decay = re.compile(r'\[?(\S+\s*\[\s*\])')

    current, index = '', 0
    for i in range(0, max_idx):
        intermediate_states[i] = []
    while lst:
        iter_lst = lst[:]  # iterate over copy in order to remove elements from list
        if index == 0:
            current = ''
        for line in iter_lst:
            if not current and '[' not in line[1]:
                final_state += line[1].split()
                if line[0] < max_idx:
                    intermediate_states[line[0]] += line[1].split()
                current, index = sub_spaces(line[1]), line[0]
                lst.remove(line)
                continue
            if not current and '[' in line[1]:
                print('[ERROR] This should not happen! Broken regex?')
                return None
            if current and current in line[1] and index-1 == line[0]:
                index -= 1
                if re.search(r'\[.*' + current + r'.*\]', line[1]):
                    removed_current_state = re.sub(r'(\[.*)' + current + r'(.*\].*)', r'\1\2', line[1], 1)
                else:
                    removed_current_state = line[1].replace(current, '', 1)
                new_val = empty_brackets.sub('', removed_current_state)
                new_val = sub_spaces(new_val)
                intmed = intermediate.findall(removed_current_state)
                if not intmed:
                    lst[lst.index(line)] = (line[0], new_val)
                    continue

                intmed = intmed[0]
                new_val = intermediate_decay.sub('', removed_current_state, count=1)
                intermediate_states[line[0]].append(intmed)
                decay = intermediate_decay.findall(removed_current_state)[0]
                insert_pos = brackets.search(decay).start()+1
                current = decay[:insert_pos] + current + decay[insert_pos:]

                if current == line[1]:
                    lst.remove(line)
                else:
                    lst[lst.index(line)] = (line[0], new_val)

    return final_state, intermediate_states

def get_decay_string(channel, max_level=0):
    """Parses a Pluto string, identifies initial, final and intermediate states
    and returns a string with sorted particles in order to be reproducible.
    The level of intermediate states can be specified by the second argument,
    the default value is max_level=0 which means all intermediate states,
    speratated by underscores, and 1 means only the first intermediate state"""

    channel = sub_spaces(channel)

    final_state, intermediate_states = get_final_intermediate_states(channel)
    initial_state = get_initial_state(channel)
    im_states = []
    for level, intermediates in intermediate_states.items():
        im_states.append([level+1, particle_list_to_string(intermediates)])

    im_string = ''
    for level, string in im_states:
        if max_level and level > max_level:
            break
        im_string += '_' + string
    is_string = particle_list_to_string(initial_state)
    fs_string = particle_list_to_string(final_state)

    return is_string + im_string + '_' + fs_string


def main():
    """Main method for testing the functions"""
    #channel = "a1[b1 [c1 [ d1  d2] c2] b2  ]  a2 [b3[ c3 c4]]"
    #channel = "eta' [pi0 [g g] pi0 [dilepton [e+ e-] g] eta [pi0 [g g] pi0 [dilepton [e+ e-] g] pi0 [g g]]]"
    channel = "eta' [g rho0 [g pi0 [g g]]]"
    print('channel:', channel)

    decay_string = get_decay_string(channel, 1)
    print('final decay string:', decay_string)


if __name__ == "__main__":
    main()

