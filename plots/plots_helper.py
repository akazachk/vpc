import pandas as pd
pd.set_option("multi_sparse", True)

import numpy as np
import matplotlib.lines as mlines
from matrix2latex import matrix2latex

import matplotlib.pyplot as plt
scale=2
DPI = 200
plt.rc('text', usetex=True)
plt.rc('font', family='serif')
plt.rc('axes.spines', **{'bottom':True, 'left':True, 'right':False, 'top':False})
plt.rc('axes', titlesize=12*scale)
plt.rc('axes', labelsize=8*scale)
plt.rc('xtick', labelsize=8*scale)
plt.rc('ytick', labelsize=8*scale)
plt.rc("legend", fontsize=8*scale)
plt.rc("figure", figsize=[6*scale,4*scale])
# plt.rc("figure", figsize=[6,4])
#plt.rc("figure", figsize=[3,2])
plt.rc("savefig", dpi=DPI)

# For LaTeX helper functions
import re

# Some columns report both floats and ints
# This is a problem for siunitx that we need to fix explicitly
# We can check for any int values in the table and apply a format to all of them
from math import floor, ceil

def tex_escape(text):
    """
        :param text: a plain text message
        :return: the message escaped to appear correctly in LaTeX
    """
    conv = {
        '&': r'\&',
        '%': r'\%',
        '$': r'\$',
        '#': r'\#',
        '_': r'\_',
        '{': r'\{',
        '}': r'\}',
        '~': r'\textasciitilde{}',
        '^': r'\^{}',
        '\\': r'\textbackslash{}',
        '<': r'\textless{}',
        '>': r'\textgreater{}',
        '≥': r'$\ge$'
    }
    regex = re.compile('|'.join(re.escape(str(key)) for key in sorted(conv.keys(), key = lambda item: - len(item))))
    return regex.sub(lambda match: conv[match.group()], text)


def remove_presolved_from_name(name:str) -> str:
    """Remove _presolved from instance names"""
    return name.removesuffix("_presolved")


def create_multirow_string(strval: str, num_rows: int = 2, alignment: str = 'c', extra_format: str = ""):
    """
    Wrap \p strval in a multirow environment for a table.
    """
    return \
        "{" + \
        "\\multirow[" + alignment + "]{"+ str(num_rows) + "}{*}{" + \
        (extra_format + "{" if extra_format != "" else "") + \
        str(strval) + \
        ("}" if extra_format != "" else "") + \
        "}" + "}"


def format_col_as_multirow(curr_series: pd.core.series.Series):
    start_val = ''
    start_row = -1
    end_row = -1
    for val in curr_series:
        end_row += 1
        is_last_row = end_row == len(curr_series)-1
        if val != start_val or is_last_row:
            num_rows = (end_row - start_row) + is_last_row
            if start_row >= 0 and num_rows > 1:
                multirow_string = create_multirow_string(str(start_val), num_rows = num_rows)
                curr_series[start_row] = multirow_string
                if is_last_row:
                    curr_series[end_row] = ""
            start_row = end_row
            start_val = val
        else:
            curr_series[end_row] = ""


def is_val(val1: float, val2: float) -> bool:
    return abs(val1 - val2) < 1e-7


def is_int(val):
    """
    Checks whether given value should be treated as an int.

    Currently treats zero as a float always which is not ideal.
    """
    if isinstance(val, str) and val == '':
        return False
    try:
        floatval = float(val)
    except ValueError:
        # print("ValueError: ", val)
        return False
    # print("DEBUG:", val, ":", type(val))
    rounds_to_int = is_val(floatval, floor(floatval)) and is_val(floatval, ceil(floatval))
    is_zero = is_val(floatval, 0.0)
    # is_float_zero = (isinstance(val,str) and val.find('.') >= 0 and is_zero)
    return rounds_to_int and (not is_zero)


# def is_int_style(col : pd.core.series.Series):
#     # return ['background-color: green' if is_int(v) else '' for v in col]
#     return ['background-color: green' if is_int(v) else '' for v in col]


def int_format(val, num_digits = 3, add_phantom = False):
    if is_int(val):
        new_str = "{\\tablenum[table-format=" + str(num_digits) +  ".0]{" + str(val) + "}"
        new_str += "\\phantom{.00}" if add_phantom else ''
        new_str += "}"
        return new_str
    else:
        return val


def enclose_in_braces(val, val_to_match : str) -> str:
    if str(val)==str(val_to_match):
        return "{" + str(val) + "}"
    else:
        return str(val)


# The styler from pandas has some limitations, particularly no way to add \midrule at arbitrary places
def add_midrule(latex: str, index: int) -> str:
    """
    Adds a midrule either `index` lines after the start or -index lines before the last line of the table

    Args:
        latex: latex table
        index: index of horizontal line insertion (in lines)
    """
    lines = latex.splitlines()
    if index >= 0:
        lines.insert(index, r'\midrule')
    else:
        index_from_bottom = -index
        lines.insert(len(lines) - index_from_bottom - 2, r'\midrule')
    return '\n'.join(lines).replace('NaN', '')


# To add adjustbox, needs to be done after LaTeX string has been generated
def add_adjustbox_environment(latex: str) -> str:
    lines = latex.splitlines()
    start_env_ind = -1
    end_env_ind = -1
    curr_ind = -1
    for line in lines:
        curr_ind += 1
        if line.startswith(r"\begin{tabular}"):
            start_env_ind = curr_ind
        if line.startswith(r"\end{tabular}"):
            end_env_ind = curr_ind+1
    if (start_env_ind >= 0 and end_env_ind >= 0):
        lines.insert(start_env_ind, r'\begin{adjustbox}{width=1\textwidth}')
        lines.insert(end_env_ind+1, r'\end{adjustbox}')
    return '\n'.join(lines).replace('NaN', '')


def add_sisetup(latex: str, table_format = "2.2") -> str:
    latex = \
"""
{
\sisetup{
    table-alignment-mode = format,
    table-number-alignment = center,
    table-format = """ + \
    table_format + ",\n" + \
        "}\n" + \
        latex + \
        "\n}"
    return latex