import xml.etree.ElementTree as ET
import re
import os

import tokenize
from io import StringIO

# Define unit conversion to mm
units = {
    'mm': 1.0,
    'cm': 10.,
    'm': 1e3,
    'um': 1e-3,
    'nm': 1e-6,
}

constants = {}

included_files = set()

def replace_units_tokenized(expr):
    tokens = []
    g = tokenize.generate_tokens(StringIO(expr).readline)

    prev_tok = None
    prev_prev_tok = None
    for toknum, tokval, *_ in g:
        # Case: number followed by unit name (e.g. 1.6 cm)
        if (
            toknum == tokenize.NAME
            and prev_tok is not None and prev_tok[0] == tokenize.OP and prev_tok[1] == "*"
            and prev_prev_tok is not None and prev_prev_tok[0] == tokenize.NUMBER
            and tokval in units
        ):
            number_val = prev_prev_tok[1]
            tokens.pop()  # remove '*'
            tokens.pop()  # remove number
            tokens.append((tokenize.OP, "("))
            tokens.append((tokenize.NUMBER, number_val))
            tokens.append((tokenize.OP, "*"))
            tokens.append((tokenize.NAME, 'units'))
            tokens.append((tokenize.OP, "["))
            tokens.append((tokenize.STRING, f'"{tokval}"'))
            tokens.append((tokenize.OP, "]"))
            tokens.append((tokenize.OP, ")"))
            prev_tok = None
            prev_prev_tok = None
        # Case: 1.6cm (no *) — handled as before
        elif (
            toknum == tokenize.NAME
            and prev_tok is not None and prev_tok[0] == tokenize.NUMBER
            and (prev_tok[1][-1].isdigit() or prev_tok[1][-1] == '.')  # last character numeric-ish
            and tokval in units
        ):
            number_val = prev_tok[1]
            tokens.pop()  # remove number
            tokens.append((tokenize.OP, "("))
            tokens.append((tokenize.NUMBER, number_val))
            tokens.append((tokenize.OP, "*"))
            tokens.append((tokenize.NAME, 'units'))
            tokens.append((tokenize.OP, "["))
            tokens.append((tokenize.STRING, f'"{tokval}"'))
            tokens.append((tokenize.OP, "]"))
            tokens.append((tokenize.OP, ")"))
            prev_prev_tok = None
            prev_tok = None
        else:
            tokens.append((toknum, tokval))
            prev_prev_tok = prev_tok
            prev_tok = (toknum, tokval)

    return tokenize.untokenize(tokens)

def parse_xml_file(filename, basepath=""):
    full_path = os.path.join(basepath, filename)
    if full_path in included_files:
        return  # Avoid duplicate includes
    included_files.add(full_path)

    tree = ET.parse(full_path)
    root = tree.getroot()

    # First handle any includes
    for inc in root.findall(".//include"):
        ref_file = os.path.expandvars(inc.attrib.get("ref"))
        if ref_file:
            parse_xml_file(ref_file, basepath=os.path.dirname(full_path))

    # Then parse constants
    for const in root.findall(".//constant"):
        name = const.attrib['name']
        value = const.attrib['value']
        constants[name] = value

def evaluate_expr(expr, seen=None):
    if seen is None:
        seen = set()

    # Replace constant names using token-by-token replacement
    tokens = []
    g = tokenize.generate_tokens(StringIO(expr).readline)
    for toknum, tokval, *_ in g:
        if toknum == tokenize.NAME and tokval in constants:
            if tokval in seen:
                raise ValueError(f"Circular reference: {tokval}")
            sub_expr = str(evaluate_expr(constants[tokval], seen | {tokval}))
            tokens.append((tokenize.NUMBER, f"({sub_expr})"))
        else:
            tokens.append((toknum, tokval))

    expr_resolved = tokenize.untokenize(tokens)
    expr_resolved = replace_units_tokenized(expr_resolved)

    try:
        return eval(expr_resolved, {"__builtins__": None}, {"units": units})
    except Exception as e:
        raise ValueError(f"Failed to evaluate: '{expr_resolved}': {e}")

# Example usage
if __name__ == "__main__":

    ## IDEA
    parse_xml_file("/cvmfs/sw.hsf.org/key4hep/releases/2025-05-29/x86_64-almalinux9-gcc14.2.0-opt/k4geo/00-22-ubhvqv/share/k4geo/FCCee/IDEA/compact/IDEA_o1_v03/IDEA_o1_v03.xml")
    parse_xml_file("/cvmfs/sw.hsf.org/key4hep/releases/2025-05-29/x86_64-almalinux9-gcc14.2.0-opt/k4geo/00-22-ubhvqv/share/k4geo/FCCee/IDEA/compact/IDEA_o1_v02/IDEA_o1_v02.xml")
    targets = ["VTXIB_z1", "VTXIB_z2", "VTXIB_z3", "VTXIB_r1", "VTXIB_r2", "VTXIB_r3", "VTXD_z1", "VTXD_z2", "VTXD_z3", "VTXD_rmax1", "VTXD_rmax2", "VTXD_rmax3", "VTXD_rmax1", "VTXD_rmin1", "VTXD_rmin2", "VTXD_rmin3", "VTXOB_z1", "VTXOB_z2", "VTXOB_r1", "VTXOB_r2"]

    ## CLD
    parse_xml_file("/cvmfs/sw.hsf.org/key4hep/releases/2025-05-29/x86_64-almalinux9-gcc14.2.0-opt/k4geo/00-22-ubhvqv/share/k4geo/FCCee/CLD/compact/CLD_o2_v07/Vertex_o4_v07_smallBP.xml")
    targets = ["VertexBarrel_zmax", "VertexBarrel_r1", "VertexBarrel_r2", "VertexBarrel_r3", "VertexBarrel_Support_Thickness", "VertexBarrel_Sensitive_Thickness", "VertexBarrel_DoubleLayer_Gap"]




    for target in targets:
        value_in_mm = evaluate_expr(constants[target])
        print(f"{target} = {value_in_mm} mm")

