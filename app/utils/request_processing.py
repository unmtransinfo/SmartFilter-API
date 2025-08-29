from typing import Union
import re
from rdkit import Chem
from typing import List
from flask import abort
from rdkit import Chem
from typing import List, Dict

def int_check(
    request,
    var_name: str,
    lower_limit: Union[None, int] = None,
    upper_limit: Union[None, int] = None,
    default_val: Union[None, int] = None,
):
    n = request.args.get(var_name, type=int) or default_val
    try:
        n = int(n)
    except:
        return abort(400, f"Invalid {var_name} provided. Expected int but got: {n}")
    if lower_limit is not None and n < lower_limit:
        return abort(400, f"Error: {var_name} must be greater than {lower_limit}")
    if upper_limit is not None and n > upper_limit:
        return abort(400, f"Error: {var_name} must be less than {upper_limit}")
    return n

def process_smiles_input(request, param_name: str, limit: int):
    value_list = request.args.get(param_name, type=str)
    if not value_list:
        return abort(400, f"No {param_name} provided")
    value_list = value_list.split(",")
    if len(value_list) > limit:
        return abort(
            400,
            f"Provided list of {param_name} exceeded limit of {limit}. Please provide <= {limit} {param_name} at a time.",
        )
    return value_list

def process_smarts_input(request, param_name: str = "smarts") -> list[str]:
    """
    Fetches and validates one or more SMARTS strings from the request.

    Args:
        request: Flask request object.
        param_name: Name of the query parameter.

    Returns:
        A list of valid SMARTS strings.

    Raises:
        400 Bad Request for missing or invalid SMARTS.
    """
    smarts_values = request.args.getlist(param_name)

    if not smarts_values:
        return abort(400, f"No {param_name} string(s) provided")

    valid_smarts = []
    for s in smarts_values:
        if not s or Chem.MolFromSmarts(s) is None:
            return abort(400, f"Invalid {param_name} string: {s}")
        valid_smarts.append(s)

    return valid_smarts 

