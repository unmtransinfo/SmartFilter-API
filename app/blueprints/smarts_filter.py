from flask import Blueprint, request, jsonify
from rdkit import Chem
from rdktools import smarts
from models import PainsRun, db
from utils.request_processing import process_smiles_input, process_smarts_input, process_multi_smarts_input
import tempfile

smarts_filter = Blueprint("smarts_filter", __name__, url_prefix="/smarts_filter")


import re

def parse_smiles_input(inputs, smiles_col=0, name_col=1, delim=None):
    parsed = []
    invalid = []

    splitter = re.compile(delim) if delim else re.compile(r"\s+")

    for line in inputs:
        parts = splitter.split(line.strip())
        if len(parts) <= smiles_col:
            invalid.append(line.strip())
            continue

        smiles = parts[smiles_col].strip()
        name = parts[name_col].strip() if len(parts) > name_col else smiles

        mol = Chem.MolFromSmiles(smiles)
        if mol:
            mol.SetProp("_Name", name)
            parsed.append((name, smiles, mol))
        else:
            invalid.append(smiles)

    return parsed, invalid


class MoleculeCollector:
    def __init__(self):
        self.accepted = []
        self.accepted_names = []

    def SetProps(self, props):
        pass

    def write(self, mol):
        canon = Chem.MolToSmiles(mol, canonical=True)
        name = mol.GetProp("_Name") if mol.HasProp("_Name") else ""
        self.accepted.append(canon)
        self.accepted_names.append((name, canon))


class MatchCountCollector:
    def __init__(self):
        self.results = []

    def SetProps(self, props):
        pass

    def write(self, mol):
        self.results.append({
            "smiles": Chem.MolToSmiles(mol, canonical=True),
            "name": mol.GetProp("_Name") if mol.HasProp("_Name") else "",
            "n_matches": mol.GetIntProp("n_matches")
        })


class MatchMultiCountCollector:
    def __init__(self):
        self.results = []

    def SetProps(self, props):
        self.props = props

    def write(self, mol):
        record = {
            "smiles": Chem.MolToSmiles(mol, canonical=True),
            "name": mol.GetProp("_Name") if mol.HasProp("_Name") else "",
        }
        for prop in self.props:
            if mol.HasProp(prop):
                record[prop] = mol.GetIntProp(prop)
            else:
                record[prop] = 0
        self.results.append(record)


class MatchFilter:
    def __init__(self):
        self.accepted = []

    def SetProps(self, props):
        pass

    def write(self, mol):
        self.accepted.append({
            "SMILES": Chem.MolToSmiles(mol, canonical=True),
            "name": mol.GetProp("_Name") if mol.HasProp("_Name") else "",
        })


@smarts_filter.route('/get_filterpains', methods=['GET'])
def get_filterpains():
    smiles_list = process_smiles_input(request, "SMILES", 1000)
    names_list = process_smiles_input(request, "Names", 1000)

    if not isinstance(smiles_list, list) or not smiles_list:
        return jsonify({"error": "No SMILES provided"}), 400

    if names_list and len(names_list) != len(smiles_list):
        return jsonify({
            "error": f"'SMILES' and 'Names' must be the same length: got {len(smiles_list)} and {len(names_list)}"
        }), 400

    if not names_list:
        names_list = smiles_list

    parsed = []
    invalid = []

    for smiles, name in zip(smiles_list, names_list):
        mol = Chem.MolFromSmiles(smiles)
        if mol:
            mol.SetProp("_Name", name)
            parsed.append((name, smiles, mol))
        else:
            invalid.append(name)

    if not parsed:
        return jsonify({"error": "No valid molecules found"}), 400

    # Collector for passing molecules
    collector = MoleculeCollector()
    print(parsed)
    mols = [mol for _, _, mol in parsed]
    exclude_mol_props = request.args.get("ExcludeMolProp", type=bool, default=False)

    smarts.FilterPAINS(mols, collector, exclude_mol_props)
    accepted_set = set(collector.accepted)

    passed = []
    failed = []

    for name, smiles, mol in parsed:
        if Chem.MolToSmiles(mol, canonical=True) in accepted_set:
            passed.append({"name": name, "smiles": smiles})
        else:
            failed.append({"name": name, "smiles": smiles})

    for inv in invalid:
        failed.append({"name": inv, "smiles": ""})

    record = PainsRun(inputs=[{"smiles": s} for s in smiles_list], passed=passed, failed=failed)
    db.session.add(record)
    db.session.commit()

    return jsonify({
        "passed": passed,
        "failed": failed
    }), 200



@smarts_filter.route('/get_matchcounts', methods=['GET'])
def get_matchcounts():
    smiles_list = process_smiles_input(request, "SMILES", 1000)
    names_list = process_smiles_input(request, "Names", 1000)

    if not isinstance(smiles_list, list) or not smiles_list:
        return jsonify({"error": "No SMILES provided"}), 400

    if names_list and len(names_list) != len(smiles_list):
        return jsonify({"error": f"'SMILES' and 'Names' must be the same length: got {len(smiles_list)} and {len(names_list)}"}), 400

    if not names_list:
        names_list = smiles_list

    smart = process_smarts_input(request)[0]
    if not smart:
        return jsonify({"error": "Smarts pattern is required"}), 400

    exclude_mol_props = request.args.get("ExcludeMolProp", type=bool, default=False)
    usa = request.args.get("usa", type=bool, default=False)

    parsed = []
    for smiles, name in zip(smiles_list, names_list):
        mol = Chem.MolFromSmiles(smiles)
        if mol:
            mol.SetProp("_Name", name)
            parsed.append(mol)

    if not parsed:
        return jsonify({"error": "No valid molecules parsed from input"}), 400

    collector = MatchCountCollector()

    smarts.MatchCounts(
        smarts=smart,
        usa=usa,
        molReader=parsed,
        molWriter=collector,
        exclude_mol_props=exclude_mol_props,
    )

    return jsonify(collector.results), 200


@smarts_filter.route('/get_matchfilter', methods=['GET'])
def get_matchfilter():
    smiles_list = process_smiles_input(request, "SMILES", 1000)
    names_list = process_smiles_input(request, "Names", 1000)

    if not isinstance(smiles_list, list) or not smiles_list:
        return jsonify({"error": "No SMILES provided"}), 400

    if not isinstance(names_list, list) or not names_list:
        return jsonify({"error": "No Names provided"}), 400

    if len(smiles_list) != len(names_list):
        return jsonify({"error": f"'SMILES' and 'Names' must be the same length: got {len(smiles_list)} and {len(names_list)}"}), 400

    smart = process_smarts_input(request)[0]
    if not smart:
        return jsonify({"error": "Smarts pattern is required"}), 400

    exclude_mol_props = request.args.get("ExcludeMolProp", type=bool, default=False)

    parsed = []
    for smiles, name in zip(smiles_list, names_list):
        mol = Chem.MolFromSmiles(smiles)
        if mol:
            mol.SetProp("_Name", name)
            parsed.append(mol)

    if not parsed:
        return jsonify({"error": "No valid molecules parsed from input"}), 400

    collector = MatchFilter()
    smarts.MatchFilter(
        smarts=smart,
        molReader=parsed,
        molWriter=collector,
        exclude_mol_props=exclude_mol_props,
    )

    return jsonify(collector.accepted), 200


@smarts_filter.route('/get_multi_matchcounts', methods=['GET'])
def get_multi_matchcounts():
    smiles_list = process_smiles_input(request, "SMILES", 1000)
    names_list = process_smiles_input(request, "Names", 1000)

    if not isinstance(smiles_list, list) or not smiles_list:
        return jsonify({"error": "No SMILES provided"}), 400

    if names_list and len(names_list) != len(smiles_list):
        return jsonify({"error": f"'SMILES' and 'Names' must be the same length: got {len(smiles_list)} and {len(names_list)}"}), 400

    if not names_list:
        names_list = smiles_list

    try:
        smarts_list = process_multi_smarts_input(request, "smarts")
    except Exception as e:
        return jsonify({"error": str(e)}), 400

    if not smarts_list:
        return jsonify({"error": "No valid SMARTS patterns provided"}), 400

    with tempfile.NamedTemporaryFile(mode="w+", delete=False) as temp_smarts_file:
        for sma in smarts_list:
            temp_smarts_file.write(sma.strip() + "\n")
        temp_smarts_file_path = temp_smarts_file.name

    exclude_mol_props = request.args.get("ExcludeMolProp", type=bool, default=False)
    usa = request.args.get("usa", type=bool, default=False)
    raiseError = request.args.get("strict", type=bool, default=False)

    parsed = []
    for smiles, name in zip(smiles_list, names_list):
        mol = Chem.MolFromSmiles(smiles)
        if mol:
            mol.SetProp("_Name", name)
            parsed.append(mol)

    if not parsed:
        return jsonify({"error": "No valid molecules parsed from input"}), 400

    collector = MatchMultiCountCollector()

    smarts.MatchCountsMulti(
        smarts_file_path=temp_smarts_file_path,
        strict_smarts=raiseError,
        usa=usa,
        molReader=parsed,
        molWriter=collector,
        exclude_mol_props=exclude_mol_props,
    )

    return jsonify(collector.results), 200


@smarts_filter.route('/get_multi_matchfilter', methods=['GET'])
def get_multi_matchfilter():
    smiles_list = process_smiles_input(request, "SMILES", 1000)
    names_list = process_smiles_input(request, "Names", 1000)

    if not isinstance(smiles_list, list) or not smiles_list:
        return jsonify({"error": "No SMILES provided"}), 400

    if names_list and len(names_list) != len(smiles_list):
        return jsonify({"error": f"'SMILES' and 'Names' must be the same length: got {len(smiles_list)} and {len(names_list)}"}), 400

    if not names_list:
        names_list = smiles_list

    try:
        smarts_list = process_multi_smarts_input(request)
    except Exception as e:
        return jsonify({"error": str(e)}), 400

    with tempfile.NamedTemporaryFile(mode="w+", delete=False) as temp_smarts_file:
        for sma in smarts_list:
            temp_smarts_file.write(sma.strip() + "\n")
        temp_smarts_file_path = temp_smarts_file.name

    exclude_mol_props = request.args.get("ExcludeMolProp", type=bool, default=False)
    raiseError = request.args.get("strict", type=bool, default=False)

    parsed = []
    for smiles, name in zip(smiles_list, names_list):
        mol = Chem.MolFromSmiles(smiles)
        if mol:
            mol.SetProp("_Name", name)
            parsed.append(mol)

    if not parsed:
        return jsonify({"error": "No valid molecules parsed from input"}), 400

    collector = MatchFilter()

    smarts.MatchFilterMulti(
        smarts_file_path=temp_smarts_file_path,
        strict_smarts=raiseError,
        molReader=parsed,
        molWriter=collector,
        exclude_mol_props=exclude_mol_props,
    )

    return jsonify(collector.accepted), 200
