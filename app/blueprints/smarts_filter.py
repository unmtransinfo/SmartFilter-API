from flask import Blueprint, request, jsonify
from rdkit import Chem
from rdktools import smarts
from rdkit import Chem
from models import PainsRun, db
from utils.request_processing import process_smiles_input, process_smarts_input, process_multi_smarts_input
import tempfile

smarts_filter = Blueprint("smarts_filter", __name__, url_prefix="/smarts_filter")


def parse_smiles_input(inputs, smiles_col=0, name_col=1, delim=" "):
    """Parses list of SMILES strings into RDKit Mol objects with optional names."""
    parsed = []
    invalid = []

    for line in inputs:
        parts = line.strip().split(delim)
        if len(parts) <= smiles_col:
            invalid.append(line)
            continue

        smiles = parts[smiles_col]
        name = parts[name_col] if len(parts) > name_col else line.strip()
        mol = Chem.MolFromSmiles(smiles)
        if mol:
            mol.SetProp("_Name", name)
            parsed.append((name, smiles, mol))
        else:
            invalid.append(name)

    return parsed, invalid


class MoleculeCollector:
    """Collector for canonical SMILES that pass filters."""
    def __init__(self):
        self.accepted = []

    def SetProps(self, props):
        pass

    def write(self, mol):
        canon = Chem.MolToSmiles(mol, canonical=True)
        self.accepted.append(canon)


class MatchCountCollector:
    """Collector for molecules and their SMARTS match counts."""
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
    """Collector for molecules and their SMARTS match counts."""
    def __init__(self):
        self.results = []

    def SetProps(self, props):
        self.props = props  # this will contain the SMARTS labels like "n_matches_01", etc.

    def write(self, mol):
        record = {
            "smiles": Chem.MolToSmiles(mol, canonical=True),
            "name": mol.GetProp("_Name") if mol.HasProp("_Name") else "",
        }

        # Collect all SMARTS match counts
        for prop in self.props:
            if mol.HasProp(prop):
                record[prop] = mol.GetIntProp(prop)
            else:
                record[prop] = 0  # fallback if missing

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
    data = process_smiles_input(request, "SMILES", 1000)
    if not isinstance(data, list) or not data:
        return jsonify({"error": "No molecules provided"}), 400

    parsed, invalid = parse_smiles_input(data, smiles_col=0, name_col=0, delim=" ")

    if not parsed:
        return jsonify({"error": "No valid molecules found"}), 400

    collector = MoleculeCollector()
    mols = [mol for _, _, mol in parsed]
    exclude_mol_props = request.args.get("ExcludeMolProp", type=bool)

    smarts.FilterPAINS(mols, collector, exclude_mol_props)
    accepted_set = set(collector.accepted)

    passed = []
    failed = []

    for name, canon, _ in parsed:
        if canon in accepted_set:
            passed.append(name)
        else:
            failed.append(name)

    failed.extend(invalid)

    record = PainsRun(inputs=[{"smiles": e} for e in data], passed=passed, failed=failed)
    db.session.add(record)
    db.session.commit()

    return jsonify({"passed": passed, "failed": failed}), 200


@smarts_filter.route('/get_matchcounts', methods=['GET'])
def get_matchcounts():
    data = process_smiles_input(request, "SMILES", 1000)
    if not isinstance(data, list) or not data:
        return jsonify({"error": "No molecules provided"}), 400

    smart = process_smarts_input(request)[0]
    print(smart)
    if not smart:
        return jsonify({"error": "Smarts pattern is required"}), 400

    # Read optional args
    smiles_col = request.args.get("smiles_column", default=0, type=int)
    name_col = request.args.get("name_column", default=1, type=int)
    delim = request.args.get("delim", default=" ")
    exclude_mol_props = request.args.get("ExcludeMolProp", type=bool, default=False)
    usa = request.args.get("usa", type=bool, default=False)

    parsed, _ = parse_smiles_input(data, smiles_col, name_col, delim)

    if not parsed:
        return jsonify({"error": "No valid molecules found"}), 400

    mols = [mol for _, _, mol in parsed]
    collector = MatchCountCollector()

    smarts.MatchCounts(
        smarts=smart,
        usa=usa,
        molReader=mols,
        molWriter=collector,
        exclude_mol_props=exclude_mol_props,
    )

    return jsonify(collector.results), 200

@smarts_filter.route('/get_matchFilter', methods=['GET'])
def get_matchfilter():
    data = process_smiles_input(request, "SMILES", 1000)
    if not isinstance(data, list) or not data:
        return jsonify({"error": "No molecules provided"}), 400

    smart = process_smarts_input(request)[0]
    if not smart:
        return jsonify({"error": "Smarts pattern is required"}), 400

    smiles_col = request.args.get("smiles_column", default=0, type=int)
    name_col = request.args.get("name_column", default=1, type=int)
    delim = request.args.get("delim", default=" ")
    exclude_mol_props = request.args.get("ExcludeMolProp", type=bool, default=False)

    parsed, _ = parse_smiles_input(data, smiles_col, name_col, delim)

    if not parsed:
        return jsonify({"error": "No valid molecules found"}), 400

    mols = [mol for _, _, mol in parsed]
    collector = MatchFilter()

    smarts.MatchFilter(
        smarts=smart,
        molReader=mols,
        molWriter=collector,
        exclude_mol_props=exclude_mol_props,
    )

    return jsonify(collector.accepted), 200


@smarts_filter.route('/get_multi_matchcounts', methods=['GET'])
def get_multi_matchcounts():
    data = process_smiles_input(request, "SMILES", 1000)
    if not isinstance(data, list) or not data:
        return jsonify({"error": "No molecules provided"}), 400

    # Get list of validated SMARTS (using helper that supports comma/space separation)
    try:
        smarts_list = process_multi_smarts_input(request, "smarts")
    except Exception as e:
        return jsonify({"error": str(e)}), 400

    if not smarts_list:
        return jsonify({"error": "No valid SMARTS patterns provided"}), 400

    # Write SMARTS to a temporary file
    with tempfile.NamedTemporaryFile(mode="w+", delete=False) as temp_smarts_file:
        for sma in smarts_list:
            temp_smarts_file.write(sma.strip() + "\n")
        temp_smarts_file_path = temp_smarts_file.name

    # Optional params
    smiles_col = request.args.get("smiles_column", default=0, type=int)
    name_col = request.args.get("name_column", default=1, type=int)
    delim = request.args.get("delim", default=" ")
    exclude_mol_props = request.args.get("ExcludeMolProp", type=bool, default=False)
    usa = request.args.get("usa", type=bool, default=False)
    raiseError = request.args.get("strict", type=bool, default=False)

    # Parse SMILES input
    parsed, _ = parse_smiles_input(data, smiles_col, name_col, delim)
    if not parsed:
        return jsonify({"error": "No valid molecules found"}), 400

    mols = [mol for _, _, mol in parsed]
    collector = MatchMultiCountCollector()

    # Call SMARTS matcher
    smarts.MatchCountsMulti(
        smarts_file_path=temp_smarts_file_path,
        strict_smarts=raiseError,
        usa=usa,
        molReader=mols,
        molWriter=collector,
        exclude_mol_props=exclude_mol_props,
    )

    return jsonify(collector.results), 200


@smarts_filter.route('/get_multi_matchfilter', methods=['GET'])
def get_multi_matchfilter():
    data = process_smiles_input(request, "SMILES", 1000)
    if not isinstance(data, list) or not data:
        return jsonify({"error": "No molecules provided"}), 400

    # Use the reusable SMARTS parsing function
    try:
        smarts_list = process_multi_smarts_input(request)
    except Exception as e:
        return jsonify({"error": str(e)}), 400

    # Write SMARTS list to a temporary file
    with tempfile.NamedTemporaryFile(mode="w+", delete=False) as temp_smarts_file:
        for sma in smarts_list:
            temp_smarts_file.write(sma.strip() + "\n")
        temp_smarts_file_path = temp_smarts_file.name

    smiles_col = request.args.get("smiles_column", default=0, type=int)
    name_col = request.args.get("name_column", default=1, type=int)
    delim = request.args.get("delim", default=" ")
    raiseError = request.args.get("strict", type=bool, default=False)
    exclude_mol_props = request.args.get("ExcludeMolProp", type=bool, default=False)

    parsed, _ = parse_smiles_input(data, smiles_col, name_col, delim)
    if not parsed:
        return jsonify({"error": "No valid molecules found"}), 400

    mols = [mol for _, _, mol in parsed]
    collector = MatchFilter()

    smarts.MatchFilterMulti(
        smarts_file_path=temp_smarts_file_path,
        strict_smarts=raiseError,
        molReader=mols,
        molWriter=collector,
        exclude_mol_props=exclude_mol_props,
    )

    return jsonify(collector.accepted), 200
