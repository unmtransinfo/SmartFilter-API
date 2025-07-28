from flask import Blueprint, request, jsonify
from rdkit import Chem
from rdktools import smarts
from rdkit.Chem import FilterCatalog
from models import PainsRun, db
from utils.request_processing import process_smiles_input, process_smarts_input, process_multi_smarts_input
import tempfile

smarts_filter = Blueprint("smarts_filter", __name__, url_prefix="/smarts_filter")


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
        self.accepted_names.append([name, canon])  # Use list, not tuple


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
    def __init__(self, named_smarts):
        self.named_smarts = named_smarts  # list of (name, smarts)
        self.results = []

    def SetProps(self, props):
        self.query_names = []
        for prop in props:
            if prop.startswith("n_matches("):
                self.query_names.append(prop)

    def write(self, mol):
        record = {
            "smiles": Chem.MolToSmiles(mol, canonical=True),
            "name": mol.GetProp("_Name") if mol.HasProp("_Name") else "",
        }

        for key in self.query_names:
            print(key)
            if mol.HasProp(key):
                record[key] = mol.GetIntProp(key)
            else:
                record[key] = 0

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

# Initialize the PAINS filter catalog once (reuse on every request)
params = FilterCatalog.FilterCatalogParams()
params.AddCatalog(FilterCatalog.FilterCatalogParams.FilterCatalogs.PAINS_A)
params.AddCatalog(FilterCatalog.FilterCatalogParams.FilterCatalogs.PAINS_B)
params.AddCatalog(FilterCatalog.FilterCatalogParams.FilterCatalogs.PAINS_C)
catalog = FilterCatalog.FilterCatalog(params)

@smarts_filter.route('/get_filterpains', methods=['GET'])
def get_filterpains():
    smiles_list = process_smiles_input(request, "SMILES", 1000)
    names_list = process_smiles_input(request, "Smile_Names", 1000)

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
            invalid.append({"name": name, "smiles": ""})

    if not parsed:
        return jsonify({"error": "No valid molecules found"}), 400

    results = []
    for name, smiles, mol in parsed:
        entry = catalog.GetFirstMatch(mol)
        if entry:
            description = entry.GetDescription()
            highlight_atom_sets = []
            for filter_match in entry.GetFilterMatches(mol):
                atom_indices = [atom_idx for _, atom_idx in filter_match.atomPairs]
                highlight_atom_sets.append(atom_indices)

            results.append({
                "name": name,
                "smiles": smiles,
                "failed": True,
                "reason": description,
                "highlight_atoms": highlight_atom_sets
            })
        else:
            results.append({
                "name": name,
                "smiles": smiles,
                "failed": False
            })

    # Add invalid molecules as failed with empty smiles
    results.extend([{"name": inv["name"], "smiles": "", "failed": True, "reason": "Invalid SMILES", "highlight_atoms": []} for inv in invalid])

    return jsonify(results), 200

@smarts_filter.route('/get_matchcounts', methods=['GET'])
def get_matchcounts():
    smiles_list = process_smiles_input(request, "SMILES", 1000)
    names_list = process_smiles_input(request, "Smile_Names", 1000)

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
    non_zero_rows = request.args.get("nonzero_rows", type=bool, default=False)
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
        nonzero_rows=non_zero_rows,
    )

    return jsonify(collector.results), 200


@smarts_filter.route('/get_matchfilter', methods=['GET'])
def get_matchfilter():
    smiles_list = process_smiles_input(request, "SMILES", 1000)
    names_list = process_smiles_input(request, "Smile_Names", 1000)

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
    invalid = []
    for smiles, name in zip(smiles_list, names_list):
        mol = Chem.MolFromSmiles(smiles)
        if mol:
            mol.SetProp("_Name", name)
            parsed.append((name, smiles, mol))
        else:
            invalid.append({"name": name, "smiles": ""})

    if not parsed and not invalid:
        return jsonify({"error": "No valid molecules parsed from input"}), 400

    collector = MatchFilter()
    smarts.MatchFilter(
        smarts=smart,
        molReader=[mol for _, _, mol in parsed],
        molWriter=collector,
        exclude_mol_props=exclude_mol_props,
    )

    accepted_smiles = set(item['SMILES'] for item in collector.accepted)

    passed = []
    failed = []

    for name, smiles, mol in parsed:
        if Chem.MolToSmiles(mol, canonical=True) in accepted_smiles:
            failed.append({"name": name, "smiles": smiles})
        else:
            passed.append({"name": name, "smiles": smiles})

    failed.extend(invalid)

    return jsonify({
        "passed": passed,
        "failed": failed
    }), 200


@smarts_filter.route('/get_multi_matchcounts', methods=['GET'])
def get_multi_matchcounts():
    smiles_list = process_smiles_input(request, "SMILES", 1000)
    names_list = process_smiles_input(request, "Smile_Names", 1000)

    if not isinstance(smiles_list, list) or not smiles_list:
        return jsonify({"error": "No SMILES provided"}), 400

    if names_list and len(names_list) != len(smiles_list):
        return jsonify({
            "error": f"'SMILES' and 'Names' must be the same length: got {len(smiles_list)} and {len(names_list)}"
        }), 400

    if not names_list:
        names_list = smiles_list

    try:
        smarts_list = request.args.getlist("smarts")
        smart_names = request.args.getlist("Smart_Names")
    except Exception as e:
        return jsonify({"error": str(e)}), 400

    if not smarts_list or not smart_names:
        return jsonify({"error": "Missing SMARTS patterns or Smart_Names"}), 400

    print(smarts_list, smart_names)
    if len(smarts_list) != len(smart_names):
        return jsonify({
            "error": f"'Smart_Names' and 'smarts' must be the same length: got {len(smart_names)} and {len(smarts_list)}"
        }), 400

    # Zip name and pattern
    named_smarts = list(zip(smarts_list, smart_names))
    # Write SMARTS only to file (not names)
    with tempfile.NamedTemporaryFile(mode="w+", delete=False) as temp_smarts_file:
        for pattern, name in named_smarts:
            temp_smarts_file.write(f"{pattern} {name}\n")
        temp_smarts_file_path = temp_smarts_file.name

    # Now read and print the file contents
    with open(temp_smarts_file_path, 'r') as f:
        print("Temp SMARTS file content:")
        print(f.read())


    exclude_mol_props = request.args.get("ExcludeMolProp", type=bool, default=False)
    usa = request.args.get("usa", type=bool, default=False)
    raiseError = request.args.get("strict", type=bool, default=False)
    non_zero_rows = request.args.get("nonzero_rows", type=bool, default=False)
    parsed = []
    for smiles, name in zip(smiles_list, names_list):
        mol = Chem.MolFromSmiles(smiles)
        if mol:
            mol.SetProp("_Name", name)
            parsed.append(mol)

    if not parsed:
        return jsonify({"error": "No valid molecules parsed from input"}), 400

    # Custom collector to tag match counts with SMARTS names
    collector = MatchMultiCountCollector(named_smarts)

    smarts.MatchCountsMulti(
        smarts_file_path=temp_smarts_file_path,
        strict_smarts=raiseError,
        usa=usa,
        molReader=parsed,
        molWriter=collector,
        exclude_mol_props=exclude_mol_props,
        nonzero_rows = non_zero_rows
    )
    return jsonify(collector.results), 200
@smarts_filter.route('/get_multi_matchfilter', methods=['GET'])
def get_multi_matchfilter():
    smiles_list = process_smiles_input(request, "SMILES", 1000)
    names_list = process_smiles_input(request, "Smile_Names", 1000)

    if not isinstance(smiles_list, list) or not smiles_list:
        return jsonify({"error": "No SMILES provided"}), 400

    if names_list and len(names_list) != len(smiles_list):
        return jsonify({"error": f"'SMILES' and 'Names' must be the same length: got {len(smiles_list)} and {len(names_list)}"}), 400

    if not names_list:
        names_list = smiles_list

    try:
        smarts_list = request.args.getlist("smarts")
    except Exception as e:
        return jsonify({"error": str(e)}), 400

    if not smarts_list:
        return jsonify({"error": "No valid SMARTS patterns provided"}), 400

    # Write SMARTS patterns to temp file
    with tempfile.NamedTemporaryFile(mode="w+", delete=False) as temp_smarts_file:
        for smarts_pattern in smarts_list:
            temp_smarts_file.write(smarts_pattern.strip() + "\n")
        temp_smarts_file_path = temp_smarts_file.name

    exclude_mol_props = request.args.get("ExcludeMolProp", type=bool, default=False)
    raiseError = request.args.get("strict", type=bool, default=False)

    parsed = []
    invalid = []
    for smiles, name in zip(smiles_list, names_list):
        mol = Chem.MolFromSmiles(smiles)
        if mol:
            mol.SetProp("_Name", name)
            parsed.append((name, smiles, mol))  # keep name and smiles for output
        else:
            invalid.append({"name": name, "smiles": ""})

    if not parsed and not invalid:
        return jsonify({"error": "No valid molecules parsed from input"}), 400

    collector = MatchFilter()
    smarts.MatchFilterMulti(
        smarts_file_path=temp_smarts_file_path,
        strict_smarts=raiseError,
        molReader=[mol for _, _, mol in parsed],
        molWriter=collector,
        exclude_mol_props=exclude_mol_props,
    )

    accepted_smiles = set(item['SMILES'] for item in collector.accepted)

    passed = []
    failed = []

    for name, smiles, mol in parsed:
        if Chem.MolToSmiles(mol, canonical=True) in accepted_smiles:
            failed.append({"name": name, "smiles": smiles})
        else:
            passed.append({"name": name, "smiles": smiles})

    failed.extend(invalid)

    return jsonify({
        "passed": passed,
        "failed": failed
    }), 200
