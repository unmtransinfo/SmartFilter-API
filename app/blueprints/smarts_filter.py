from flask import Blueprint, request, jsonify
from rdkit import Chem
from rdktools import smarts
from rdkit.Chem import FilterCatalog
from models import PainsRun, db
from utils.request_processing import process_smiles_input, process_smarts_input, process_multi_smarts_input
import tempfile
import re

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
        self.named_smarts = named_smarts
        self.results = {}

    def SetProps(self, props):
        self.query_names = [p for p in props if p.startswith("n_matches(")]

    def write(self, mol):
        mol_name = mol.GetProp("_Name") if mol.HasProp("_Name") else ""
        if mol_name not in self.results:
            self.results[mol_name] = {}

        for prop_name in self.query_names:
            count = mol.GetIntProp(prop_name) if mol.HasProp(prop_name) else 0

            # Extract what's inside n_matches(...)
            raw_query = prop_name[len("n_matches("):-1]  # e.g. query_02 = "[OH]"

            # Use regex to extract the SMARTS part inside quotes
            match = re.search(r'=\s*"(.+)"', raw_query)
            if match:
                smarts_pattern = match.group(1)  # This gives just [OH]
            else:
                smarts_pattern = raw_query  # fallback

            # Map SMARTS pattern to smart name if it exists
            smart_name = None
            for pattern, name in self.named_smarts:
                if pattern == smarts_pattern:
                    smart_name = name
                    break
            if smart_name is None:
                smart_name = smarts_pattern

            self.results[mol_name][smart_name] = count





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
            invalid.append({"name": name, "smiles": smiles})

    if not parsed:
        return jsonify({"error": "No valid molecules found"}), 400

    # Get all PAINS filters
    all_pains_filters = []
    num_entries = catalog.GetNumEntries()
    for i in range(num_entries):
        entry = catalog.GetEntry(i)
        all_pains_filters.append(entry.GetDescription())

    results = []

    for name, smiles, mol in parsed:
        reasons = []
        highlight_atom_sets = []

        try:
            matched_entries = catalog.GetMatches(mol)
            for entry in matched_entries:
                reasons.append(entry.GetDescription())

                try:
                    matches = entry.GetFilterMatches(mol)
                    for match in matches:
                        atom_indices = [atom_idx for _, atom_idx in match.atomPairs]
                        if atom_indices:
                            highlight_atom_sets.append(atom_indices)
                except Exception as match_err:
                    print(f"[Match Error] {entry.GetDescription()} on '{name}': {match_err}")
        except Exception as e:
            print(f"[Catalog Match Error] Failed on molecule '{name}': {e}")

        results.append({
            "name": name,
            "smiles": smiles,
            "failed": bool(reasons),
            "reasons": reasons,
            "highlight_atoms": highlight_atom_sets
        })

    return jsonify({
        "results": results,
        "all_pains_filters": all_pains_filters,
        "invalid": invalid
    }), 200


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
    smart = process_smarts_input(request)[0]
    smart_name = request.args.get("Smart_Names", smart)
    exclude_mol_props = request.args.get("exclude_mol_props", type=bool, default=False)

    if len(smiles_list) != len(names_list):
        return jsonify({"error": "SMILES and Names must be same length"}), 400

    parsed = []
    invalid = []
    for smiles, name in zip(smiles_list, names_list):
        mol = Chem.MolFromSmiles(smiles)
        if mol:
            mol.SetProp("_Name", name)
            parsed.append((name, smiles, mol))
        else:
            invalid.append({"name": name, "smiles": ""})

    collector = MatchFilter()
    smarts.MatchFilter(
        smarts=smart,
        molReader=[mol for _, _, mol in parsed],
        molWriter=collector,
        exclude_mol_props = exclude_mol_props
    )

    accepted_smiles = set(item['SMILES'] for item in collector.accepted)
    failed = []
    passed = []
    qmol = Chem.MolFromSmarts(smart)

    for name, smiles, mol in parsed:
        canon = Chem.MolToSmiles(mol, canonical=True)
        matches = mol.GetSubstructMatches(qmol)
        if canon in accepted_smiles:
            failed.append({
                "name": name,
                "smiles": smiles,
                "failed": True,
                "reason": smart_name,
                "highlight_atoms": [list(m) for m in matches]
            })
        else:
            passed.append({
                "name": name,
                "smiles": smiles,
                "failed": False
            })

    failed.extend(invalid)
    return jsonify({"passed": passed, "failed": failed}), 200

@smarts_filter.route('/get_multi_matchcounts', methods=['GET'])
def get_multi_matchcounts():
    smiles_list = process_smiles_input(request, "SMILES", 1000)
    names_list = process_smiles_input(request, "Smile_Names", 1000)

    if not isinstance(smiles_list, list) or not smiles_list:
        return jsonify({"error": "No SMILES provided"}), 400

    if names_list and len(names_list) != len(smiles_list):
        return jsonify({
            "error": f"'SMILES' and 'Smile_Names' must be the same length: got {len(smiles_list)} and {len(names_list)}"
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

    if len(smarts_list) != len(smart_names):
        return jsonify({
            "error": f"'Smart_Names' and 'smarts' must be the same length: got {len(smart_names)} and {len(smarts_list)}"
        }), 400

    named_smarts = list(zip(smarts_list, smart_names))

    # Write SMARTS and names to temporary file
    with tempfile.NamedTemporaryFile(mode="w+", delete=False) as temp_smarts_file:
        for pattern, name in named_smarts:
            temp_smarts_file.write(f"{pattern} {name}\n")
        temp_smarts_file_path = temp_smarts_file.name

    exclude_mol_props = request.args.get("ExcludeMolProp", type=bool, default=False)
    usa = request.args.get("usa", type=bool, default=False)
    raiseError = request.args.get("strict", type=bool, default=False)
    nonzero_rows = request.args.get("nonzero_rows", type=bool, default=False)
    # Parse molecules
    parsed = []
    for smiles, name in zip(smiles_list, names_list):
        mol = Chem.MolFromSmiles(smiles)
        if mol:
            mol.SetProp("_Name", name)
            parsed.append(mol)

    if not parsed:
        return jsonify({"error": "No valid molecules parsed from input"}), 400

    # Run the match counts
    collector = MatchMultiCountCollector(named_smarts)
    
    smarts.MatchCountsMulti(
        smarts_file_path=temp_smarts_file_path,
        strict_smarts=raiseError,
        usa=usa,
        molReader=parsed,
        molWriter=collector,
        exclude_mol_props=exclude_mol_props,
        nonzero_rows=nonzero_rows
    )
    results = []
    for i, mol in enumerate(parsed):
        mol_name = mol.GetProp("_Name")
        result = {
            "name": mol_name,
            "smiles": smiles_list[i],
            "matches": []
        }

        for pattern, name in named_smarts:
            count = collector.results.get(mol_name, {}).get(name, 0)
            result["matches"].append({
                "smarts": pattern,
                "name": name,
                "count": count
            })

        results.append(result)


    return jsonify(results), 200


@smarts_filter.route('/get_multi_matchfilter', methods=['GET'])
def get_multi_matchfilter():
    smiles_list = process_smiles_input(request, "SMILES", 1000)
    names_list = process_smiles_input(request, "Smile_Names", 1000)

    if not isinstance(smiles_list, list) or not smiles_list:
        return jsonify({"error": "No SMILES provided"}), 400

    if names_list and len(names_list) != len(smiles_list):
        return jsonify({"error": "'SMILES' and 'Names' must be the same length"}), 400

    if not names_list:
        names_list = smiles_list

    smarts_list = request.args.getlist("smarts")
    smart_names = request.args.getlist("Smart_Names")

    if not smarts_list or not smart_names:
        return jsonify({"error": "Missing SMARTS or Smart_Names"}), 400

    if len(smarts_list) != len(smart_names):
        return jsonify({"error": "'Smart_Names' and 'smarts' must be the same length"}), 400

    exclude_mol_props = request.args.get("exclude_mol_props", type=bool, default=False)
    strict = request.args.get("strict", type=bool, default=False)

    # Write SMARTS file for MatchFilterMulti
    with tempfile.NamedTemporaryFile(mode="w+", delete=False) as temp_smarts_file:
        for smarts_pattern in smarts_list:
            temp_smarts_file.write(smarts_pattern.strip() + "\n")
        temp_smarts_file_path = temp_smarts_file.name

    parsed = []
    invalid = []
    for smiles, name in zip(smiles_list, names_list):
        mol = Chem.MolFromSmiles(smiles)
        if mol:
            mol.SetProp("_Name", name)
            parsed.append((name, smiles, mol))
        else:
            invalid.append({
                "name": name,
                "smiles": "",
                "failed": True,
                "reason": "Invalid SMILES",
                "highlight_atoms": []
            })

    if not parsed and not invalid:
        return jsonify({"error": "No valid molecules parsed"}), 400

    collector = MatchFilter()
    smarts.MatchFilterMulti(
        smarts_file_path=temp_smarts_file_path,
        strict_smarts=strict,
        molReader=[mol for _, _, mol in parsed],
        molWriter=collector,
        exclude_mol_props=exclude_mol_props,
    )

    accepted_smiles = set(item['SMILES'] for item in collector.accepted)

    failed = []
    passed = []

    # For each molecule, test which SMARTS it failed
    for name, smiles, mol in parsed:
        canon = Chem.MolToSmiles(mol, canonical=True)

        if canon in accepted_smiles:
            # Failed: find which SMARTS matched (may be more than one)
            for smarts_pattern, smart_name in zip(smarts_list, smart_names):
                qmol = Chem.MolFromSmarts(smarts_pattern)
                if not qmol:
                    continue
                matches = mol.GetSubstructMatches(qmol)
                if matches:
                    failed.append({
                        "name": name,
                        "smiles": smiles,
                        "failed": True,
                        "reason": smart_name,
                        "highlight_atoms": [list(m) for m in matches]
                    })
                    break  # stop after first match for same mol
        else:
            passed.append({
                "name": name,
                "smiles": smiles,
                "failed": False
            })

    failed.extend(invalid)

    return jsonify({"passed": passed, "failed": failed}), 200
