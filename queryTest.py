import os, json, argparse, ast
from pymatgen.core import Structure
from config import mpAPI
from mp_api.client import MPRester

FIELDS  = [
    # Identifiers
    "material_id", "formula_pretty", "formula_anonymous",
    "chemsys", "nelements", "elements",
    "deprecated", "theoretical",
    # Structural
    "symmetry",
    "composition",
    "nsites",
    "structure",
    # Stability
    "formation_energy_per_atom",
    "energy_above_hull",
    "is_stable",
    # Electronic
    "band_gap", "is_gap_direct", "is_metal",
    "cbm", "vbm", "efermi",
    # Oxidation states
    "possible_species",
    "composition_reduced"
]

class queryStructure():
    def __init__(self, structural_ID, fields=FIELDS, CIF_DIRS = "TEST_CIF", JSON_PATH = "TEST_QUERY"):
        """
        Query exactly one structure from Materials Project
        structural_ID: str
        fields: list
        """
        self.IDs = structural_ID if isinstance(structural_ID, list) else [structural_ID]
        self.fields = fields
        self.cif_dirs = CIF_DIRS
        self.json_path = os.path.join(JSON_PATH+".json")

    def query(self):
        with MPRester(mpAPI) as mpr:
            docs = mpr.materials.summary.search(
                material_ids=self.IDs,
                fields = self.fields,
            )
        if len(docs) == 0:
            raise ValueError(f"No material found for {self.IDs}")
        
        for doc in docs:
        
            data = {
                k: v for k, v in doc.model_dump().items()
                if k!= "fields_not_requested"
            }
            self._save_cif(data)
            self._save_json(data)
            print(f"Structure {doc['material_id']} succesfully queried")

    def _save_cif(self, data):
        """
        Save structure to CIF files
        """
        os.makedirs(self.cif_dirs, exist_ok=True)
        structure = data.get("structure")
        if structure is None:
            raise ValueError("No structure available")
        if isinstance(structure, dict):
            structure = Structure.from_dict(structure)
        
        cif_path = os.path.join(self.cif_dirs,f"{data['material_id']}.cif")
        with open(cif_path, "w") as f:
            f.write(structure.to(fmt="cif"))
    
    def _save_json(self, data):
        """
        Save metadata to JSON
        """
        metadata = data.copy()
        metadata.pop("structure", None)

        if os.path.exists(self.json_path):
            with open(self.json_path, "r") as f:
                try:
                    results = json.load(f)
                    if not isinstance(results, list):
                        results = [results]
                except json.JSONDecodeError:
                    results = []
        else:
            results = []
        
        results.append(metadata)
        with open(self.json_path, "w") as f:
            json.dump(results, f, indent=2)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Query one structure from Materials Project"
    )

    parser.add_argument(
        "-q", "--query",
        nargs="+",
        required=True,
    )
    args = parser.parse_args()
    structural_ID = args.query
    query = queryStructure(structural_ID)
    query.query()

#["mp-4651", "mp-2998", "mp-5827", "mp-4019", "mp-1079517", "mp-19127", "mp-18857"]