import os, json, argparse, ast, logging
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
    def __init__(self, fields=FIELDS, CIF_DIRS = "TEST_CIF", JSON_PATH = "TEST_QUERY.json",
                 logger = None):
        """
        Query structures from Materials Project based on its MPIDs
        structural_ID   :   str or list(str)
            MPIDs for a given structure, enter using argument parser
        fields          :   list(str), optional
            fields of the projected output for MP queries on each structure,
            default is listed
        CIF_DIRS        :   str, optional
            directory path for cif files of the queried structure
        JSON_PATH       :   str, optional
            file path for the json files. 
        """

        self.fields = fields
        self.cif_dirs = CIF_DIRS
        self.json_path = JSON_PATH
        self.logger  = logger if logger is not None else logging.getLogger(__name__)

    def query(self, ID):
        """
        Main function for querying MPIDs
        """
        IDs = ID if isinstance(ID, list) else [ID]
        with MPRester(mpAPI) as mpr:
            docs = mpr.materials.summary.search(
                material_ids=IDs,
                fields = self.fields,
            )
        if len(docs) == 0:
            self.logger.error(f"No material found for {IDs}")

        return docs
    
    def processing_query(self, docs):
        for doc in docs:
            data = {
                k: v for k, v in doc.model_dump().items()
                if k!= "fields_not_requested"
            }
            data['material_id'] = doc['material_id']
            self.save_cif(data)
            self.save_json(data)

    def save_cif(self, data):
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
    
    def save_json(self, data):
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
