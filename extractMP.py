from config import mpAPI
from mp_api.client import MPRester
from pymatgen.core import Composition, Species, Structure
import json, os
from datetime import datetime

x_site_elements = ["O", "F", "Cl", "Br", "I"]
fields = [
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

class mpExtraction:
    def __init__(self, x_site_elements=x_site_elements, fields=fields, api_key=mpAPI, verbose=True):
        self.x_site_elements = x_site_elements if isinstance(x_site_elements, list) else [x_site_elements]
        self.fields = fields
        self.api_key = api_key
        self.number_inquiries = 0
        self.results = []
        self.cif_dir = "cif_files"
        self.metadata = "metadata.json"
        self.verbose = verbose
        os.makedirs(self.cif_dir, exist_ok=True)

        # Counters for debugging
        self._count_total = 0
        self._count_no_species = 0
        self._count_unbalanced = 0
        self._count_no_structure = 0
        self._count_passed = 0

    def _log(self, msg):
        #logging purpose, use verbose=False to not update the output.log
        filename = "output.log"
        if self.verbose:
            if not os.path.isfile(filename):
                with open(filename, "w", encoding="utf-8") as file:
                    now = datetime.now()
                    formatted_time = now.strftime("%d:%m:%Y %H:%M:%S")
                    file.write(f"{formatted_time}       Materials Project Data Extraction Starts")
            else:
                with open(filename, "a", encoding="utf-8") as file:
                    now = datetime.now()
                    formatted_time = now.strftime("%d:%m:%Y %H:%M:%S")
                    file.write(f"\n{formatted_time}      "+msg)

    def find_perovskite_material(self):
        seen = set()

        with MPRester(self.api_key) as mpr:
            for x_el in self.x_site_elements:
                for n_wildcards in (2, 3, 4): 
                    #X_EL-*-*-* (3,4,5 configurations)
                    pattern = f"{x_el}-" + "-".join(["*"] * n_wildcards)
                    self._log(f"[Query {self.number_inquiries + 1}] chemsys pattern: {pattern}")

                    docs = mpr.materials.summary.search( 
                        chemsys=pattern,
                        energy_above_hull=(0, 0.4),
                        fields=self.fields,
                        chunk_size=1, #set None for production, 1 for testing
                        num_chunks=1, #set None for production, 1 for testing
                    )
                    self.number_inquiries += 1
                    self._log(f"  → {len(docs)} docs returned")

                    for doc in docs:
                        self._count_total += 1
                        self._log(f"[{doc.material_id}] {doc.formula_pretty}")

                        # Check structure availability
                        if doc.structure is None:
                            self._count_no_structure += 1
                            self._log(f"    ✗ No structure available")

                        # Stage 1: oxidation balance
                        balanced, charge_sum = self._check_oxidation(doc)
                        if not balanced:
                            self._log(f"    ✗ Oxidation unbalanced (charge sum: {charge_sum})")
                            self._count_unbalanced += 1
                            continue

                        self._log(f"    ✓ Oxidation balanced (charge sum: {charge_sum})")

                        # Stage 2: 
                        
                        if doc.material_id not in seen:
                            seen.add(doc.material_id)
                            self._count_passed += 1
                            self.results.append(
                                {k: v for k, v in doc.model_dump().items()
                                 if k != "fields_not_requested"}
                            )

        self._print_summary()
        return self.results

    def _print_summary(self):
        print(f"\n{'='*40}")
        print(f"Total API calls       : {self.number_inquiries}")
        print(f"Total docs fetched    : {self._count_total}")
        print(f"No possible_species   : {self._count_no_species}")
        print(f"Unbalanced charge     : {self._count_unbalanced}")
        print(f"No structure          : {self._count_no_structure}")
        print(f"Passed all filters    : {self._count_passed}")
        print(f"{'='*40}\n")

    def _check_oxidation(self, doc):
        """Returns (bool, charge_sum)"""
        if not doc.possible_species:
            self._count_no_species += 1
            return False, None

        charge_sum = None
        try:
            elem_amts = {str(el): float(amt) for el, amt in doc.composition_reduced.items()}
            ox_states = {}
            for species in doc.possible_species:
                s = Species.from_str(species)
                ox_states[str(s.element)] = s.oxi_state

            if not all(el in ox_states for el in elem_amts):
                return False, None

            charge_sum = round(sum(ox_states[el] * amt for el, amt in elem_amts.items()), 4)
            return abs(charge_sum) < 0.01, charge_sum

        except Exception as e:
            self._log(f"    ! Exception in oxidation check: {e}")
            return False, charge_sum

    def to_json(self):
        results = []
        for entry in self.results:
            structure = entry.pop("structure", None)
            if structure is not None:
                try:
                    from pymatgen.core import Structure
                    if isinstance(structure, dict):
                        structure = Structure.from_dict(structure)
                    cif_str = structure.to(fmt="cif")
                    cif_path = os.path.join(self.cif_dir, f"{entry['material_id']}.cif")
                    with open(cif_path, "w") as f:
                        f.write(cif_str)
                    entry["cif_file"] = cif_path
                    self._log(f"  CIF saved: {cif_path}")
                except Exception as e:
                    entry["cif_file"] = None
                    entry["cif_error"] = str(e)
                    self._log(f"  ! CIF error for {entry['material_id']}: {e}")
            else:
                entry["cif_file"] = None
                self._count_no_structure += 1

            entry.pop("fields_not_requested", None)
            results.append(entry)

        with open(self.metadata, "w") as f:
            json.dump(results, f, indent=2, default=str)

        print(f"Saved {len(results)} entries to {self.metadata}")
        print(f"CIF files saved to ./{self.cif_dir}/")

mpEx = mpExtraction(verbose=True)
mpEx.find_perovskite_material()
mpEx.to_json()