from itertools import product
import json, os
from pymatgen.core import Species 

class checkOxidation():
    def __init__(self, tolerance = 0.01, TEST_JSON_PATH = "TEST_QUERY", OUTPUT_JSON_PATH = "OXIDATION_QUERY"):
        """
        Checks whether a material's composition can satisfy charge
        neutrality, by trying combinations of oxidation states pulled
        from Materials Project's `possible_species` field.

        tolerance : float, optional
            How far off zero total charge is still counted as "neutral".
            Default 0.01.

        TEST_JSON_PATH : str, optional
            Filename (without ".json") of the input query file to read from.

        OUTPUT_JSON_PATH : str, optional
            Filename (without ".json") to write oxidation-state results to.
        """

        self.tolerance = tolerance 
        self.json_path = TEST_JSON_PATH + ".json"
        self.output_json = OUTPUT_JSON_PATH + ".json"
    
    def load_json(self):
        with open(self.json_path, "r", encoding="utf-8") as f:
            return json.load(f)
    
    def save_json(self, query, charge, assignment):
        """
        Save metadata to JSON
        """
        metadata = {"formula_pretty": query["formula_pretty"],
                "material_id": query["material_id"],
                "total_charge": charge,
                "ion_assignment": assignment
                }

        if os.path.exists(self.output_json):
            with open(self.output_json, "r") as f:
                try:
                    results = json.load(f)
                    if not isinstance(results, list):
                        results = [results]
                except json.JSONDecodeError:
                    results = []
        else:
            results = []
        
        results.append(metadata)
        with open(self.output_json, "w") as f:
            json.dump(results, f, indent=2)

    def _get_composition(self, doc):
        """
        Extract reduced composition
        """
        return {
            str(el): float(amount) 
            for el, amount in doc["composition_reduced"].items()
        }
    
    def _get_oxidation_states(self, doc):
        """
        Reads doc["possible_species"] (e.g. ["Ca2+", "Ti4+", "O2-"]) and
        groups the oxidation states by element.
        """

        candidates = {}
        for species in doc.get("possible_species", []):
            sp = Species.from_str(species)
            element = str(sp.element)
            candidates.setdefault(element, set()).add(sp.oxi_state)
        
        return candidates

    def _evaluate_combination(self, composition, elements, oxidation_states):
        """
        Adds up total charge for one specific choice of oxidation state
        per element: sum(oxidation_state * how_many_atoms) across all
        elements.
        """
        return round(sum(ox*composition[element]
                for element, ox in zip(elements, oxidation_states))
                ,4)
    
    def _find_best_assigment(self, composition, candidates):
        """
        Tries every combination of oxidation states (one per element, all
        combos from `candidates`) and looks for one that adds up to ~0
        total charge. 
        """




        elements = list(composition.keys())
        option_lists = [sorted(candidates[e]) for e in elements]

        best_charge = None 
        best_assignment = None 

        for combo in product(*option_lists):
            charge = self._evaluate_combination(
                composition, elements, combo)
            assignment = dict(zip(elements, combo))

            if best_charge is None or abs(charge)<abs(best_charge):
                best_charge = charge
                best_assignment = assignment
            
            if abs(charge) <= self.tolerance:
                return True, charge, assignment

        return False, best_charge, best_assignment

    def check_charge(self, doc):
        """
        Main entry point: given one Materials Project doc, checks if it
        can be assigned oxidation states that balance to zero charge.

        Returns (is_neutral, charge, assignment):
          - is_neutral: True if a neutral (or near-neutral, within
            tolerance) assignment was found.

          - charge: the total charge for the best assignment found.

          - assignment: dict of {element: oxidation_state} for that
            assignment.

        Returns (False, None, None) if there's no possible_species data,
        or if some element in the composition has no candidate oxidation
        state to try at all.
        """

        if not doc.get("possible_species"):
            return False, None, None
        composition = self._get_composition(doc)
        candidates = self._get_oxidation_states(doc)

        if not all(element in candidates for element in composition):
            return False, None, None
        return self._find_best_assigment(
            composition, candidates
        )
    
def write_log(log, filename = "oxidation.log"):
    if not os.path.isfile(filename):
        with open(filename, "w", encoding = "utf-8") as file:
            file.write(log)
    else:
        with open(filename, "a", encoding = "utf-8") as file:
            file.write("\n"+log)

if __name__ == "__main__":
    cO = checkOxidation()
    if os.path.exists("oxidation.log"):
        os.remove("oxidation.log")
        os.remove("OXIDATION_QUERY.json")
    queries = cO.load_json()
    for query in queries:
        isNeutral, charge, assignment = cO.check_charge(query)
        log = f"Structure {query['formula_pretty']} ({query['material_id']}) -> Total Charge: {charge} with ion assignment {assignment}"
        write_log(log)
        cO.save_json(query, charge, assignment)
        print(f"{query['material_id']} is logged")