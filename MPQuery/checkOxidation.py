from itertools import product
import json, os, logging
from pymatgen.core import Species, Element

class checkOxidation():
    def __init__(self, tolerance = 0.01, TEST_JSON_PATH = "TEST_QUERY.json", OUTPUT_JSON_PATH = "OXIDATION_QUERY.json",
                 logger = None):
        """
        Checks whether a material's composition can satisfy charge
        neutrality, by trying combinations of oxidation states pulled
        from Materials Project's `possible_species` field.

        tolerance : float, optional
            How far off zero total charge is still counted as "neutral".
            Default 0.01.

        TEST_JSON_PATH : str, optional
            Filename of the input query file to read from.

        OUTPUT_JSON_PATH : str, optional
            Filename to write oxidation-state results to.
        """

        self.tolerance = tolerance 
        self.json_path = TEST_JSON_PATH
        self.output_json = OUTPUT_JSON_PATH
        self.logger  = logger if logger is not None else logging.getLogger(__name__)
        self._TIER_SEVERITY = {"common": 0,
                          "icsd": 1,
                          "exotic": 2}
    def load_json(self):
        with open(self.json_path, "r", encoding="utf-8") as f:
            return json.load(f)
    
    def save_json(self, metadata):
        """
        Save metadata to JSON
        """

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
        Reads doc["possible_species"] 
        groups the oxidation states by element. 
        Non-integer oxidation states are dropped
        """

        candidates = {}
        dropped = []
        for species in doc.get("possible_species", []):
            sp = Species.from_str(species)
            if not float(sp.oxi_state).is_integer():
                dropped.append(species)
                continue
            element = str(sp.element)
            candidates.setdefault(element, set()).add(sp.oxi_state)
        
        if dropped:
            self.logger.warning(
                f"{doc.get('material_id')}: dropped fractional oxidation "
                f"state candidate(s) {dropped}"
            )
        return candidates

    def _rank_states(self, element_symbol, states):
        """
        Order one element's oxidation states based on its chemical plausibility tier
        and tags with their label.
        common (1): own preferred order 
        icsd (2): covered in icsd
        exotic (3): else
        """
        try:
            el = Element(element_symbol)
            common = [s for s in el.common_oxidation_states if s in states]
            icsd = [s for s in getattr(el, "icsd_oxidation_states", ())]
        except ValueError:
            common, icsd = [], []
        
        rest = sorted(s for s in states if s not in common and s not in icsd)
        return([(s, "common") for s in common] +
               [(s, "icsd") for s in icsd] +
               [(s, "exotic") for s in rest])
    
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
        combos from `candidates`, lookup from rank states) and looks for one that adds up to ~0
        total charge. 
        """

        elements = list(composition.keys())
        option_lists = [self._rank_states(e, candidates[e]) for e in elements]

        best_charge = None 
        best_assignment = None 
        best_tier = None

        for combo in product(*option_lists):
            states = tuple(state for state, _tier in combo)
            tiers = [tier for _state, tier in combo]
            overall_tier = max(tiers, key=lambda t: self._TIER_SEVERITY[t])

            charge = self._evaluate_combination(
                composition, elements, states)
            assignment = dict(zip(elements, states))

            if best_charge is None or abs(charge)<abs(best_charge):
                best_charge = charge
                best_assignment = assignment
                best_tier = overall_tier
            
            if abs(charge) <= self.tolerance:
                return True, charge, assignment, overall_tier

        return False, best_charge, best_assignment, best_tier

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

        material_id = doc.get("material_id")
        formula_pretty = doc.get("formula_pretty")
        state_tier = None

        try:

            if not doc.get("possible_species"):
                self.logger.error(f"Possible species of {doc['material_id']} is not found")
                isNeutral, charge, assignment = False, None, None
            else:
                composition = self._get_composition(doc)
                candidates = self._get_oxidation_states(doc)
                if not all(element in candidates for element in composition):
                    self.logger.error(f"Inconsistency in {material_id} oxidation states and composition")
                    isNeutral, charge, assignment = False, None, None
                else:        
                    isNeutral, charge, assignment, state_tier = self._find_best_assigment(composition, candidates)
        except (KeyError, TypeError, ValueError) as e:
            self.logger.error(f"Oxidation-state check failed ({material_id}): {e!r}")
            isNeutral, charge, assignment = False, None, None
        
        return  {"formula_pretty": formula_pretty,
                "material_id": material_id,
                "is_neutral": isNeutral,
                "total_charge": charge,
                "ion_assignment": assignment,
                "state_tier": state_tier
                }
    
