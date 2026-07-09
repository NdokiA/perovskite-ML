import json, math, os
from pymatgen.core import Species
import traceback
import logging
_CN_ROMAN = {1: "I", 2: "II", 3: "III", 4: "IV", 5: "V", 6: "VI",
             7: "VII", 8: "VIII", 9: "IX", 10: "X", 11: "XI", 12: "XII"}
 
class checkTolerance():
    def __init__(self, JSON_PATH = "TEST_QUERY.json", OXI_PATH="OXIDATION_QUERY.json", NEIGH_PATH="NEIGHBOR_QUERY.json",
                OUTPUT_JSON = "TOLERANCE_QUERY.json", logger = None):
        """
        Computes perovskite tolerance factors (octahedral factor mu,
        Goldschmidt t, Bartel tau) from checkOxidation.py and
        checkNeighbor.py results. 
        Every radius/oxidation-state value is a lookup by element, not a
        measurement off the actual structure.
 
        OXI_PATH   : str, filename of checkOxidation.py's output
        NEIGH_PATH : str, filename of checkNeighbor.py's output
        """
        self.OXI_PATH = OXI_PATH
        self.NEIGH_PATH = NEIGH_PATH
        self.JSON_PATH = JSON_PATH
        self.OUTPUT_PATH = OUTPUT_JSON
        self.logger  = logger if logger is not None else logging.getLogger(__name__)
    
    def load_json(self):
        with open(self.JSON_PATH, "r", encoding="utf-8") as f:
            return json.load(f)

    def save_json(self, metadata):
        """
        Save metadata to JSON
        """

        if os.path.exists(self.OUTPUT_PATH):
            with open(self.OUTPUT_PATH, "r") as f:
                try:
                    results = json.load(f)
                    if not isinstance(results, list):
                        results = [results]
                except json.JSONDecodeError:
                    results = []
        else:
            results = []
        
        results.append(metadata)
        with open(self.OUTPUT_PATH, "w") as f:
            json.dump(results, f, indent=2)

    def _load_oxi_index(self):
        """index checkOxidation.py results by material_id"""
        with open(self.OXI_PATH, "r", encoding="utf-8") as f:
            data = json.load(f)
        return {d["material_id"]: d.get("ion_assignment", {}) for d in data}
 
    def _get_oxi_map(self, doc):
        """Oxidation-state dict for this material (element -> charge)."""
        return self.oxi_index.get(doc.get("material_id"))
 
    def _load_ion_assigments(self):
        """
        Loads checkNeighbor.py's verify_bx6 results, indexed by material_id.
        Keeps the whole dict per material -- need is_perovskite, x_elements,
        b_candidates, AND b_true all together, not just a few fields.
        """
        with open(self.NEIGH_PATH, "r", encoding="utf-8") as f:
            data = json.load(f)
        return {d["material_id"]: d for d in data}
 
    def _get_ion_assignments(self, doc):
        """
        A/B/X element roles for this material, gated on is_perovskite.
 
        Returns (a_elements, b_elements, x_elements) as lists -- lists,
        not single strings, because double perovskites can have >1
        element per role. Returns (None, None, None) if this material
        wasn't confirmed as perovskite by checkNeighbor.
        """
        entry = self.ion_index.get(doc.get("material_id"))
        if not entry or not entry.get("is_perovskite"):
            return None, None, None
        x_elements = entry["x_elements"]
        b_elements = entry["b_true"]
        a_elements = sorted(set(entry["b_candidates"]) - set(b_elements))
        return a_elements, b_elements, x_elements
 
    def _snap_to_shannon(self, actual_cn):
        """
        Snaps a float/measured CN to the nearest available Shannon CN integer.
        """
        valid_cns = sorted(_CN_ROMAN.keys())
        return min(valid_cns, key=lambda x: abs(x - actual_cn))

    def _get_radius(self, element, oxi_map, cn):
        """
        Shannon ionic radius for one element, at its assigned oxidation
        state (from oxi_map) and the given coordination number.
        cn: int, 6 or 12 -- converted internally to the roman-numeral
        string pymatgen's table expects.
        """
        oxi_state = oxi_map[element]
        sp = Species(element, oxi_state)
        cn_roman = _CN_ROMAN[cn]

        try:
            return sp.get_shannon_radius(cn_roman)
        except KeyError:
            try:
                return sp.get_shannon_radius(cn_roman, spin="High Spin")
            except KeyError:
                if sp.ionic_radius is not None:
                    return sp.ionic_radius
                raise ValueError(f"No radius data found for {sp} (CN={cn_roman}).")

 
    def _avg_radius(self, elements, oxi_map, cn):
        """
        Mean ionic radius across one or more elements sharing a site
        (e.g. B and B' in a double perovskite). Averaging a single
        element is a no-op, so ABX3 and A2BB'X6/AA'BB'X6 use the exact
        same code path -- this IS the "modified Goldschmidt" mean-radius
        variant, not a separate formula.
        """
        radii = [self._get_radius(el, oxi_map, cn) for el in elements]
        return sum(radii) / len(radii)
 
    def _avg_oxi_state(self, elements, oxi_map):
        """Mean oxidation state across one or more elements (needed for
        Bartel's n_A term when the A-site is mixed)."""
        states = [oxi_map[el] for el in elements]
        return sum(states) / len(states)

    def _get_a_site_params(self, doc, a_elements, oxi_map):
        """
        Helper to resolve A-site coordination from the neighbor data and 
        compute the radius and oxidation state.
        """
        entry = self.ion_index.get(doc.get("material_id"), {})
        cn_data = entry.get("cn", {})

        # Extract measured CNs for A-site elements
        a_cns = [
            data["CN"] for site_idx, data in cn_data.items() 
            if data.get("element") in a_elements
        ]
        
        # Determine target CN
        if a_cns:
            avg_cn = self._snap_to_shannon(sum(a_cns) / len(a_cns))
        else:
            avg_cn = 12 # Default fallback

        r_A = self._avg_radius(a_elements, oxi_map, avg_cn)
        n_A = self._avg_oxi_state(a_elements, oxi_map)
        
        return r_A, n_A

    def get_tolerance_factors(self, doc):
        """
        Computes octahedral factor (mu), Goldschmidt tolerance factor (t),
        and Bartel factor (tau) for one material.
 
        Idealized coordination numbers are used for radius lookup
        (A: CN=12, B/X: CN=6) regardless of what CrystalNN measured on
        the real structure
 
        """
        self.oxi_index = self._load_oxi_index()
        self.ion_index = self._load_ion_assigments()
        oxi_map = self._get_oxi_map(doc)
        a_elements, b_elements, x_elements = self._get_ion_assignments(doc)
        if a_elements is None:
            print("Tolerance None")
            self.logger.error(f"{doc.get('material_id')} is not perovskite, returning None")
            return {
            "material_id": doc.get("material_id"),
            "a_elements": a_elements, "b_elements": b_elements, "x_elements": x_elements,
            "r_A": None, "r_B": None, "r_X": None, "n_A": None,
            "octahedral_factor": None,
            "goldschmidt_t": None,
            "bartel_tau": None,
        }

        try:
            # Modularized calculation steps
            r_A, n_A = self._get_a_site_params(doc, a_elements, oxi_map)
            r_B = self._avg_radius(b_elements, oxi_map, 6)
            r_X = self._avg_radius(x_elements, oxi_map, 6)
    
            # Compute Factors
            mu = r_B / r_X
            t = (r_A + r_X) / (math.sqrt(2) * (r_B + r_X))
            ratio = r_A / r_B
            tau = r_X / r_B - n_A * (n_A - ratio / math.log(ratio))

        except (KeyError, TypeError, ValueError, ZeroDivisionError) as e:
            print(f"Tolerance error: {e}")
            self.logger.error(f"Tolerance factor computation failed ({doc.get('material_id')}): {e!r}")
            return {
                "material_id": doc.get("material_id"),
                "a_elements": a_elements, "b_elements": b_elements, "x_elements": x_elements,
                "r_A": None, "r_B": None, "r_X": None, "n_A": None,
                "octahedral_factor": None,
                "goldschmidt_t": None,
                "bartel_tau": None,
            }

 
        return {
            "material_id": doc.get("material_id"),
            "a_elements": a_elements, "b_elements": b_elements, "x_elements": x_elements,
            "r_A": r_A, "r_B": r_B, "r_X": r_X, "n_A": n_A,
            "octahedral_factor": mu,
            "goldschmidt_t": t,
            "bartel_tau": tau,
        }