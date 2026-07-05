from itertools import combinations
from pymatgen.core import Structure, Species
from pymatgen.analysis.local_env import CrystalNN
from pymatgen.analysis.bond_valence import BVAnalyzer
import os, json

class checkNeighbor():
    def __init__(self, JSON_PATH = "TEST_QUERY", OXI_PATH = "OXIDATION_QUERY", CIF_DIR = "TEST_CIF"):
        """
        Check neighboring BX6 structure in .cif structure to determine perovskite order
        cif_path :   .cif document
        """
        self.JSON_PATH = JSON_PATH + ".json"
        self.OXI_PATH = OXI_PATH + ".json" if OXI_PATH else None
        self.CIF_DIR = CIF_DIR
        self.oxi_index = self._load_oxi_index() if OXI_PATH else {}

    def _load_json(self):
        with open(self.JSON_PATH, "r", encoding="utf-8") as f:
            return json.load(f)
    
    def _load_oxi_index(self):
        """index checkoxidation results by material_id"""
        with open(self.OXI_PATH, "r", encoding="utf-8") as f:
            data = json.load(f)
        return {d["material_id"]: d.get("ion_assignment", {}) for d in data}

    def _get_oxi_map(self, doc):
        """oxidation dict for this material, or None if oxidation_query available"""
        return self.oxi_index.get(doc.get("material_id"))

    def _extract_cif(self,doc, oxi_map = None):
        """
        doc is extracted from self._load_json()
        """
        if doc.get("cif_file"):
            cif_path = os.path.join(self.CIF_DIR, os.path.basename(doc["cif_file"]))
        else:
            cif_path = os.path.join(self.CIF_DIR, doc["material_id"] + ".cif")
        
        structure = Structure.from_file(cif_path)
        if oxi_map:
            structure.add_oxidation_state_by_element(oxi_map)
        else:
            #fallback using BVAnalyzer
            try:
                valences = BVAnalyzer().get_valences(structure)
                structure.add_oxidation_state_by_site(valences)
            except Exception:
                pass
        return structure

    def _get_candidates(self, doc):
        """
        From possible_species oxidation states:
          x_elements   = all anions  (oxi < 0)
          b_candidates = all cations (oxi > 0)  <-- every cation is a B candidate
        CN=6 test later decides which are real octahedral B sites.
        Fallback if oximap unavailable
        """
        parsed = [Species.from_str(s) for s in doc.get("possible_species", [])]
        x_elements = sorted({sp.symbol for sp in parsed if sp.oxi_state < 0})
        b_candidates = sorted({sp.symbol for sp in parsed if sp.oxi_state > 0})
        return x_elements, b_candidates

    def _get_candidates_oximap(self, oxi_map):
        """X = anions (oxi<0), B-candidates = all cations (oxi>0)"""
        x_elements = sorted(el for el, ox in oxi_map.items() if ox < 0)
        b_candidates = sorted(el for el, ox in oxi_map.items() if ox > 0)
        return x_elements, b_candidates

    def _get_element_indices(self, structure, elements):
        """
        obtain site indices from element symbols 
        (ordered + disordered aware)
        """
        elements = set(elements)
        indices = []
        for i, site in enumerate(structure):
            if site.is_ordered:
                if site.specie.symbol in elements:
                    indices.append(i)
            else:
                if any(sp.symbol in elements for sp in site.species):
                    indices.append(i)
        return indices
    
    def _get_site_cn(self, structure, site_idx, anion_indices, cnn=None):
        """
        CN for one site (site B, anion-restricted), periodic-image safe
        """
        cnn = cnn or CrystalNN()
        nn = cnn.get_nn_info(structure, site_idx)
        #key by (site_index, image) so two different periodic copies of the 
        # same crystallographic site are not collapsed into one neighbor
        anion_nn = {
            (n['site_index'], n['image'])
            for n in nn if n['site_index'] in anion_indices
        }
        return anion_nn
    
    def _build_polyhedra_map(self, structure, cation_indices, anion_indices, cnn = None):
        cnn = cnn or CrystalNN()
        return {ci: self._get_site_cn(structure, ci, anion_indices, cnn)
                for ci in cation_indices}
    
    def _check_cn(self, poly_map, cn_target=6, tol=0):
        """
        CN pass/fail check (default tol=0, strict for BX6)
        """
        return {ci: {"CN": len(anions), 
                     "CN_ok": abs(len(anions) - cn_target)<=tol
                     }
                     for ci, anions in poly_map.items()
                     }
    def _get_shared_anions(self, poly_map, i, j):
        return poly_map[i] & poly_map[j]
    
    def _classify_share_mode(self, n_shared):
        return {0: 'none', 1: 'corner', 2: 'edge', 3: 'face'
                }.get(n_shared, 'over-face')
    
    def _site_label(self, structure, idx):
        """
        helper: species label for a site, ordered or disordered
        """
        site = structure[idx]
        if site.is_ordered:
            return site.specie.symbol
        return max(site.species.items(), key=lambda kv:kv[1])[0].symbol
    
    def _check_pairwise_connectivity(self, structure, poly_map):
        pairs = []
        for i, j in combinations(poly_map.keys(), 2):
            shared = self._get_shared_anions(poly_map, i, j)
            n = len(shared)
            if n == 0:
                continue
            pairs.append({
                "i": i, "j": j, "n_shared": n,
                "mode": self._classify_share_mode(n),
                "hetero": self._site_label(structure, i) != self._site_label(structure, j)
            })
        return pairs

    def is_perovskite(self, cn_results, conn_results, require_corner=False):
        """
        Perovskite iff >=1 site has CN==6 AND those octahedra form a
        connected network. require_corner -> demand corner-sharing link.
        """
        b_sites = [ci for ci, r in cn_results.items() if r["CN_ok"]]
        if not b_sites:
            return False, b_sites
        if require_corner:
            has_net = any(p["mode"] == "corner"
                          and p["i"] in b_sites and p["j"] in b_sites
                          for p in conn_results)
        else:
            has_net = any(p["i"] in b_sites and p["j"] in b_sites
                          for p in conn_results)
        return has_net, b_sites

    
    def verify_bx6(self, doc, x_elements = None, b_elements = None, cn_target=6, cn_tol=0,
                   require_corner = False):
        oxi_map = self._get_oxi_map(doc)
        structure = self._extract_cif(doc, oxi_map)

        oxi_map = self._get_oxi_map(doc)
        # candidate elements: prefer oxi_map, else possible_species
        if x_elements is None or b_elements is None:
            if oxi_map:
                x_elements, b_elements = self._get_candidates_oximap(oxi_map)
            else:
                x_elements, b_elements = self._get_candidates(doc)

        x_idx = self._get_element_indices(structure, x_elements)
        b_idx = self._get_element_indices(structure, b_elements)
        cnn = CrystalNN()
        poly_map = self._build_polyhedra_map(structure, b_idx, x_idx, cnn)
        cn_results = self._check_cn(poly_map, cn_target, cn_tol)
        conn_results = self._check_pairwise_connectivity(structure, poly_map)

        perov, b_sites = self.is_perovskite(cn_results, conn_results, require_corner)
 
        # true B elements = elements actually at CN=6 sites
        b_true = sorted({self._site_label(structure, ci) for ci in b_sites})
 
        corner_frac = (sum(p['mode'] == 'corner' for p in conn_results)
                       / len(conn_results) if conn_results else 0.0)
        n_corner = sum(p['mode'] == "corner" for p in conn_results)
        hetero_corner_frac = (sum(p["mode"] == "corner" and p["hetero"] for p in conn_results)
                              / n_corner if n_corner else None)
 
        return {
            "material_id": doc.get("material_id"),
            "is_perovskite": perov,
            "x_elements": x_elements,
            "b_candidates": b_elements,
            "b_true": b_true,          # elements at real CN=6 octahedral sites
            "n_b_sites": len(b_sites),
            "cn": cn_results,
            "connectivity": conn_results,
            "corner_share_fraction": corner_frac,
            "hetero_corner_fraction": hetero_corner_frac,
        }


def write_log(log, filename = "neighbor.log"):
    mode = "w" if not os.path.isfile(filename) else "a"
    prefix = "" if mode == "w" else "\n"
    with open(filename, mode, encoding="utf-8") as file:
        file.write(prefix + log)

if __name__ == "__main__":
    cN = checkNeighbor(JSON_PATH="TEST_QUERY", CIF_DIR="TEST_CIF")
    if os.path.exists("neigbor.log"):
        os.remove("neighbor.log")
    queries = cN._load_json()
    for query in queries:
        try:
            result = cN.verify_bx6(query)
        except Exception as e:
            print(f"SKIP {query.get('material_id')} -> {e}")
            continue
        log = (f"{result['material_id']} | perovskite={result['is_perovskite']} "
               f"| B={result['b_true']} | n_B={result['n_b_sites']} "
               f"| corner_frac={result['corner_share_fraction']:.2f}")
        write_log(log)
        print(log)
