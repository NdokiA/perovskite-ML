from itertools import combinations
from pymatgen.core import Structure, Species
from pymatgen.analysis.local_env import CrystalNN
from pymatgen.analysis.bond_valence import BVAnalyzer
import os, json

class checkNeighbor():
    def __init__(self, JSON_PATH = "TEST_QUERY.json", OXI_PATH = "OXIDATION_QUERY.json", CIF_DIR = "TEST_CIF",
                 OUTPUT_JSON = "NEIGHBOR_QUERY.json"):
        """
        Given a crystal structure (.cif), checks whether its cation sites
        form BX6 octahedra (6 anion neighbors) and whether those octahedra
        connect to each other by sharing corners — the defining feature of
        a perovskite-type structure.

        JSON_PATH : str, optional
            Filename of the query file listing materials
            to check (each entry needs at least "material_id", and either
            "cif_file" or a CIF named after its material_id).

        OXI_PATH : str, optional
            Filename of oxidation-state results from
            checkOxidation.py, keyed by material_id. Pass None to skip
            loading this and rely on possible_species / BVAnalyzer instead.

        CIF_DIR : str, optional
            Folder containing the .cif files referenced by JSON_PATH.
        """

        self.JSON_PATH = JSON_PATH
        self.OXI_PATH = OXI_PATH if OXI_PATH else None
        self.CIF_DIR = CIF_DIR
        self.OUTPUT_PATH = OUTPUT_JSON

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
        """index checkoxidation results by material_id"""
        with open(self.OXI_PATH, "r", encoding="utf-8") as f:
            data = json.load(f)
        return {d["material_id"]: d.get("ion_assignment", {}) for d in data}

    def _get_oxi_map(self, doc):
        """Oxidation-state dict for this material (element -> charge), or
        None if this material_id has no entry in the loaded oxi_index."""

        return self.oxi_index.get(doc.get("material_id"))

    def _extract_cif(self,doc, oxi_map = None):
        """
        Loads the structure from its .cif file and, if possible, tags each
        atom with an oxidation state

        Priority: use oxi_map if given (from checkOxidation.py results).
        Otherwise, try guessing oxidation states automatically via
        BVAnalyzer (bond-valence method)
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
        Fallback if oximap is unavailable
        """
        parsed = [Species.from_str(s) for s in doc.get("possible_species", [])]
        x_elements = sorted({sp.symbol for sp in parsed if sp.oxi_state < 0})
        b_candidates = sorted({sp.symbol for sp in parsed if sp.oxi_state > 0})
        return x_elements, b_candidates

    def _get_candidates_oximap(self, oxi_map):
        """
        reads from an oxi_map (element-> oxidation state) 
        Anions become x_elements; cations become b_candidates.
        """

        x_elements = sorted(el for el, ox in oxi_map.items() if ox < 0)
        b_candidates = sorted(el for el, ox in oxi_map.items() if ox > 0)
        return x_elements, b_candidates

    def _get_element_indices(self, structure, elements):
        """
        Finds which sites in the structure hold any of the given elements
        and returns their index numbers.

        Handles both ordered sites (one element per site) and disordered
        sites (partial occupancy, mixed elements) for a disordered
        site
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
        Finds the anion neighbors of one cation site using CrystalNN, and
        returns them as a set of (site_index, image) pairs. The count of
        this set is the site's coordination number (CN).

        Using (site_index, image) instead of just site_index because
        two neighbors can be the SAME crystallographic site but different
        periodic copies
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
        """
        Runs _get_site_cn for every candidate cation site and collects the
        results into one dict: {cation_site_index: set of anion neighbors}.
        """

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
    
    def _classify_share_mode(self, n_shared):
        """

        Turns a shared-anion count into the usual crystallography name:
        1 shared anion = corner-sharing, 
        2 = edge-sharing, 
        3 = face-sharing,
        0 = not touching
        >3 = unusual — labeled "over-face", maybe geometry issue, double-check
        """

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
    
    def _find_edges(self, poly_map, i, j):

        """
        Find all distinct periodic 'edges' (specific corner/edge/face-sharing
        contacts) between the home-cell polyhedron at site i and translated
        copies of the polyhedron at site j.

        Returns dict: T (image offset) -> shared anion set for that specific
        instance pair. Necessary because a raw poly_map[i] & poly_map[j]
        only checks the home (T=0) instances against each other, missing:
        self-sharing and cross-site sharing across different unit cell
        """

        by_anion_j = {}

        for (a, img) in poly_map[j]:
            by_anion_j.setdefault(a, []).append(img)

        edges = {}
        for (a, imgA) in poly_map[i]:
            for imgB in by_anion_j.get(a, []):
                T = tuple(round(x - y) for x, y in zip(imgA, imgB))
                if i == j:
                    if T == (0, 0, 0):
                        continue  # same instance, not a neighbor

                    if T > tuple(-x for x in T):
                        continue  # dedupe: T and -T are the same physical bond

                edges.setdefault(T, set()).add((a, imgA))
        return edges



    def _check_pairwise_connectivity(self, structure, poly_map):
        """
        Goes through every pair of cation sites (including a site paired
        with itself) and records each touching point found by
        _find_edges as one entry: which two sites, how many anions they
        share, what kind of sharing it is (corner/edge/face), and
        whether the two sites are different elements ("hetero").
        """

        pairs = []
        keys = list(poly_map.keys())

        for idx_a in range(len(keys)):
            for idx_b in range(idx_a, len(keys)):
                i, j = keys[idx_a], keys[idx_b]
                edges = self._find_edges(poly_map, i, j)

                for T, shared in edges.items():

                    n = len(shared)
                    if n == 0:
                        continue

                    pairs.append({
                        "i": i, "j": j, "image": T, "n_shared": n,
                        "mode": self._classify_share_mode(n),
                        "hetero": self._site_label(structure, i) != self._site_label(structure, j)
                    })
        return pairs

    def is_perovskite(self, cn_results, conn_results,
                    require_corner=True, corner_frac_threshold=0.9):
        """
        True/False call: a structure counts as perovskite if at least one
        site has CN==6 (a real BX6 octahedron) AND the octahedral network
        connecting those confirmed B sites is corner-sharing OVERWHELMINGLY
        (>= corner_frac_threshold of B-B contacts), not just partially.
        """
        b_sites = [ci for ci, r in cn_results.items() if r["CN_ok"]]
        if not b_sites:
            return False, b_sites, 0.0

        b_conn = [p for p in conn_results if p["i"] in b_sites and p["j"] in b_sites]
        if not b_conn:
            return False, b_sites, 0.0

        corner_frac = sum(p["mode"] == "corner" for p in b_conn) / len(b_conn)

        if require_corner:
            has_net = corner_frac >= corner_frac_threshold
        else:
            has_net = True

        return has_net, b_sites, corner_frac

    
    def verify_bx6(self, doc, x_elements = None, b_elements = None, cn_target=6, cn_tol=1,
                   require_corner = True):
        """

        Main entry point. Runs the full check on one material:

          1. Load structure, tag it with oxidation states.

          2. Figure out which elements are anions (X) and which are
             cation candidates (B) — every cation counts as a candidate
             at this stage, not just the "real" B-site ones.

          3. Build each cation's anion-neighbor set (poly_map) and check
             which ones actually hit CN==6 (the real octahedral sites).

          4. Check how those octahedra connect to each other (corner /
             edge / face sharing), including through periodic images.

          5. Decide is_perovskite, and report sharing stats restricted to
             the confirmed CN==6 sites only (so A-site cations, which
             aren't 6-coordinate, don't distort the numbers).

        require_corner: passed straight to is_perovskite 

        Returns a dict with the verdict plus all the intermediate data
        (cn, connectivity, corner_share_fraction, hetero_corner_fraction)
        for inspection.
        """       
        self.oxi_index = self._load_oxi_index() if self.OXI_PATH else {}

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

        perov, b_sites, corner_frac = self.is_perovskite(cn_results, conn_results, require_corner)
 
        # true B elements = elements actually at CN=6 sites
        b_true = sorted({self._site_label(structure, ci) for ci in b_sites})
 
        # restrict the sharing stats to edges between CONFIRMED B sites only.
        # conn_results still includes every cation candidate (A-site included,
        # since b_elements/b_idx cover all cations by design). Without this
        # filter, A-site cations (CN~12, not octahedral) leak into the
        # corner-sharing stats and dilute them.

        b_conn = [p for p in conn_results if p['i'] in b_sites and p['j'] in b_sites]

        n_corner = sum(p['mode'] == "corner" for p in b_conn)
        hetero_corner_frac = (sum(p["mode"] == "corner" and p["hetero"] for p in b_conn)
                              / n_corner if n_corner else None)

        return {
            "material_id": doc.get("material_id"),
            "formula_pretty": doc.get("formula_pretty"),
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