import os, logging
from config import mpAPI
from mp_api.client import MPRester
from emmet.core.symmetry import CrystalSystem
from queryTest import queryStructure as qS
from checkNeighbor import checkNeighbor as cN
from checkOxidation import checkOxidation as cO 
from checkTolerance import checkTolerance as cT
from pathlib import Path

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

class queryPerovskite(qS):
    def __init__(self, fields=FIELDS, CIF_DIRS="QUERY/CIF", JSON_PATH="QUERY/QUERY.json",
                 OXI_PATH="QUERY/OXIDATION_QUERY.json", NEIGH_PATH="QUERY/NEIGHBOR_QUERY.json",
                 TOL_PATH="QUERY/TOLERANCE_QUERY.json", LOG_PATH="QUERY/LOG_QUERY/query.log"):
        """
        Query and structurally verify candidate perovskite structures from Materials Project.

        Candidates are first identified via a robocrystallographer keyword search unioned
        with a provenance tags/remarks search (see `obtain_ID`), then each candidate is
        verified via oxidation-state, CrystalNN connectivity, and tolerance-factor checks

        Parameters
        ----------
        fields : list(str), optional
            Fields of the projected output for MP queries on each structure, default is FIELDS.
        CIF_DIRS : str, optional
            Directory path for CIF files of the queried structure.
        JSON_PATH : str, optional
            File path for the general query JSON output.
        OXI_PATH : str, optional
            File path for the oxidation-state check JSON output.
        NEIGH_PATH : str, optional
            File path for the CrystalNN neighbor/connectivity JSON output.
        TOL_PATH : str, optional
            File path for the tolerance-factor JSON output.
        LOG_PATH : str, optional
            File path for the run log (parent directory is created automatically if
            it doesn't exist).
        """

        CIF_DIRS = Path(CIF_DIRS)
        JSON_PATH = Path(JSON_PATH)
        OXI_PATH = Path(OXI_PATH)
        NEIGH_PATH = Path(NEIGH_PATH)
        TOL_PATH = Path(TOL_PATH)
        LOG_PATH = Path(LOG_PATH)

        for path in [CIF_DIRS, JSON_PATH, OXI_PATH, NEIGH_PATH, TOL_PATH, LOG_PATH]:
            path.parent.mkdir(parents=True, exist_ok=True)
        
        super().__init__(fields, CIF_DIRS, JSON_PATH)
        self.cN = cN(JSON_PATH, OXI_PATH, CIF_DIRS, NEIGH_PATH)
        self.cO = cO(TEST_JSON_PATH=JSON_PATH,
                     OUTPUT_JSON_PATH=OXI_PATH)
        self.cT = cT(JSON_PATH, OXI_PATH, NEIGH_PATH, TOL_PATH)
        self.log_path = LOG_PATH
        self._setup_logger()

    def _setup_logger(self):
        """
        Configure a logger that writes timestamped action logs to `self.log_path` and
        echoes them to console.

        Uses a unique logger name per instance so that multiple `queryPerovskite`
        instances pointed at different `LOG_PATH`s don't accumulate each other's
        file handlers.
        """
        log_dir = os.path.dirname(self.log_path)
        if log_dir:
            os.makedirs(log_dir, exist_ok=True)

        self.logger = logging.getLogger(f"{__name__}.queryPerovskite.{id(self)}")
        self.logger.setLevel(logging.INFO)
        self.logger.propagate = False

        if not self.logger.handlers:
            formatter = logging.Formatter("%(asctime)s\t\t%(message)s", datefmt="%Y-%m-%d %H:%M:%S")

            file_handler = logging.FileHandler(self.log_path)
            file_handler.setFormatter(formatter)
            self.logger.addHandler(file_handler)

            stream_handler = logging.StreamHandler()
            stream_handler.setFormatter(formatter)
            self.logger.addHandler(stream_handler)

    def _log(self, action):
        """
        Log a single action with a timestamp, formatted as `time<TAB><TAB>action`.

        Parameters
        ----------
        action : str
            Description of the action being logged.
        """
        self.logger.info(action)

    def obtain_ID(self):
        """
        Query candidate perovskite MPIDs from Materials Project.
 
        Combines two independent candidate sources so that structures missed by one
        are still caught by the other:
          - robocrystallographer descriptions containing the keyword "perovskite"
          - provenance tags/remarks containing the keyword "perovskite"
 
        Returns
        -------
        list(str)
            Union of MPIDs from both sources.
        """
        self._log("Starting robocrys keyword search for 'perovskite'")

        with MPRester(mpAPI) as mpr:
            robocrys_docs = mpr.materials.robocrys.search(
                keywords=["perovskite"],
                chunk_size=1, num_chunks=1,  # set None for production, 1 for testing
            )
            robo_ids = {d.material_id for d in robocrys_docs}
        self._log(f"Robocrys search returned {len(robo_ids)} candidate MPIDs")
 
        self._log("Narrowing candidates by crystal system before provenance search")
        candidate_ids = set()
        with MPRester(mpAPI) as mpr_summary:
            summary_docs = mpr_summary.materials.summary.search(
                num_elements=(3, 5),
                energy_above_hull=(0, 0.4),
                fields=["material_id"],
                chunk_size=1, num_chunks=1,  # set None for production, 1 for testing
            )
            candidate_ids.update(d.material_id for d in summary_docs)

        self._log(f"Summary pre-filter narrowed to {len(candidate_ids)} material IDs ")
 
        self._log("Starting provenance tags/remarks search for 'perovskite'")
        with MPRester(mpAPI, use_document_model=False) as mpr2:
            prov_docs = mpr2.materials.provenance.search(
                material_ids=list(candidate_ids),
                fields=["material_id", "remarks", "tags"],
                chunk_size=1, num_chunks=1,  # set None for production, 1 for testing
            )

        prov_ids = {
            doc.get("material_id") for doc in prov_docs
            if any("perovskite" in t.lower() for t in (doc.get("tags", []) + doc.get("remarks", [])))
        }
        self._log(f"Provenance search returned {len(prov_ids)} candidate MPIDs")
 
        candidate_mpids = list(robo_ids | prov_ids)
        self._log(f"Total {len(candidate_mpids)} candidate MPIDs successfully obtained (union of both sources)")
        return candidate_mpids


    def query_ID(self, mpids):
        """
        Run the full verification pipeline (oxidation-state check, CrystalNN BX6
        connectivity check, tolerance-factor calculation) on each candidate MPID and
        save each result to its respective JSON output.

        Parameters
        ----------
        mpids : list(str)
            Candidate MPIDs to verify, typically the output of `obtain_ID`.
        """
        self._log(f"Starting verification pipeline for {len(mpids)} MPIDs")
        for id in mpids:
            self._log(f"Querying structure for {id}")
            try:
                doc = self.query(id)[0]
            except Exception as e:
                self._log(f"ERROR querying {id}: {e}")
                continue

            self._log(f"Running oxidation-state check for {id}")
            results_cO = self.cO.check_charge(doc)
            self.cO.save_json(results_cO)

            self._log(f"Running CrystalNN BX6 verification for {id}")
            results_cN = self.cN.verify_bx6(doc)
            self.cN.save_json(results_cN)

            self._log(f"Computing tolerance factors for {id}")
            results_cT = self.cT.get_tolerance_factors(doc)
            self.cT.save_json(results_cT)

            self._log(f"Completed verification pipeline for {id}")
        self._log(f"Finished verification pipeline for all {len(mpids)} MPIDs")

if __name__ == "__main__":
    query = queryPerovskite()
    IDs = query.obtain_ID()
    print(IDs)


    
