import os, logging, json, time
from config import mpAPI
from mp_api.client import MPRester
from queryStructure import queryStructure as qS
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
                 TOL_PATH="QUERY/TOLERANCE_QUERY.json", LOG_PATH="QUERY/LOG_QUERY/queryRobo.log",
                 BATCH_QUERY = 9999):
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

        for path in [JSON_PATH, OXI_PATH, NEIGH_PATH, TOL_PATH, LOG_PATH]:
            path.parent.mkdir(parents=True, exist_ok=True)
        os.makedirs(CIF_DIRS, exist_ok=True)
        
        self.log_path = LOG_PATH
        self._setup_logger()

        super().__init__(fields, CIF_DIRS, JSON_PATH, logger=self.logger)
        self.cN = cN(JSON_PATH, OXI_PATH, CIF_DIRS, NEIGH_PATH, logger = self.logger)
        self.cO = cO(TEST_JSON_PATH=JSON_PATH,
                     OUTPUT_JSON_PATH=OXI_PATH, logger = self.logger)
        self.cT = cT(JSON_PATH, OXI_PATH, NEIGH_PATH, TOL_PATH, logger = self.logger)

        self.NEIGH_PATH = NEIGH_PATH

        self.BATCH_QUERY = BATCH_QUERY
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
                chunk_size = 1000,
                num_chunks = None, # set None for production, 1 for testing
            )
            candidate_mpids = {d.material_id.string for d in robocrys_docs}
        
        self._log(f"Robocrys search returned {len(candidate_mpids)} candidate MPIDs")
 
        self._log(f"Total {len(candidate_mpids)} candidate MPIDs successfully obtained")
        return list(candidate_mpids)

    def batch_ID(self, mpids):
        batches = [mpids[i:i+self.BATCH_QUERY] 
                   for i in range(0, len(mpids), self.BATCH_QUERY)
                    ]
        return batches
    
    def query_ID(self, mpids):
        """
        Run the full verification pipeline (oxidation-state check, CrystalNN BX6
        connectivity check, tolerance-factor calculation) on each candidate MPID and
        save each result to its respective JSON output.
        """
        total_start = time.perf_counter()
        if not mpids:
            self._log("No candidate MPIDs provided. Aborting pipeline")
            return
        
        batches = self._timed(f"Batching queried structures", self.batch_ID, mpids)
        self._log(f"Querying in {len(batches)} batch(es) of up to {self.BATCH_QUERY} MPIDs")

        docs = []
        self._log(f"Starting verification pipeline for {len(mpids)} MPIDs")
        for batch_num, batch in enumerate(batches, start = 1):
            batch_docs = self._timed(
                f"Querying MPIDs (batch {batch_num}/{len(batches)})",
                self.query,
                batch
            )

            if batch_docs is None:
                self._log(
                    f"Batch {batch_num}/{len(batches)} failed -- "
                    f"skipping its {len(batch)} MPIDs"
                )
                continue
            docs.extend(batch_docs)


        self._timed("Processing queried structures", self.processing_query, docs)
        self._log(f"Finished Querying. Total downloaded files: {len(docs)} of requested {len(mpids)}")

        for doc in docs:
            mpid = doc["material_id"]
            struct_start = time.perf_counter()

            self._log(f"Starting verification pipeline for {mpid}")

            results_cO = self._timed(
                f"Oxidation-state check ({mpid})",
                self.cO.check_charge,
                doc,
            )
            if results_cO is not None:
                self.cO.save_json(results_cO)

            results_cN = self._timed(
                f"CrystalNN BX6 verification ({mpid})",
                self.cN.verify_bx6,
                doc,
            )
            if results_cN is not None:
                self.cN.save_json(results_cN)

            results_cT = self._timed(
                f"Tolerance factor computation ({mpid})",
                self.cT.get_tolerance_factors,
                doc,
            )
            if results_cT is not None:
                self.cT.save_json(results_cT)

            self._log(
                f"Completed verification pipeline for {mpid} "
                f"in {time.perf_counter() - struct_start:.3f} s"
            )

        total_elapsed = time.perf_counter() - total_start
        self._log(f"Finished verification pipeline for all {len(mpids)} MPIDs")
        self._log(f"Total Discovered Perovskite: {self.count_perovskite()}")
        self._log(f"Total execution time: {total_elapsed:.3f} s")

    def count_perovskite(self):
        """
        Report the number of succesfully acquired perovskite
        """
        with open(self.NEIGH_PATH, 'r') as f:
            data = json.load(f)

        count = sum(1 for item in data if item.get("is_perovskite") is True)
        return count

    def _timed(self, name, func, *args, **kwargs):
        """
        Execute a function, log its execution time, and return the result.
        """
        start = time.perf_counter()
        try:
            result = func(*args, **kwargs)
            elapsed = time.perf_counter() - start
            self._log(f"{name} completed in {elapsed:.3f} s")
            return result
        
        except Exception as e:
            elapsed = time.perf_counter() - start
            self._log(f"ERROR in {name} after {elapsed:.3f} s: {e}")
            return None
        
if __name__ == "__main__":
    query = queryPerovskite()
    IDs = query.obtain_ID()
    print(IDs)
    query.query_ID(IDs)


    
