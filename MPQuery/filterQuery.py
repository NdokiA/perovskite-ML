import json, os, logging, time
from .checkTolerance import checkTolerance
from .checkNeighbor import checkNeighbor
from .checkOxidation import checkOxidation

class filterQuery():
    def __init__(self, QUERY_DIRS = "QUERY", LOG_DIRS="QUERY/LOG_QUERY/", DUMP_DIRS = "QUERY/DUMP"):
        
        self.cif_dirs = os.path.join(QUERY_DIRS, "CIF")
        self.query_path = os.path.join(QUERY_DIRS, "QUERY.json")
        self.oxi_path = os.path.join(QUERY_DIRS, "OXIDATION_QUERY.json")
        self.neigh_path = os.path.join(QUERY_DIRS, "NEIGHBOR_QUERY.json")
        self.tol_path = os.path.join(QUERY_DIRS, "TOLERANCE_QUERY.json")

        #logging directories
        os.makedirs(LOG_DIRS, exist_ok=True)
        self.log_path = os.path.join(LOG_DIRS, "queryFilter.log")

        #Dump directories
        self.cif_dump_dirs = os.path.join(DUMP_DIRS, "CIF")
        os.makedirs(self.cif_dump_dirs, exist_ok=True)

        self.query_dump_path = os.path.join(DUMP_DIRS, "QUERY.json")
        self.oxi_dump_path = os.path.join(DUMP_DIRS, "OXIDATION_QUERY.json")
        self.neigh_dump_path = os.path.join(DUMP_DIRS, "NEIGHBOR_QUERY.json")
        self.tol_dump_path = os.path.join(DUMP_DIRS, "TOLERANCE_QUERY.json")
        
        
        docs = self.load_json(self.query_path)
        self.docs = {doc["material_id"]: doc for doc in docs}
        self.oxi = {}
        self.neigh = {}
        self.tol = {}

        self.removed_docs = {}
        self.removed_oxi = {}
        self.removed_neigh = {}
        self.removed_tol = {}
        
        self._setup_logger()
            
    def load_json(self, path):
        with open(path, "r", encoding="utf-8") as f:
            return json.load(f)

    def save_json(self, metadata, path, update=True):

        if update and os.path.exists(path):
            with open(path, "r") as f:
                try:
                    results = json.load(f)
                    if not isinstance(results, list):
                        results = [results]
                except json.JSONDecodeError:
                    results = []
        else:
            results = []

        if isinstance(metadata, list):
            results.extend(metadata)
        else:
            results.append(metadata)

        with open(path, "w") as f:
            json.dump(results, f, indent=2)

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
    
    def _separate_doc(self, mpid, state, results_CO = None, results_CN = None, results_CT = None):
        match state:
            case "ehull":
                                
                #remove mpid doc
                removed_doc = self.docs.pop(mpid, None)
                removed_doc["excluded_reason"] = {"stage": state,
                                                  "reason": "Energy Above Hull Above 0.4 eV/atom"}
                self.removed_docs[mpid] = removed_doc
            case "oxidation":

                #remove mpid doc
                removed_doc = self.docs.pop(mpid, None)
                removed_doc["excluded_reason"] = {"stage": state,
                                                  "reason": "Oxidation state is not neutral"}
                self.removed_docs[mpid] = removed_doc
                self.removed_oxi[mpid] = results_CO

            case "neighbor":

                #remove mpid doc
                removed_doc = self.docs.pop(mpid, None)
                removed_doc["excluded_reason"] = {"stage": state,
                                                  "reason": "Neighboring CN state indicates no perovskite structure"}
                self.removed_docs[mpid] = removed_doc 
                self.removed_neigh[mpid] = results_CN

                #remove oxi doc
                removed_oxi = self.oxi.pop(mpid, None)
                self.removed_oxi[mpid] = removed_oxi      

            case "tolerance":

                #remove mpid doc
                removed_doc = self.docs.pop(mpid, None)
                removed_doc["excluded_reason"] = {"stage": state,
                                                  "reason": "Tolerance Factor Calculation Error!"}
                self.removed_docs[mpid] = removed_doc 
                self.removed_tol[mpid] = results_CT

                #remove oxi doc
                removed_oxi = self.oxi.pop(mpid, None)
                self.removed_oxi[mpid] = removed_oxi  

                #remove neigh doc 
                removed_neigh = self.neigh.pop(mpid, None)
                self.removed_neigh[mpid] = removed_neigh 

    def filter_Ehull(self):
        removed = 0
        self._log(f"Starting E_hull check")
        struct_start = time.perf_counter()

        for i, (mpid, doc) in enumerate(list(self.docs.items())):
            eHull = doc["energy_above_hull"] 
            if eHull > 0.4:
                self._log(f"{mpid} is NOT thermodynamicaly stable (Ehull = {eHull} eV/Atom)! Dumping data...")
                self._separate_doc(mpid, state="ehull")
                removed += 1

        total_elapsed = time.perf_counter() - struct_start
        remaining = len(self.docs)
        total = remaining + removed

        self._log(f"Finished Energy Above Hull check for all {total} MPIDs")
        self._log(
            f"Ehull filter: kept {remaining}/{total} "
            f"({remaining/total:.1%}), removed {removed}"
        )
        self._log(f"Total execution time: {total_elapsed:.3f} s")

    def filter_oxidation(self):
        
        removed = 0
        self._log(f"Starting oxidation-state")
        self.cO = checkOxidation(TEST_JSON_PATH=self.query_path,
                                 OUTPUT_JSON_PATH=self.oxi_path,
                                 logger = self.logger)
        
        struct_start = time.perf_counter()

        for i, (mpid, doc) in enumerate(list(self.docs.items())):

            results_cO = self._timed(
                f"Oxidation-state check ({mpid})",
                self.cO.check_charge,
                doc)

            if results_cO["is_neutral"]:
                self.oxi[mpid] = results_cO
            else:
                self._log(f"{mpid} is NOT neutral!. Dumping data...")
                self._separate_doc(mpid, state="oxidation", results_CO=results_cO)
                removed += 1
            
        total_elapsed = time.perf_counter() - struct_start
        remaining = len(self.docs)
        total = remaining + removed

        self._log(f"Finished oxidation-state calculation for all {total} MPIDs")
        self._log(
            f"Oxidation filter: kept {remaining}/{total} "
            f"({remaining/total:.1%}), removed {removed}"
        )
        self._log(f"Total execution time: {total_elapsed:.3f} s")
    
    def filter_neighbor(self):

        removed = 0
        self._log(f"Starting neighbor-state")

        oxidation_index = {d["material_id"]: d.get("ion_assignment", {}) for d in list(self.oxi.values())}
        self.cN = checkNeighbor(JSON_PATH=self.query_path,
                                OXI_PATH=self.oxi_path,
                                OUTPUT_JSON=self.neigh_path,
                                CIF_DIR=self.cif_dirs,
                                logger = self.logger,
                                OXI_INDEX=oxidation_index)

        self._log(f"Starting Neighbor Check")
        struct_start = time.perf_counter()

        for i, (mpid, doc) in enumerate(list(self.docs.items())):

            results_CN = self._timed(
                f"Neighbor (BX6) check ({mpid})",
                self.cN.verify_bx6,
                doc)

            if results_CN["is_perovskite"]:
                self.neigh[mpid] = results_CN
            else:
                self._log(f"{mpid} is NOT perovskite!. Dumping data...")
                self._separate_doc(mpid, state="neighbor", results_CN=results_CN)
                removed += 1
            
        total_elapsed = time.perf_counter() - struct_start
        remaining = len(self.docs)
        total = remaining + removed

        self._log(f"Finished neighbor (BX6) check calculation for all {total} MPIDs")
        self._log(
            f"Neighbor filter: kept {remaining}/{total} "
            f"({remaining/total:.1%}), removed {removed}"
        )
        self._log(f"Total execution time: {total_elapsed:.3f} s")
    
    def filter_tolerance(self):

        removed = 0
        self._log(f"Starting Tolerance Check")

        oxidation_index = list(self.oxi.values())
        ionic_index = list(self.neigh.values())
        self.cT = checkTolerance(JSON_PATH=self.query_path,
                                 OXI_PATH=self.oxi_path,
                                 NEIGH_PATH=self.neigh_path,
                                 OUTPUT_JSON=self.tol_path,
                                 logger = self.logger,
                                 OXI_INDEX=oxidation_index,
                                 ION_INDEX=ionic_index)

        
        struct_start = time.perf_counter()

        for i, (mpid, doc) in enumerate(list(self.docs.items())):

            results_CT = self._timed(
                f"Tolerance Check for ({mpid})",
                self.cT.get_tolerance_factors,
                doc)

            if results_CT["octahedral_factor"] is not None:
                self.tol[mpid] = results_CT
            else:
                self._log(f"Error occured at {mpid}. Aborting data...")
                self._separate_doc(mpid, "tolerance", results_CT = results_CT)
            
        total_elapsed = time.perf_counter() - struct_start
        remaining = len(self.docs)
        total = remaining + removed

        self._log(f"Finished tolerance check calculation for all {total} MPIDs")
        self._log(
            f"Tolerance filter: kept {remaining}/{total} "
            f"({remaining/total:.1%}), removed {removed}"
        )
        self._log(f"Total execution time: {total_elapsed:.3f} s")

    def save_all(self):
        self._log("Saving All Results...")

        self.save_json(list(self.docs.values()), self.query_path, update=False)
        self.save_json(list(self.oxi.values()), self.oxi_path, update=False)
        self.save_json(list(self.neigh.values()), self.neigh_path, update=False)
        self.save_json(list(self.tol.values()), self.tol_path, update=False)

        self.save_json(list(self.removed_docs.values()), self.query_dump_path, update=False)
        self.save_json(list(self.removed_oxi.values()), self.oxi_dump_path, update=False)
        self.save_json(list(self.removed_neigh.values()), self.neigh_dump_path, update=False)
        self.save_json(list(self.removed_tol.values()), self.tol_dump_path, update=False)

        self._log("Moving .cif Files...")
        for mpid in list(self.removed_docs.keys()):
            os.rename(os.path.join(self.cif_dirs, mpid+".cif"), 
                os.path.join(self.cif_dump_dirs, mpid+".cif"))
                
        self._log("Results Saved")
                

if __name__ == "__main__":
    fQ = filterQuery()
    fQ.filter_Ehull()
    fQ.filter_oxidation()
    fQ.filter_neighbor()
    fQ.filter_tolerance()
    fQ.save_all()