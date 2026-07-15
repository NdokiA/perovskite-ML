import json, os, logging, time
from .checkTolerance import checkTolerance
from .checkNeighbor import checkNeighbor
from .checkOxidation import checkOxidation

class filterQuery():
    def __init__(self, CIF_DIRS="QUERY/CIF", JSON_PATH="QUERY/QUERY.json",
                 OXI_PATH="QUERY/OXIDATION_QUERY.json", NEIGH_PATH="QUERY/NEIGHBOR_QUERY.json",
                 TOL_PATH="QUERY/TOLERANCE_QUERY.json", LOG_PATH="QUERY/LOG_QUERY/queryFilter.log",
                 DUMP_DIRS = "QUERY/DUMP", DUMP_JSON="QUERY/DUMP.json"):
        
        self.cif_dirs = CIF_DIRS
        self.json_path = JSON_PATH
        self.oxi_path = OXI_PATH 
        self.neigh_path = NEIGH_PATH
        self.tol_path = TOL_PATH
        self.log_path = LOG_PATH
        self.dump_dirs = DUMP_DIRS
        self.json_dump = DUMP_JSON
        
        os.makedirs(self.dump_dirs, exist_ok=True)
        
        docs = self.load_json(self.json_path)
        self.docs = {doc["material_id"]: doc for doc in docs}
        
        self._setup_logger()
        #Disregard these calls below if you're sure that tolerance query is filled
        #I need to rebuild the neighbor query because of a bug in tolerance query assumption
        #Not to mention I have filled the data
        #Erase this part later
            
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
    
    def filter_Ehull(self):
        removed_docs = []
        self._log(f"Starting E_hull check")
        struct_start = time.perf_counter()

        for i, (mpid, doc) in enumerate(list(self.docs.items())):
            eHull = doc["energy_above_hull"] 
            if eHull > 0.4:
                self._log(f"{mpid} is NOT thermodynamicaly stable (Ehull = {eHull} eV/Atom)! Dumping data...")
                os.rename(os.path.join(self.cif_dirs, mpid+".cif"), 
                          os.path.join(self.dump_dirs, mpid+".cif"))
                removed_doc = self.docs.pop(mpid, None)
                removed_docs.append(removed_doc)

        self._log("Saving JSON Files...")
        self.save_json(removed_docs, self.json_dump)
        self.save_json(list(self.docs.values()), self.json_path, update=False)

        total_elapsed = time.perf_counter() - struct_start
        remaining = len(self.docs)
        removed = len(removed_docs)
        total = remaining + removed

        self._log(f"Finished Energy Above Hull check for all {total} MPIDs")
        self._log(
            f"Ehull filter: kept {remaining}/{total} "
            f"({remaining/total:.1%}), removed {removed}"
        )
        self._log(f"Total execution time: {total_elapsed:.3f} s")

    def filter_oxidation(self):
        
        self.cO = checkOxidation(TEST_JSON_PATH=self.json_path,
                                 OUTPUT_JSON_PATH=self.oxi_path,
                                 logger = self.logger)
        oxidation_docs = []
        removed_docs = []
        self._log(f"Starting oxidation-state")
        struct_start = time.perf_counter()

        for i, (mpid, doc) in enumerate(list(self.docs.items())):

            results_cO = self._timed(
                f"Oxidation-state check ({mpid})",
                self.cO.check_charge,
                doc)

            if results_cO["is_neutral"]:
                oxidation_docs.append(results_cO)
            else:
                self._log(f"{mpid} is NOT neutral!. Dumping data...")
                os.rename(os.path.join(self.cif_dirs, mpid+".cif"), 
                          os.path.join(self.dump_dirs, mpid+".cif"))
                removed_doc = self.docs.pop(mpid, None)
                removed_docs.append(removed_doc)
        
        self._log("Saving JSON Files...")
        self.save_json(removed_docs, self.json_dump)
        self.save_json(list(self.docs.values()), self.json_path, update=False)
        self.save_json(oxidation_docs, self.oxi_path, update=False)
            
        total_elapsed = time.perf_counter() - struct_start
        remaining = len(self.docs)
        removed = len(removed_docs)
        total = remaining + removed

        self._log(f"Finished oxidation-state calculation for all {total} MPIDs")
        self._log(
            f"Oxidation filter: kept {remaining}/{total} "
            f"({remaining/total:.1%}), removed {removed}"
        )
        self._log(f"Total execution time: {total_elapsed:.3f} s")
    
    def filter_neighbor(self):

        self.cN = checkNeighbor(JSON_PATH=self.json_path,
                                OXI_PATH=self.oxi_path,
                                OUTPUT_JSON=self.neigh_path,
                                CIF_DIR=self.cif_dirs,
                                logger = self.logger)
        removed_docs = []
        neighbor_docs = []

        oxidation_docs = {doc["material_id"]: doc for doc in self.load_json(self.oxi_path)}

        self._log(f"Starting Neighbor Check")
        struct_start = time.perf_counter()

        for i, (mpid, doc) in enumerate(list(self.docs.items())):

            results_CN = self._timed(
                f"Neighbor (BX6) check ({mpid})",
                self.cN.verify_bx6,
                doc)

            if results_CN["is_perovskite"]:
                neighbor_docs.append(results_CN)
            else:
                self._log(f"{mpid} is NOT perovskite!. Dumping data...")
                os.rename(os.path.join(self.cif_dirs, mpid+".cif"), 
                          os.path.join(self.dump_dirs, mpid+".cif"))
                oxidation_docs.pop(mpid, None)
                removed_doc = self.docs.pop(mpid, None)
                removed_docs.append(removed_doc)
        
        self._log("Saving JSON Files...")
        self.save_json(removed_docs, self.json_dump)
        self.save_json(list(oxidation_docs.values()), self.oxi_path, update=False)
        self.save_json(list(self.docs.values()), self.json_path, update=False)
        self.save_json(neighbor_docs, self.neigh_path, update=False)
            
        total_elapsed = time.perf_counter() - struct_start
        remaining = len(self.docs)
        removed = len(removed_docs)
        total = remaining + removed

        self._log(f"Finished neighbor (BX6) check calculation for all {total} MPIDs")
        self._log(
            f"Neighbor filter: kept {remaining}/{total} "
            f"({remaining/total:.1%}), removed {removed}"
        )
        self._log(f"Total execution time: {total_elapsed:.3f} s")
    
    def filter_tolerance(self):

        self.cT = checkTolerance(JSON_PATH=self.json_path,
                                 OXI_PATH=self.oxi_path,
                                 NEIGH_PATH=self.neigh_path,
                                 OUTPUT_JSON=self.tol_path,
                                 logger = self.logger)
        removed_docs = []
        tolerance_docs = []

        oxidation_docs = {doc["material_id"]: doc for doc in self.load_json(self.oxi_path)}
        neigh_docs = {doc["material_id"]: doc for doc in self.load_json(self.neigh_path)}

        self._log(f"Starting Tolerance Check")
        struct_start = time.perf_counter()

        for i, (mpid, doc) in enumerate(list(self.docs.items())):

            results_CT = self._timed(
                f"Tolerance Check for ({mpid})",
                self.cT.get_tolerance_factors,
                doc)

            if results_CT["octahedral_factor"] is not None:
                tolerance_docs.append(results_CT)
            else:
                self._log(f"Error occured at {mpid}. Aborting data...")
                os.rename(os.path.join(self.cif_dirs, mpid+".cif"), 
                          os.path.join(self.dump_dirs, mpid+".cif"))
                oxidation_docs.pop(mpid, None)
                neigh_docs.pop(mpid, None)
                removed_doc = self.docs.pop(mpid, None)
                removed_docs.append(removed_doc)
        
        self._log("Saving JSON Files...")
        self.save_json(removed_docs, self.json_dump)
        self.save_json(list(oxidation_docs.values()), self.oxi_path, update=False)
        self.save_json(list(neigh_docs.values()), self.neigh_path, update=False)
        self.save_json(list(self.docs.values()), self.json_path, update=False)
        self.save_json(tolerance_docs, self.tol_path, update=False)
            
        total_elapsed = time.perf_counter() - struct_start
        remaining = len(self.docs)
        removed = len(removed_docs)
        total = remaining + removed

        self._log(f"Finished tolerance check calculation for all {total} MPIDs")
        self._log(
            f"Tolerance filter: kept {remaining}/{total} "
            f"({remaining/total:.1%}), removed {removed}"
        )
        self._log(f"Total execution time: {total_elapsed:.3f} s")

if __name__ == "__main__":
    fQ = filterQuery()
    fQ.filter_Ehull()
    fQ.filter_oxidation()
    fQ.filter_neighbor()
    fQ.filter_tolerance()