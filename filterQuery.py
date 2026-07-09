import json, os, logging, shutil
from checkTolerance import checkTolerance
from checkNeighbor import checkNeighbor

class filterQuery():
    def __init__(self, CIF_DIRS="QUERY/CIF", JSON_PATH="QUERY/QUERY.json",
                 OXI_PATH="QUERY/OXIDATION_QUERY.json", NEIGH_PATH="QUERY/NEIGHBOR_QUERY.json",
                 TOL_PATH="QUERY/TOLERANCE_QUERY.json", LOG_PATH="QUERY_FILTERED/LOG_QUERY/queryFilter.log"):
        
        self.cif_dirs = CIF_DIRS
        self.json_path = JSON_PATH
        self.oxi_path = OXI_PATH 
        self.neigh_path = NEIGH_PATH
        self.tol_path = TOL_PATH
        self.log_path = LOG_PATH

        self.make_output_paths()
        self.load_query()
        self._setup_logger()
        #Disregard these calls below if you're sure that tolerance query is filled
        #I need to rebuild the neighbor query because of a bug in tolerance query assumption
        #Not to mention I have filled the data
        #Erase this part later

        self.cT = self.build_tolerance_query()
        self.cN = self.build_neighbor_query()
            
    def load_json(self, path):
        with open(path, "r", encoding="utf-8") as f:
            return json.load(f)

    def save_json(self, metadata, path):
        """
        Save metadata to JSON
        """

        if os.path.exists(path):
            with open(path, "r") as f:
                try:
                    results = json.load(f)
                    if not isinstance(results, list):
                        results = [results]
                except json.JSONDecodeError:
                    results = []
        else:
            results = []
        
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

    def _convert(self, path):
        new_path = path.replace("QUERY/", "QUERY_FILTERED/", 1)

        if os.path.splitext(new_path)[1]:
            os.makedirs(os.path.dirname(new_path), exist_ok = True)
        else:
            os.makedirs(new_path, exist_ok = True)
        
        return new_path
    
    def make_output_paths(self):
        setattr(self, "cif_dirs_output", self._convert(self.cif_dirs))
        setattr(self, "json_path_output", self._convert(self.json_path))
        setattr(self, "oxi_path_output", self._convert(self.oxi_path))
        setattr(self, "neigh_path_output", self._convert(self.neigh_path))
        setattr(self, "tol_path_output", self._convert(self.tol_path))
    
    def load_query(self):
        json_queries = self.load_json(self.json_path)
        oxi_queries = self.load_json(self.oxi_path)
        neigh_queries = self.load_json(self.neigh_path)
        tol_queries = self.load_json(self.tol_path)
        
        self.query_IDs = [query["material_id"] for query in json_queries]

        self.json_lookup = {q["material_id"]: q for q in json_queries}
        self.oxi_lookup = {q["material_id"]: q for q in oxi_queries}
        self.neigh_lookup = {q["material_id"]: q for q in neigh_queries}
        self.tol_lookup = {q["material_id"]: q for q in tol_queries}
        
    def build_neighbor_query(self):
        cN = checkNeighbor(self.json_path, 
                           self.oxi_path, 
                           self.cif_dirs,
                           logger=self.logger)
        return cN
    
    def build_tolerance_query(self):
        cT = checkTolerance(self.json_path,
                            self.oxi_path,
                            self.neigh_path,
                            logger = self.logger)
        return cT

    def filter_queries(self):
        filtered_query = []
        filtered_oxi = []
        filtered_neigh = []
        filtered_tol = []
        
        print("Start Filtering")
        for id in self.query_IDs:
            json_query = self.json_lookup.get(id)
            oxi_query = self.oxi_lookup.get(id)
            neigh_query = self.neigh_lookup.get(id)
            tol_query = self.tol_lookup.get(id)
        
            if json_query["energy_above_hull"]<0.4 and neigh_query["is_perovskite"]:
                filtered_query.append(json_query)
                filtered_oxi.append(oxi_query)
                shutil.copy(os.path.join(self.cif_dirs, id+".cif"), 
                            os.path.join(self.cif_dirs_output, id+".cif"))

                neigh_query = self.cN.verify_bx6(json_query) #erase this once the process is completed
                filtered_neigh.append(neigh_query)

                if tol_query is None or tol_query.get("octahedral_factor") is None:
                    tol_query = self.cT.get_tolerance_factors(json_query)
                filtered_tol.append(tol_query)
        
        self.save_json(filtered_query, self.json_path_output)
        self.save_json(filtered_oxi, self.oxi_path_output)
        self.save_json(filtered_tol, self.tol_path_output)
        self.save_json(filtered_neigh, self.neigh_path_output)

if __name__ == "__main__":
    fQ = filterQuery()
    fQ.filter_queries()