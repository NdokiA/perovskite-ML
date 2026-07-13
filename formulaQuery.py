import os, logging, json, time
from config import mpAPI
from robocrystQuery import queryPerovskite as qP
from mp_api.client import MPRester
from queryStructure import queryStructure as qS
from checkNeighbor import checkNeighbor as cN
from checkOxidation import checkOxidation as cO 
from checkTolerance import checkTolerance as cT
from pathlib import Path

def load_json(path):
    with open(path, "r", encoding="utf-8") as f:
        return json.load(f)
    
class queryPerovskite_Formula(qP):
    def __init__(self, formula = None, LOG_PATH="QUERY/LOG_QUERY/queryFormula.log",
                 BATCH_QUERY = 9999):
        """
        Query and structurally verify candidate perovskite structures from Materials Project.

        Candidates are first identified based on formula_anonymous keyword and various others
        """
        
        super().__init__(LOG_PATH=LOG_PATH, BATCH_QUERY=BATCH_QUERY)
        files = load_json(self.json_path) if os.path.exists(self.json_path) else None
        self.formula = formula 
        self.queried = {}
        if files is not None: #prioritize if previous robocrystQuery is ran
            self.formula = set([file["formula_anonymous"] for file in files])
            self.queried = set([file["material_id"] for file in files])
        

    def obtain_ID(self):
        """
        Query candidate perovskite MPIDs manually based on the formula_anonymous
        from previous run (robocryst Query) or initialized in class

        Returns
        -------
        list(str)
            Union of MPIDs
        """
        self._log(f"Starting formula keyword search for formula: {self.formula}")

        with MPRester(mpAPI) as mpr:
            formula_docs = mpr.materials.summary.search(
                formula  = list(self.formula),
                energy_above_hull = (0,0.4),
                fields = ["material_id"],
                chunk_size = 10, #set to 100 for prod, 1 for test
                num_chunks = 1, #set to None for prod, 1 for test

            )

            candidate_mpids = {d.material_id.string for d in formula_docs} - self.queried
        
        self._log(f"Formula search returned {len(candidate_mpids)} candidate MPIDs")
 
        self._log(f"Total {len(candidate_mpids)} candidate MPIDs successfully obtained")
        return list(candidate_mpids)


if __name__ == "__main__":
    query = queryPerovskite_Formula()
    IDs = query.obtain_ID()
    if not IDs:
        raise ValueError("query() called with an empty list of material IDs. Nothing to query.")
    print(IDs)
    query.query_ID(IDs)
