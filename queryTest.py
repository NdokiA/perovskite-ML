from MPQuery.testingQuery import queryPerovskite_Formula
from MPQuery.filterQuery import filterQuery
from lineGrab import load_json, save_json

formula = [
  "SrTiO3",
  "CaTiO3",
  "GdFeO3",
  "LaAlO3",
  "BaTiO3",
  "La2NiMnO6",
  "Sr2FeMoO6",
  "Ba2NdSbO6",
  "La2MnCoO6",
  "Nd2CoMnO6",
  "SrNdCoMnO6"
  "NdSbO3",
  "CsPbI3",
  "NaMgF3",
  "BaNiO3",
  "BaMnO3",
  "FeTiO3",
  "Ca3PbO",
  "MgSiO3",
  "ReO3",
  "WO3",
  "Sr2TiO4",
  "Sr3Ti2O7"
]

tQ = queryPerovskite_Formula(formula=formula,
                             CIF_DIRS = "QUERY_TRY/CIF",
                             JSON_PATH="QUERY_TRY/QUERY.json",
                             LOG_PATH = "QUERY_TRY/LOG_QUERY/queryFormula.log")

mpids = tQ.obtain_ID()
tQ.query_ID(mpids)

"""
fQ = filterQuery(QUERY_DIRS="QUERY_TRY",
                 LOG_DIRS="QUERY_TRY/LOG_QUERY",
                 DUMP_DIRS="QUERY_TRY/DUMP")
fQ.filter_Ehull()
fQ.filter_oxidation()
fQ.filter_neighbor()
fQ.save_all()
"""


