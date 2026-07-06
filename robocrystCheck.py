from config import mpAPI
from mp_api.client import MPRester

with MPRester(mpAPI) as mpr:
    candidate_ids = {
        d.material_id for d in mpr.materials.summary.search(
            num_elements=(3, 5),
            energy_above_hull=(0, 0.4),
            fields=["material_id"],
        )
    }
print(len(candidate_ids))