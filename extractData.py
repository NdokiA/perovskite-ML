from config import mpAPI
from mp_api.client import MPRester
from emmet.core.summary import HasProps
import json, os

def find_full_metadata(chunk_size=1):    
    """
    Retrieve the complete summary metadata for all materials using
    `mpr.materials.summary.search()`.

    Parameters
    ----------
    chunk_size : int, optional
        Number of materials to process in each batch.  Default is 1.

    Returns
    -------
    Exported: cif_files/cif_file.cif
            : perovskite_metadata.json
    """
    with MPRester(mpAPI) as mpr:
        docs = mpr.materials.summary.search(
            chemsys="O-*-*",
            formula="ABC3",
            fields=[
                # Identifiers
                "material_id", "formula_pretty", "formula_anonymous", "chemsys",
                "elements", "nelements", "nsites",
                "deprecated", "theoretical", "database_IDs",

                # Structural / composition
                "composition", "composition_reduced", "symmetry",
                "volume", "density", "density_atomic", "structure",

                # Stability / energy
                "energy_per_atom", "uncorrected_energy_per_atom",
                "formation_energy_per_atom", "energy_above_hull",
                "is_stable", "equilibrium_reaction_energy_per_atom", "decomposes_to",

                # Electronic
                "band_gap", "is_gap_direct", "is_metal", "cbm", "vbm", "efermi",

                # Magnetic
                "is_magnetic", "ordering", "total_magnetization",

                # Mechanical
                "bulk_modulus", "shear_modulus", "universal_anisotropy", "homogeneous_poisson",
            ],
            chunk_size=chunk_size,
            num_chunks=None,
        )

    results = []
    cif_dir = "cif_files"
    os.makedirs(cif_dir, exist_ok=True)

    for doc in docs:
        entry = doc.model_dump()

        # Convert structure -> CIF, save separately, drop raw structure from JSON
        structure = entry.pop("structure", None)
        if structure is not None:
            try:
                cif_str = doc.structure.to(fmt="cif")
                cif_path = os.path.join(cif_dir, f"{doc.material_id}.cif")
                with open(cif_path, "w") as f:
                    f.write(cif_str)
                entry["cif_file"] = cif_path
            except Exception as e:
                entry["cif_file"] = None
                entry["cif_error"] = str(e)

        results.append(entry)

    with open("perovskite_metadata.json", "w") as f:
        json.dump(results, f, indent=2, default=str)

    print(f"Saved {len(results)} entries to perovskite_metadata.json")
    print(f"CIF files saved to ./{cif_dir}/")
    return results


if __name__ == "__main__":
    find_full_metadata(chunk_size=1)