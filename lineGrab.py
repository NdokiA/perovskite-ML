import re, os, json

def open_log(path):
    with open(path, "r", encoding="utf-8") as file:
        return [
            line
            for line in file
            if "Tolerance factor computation failed" in line
        ]

def write_log(log, path):
    with open(path, "w", encoding="utf-8") as file:
        file.write("".join(log))

def obtain_mpids(log):

    mpids = []
    for line in log:
        match = re.search(r"mp-\d+", line)
        if match:
            mpids.append(match.group())
    return mpids

def load_json(path):
        with open(path, "r", encoding="utf-8") as f:
            return json.load(f)

def save_json(metadata, path, update=True):

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

if __name__ == "__main__":
    os.makedirs("REQUERY", exist_ok=True)
    os.makedirs("REQUERY/CIF", exist_ok=True)

    log = open_log("requiresAttention.txt")
    #print("Writing Tolerance failure on requiresAttention.txt")
    #write_log(log, "requiresAttention.txt")

    print("Obtaining MPIDs")
    mpids = obtain_mpids(log)
    print(f"Number of Error MPIDS: {len(mpids)}")
    queries = {query["material_id"]:query for query in load_json("QUERY/DUMP.json")}
    print(f"Number of queries: {len(queries)}")
    discarded_queries = []

    list_cifs = {dir_[:-4] for dir_ in os.listdir("QUERY/DUMP") if dir_.endswith(".cif")}

    for mpid in mpids:
        if mpid in queries:
            discarded_queries.append(queries.pop(mpid))
            os.rename(f"QUERY/DUMP/{mpid}.cif", f"REQUERY/CIF/{mpid}.cif")

    print(len(discarded_queries))
    save_json(list(queries.values()), "QUERY/DUMP.json", update=False)
    save_json(list(discarded_queries), "REQUERY/QUERY.json")


