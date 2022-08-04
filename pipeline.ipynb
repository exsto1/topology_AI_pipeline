import os
from urllib import request
from tqdm import tqdm
import multiprocessing
from topoly import alexander
import argparse


def download_Pfam(family, path):
    url1 = f"https://pfam.xfam.org/family/{family}/alignment/full/format?format=fasta&alnType=full&order=t&case=u&gaps=none&download=0"
    url2 = f"https://pfam.xfam.org/family/{family}/alignment/rp75/format?format=fasta&alnType=rp75&order=t&case=u&gaps=none&download=0"
    url3 = f"https://pfam.xfam.org/family/{family}/alignment/rp55/format?format=fasta&alnType=rp55&order=t&case=u&gaps=none&download=0"
    url4 = f"https://pfam.xfam.org/family/{family}/alignment/rp35/format?format=fasta&alnType=rp35&order=t&case=u&gaps=none&download=0"
    url5 = f"https://pfam.xfam.org/family/{family}/alignment/rp15/format?format=fasta&alnType=rp15&order=t&case=u&gaps=none&download=0"
    urls = [url1, url2, url3, url4, url5]
    id = ["full", "rp75", "rp55", "rp35", "rp15"]

    if not os.path.exists(path):
        os.makedirs(path)

    state = False
    for url in range(len(urls)):
        if not state:
            try:
                request.urlretrieve(urls[url], f"{path}/{family}_{id[url]}.fasta")
                state = True
            except:
                continue


def CD_HIT_wrapper(path):
    newpath = path.split(".")[0] + "_cluster.fasta"
    os.system(f"cdhit -i {path} -o {newpath} -c 0.95 -d 0")


def get_IDs_from_fasta(fasta):
    file_h = open(fasta)
    file = file_h.readlines()
    file_h.close()

    IDs = [i.lstrip(">").split("/")[0] for i in file if ">" in i]
    return IDs


def get_data_from_Uniprot(IDs):
    id_string = ""

    if "_" in IDs[0]:
        mode = "entry"
    else:
        mode = "accession"

    if mode == "entry":
        id_string = "%20OR%20".join([f"%28id%3A{i}%29" for i in IDs])
    elif mode == "accession":
        id_string = "%20OR%20".join([f"%28accession%3A{i}%29" for i in IDs])

    if id_string:
        url = f"https://rest.uniprot.org/uniprotkb/stream?fields=accession%2Cxref_pdb%2Cxref_alphafolddb&format=tsv&query={id_string}"
    else:
        print(f"incorrect id_string {IDs} / {mode}")
        return []

    data = request.urlopen(url).readlines()
    data = [i.decode("utf-8") for i in data]

    data = [i.rstrip().split("\t") for i in data[1:]]
    return data


def alpha_fold_cif_download(valid_alphafold_ids, FOLDER_PATH):
    if not os.path.exists(FOLDER_PATH):
        os.makedirs(FOLDER_PATH)

    for id in valid_alphafold_ids:
        filepath = f"{FOLDER_PATH}/{id}.cif"
        url = f"https://alphafold.ebi.ac.uk/files/AF-{id}-F1-model_v3.cif"
        request.urlretrieve(url, filepath)


def download_job_wrapper(batch):
    path, IDs = batch
    Alpha_IDs = get_data_from_Uniprot(IDs)
    Alpha_IDs = [i[2].rstrip(";") for i in Alpha_IDs if len(i) == 3]
    alpha_fold_cif_download(Alpha_IDs, path)


def simple_knot_calculation(file):
    """
    :param source: path to folder with PDB / cif files
    :param description: entry description (crystal / predicted / alphafold / other)
    """

    id = file.split("/")[-1].split(".")[0]
    try:
        knot = alexander(file)
        if "0_1" in knot:
            if knot["0_1"] > 0.5:
                res = [id, "0", f"Unknot-Alpha", str(knot)]
            else:
                tops = list(knot.keys())
                tops.sort(key=lambda i: knot[i], reverse=True)
                res = [id, tops[0], f"Knot-Alpha", str(knot)]
        else:
            tops = list(knot.keys())
            tops.sort(key=lambda i: knot[i], reverse=True)
            res = [id, tops[0], f"Knot-Alpha", str(knot)]
    except:
        res = [id, "?", f"?-Alpha", "Topoly-Error"]

    return res


def pipeline_families(families, mode="families"):
    batch_size = 100
    major_batch_size = 10
    major_batch_size_files = 10
    workers_n = 4
    file_workers_n = 3

    if mode == "clans":
        clans_map_h = open("family_clan.tsv")
        clans_map = clans_map_h.readlines()
        clans_map_h.close()

        clans_map = [i.rstrip().split("\t") for i in clans_map]

    for family in tqdm(families, desc="MAIN LOOP"):
        if mode == "families":
            base_path1 = f"DATA/{family}/family_seq"
            base_path2 = f"DATA/{family}/AlphaFold"
        elif mode == "clans":
            clan = ""
            for i in range(len(clans_map)):
                if clans_map[i][0] == family:
                    clan = clans_map[i][1]
            if not clan:
                clan = "No_Clan"

            base_path1 = f"DATA/{clan}/{family}/family_seq"
            base_path2 = f"DATA/{clan}/{family}/AlphaFold"

        download_Pfam(family, f"{base_path1}")
        CD_HIT_wrapper(f"{base_path1}/{family}_full.fasta")
        IDs = get_IDs_from_fasta(f"{base_path1}/{family}_full_cluster.fasta")
        batches = [[base_path2, IDs[i:i+batch_size]] for i in range(0, len(IDs), batch_size)]
        major_batches = [batches[i:i+major_batch_size] for i in range(0, len(batches), major_batch_size)]
        workers = multiprocessing.Pool(workers_n)

        for major_batch in tqdm(major_batches):
            workers.map(download_job_wrapper, major_batch)
        workers.close()
        workers.join()

        workers = multiprocessing.pool.ThreadPool(file_workers_n)
        files = os.listdir(f"{base_path2}")
        files = [f"{base_path2}/{i}" for i in files]
        file_batches = [files[i: i+major_batch_size_files] for i in range(0, len(files), major_batch_size_files)]

        results = []

        for file_batch in tqdm(file_batches):
            t = workers.map(simple_knot_calculation, file_batch)
            results.extend(t)

        workers.close()
        workers.join()

        print(f"FAMILY {family}")
        for i in results:
            print(i)

        outfile = open(f"{base_path1}/{family}_summary.tsv", "w")
        outfile.write("\n".join(["\t".join(i) for i in results]))
        outfile.close()


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", help="File with familiy IDs - 1 ID per line")
    args = parser.parse_args()

    input_data = args.input
    try:
        file_h = open(input_data)
    except:
        print("Wrong path!")
        parser.print_help()
        exit()

    file = file_h.readlines()
    file_h.close()

    families = [i.rstrip().upper() for i in file]
    pipeline_families(families, "clans")

