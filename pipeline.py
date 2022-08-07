import os
from urllib import request
from tqdm.auto import tqdm
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
    newpath = "./" + path.split(".")[1] + "_cluster.fasta"
    os.system(f"./cd-hit -i {path} -o {newpath} -c 0.95 -d 0 > /dev/null")


def get_IDs_from_cluster(clstr):
    file_h = open(clstr)
    file = file_h.readlines()
    file_h.close()
    IDs = []

    temp = []
    for i in file:
        if ">Cluster" in i:
            if temp:
                IDs.append(temp)
                temp = []
        else:
            temp.append(i.split(">")[1].split("/")[0])

    return IDs


def get_data_from_Uniprot(IDs, mode):
    id_string = ""
    string_list = []
    data_full = []

    for i in IDs:
        for i1 in i:
            if "_" in i1:
                string_list.append(f"%28id%3A{i1}%29")
            else:
                string_list.append(f"%28accession%3A{i1}%29")

    if string_list:
        for i in range(0, len(string_list), 100):
            id_string = "%20OR%20".join(string_list[i:i+100])
            url = f"https://rest.uniprot.org/uniprotkb/stream?fields=accession%2Cid%2Cxref_pdb%2Cxref_alphafolddb&format=tsv&query={id_string}"
            data = request.urlopen(url).readlines()
            data = [i.decode("utf-8") for i in data]
            data = [i.rstrip().split("\t") for i in data[1:]]
            data_full.extend(data)
    else:
        print(f"incorrect id_string {IDs} / {mode}")
        return []

    return data_full


def alpha_fold_cif_download(IDs_clusters, FOLDER_PATH):
    def quality_check(path):
        file_h = open(path)
        file = file_h.readlines()
        file_h.close()

        val = [i for i in file if "_ma_qa_metric_global.metric_value" in i]
        if not val or float(val[0].split(" ")[1]) < 50:
            print("LOW QUALITY - DELETING: ", path)
            return False
        else:
            return True

    if not os.path.exists(FOLDER_PATH):
        os.makedirs(FOLDER_PATH)

    for ids in IDs_clusters:
        for id in ids:
            filepath = f"{FOLDER_PATH}/{id}.cif"
            url = f"https://alphafold.ebi.ac.uk/files/AF-{id}-F1-model_v3.cif"
            request.urlretrieve(url, filepath)
            if quality_check(filepath):
                break
            else:
                os.remove(filepath)


def download_job_wrapper(batch):
    path, IDs = batch

    if "_" in IDs[0][0]:
        mode = "entry"
    else:
        mode = "accession"

    Alpha_IDs = get_data_from_Uniprot(IDs, mode)
    Alpha_IDs = [i for i in Alpha_IDs if len(i) == 4]

    clusters = []
    for cluster in IDs:
        temp = []
        for id in cluster:
            for ref in Alpha_IDs:
                if mode == "entry":
                    if ref[1] == id:
                        temp.append(ref[3].rstrip(";"))
                elif mode == "accession":
                    if ref[0] == id:
                        temp.append(ref[3].rstrip(";"))
        if temp:
            clusters.append(temp)
            temp = []


    alpha_fold_cif_download(clusters, path)


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


def get_sequences(path):
    IDs = os.listdir(path)
    IDs = [i.rstrip(".cif") for i in IDs]
    id_strings = [f"%28accession%3A{i}%29" for i in IDs]
    data_full = []

    for i in range(0, len(id_strings), 50):
        url = "https://rest.uniprot.org/uniprotkb/stream?fields=accession%2Csequence%2Cxref_pdb&format=tsv&query="
        url += "%20OR%20".join(id_strings[i:i+50])

        data = request.urlopen(url).readlines()
        data = [i.decode("utf-8") for i in data]
        data = [i.rstrip().split("\t") for i in data[1:]]

        data_full.extend(data)

    return data_full


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
            base_path1 = f"./DATA/{family}/family_seq"
            base_path2 = f"./DATA/{family}/AlphaFold"
        elif mode == "clans":
            clan = ""
            for i in range(len(clans_map)):
                if clans_map[i][0] == family:
                    if len(clans_map[i]) == 2:
                      clan = clans_map[i][1]
            if not clan:
                clan = "No_Clan"

            base_path1 = f"./DATA/{clan}/{family}/family_seq"
            base_path2 = f"./DATA/{clan}/{family}/AlphaFold"

        download_Pfam(family, f"{base_path1}")
        CD_HIT_wrapper(f"{base_path1}/{family}_full.fasta")
        IDs = get_IDs_from_cluster(f"{base_path1}/{family}_full_cluster.fasta.clstr")


        batches = [[base_path2, IDs[i:i+batch_size]] for i in range(0, len(IDs), batch_size)]
        major_batches = [batches[i:i+major_batch_size] for i in range(0, len(batches), major_batch_size)]
        workers = multiprocessing.Pool(workers_n)

        for major_batch in tqdm(major_batches, leave=False, desc="DOWNLOAD"):
            workers.map(download_job_wrapper, major_batch)
        workers.close()
        workers.join()

        workers = multiprocessing.pool.ThreadPool(file_workers_n)
        files = os.listdir(f"{base_path2}")
        files = [f"{base_path2}/{i}" for i in files]
        file_batches = [files[i: i+major_batch_size_files] for i in range(0, len(files), major_batch_size_files)]

        results = []

        for file_batch in tqdm(file_batches, leave=False, desc="KNOT CAL"):
            t = workers.map(simple_knot_calculation, file_batch)
            results.extend(t)

        workers.close()
        workers.join()

        seq_data = get_sequences(base_path2)
        PDB = [i for i in seq_data if len(i) == 3]
        for i in PDB:
            print(i[0], i[2])
        print(len(PDB))

        outfile = open(f"{base_path1}/{family}_summary.tsv", "w")
        for i in results:
            for i1 in seq_data:
                if i[0] == i1[0]:
                    outfile.write(f"{i1[0]}\t{i1[1]}\t{''.join(i[1].split('_'))}\n")
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

