file_h = open("_all_chains_knotted.txt")
file = file_h.readlines()
file_h.close()

data_families = [i.split(";")[2] for i in file]

unique_families = []
for i in data_families:
    if i not in unique_families:
        unique_families.append(i)

file_h = open("family_clan.tsv")
file = file_h.readlines()
file_h.close()

data_clans = [i.rstrip().split("\t") for i in file]

clans = []
all_families = []
for i in unique_families:
    for i1 in data_clans:
        if i == i1[0]:
            if len(i1) == 2 and i1[1] not in clans:
                clans.append(i1[1])
            else:
                all_families.append(i)

for i in clans:
    for i1 in data_clans:
        if len(i1) == 2 and i1[1] == i:
            all_families.append(i1[0])

newfile = open("test_full.txt", "w")
newfile.write("\n".join(all_families))
newfile.close()


