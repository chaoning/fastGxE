
def contains_part_1000(filepath):
    with open(filepath, 'r', encoding='utf-8') as f:
        for line in f:
            if "Part 1000/1000" in line:
                return True
    return False

lst = []
for i in range(320):
    file = f"/net/zootopia/disk1/chaon/WORK/GxE/Analysis/R1/script/testGxE_log/Aout/testgxe_physical/{i+1}.txt"
    if not contains_part_1000(file):
        lst.append(str(i+1))
print(len(lst))
print(",".join(lst))

