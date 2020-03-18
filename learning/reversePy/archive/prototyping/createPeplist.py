f = open("pepList.txt", "r")
out = open("pepList.peptide.csv", "w+")

valid = 'ARNDCEQGHILKMFPSTWYV'

proteinList = []
peptideList = []



for line in f:

    if len(line.split()) != 2:
        continue
    line = line.split()
    try:
        #protein = line[1].split("|")[1]
        protein = line[1]
        peptide = line[0]
    except:
        print(line)

    found = False
    for letter in peptide:
        if letter not in valid:
            print(peptide)
            found = True
            continue
    if found:
        continue

    """

    if peptide in peptideList:
        continue
    elif '(' in peptide or 'Z' in peptide or 'B' in peptide or 'U' in peptide:
        continue
    else:
        proteinList.append(protein)
        peptideList.append(peptide)
    """

    proteinList.append(protein)
    peptideList.append(peptide)

out.write('"protein","sequence"\n')
for ele in zip(proteinList, peptideList):
    out.write('"' + ele[0] + '","' + ele[1] + '"\n')
