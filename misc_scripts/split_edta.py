import sys

input_file = sys.argv[1]
prefix = sys.argv[2]

features = {
    "copia": ["Copia_LTR_retrotransposon"],
    "gypsy": ["Gypsy_LTR_retrotransposon"],
    "helitron": ["helitron"],
    "tir": ["hAT_TIR_transposon", "Mutator_TIR_transposon", "PIF_Harbinger_TIR_transposon", "Tc1_Mariner_TIR_transposon", "CACTA_TIR_transposon"],
}

for feature in features.keys():
    with open(input_file) as infile, open(f"{prefix}{feature}.bed", "w") as outfile:
        for line in infile:
            if line.startswith("#"):
                continue
            line = line.split("\t")
            if line[2] in features[feature]:
                outfile.write(f"{line[0]}\t{line[3]}\t{line[4]}\n")