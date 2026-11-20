# Manhattan plot for late-onset Alzheimer's disease

# The code reads GWAS summary statistics from a text file and plots a Manhattan plot.
# The file must be downloaded first and saved in the SAME folder as this script (about 543 MB):
# https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST007001-GCST008000/GCST007511/Kunkle_etal_Stage1_results.txt

# Import libraries
# Import the plotting library, used to draw the graph
import matplotlib.pyplot as plt
# Import the maths library, used to calculate -log10 of the p-values
import math
# Import the json library, used to read the answer from the Ensembl website
import json
# Import the time library, used to pause between lookup tries, to avoid overloading
import time
# Import the library used to connect to the Ensembl website
import urllib.request

# The path of the file that was downloaded
file_name = "Kunkle_etal_Stage1_results.txt"


# 1. Python - read file name and check data  
# The code reads the file line by line
# Only SNPs (single nucleotide polymorphisms) with a p-value below 0.01 are kept for plotting, because the file is very large
# The data is stored using a dictionary: for each chromosome, a list of positions, -log10 p-values and SNP names
snp_data = {}

# These counters are used to check the data
total_rows = 0
skipped_rows = 0
smallest_p = 1.0
smallest_p_snp = ""

with open(file_name) as f:
    # The first line is the header, so it is read and skipped
    header = f.readline().split()
    print("Columns:", header)

    for line in f:
        # The function .split() splits the line into a list of items using the spaces
        items = line.split()
        total_rows += 1

        # The code extracts the chromosome, position, SNP name and p-value
        # If a value cannot be converted, the row is skipped and counted
        try:
            chromosome = int(items[0])
            position = int(items[1])
            snp_name = items[2]
            p_value = float(items[7])
        except (ValueError, IndexError):
            skipped_rows += 1
            continue

        # The code ignores p-values of 0 or below because log10 cannot be calculated for them
        if p_value <= 0:
            skipped_rows += 1
            continue

        # The code keeps track of the smallest p-value (the strongest association)
        if p_value < smallest_p:
            smallest_p = p_value
            smallest_p_snp = snp_name

        # Only SNPs with p below 0.01 are saved
        if p_value < 0.01:
            if chromosome not in snp_data:
                snp_data[chromosome] = {"Positions": [], "Log_p": [], "Names": []}
            snp_data[chromosome]["Positions"].append(position)
            snp_data[chromosome]["Log_p"].append(-math.log10(p_value))
            snp_data[chromosome]["Names"].append(snp_name)

# The code prints a summary so the data can be checked
print("Total rows read:", total_rows)
print("Rows skipped:", skipped_rows)
print("Chromosomes found:", sorted(snp_data.keys()))
print("Smallest p-value:", smallest_p, "(", smallest_p_snp, ")")



# 2. Creating the Manhattan plot
# The significance used in GWAS is p = 5 x 10^-8
significance_line = -math.log10(5e-8)

# This dictionary stores the x position of each SNP on the graph
# The code puts the chromosomes side by side on one x-axis
# The offset moves each chromosome to the right by the length of the previous ones
snp_x = {}
offset = 0
tick_positions = []
tick_labels = []
colours = ["blue", "grey"]

plt.figure(figsize=(13, 5))

for chromosome in sorted(snp_data.keys()):
    data = snp_data[chromosome]

    # The positions are shifted by the offset so the chromosomes do not overlap
    x_values = [pos + offset for pos in data["Positions"]]
    y_values = data["Log_p"]

    # The x position of every SNP is saved so the top SNPs can be circled later
    for i in range(len(data["Positions"])):
        snp_x[(chromosome, data["Positions"][i])] = x_values[i]

    # The colours alternate between chromosomes
    plt.scatter(x_values, y_values, s=4, color=colours[chromosome % 2])

    # The label for this chromosome is placed in the middle of it
    tick_positions.append(offset + max(data["Positions"]) / 2)
    tick_labels.append(chromosome)

    # The offset is updated for the next chromosome
    offset = offset + max(data["Positions"])

# The significance line is added
plt.axhline(significance_line, color="red", linestyle="--", label=r"$p = 5 \times 10^{-8}$")

# A title, axis labels and chromosome labels for the plot
plt.title("Manhattan plot: late-onset Alzheimer's disease")
plt.xlabel("Chromosome")
plt.ylabel(r"$-\log_{10}(p)$")
plt.xticks(tick_positions, tick_labels, fontsize=8)
plt.legend()



# 3. Selecting the top 10 SNPs 
# The code collects every SNP above the significance line into one list
significant = []
for chromosome in snp_data:
    data = snp_data[chromosome]
    for i in range(len(data["Names"])):
        if data["Log_p"][i] > significance_line:
            significant.append((data["Log_p"][i], chromosome, data["Positions"][i], data["Names"][i]))

# The list is sorted so the strongest associations are first
significant.sort(reverse=True)
print("SNPs above the significance line:", len(significant))

# The top 10 are printed
for item in significant[:10]:
    print("Chromosome", item[1], "| position", item[2], "|", item[3], "| -log10(p) =", round(item[0], 1))


# 4. Finding the strongest SNP on each chromosome
# The top 10 SNPs are all in the same peak (chromosome 19), so one SNP per chromosome is used instead
# The list "significant" is already sorted from strongest to weakest, so the first SNP seen for each chromosome is the strongest
best_per_chromosome = []
chromosomes_seen = []
for item in significant:
    if item[1] not in chromosomes_seen:
        chromosomes_seen.append(item[1])
        best_per_chromosome.append(item)

# The list is sorted by chromosome so the numbers on the graph go from left to right
best_per_chromosome.sort(key=lambda item: item[1])


# 5. Looking up the gene for each SNP using the Ensembl website (needs internet)
# The function sends the rs number to Ensembl and reads the gene name from the answer
def find_gene(rs_number):
    # Only names starting with "rs" can be looked up (some names in the file are insertions or deletions)
    if not rs_number.startswith("rs"):
        return "no rs number"
    url = "https://rest.ensembl.org/vep/human/id/" + rs_number + "?content-type=application/json"
    # The code tries up to 3 times, with a pause between tries, because the website sometimes refuses or is slow
    answer = None
    error_text = ""
    for attempt in range(3):
        try:
            with urllib.request.urlopen(url, timeout=30) as response:
                answer = json.loads(response.read())
            break
        except Exception as error:
            error_text = str(error)[:40]
            time.sleep(2)
    if answer is None:
        return "lookup failed (" + error_text + ")"

    # The code looks through the results for the closest gene
    # A SNP inside a gene has no distance, so it counts as 0
    best_gene = "no gene found"
    best_distance = 10**9
    for result in answer:
        for consequence in result.get("transcript_consequences", []):
            if "gene_symbol" in consequence:
                distance = consequence.get("distance", 0)
                if distance < best_distance:
                    best_distance = distance
                    best_gene = consequence["gene_symbol"]
    return best_gene

# The gene is looked up for each SNP, with a short pause so the website is not overloaded
genes = []
for item in best_per_chromosome:
    genes.append(find_gene(item[3]))
    time.sleep(0.5)


# 6. Circling the SNPs on the graph with numbers
# The code draws a red circle around each SNP and writes its number next to it
for i in range(len(best_per_chromosome)):
    item = best_per_chromosome[i]
    x = snp_x[(item[1], item[2])]
    y = item[0]
    plt.scatter([x], [y], s=150, facecolors="none", edgecolors="red", linewidths=1.5)
    # The highest SNP is labelled on the right-hand side of its circle because there is no room above it
    if y > 100:
        plt.annotate(str(i + 1), (x, y), xytext=(16, -4), textcoords="offset points", fontsize=9, color="red", fontweight="bold", ha="left")
    else:
        # Every second number is placed higher so numbers on neighbouring chromosomes do not overlap
        height = 14 + 12 * (i % 2)
        plt.annotate(str(i + 1), (x, y), xytext=(0, height), textcoords="offset points", fontsize=9, color="red", fontweight="bold", ha="center")

# The key is printed so each number can be matched to a SNP and gene
print("Key for the numbered circles:")
for i in range(len(best_per_chromosome)):
    item = best_per_chromosome[i]
    # The p-value is calculated back from -log10(p)
    p_value = 10 ** (-item[0])
    print(i + 1, "| chromosome", item[1], "|", item[3], "| gene:", genes[i], "| p =", "%.1e" % p_value)

# The plot is saved as a picture and displayed
plt.savefig("manhattan_alzheimers.png", dpi=150)
plt.show()

# 7. Interpretation
    # 9 chromosomes have a SNP with p < 5 x 10^-8, therefore, these regions are linked to late-onset Alzheimer's disease
    # The strongest signal is on chromosome 19 (p = 2.7e-238), the region of the APOE gene, the best known genetic risk factor for Alzheimer's disease
    # The other circled SNPs are weaker, which indicates they have smaller effects
    # Limitations:
        # 1. GWAS shows association, but not cause. The gene listed is only the closest one to the SNP, and the effect could be on a different gene.
        # 2. The SNPs on chromosome 2 and chromosome 10 (circles 2 and 6) were found, but no gene was found for them. This could be because the SNPs are in non-coding regions further from a gene than the distance Ensembl searches. Will investigate further.


# 8. Reference : Kunkle BW et al. (2019) Genetic meta-analysis of diagnosed Alzheimer's disease identifies new risk loci. Nature Genetics. Data from the NHGRI-EBI GWAS Catalog (GCST007511).
# 9. Acknowledgement of AI : Claude AI was used to review errors in the code. 
