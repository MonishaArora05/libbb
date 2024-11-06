import re

# Reading and formatting the sequence
sequence = ""
with open("sequence.fasta") as f:
    sequence = ''.join(f.read().splitlines()[1:]).upper()

print(sequence)

# Calculating GC content
g_c = sequence.count('G')
c_c = sequence.count('C')
total_c = len(sequence)
gc_per = ((g_c + c_c) / total_c) * 100

print("G count: ", g_c)
print("C count: ", c_c)
print("GC content percentage: ", gc_per)

# Calculating AT content
a_c = sequence.count('A')
t_c = sequence.count('T')
at_per = ((a_c + t_c) / total_c) * 100

print("A count: ", a_c)
print("T count: ", t_c)
print("AT content percentage: ", at_per)

# Calculating AT/GC Ratio
ratio = (a_c + t_c) / (g_c + c_c)
print(f"AT/GC Ratio: {ratio:.2f}")

# Finding coding regions
start_codon = 'ATG'
stop_codons = ['TAA', 'TAG', 'TGA']
coding_regions = []  # Initialize an empty list to store coding regions

# Find the start codon and search for coding regions
start_index = sequence.find(start_codon)

while start_index != -1:
    for stop_codon in stop_codons:
        stop_index = sequence.find(stop_codon, start_index + 3)
        if stop_index != -1 and (stop_index - start_index) % 3 == 0:
            coding_region = sequence[start_index:stop_index + 3]
            coding_regions.append(coding_region)
            break
    start_index = sequence.find(start_codon, start_index + 1)

# Output results for coding regions
if coding_regions:
    print("Coding Regions Found")
    for i, coding_region in enumerate(coding_regions, start=1):
        print(f"\nRegion {i}: {coding_region}\nLength: {len(coding_region)}")
else:
    print("No coding regions found")

def find_motifs(sequence, motif):
    print(f"\nSearching for motif: {motif}")
    matches = [match.start() for match in re.finditer(motif, sequence)]
    if matches:
        print(f"Motif '{motif}' found at positions: {matches}")
    else:
        print(f"Motif '{motif}' not found")

# Define motif
motif = "TATAA"
find_motifs(sequence, motif)