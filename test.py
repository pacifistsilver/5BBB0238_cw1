frame_number = 1
start_codon_indices = []
seq = "CTATGGGCACGATGGTCCCATGACC"
seq_frame = seq[frame_number:]

for index in range(frame_number, len(seq)):
    if seq_frame[index: index + 3] == "ATG":
        start_codon_indices.append(index)

print(start_codon_indices)