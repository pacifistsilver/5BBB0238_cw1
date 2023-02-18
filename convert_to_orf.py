filename = "dnasequence.txt"
file_path = "C:\\Users\\danie\\Downloads\\dnasequence.txt"
dna_set = {'T', 'C', 'A', 'G'}

with open(file_path, "r") as f:
    data = f.read().strip("\n")

class ProcessORF():
    def __init__(self) -> None:
        pass
def check_error_cases(seq, check_set):
    """
    Return errors if unexepcted nucleotides are encountered in the given sequence.

    :param seq_set: nucleotide sequence
    :param check_set: set to check against for membership testing
    :type seq_set: array
    :type check_set: set
    :return: error message
    :rtype: str
    """
    error_message_dict = {0: "DNA Sequence TypeCheck: ok", 1: "Error: Unexpected Type encountered in dna sequence", 2: "Error: uracil found in dna sequence"}
    error_key = None
    seq = set(seq)

    if seq.issubset(check_set) == True:
        error_key = 0
    elif "U" in seq: error_key = 2
    else: error_key = 1

    return error_message_dict[error_key]

def get_start_codon_indices(seq):
    start_codon_indices = [index for index in range(len(seq)) if seq[index: index + 3] == "ATG"]
    return start_codon_indices

def get_terminal_codons(seq):
    tga_codon_indices = [index for index in range(len(seq)) if seq[index: index + 3] == "TGA"]
    taa_stop_codon_indices = [index for index in range(len(seq)) if seq[index: index + 3] == "TAA" or "TAG"]
    tag_stop_codon_indicies = [index for index in range(len(seq)) if seq[index: index + 3] == "TAG"]
    return (tga_codon_indices + taa_stop_codon_indices + tag_stop_codon_indicies)

def get_reading_frame(seq, frame_number = 0):
    return seq[frame_number:]

if __name__ == "__main__": 
    dnaseq_code = check_error_cases(data, dna_set)
    #if 'ATG' in seq:  print("ATG Codons found in sequence")
    for x in range(3):
        current_frame = get_reading_frame(x)
        print(current_frame)
        print(get_start_codon_indices(current_frame), get_terminal_codons(current_frame))
        
