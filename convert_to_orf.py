file_path = "data\dnasequence_test.txt"

with open(file_path, "r") as f:
    data = f.read().strip("\n")

class SequenceHandler():
    def __init__(self, seq) -> None:
        self.seq = seq
        self.check_set = {'T', 'C', 'A', 'G'}

    def validate_dna(self):
        """
        Return errors if unexepcted nucleotides are encountered in the given sequence.

        :param seq_set: nucleotide sequence
        :param check_set: set to check against for membership testing
        :type seq_set: array
        :type check_set: set
        :return: error message
        :rtype: str
        """
        error_messages = [
        "DNA Validation: Ok",
        "Error: DNA Sequence failed to Validate!",
        "Error: RNA Sequence found!"
        ]
        error_key = None
        set_seq = set(self.seq)

        if set_seq.issubset(self.check_set) == True:
            error_key = 0
        elif "U" in set_seq: error_key = 2
        else: error_key = 1

        return error_messages[error_key]

    def get_start_codon_indices(self, frame_number):
        seq_frame = self.seq[frame_number:]
        start_codon_indices = [index for index in range(frame_number, len(seq_frame)) if seq_frame[index: index + 3] == "ATG"]
        print(frame_number, self.seq[frame_number:])
        return start_codon_indices

    def get_end_codon_indicies(self, frame_number):
        seq_frame = self.seq[frame_number:]
        stop_codon_indicies = [index for index in range(frame_number, len(seq_frame)) if seq_frame[index: index + 3] in ["TGA", "TAA", "TAG"]]
        return stop_codon_indicies
    
    def get_orfs(self, frame_number = 0):
        seq_frame = self.seq[frame_number:]
        start_codon_indices = self.get_start_codon_indices(frame_number) 
        end_codon_indicies = self.get_end_codon_indicies(frame_number)
        print(start_codon_indices, end_codon_indicies)

        orfs = [] 
        for start_index in start_codon_indices:
            end_index = None
            for stop_index in end_codon_indicies:
                if stop_index > start_index and (stop_index - start_index) % 3 == 0:
                    end_index = stop_index
                    break

            if end_index is not None:
                orf = seq_frame[start_index:end_index]
                orfs.append((start_index, orf))
        
        return orfs
    
handler = SequenceHandler(data)

if __name__ == "__main__": 
    orf = []
    for x in range(3):
        orf.append(handler.get_orfs(x))
    print(orf)
