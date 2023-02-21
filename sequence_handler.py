with open("data\dnasequence.txt", "r") as a, open("data\dnasequence2.txt", "r") as b:
    data = a.read().strip("\n")
    data2 = b.read().strip("\n")

class SequenceHandler():
    """
    Contains methods to validate DNA, and extract ORF data. 
    
    :param seq: DNA sequence
    :type seq: str
    :ivar check_set: set containing T, C, A, G
    :ivartype check_set:set
    """
    def __init__(self, seq) -> None:
        self.seq = seq
        self.check_set = {'T', 'C', 'A', 'G'}

    def validate_dna(self):
        """
        Return errors if unexpected nucleotides are encountered in the given sequence.

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

        if set_seq.issubset(self.check_set) == True: # checks if set_seq is equivalent to check_set
            error_key = 0
        elif "U" in set_seq: # uracil was found
            error_key = 2
        else: # generic error
            error_key = 1

        return [error_messages[error_key], error_key]

    def get_start_codon_indices(self, frame_number):
        """
        Return a list of start codon indices in the given reading frame.
        
        :param frame_number: The reading frame
        :type frame_number: int
        :return: A list of start codon indices.
        :rtype: list[int]
        """
        seq_frame = self.seq[frame_number:] 
        start_codon_indices = [index for index in range(0, len(seq_frame), 3) if seq_frame[index: index + 3] == "ATG"] 
        # using list comprehension to get all indexes related to ATG
        return start_codon_indices

    def get_end_codon_indices(self, frame_number):
        """
        Return a list of end codon indices in the given reading frame.
        
        :param frame_number: The reading frame
        :type frame_number: int
        :return: A list of end codon indices.
        :rtype: list[int]
        """
        seq_frame = self.seq[frame_number:]
        stop_codon_indicies = [index for index in range(0, len(seq_frame), 3) if seq_frame[index: index + 3] in ["TGA", "TAA", "TAG"]]
        return stop_codon_indicies
    # Could probably combine these. 
    def get_orfs(self, frame_number = 0):
        """
        Return a list of Open Reading Frames (ORFs) in the given reading frame.
        
        :param frame_number: The reading frame
        :type frame_number: int
        :return: A list of ORFs.
        :rtype: list[tuple(int, int, str, int)]
        """
        seq_frame = self.seq[frame_number:]
        start_codon_indices = self.get_start_codon_indices(frame_number) 
        end_codon_indicies = self.get_end_codon_indices(frame_number)
        orfs = [] 
        # we don't want to search the entire sequence again, so we just look at the positions we already have.
        for start_index in start_codon_indices:
            end_index = None
            for stop_index in end_codon_indicies:
                if stop_index > start_index and (stop_index - start_index) % 3 == 0:
                    # we don't want an index to stop before it starts. 
                    # orfs must be a multiple of 3
                    end_index = stop_index # valid stop was found so leave
                    break
            
            if end_index is not None:
                orf = seq_frame[start_index:end_index] # the orf between ATG and valid stop codon.
                orfs.append((frame_number + 1, (start_index + frame_number) + 1, orf, len(orf))) 
                # must add frame number to start_index as otherwise, the base position would be wrong.
        return orfs
    
    def print_orf_table(self, paired_list):
        """
        Prints given orf data in tabel format

        :param paired_list: A list of tuples containing orf data
        :type paired_list: list[tuples]
        :return: None
        :rtype: None
        """
        headers = ("Frame", "Base Number", "Sequence", "Sequence Length")
        col_widths = []
        for i in range(len(headers)):
            max_len = max(len(str(row[i])) for row in paired_list + [headers]) 
            # spacings dependant on greatest length in paired_list
            col_widths.append(max_len + 2) # 2 for extra spacing

        row_format = "".join(["{:<" + str(width) + "}" for width in col_widths])
        # {:< int} use f string to dynamically change width of each column
        print(row_format.format(*headers))

        for row in paired_list:
            print(row_format.format(*row))
            # * used to unpack paired_list tuple. 
    
    def return_single_orf(self, frame, longest_only = False):
        """
        Returns longest ORF given the reading frame

        :param frame: The frame (0, 1, or 2) to search for ORFs in.
        :type frame: int
        :param longest_only: If True, only the longest ORF will be printed. Defaults to False.
        :type longest_only: bool
        :return: A list containing the ORFs to be printed.
        :rtype: list[tuple]
        """
        orf_data = sorted(self.get_orfs(frame), key=lambda tup: tup[3], reverse=True)
        orfs_to_print = [orf_data[0] if orf_data and longest_only else orf_data]
        if orfs_to_print:
            return orfs_to_print

    def determine_mutation(self, alt_seq):
        # 1 - sub, 2 - del, 3 - insert 
        mutation_type = None

        return mutation_type

handler = SequenceHandler(data2)


if __name__ == "__main__": 
    result, error_key = handler.validate_dna()
    if error_key != 0: 
        print(result)
    else: 
        for frame in range(3):
            longest_orf = handler.return_single_orf(frame, True)
            longest_orf_tuples = [(*tup, ) for tup in longest_orf]   
            handler.print_orf_table(longest_orf_tuples), print()
           
        for frame in range(3):
            orf_data = handler.get_orfs(frame)
            if orf_data:
                handler.print_orf_table(orf_data)
                print()



 
            


