dnaseq = "TCTGACTGAGCTACCCTCGGCTCGATCGCTCTCGATCGTCTCGGCTCTATATAATATGCTGATTATATATAGCCTCTCGCCTAGCTGCTATTATAGAATGCTGCTCGCTCGATGCTGCCGCTATGAGATATGATGCTGCTCGTAGTAGCCTGCTCGATCGTGCTAGCGACTGACTGTCGTAGCTGCTAGCTGCTAGTCGACGATCGACTGACTTCTCTCTCGCGATGAAGTAGCTGCTAGCTAGCTAGCTACGTAGCTATCGATCGATCGACTGATGCATGCAGTCGTAGCATGCTTAATTAGAGAGAG"
frame = 0

def orf(seq):
    orf = ''
    i = 0
    while i <= len(seq):
        for j in range (i,len(seq) - 2, 3):
            start = seq[i:i+3]
            if start == "ATG":
                starthere = i

                for k in range(starthere, len(seq) - 2, 3):
                    stophere = seq[k:k+3]
                    if stophere in ("TAG", "TGA", "TAA"):
                        break
                    orf += stophere
                print(orf)
                orf = ''
                i += 3
                break
        else:
            i+=3

orf(dnaseq[frame:])