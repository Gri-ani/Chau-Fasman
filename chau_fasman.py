alpha_props = {
    'E': 1.53, 'A': 1.45, 'L': 1.34, 'H': 1.24, 'M': 1.20,
    'Q': 1.17, 'W': 1.14, 'V': 1.14, 'F': 1.12, 'K': 1.07,
    'I': 1.00, 'D': 0.98, 'T': 0.82, 'S': 0.79, 'R': 0.79,
    'C': 0.77, 'N': 0.73, 'Y': 0.61, 'P': 0.59, 'G': 0.53
}

beta_props = {
    'M': 1.67,
    'V': 1.65,
    'I': 1.60,
    'C': 1.30,
    'Y': 1.29,
    'F': 1.28,
    'Q': 1.23,
    'L': 1.22,
    'T': 1.20,
    'W': 1.19,
    'A': 0.97,
    'R': 0.90,
    'G': 0.81,
    'D': 0.80,
    'K': 0.74,
    'S': 0.72,
    'H': 0.71,
    'N': 0.65,
    'P': 0.62,
    'E': 0.26
}


def find_alpha_helices(sequence):
    """Find all likely alpha helices in the sequence. Returns a list
    of [start, end] pairs for the alpha helices."""
    start = 0
    helices = []
    # Try each window
    while (start + 6 < len(sequence)):
        # Count the number of "good" amino acids (those likely to be
        # in an alpha helix).
        num_good = 0
        found_four_good = True
        # Taking a set of 6 amino acids as given in the example
        for i in range(start, start+6):
            if i < len(sequence) and alpha_props[sequence[i]] >= 1:
                num_good += 1
                if num_good == 4:
                    found_four_good = False
        # If we have at least 4 residues with P[H] >= 1
        if not found_four_good:
            # Extend the residues if the residue is on the right corner
            # Extension only in the left and vice versa
            [extended_start, extended_end] = extend(sequence, start, start+6)
            if [extended_start, extended_end] not in helices:
                helices.append([extended_start, extended_end])
        # Go on to the next frame
        start += 1

    return helices



def CF_find_beta(seq):
    """Find all likely beta strands in sequence. Returns a list
    of [start, end] pairs for the beta strands."""
    start = 0
    beta_strands = []
    # Try each window
    while (start + 5 < len(seq)):
        # Count the number of "good" amino acids (those likely to be
        # in a beta strand).
        num_good = 0
        found_three_good = False
        # Taking a set of 5 amino acids as given in the example
        for i in range(start, start+5):
            if i < len(seq) and beta_props[seq[i]] >= 1:
                num_good += 1
                if num_good == 3:
                    found_three_good = True
        # If we have at least 3 residues with P(S) >= 1
        if found_three_good:
            # Extend the residues if the residue is on the right corner
            [extended_start, extended_end] = extend(seq, start, start+5)
            if [extended_start, extended_end] not in beta_strands:
                beta_strands.append([extended_start, extended_end])
        # Go on to the next frame
        start += 1

    return beta_strands


def extend(seq, start, end):

    # We extend the region in both directions until the sum of propensity < 4

    # Extend the sequence to the right
    while end <= len(seq)-1:
        sum_right = 0
        for x in seq[end-3:end+1]:
            sum_right += alpha_props[x]
        if sum_right >= 4:
            end += 1
        else:
            break

    # Extend the sequence to the left
    while start > 0:
        sum_left = 0
        for x in seq[start-1:start+3]:
            sum_left += alpha_props[x]
        if sum_left >= 4:
            start -= 1
        else:
            break

    return [start, end]







protein_sequence = "MNASSEGESFAGSVQIPGGTTVLVELTPDIHICGICKQQFNNLDAFVAHKQSGCQLTGTSAAAP\
STVQFVSEETVPATQTQTTTRTITSETQTITVSAPEFVFEHGYQTYLPTESNENQTATVISLPA\
KSRTKKPTTPPAQKRLNCCYPGCQFKTAYGMKDMERHLKIHTGDKPHKCEVCGKCFSRKDKLKT\
HMRCHTGVKPYKCKTCDYAAADSSSLNKHLRIHSDERPFKCQICPYASRNSSQLTVHLRSHTAS\
ELDDDVPKANCLSTESTDTPKAPVITLPSEAREQMATLGERTFNCCYPGCHFKTVHGMKDLDRH\
LRIHTGDKPHKCEFCDKCFSRKDNLTMHMRCHTSVKPHKCHLCDYAAVDSSSLKKHLRIHSDER\
PYKCQLCPYASRNSSQLTVHLRSHTGDTPFQCWLCSAKFKISSDLKRHMIVHSGEKPFKCEFCD\
VRCTMKANLKSHIRIKHTFKCLHCAFQGRDRADLLEHSRLHQADHPEKCPECSYSCSSAAALRV\
HSRVHCKDRPFKCDFCSFDTKRPSSLAKHVDKVHRDEAKTENRAPLGKEGLREGSSQHVAKIVT\
QRAFRCETCGASFVRDDSLRCHKKQHSDQSENKNSDLVTFPPESGASGQLSTLVSVGQLEAPLE\
PSQDL"


def print_protein_sequence(protein_sequence):
    print("Protein Sequence:")
    print(protein_sequence)

# Example usage:

print()
print_protein_sequence(protein_sequence)

print()




def print_helix_regions(seq, helix_regions):
    result = ''
    helix_sequence = ''
    for i in range(len(seq)):
        is_helix = False
        for start, end in helix_regions:
            if i >= start and i < end:
                is_helix = True
                break
        if is_helix:
            result += 'H'
            helix_sequence += 'H'
        else:
            result += '-'
            helix_sequence += '_'
    print(result)
    return helix_sequence


def print_beta_regions(seq, beta_strands):
    result = ''
    sheet_sequence = ''
    for i in range(len(seq)):
        is_strand = False
        for start, end in beta_strands:
            if i >= start and i < end:
                is_strand = True
                break
        if is_strand:
            result += 'S'
            sheet_sequence += 'S'
        else:
            result += '-'
            sheet_sequence += '_'
    print(result)
    return sheet_sequence

#
# def compile_sec_struct():
#     secondary_structure = ['-'] * len(helix_sequence)
#     for aa in range(len(secondary_structure)):
#         if helix_sequence[aa] == 'H' and sheet_sequence[aa] == '-':  # only H
#             secondary_structure[aa] = 'H'
#         elif sheet_sequence[aa] == 'S':  # only S
#             secondary_structure[aa] = 'S'
#         elif helix_sequence[aa] == 'H' and sheet_sequence[aa] == 'S':  # H or S
#             if helix_value_list[aa] > sheet_value_list[aa]:
#                 secondary_structure[aa] = 'H'
#             else:
#                 secondary_structure[aa] = 'S'
#         elif helix_sequence[aa] == '-' and sheet_sequence[aa] == '-':  # none -> C
#             secondary_structure[aa] = 'C'
#         elif helix_sequence[aa] == 'H' and sheet_sequence[aa] == 'S':  # H or S
#             if helix_value_list[aa] > sheet_value_list[aa]:
#                 secondary_structure[aa] = 'H'
#             else:
#                 secondary_structure[aa] = 'S'
#     sec_struct_as_string = ''.join(secondary_structure)
#     return sec_struct_as_string
#

# Example usage
# protein_sequence = "MNASSEGESFAGSVQIPGGTTVLVELTPDIHICGICKQQFNNLDAFVAHKQSGCQLTGTSAAAPSTVQFVSEETVPATQTQTTTRTITSETQTITVSAPEFVFEHGYQTYLPTESNENQTATVISLPAKSRTKKPTTPPAQKRLNCCYPGCQFKTAYGMKDMERHLKIHTGDKPHKCEVCGKCFSRKDKLKT...
helix_regions = find_alpha_helices(protein_sequence)
beta_strands = CF_find_beta(protein_sequence)

print("------------------------alpha sequence--------------------------------")
print()
helix_sequence = print_helix_regions(protein_sequence, helix_regions)


print()
# print(find_alpha_helices(protein_sequence))
print()
print()


print("-----------------------------beta sequence---------------------------------")
print()
sheet_sequence = print_beta_regions(protein_sequence, beta_strands)
print()
print()
# print(CF_find_beta(protein_sequence))

#find the code for conflicting regions








c_sequence = ''
for helix_residue, sheet_residue in zip(helix_sequence, sheet_sequence):
    if helix_residue == '-' or sheet_residue == '-':
        c_sequence += '-'
    elif helix_residue == 'H' and sheet_residue == 'S':
        c_sequence += 'C'
    else:
        c_sequence += '-'

print("--------------------------------------conflicting sequence------------------------------------------------:")
print()
print(c_sequence)

def secondary_structure(helix_sequence, sheet_sequence, seq):
    result = ['' for _ in range(len(seq))]
    sequence_length = len(seq)
    i = 0
    while i < sequence_length:
        if (helix_sequence[i] == 'H' and sheet_sequence[i] == '_') or (helix_sequence[i] == '_' and sheet_sequence[i] == 'H'):
            result[i] = 'H'
            i += 1
        elif helix_sequence[i] == '_' and sheet_sequence[i] == '_':
            result[i] = '-'
            i += 1
        elif (helix_sequence[i] == '_' and sheet_sequence[i] == 'S') or (helix_sequence[i] == 'S' and sheet_sequence[i] == '_'):
            result[i] = 'S'
            i += 1
        else:
            n = 0
            while i < sequence_length and ((helix_sequence[i] == 'H' and sheet_sequence[i] == 'S') or (helix_sequence[i] == 'S' and sheet_sequence[i] == 'H')):
                n += 1
                i += 1
            p1 = sum(alpha_props[residue] for residue in seq[i-n:i])
            p2 = sum(beta_props[residue] for residue in seq[i-n:i])
            if p1 > p2:
                for k in range(i-n, i):
                    result[k] = 'H'
            else:
                for k in range(i-n, i):
                    result[k] = 'S'
    return ''.join(result)


print()
print()
print("---------------------------------------------------Secondary Structure------------------------------------------------------")

print()
print(secondary_structure(helix_sequence, sheet_sequence, protein_sequence))