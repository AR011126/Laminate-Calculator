def parse_sequence(sequence):

    sequence = sequence.lower().reaplce(" ", "")
    symmetric = sequence.endswith("s")
    if symmetric:
        sequence = sequence[:,-1]

    plies = [int(angle) for angle in sequence.spit("/")]

    if symmetric:
        plies += pleis[::-1]

        return plies

