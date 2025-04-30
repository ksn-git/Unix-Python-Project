#!/bin/env python3

#check deviation
def check_deviation(window_elem,motif_elem,penalty_elem,deviation,max_deviation):
    """Compare the element from the window with the motif element to return the updated deviation
    and boolean indicating success or not (exceeding max)"""

    #if multiple charactors are found
    if isinstance(motif_elem,set):
        if window_elem not in motif_elem:
            deviation += penalty_elem
            if deviation > max_deviation:
                return deviation, True,'X'      #signal stop now
            return deviation, False, window_elem
    #if a single charactor is found
    else: 
        if window_elem != motif_elem:
            deviation += penalty_elem
            if deviation > max_deviation:
                return deviation, True, 'X'
            return deviation, False, window_elem
    #return False if match continues
    return deviation, False, window_elem         #continue match

def find_motif(sequence, motif_list, penalty_list, max_deviation, minimum_gap, maximum_gap):
    """Generator that searches for motif.
    Yields (position, deviation, matched_sequence) when a match is found."""

    try:
        star_index = motif_list.index("*")
        len_part_2 = len(motif_list) - star_index - 1
    except ValueError:
        star_index = None

    #define motif length (gap or no gap)
    motif_core_length = len(motif_list) - (1 if star_index is not None else 0)

    for i in range(0, len(sequence) - motif_core_length - maximum_gap + 1):
        deviation = 0
        success = False

        if star_index is None:
            # No gaps in motif
            window = sequence[i:i + len(motif_list)]
            match = []

            for j in range(len(motif_list)):
                w = window[j]
                m = motif_list[j]
                p = penalty_list[j]
                deviation, exceeded, char = check_deviation(w, m, p, deviation, max_deviation)
                match.append(char)
                if exceeded:
                    break
            else:
                success = True
                matched_seq = ''.join(match)

        else:
            # Gapped motif
            for gap_size in range(minimum_gap, maximum_gap + 1):
                window = sequence[i:i + len(motif_list) - 1 + gap_size]
                match = []
                deviation = 0
                exceeded = False

                # Before the gap
                for j in range(star_index):
                    w = window[j]
                    m = motif_list[j]
                    p = penalty_list[j]
                    deviation, exceeded, char = check_deviation(w, m, p, deviation, max_deviation)
                    match.append(char)
                    if exceeded:
                        break

                if exceeded:
                    continue

                match.append('*')  # gap symbol

                # After the gap
                for k in range(len_part_2):
                    w_idx = star_index + gap_size + k
                    if w_idx >= len(window):
                        exceeded = True
                        break
                    w = window[w_idx]
                    m = motif_list[star_index + 1 + k]
                    p = penalty_list[star_index + 1 + k]
                    deviation, exceeded, char = check_deviation(w, m, p, deviation, max_deviation)
                    match.append(char)
                    if exceeded:
                        break

                if not exceeded:
                    success = True
                    matched_seq = ''.join(match)
                    break

        if success and deviation <= max_deviation:
            print(f"Yielding match at {i}: {matched_seq}, deviation: {deviation}")
            yield (i, deviation, matched_seq)
        else:
            print(f"No match at {i}, deviation too high: {deviation}")

