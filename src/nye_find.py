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
                return deviation, True      #signal stop now
    #if a single charactor is found
    else: 
        if window_elem != motif_elem:
            deviation += penalty_elem
            if deviation > max_deviation:
                return deviation, True
    #return False if match continues
    return deviation, False         #continue match

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

            j = 0
            while j < len(motif_list):
                w = window[j]
                m = motif_list[j]
                p = penalty_list[j]
                deviation, exceeded = check_deviation(w, m, p, deviation, max_deviation)
                if exceeded:
                    break
                j += 1
            else:
                # Only if loop finishes without break
                success = True
                matched_seq = window

        else:
            # Gapped motif
            window = sequence[i:i + len(motif_list) - 1 + maximum_gap]

            j = 0
            while j < len(motif_list):
                if j == star_index:
                    # Handle gap
                    for gap_size in range(minimum_gap, maximum_gap + 1):
                        local_deviation = deviation
                        exceeded = False
                        m = 0
                        while m < len_part_2:
                            w_idx = j + gap_size + m
                            motif_idx = j + m + 1
                            if w_idx >= len(window):
                                exceeded = True
                                break
                            w = window[w_idx]
                            m_elem = motif_list[motif_idx]
                            p = penalty_list[motif_idx]
                            local_deviation, exceeded = check_deviation(w, m_elem, p, local_deviation, max_deviation)
                            if exceeded:
                                break
                            m += 1
                        if not exceeded:
                            deviation = local_deviation
                            matched_seq = window[:star_index] + window[star_index + gap_size: star_index + gap_size + len_part_2]
                            success = True
                            break
                    break  # after handling the gap
                else:
                    w = window[j]
                    m_elem = motif_list[j]
                    p = penalty_list[j]
                    deviation, exceeded = check_deviation(w, m_elem, p, deviation, max_deviation)
                    if exceeded:
                        break
                j += 1

        if success and deviation <= max_deviation:
            print(f"Yielding match at {i}: {matched_seq}, deviation: {deviation}")
            yield (i, deviation, matched_seq)
        else:
            print(f"No match at {i}, deviation too high: {deviation}")

