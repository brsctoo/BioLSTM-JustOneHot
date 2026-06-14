def clean_location_string(location):
    """
    Clean the location string by removing unwanted characters.

    ex. location: join{[<133:164](+), [344:400](+), [541:572](+)}
    returns: ['133:164', '344:400', '541:572']
    """

    # Remove 'join{' and '}' from the location string
    location = str(location)
    cleaned_location = location.split("{")[-1]
    
    # Remove unwanted characters
    for char in "}[]<>+-(){ ":
        cleaned_location = cleaned_location.replace(char, "")
    
    # Separete by commas to get individual intervals
    cleaned_location = cleaned_location.split(",")
    
    return cleaned_location

def make_exons_intervals_list(location):
    """
    Create a list of exon sequences from the split sequences.

    - location: A string representing the sequence with exon markers.
    ex. location: join{[<133:164](+), [344:400](+), [541:572](+)}
    returns: A list of exon sequences.
    ex. returns: [[133, 164], [344, 400], [541, 572]]
    """

    exons_intervals = []
    # Split the location string to extract exon parts
    intervals = clean_location_string(location)
    
    for interval in intervals:
        bonds = interval.split(":")
        start = int(bonds[0])
        end = int(bonds[1]) - 1
        exons_intervals.append([start, end])
        
    return exons_intervals

def make_introns_intervals_list(exons_intervals, seq_length):
    """
    Create a list of intron sequences from the split sequences.

    - exons: A list of exon sequences.
    ex. exons: [[133, 164], [344, 400], [541, 572]]
    returns: A list of intron sequences.
    ex. returns: [[0, 132], [165, 343], [401, 540], [573, seq_length-1]]
    """

    introns_intervals = []
        
    # If the first exon does not start at 0
    if exons_intervals[0][0] > 0: 
        introns_intervals.append([0, exons_intervals[0][0] - 1])
    
    for i in range(len(exons_intervals) - 1):   
        intron_start = exons_intervals[i][1] + 1
        intron_end = exons_intervals[i + 1][0] - 1
        introns_intervals.append([intron_start, intron_end])
    
    # If the last exon does not end at the last position of the sequence
    if exons_intervals[-1][1] < seq_length:
        introns_intervals.append([exons_intervals[-1][1] + 1, seq_length - 1])
    
    return introns_intervals

def make_exons_list(exons_intervals, seq):
    """
    Create a list of intron sequences from the split sequences.
    
    - exons_intervals: A list of exon intervals.
    ex. exons_intervals: [[133, 164], [344, 400], [541, 572]]
    returns: A list of exon sequences.
    ex. returns: ['ATG...TAA', 'GGC...TGA', 'CCT...TAG']
    """

    exons = []
    for exon_interval in exons_intervals:
        exons.append(seq[exon_interval[0]:exon_interval[1]+1])
        
    return exons

def make_introns_list(introns_intervals, seq):
    """
    Create a list of intron sequences from the split sequences. 
    For verification, we see if it starts with 'GT' and ends with 'AG'.
    
    - introns_intervals: A list of intron intervals.
    ex. introns_intervals: [[0, 132], [165, 343], [401, 540], [573, seq_length-1]]
    returns: A list of intron sequences.
    ex. returns: ['GTA...CAG', 'TTC...GGA', 'AAG...TTC', 'GGC...AAT']
    """

    introns = []

    # Verify the intron sequences based on 'GT' and 'AG' rules
    for intron_interval in introns_intervals:
        start = seq[intron_interval[0]:intron_interval[0]+2]
        end = seq[intron_interval[1]-1:intron_interval[1]+1] 

        # Caso padrão GT...AG
        if start == "GT" and end == "AG":
           introns.append(seq[intron_interval[0]:intron_interval[1]+1])
        # Caso especial ...AG
        elif start == 0 and end == "AG":
           introns.append(seq[intron_interval[0]:intron_interval[1]+1])
        # Caso especial GT...
        elif start == "GT" and end == 0:
           introns.append(seq[intron_interval[0]:intron_interval[1]+1])
        
    return introns


