# codon_usage_homosapiens = [
#     ['UUU', 'F', 0.49, 17.08,  1341983],  ['UCU', 'S', 0.20, 16.88,  1326337], ['UAU', 'Y', 0.47, 12.06,   947642],   ['UGU', 'C', 0.49, 10.48, 823726],
#     ['UUC', 'F', 0.51, 17.53,  1377246],  ['UCC', 'S', 0.20, 17.35,  1363614], ['UAC', 'Y', 0.53, 13.46,   1057591],  ['UGC', 'C', 0.51, 10.9, 856646],
#     ['UUA', 'L', 0.09,  8.78,  690119],   ['UCA', 'S', 0.16, 14.15,  1111829], ['UAA', '*', 0.28,  0.43,    33559],   ['UGA', '*', 0.50,  0.77,  60621],
#     ['UUG', 'L', 0.14, 13.42,  1054214],  ['UCG', 'S', 0.05,  4.08,  320758],  ['UAG', '*', 0.22,  0.34,    26407],   ['UGG', 'W', 1.00, 11.67, 916734],
    
#     ['CUU', 'L', 0.14, 14.11,  1108931],  ['CCU', 'P', 0.30, 19.07,  1498271],  ['CAU', 'H', 0.45, 11.87,  932958],   ['CGU', 'R', 0.08,  4.56, 358155],
#     ['CUC', 'L', 0.18, 17.85,  1402940],  ['CCC', 'P', 0.30, 18.98,  1491270],  ['CAC', 'H', 0.55, 14.65,  1151545],  ['CGC', 'R', 0.16, 8.78, 690020],
#     ['CUA', 'L', 0.08,  7.49,  588397],   ['CCA', 'P', 0.30, 18.67,  1466930],  ['CAA', 'Q', 0.29, 14.15,  1111835],  ['CGA', 'R', 0.11,  6.39, 502451],
#     ['CUG', 'L', 0.37, 36.24, 2848162],   ['CCG', 'P', 0.10,  6.19,  486595],   ['CAG', 'Q', 0.71, 35.42, 2783737],   ['CGG', 'R', 0.19, 10.70, 841125],
    
#     ['AUU', 'I', 0.38, 16.50,  1296263],  ['ACU', 'T', 0.26, 14.26,  1120650], ['AAU', 'N', 0.50, 18.40,   1445838],  ['AGU', 'S', 0.16, 13.98, 1098334],
#     ['AUC', 'I', 0.43, 18.67,  1467275],  ['ACC', 'T', 0.33, 17.80,  1398865], ['AAC', 'N', 0.50, 18.22,   1431977],  ['AGC', 'S', 0.23, 19.70, 1547892],
#     ['AUA', 'I', 0.19,  8.14,  639457],   ['ACA', 'T', 0.30, 16.54,  1299919], ['AAA', 'K', 0.46, 27.65,   2172879],  ['AGA', 'R', 0.24, 13.32, 1046695],
#     ['AUG', 'M', 1.00, 21.39,  1680759],  ['ACG', 'T', 0.10,  5.65,  443078],  ['AAG', 'K', 0.54, 31.85,  2503111],   ['AGG', 'R', 0.22, 12.19, 958126],
    
#     ['GUU', 'V', 0.20, 11.72,  920619],   ['GCU', 'A', 0.28, 18.85,  1481387],  ['GAU', 'D', 0.50, 23.99,  1885552],  ['GGU', 'G', 0.16, 10.77, 846406],
#     ['GUC', 'V', 0.23, 13.44,  1055851],  ['GCC', 'A', 0.38, 25.84, 2030499],   ['GAC', 'D', 0.50, 24.29, 1908799],   ['GGC', 'G', 0.30, 19.87, 1558404],
#     ['GUA', 'V', 0.13,  7.67,  603076],   ['GCA', 'A', 0.25, 17.05,  1339676],  ['GAA', 'E', 0.46, 33.72, 2649394],   ['GGA', 'G', 0.31, 17.07, 1641033],
#     ['GUG', 'V', 0.44, 25.90, 2035509],   ['GCG', 'A', 0.09,  6.03,  473616],   ['GAG', 'E', 0.54, 39.75, 3123257],   ['GGG', 'G', 0.23, 15.33, 1204755],
# ]

codon_usage_homosapiens = [
    #https://www.biologicscorp.com/tools/CodonUsage
    ['UUU', 'F', 0.46, 17.6, 714298],  ['UCU', 'S', 0.19, 15.2, 618711],  ['UAU', 'Y', 0.44, 12.2, 495699],  ['UGU', 'C', 0.46, 10.6, 430311],
    ['UUC', 'F', 0.54, 20.3, 824692],  ['UCC', 'S', 0.22, 17.7, 718892],  ['UAC', 'Y', 0.56, 15.3, 622407],  ['UGC', 'C', 0.54, 12.6, 513028],
    ['UUA', 'L', 0.08,  7.7, 311881],  ['UCA', 'S', 0.15, 12.2, 496448],  ['UAA', '*', 0.30,  1.0,  40285],  ['UGA', '*', 0.47,  1.6,  63237],
    ['UUG', 'L', 0.13, 12.9, 525688],  ['UCG', 'S', 0.05,  4.4, 179419],  ['UAG', '*', 0.24,  0.8,  32109],  ['UGG', 'W', 1.00, 13.2, 535595],
    
    ['CUU', 'L', 0.13, 13.2, 536515],  ['CCU', 'P', 0.29, 17.5, 713233],  ['CAU', 'H', 0.42, 10.9, 441711],  ['CGU', 'R', 0.08,  4.5, 184609],
    ['CUC', 'L', 0.20, 19.6, 796638],  ['CCC', 'P', 0.32, 19.8, 804620],  ['CAC', 'H', 0.58, 15.1, 613713],  ['CGC', 'R', 0.18, 10.4, 423516],
    ['CUA', 'L', 0.07,  7.2, 290751],  ['CCA', 'P', 0.28, 16.9, 688038],  ['CAA', 'Q', 0.27, 12.3, 501911],  ['CGA', 'R', 0.11,  6.2, 250760],
    ['CUG', 'L', 0.40, 39.6, 1611801], ['CCG', 'P', 0.11,  6.9, 281570],  ['CAG', 'Q', 0.73, 34.2, 1391973], ['CGG', 'R', 0.20, 11.4, 464485],
    
    ['AUU', 'I', 0.36, 16.0, 650473],  ['ACU', 'T', 0.25, 13.1, 533609],  ['AAU', 'N', 0.47, 17.0, 689701],  ['AGU', 'S', 0.15, 12.1, 493429],
    ['AUC', 'I', 0.47, 20.8, 846466],  ['ACC', 'T', 0.36, 18.9, 768147],  ['AAC', 'N', 0.53, 19.1, 776603],  ['AGC', 'S', 0.24, 19.5, 791383],
    ['AUA', 'I', 0.17,  7.5, 304565],  ['ACA', 'T', 0.28, 15.1, 614523],  ['AAA', 'K', 0.43, 24.4, 993621],  ['AGA', 'R', 0.21, 12.2, 494682],
    ['AUG', 'M', 1.00, 22.0, 896005],  ['ACG', 'T', 0.11,  6.1, 246105],  ['AAG', 'K', 0.57, 31.9, 1295568], ['AGG', 'R', 0.21, 12.0, 486463],
    
    ['GUU', 'V', 0.18, 11.0, 448607],  ['GCU', 'A', 0.27, 18.4, 750096],  ['GAU', 'D', 0.46, 21.8, 885429],  ['GGU', 'G', 0.16, 10.8, 437126],
    ['GUC', 'V', 0.24, 14.5, 588138],  ['GCC', 'A', 0.40, 27.7, 1127679], ['GAC', 'D', 0.54, 25.1, 1020595], ['GGC', 'G', 0.34, 22.2, 903565],
    ['GUA', 'V', 0.12,  7.1, 287712],  ['GCA', 'A', 0.23, 15.8, 643471],  ['GAA', 'E', 0.42, 29.0, 1177632], ['GGA', 'G', 0.25, 16.5, 669873],
    ['GUG', 'V', 0.46, 28.1, 1143534], ['GCG', 'A', 0.11,  7.4, 299495],  ['GAG', 'E', 0.58, 39.6, 1609975], ['GGG', 'G', 0.25, 16.5, 669768],
]

codon_usage_ecoli = [
    ['UUU', 'F', 0.58, 22.40,  88573],  ['UCU', 'S', 0.15,  8.46,  33454], ['UAU', 'Y', 0.57, 16.17,   63949],  ['UGU', 'C', 0.45, 5.22, 20643],
    ['UUC', 'F', 0.42, 16.54,  65392],  ['UCC', 'S', 0.15,  8.64,  34157], ['UAC', 'Y', 0.43, 12.16,   48059],  ['UGC', 'C', 0.55, 6.48, 25625],
    ['UUA', 'L', 0.13, 13.95,  55171],  ['UCA', 'S', 0.12, 7.17,  28362],  ['UAA', '*', 0.64,  2.10,    8296],  ['UGA', '*', 0.29,  0.97,  3849],
    ['UUG', 'L', 0.13, 13.77,  54440],  ['UCG', 'S', 0.15,  8.89,  35153], ['UAG', '*', 0.07,  0.23,    903],   ['UGG', 'W', 1.00, 15.26, 60314],
    
    ['CUU', 'L', 0.10, 11.11,  43907],  ['CCU', 'P', 0.16, 7.01,  27703],  ['CAU', 'H', 0.57, 12.84,  50765],    ['CGU', 'R', 0.38,  20.94, 82796],
    ['CUC', 'L', 0.10, 11.06,  43733],  ['CCC', 'P', 0.12, 5.44,  21522],  ['CAC', 'H', 0.43, 9.65,  38150],     ['CGC', 'R', 0.40, 21.99, 86947],
    ['CUA', 'L', 0.04,  3.93,  15553],  ['CCA', 'P', 0.19, 8.49,  33554],  ['CAA', 'Q', 0.35, 15.37,  60750],    ['CGA', 'R', 0.07,  3.60, 14242],
    ['CUG', 'L', 0.49, 52.75, 208539],  ['CCG', 'P', 0.53,  23.35,  92332],  ['CAG', 'Q', 0.65, 28.70, 113461],  ['CGG', 'R', 0.09, 5.21, 20580],
    
    ['AUU', 'I', 0.51, 30.64,  121134],  ['ACU', 'T', 0.17, 8.98,  35521],   ['AAU', 'N', 0.45, 17.66,   69838],   ['AGU', 'S', 0.15, 8.70, 34392],
    ['AUC', 'I', 0.42, 25.14,  99402],   ['ACC', 'T', 0.43, 23.25,  91927],  ['AAC', 'N', 0.55, 21.53,   85123],   ['AGC', 'S', 0.28, 16.00, 63251],
    ['AUA', 'I', 0.07,  4.37,  17287],   ['ACA', 'T', 0.13, 7.03,  27804],   ['AAA', 'K', 0.77, 33.68,   133176],  ['AGA', 'R', 0.04, 2.04, 8050],
    ['AUG', 'M', 1.00, 27.89,  110250],  ['ACG', 'T', 0.27,  14.38,  56860], ['AAG', 'K', 0.23, 10.23,  40444],    ['AGG', 'R', 0.02, 1.15, 4562],
    
    ['GUU', 'V', 0.26, 18.40,  72758],  ['GCU', 'A', 0.16, 15.36,  60728],    ['GAU', 'D', 0.63, 31.96,  126351],  ['GGU', 'G', 0.34, 24.85, 98251],
    ['GUC', 'V', 0.21, 15.17,  59967],  ['GCC', 'A', 0.27, 25.51, 100864],    ['GAC', 'D', 0.37, 18.97, 74987],    ['GGC', 'G', 0.40, 29.65, 117205],
    ['GUA', 'V', 0.16, 11.03,  43593],  ['GCA', 'A', 0.21, 20.35,  80448],    ['GAA', 'E', 0.69, 39.50, 156171],   ['GGA', 'G', 0.11, 7.93, 31353],
    ['GUG', 'V', 0.37, 26.32, 104043],  ['GCG', 'A', 0.36,  33.86,  133855],  ['GAG', 'E', 0.31, 17.63, 69702],    ['GGG', 'G', 0.15, 10.98, 43422],
]

codon_usages = {
    'homosapiens': codon_usage_homosapiens,
    'ecoli': codon_usage_ecoli,
}
list_seqs = {
    'list_codon': ['GCU', 'GCC', 'GCA', 'GCG', 'UGU', 'UGC', 'GAU', 'GAC',
                   'GAA', 'GAG', 'UUU', 'UUC', 'GGU', 'GGC', 'GGA', 'GGG',
                   'CAU', 'CAC', 'AUU', 'AUC', 'AUA', 'AAA', 'AAG', 'UUA',
                   'UUG', 'CUU', 'CUC', 'CUA', 'CUG', 'AUG', 'AAU', 'AAC',
                   'CCU', 'CCC', 'CCA', 'CCG', 'CAA', 'CAG', 'CGU', 'CGC',
                   'CGA', 'CGG', 'AGA', 'AGG', 'UCU', 'UCC', 'UCA', 'UCG',
                   'AGU', 'AGC', 'ACU', 'ACC', 'ACA', 'ACG', 'GUU', 'GUC',
                   'GUA', 'GUG', 'UGG', 'UAU', 'UAC', 'UAA', 'UAG', 'UGA'],
    'list_pad' : ['PAD'],
    'list_aminoacid' : ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 
                  'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y', '*'],
}

codon_to_aa = {
        # RNA to Amino Acid
    'GCU': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
    'UGU': 'C', 'UGC': 'C',
    'GAU': 'D', 'GAC': 'D',
    'GAA': 'E', 'GAG': 'E',
    'UUU': 'F', 'UUC': 'F',
    'GGU': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G',
    'CAU': 'H', 'CAC': 'H',
    'AUU': 'I', 'AUC': 'I', 'AUA': 'I',
    'AAA': 'K', 'AAG': 'K',
    'UUA': 'L', 'UUG': 'L', 'CUU': 'L', 'CUC': 'L', 'CUA': 'L', 'CUG': 'L',
    'AUG': 'M',
    'AAU': 'N', 'AAC': 'N',
    'CCU': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
    'CAA': 'Q', 'CAG': 'Q',
    'CGU': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R', 'AGA': 'R', 'AGG': 'R',
    'UCU': 'S', 'UCC': 'S', 'UCA': 'S', 'UCG': 'S', 'AGU': 'S', 'AGC': 'S',
    'ACU': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
    'GUU': 'V', 'GUC': 'V', 'GUA': 'V', 'GUG': 'V',
    'UGG': 'W',
    'UAU': 'Y', 'UAC': 'Y',
    'UAA': '*', 'UAG': '*', 'UGA': '*',
}

aa_to_codon = {
      'A': ['GCU', 'GCC', 'GCA', 'GCG'],
      'C': ['UGU', 'UGC'],
      'D': ['GAU', 'GAC'],
      'E': ['GAA', 'GAG'],
      'F': ['UUU', 'UUC'],
      'G': ['GGU', 'GGC', 'GGA', 'GGG'],
      'H': ['CAU', 'CAC'],
      'I': ['AUU', 'AUC', 'AUA'],
      'K': ['AAA', 'AAG'],
      'L': ['UUA', 'UUG', 'CUU', 'CUC', 'CUA', 'CUG'],
      'M': ['AUG'],
      'N': ['AAU', 'AAC'],
      'P': ['CCU', 'CCC', 'CCA', 'CCG'],
      'Q': ['CAA', 'CAG'],
      'R': ['CGU', 'CGC', 'CGA', 'CGG', 'AGA', 'AGG'],
      'S': ['UCU', 'UCC', 'UCA', 'UCG', 'AGU', 'AGC'],
      'T': ['ACU', 'ACC', 'ACA', 'ACG'],
      'V': ['GUU', 'GUC', 'GUA', 'GUG'],
      'W': ['UGG'],
      'Y': ['UAU', 'UAC'],
      '*': ['UAA', 'UAG', 'UGA'],
}


genetic_code = {
    # RNA to Amino Acid
    'GCU': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
    'UGU': 'C', 'UGC': 'C',
    'GAU': 'D', 'GAC': 'D',
    'GAA': 'E', 'GAG': 'E',
    'UUU': 'F', 'UUC': 'F',
    'GGU': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G',
    'CAU': 'H', 'CAC': 'H',
    'AUU': 'I', 'AUC': 'I', 'AUA': 'I',
    'AAA': 'K', 'AAG': 'K',
    'UUA': 'L', 'UUG': 'L', 'CUU': 'L', 'CUC': 'L', 'CUA': 'L', 'CUG': 'L',
    'AUG': 'M',
    'AAU': 'N', 'AAC': 'N',
    'CCU': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
    'CAA': 'Q', 'CAG': 'Q',
    'CGU': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R', 'AGA': 'R', 'AGG': 'R',
    'UCU': 'S', 'UCC': 'S', 'UCA': 'S', 'UCG': 'S', 'AGU': 'S', 'AGC': 'S',
    'ACU': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
    'GUU': 'V', 'GUC': 'V', 'GUA': 'V', 'GUG': 'V',
    'UGG': 'W',
    'UAU': 'Y', 'UAC': 'Y',
    'UAA': '*', 'UAG': '*', 'UGA': '*',

    # PAD
    'PAD': 'PAD',

    # DNA to Amino Acid
    'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
    'TGT': 'C', 'TGC': 'C',
    'GAT': 'D', 'GAC': 'D',
    'GAA': 'E', 'GAG': 'E',
    'TTT': 'F', 'TTC': 'F',
    'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G',
    'CAT': 'H', 'CAC': 'H',
    'ATT': 'I', 'ATC': 'I', 'ATA': 'I',
    'AAA': 'K', 'AAG': 'K',
    'TTA': 'L', 'TTG': 'L', 'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',
    'ATG': 'M',
    'AAT': 'N', 'AAC': 'N',
    'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
    'CAA': 'Q', 'CAG': 'Q',
    'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R', 'AGA': 'R', 'AGG': 'R',
    'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S', 'AGT': 'S', 'AGC': 'S',
    'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
    'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',
    'TGG': 'W',
    'TAT': 'Y', 'TAC': 'Y',
    'TAA': '*', 'TAG': '*', 'TGA': '*',

    'F': ['UUU', 'UUC'], 'L': ['UUA', 'UUG', 'CUU', 'CUC', 'CUA', 'CUG'],
    'I': ['AUU', 'AUC', 'AUA'], 'M': ['AUG'], 'V': ['GUU', 'GUC', 'GUA', 'GUG'],
    'S': ['UCU', 'UCC', 'UCA', 'UCG', 'AGU', 'AGC'], 'P': ['CCU', 'CCC', 'CCA', 'CCG'],
    'T': ['ACU', 'ACC', 'ACA', 'ACG'], 'A': ['GCU', 'GCC', 'GCA', 'GCG'],
    'Y': ['UAU', 'UAC'], 'H': ['CAU', 'CAC'], 'Q': ['CAA', 'CAG'],
    'N': ['AAU', 'AAC'], 'K': ['AAA', 'AAG'], 'D': ['GAU', 'GAC'],
    'E': ['GAA', 'GAG'], 'C': ['UGU', 'UGC'], 'W': ['UGG'],
    'R': ['CGU', 'CGC', 'CGA', 'CGG', 'AGA', 'AGG'], 'G': ['GGU', 'GGC', 'GGA', 'GGG'],
    '*': ['UAA', 'UAG', 'UGA']
}


label2index = {
    # A 
    'GCU': 0, 'GCC': 1, 'GCA': 2, 'GCG': 3,
    # C
    'UGU': 0, 'UGC': 1,
    # D
    'GAU': 0, 'GAC': 1,
    # E
    'GAA': 0, 'GAG': 1,
    # F
    'UUU': 0, 'UUC': 1,
    # G
    'GGU': 0, 'GGC': 1, 'GGA': 2, 'GGG': 3,
    # H
    'CAU': 0, 'CAC': 1,
    # I
    'AUU': 0, 'AUC': 1, 'AUA': 2,
    # K
    'AAA': 0, 'AAG': 1,
    # L
    'UUA': 0, 'UUG': 1, 'CUU': 2, 'CUC': 3, 'CUA': 4, 'CUG': 5,
    # M
    'AUG': 0,
    # N
    'AAU': 0, 'AAC': 1,
    # P
    'CCU': 0, 'CCC': 1, 'CCA': 2, 'CCG': 3,
    # Q
    'CAA': 0, 'CAG': 1,
    # R
    'CGU': 0, 'CGC': 1, 'CGA': 2, 'CGG': 3, 'AGA': 4, 'AGG': 5,
    # S
    'UCU': 0, 'UCC': 1, 'UCA': 2, 'UCG': 3, 'AGU': 4, 'AGC': 5,
    # T
    'ACU': 0, 'ACC': 1, 'ACA': 2, 'ACG': 3,
    # V
    'GUU': 0, 'GUC': 1, 'GUA': 2, 'GUG': 3,
    # W
    'UGG': 0,
    # Y
    'UAU': 0, 'UAC': 1,
    # *
    'UAA': 0, 'UAG': 1, 'UGA': 2,
    # PAD
    'PAD': 6,
}


seq2index = {
    
    'A': 0, 'C': 1, 'D': 2, 'E': 3, 'F': 4, 
    'G': 5, 'H': 6, 'I': 7, 'K': 8, 'L': 9, 
    'M': 10, 'N': 11, 'P': 12, 'Q': 13, 'R': 14, 
    'S': 15, 'T': 16, 'V': 17, 'W': 18, 'Y': 19, '*': 20,
    'PAD': 21,

    'GCU': 22, 'GCC': 23, 'GCA': 24, 'GCG': 25,
    'UGU': 26, 'UGC': 27,
    'GAU': 28, 'GAC': 29,
    'GAA': 30, 'GAG': 31,
    'UUU': 32, 'UUC': 33,
    'GGU': 34, 'GGC': 35, 'GGA': 36, 'GGG': 37,
    'CAU': 38, 'CAC': 39,
    'AUU': 40, 'AUC': 41, 'AUA': 42,
    'AAA': 43, 'AAG': 44,
    'UUA': 45, 'UUG': 46, 'CUU': 47, 'CUC': 48, 'CUA': 49, 'CUG': 50,
    'AUG': 51,
    'AAU': 52, 'AAC': 53,
    'CCU': 54, 'CCC': 55, 'CCA': 56, 'CCG': 57,
    'CAA': 58, 'CAG': 59,
    'CGU': 60, 'CGC': 61, 'CGA': 62, 'CGG': 63, 'AGA': 64, 'AGG': 65,
    'UCU': 66, 'UCC': 67, 'UCA': 68, 'UCG': 69, 'AGU': 70, 'AGC': 71,
    'ACU': 72, 'ACC': 73, 'ACA': 74, 'ACG': 75,
    'GUU': 76, 'GUC': 77, 'GUA': 78, 'GUG': 79,
    'UGG': 80,
    'UAU': 81, 'UAC': 82,
    'UAA': 83, 'UAG': 84, 'UGA': 85,
}


# 코돈 vocab id → 코돈 문자열
codon_id2triplet = {
    22: 'GCU', 23: 'GCC', 24: 'GCA', 25: 'GCG',
    26: 'UGU', 27: 'UGC',
    28: 'GAU', 29: 'GAC',
    30: 'GAA', 31: 'GAG',
    32: 'UUU', 33: 'UUC',
    34: 'GGU', 35: 'GGC', 36: 'GGA', 37: 'GGG',
    38: 'CAU', 39: 'CAC',
    40: 'AUU', 41: 'AUC', 42: 'AUA',
    43: 'AAA', 44: 'AAG',
    45: 'UUA', 46: 'UUG', 47: 'CUU', 48: 'CUC', 49: 'CUA', 50: 'CUG',
    51: 'AUG',
    52: 'AAU', 53: 'AAC',
    54: 'CCU', 55: 'CCC', 56: 'CCA', 57: 'CCG',
    58: 'CAA', 59: 'CAG',
    60: 'CGU', 61: 'CGC', 62: 'CGA', 63: 'CGG', 64: 'AGA', 65: 'AGG',
    66: 'UCU', 67: 'UCC', 68: 'UCA', 69: 'UCG', 70: 'AGU', 71: 'AGC',
    72: 'ACU', 73: 'ACC', 74: 'ACA', 75: 'ACG',
    76: 'GUU', 77: 'GUC', 78: 'GUA', 79: 'GUG',
    80: 'UGG',
    81: 'UAU', 82: 'UAC',
    83: 'UAA', 84: 'UAG', 85: 'UGA'
}




index2codon = {
    'A': {0: 'GCU', 
          1: 'GCC',
          2: 'GCA',
          3: 'GCG',
          4: 'PAD', 5: 'PAD', 6: 'PAD'},
    'C': {0: 'UGU',
          1: 'UGC',
          2: 'PAD', 3: 'PAD', 4: 'PAD', 5: 'PAD', 6: 'PAD'},
    'D': {0: 'GAU',
          1: 'GAC',
          2: 'PAD', 3: 'PAD', 4: 'PAD', 5: 'PAD', 6: 'PAD'},
    'E': {0: 'GAA',
          1: 'GAG',
          2: 'PAD', 3: 'PAD', 4: 'PAD', 5: 'PAD', 6: 'PAD'},
    'F': {0: 'UUU',
          1: 'UUC',
          2: 'PAD', 3: 'PAD', 4: 'PAD', 5: 'PAD', 6: 'PAD'},
    'G': {0: 'GGU',
          1: 'GGC',
          2: 'GGA',
          3: 'GGG',
          4: 'PAD', 5: 'PAD', 6: 'PAD'},
    'H': {0: 'CAU',
          1: 'CAC',
          2: 'PAD', 3: 'PAD', 4: 'PAD', 5: 'PAD', 6: 'PAD'},
    'I': {0: 'AUU',
          1: 'AUC',
          2: 'AUA',
          3: 'PAD', 4: 'PAD', 5: 'PAD', 6: 'PAD'},
    'K': {0: 'AAA',
          1: 'AAG',
          2: 'PAD', 3: 'PAD', 4: 'PAD', 5: 'PAD', 6: 'PAD'},
    'L': {0: 'UUA',
          1: 'UUG',
          2: 'CUU',
          3: 'CUC',
          4: 'CUA',
          5: 'CUG',
          6: 'PAD'},
    'M': {0: 'AUG',
          1: 'AUG',
          2: 'PAD', 3: 'PAD', 4: 'PAD', 5: 'PAD', 6: 'PAD'},
    'N': {0: 'AAU',
          1: 'AAC',
          2: 'PAD', 3: 'PAD', 4: 'PAD', 5: 'PAD', 6: 'PAD'},
    'P': {0: 'CCU',
          1: 'CCC',
          2: 'CCA',
          3: 'CCG', 
          4: 'PAD', 5: 'PAD', 6: 'PAD'},
    'Q': {0: 'CAA',
          1: 'CAG', 
          2: 'PAD', 3: 'PAD', 4: 'PAD', 5: 'PAD', 6: 'PAD'},
    'R': {0: 'CGU', 
          1: 'CGC', 
          2: 'CGA', 
          3: 'CGG', 
          4: 'AGA', 
          5: 'AGG',
          6: 'PAD'},
    'S': {0: 'UCU',
          1: 'UCC',
          2: 'UCA',
          3: 'UCG',
          4: 'AGU',
          5: 'AGC',
          6: 'PAD'},
    'T': {0: 'ACU',
          1: 'ACC',
          2: 'ACA',
          3: 'ACG',
          4: 'PAD', 5: 'PAD', 6: 'PAD'},
    'V': {0: 'GUU',
          1: 'GUC',
          2: 'GUA',
          3: 'GUG',
          4: 'PAD', 5: 'PAD', 6: 'PAD'},
    'W': {0: 'UGG',
          1: 'UGG',
          2: 'PAD', 3: 'PAD', 4: 'PAD', 5: 'PAD', 6: 'PAD'},
    'Y': {0: 'UAU',
          1: 'UAC', 
          2: 'PAD', 3: 'PAD', 4: 'PAD', 5: 'PAD', 6: 'PAD'},
    '*': {0: 'UAA',
          1: 'UAG',
          2: 'UGA',
          3: 'PAD', 4: 'PAD', 5: 'PAD', 6: 'PAD'},
    'PAD': {0: 'PAD', 1: 'PAD', 2: 'PAD', 3: 'PAD', 4: 'PAD', 5: 'PAD', 6: 'PAD'},
}


index2aminoacid = {
    0: 'A', 1: 'C', 2: 'D', 3: 'E', 4: 'F',
    5: 'G', 6: 'H', 7: 'I', 8: 'K', 9: 'L', 
    10: 'M', 11: 'N', 12: 'P', 13: 'Q', 14: 'R', 
    15: 'S', 16: 'T', 17: 'V', 18: 'W', 19: 'Y', 
    20: '*', 21: 'PAD', 
}


genetic_dictionary_7 = {
    'usages': codon_usages,
    'list_seqs': list_seqs,
    'codon_to_aa': codon_to_aa,
    'aa_to_codon': aa_to_codon,
    'genetic_code': genetic_code,
    'seq2index': seq2index,
    'label2index': label2index,
    'index2codon': index2codon,
    'index2aminoacid': index2aminoacid,
}


model_dict = {
    'dim_in': 32,
    'dim_attn': 16,
    'dim_out': 7,
    'n_heads': 4,
    'dropout_attn': 0.1,
    'dropout_ff': 0.1,
    'mult_ff': 2,
    'depth': 3,
    'name': 'coformer'
}

train_dict = {
    'batchsize': 512,
    'epoch': 200,
    'optimizer': 'adamw',
    'lr': 1e-3,
    'wo_pad': True,
    'verbose': True,
}

pretrain_dict = {
    'batchsize': 512,
    'epochs': 200,
    'optimizer': 'adamw',
    'lr': 1e-3,
    'mlm_codon': True,
    'mlm_aa': True,
    'contrastive': True,
    'wo_pad': True,
    'verbose': True,
}

data_dict = {
    'data_name':'homosapiens',
    'clip_size': 128,
    'batchsize': 512,
    'species': 'homosapiens',
    'size_info': [
        [0, 128, 128, 256],
        [256, 512, 512, 64],
        [512, 100000, 1024, 16]
    ]
}

init_dict = {
    'file_name': 'best_model',
    'seed':0,
    'device':-1,
}

param_dict = {
    'model': model_dict,
    'train': train_dict,
    'pretrain': pretrain_dict,
    'data': data_dict,
    'init': init_dict,
}





index2codon_global = {
    0: {'amino_acid': 'A', 'codon': 'GCU', 'local_index': 0},
    1: {'amino_acid': 'A', 'codon': 'GCC', 'local_index': 1},
    2: {'amino_acid': 'A', 'codon': 'GCA', 'local_index': 2},
    3: {'amino_acid': 'A', 'codon': 'GCG', 'local_index': 3},
    4: {'amino_acid': 'C', 'codon': 'UGU', 'local_index': 0},
    5: {'amino_acid': 'C', 'codon': 'UGC', 'local_index': 1},
    6: {'amino_acid': 'D', 'codon': 'GAU', 'local_index': 0},
    7: {'amino_acid': 'D', 'codon': 'GAC', 'local_index': 1},
    8: {'amino_acid': 'E', 'codon': 'GAA', 'local_index': 0},
    9: {'amino_acid': 'E', 'codon': 'GAG', 'local_index': 1},
    10: {'amino_acid': 'F', 'codon': 'UUU', 'local_index': 0},
    11: {'amino_acid': 'F', 'codon': 'UUC', 'local_index': 1},
    12: {'amino_acid': 'G', 'codon': 'GGU', 'local_index': 0},
    13: {'amino_acid': 'G', 'codon': 'GGC', 'local_index': 1},
    14: {'amino_acid': 'G', 'codon': 'GGA', 'local_index': 2},
    15: {'amino_acid': 'G', 'codon': 'GGG', 'local_index': 3},
    16: {'amino_acid': 'H', 'codon': 'CAU', 'local_index': 0},
    17: {'amino_acid': 'H', 'codon': 'CAC', 'local_index': 1},
    18: {'amino_acid': 'I', 'codon': 'AUU', 'local_index': 0},
    19: {'amino_acid': 'I', 'codon': 'AUC', 'local_index': 1},
    20: {'amino_acid': 'I', 'codon': 'AUA', 'local_index': 2},
    21: {'amino_acid': 'K', 'codon': 'AAA', 'local_index': 0},
    22: {'amino_acid': 'K', 'codon': 'AAG', 'local_index': 1},
    23: {'amino_acid': 'L', 'codon': 'UUA', 'local_index': 0},
    24: {'amino_acid': 'L', 'codon': 'UUG', 'local_index': 1},
    25: {'amino_acid': 'L', 'codon': 'CUU', 'local_index': 2},
    26: {'amino_acid': 'L', 'codon': 'CUC', 'local_index': 3},
    27: {'amino_acid': 'L', 'codon': 'CUA', 'local_index': 4},
    28: {'amino_acid': 'L', 'codon': 'CUG', 'local_index': 5},
    29: {'amino_acid': 'M', 'codon': 'AUG', 'local_index': 0},
    30: {'amino_acid': 'N', 'codon': 'AAU', 'local_index': 0},
    31: {'amino_acid': 'N', 'codon': 'AAC', 'local_index': 1},
    32: {'amino_acid': 'P', 'codon': 'CCU', 'local_index': 0},
    33: {'amino_acid': 'P', 'codon': 'CCC', 'local_index': 1},
    34: {'amino_acid': 'P', 'codon': 'CCA', 'local_index': 2},
    35: {'amino_acid': 'P', 'codon': 'CCG', 'local_index': 3},
    36: {'amino_acid': 'Q', 'codon': 'CAA', 'local_index': 0},
    37: {'amino_acid': 'Q', 'codon': 'CAG', 'local_index': 1},
    38: {'amino_acid': 'R', 'codon': 'CGU', 'local_index': 0},
    39: {'amino_acid': 'R', 'codon': 'CGC', 'local_index': 1},
    40: {'amino_acid': 'R', 'codon': 'CGA', 'local_index': 2},
    41: {'amino_acid': 'R', 'codon': 'CGG', 'local_index': 3},
    42: {'amino_acid': 'R', 'codon': 'AGA', 'local_index': 4},
    43: {'amino_acid': 'R', 'codon': 'AGG', 'local_index': 5},
    44: {'amino_acid': 'S', 'codon': 'UCU', 'local_index': 0},
    45: {'amino_acid': 'S', 'codon': 'UCC', 'local_index': 1},
    46: {'amino_acid': 'S', 'codon': 'UCA', 'local_index': 2},
    47: {'amino_acid': 'S', 'codon': 'UCG', 'local_index': 3},
    48: {'amino_acid': 'S', 'codon': 'AGU', 'local_index': 4},
    49: {'amino_acid': 'S', 'codon': 'AGC', 'local_index': 5},
    50: {'amino_acid': 'T', 'codon': 'ACU', 'local_index': 0},
    51: {'amino_acid': 'T', 'codon': 'ACC', 'local_index': 1},
    52: {'amino_acid': 'T', 'codon': 'ACA', 'local_index': 2},
    53: {'amino_acid': 'T', 'codon': 'ACG', 'local_index': 3},
    54: {'amino_acid': 'V', 'codon': 'GUU', 'local_index': 0},
    55: {'amino_acid': 'V', 'codon': 'GUC', 'local_index': 1},
    56: {'amino_acid': 'V', 'codon': 'GUA', 'local_index': 2},
    57: {'amino_acid': 'V', 'codon': 'GUG', 'local_index': 3},
    58: {'amino_acid': 'W', 'codon': 'UGG', 'local_index': 0},
    59: {'amino_acid': 'Y', 'codon': 'UAU', 'local_index': 0},
    60: {'amino_acid': 'Y', 'codon': 'UAC', 'local_index': 1},
    61: {'amino_acid': '*', 'codon': 'UAA', 'local_index': 0},
    62: {'amino_acid': '*', 'codon': 'UAG', 'local_index': 1},
    63: {'amino_acid': '*', 'codon': 'UGA', 'local_index': 2},
}




aa2global_index = {
    ('A', 0): 0,
    ('A', 1): 1,
    ('A', 2): 2,
    ('A', 3): 3,
    ('C', 0): 4,
    ('C', 1): 5,
    ('D', 0): 6,
    ('D', 1): 7,
    ('E', 0): 8,
    ('E', 1): 9,
    ('F', 0): 10,
    ('F', 1): 11,
    ('G', 0): 12,
    ('G', 1): 13,
    ('G', 2): 14,
    ('G', 3): 15,
    ('H', 0): 16,
    ('H', 1): 17,
    ('I', 0): 18,
    ('I', 1): 19,
    ('I', 2): 20,
    ('K', 0): 21,
    ('K', 1): 22,
    ('L', 0): 23,
    ('L', 1): 24,
    ('L', 2): 25,
    ('L', 3): 26,
    ('L', 4): 27,
    ('L', 5): 28,
    ('M', 0): 29,
    ('N', 0): 30,
    ('N', 1): 31,
    ('P', 0): 32,
    ('P', 1): 33,
    ('P', 2): 34,
    ('P', 3): 35,
    ('Q', 0): 36,
    ('Q', 1): 37,
    ('R', 0): 38,
    ('R', 1): 39,
    ('R', 2): 40,
    ('R', 3): 41,
    ('R', 4): 42,
    ('R', 5): 43,
    ('S', 0): 44,
    ('S', 1): 45,
    ('S', 2): 46,
    ('S', 3): 47,
    ('S', 4): 48,
    ('S', 5): 49,
    ('T', 0): 50,
    ('T', 1): 51,
    ('T', 2): 52,
    ('T', 3): 53,
    ('V', 0): 54,
    ('V', 1): 55,
    ('V', 2): 56,
    ('V', 3): 57,
    ('W', 0): 58,
    ('Y', 0): 59,
    ('Y', 1): 60,
    ('*', 0): 61,
    ('*', 1): 62,
    ('*', 2): 63
    }

