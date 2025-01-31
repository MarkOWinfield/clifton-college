# -*- coding: utf-8 -*-
"""
Created on Tue Dec  3 13:06:03 2024

@author: Mark Winfield
"""

import streamlit as st
from random import choice

st.set_page_config(layout="wide")

st.markdown(""" <style> .font1 {
    font-size: 45px; font-family: 'Copper Black'; color: #FF9633}
    </style> """, unsafe_allow_html=True)


###############################################################################
################################ FUNCTIONS ####################################
###############################################################################

def chunkstring(string, length):
    return (string[0+i:length+i] for i in range(0, len(string), length))


def translate(codons):
    
    converter = {'TTT': 'Phe',
                 'TTC': 'Phe',
                 'TTA': 'Leu',
                 'TTG': 'Leu',
                 'CTT': 'Leu',
                 'CTC': 'Leu',
                 'CTA': 'Leu',
                 'CTG': 'Leu',
                 'ATT': 'Ile',
                 'ATC': 'Ile',
                 'ATA': 'Ile',
                 'ATG': 'Met',
                 'GTT': 'Val',
                 'GTC': 'Val',
                 'GTA': 'Val',
                 'GTG': 'Val',
                 'TCT': 'Ser',
                 'TCC': 'Ser',
                 'TCA': 'Ser',
                 'TCG': 'Ser',
                 'CCT': 'Pro',
                 'CCC': 'Pro',
                 'CCA': 'Pro',
                 'CCG': 'Pro',
                 'ACT': 'Thr',
                 'ACC': 'Thr',
                 'ACA': 'Thr',
                 'ACG': 'Thr',
                 'GCT': 'Ala',
                 'GCC': 'Ala',
                 'GCA': 'Ala',
                 'GCG': 'Ala',
                 'TAT': 'Tyr',
                 'TAC': 'Tyr',
                 'TAA': 'STOP',
                 'TAG': 'STOP',
                 'CAT': 'His',
                 'CAC': 'His',
                 'CAA': 'Glu',
                 'CAG': 'Glu',
                 'AAT': 'Asn',
                 'AAC': 'Asn',
                 'AAA': 'Lys',
                 'AAG': 'Lys',
                 'GAT': 'Asp',
                 'GAC': 'Asp',
                 'GAA': 'Glu',
                 'GAG': 'Glu',
                 'TGT': 'Cys',
                 'TGC': 'Cys',
                 'TGA': 'STOP',
                 'TGG': 'Trp',
                 'CGT': 'Arg',
                 'CGC': 'Arg',
                 'CGA': 'Arg',
                 'CGG': 'Arg',
                 'AGT': 'Ser',
                 'AGC': 'Ser',
                 'AGA': 'Arg',
                 'AGG': 'Arg',
                 'GGT': 'Gly',
                 'GGC': 'Gly',
                 'GGA': 'Gly',
                 'GGG': 'Gly'             
                 }

    aa = []

    for i in codons:
        if i in ['TGA', 'TAA', 'TAG']:
            aa.append(' :red[STOP]')
            break
        if len(i) == 3:
            aa.append(converter[i])
        else:
            continue
        
    return(aa)


def remove_short_codon(codons):
    
    for i in codons:
        if len(i) < 3:
            codons.remove(i)
            
    codon_string = ' '.join(codons)
            
    return(codon_string)


###############################################################################
###############################################################################


#st.markdown('<h1 class="font1">Quick Check Aliens vs Reference</h1>', unsafe_allow_html=True) 

st.title(' :blue[Transcribe the DNA into mRNA]')

# The random sequence created on by the 'Sequence Genertor' is passed to this
# page so that open reading frames can be found.
sequence = st.session_state['sequence']
st.write(sequence)

st.markdown('''
            Identify open reading frames in the DNA sequence:
            1. Locate a sequence corresponding to a start codon (usually ATG) - this establishes 
            the reading frame.
            2. Read this sequence in base triplets until a stop codon (TGA, TAG or TAA) is reached.
            3. The longer the sequence, the greater the likelihood that it represents
            a genuine open reading frame and, therefore, a gene.

            Gene sequences are largely conserved – so if an ORF sequence is present 
            in multiple genomes, it likely represents a gene.
            '''
            )



st.markdown('#### Find the ORFs')

# st.write(f'The number of Open Reading Frames (ORFs) in this sequence is {sequence.count('ATG')}')

ORF_count = sequence.count('ATG')
met1 = sequence.find('ATG')

if ORF_count == 1: 
    
    st.write(f'There is only one ORF in this sequence; it is a position {met1 + 1}.  Using this reading frame, the codons are as follows:')
    
    ORF = sequence[met1:]
    
    codons = list(chunkstring(ORF, 3))
    
    codon_string = remove_short_codon(codons)
    
    st.write(codon_string)
    
    RNA_string = codon_string.replace('T', 'U')
    
    st.markdown('#### Transcribe the DNA into mRNA')

    st.write(RNA_string)
   
    
elif ORF_count > 1:
    
    st.write(f'There are {ORF_count} ORFs in this sequence; the first is at position {met1 + 1}.  Using this reading frame, the codons are as follows:')
    
    ORF = sequence[met1:]
    # st.write(ORF)
    
    codons = list(chunkstring(ORF, 3))
    
    codon_string = remove_short_codon(codons)
    
    st.write(codon_string)
    
    RNA_string = codon_string.replace('T', 'U')
    
    st.markdown('#### Transcribe the DNA into mRNA')

    st.write(RNA_string)
   
    aa = translate(codons)
    

else:

    st.write('The are no ORFs in this sequence.  Try generating a longer sequence')
    
    
