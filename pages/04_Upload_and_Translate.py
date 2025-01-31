# -*- coding: utf-8 -*-
"""
Created on Tue Dec  3 13:06:03 2024

@author: Mark Winfield
"""

import streamlit as st
import pandas as pd
from io import StringIO

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
        if '\r' in i:
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

st.title(' :blue[Upload Fasta file and Translate]')


st.markdown('''
            Identify an open reading frame in an uploaded sequence:
            1. Locate a sequence corresponding to a start codon (usually ATG) - this establishes 
            the reading frame.
            2. Read this sequence in base triplets until a stop codon (TGA, TAG or TAA) is reached.
            3. The longer the sequence, the greater the likelihood that it represents
            a genuine open reading frame and, therefore, a gene.
            4. Until further proofs of its veracity, the '*gene*' may be labelled as a
            '*:orange[hypothetical protein coding]*'' sequence.

            Certain bioinformatic programs can automatically identify potential ORFs 
            when provided with a candidate sequence.

            Gene sequences are largely conserved – so if an ORF sequence is present 
            in multiple genomes, it likely represents a gene.
            '''
            )  

    


uploaded_file = st.file_uploader("Choose a file")
if uploaded_file is not None:
    # # To read file as bytes:
    # bytes_data = uploaded_file.getvalue()
    # st.write(bytes_data)

    # To convert to a string based IO:
    stringio = StringIO(uploaded_file.getvalue().decode("utf-8"))
    # st.write(stringio)

    # To read file as string:
    string_data = stringio.read()
    # st.write(string_data)

    # # Can be used wherever a "file-like" object is accepted:
    # dataframe = pd.read_csv(uploaded_file)
    # st.write(dataframe)    
    


if uploaded_file is not None:    
    sequence = string_data

   
    st.text_area(label = 'Generated Random Sequence', value = sequence)
# st.write(sequence2[0])
# st.write(sequence2.find('>'))
# st.write(sequence2.find(':'))
# st.write(sequence2.find('\r'))
    st.write(len(sequence))
    index = sequence.find('\r')
    sequence1 = sequence[index:] #returns the chars after the seen char or substring
    sequence1 = sequence1.lstrip().rstrip()
    sequence1 = sequence1.replace('\r\n', '')
    st.write(len(sequence1))
    st.write(len(sequence) - len(sequence1))
    st.markdown('#### Find the ORFs')

# st.write(f'The number of Open Reading Frames (ORFs) in this sequence is {sequence.count('ATG')}')

    ORF_count = sequence1.count('ATG')
    met1 = sequence1.find('ATG')

    if ORF_count == 1: 
    
        st.write(f'There is only one ORF in this sequence; it is a position {met1 + 1}.  Using this reading frame, the codons are as follows:')
    
        ORF = sequence1[met1:]
    
        codons = list(chunkstring(ORF, 3))
    
        codon_string = remove_short_codon(codons)
    
        st.write(codon_string)
    
        RNA_string = codon_string.replace('T', 'U')
    
        st.markdown('#### Transcribe the DNA into mRNA')

        st.write(RNA_string)
   
        aa = translate(codons)
    
        st.markdown('#### Translate mRNA into protein sequence')  

    
        col1, col2 = st.columns([2,4], gap='medium')
        with col1:
        
            st.image('./images/Aminoacids_table.png')

        with col2:
            protein_length = len(aa)
            # if 'STOP' in aa:
                #     protein_length = protein_length -1
                # else:
                    #     protein_length = protein_length
            
            st.write(f'The candidate protein is {protein_length} amino acids long.')
            aa_string = '.'.join(aa)
            st.write(aa_string)
        
            if protein_length <= 10:
                st.markdown('# :red[This is too short to be a protein candidate!]')
            else:
                st.markdown(('# :green[This looks interesting!  Requires further investigation]'))
        
    elif ORF_count > 1:
    
        st.write(f'There are {ORF_count} ORFs in this sequence; the first is at position {met1 + 1}.  Using this reading frame, the codons are as follows:')
    
        ORF = sequence1[met1:]
        # st.write(ORF)
        
        codons = list(chunkstring(ORF, 3))
    
        codon_string = remove_short_codon(codons)
    
        st.write(codon_string)
    
        RNA_string = codon_string.replace('T', 'U')
    
        st.markdown('#### Transcribe the DNA into mRNA')

        st.write(RNA_string)
   
        aa = translate(codons)
    
        st.markdown('#### Translate mRNA into protein sequence')  

    
        col1, col2 = st.columns([2,4], gap='medium')
        with col1:
        
            st.image('./images/Aminoacids_table.png')

        with col2:
            
            protein_length = len(aa)
            # st.write(type(protein_length))
            # if 'STOP' in aa:
                #     protein_length = protein_length -1
                #     st.write(protein_length)
                # else:
                    #     protein_length = protein_length
        
            st.write(f'The candidate protein is {protein_length} amino acids long.')
            aa_string = '.'.join(aa)
            st.write(aa_string)
        
            if protein_length <= 10:
                st.markdown('# :red[This is too short to be a protein candidate!]')
            else:
                st.markdown(('# :green[This looks interesting!  Requires further investigation]'))


    else:

        st.write('The are no ORFs in this sequence.  Try generating a longer sequence')
    
    
