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
###############################################################################


st.title(' :blue[How to find Hypothetical Genes]')

st.markdown('''
           DNA sequencing has experienced a boom over the last couple of decades since
           the Human Genome Project was completed in the year 2000.  Since that 
           time, many species have had their genome sequenced, and research scientist
           are now faced with the problem of making sence of the sequences.  One of
           their tasks is to identify coding sequences (genes).  One approach to this
           problem is to search for Open Reading Frames (ORFs) - spans of DNA
           between a START (ATG) and a STOP (TAA, TAG, TGA) codon - and predict the
           hypothetical proteins that would be generated from them.
           
           Since DNA is interpreted in groups of three nucleotides (codons), a 
           DNA strand (not the double strand! )has three possible reading frames (see image below).  Given
           that DNA is a double helix made up of two strands each with three possible
           reading frame, any molecule of DNA has six possible reading frames.
           '''
           )

st.image('./images/reading-frames.jpeg')          
           
st.markdown('''
           The average bacterial protein, I believe, is in the range 350 - 400
           amino acids.  Given that bacterial genes do not have introns then the
           average gene coding for these proteins must be in the range 1050 -
           1200 nucleotides.  This simple program generates random DNA sequences
           that span this range and then performs the following analysis of the
           sequence:
           - Looks for all Open Reading Frames (ATG)
           - Highlights the codons reading from the first start codon (ATG)
           - Transcribes the DNA sequence to its mRNA equivalent (copied from the -ve strand)
           - Translates the mRNA into its hypothetical protein sequence
            
           '''
            )
    

    
st.markdown('''
            To identify an open reading frame:
            Locate a sequence corresponding to a start codon in order to determine the reading frame – this will be ATG (sense strand)
            Read this sequence in base triplets until a stop codon is reached (TGA, TAG or TAA)
            The longer the sequence, the more significant the likelihood that the sequence corresponds to an open reading frame


            Certain bioinformatic programs can automatically identify potential ORFs when provided with a candidate sequence

            Gene sequences are largely conserved – so if an ORF sequence is present in multiple genomes, it likely represents a gene
            '''
             )
