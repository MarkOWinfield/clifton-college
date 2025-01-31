# -*- coding: utf-8 -*-
"""
Created on Tue Dec  3 13:06:03 2024

@author: Mark Winfield
"""

import streamlit as st
from random import choice
import pandas as pd
import plotly.express as px
import matplotlib.pyplot as plt
import numpy as np

st.set_page_config(layout="wide")

st.markdown(""" <style> .font1 {
    font-size: 45px; font-family: 'Copper Black'; color: #FF9633}
    </style> """, unsafe_allow_html=True)

###############################################################################
################################ FUNCTIONS ####################################
###############################################################################

def terminator(unknown_seq):
    
    seq_len = len(unknown_seq)

    adenine = ['A', 'A', 'A', 'A', 'A', 'A', 'A', 'A','a']
    cytosine  = ['C', 'C', 'C', 'C', 'C', 'C', 'C', 'C','c']
    guanine = ['G', 'G', 'G', 'G', 'G', 'G', 'G', 'G','g']
    thymine = ['T', 'T', 'T', 'T', 'T', 'T', 'T', 'T','t']
        
    sequence_list = []
        
    for k in range(2000):
        sequence = ''
        for i in range(seq_len):
     
                
            if unknown_seq[i] == 'A':
                base=choice(thymine)
                sequence += base
                if base == 't':
                    break
                else:
                    continue            

            if unknown_seq[i] == 'C':
                base=choice(guanine)
                sequence += base
                if base == 'g':
                    break
                else:
                    continue            

            if unknown_seq[i] == 'G':
                 base=choice(cytosine)
                 sequence += base
                 if base == 'c':
                     break
                 else:
                     continue           

            if unknown_seq[i] == 'T':
                base=choice(adenine)
                sequence += base
                if base == 'a':
                    break
                else:
                    continue
               
            sequence += unknown_seq[i]

        sequence_list.append(sequence)
        sequence_list.sort(key=len)
        
        # Convert list to a dictionary and back in order to eliminate duplicate values
        sequence_list_sorted = list(dict.fromkeys(sequence_list))
        
        if unknown_seq in sequence_list:
            sequence_list.remove(unknown_seq)
            sequence_list_sorted.remove(unknown_seq)
        else:
            continue
        
        # for i in sequence_list_sorted:
        #     freq = sequence_list.count(i)
#            st.write(f'{i} = {freq}')
        
    return(sequence_list_sorted)


###############################################################################
###############################################################################

st.title(' :blue[Sanger Sequencing Simulation]')

col1, col2 = st.columns([2,1], gap='large')

with col1:

    st.markdown(

'''
## Introduction

Developed by Fred Sanger (hence the name) in 1975, Sanger sequencing was the first
method of DNA sequencing. It was the method used for the ground-breaking Human 
Genome Project, completed in 2003. Other names for this technique are 
‘chain-termination sequencing’ and ‘dideoxy sequencing’.

## How does Sanger sequencing work?

DNA is used as a template in a polymerase chain reaction (PCR).  In four
separate Eppendorf tube one places a mix of normal nucleotides (dNTPs) and 
chain-terminating nucleotides (ddNTPs) is used in the PCR reaction.  The ddNTPs
are labelled with either a radioisotope of a fluorochrome.  Four 
reactions are set up in four different Eppendorf tubes: one with ddAdenine, 
one with ddCytosine one with ddGuanine and one with ddThymine.
When a chain-terminating base is randomly incorporated into a growing DNA 
chain, it cannot grow any further. This means that DNA fragments of different 
lengths are generated. Each fragment ends in a chain-terminating base.
The DNA fragments are then separated by size using electrophoresis - this 
could be through a gel or through a capillary tube.  If a gel is used then the
dideoxynucloetides will be radiactively labelled and, once the products of the
PCR reaction have been run through the gel, will be visualised by exposure of the 
gel to a film.  If fluorochomes are used, each of the four chain-terminating bases (A/T/C/G) 
will have a different fluorescent label. As the products of the PCR reaction
pass along the capillary tube a laser is used to excite the fluorescently labelled bases at the end of each fragment.
Shorter fragments come first in the sequence followed by increasingly longer fragments.
The fluorescence of the base that terminated each length of fragment is recorded, and a chromatograph is generated showing which base is present at which position along the DNA fragment.
The chromatogram is compared with a reference file to identify any variants.

### Advantages and limitations of Sanger sequencing
#### Advantages

- Gold standard method for accurate detection of single nucleotide variants and small insertions/deletions.
- More flexible for testing for a specific familial variant than massively parallel sequencing.
- Cost effective where single samples need to be tested very urgently, so they cannot be batched up (for example, in prenatal testing or parental carrier testing during a pregnancy).
- Less reliant on computational tools than massively parallel sequencing.
- In some cases, longer fragments (up to approximately 1,000 base pairs) can be sequenced than in massively parallel sequencing.

#### Limitations

- Limited throughput.
- Not cost effective for sequencing many genes in parallel, or for sequencing 
the same region in many samples.
- May not detect mosaicism.
- Can require a larger amount of input DNA than massively parallel sequencing.


'''
)


with col2: 
    
    st.image('./images/deoxy.jpg')
    st.markdown('''
                - At carbon 2 there is an hydrogen (no oxygen}
                - At carbon 3 there is a hydroxyl group (OH)
                - That is, the sugar is 2-deoxyribose
                
                Deoxynucleotides can form bonds with the phosphate group of
                other deoxynucleotides so the DNA fragment can continue to grow.
                ''')
    
    st.image('./images/dideoxy.jpg')
    st.markdown('''
                - At carbon 2 there is a hydrogen (no oxygen}
                - At carbon 3 there is a hydrogen (no oxygen)
                - That is, the sugar is 2,3-dideoxyribose
                
                Dideoxynucleotides can't form bonds so the DNA fragement can't be extended once
                this is incorporated - that is, it TERMINATES elongation!
                
                ''')
    


st.header('Begin the Experiment')

st.image('./images/Eppendorfs_1.png', caption='Load the Eppendorfs', width = 750)


st.header('Add the DNA to the mix')

unknown_seq = st.text_input(
        "Enter a sequence of about 30 bases (A, C, G, T) 👇", ''

        )
unknown_seq = unknown_seq.upper()

st.header('Run the PCR reaction')

st.image('./images/Eppendorfs_2.png', width = 1000)


sequence_list = terminator(unknown_seq)

a_list = []
c_list = []
g_list = []
t_list = []

fragment_vis = []

df = pd.DataFrame(columns=['nucleotide', 'fragment length', 'track']) 
df_fragments = pd.DataFrame(columns=['Fragment', 'Length (bases)'])

    
for fragment in sequence_list: 
#    print(f'{fragment} = {len(fragment)} bases long')
    fragment_vis.append(fragment)
    df_fragments.loc[len(df_fragments)] = {'Fragment': fragment, 'Length (bases)': len(fragment)}
    if fragment[-1] == 'a':
        a_list.append(len(fragment))
        new_value = {'nucleotide': 'A', 'fragment length': len(fragment), 'track': 1}
        df.loc[len(df)] = new_value
    if fragment[-1] == 'c':
        c_list.append(len(fragment))
        new_value = {'nucleotide': 'C', 'fragment length': len(fragment), 'track': 2}
        df.loc[len(df)] = new_value
    if fragment[-1] == 'g':
        g_list.append(len(fragment))
        new_value = {'nucleotide': 'G', 'fragment length': len(fragment), 'track': 3}
        df.loc[len(df)] = new_value
    if fragment[-1] == 't':
        t_list.append(len(fragment))
        new_value = {'nucleotide': 'T', 'fragment length': len(fragment), 'track': 4}
        df.loc[len(df)] = new_value

st.header('Polyacrylamide Gel Electrophoresis')
        
 
col1, col2, col3 = st.columns([1,2,1], gap='large')

with col1:
    
    st.markdown("""
                Because of the inclusion of radio-labelled dideoxynucleotide 
                (represented by lower case letters) in the PCR mix, extending
                fragments will be terminated.
                
                In order to'read' the sequence, start from the bottom of the gel 
                (remember, the smaller the fragment the further if migrates through
                 the gel) and record the first band (from left to right the nucleotides A, C G T)
                and then move up the gel recording each band in turn.  You now 
                have your sequence!  However, be aware that there amy be reading
                errors: a nucleotide not recorded. 
        
                """
                )

with col2:

    fig = px.scatter(df, x="track", y="fragment length", width = 500, height = 850)

    fig.update_traces(
        marker=dict(size=50, symbol="line-ew", line=dict(width=3, color="DarkSlateGrey")),
        selector=dict(mode="markers")
        )
    
    fig.update_layout(
        xaxis = dict(
            tickmode = 'array',
            tickvals = [1, 2, 3, 4],
            ticktext = ['A', 'C', 'G', 'T'], color = 'blue'
        )
    )

    fig.update_layout(yaxis_title=None)
    fig.update_layout(xaxis_title=None)
    fig.update_xaxes(showgrid=False)
    fig.update_yaxes(showgrid=False)

    fig.update_xaxes(title_font=dict(size=18, family='Courier', color='crimson'))
    fig.update_yaxes(showticklabels=False)
    fig.update_layout(
        xaxis = dict(
        tickfont = dict(size=40)))
    fig.update_xaxes(title_font_color="blue")
  
    fig.update_layout({
       'plot_bgcolor': 'rgba(189, 195, 199, 1)',
       'paper_bgcolor': 'rgba(0, 0, 0, 0)',
        })

    st.plotly_chart(fig)

with col3:

    st.write("Amplified fragments produced in the test tubes.")

    # Reverse row order and reset index

    df_fragments = df_fragments.loc[::-1].reset_index(drop = True)

    st.table(df_fragments.Fragment)
