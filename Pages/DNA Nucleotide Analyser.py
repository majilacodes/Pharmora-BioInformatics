import streamlit as st
import pandas as pd
import altair as alt

st.sidebar.title('Input DNA Sequence')
st.sidebar.write('DNA Sequence Analyser, a comprehensive tool to analyse DNA Nucleotides of a given DNA sequence.')
st.write("")

default = "GAACACGTGGAGGCAAACAGGAAGGTGAAGAAGAACTTATCCTATCAGGACGGAAGGTCCTGTGCTCGGGATCTTCCAGACGTCGCGACTCTAAATTGCCCCCTCTGAGGTCAAGGAACACAAGATGGTTTTGGAAATGCTGAACCCGATACATTATAACATCACCAGCATCGTGCCTGAAGCCATGCCTGCTGCCACCATGCCAGTCCT"

sequence = st.sidebar.text_area('User Input', default, height=275)

st.write("""
# DNA Nucleotide Analyser
""")

def DNA_nucleotide_count(seq):
    return dict([
        ('A', seq.count('A')),
        ('T', seq.count('T')),
        ('G', seq.count('G')),
        ('C', seq.count('C'))
    ])

def update():
    cleaned_sequence = ''.join(sequence.split()).upper()

    # 1. Print dictionary
    st.subheader('Dictionary')
    X = DNA_nucleotide_count(cleaned_sequence)
    st.write(X)

    # 2. Print text
    st.subheader('Text output')
    st.write('There are  ' + str(X['A']) + ' adenine (A)')
    st.write('There are  ' + str(X['T']) + ' thymine (T)')
    st.write('There are  ' + str(X['G']) + ' guanine (G)')
    st.write('There are  ' + str(X['C']) + ' cytosine (C)')

    # 3. Display DataFrame
    st.subheader('DataFrame')
    df = pd.DataFrame.from_dict(X, orient='index', columns=['count'])
    df.reset_index(inplace=True)
    df = df.rename(columns={'index': 'nucleotide'})
    st.write(df)

    # 4. Display Bar Chart using Altair
    st.subheader('Bar chart')
    p = alt.Chart(df).mark_bar().encode(
        x='nucleotide',
        y='count'
    )
    p = p.properties(
        width=alt.Step(80)  # controls width of bar.
    )
    st.write(p)

if st.sidebar.button('Analyze'):
    update()
else:
    st.write("")
    st.warning('👈 Enter DNA sequence data!')