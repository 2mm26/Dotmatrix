'''
23261089
毛利仁香
'''

import streamlit as st

from Bio import SeqIO
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
matplotlib.rcParams["font.family"] = "BIZ UDGothic"


def dotmatrix(f1,f2,win):
    record1 = next( SeqIO.parse(f1,"fasta"))
    record2 = next( SeqIO.parse(f2,"fasta"))
    seq1 = record1.seq
    seq2 = record2.seq

    #win = 10

    a = len(seq1)-win+1
    b = len(seq2)-win+1

    width = 500
    height = 500

    image = np.zeros( (height,width) ) 

    hash = {}

    for x in range(a):
        subseq1 = seq1[x:x+win]
        if subseq1 not in hash:
            hash[subseq1] = []
        hash[subseq1].append(x)

    for y in range(b):
        subseq2 = seq2[y:y+win]
        py = int(y/b*height)
        if subseq2 in hash:
            for x in hash[subseq2]:
                px = int(x/a*width)
                image[py,px] = 1

    plt.imshow(image,extent=(1,a,b,1),cmap="Grays")
    #plt.show()

    st.pyplot(plt)


st.title("Dot matrix")

file1 = st.sidebar.file_uploader("Sequence file 1:")
file2 = st.sidebar.file_uploader("Sequence file 2:")

win = st.sidebar.slider("Window size:",4,100,10)

from io import StringIO

if file1 and file2:
    with StringIO(file1.getvalue().decode("utf-8")) as f1,\
         StringIO(file2.getvalue().decode("utf-8")) as f2:
        dotmatrix(f1,f2,win)