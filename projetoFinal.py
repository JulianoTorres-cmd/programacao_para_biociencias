#!/usr/bin/env python3
#sempre fazer isso para executar o programa no terminal: "chmod +x programa.py"
#para ver a versao de um import ou se ele realmente está no pc: "print(SeuImport.__version__)"

import sys
import pandas as pd
from Bio import SeqIO
from Bio.Blast.Applications import NcbiblastxCommandline

# c) cpm = num de leitura mapeadas * 10^6 / total de mapeados = 15115 para cada celula
# d) (rep1_a_cpm + rep2_a_cpm / 2) para cada celula
# e) 5 maiores de cada letra (A e B [cond_x_cpm_media])
# f) procurar uma pesquisa blast dos 10 genes selecionados (5 A e 5 B)
# g) imprimir o melhor hit dos 10 genes com o maior bitscore.
# h) deve estar no formato:	gene_id     Cond_A_CPM_media     Cond_B_CPM_media      id_proteína_encontrada 

def createCPM(dataframe):

    # lista com os valores a serem adicionados
    rep1_A_CPM = []
    rep2_A_CPM = []
    rep1_B_CPM = []
    rep2_B_CPM = []

    # para cada celula do data frame, normalizo e salvo em sua respectiva lista. no final, eu insiro o valor
    for i in range(dataframe.index.stop):
        value1_A = int( (dataframe['Rep1_A'][i] ) * 10 ** 6) / float(dataframe.index.stop)
        rep1_A_CPM.append(value1_A)

        value2_A = int( ( dataframe['Rep2_A'][i] ) * 10 ** 6  ) / float(dataframe.index.stop)
        rep2_A_CPM.append(value2_A)

        value1_B = int( ( dataframe['Rep1_B'][i] ) * 10 ** 6  ) / float(dataframe.index.stop)
        rep1_B_CPM.append(value1_B)

        value2_B = int( ( dataframe['Rep2_B'][i] ) * 10 ** 6  ) / float(dataframe.index.stop)
        rep2_B_CPM.append(value2_B)
    
    dataframe.insert(loc = 5, column = 'Rep1_A_CPM', value = rep1_A_CPM)
    dataframe.insert(loc = 6, column = 'Rep2_A_CPM', value = rep2_A_CPM)
    dataframe.insert(loc = 7, column = 'Rep1_B_CPM', value = rep1_B_CPM)
    dataframe.insert(loc = 8, column = 'Rep2_B_CPM', value = rep2_B_CPM)

    return dataframe

def createCPM_media(dataframe):

    cond_A_CPM_media = []
    cond_B_CPM_media = []

    for i in range(dataframe.index.stop):
        value_A_media = ( int( dataframe['Rep1_A_CPM'][i] ) + int( dataframe['Rep2_A_CPM'][i] ) ) / float(2)
        cond_A_CPM_media.append(value_A_media)

        value_B_media = ( int( dataframe['Rep1_B_CPM'][i] ) + int( dataframe['Rep2_B_CPM'][i] ) ) / float(2)
        cond_B_CPM_media.append(value_B_media)
    
    dataframe.insert(loc = 9, column = 'Cond_A_CPM_media', value = cond_A_CPM_media)
    dataframe.insert(loc = 10, column = 'Cond_B_CPM_media', value = cond_B_CPM_media)

def expressos_A_e_B(dataframe):
    # retorno uma lista com os 10 maiores valores encontrados. 5 primeiro referentes a A_cpm_medio e os ultimos 5 ao B_cpm_medio

    maioresValoresA_B = []
    exp_a = dataframe.nlargest(5, 'Cond_A_CPM_media')['gene_id']
    exp_b = dataframe.nlargest(5, 'Cond_B_CPM_media')['gene_id']
    for i in exp_a:
        maioresValoresA_B.append(i)
    for i in exp_b:
        maioresValoresA_B.append(i)
    return maioresValoresA_B

# lendo por linha de comando
tabela = sys.argv[1]
desconhecido = sys.argv[2]
prolixus = sys.argv[3]

# tabela excel
df = pd.read_excel(tabela)

# fluxo do programa
createCPM(df)
createCPM_media(df)
maioresValores = expressos_A_e_B(df)
# o programa gera a tabela com todas as coisas que o senhor pediu no codigo abaixo (tera print com o resultado na tabela): 
#df.to_excel('tabelaTESTE.xlsx')

# repartindo o fasta "Rdesconhecidus.fasta" e armazenando as sequencias mais expressas, OBS: só uso isso aqui para confirmar quais eram as sequencias mais expressas 
seq_maiores_valores = []
fasta = SeqIO.parse(open(desconhecido, 'r'), 'fasta')
#fasta = SeqIO.parse(open('Rdesconhecidus.fasta', 'r'), 'fasta')
for line in fasta:
    if(line.id in maioresValores):
        #seq_maiores_valores.append(line.id)
        seq_maiores_valores.append(line.seq)

# esses print confirmam quais são os 10 genes mais expressos, OBS: os genes "4174" e "11244" estao duplicados, ou seja, sao os mais espressos tanto no Cond_A_CPM_media quanto em Cond_B_CPM_media
#print(maioresValores)
#print(seq_maiores_valores)
        
 # abaixo segue o blast individual que fiz para cada gene (tera um print com os resultados), peguei a sequencia de cada "geneid" no arquivo "Rdesconhecidus.fasta" segundo a minha lista "maioresValores"       
blastx = "/home/juliano/anaconda3/bin/blastx"
blast = r"/home/juliano/Documentos/projetoFinalprog2/blast_projeto_final.txt"

meu_blast = NcbiblastxCommandline(cmd = blastx ,query = desconhecido, subject = prolixus, evalue = 0.05, outfmt = 6, out = blast)
stdout, stdeer = meu_blast()
result = pd.read_csv("/home/juliano/Documentos/projetoFinalprog2/blast_projeto_final.txt", sep='\t', names=["qseqid","sseqid","pident","length","mismatch","gapopen","qstart","qend","sstart","send","evalue","bitscore"])
max_score = result.sort_values('bitscore')
print(max_score.iloc[[-1]])








