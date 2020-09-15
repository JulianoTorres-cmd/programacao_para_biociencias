import csv

#TAC1 PROG2: Juliano Torres
#resposta da questão 1:

#refArquivoEntrada1 = open("/home/juliano/Downloads/TcCLB.506717.80mRNA-p1.fasta")
#refArquivoSaida = open("/home/juliano/Downloads/TcCLB.506717.80mRNA-p1.tsv", "w")
#cabecalho = refArquivoEntrada1.readline()[1:-346]
#sequencia = ""
#for linha in refArquivoEntrada1:
    #sequencia += linha.replace('\n', '')
#refArquivoSaida.write("{}\t{}".format(cabecalho, sequencia))
#refArquivoEntrada1.close()
#refArquivoSaida.close()

#respota da questão 2:
#muitas dúvidas nessa questão!

refArquivo = open("/home/juliano/Downloads/Tcruzi_AnnotatedProteins.fasta")
cabecalho = refArquivo.readline()[1:-1]
sequencia = ""
for linha in refArquivo:
    sequencia += linha.replace('\n','')
print("Cabecalho: {}".format(cabecalho))
print("Sequencia: {}".format(sequencia))
refArquivo.close()



#respota da questão 3:

#with open("/home/juliano/Downloads/species.csv") as species:
    #csv = csv.reader(species, delimiter=',')
    #for row in csv:
        #print(row)




