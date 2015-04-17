
# coding: utf-8

# In[ ]:

#####Part 1 - BLAST eukaryotic queries against a local human and bacterial database
#Module imported
from Bio.Blast.Applications import NcbiblastpCommandline #Needed to run the BLAST search using the queries as the input
#Directories
name = "Xenoturbella_24JAN_preclast.fa"
data_base = "bacterial_and_humanDB.faa"
out_file = "Xenoturbellavs2and9606.xml"
#BLAST search - 32 threads, max of 3 hits per query
comline = NcbiblastpCommandline(query=name, db=data_base, evalue= 1e-5, max_target_seqs=3, outfmt=5, out=out_file, num_threads=32)
stdout, stderr = comline()
#####Part 2 - Parse BLAST Results and Queries and create empty lists to append queries to based on their response to criteria
#Modules imported
from Bio.Blast import NCBIXML #Needed to parse the BLAST output into a usable format
from Bio import SeqIO #Needed to parse the queries
#Empty Lists to append query sequences too based on how their BLAST results respond to contaminated criteria
contaminated_sequences = []
uncontaminated_sequences = []
check_sequences = []
#Query record lists and dictionaries
query_records = list(SeqIO.parse(name, "fasta")) #list created to allow iteration through query records
query_record_dict = SeqIO.to_dict(query_records)
#Blast record lists and dictionaries
result_handle = open(out_file)
blast_records = list(NCBIXML.parse(result_handle))
query_record_dict = SeqIO.to_dict(query_records)
#####Part 3 - Filter out Contaminated Queries using a For-Loop to iterate through parsed BLAST Results and Queries
for blast_record in blast_records: #for each blast record in the BLAST results
    split_query_title = list((blast_record.query).split()) # making a list out of the words in the query title for the blast record
    if len(blast_record.alignments) == 0: #if there are no results given
        uncontaminated_sequences.append(query_record_dict[split_query_title[0]]) #appends the query record from the dictionary with the key given from the first word of the query title for this blast record
        continue # continue iterating through the for-loop
    concluded = False # at this point in the for loop the contamination identity is not yet concluded
    for alignment in blast_record.alignments: 
        if concluded: #if this query has been appended to one of the lists
            break #skip to the next blast record
        for hsp in alignment.hsps: #for each high scoring pair in this blast record
            if concluded: #if this query has been appended to one of the lists
                break #skip to the next blast record
            A=hsp.identities*100.0/(len(hsp.sbjct)) >= 50.0 #identity>=50
            B=hsp.expect == 0.0 #evalue=0
            C=hsp.query == hsp.sbjct #'draft' query sequence=hit sequence. Therefore it is not a draft sequence
            D=len(blast_record.alignments) > 1 #there are more alignments
            if (A) and ((not (B)) or (not (C))) or (A) and (not((B) and (C))): #if AnBn!C or An!B
                ''' If % identity>=50, evalue=0 but query sequence different from hit sequence:
                or % identity>=50 and evalue is not 0:'''
                contaminated_sequences.append(query_record_dict[split_query_title[0]]) #appends the query record from the dictionary with the key given from the first word of the query title for this blast record
                concluded = True #at this point in the for loop the contamination identity is concluded
            elif (not (A)) or (not(D)): #if A! or AnBnCn!D
                ''' If % identity<50:
                or % identity >=50, evalue=0, query sequence=hit sequence and there are no more alignments:'''
                uncontaminated_sequences.append(query_record_dict[split_query_title[0]]) #appends the query record from the dictionary with the key given from the first word of the query title for this blast record
                concluded = True #at this point in the for loop the contamination identity is concluded
        if concluded: #if this query has been appended to one of the lists
            break #skip to the next blast record
    if not concluded:  #if a query has not been appended to one of the lists after iterating through all hsps given, assume it is a bacterial or human sequence
        check_sequences.append(query_record_dict[split_query_title[0]])  #appends the query record from the dictionary with the key given from the first word of the query title for this blast record
#####Part 4 - Write lists created to three distinct new fasta files seperating queries according to contamination status
SeqIO.write(contaminated_sequences, "contaminated_xenoturbella_queries.fa", "fasta")
SeqIO.write(uncontaminated_sequences, "xenoturbella_queries.fa", "fasta")
SeqIO.write(check_sequences, "check_xenoturbella_queries.fa", "fasta")

