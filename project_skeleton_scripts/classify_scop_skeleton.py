#!/usr/bin/python

import itertools
import argparse
import csv

def retrieve_scop_data(scop_file):
    """
    Reads a database file from SCOP and returns a dictionary with the protein IDS mapping.
    :param scop_file: database file containing mapping of PDB's to SCOP ID's.
    :return: dictionary containing the mapping of PDBs (keys) and their corresponding SCOP data (values).
    """
    scop_data = {}

    ##########################
    ### START CODING HERE ####
    ##########################
    # You can parse SCOP data in various ways. E.g. you can use dictionary of dictionaries
    # {proteinID: {"class": class, "fold": fold, "superfamily": superfamily, 'family': family}}
    with open(scop_file) as sf:
        reader = csv.reader(sf, delimiter="\t")
        for row in reader:
            if(row[0][0] != "#"):
                ## Get the PDB (key)
                pdb = row[1]
                
                ## generate the scop hierarchy dic (value)
                scop_dic = {}
                scop_hierarchy = [x.strip() for x in row[5].split(',')]
                scop_class = scop_hierarchy[1].split("=", 1)[1]
                scop_fold = scop_hierarchy[1].split("=", 1)[1]
                scop_superfamily = scop_hierarchy[1].split("=", 1)[1]
                scop_family = scop_hierarchy[1].split("=", 1)[1]

                scop_dic["class:"] = scop_class
                scop_dic["fold:"] = scop_fold
                scop_dic["superfamily:"] = scop_superfamily
                scop_dic["family:"] = scop_family

                ## add an element to the scop_data dictionary
                scop_data[pdb] = scop_dic
    
    ########################
    ### END CODING HERE ####
    ########################

    return scop_data


def compute_similarity_score(prot1_scop, prot2_scop):
    """
    Computes the score for two proteins on the basis of the data from SCOP database.
    :param prot1_scop: data for protein 1 from SCOP database.
    :param prot2_scop: data for protein 2 from SCOP database.
    :return: similarity score (float)
    """
    ##########################
    ### START CODING HERE ####
    ##########################
    # You need to decide whether you need this function for SCOP database.
    pass


    ########################
    ### END CODING HERE ####
    ########################


def check_similarity_for_protein_pair(prot1_scop, prot2_scop):
    """
    Returns the similarity score between two proteins.
    :param prot1_scop: data for protein 1 from SCOP database.
    :param prot2_scop: data for protein 2 from SCOP database.
    :param pair: a tuple with the UniProt IDs of the two proteins to compare.
    :return: "different", "similar" or "ambiguous".
    """
    ##########################
    ### START CODING HERE ####
    ##########################
    pass
    
    
    ########################
    ### END CODING HERE ####
    ########################

# If you will use the numeric score for SCOP (similar to GO), you may want to use check_similarity_for_protein_pair
# with other arguments. See the example below.
# def check_similarity_for_protein_pair(score, threshold):
#    pass

def generate_all_possible_protein_pairs(protein_ids):
    """
    Returns a list containing all unique protein pairs.
    :param protein_ids: list of all proteins IDs.
    :return: list of possible unique protein pairs.
    """
    pairs = list()
    ##########################
    ### START CODING HERE ####
    ##########################
    # You can add a pair of proteins to the list using the following code:
    # pairs.append((protein1, protein2))
    # Generate all possible combinations of IDs



    ########################
    ### END CODING HERE ####
    ########################
    return pairs


def assign_homology(scop_dict, protein_ids_pdbs, pairs):
    """
    Computes the similarity score between all protein pairs from the list, and decides if two proteins are homologs
    (different, ambiguous or similar).
    :param scop_dict: dictionary containing the mapping of PDBs (keys) and their corresponding SCOP data (values).
    :param protein_ids_pdbs: dictionary with UniprotID as key and PDB ID as a value.
    :param pairs: list of all possible unique protein pairs.
    :return: dictionary with UniProt ID (key), similarity(different, ambiguous or similar).
    """
    scop_homology = {}

    ##########################
    ### START CODING HERE ####
    ##########################
    # You should remember to take care about the proteins that are not in the SCOP database.



    ########################
    ### END CODING HERE ####
    ########################
    return scop_homology


def write_results(filename, scop_homology):
    """
    Writes in an output file the all of the protein pairs and their similarity/dissimilarity.
    :param output_filename: the name of the output file.
    :param scop_homology: dictionary (keys: protein pairs as tuples; values: one of the value - different/similar/ambiguous)
    """
    with open(filename, "w") as f:
        for (p1, p2), value in scop_homology.items():
            f.write("\t".join([p1, p2, value]) + "\n")


def read_protein_ids_file(filename):
    """
    Returns the list of UniProt IDs contained in the input file.
    :param filename:
    :return: list of UniProt IDs in the input file.
    """
    ##########################
    ### START CODING HERE ####
    ##########################
    uniprot_id_list = []

    uniprot_ids = open(filename)
    for line in uniprot_ids:
        protein_id = line.strip()
        uniprot_id_list.append(protein_id)
    
    uniprot_ids.close()
    return uniprot_ids

    #######################
    ### END CODING HERE ###
    #######################



def read_lookup_table(filename):
    """
    Reads the specified file and returns the dictionary with UniprotID as key and PDB ID as a value.
    :param filename: file with the mapping between Uniprot ids and PDB ids.
    :return: dictionary with UniprotID as key and PDB ID as a value.
    """
    ##########################
    ### START CODING HERE ####
    ##########################
    result_dict = {}

    with open(filename) as csv_file:
        reader = csv.reader(csv_file, delimiter="\t")
        for row in reader:
            uniprotId = row[0]
            PDB_ID = row[1]
            result_dict[uniprotId] = PDB_ID
    
    return result_dict
    #######################
    ### END CODING HERE ###
    #######################



def main(uniprot_filename, output_filename, pdb_id_file, scop_file):
    ##########################
    ### START CODING HERE ####
    ##########################
    pass

    uniprot_PDB_mapping = read_lookup_table(pdb_id_file)
    retrieve_scop_data(uniprot_PDB_mapping)


    #######################
    ### END CODING HERE ###
    #######################


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='The script retrieves data from SCOP database (from local file)'
                                                 ' and provides an output file'
                                                 ' with the strings "<id1>   <id2>   <similarity>", where'
                                                 ' <similarity> is a string with one of the values from '
                                                 ' different/ambiguous/similar.')

    parser.add_argument("-ids", "--uniprot_filename", help="File with the protein Uniprot IDs", required=True)
    parser.add_argument("-pdb", "--pdb_id_file", help="File with the mapping between Uniprot ids and PDB ids", required=True)
    parser.add_argument("-s", "--scop_file", help="SCOP database file", required=True)
    parser.add_argument("-o", "--output_filename", help="Output file name", required=True)

    args = parser.parse_args()

    uniprot_filename = args.uniprot_filename
    output_filename = args.output_filename
    pdb_id_file = args.pdb_id_file
    scop_file = args.scop_file

    main(uniprot_filename, output_filename, pdb_id_file, scop_file)