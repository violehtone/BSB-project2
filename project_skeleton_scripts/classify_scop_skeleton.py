#!/usr/bin/python

import itertools
import argparse
import csv
from itertools import permutations

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
                ## Get the PDB (dic key)
                pdb = row[1]

                ## generate the scop hierarchy dic (dic value)
                scop_dic = {}
                scop_hierarchy = [x.strip() for x in row[5].split(',')] # creates a list of SCOP hierarchy data
                scop_class = scop_hierarchy[0].split("=", 1)[1]
                scop_fold = scop_hierarchy[1].split("=", 1)[1]
                scop_superfamily = scop_hierarchy[2].split("=", 1)[1]
                scop_family = scop_hierarchy[3].split("=", 1)[1]
                scop_dm = scop_hierarchy[4].split("=", 1)[1]
                scop_sp = scop_hierarchy[5].split("=", 1)[1]
                scop_px = scop_hierarchy[6].split("=", 1)[1]

                scop_dic["class"] = scop_class
                scop_dic["fold"] = scop_fold
                scop_dic["superfamily"] = scop_superfamily
                scop_dic["family"] = scop_family
                scop_dic["domain"] = scop_dm
                scop_dic["species"] = scop_sp # similar
                scop_dic["protein"] = scop_px

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
    if (prot1_scop["superfamily"] == prot2_scop["superfamily"]):
        return "similar"
    elif (prot1_scop["family"] == prot2_scop["family"]):
        return "ambiguous"
    else:
        return "different"

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
    pairs = []
    ##########################
    ### START CODING HERE ####
    ##########################
    # You can add a pair of proteins to the list using the following code:
    # pairs.append((protein1, protein2))
    # Generate all possible combinations of IDs

    # Get all permutations of length 2 from protein_ids
    perm = permutations(protein_ids, 2)
    permList = list(perm)

    temp = []
    for pair in permList:
        pairs.append(pair)

    for a, b in pairs:
        if (a,b) not in temp and (b,a) not in temp:
            temp.append((a,b))
    
    pairs = temp * 1 # copy temp to permList

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
    scop_homology = {} ## {{:UniprotID, :similarity(i.e. similar}}

    ##########################
    ### START CODING HERE ####
    ##########################
    # You should remember to take care about the proteins that are not in the SCOP database.
    for pair in pairs:
        protein1_id = pair[0] ## uniprot id
        protein2_id = pair[1] ## uniprot id
        
        protein1_pdb = protein_ids_pdbs[protein1_id].lower() #pdb
        protein2_pdb = protein_ids_pdbs[protein2_id].lower() #pdb

        try:
            # find scop-data for proteins:
            protein1_scop_data = scop_dict[protein1_pdb]
            protein2_scop_data = scop_dict[protein2_pdb]

            # find similarity for the protein pair
            similarity = check_similarity_for_protein_pair(protein1_scop_data, protein2_scop_data)

            # add uniprot id : similarity -pair into scop_homology dict
            scop_homology[(protein1_id, protein2_id)] = similarity
        except KeyError:
            print("Protein pair not found from SCOP data [%s, %s]", protein1_pdb, protein2_pdb)

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

    return uniprot_id_list

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

    ## Uniprot id i.e. P16099 / Q92547..
    ## PDB id i.e. 2D5Z, 1FSL, 1LH3..
    ## SCOP-ID i.e. 15225, 15230..

    # fetch a dict of {Uniprot_id : PDB}
    uniprot_PDB_mapping = read_lookup_table(pdb_id_file)

    # fetch a list of Uniprot ids (id, id, id)
    uniprot_ids = read_protein_ids_file(uniprot_filename)

    # retrieve scop data in dictionary {pdb: {family : "", ... "species" : ""}}
    pdb_scop_data_dic = retrieve_scop_data(scop_file)

    # Get all possible permutations of protein pairs
    all_possible_pairs = generate_all_possible_protein_pairs(uniprot_ids)

    # calculate scop homology
    scop_homology = assign_homology(pdb_scop_data_dic, uniprot_PDB_mapping, all_possible_pairs)

    # write the results
    write_results(output_filename, scop_homology)

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