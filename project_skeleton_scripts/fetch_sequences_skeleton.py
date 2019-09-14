#!/usr/bin/python

# This script downloads all the sequences in the list provided by the user, and puts them into a single file to be used
# as a database.
# It also stores the fetched fasta sequences into individual files which will be used as queries for (PSI-)BLAST search.

import urllib.request
import argparse
import os.path

def fetch_one_fasta(uniprot_id):
    """
    Fetch the fasta formatted sequence of the uniProtID supplied in the argument.
    :param uniProtID: id of the protein to fetch in uniprot database.
    :return: sequence data in fasta format.
    """
    ##########################
    ### START CODING HERE ####
    ##########################
    # Define the variable 'url' as a string with the URL to the Fasta formatted sequence of the uniProtID
    # in the uniprot website (www.uniprot.org).
    url = "https://uniprot.org/uniprot/" + str(uniprot_id) + ".fasta"
    ##########################
    ###  END CODING HERE  ####
    ##########################

    fh = urllib.request.urlopen(url)
    result = fh.read().decode("utf8")
    fh.close()

    print(result)

    return result


def check_query_folder(query_folder):
    if not os.path.exists(query_folder):
        return False
    return True


def fetch_all_sequences(query_folder, uniprot_filename, database_filename):
    """
    Fetch the fasta formatted sequence for each uniProtID in the list.
    :param query_folder: folder containing fasta files.
    :param uniprot_filename: path to the file containing the uniprot ids
    :param database_filename: path to the file saving the database
    """
    if not query_folder.endswith("/"):
        query_folder += "/"

    print("Processing the list of ids...")

    uniprot_file = open(uniprot_filename, "r")
    database_file = open(database_filename, "w")

    for line in uniprot_file:
        uniprot_id = line.strip()
        ##########################
        ### START CODING HERE ####
        ##########################
        # Fetch the fasta formatted sequence for each uniProtID.
        protein_fasta = fetch_one_fasta(uniprot_id)

        # Store the fasta sequences as individual fasta file in your query directory.
        with open(os.path.join(query_folder, str(uniprot_id) + ".fasta"), "w") as f:
            f.write(protein_fasta)

        # Store all the fasta sequences in one single fasta file as well. These individual files will be used
        # as (PSI-)BLAST queries later on.
        database_file.write(protein_fasta)
        database_file.write("\n")

        ##########################
        ###  END CODING HERE  ####
        ##########################

    print("Processing finished.")
    uniprot_file.close()
    database_file.close()


def main(query_folder, uniprot_file, database_file):
    if check_query_folder(query_folder):
        fetch_all_sequences(query_folder, uniprot_file, database_file)
    else:
        print("Query folder does not exist. "
              "Please make sure that the specified folder exists before you run this script.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='This script downloads all the sequences in the list provided'
                                                 'by the user, and puts them into a single file to be used'
                                                 'as a database. It also stores the fetched fasta sequences into'
                                                 'individual files which will be used as queries for (PSI-)BLAST search')

    parser.add_argument("-ids", "--uniprot_filename", help="the list of UniProt IDs", required=True)
    parser.add_argument("-db", "--database_filename", help="Output file containing all fetched sequences", required=True)
    parser.add_argument("-q", "--qfolder", help="Output folder for your queries", required=True)
    args = parser.parse_args()

    query_folder = args.qfolder
    uniprot_filename = args.uniprot_filename
    database_filename = args.database_filename

    main(query_folder, uniprot_filename, database_filename)
