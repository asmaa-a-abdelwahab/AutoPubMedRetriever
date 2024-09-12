import os
import time
import csv
import requests
import requests
import pandas as pd
from Bio import Entrez
from metapub import FindIt, PubMedFetcher

pd.set_option('display.max_colwidth', 1)

class CompoundSynonymChecker:

    def __init__(self, email):
        """Initialize with an email address to use with Entrez."""
        self.email = email
        Entrez.email = email
        self.synonym_data = {}  # Cache for storing synonym data

    def get_compound_synonyms(self, compound_name):
        """Retrieve compound synonyms from PubChem."""
        try:
            cid_url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{compound_name}/cids/JSON"
            cid_response = requests.get(cid_url)
            cid_response.raise_for_status()

            cid_data = cid_response.json()
            cid = cid_data.get('IdentifierList', {}).get('CID', [None])[0]
            if cid is None:
                raise ValueError(f"No CID found for {compound_name}")

            synonyms_url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/synonyms/JSON"
            synonyms_response = requests.get(synonyms_url)
            synonyms_response.raise_for_status()

            synonyms_data = synonyms_response.json()
            all_synonyms = synonyms_data.get('InformationList', {}).get('Information', [{}])[0].get('Synonym', [])
            if not all_synonyms:
                raise ValueError(f"No synonyms found for {compound_name}")

            return all_synonyms

        except requests.HTTPError as http_err:
            print(f"HTTP error occurred: {http_err}")
        except Exception as err:
            print(f"An error occurred: {err}")
        return []

    def get_synonym_frequency(self, synonyms):
        """Check the frequency of each synonym in PubMed literature."""
        synonym_frequency = {}
        for synonym in synonyms:
            try:
                handle = Entrez.esearch(db="pubmed", term=synonym, retmax=0)
                record = Entrez.read(handle)
                handle.close()
                synonym_frequency[synonym] = int(record.get("Count", 0))
                time.sleep(0.3)  # Sleep to respect PubMed's rate limit policy
            except Exception as e:
                print(f"Error fetching frequency for {synonym}: {e}")
        return synonym_frequency

    def most_common_synonyms(self, compound_name, top_n=5):
        """Get the most commonly used synonyms based on literature frequency."""
        if compound_name in self.synonym_data:
            # Return the top N synonyms from the cache
            return [synonym for synonym, _ in self.synonym_data[compound_name][:top_n]]

        all_synonyms = self.get_compound_synonyms(compound_name)
        if not all_synonyms:
            return []

        frequency = self.get_synonym_frequency(all_synonyms)
        sorted_synonyms = sorted(frequency.items(), key=lambda item: item[1], reverse=True)
        self.synonym_data[compound_name] = sorted_synonyms  # Cache the full list
        return [synonym for synonym, _ in sorted_synonyms[:top_n]]

    def save_to_csv(self, compound_name, top_n=5, file_name='synonyms.csv'):
        """Save the most common synonyms to a CSV file."""
        if compound_name not in self.synonym_data:
            # Populate the cache if it hasn't been done
            self.most_common_synonyms(compound_name, top_n)

        sorted_synonyms = self.synonym_data[compound_name][:top_n]

        # Write to CSV using the cached data
        with open(file_name, mode='w', newline='', encoding='utf-8') as file:
            writer = csv.writer(file)
            writer.writerow(['Compound', 'Synonym', 'Frequency'])
            for synonym, freq in sorted_synonyms:
                writer.writerow([compound_name, synonym, freq])
        print(f"Saved the top {top_n} synonyms for {compound_name} to {file_name}")



class ChemicalInteractionsLiteratueCollector:
    """
    A class used to fetch and download articles from PubMed based on specific search terms.
    """
    def __init__(self):
        """
        Initializes a new instance of the PubMedArticlesFetcher class.
        """
        self.pubmed = PubMedFetcher()
        self.articleList = []


    def element_in_text(self, elements, text):
        """
        Checks if any element of a list is present in the specified text.

        Parameters:
        elements (list[str]): List of substrings to check.
        text (str): Text in which to check for substrings.

        Returns:
        bool: True if any element is found in text, otherwise False.
        """
        text = text.lower()
        return any(element.lower() in text for element in elements)


    def fetch_articles(self, search_term: str, synonyms, retmax: int = 1000):
        """
        Fetch articles' metadata from PubMed based on a given search term.

        Parameters:
        search_term (str): The term to search articles for.
        retmax (int): The maximum number of articles to retrieve.
        """
        self.search_term = search_term

        # Try to fetch pmids, log error if not successful
        try:
            pmids = self.pubmed.pmids_for_query(search_term, retmax=retmax, sort='relevance')
        except Exception as e:
            logging.error(f"Error fetching pmids for {search_term}: {e}")
            return

        # Iterate over each pmid
        for pmid in pmids:
            try:
                article = self.pubmed.article_by_pmid(pmid)
                articleDict = article.to_dict()
                # Only keep necessary keys
                keys_to_keep = ['title','pmid','url','authors','doi','pmc','issn', 'mesh','chemicals','journal','abstract','year']
                articleDict = {k: articleDict[k] for k in keys_to_keep}

                try:
                    time.sleep(1)  # Consider more intelligent rate-limiting
                    src = FindIt(pmid, retry_errors=True)
                    if src and src.url:
                        articleDict['url'] = src.url
                except Exception as inner_e:
                    logging.warning(f"Warning while trying to find URL for {pmid}: {inner_e}")

                # If any synonym is in the article abstract, append it
                # if self.element_in_text(synonyms, articleDict['abstract']):
                self.articleList.append(articleDict)
            except Exception as e:
                logging.error(f"Error processing article with pmid {pmid}: {e}")

    def save_to_csv(self, filename: str):
        """
        Save the fetched articles' metadata into a CSV file.

        Parameters:
        filename (str): The name of the CSV file to save data into.
        """
        self.df = pd.DataFrame(self.articleList)
        self.df.reset_index(inplace=True, drop=True)  # Reset index and drop the old one immediately
        self.df = self.df.astype(str).drop_duplicates()  # Convert to string and then drop duplicates

        # Save the CSV to the directory named after the search term
        self.df.to_csv(filename, index=False)


# Usage: (Time  ~ 1h)
email = "asmaa@edelweissconnect.com"
checker = CompoundSynonymChecker(email)

# Initialize an empty dictionary to store compound synonyms
compound_synonyms_dict = {}

# # Loop over each compound, cache the results, and add them to the dictionary
# for compound in compounds:

#     top_synonyms = checker.most_common_synonyms(compound, top_n=10)
#     compound_synonyms_dict[compound] = top_synonyms

#     # Optionally save each to CSV as well
#     checker.save_to_csv(compound, top_n=100, file_name=f'Synonyms/{compound}_synonyms_frequency.csv')

# # At this point, compound_synonyms_dict contains each compound and their top synonyms.
# print(compound_synonyms_dict)