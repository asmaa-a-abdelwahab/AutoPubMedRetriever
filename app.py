import os
import time
import logging
import requests
import datetime
import pandas as pd
import streamlit as st
from Bio import Entrez
from metapub import FindIt, PubMedFetcher
from itertools import product
import smtplib
from email.mime.multipart import MIMEMultipart
from email.mime.base import MIMEBase
from email import encoders
from email.mime.text import MIMEText
from io import BytesIO

pd.set_option('display.max_colwidth', 1)

# Streamlit App Setup
st.markdown("""
    <div style="text-align: center;">
        <h1>🧬 --- PubMed Retriever --- 🧬\nCompound-Gene Relationships</h1>
    </div>
    """, unsafe_allow_html=True)

st.markdown("""
    This tool allows you to search PubMed for combinations of compounds and genes. 
    Results will be processed and emailed to you.
""")

# Streamlit Input for User (Main Section)
email = st.text_input("📧 Enter your email", value="example@example.com")
compounds_input = st.text_area("🧪 Enter Compounds (One Compound per Line)")
compounds_list = [compound.strip() for compound in compounds_input.split("\n")]

genes_input = st.text_area("🧬 Enter Genes (One Gene per Line)")
genes = [gene.strip() for gene in genes_input.split("\n")]

# New input for additional keywords from user
additional_keywords_input = st.text_area("🔗 Enter Relationship Keywords (One Keyword per Line)")

# Convert the additional keywords into a PubMed-compatible OR query
additional_keywords_list = [keyword.strip() for keyword in additional_keywords_input.split("\n")]
additional_condition = f"AND ({' OR '.join([f'{kw}[Title/Abstract]' for kw in additional_keywords_list])})"

# Sidebar for Configuration
st.sidebar.markdown(
    """
    <div style="float: left;">
        <img src="https://raw.githubusercontent.com/asmaa-a-abdelwahab/AIGraphQuery-/main/EwC%20full%20logo.png" alt="Logo" width="100" style="border-radius: 1px; float: left;">
    </div>
    """, 
    unsafe_allow_html=True
)

# Sidebar for input fields
st.sidebar.write("\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n")
st.sidebar.title("Configuration")
top_n = st.sidebar.slider("Select number of top synonyms per compound", 1, 10, 5)
retmax = st.sidebar.slider("Maximum number of PubMed articles to retrieve", 50, 100, 1000, 500)
start_year = st.sidebar.number_input("Start Year", value=2000, step=1, min_value=1900, max_value=datetime.datetime.now().year)
end_year = st.sidebar.number_input("End Year", value=datetime.datetime.now().year, step=1, min_value=1900, max_value=datetime.datetime.now().year)
st.sidebar.write("\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n")
st.sidebar.markdown(
    """
    <div style="position: absolute; bottom: 10; width: 100%; text-align: left;">
        <p style="font-size:14px; color:black;">BY: 
            <a href="https://github.com/asmaa-a-abdelwahab" target="_blank" style="font-size:14px; color:black;">Asmaa A. Abdelwahab</a>
        </p>
    </div>
    """, 
    unsafe_allow_html=True
)

# Cache to prevent redundant queries
@st.cache_data(show_spinner=False)
def cache_synonyms(compound):
    return CompoundResearchHelper(email).most_common_synonyms(compound, top_n)

class CompoundResearchHelper:
    """A class to fetch compound synonyms and retrieve relevant articles from PubMed."""
    
    def __init__(self, email, retmax=1000):
        self.email = email
        Entrez.email = email
        self.pubmed = PubMedFetcher()
        self.retmax = retmax
        self.synonym_data = {}
        self.articleList = []
        self.current_year = datetime.datetime.now().year

    def get_compound_synonyms(self, compound_name):
        """Retrieve compound synonyms from PubChem."""
        try:
            cid_url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{compound_name}/cids/JSON"
            cid_response = requests.get(cid_url)
            cid_response.raise_for_status()
            
            cid_data = cid_response.json()
            cid = cid_data.get('IdentifierList', {}).get('CID', [None])[0]
            if not cid:
                logging.warning(f"No CID found for {compound_name}")
                return []

            synonyms_url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/synonyms/JSON"
            synonyms_response = requests.get(synonyms_url)
            synonyms_response.raise_for_status()

            synonyms_data = synonyms_response.json()
            all_synonyms = synonyms_data.get('InformationList', {}).get('Information', [{}])[0].get('Synonym', [])
            return all_synonyms or []

        except requests.HTTPError as http_err:
            logging.error(f"HTTP error occurred: {http_err}")
        except Exception as err:
            logging.error(f"An error occurred: {err}")
        return []

    def get_synonym_frequency(self, synonyms):
        """Check the frequency of each synonym in PubMed literature."""
        synonym_frequency = {}
        for synonym in synonyms:
            try:
                handle = Entrez.esearch(db="pubmed", term=synonym, retmax=0)
                record = Entrez.read(handle)
                synonym_frequency[synonym] = int(record.get("Count", 0))
                handle.close()
                time.sleep(0.3)  # Respect PubMed rate limit
            except Exception as e:
                logging.error(f"Error fetching frequency for {synonym}: {e}")
        return synonym_frequency

    def most_common_synonyms(self, compound_name, top_n=5):
        """Get the most commonly used synonyms based on literature frequency."""
        if compound_name in self.synonym_data:
            return [synonym for synonym, _ in self.synonym_data[compound_name][:top_n]]

        all_synonyms = self.get_compound_synonyms(compound_name)
        if not all_synonyms:
            return []

        frequency = self.get_synonym_frequency(all_synonyms)
        sorted_synonyms = sorted(frequency.items(), key=lambda item: item[1], reverse=True)
        self.synonym_data[compound_name] = sorted_synonyms
        return [compound_name] + [synonym for synonym, _ in sorted_synonyms[:top_n]]

    def fetch_articles(self, search_term, retmax=1000, start_year=2000, end_year=None):
        """Fetch articles' metadata from PubMed based on a given search term and date range."""
        date_range = f" AND (\"{start_year}/01/01\"[PDat] : \"{end_year}/12/31\"[PDat])" if end_year else ""
        full_search_term = f"{search_term}{date_range}"
        
        try:
            pmids = self.pubmed.pmids_for_query(full_search_term, retmax=retmax, sort='relevance')
        except Exception as e:
            logging.error(f"Error fetching pmids for {search_term}: {e}")
            return []
        
        articles = []
        for pmid in pmids:
            try:
                article = self.pubmed.article_by_pmid(pmid)
                article_dict = article.to_dict()
                keys_to_keep = ['title', 'pmid', 'url', 'authors', 'doi', 'pmc', 'issn', 'mesh', 'chemicals', 'journal', 'abstract', 'year']
                article_dict = {k: article_dict.get(k, None) for k in keys_to_keep}

                try:
                    time.sleep(1)
                    src = FindIt(pmid, retry_errors=True)
                    if src and src.url:
                        article_dict['url'] = src.url
                except Exception as inner_e:
                    logging.warning(f"Warning while trying to find URL for {pmid}: {inner_e}")

                articles.append(article_dict)
            except Exception as e:
                logging.error(f"Error processing article with pmid {pmid}: {e}")
        return articles

    def process_compound_and_genes(self, compound, genes, start_year, end_year, additional_condition):
        """Process a single compound and multiple genes for synonym lookup and article fetching."""
        logging.info(f"Processing compound: {compound}")
        top_synonyms = cache_synonyms(compound)

        # Combine each synonym with every gene
        queries = [f"({synonym}[Title/Abstract]) AND ({gene}[Title/Abstract]) {additional_condition}" 
                for synonym, gene in product(top_synonyms, genes)]

        # Fetch articles for each combination
        for query in queries:
            self.articleList.extend(self.fetch_articles(query, retmax=self.retmax, start_year=start_year, end_year=end_year))

        # Return articles for the compound-gene combination
        return self.articleList

# Run when the user clicks the "Search" button
if st.button("🚀 Launch Search"):
    if not email or not compounds_list or not genes:
        st.error("Please fill out all fields.")
    else:
        st.info("⏳ Starting article retrieval process...")
        helper = CompoundResearchHelper(email)

        for compound in compounds_list:
            st.info(f"🔍 Processing compound: {compound}")
            articles = helper.process_compound_and_genes(compound, genes, start_year, end_year, additional_condition)

            # Prepare DataFrame to send as CSV attachment
            if articles:
                df = pd.DataFrame(articles)

                # Convert any unhashable columns (e.g., lists, dictionaries) to strings
                df = df.map(lambda x: str(x) if isinstance(x, (list, dict)) else x)

                # Drop duplicates based on specific columns that are less likely to have complex data
                df = df.drop_duplicates(subset=['title', 'pmid'], keep='first')

                # Convert 'year' column to numeric, ignoring errors
                if 'year' in df.columns:
                    df['year'] = pd.to_numeric(df['year'], errors='coerce')

                df.reset_index(drop=True, inplace=True)

                # Save to a buffer (in-memory file)
                output = BytesIO()
                df.to_csv(output, index=False)
                output.seek(0)

                # Save the CSV temporarily
                filename = f"{compound}_articles.csv"
                with open(filename, "wb") as f:
                    f.write(output.read())

                # Show success for compound
                st.success(f"✔️ Articles processed for compound: {compound}")

                # Optionally display the DataFrame
                st.dataframe(df)

                # Cleanup: Remove the file from the server
                if os.path.exists(filename):
                    os.remove(filename)

            else:
                st.warning(f"No articles found for compound: {compound}")

        # Show the articles on the Streamlit app (optional if you want real-time display)
        if len(helper.articleList) > 0:
            st.subheader("Fetched Articles for Compounds and Genes")

            # Convert article list to a DataFrame
            df = pd.DataFrame(helper.articleList)

            # Convert unhashable columns to strings
            df = df.map(lambda x: str(x) if isinstance(x, (list, dict)) else x)

            # Drop duplicates based on specific columns
            df = df.drop_duplicates(subset=['title', 'pmid'], keep='first')

            # Convert 'year' column to numeric
            if 'year' in df.columns:
                df['year'] = pd.to_numeric(df['year'], errors='coerce')

            df.reset_index(drop=True, inplace=True)

            # Save the cleaned DataFrame to a CSV file
            df.to_csv('pubmed_articles.csv', index=False)

            # Display the cleaned DataFrame in Streamlit
            st.dataframe(df)

        else:
            st.warning("No articles found for the given combinations.")