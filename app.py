import datetime
import logging
import time
from io import BytesIO
import base64
from itertools import product
from urllib.parse import quote

import pandas as pd
import requests
import streamlit as st
from Bio import Entrez
from metapub import FindIt, PubMedFetcher

pd.set_option("display.max_colwidth", 1)

# Streamlit App Setup
st.markdown(
    """
    <div style="text-align: center;">
        <h1>🧬 --- PubMed Retriever --- 🧬\nInvestigate Compound Interactions</h1>
    </div>
    """,
    unsafe_allow_html=True,
)

st.markdown("""
    This tool allows you to search PubMed for combinations of compounds and genes. 
    Results will be processed and displayed here.
""")

# Streamlit Input for User (Main Section)
compounds_input = st.text_area("🧪 Enter Compounds (One Compound per Line)")
compounds_list = [
    compound.strip() for compound in compounds_input.split("\n") if compound.strip()
]

genes_input = st.text_area("🧬 Enter Genes (One Gene per Line) (Optional)")
genes = [gene.strip() for gene in genes_input.split("\n") if gene.strip()]

# New input for additional keywords from user
additional_keywords_input = st.text_area(
    "🔗 Enter Relationship Keywords (One Keyword per Line) (Optional)"
)
additional_keywords_list = [
    keyword.strip()
    for keyword in additional_keywords_input.split("\n")
    if keyword.strip()
]

# Modify the additional condition only if additional keywords are provided
additional_condition = (
    f"AND ({' OR '.join([f'{kw}[Title/Abstract]' for kw in additional_keywords_list])})"
    if additional_keywords_list
    else ""
)

file_ = open("images/logo.png", "rb").read()
base64_image = base64.b64encode(file_).decode("utf-8")

st.sidebar.markdown(
    f"""
    <div style="float: left;">
        <img src="data:image/png;base64,{base64_image}" alt="Logo" width="200" style="border-radius: 5px;">
    </div>
    """,
    unsafe_allow_html=True,
)


# Sidebar for input fields
st.sidebar.write("\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n")
st.sidebar.title("Configuration")
# Ask for Email and NCBI API Key
# email = st.sidebar.text_input("📧 Enter your email")
# api_key = st.sidebar.text_input("🔑 Enter your NCBI API key", type="password")
st.sidebar.markdown(
    """
    <div style="text-align: left;">
        <p style="font-size:14px; color:black;">
            <a href="https://support.nlm.nih.gov/knowledgebase/article/KA-05317/en-us" target="_blank" style="font-size:14px; color:black;">Instructions for creating NCBI API key</a>
        </p>
    </div>
    """,
    unsafe_allow_html=True,
)
top_n = st.sidebar.slider("Select number of top synonyms per compound", 1, 10, 2)
retmax = st.sidebar.slider(
    "Maximum number of PubMed articles to retrieve", 100, 1000, 100
)
start_year = st.sidebar.number_input(
    "Start Year",
    value=2000,
    step=1,
    min_value=1900,
    max_value=datetime.datetime.now().year,
)
end_year = st.sidebar.number_input(
    "End Year",
    value=datetime.datetime.now().year,
    step=1,
    min_value=1900,
    max_value=datetime.datetime.now().year,
)
st.sidebar.write("\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n")
# Sidebar for Configuration
st.sidebar.markdown(
    """
    <div style="display: flex; align-items: center; justify-content: left;">
        <a href="https://github.com/asmaa-a-abdelwahab" target="_blank" style="text-decoration: none;">
            <img src="https://github.githubassets.com/images/modules/logos_page/GitHub-Mark.png" alt="GitHub Logo" style="width:40px; height:40px; margin-right: 10px;">
        </a>
        <a href="https://github.com/asmaa-a-abdelwahab" target="_blank" style="text-decoration: none;">
            <p style="font-size: 16px; font-weight: bold; color: black; margin: 0;">@asmaa-a-abdelwahab</p>
        </a>
    </div>
    """,
    unsafe_allow_html=True,
)


# Cache to prevent redundant queries
@st.cache_data(show_spinner=False)
def cache_synonyms(compound):
    return CompoundResearchHelper().most_common_synonyms(
        compound, top_n
    )  # email, api_key


class CompoundResearchHelper:
    """A class to fetch compound synonyms and retrieve relevant articles from PubMed."""

    def __init__(self, retmax=1000):  # email, api_key,
        # Entrez.email = email  # Set user-provided email
        # Entrez.api_key = api_key  # Set user-provided NCBI API key
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
            cid = cid_data.get("IdentifierList", {}).get("CID", [None])[0]
            if not cid:
                logging.warning(f"No CID found for {compound_name}")
                return []

            synonyms_url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/synonyms/JSON"
            synonyms_response = requests.get(synonyms_url)
            synonyms_response.raise_for_status()

            synonyms_data = synonyms_response.json()
            all_synonyms = (
                synonyms_data.get("InformationList", {})
                .get("Information", [{}])[0]
                .get("Synonym", [])
            )
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
                # Encode the synonym to make it URL safe and wrap it in quotes for PubMed
                encoded_synonym = quote(f'"{synonym}"')
                handle = Entrez.esearch(
                    db="pubmed",
                    term=encoded_synonym,
                    retmax=0,
                    api_key=Entrez.api_key,
                    email=Entrez.email,
                )
                record = Entrez.read(handle)
                synonym_frequency[synonym] = int(record.get("Count", 0))
                handle.close()
                time.sleep(0.3)  # Respect PubMed rate limit
            except Exception as e:
                logging.error(f"Error fetching frequency for {synonym}: {e}")
                continue  # Skip this synonym and move to the next one
        return synonym_frequency

    def most_common_synonyms(self, compound_name, top_n=5):
        """Get the most commonly used synonyms based on literature frequency."""
        if compound_name in self.synonym_data:
            return [synonym for synonym, _ in self.synonym_data[compound_name][:top_n]]

        all_synonyms = self.get_compound_synonyms(compound_name)
        if not all_synonyms:
            return []

        frequency = self.get_synonym_frequency(all_synonyms)
        sorted_synonyms = sorted(
            frequency.items(), key=lambda item: item[1], reverse=True
        )
        self.synonym_data[compound_name] = sorted_synonyms
        return [compound_name] + [synonym for synonym, _ in sorted_synonyms[:top_n]]

    def fetch_articles(self, search_term, retmax=1000, start_year=2000, end_year=None):
        """Fetch articles' metadata from PubMed based on a given search term and date range."""
        date_range = (
            f' AND ("{start_year}/01/01"[PDat] : "{end_year}/12/31"[PDat])'
            if end_year
            else ""
        )
        full_search_term = f"{search_term}{date_range}"

        try:
            pmids = self.pubmed.pmids_for_query(
                full_search_term, retmax=retmax, sort="relevance"
            )
        except Exception as e:
            logging.error(f"Error fetching pmids for {search_term}: {e}")
            return []

        articles = []
        for pmid in pmids:
            try:
                article = self.pubmed.article_by_pmid(pmid)
                article_dict = article.to_dict()
                keys_to_keep = [
                    "title",
                    "pmid",
                    "url",
                    "authors",
                    "doi",
                    "pmc",
                    "issn",
                    "mesh",
                    "chemicals",
                    "journal",
                    "abstract",
                    "year",
                ]
                article_dict = {k: article_dict.get(k, None) for k in keys_to_keep}

                try:
                    time.sleep(1)
                    src = FindIt(pmid, retry_errors=True)
                    if src and src.url:
                        article_dict["url"] = src.url
                except Exception as inner_e:
                    logging.warning(
                        f"Warning while trying to find URL for {pmid}: {inner_e}"
                    )

                articles.append(article_dict)
            except Exception as e:
                logging.error(f"Error processing article with pmid {pmid}: {e}")
        return articles

    def process_compound_and_genes(
        self, compound, genes, start_year, end_year, additional_condition
    ):
        """Process a single compound and multiple genes for synonym lookup and article fetching."""
        logging.info(f"Processing compound: {compound}")
        top_synonyms = cache_synonyms(compound)

        # Handle queries with or without genes
        if genes:
            queries = [
                f"({synonym}[Title/Abstract]) AND ({gene}[Title/Abstract]) {additional_condition}"
                for synonym, gene in product(top_synonyms, genes)
            ]
        else:
            queries = [
                f"({synonym}[Title/Abstract]) {additional_condition}"
                for synonym in top_synonyms
            ]

        # Fetch articles for each combination
        for query in queries:
            self.articleList.extend(
                self.fetch_articles(
                    query, retmax=self.retmax, start_year=start_year, end_year=end_year
                )
            )

        # Return articles for the compound-gene combination
        return self.articleList


# Run when the user clicks the "Search" button
if st.button("🚀 Launch Search"):
    # if not email or not api_key:
    #     st.error("Please provide both email and NCBI API key.")
    if not compounds_list:
        st.error("Please fill out the compound field.")
    else:
        st.info("⏳ Starting article retrieval process...")
        helper = CompoundResearchHelper()  # email, api_key

        # List to store DataFrames for each compound
        all_articles = []

        for compound in compounds_list:
            st.info(f"🔍 Processing compound: {compound}")

            # Ensure that the helper.articleList is reset before processing each compound
            helper.articleList = []  # Clear the article list before processing a new compound

            # Retrieve articles for the current compound
            articles = helper.process_compound_and_genes(
                compound, genes, start_year, end_year, additional_condition
            )

            # Prepare DataFrame to save and display for each compound
            if articles:
                df = pd.DataFrame(articles)

                # Convert any unhashable columns (e.g., lists, dictionaries) to strings
                df = df.applymap(lambda x: str(x) if isinstance(x, (list, dict)) else x)

                # Drop duplicates based on specific columns
                df = df.drop_duplicates(subset=["title", "pmid"], keep="first")

                # Convert 'year' column to numeric, ignoring errors
                if "year" in df.columns:
                    df["year"] = pd.to_numeric(df["year"], errors="coerce")

                df.reset_index(drop=True, inplace=True)

                # Save to a buffer (in-memory file)
                output = BytesIO()
                df.to_csv(output, index=False)
                output.seek(0)

                # Save each compound's articles to a CSV
                filename = f"{compound}_articles.csv"
                with open(filename, "wb") as f:
                    f.write(output.read())

                # Show success for compound
                st.success(f"✔️ Articles processed for compound: {compound}")

                # Optionally display the DataFrame for the compound
                st.dataframe(df)

                # Append the DataFrame to the list
                all_articles.append(df)

            else:
                st.warning(f"No articles found for compound: {compound}")

        # After processing all compounds, combine all DataFrames
        if all_articles:
            combined_df = pd.concat(all_articles, ignore_index=True)

            # Drop duplicates across all combined articles
            combined_df = combined_df.drop_duplicates(
                subset=["title", "pmid"], keep="first"
            )

            # Save the combined DataFrame to a CSV
            combined_df.to_csv("combined_pubmed_articles.csv", index=False)

            # Show the combined DataFrame in Streamlit
            st.subheader("Combined Articles for All Compounds")
            st.dataframe(combined_df)

            st.success("✔️ Combined articles saved to 'combined_pubmed_articles.csv'")
        else:
            st.warning("No articles found for any of the compounds.")
