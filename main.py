import smtplib
import os
import time
from email.mime.text import MIMEText
from Bio import Entrez
import google.generativeai as genai
from datetime import datetime

# ==========================================
# [설정] 검색어 목록
# ==========================================
SEARCH_KEYWORDS = [
    "(Total Knee Replacement) AND (Robotic)",
    "Total Knee Arthroplasty",
    "Foot",
    "Ankle Instability",
    "(Ankle) AND (Arthroscopy)",
    "(Knee) AND (Arthroscopy)",
    "Arthroscopy"
]

MY_EMAIL = os.getenv("MY_EMAIL")
MY_APP_PASSWORD = os.getenv("MY_APP_PASSWORD")
GEMINI_API_KEY = os.getenv("GEMINI_API_KEY")
RECEIVER_EMAIL = os.getenv("RECEIVER_EMAIL", MY_EMAIL)

if GEMINI_API_KEY:
    genai.configure(api_key=GEMINI_API_KEY.strip())
Entrez.email = MY_EMAIL 

# ==========================================
# [기능 1] 논문 검색
# ==========================================
def fetch_pubmed_papers(keyword, max_results=5):
    try:
        handle = Entrez.esearch(db="pubmed", term=keyword, retmax=max_results, 
                                sort="relevance", reldate=2, datetype="pdat")
        record = Entrez.read(handle)
        handle.close()
    except:
        return []
    
    id_list = record["IdList"]
    papers = []
    if not id_list: return []

    try:
        handle = Entrez.efetch(db="pubmed", id=id_list, retmode="xml")
        records = Entrez.read(handle)
        handle.close()
    except:
        return []

    for article in records['PubmedArticle']:
        try:
            citation = article['MedlineCitation']
            article_data = citation['Article']
            title = article_data['ArticleTitle']
            journal = article_data['Journal'].get('Title', 'Unknown Journal')
            pmid = citation['PMID']
            link = f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/"

            if 'Abstract' in article_data:
                abstract_list = article_data['Abstract']['AbstractText']
                abstract = " ".join(abstract_list) if isinstance(abstract_list, list) else str(abstract_list)
            else:
                abstract = "초록 없음"
            
            papers.append({"title": title, "abstract": abstract, "journal": journal, "link": link})
        except:
            continue
    return papers
