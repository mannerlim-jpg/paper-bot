import smtplib
import os
from email.mime.text import MIMEText
from Bio import Entrez
import google.generativeai as genai
from datetime import datetime

# ==========================================
# [설정 구역] 요청하신 검색어 추가 완료
# ==========================================

SEARCH_KEYWORDS = [
    # 1. 기존 핵심 주제
    "(Total Knee Replacement) AND (Robotic)",
    
    # 2. 새로 추가하신 주제
    "Total Knee Arthroplasty",  # 인공슬관절 전치환술 (일반)
    "Foot",                     # 발 (족부)
    
    # 3. 지난번 추가했던 세부 주제 (계속 유지)
    "Ankle Instability",
    "(Ankle) AND (Arthroscopy)",
    "(Knee) AND (Arthroscopy)",
    "Arthroscopy"
]

# 환경변수 로드
MY_EMAIL = os.getenv("MY_EMAIL")
MY_APP_PASSWORD = os.getenv("MY_APP_PASSWORD")
GEMINI_API_KEY = os.getenv("GEMINI_API_KEY")
RECEIVER_EMAIL = os.getenv("RECEIVER_EMAIL", MY_EMAIL)

# [핵심] API 키 공백 제거 (에러 방지용)
if GEMINI_API_KEY:
    clean_key = GEMINI_API_KEY.strip()
    genai.configure(api_key=clean_key)
Entrez.email = MY_EMAIL 

# ==========================================
# [기능 1] 논문 검색
# ==========================================
def fetch_pubmed_papers(keyword, max_results=5):
    try:
        # 최근 48시간(2일) 이내 논문만 검색 (중복 방지)
        handle = Entrez.esearch(db="pubmed", term=keyword, retmax=max_results, 
                                sort="relevance", reldate=2, datetype="pdat")
        record = Entrez.read(handle)
        handle.close()
    except Exception:
        return []
    
    id_list = record["IdList"]
    papers = []
    
    if not id_list:
        return []

    try:
        handle = Entrez.efetch(db="pubmed", id=id_list, retmode="xml")
        records = Entrez.read(handle)
        handle.close()
    except Exception:
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
                
            papers.append({
                "title": title, 
                "abstract": abstract, 
                "journal": journal, 
                "link": link
            })
        except Exception:
            continue
    return papers

# ==========================================
# [기능 2] Gemini 요약 (2.5 버전 적용)
# ==========================================
def summarize_paper(title, abstract):
    if not GEMINI_API_KEY:
        return "오류: API 키가 없습니다."

    # [수정됨] 선생님 키에서 사용 가능한 최신 모델 'gemini-2.5-flash' 적용
    try:
        model = genai.GenerativeModel('gemini-2.5-flash')
        
        prompt = f"""
        아래 정형외과 논문을 한국어로 핵심만 3줄 요약하세요.
        제목: {title}
        초록: {abstract}
        형식:
        1. 목적:
        2. 결과:
        3. 의의:
        """
        response = model.generate_content(prompt)
        return response.text
        
    except Exception as e:
        # 혹시 몰라 2.0 버전으로 재시도하는 안전장치 추가
        try:
            model = genai.GenerativeModel('gemini-2.0-flash')
            response = model.generate_content(prompt)
            return response.text
        except:
            return f"요약 실패 (에러: {e})"

# ==========================================
# [기능 3] 이메일 전송
# ==========================================
def send_email(content_html):
    if not MY_EMAIL or not MY_APP_PASSWORD:
        return

    msg = MIMEText(content_html, 'html')
    today = datetime.now().strftime('%Y-%m-%d')
    msg['Subject'] = f"📢 [매일 아침] {today} 정형외과 최신 논
