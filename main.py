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

print("🚀 [진단 시작] 봇 가동 중...")

# 환경변수 확인 (비밀번호는 숨기고, 존재 여부만 확인)
MY_EMAIL = os.getenv("MY_EMAIL")
MY_APP_PASSWORD = os.getenv("MY_APP_PASSWORD")
GEMINI_API_KEY = os.getenv("GEMINI_API_KEY")
RECEIVER_EMAIL = os.getenv("RECEIVER_EMAIL", MY_EMAIL)

if not MY_EMAIL: print("❌ 경고: MY_EMAIL 환경변수가 없습니다!")
if not MY_APP_PASSWORD: print("❌ 경고: MY_APP_PASSWORD 환경변수가 없습니다!")
if not GEMINI_API_KEY: print("❌ 경고: GEMINI_API_KEY 환경변수가 없습니다!")
else:
    print("✅ 환경변수 확인 완료.")
    genai.configure(api_key=GEMINI_API_KEY.strip())

# Entrez 이메일 설정 (필수)
if MY_EMAIL:
    Entrez.email = MY_EMAIL
else:
    Entrez.email = "test@test.com" # 비상용 더미 이메일

# ==========================================
# [기능 1] 논문 검색 (디버깅 모드)
# ==========================================
def fetch_pubmed_papers(keyword, max_results=5):
    print(f"🔎 검색 시도: '{keyword}' (최근 30일)")
    
    try:
        # 30일 검색으로 고정
        handle = Entrez.esearch(db="pubmed", term=keyword, retmax=max_results, 
                                sort="relevance", reldate=30, datetype="pdat")
        record = Entrez.read(handle)
        handle.close()
    except Exception as e:
        print(f"🚨 [검색 에러] PubMed 접속 실패: {e}")
        return []
    
    id_list = record["IdList"]
    print(f"   👉 발견된 논문 ID: {id_list}")
    
    if not id_list:
        return []

    papers = []
    try:
        handle = Entrez.efetch(db="pubmed", id=id_list, retmode="xml")
        records = Entrez.read(handle)
        handle.close()
    except Exception as e:
        print(f"🚨 [다운로드 에러] 논문 상세 정보 가져오기 실패: {e}")
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
        except Exception as e:
            print(f"⚠️ 논문 파싱 중 건너뜀: {e}")
            continue
            
    print(f"   ✅ 처리 완료: {len(papers)}건")
    return papers

# ==========================================
# [기능 2] Gemini 요약
# ==========================================
def summarize_paper(title, abstract):
    if not GEMINI_API_KEY: return "N", "키 없음", "내용 없음"

    prompt = f"""
    당신은 정형외과 전문의입니다.
    [분석 지침]
    1. 이 논문이 기존 임상 관행을 바꾸거나 도전합니까? (Yes/No)
    2. 환자에게 30초 설명 (구어체 한 문장)
    3. 상세 리뷰 (배경/결과/비판)
    [출력 형식] (Y/N) @ (한문장) @ (상세리뷰)
    [논문] {title} / {abstract}
    """

    try:
        model = genai.GenerativeModel('gemini-2.5-flash')
        response = model.generate_content(prompt)
        text = response.text
        parts = text.split("@")
        if len(parts) >= 3:
            return parts[0].strip().upper(), parts[1].strip(), parts[2].strip()
        else:
            return "N", "형식 에러", text
    except Exception as e:
        print(f"⚠️ [AI 에러] : {e}")
        return "N", "요약 실패", str(e)

# ==========================================
# [기능 3] 이메일 전송
# ==========================================
def send_email(content_html):
    if not MY_EMAIL or not MY_APP_PASSWORD: 
        print("❌ 이메일 전송 실패: 계정 정보가 없습니다.")
        return

    print("📧 이메일 발송 시도 중...")
    msg = MIMEText(content_html, 'html')
    today = datetime.now().strftime('%Y-%m-%d')
    msg['Subject'] = f"📢 [Dr.AI] {today} 정형외과 논문 브리핑 (TEST)"
    msg['From'] = MY_EMAIL
    msg['To'] = RECEIVER_EMAIL

    try:
        with smtplib.SMTP_SSL('smtp.gmail.com', 465) as server:
            server.login(MY_EMAIL, MY_APP_PASSWORD)
            server.send_message(msg)
        print("✅ 이메일 발송 성공!")
    except Exception as e:
        print(f"❌ 이메일 발송 에러: {e}")

# ==========================================
# [실행]
# ==========================================
def main():
    html_body = "<h2>👨‍⚕️ 진단 모드 리포트</h2><hr>"
    total_papers = 0

    for keyword in SEARCH_KEYWORDS:
        papers = fetch_pubmed_papers(keyword, max_results=2
