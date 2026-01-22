import smtplib
import os
from email.mime.text import MIMEText
from Bio import Entrez
import google.generativeai as genai
from datetime import datetime

# ==========================================
# [설정 구역]
# ==========================================

# 1. 검색어 목록 (자유롭게 추가/삭제)
SEARCH_KEYWORDS = [
    "(Total Knee Replacement) AND (Robotic)",
    "(Total Hip Replacement) AND (Direct Anterior Approach)",
    "(Ankle Instability) AND (Reconstruction)"
]

# 2. 내 정보 가져오기 (Github Secrets)
MY_EMAIL = os.getenv("MY_EMAIL")
MY_APP_PASSWORD = os.getenv("MY_APP_PASSWORD")
GEMINI_API_KEY = os.getenv("GEMINI_API_KEY")
RECEIVER_EMAIL = os.getenv("RECEIVER_EMAIL", MY_EMAIL)

# 3. 기본 설정
if GEMINI_API_KEY:
    genai.configure(api_key=GEMINI_API_KEY)
Entrez.email = MY_EMAIL 

# ==========================================
# [기능 1] 논문 검색 (링크, 저널명 포함)
# ==========================================
def fetch_pubmed_papers(keyword, max_results=3):
    try:
        # reldate=2 (최근 2일)로 설정하여 매일 새로운 논문 위주로 검색
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
            # 기본 정보 추출
            citation = article['MedlineCitation']
            article_data = citation['Article']
            
            title = article_data['ArticleTitle']
            
            # [추가됨] 저널 이름 가져오기
            journal = article_data['Journal'].get('Title', 'Unknown Journal')
            
            # [추가됨] PMID로 링크 만들기
            pmid = citation['PMID']
            link = f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/"

            # 초록 추출
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
# [기능 2] Gemini 요약
# ==========================================
def summarize_paper(title, abstract):
    if not GEMINI_API_KEY:
        return "API 키가 없습니다."

    model = genai.GenerativeModel('gemini-pro')
    prompt = f"""
    아래 의학 논문을 한국어로 핵심만 3줄 요약하세요.
    제목: {title}
    초록: {abstract}
    
    형식:
    1. 목적:
    2. 결과:
    3. 의의:
    """
    try:
        response = model.generate_content(prompt)
        return response.text
    except Exception:
        return "요약 실패 (내용이 너무 길거나 API 오류)"

# ==========================================
# [기능 3] 이메일 전송 (HTML 디자인 적용)
# ==========================================
def send_email(content_html):
    if not MY_EMAIL or not MY_APP_PASSWORD:
        return

    # 일반 텍스트가 아니라 HTML 형식을 사용한다고 선언
    msg = MIMEText(content_html, 'html')
    
    today = datetime.now().strftime('%Y-%m-%d')
    msg['Subject'] = f"📢 [매일 논문] {today} 분야별 최신 의학 리포트"
    msg['From'] = MY_EMAIL
    msg['To'] = RECEIVER_EMAIL

    try:
        with smtplib.SMTP_SSL('smtp.gmail.com', 465) as server:
            server.login(MY_EMAIL, MY_APP_PASSWORD)
            server.send_message(msg)
        print("✅ 이메일 발송 완료")
    except Exception as e:
        print(f"❌ 발송 실패: {e}")

# ==========================================
# [실행] 메인 로직
# ==========================================
def main():
    # 이메일 본문을 HTML로 예쁘게 작성
    html_body = "<h2>📅 최근 48시간 이내 주요 논문 리포트</h2><hr>"
    total_papers_found = 0

    for keyword in SEARCH_KEYWORDS:
        papers = fetch_pubmed_papers(keyword, max_results=3)
        
        # 주제별 헤더
        html_body += f"<h3 style='color: #2E86C1;'>🔍 주제: {keyword}</h3>"

        if not papers:
            html_body += "<p>Running... (최근 48시간 내 새로운 논문 없음)</p><br>"
            continue

        for i, paper in enumerate(papers, 1):
            summary = summarize_paper(paper['title'], paper['abstract'])
            
            # 논문 1개 박스 디자인
            html_body += f"""
            <div style='background-color: #f9f9f9; padding: 15px; border-radius: 5px; margin-bottom: 20px;'>
                <p style='font-size: 16px; font-weight: bold;'>
                    <a href='{paper['link']}' style='text-decoration: none; color: #1e4f91;'>
                        [{i}] {paper['title']}
                    </a>
                </p>
                <p style='color: gray; font-size: 12px;'>📖 Journal: {paper['journal']}</p>
                <div style='background-color: #ffffff; padding: 10px; border: 1px solid #ddd;'>
                    {summary.replace('\n', '<br>')}
                </div>
            </div>
            """
            total_papers_found += 1
        
        html_body += "<br>"

    # 논문이 하나라도 있으면 메일 전송
    if total_papers_found > 0:
        send_email(html_body)
    else:
        print("새로운 논문이 없어 메일을 보내지 않았습니다.")

if __name__ == "__main__":
    main()
