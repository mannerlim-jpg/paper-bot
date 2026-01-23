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

# ==========================================
# [기능 2] Gemini 임상 판단 (속도 최적화)
# ==========================================
def summarize_paper(title, abstract):
    if not GEMINI_API_KEY: return "N", "API 키 없음", ""

    prompt = f"""
    당신은 정형외과 전문의입니다.
    
    [분석 지침]
    1. 이 논문이 기존 임상 관행(Practice)을 바꾸거나 도전하는 내용입니까? (Yes/No)
    2. 외래 진료 중 환자에게 30초 안에 설명한다면? (구어체 한 문장)
    3. 상세 리뷰 작성 (배경/결과/비판)

    [출력 형식] (구분자 @ 사용)
    (Y/N) @ (환자용 한 문장) @ (상세 리뷰)

    [논문] {title} / {abstract}
    """

    for attempt in range(2): 
        try:
            model = genai.GenerativeModel('gemini-2.5-flash')
            response = model.generate_content(prompt)
            text = response.text
            
            parts = text.split("@")
            if len(parts) >= 3:
                return parts[0].strip().upper(), parts[1].strip(), parts[2].strip()
            else:
                return "N", "형식 오류", text
        except Exception as e:
            print(f"⚠️ 에러 발생. 재시도 중... ({e})")
            time.sleep(5) # 재시도 대기는 5초로 단축
            if attempt == 1: return "N", "요약 실패", f"에러: {e}"

# ==========================================
# [기능 3] 이메일 전송
# ==========================================
def send_email(content_html):
    if not MY_EMAIL or not MY_APP_PASSWORD: return

    msg = MIMEText(content_html, 'html')
    today = datetime.now().strftime('%Y-%m-%d')
    msg['Subject'] = f"📢 [Dr.AI] {today} 정형외과 논문 브리핑 (06:00)"
    msg['From'] = MY_EMAIL
    msg['To'] = RECEIVER_EMAIL

    with smtplib.SMTP_SSL('smtp.gmail.com', 465) as server:
        server.login(MY_EMAIL, MY_APP_PASSWORD)
        server.send_message(msg)
    print("✅ 이메일 발송 완료")

# ==========================================
# [실행] 메인 컨트롤러
# ==========================================
def main():
    html_body = """
    <h2 style='color:#2c3e50;'>👨‍⚕️ 오늘의 정형외과 논문 인사이트</h2>
    <p style='color:gray; font-size:12px;'>* 🚨: 임상 관행 변화 가능성 있음</p>
    <hr>
    """
    total_papers = 0

    for keyword in SEARCH_KEYWORDS:
        papers = fetch_pubmed_papers(keyword, max_results=5)
        html_body += f"<h3 style='color: #2980b9; margin-top:30px;'>🔍 {keyword}</h3>"

        if not papers:
            html_body += "<p style='color: #bdc3c7;'>신규 논문 없음</p>"
            continue

        for i, paper in enumerate(papers, 1):
            # [속도 개선] 30초 -> 10초로 단축!
            if total_papers > 0: 
                print(f"[{total_papers}번째 완료] 10초 대기 중...")
                time.sleep(10)

            impact, one_liner, deep_rev = summarize_paper(paper['title'], paper['abstract'])
            
            # 배지 디자인
            if "Y" in impact:
                badge = "<span style='background-color:#e74c3c; color:white; padding:3px 8px; border-radius:4px; font-size:12px;'>🚨 주목</span>"
                bg_color = "#fff5f5"
            else:
                badge = "<span style='background-color:#95a5a6; color:white; padding:3px 8px; border-radius:4px; font-size:12px;'>✅ 일반</span>"
                bg_color = "#ffffff"

            html_body += f"""
            <div style='background-color: {bg_color}; border: 1px solid #dfe6e9; border-radius: 8px; margin-bottom: 15px; padding: 15px;'>
                <p style='font-size: 16px; font-weight: bold; margin: 0 0 5px 0;'>
                    {badge} <a href='{paper['link']}' target='_blank' style='text-decoration: none; color: #2c3e50;'>[{i}] {paper['title']}</a>
                </p>
                <div style='background-color: #ecf0f1; padding: 10px; border-left: 4px solid #3498db; font-size: 14px; margin: 10px 0;'>
                    <strong>🗣️ 환자용:</strong> "{one_liner}"
                </div>
                <details>
                    <summary style='cursor: pointer; color: #7f8c8d; font-size: 13px;'>🔽 상세 리뷰 보기</summary>
                    <div style='padding-top: 10px; font-size: 14px; line-height: 1.6; color: #444;'>{deep_rev.replace('\n', '<br>')}</div>
                </details>
            </div>
            """
            total_papers += 1

    if total_papers > 0: send_email(html_body)

if __name__ == "__main__":
    main()
