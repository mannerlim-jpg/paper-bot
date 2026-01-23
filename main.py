import smtplib
import os
from email.mime.text import MIMEText
from Bio import Entrez
import google.generativeai as genai
from datetime import datetime

# ==========================================
# [설정 구역]
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

# API 키 공백 제거
if GEMINI_API_KEY:
    clean_key = GEMINI_API_KEY.strip()
    genai.configure(api_key=clean_key)
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
# [기능 2] Gemini 요약 (2.5 버전)
# ==========================================
def summarize_paper(title, abstract):
    if not GEMINI_API_KEY:
        return "오류: API 키가 없습니다."

    try:
        # 선생님 키에 맞는 최신 모델 (2.5-flash)
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
        # 혹시 실패하면 2.0으로 재시도
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
    
    # [수정 포인트] 이 부분이 잘리지 않도록 주의하세요!
    msg['Subject'] = f"📢 [매일 아침] {today} 정형외과 최신 논문 리포트"
    
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
# [실행] 메인 컨트롤러
# ==========================================
def main():
    html_body = "<h2>📅 최근 48시간 내 발표된 주요 논문</h2><hr>"
    total_papers_found = 0

    for keyword in SEARCH_KEYWORDS:
        papers = fetch_pubmed_papers(keyword, max_results=5)
        html_body += f"<h3 style='color: #2E86C1; border-left: 5px solid #2E86C1; padding-left: 10px;'>🔍 {keyword}</h3>"

        if not papers:
            html_body += "<p style='color: gray; font-style: italic;'>최근 2일간 새로 등록된 논문이 없습니다.</p><br>"
            continue

        for i, paper in enumerate(papers, 1):
            summary = summarize_paper(paper['title'], paper['abstract'])
            summary_html = summary.replace('\n', '<br>')
            
            html_body += f"""
            <div style='background-color: #f9f9f9; padding: 15px; border-radius: 8px; margin-bottom: 20px; box-shadow: 0 2px 4px rgba(0,0,0,0.1);'>
                <p style='font-size: 16px; font-weight: bold; margin-top: 0;'>
                    <a href='{paper['link']}' target='_blank' style='text-decoration: none; color: #1e4f91;'>
                        [{i}] {paper['title']}
                    </a>
                </p>
                <p style='color: #666; font-size: 12px; margin-bottom: 10px;'>📖 Journal: {paper['journal']}</p>
                <div style='background-color: #ffffff; padding: 12px; border: 1px solid #eee; border-radius: 4px; line-height: 1.6;'>
                    {summary_html}
                </div>
            </div>
            """
            total_papers_found += 1
        
        html_body += "<br>"

    if total_papers_found > 0:
        send_email(html_body)
    else:
        print("새로운 논문이 없어 메일을 보내지 않았습니다.")

if __name__ == "__main__":
    main()
