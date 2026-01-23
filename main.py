import smtplib
import os
import time
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

# API 키 공백 제거 및 설정
if GEMINI_API_KEY:
    clean_key = GEMINI_API_KEY.strip()
    genai.configure(api_key=clean_key)
Entrez.email = MY_EMAIL 

# ==========================================
# [기능 1] 논문 검색 (중복 방지 유지)
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
# [기능 2] Gemini 심층 리뷰 (업그레이드 됨!)
# ==========================================
def summarize_paper(title, abstract):
    if not GEMINI_API_KEY:
        return "오류: API 키가 없습니다."

    try:
        model = genai.GenerativeModel('gemini-2.5-flash')
        
        # [핵심] AI에게 주는 명령서(Prompt)를 대폭 강화했습니다.
        prompt = f"""
        당신은 숙련된 정형외과 전문의(Orthopedic Surgeon)입니다.
        아래 논문의 초록(Abstract)을 읽고, 동료 의사에게 브리핑하듯이 전문적으로 리뷰해주세요.
        
        [지침]
        1. 단순 번역이 아닌, '대한정형외과학회'에서 통용되는 전문 의학 용어를 사용하세요.
        2. 문장은 간결하고 명확하게 작성하세요.
        3. 아래 출력 형식을 엄격히 지키세요.

        [논문 정보]
        제목: {title}
        초록: {abstract}
        
        [출력 형식]
        1. **연구 배경 및 목적**: (왜 이 연구를 했는지 1~2문장으로 요약)
        2. **주요 결과 및 수치**: (P-value, 환자 수 등 구체적 수치를 포함하여 핵심 결과 요약)
        3. **비판적 고찰(Critical Review)**: (이 연구의 장점, 또는 한계점이나 주의할 점을 분석)
        4. **임상적 제언(Take-home Message)**: (임상 현장에서 적용할 점 한 줄 요약)
        """
        response = model.generate_content(prompt)
        return response.text
        
    except Exception as e:
        try:
            model = genai.GenerativeModel('gemini-2.0-flash')
            response = model.generate_content(prompt)
            return response.text
        except:
            return f"리뷰 생성 실패 (에러: {e})"

# ==========================================
# [기능 3] 이메일 전송
# ==========================================
def send_email(content_html):
    if not MY_EMAIL or not MY_APP_PASSWORD:
        return

    msg = MIMEText(content_html, 'html')
    today = datetime.now().strftime('%Y-%m-%d')
    msg['Subject'] = f"📢 [Deep Review] {today} 정형외과 최신 논문 분석"
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
    html_body = "<h2>📅 최근 48시간 내 발표된 논문 심층 리뷰</h2><hr>"
    total_papers_found = 0

    for keyword in SEARCH_KEYWORDS:
        papers = fetch_pubmed_papers(keyword, max_results=5)
        html_body += f"<h3 style='color: #2E86C1; border-left: 5px solid #2E86C1; padding-left: 10px;'>🔍 {keyword}</h3>"

        if not papers:
            html_body += "<p style='color: gray; font-style: italic;'>최근 2일간 새로 등록된 논문이 없습니다.</p><br>"
            continue

        for i, paper in enumerate(papers, 1):
            # API 과부하 방지 (15초 대기)
            if i > 1: 
                print("심층 분석을 위해 잠시 생각 정리 중... (15초 대기)")
                time.sleep(15)

            summary = summarize_paper(paper['title'], paper['abstract'])
            summary_html = summary.replace('\n', '<br>')
            
            # [디자인 업그레이드] 가독성을 높였습니다
            html_body += f"""
            <div style='background-color: #f8f9fa; padding: 20px; border-radius: 10px; margin-bottom: 30px; border: 1px solid #e9ecef; box-shadow: 0 4px 6px rgba(0,0,0,0.05);'>
                <p style='font-size: 18px; font-weight: bold; margin-top: 0; color: #1e4f91;'>
                    <a href='{paper['link']}' target='_blank' style='text-decoration: none; color: #1e4f91;'>
                        [{i}] {paper['title']}
                    </a>
                </p>
                <p style='color: #666; font-size: 13px; margin-bottom: 15px; border-bottom: 1px solid #ddd; padding-bottom: 10px;'>
                    📖 Journal: <strong>{paper['journal']}</strong>
                </p>
                <div style='background-color: #ffffff; padding: 15px; border-radius: 5px; line-height: 1.8; font-size: 15px; color: #333;'>
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
