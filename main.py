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

# API 키 설정
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
# [기능 2] Gemini 임상 판단 (GPT 아이디어 적용)
# ==========================================
def summarize_paper(title, abstract):
    if not GEMINI_API_KEY:
        return "N", "오류: 키 없음", ""

    # 프롬프트: '정보'가 아니라 '판단'을 요구함
    prompt = f"""
    당신은 20년 차 정형외과 전문의입니다. 후배 의사에게 이 논문을 브리핑합니다.
    
    [분석 지침]
    1. 이 논문이 기존 임상 관행(Practice)을 바꾸거나 도전하는 내용입니까? (Yes/No)
    2. 외래 진료 중 환자에게 30초 안에 설명한다면 어떻게 말하겠습니까? (구어체 한 문장)
    3. 상세 리뷰를 작성하세요.

    [출력 형식]
    반드시 아래 구분자(@)를 사용하여 3부분으로 나누세요.
    
    (Y 또는 N) @ (환자용 한 문장 설명) @ (상세 전문 리뷰: 배경/결과/비판)

    [논문]
    제목: {title}
    초록: {abstract}
    """

    for attempt in range(2): 
        try:
            model = genai.GenerativeModel('gemini-2.5-flash')
            response = model.generate_content(prompt)
            text = response.text
            
            # 결과 파싱 (Impact / One-Liner / Detail)
            parts = text.split("@")
            if len(parts) >= 3:
                impact = parts[0].strip().upper() # Y or N
                one_liner = parts[1].strip()
                detail = parts[2].strip()
                return impact, one_liner, detail
            else:
                return "N", "요약 형식이 올바르지 않음", text
            
        except Exception as e:
            print(f"⚠️ API 에러 ({e}). 60초 대기 후 재시도...")
            time.sleep(60)
            if attempt == 1:
                return "N", "요약 실패", f"에러: {e}"

# ==========================================
# [기능 3] 이메일 전송 (배지 디자인 적용)
# ==========================================
def send_email(content_html):
    if not MY_EMAIL or not MY_APP_PASSWORD:
        return

    msg = MIMEText(content_html, 'html')
    today = datetime.now().strftime('%Y-%m-%d')
    msg['Subject'] = f"📢 [Dr.AI] {today} 정형외과 논문 브리핑"
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
    html_body = """
    <h2 style='color:#2c3e50;'>👨‍⚕️ 오늘의 정형외과 논문 인사이트</h2>
    <p style='color:gray; font-size:12px;'>* 🚨 표시는 임상 관행에 영향을 줄 수 있는 중요 논문입니다.</p>
    <hr>
    """
    total_papers_found = 0

    for keyword in SEARCH_KEYWORDS:
        papers = fetch_pubmed_papers(keyword, max_results=5)
        html_body += f"<h3 style='color: #2980b9; margin-top:30px;'>🔍 {keyword}</h3>"

        if not papers:
            html_body += "<p style='color: #bdc3c7;'>최근 48시간 내 새로운 논문 없음</p>"
            continue

        for i, paper in enumerate(papers, 1):
            if total_papers_found > 0: 
                time.sleep(30) # 30초 대기 (안전)

            # Gemini가 판단한 결과 받기
            impact_flag, one_liner, deep_rev = summarize_paper(paper['title'], paper['abstract'])
            
            # [배지 로직] Y면 빨간색 주목 배지, N이면 평범한 배지
            if "Y" in impact_flag:
                badge = "<span style='background-color:#e74c3c; color:white; padding:3px 8px; border-radius:4px; font-size:12px; font-weight:bold;'>🚨 주목! 임상 변화 가능</span>"
                box_color = "#fff5f5" # 연한 빨간 배경
            else:
                badge = "<span style='background-color:#95a5a6; color:white; padding:3px 8px; border-radius:4px; font-size:12px; font-weight:bold;'>✅ 기존 견해 유지</span>"
                box_color = "#ffffff" # 흰색 배경

            deep_rev_html = deep_rev.replace('\n', '<br>')
            
            html_body += f"""
            <div style='background-color: {box_color}; border: 1px solid #dfe6e9; border-radius: 8px; margin-bottom: 15px; padding: 15px; box-shadow: 0 1px 3px rgba(0,0,0,0.1);'>
                
                <p style='font-size: 16px; font-weight: bold; margin: 0 0 5px 0;'>
                    {badge} <a href='{paper['link']}' target='_blank' style='text-decoration: none; color: #2c3e50;'>[{i}] {paper['title']}</a>
                </p>
                <p style='font-size: 11px; color: #7f8c8d; margin: 0 0 10px 0;'>📖 {paper['journal']}</p>
                
                <div style='background-color: #ecf0f1; color: #2c3e50; padding: 10px; border-left: 4px solid #3498db; font-size: 14px; margin-bottom: 10px;'>
                    <strong>🗣️ 환자 설명용:</strong> "{one_liner}"
                </div>

                <details>
                    <summary style='cursor: pointer; color: #2980b9; font-size: 13px; font-weight: bold;'>
                        🔽 전문의 상세 리뷰 보기
                    </summary>
                    <div style='background-color: #f8f9fa; padding: 15px; margin-top: 10px; border-radius: 5px; font-size: 14px; line-height: 1.6; color: #444;'>
                        {deep_rev_html}
                    </div>
                </details>
            </div>
            """
            total_papers_found += 1

    if total_papers_found > 0:
        send_email(html_body)
    else:
        print("보낼 논문이 없어 종료합니다.")

if __name__ == "__main__":
    main()
