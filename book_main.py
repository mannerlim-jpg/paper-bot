import os
import json
import smtplib
import ssl
import urllib.request
import urllib.error
from email.mime.text import MIMEText
from email.mime.multipart import MIMEMultipart
from datetime import datetime

def get_book_recommendation():
    """
    urllib을 사용하여 Google Gemini API (REST)를 직접 호출합니다.
    SDK(google.generativeai)를 사용하지 않아 버전 호환성 문제에서 자유롭습니다.
    """
    api_key = os.environ.get("GEMINI_API_KEY")
    if not api_key:
        raise ValueError("GEMINI_API_KEY가 설정되지 않았습니다.")

    # 1. API 엔드포인트 설정 (gemini-1.5-flash 모델 사용)
    url = f"https://generativelanguage.googleapis.com/v1beta/models/gemini-1.5-flash:generateContent?key={api_key}"

    # 2. 정형외과 의사를 위한 프롬프트 구성
    prompt = """
    당신은 20년 차 베테랑 정형외과 전문의를 위한 지적인 비서입니다.
    오늘 그에게 영감을 줄 수 있는 책을 한 권 추천해주세요.
    
    [추천 카테고리 예시]
    1. 정형외과 최신 트렌드나 수술 기법에 관한 전문 서적
    2. 의학의 역사나 의료 윤리를 다룬 인문학 서적
    3. 리더십이나 마음의 평화를 주는 철학 서적 (예: 니체, 스토아 학파 등)
    
    [출력 형식]
    - 책 제목:
    - 저자:
    - 추천 이유 (3문장 이내, 의사의 관점에서):
    """

    # 3. 요청 데이터 구성 (JSON)
    payload = {
        "contents": [{
            "parts": [{"text": prompt}]
        }]
    }
    data = json.dumps(payload).encode("utf-8")

    # 4. HTTP 요청 헤더 설정
    headers = {
        "Content-Type": "application/json"
    }

    # 5. REST API 호출 (urllib 사용)
    try:
        req = urllib.request.Request(url, data=data, headers=headers, method="POST")
        with urllib.request.urlopen(req) as response:
            response_body = response.read().decode("utf-8")
            response_json = json.loads(response_body)
            
            # 응답 파싱
            try:
                text = response_json["candidates"][0]["content"]["parts"][0]["text"]
                return text
            except (KeyError, IndexError) as e:
                return f"API 응답 형식이 예상과 다릅니다: {e}"
                
    except urllib.error.HTTPError as e:
        return f"HTTP 에러 발생: {e.code} - {e.reason}"
    except urllib.error.URLError as e:
        return f"URL 에러 발생: {e.reason}"

def send_email(content):
    """
    SMTP를 사용하여 이메일을 발송합니다. (UTF-8 인코딩 적용)
    """
    sender_email = os.environ.get("MY_EMAIL")
    sender_password = os.environ.get("MY_APP_PASSWORD")
    
    # 수신자도 본인으로 설정 (필요 시 변경 가능)
    receiver_email = sender_email

    if not sender_email or not sender_password:
        raise ValueError("이메일 설정(MY_EMAIL, MY_APP_PASSWORD)이 누락되었습니다.")

    # 이메일 메시지 객체 생성
    msg = MIMEMultipart()
    msg["From"] = sender_email
    msg["To"] = receiver_email
    
    # 제목에 날짜 포함
    today = datetime.now().strftime("%Y-%m-%d")
    msg["Subject"] = f"[{today}] 오늘의 정형외과 의사 추천 도서"

    # 본문 추가 (UTF-8 명시)
    body = f"""
    안녕하세요, 원장님.
    오늘의 추천 도서를 보내드립니다.
    
    --------------------------------------------------
    {content}
    --------------------------------------------------
    
    오늘도 평안한 진료 되시길 바랍니다.
    """
    msg.attach(MIMEText(body, "plain", "utf-8"))

    # SMTP 서버 연결 및 발송 (Gmail 기준)
    try:
        context = ssl.create_default_context()
        with smtplib.SMTP_SSL("smtp.gmail.com", 465, context=context) as server:
            server.login(sender_email, sender_password)
            server.sendmail(sender_email, receiver_email, msg.as_string())
        print("이메일 발송 성공!")
    except Exception as e:
        print(f"이메일 발송 실패: {e}")

if __name__ == "__main__":
    print("Gemini API 호출 중...")
    recommendation = get_book_recommendation()
    
    print("추천 내용:")
    print(recommendation)
    
    print("\n이메일 전송 중...")
    send_email(recommendation)
