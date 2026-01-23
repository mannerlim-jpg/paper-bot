import os
import json
import smtplib
import ssl
import urllib.request
import urllib.error
import random
import time
from email.mime.text import MIMEText
from email.mime.multipart import MIMEMultipart
from datetime import datetime

def get_book_recommendation():
    """
    브라우저로 직접 확인한 'models/gemini-1.5-flash' 경로를 그대로 사용합니다.
    """
    api_key = os.environ.get("GEMINI_API_KEY", "").strip()
    if not api_key:
        return "오류: GEMINI_API_KEY가 설정되지 않았습니다."

    # 1. 테마 선택
    themes = ["철학", "역사", "불교", "빅 히스토리", "실존주의", "인공지능", "에세이"]
    today_theme = random.choice(themes)

    # 2. 검증된 v1beta 엔드포인트 및 모델명
    model_name = "gemini-1.5-flash"
    url = f"https://generativelanguage.googleapis.com/v1beta/models/{model_name}:generateContent?key={api_key}"

    # 3. 구글 API 표준 Payload 구조
    prompt = f"당신은 지적인 독서 큐레이터입니다. '{today_theme}' 주제의 책 1권을 한국어로 추천하고, 추천 이유와 인상 깊은 구절을 적어주세요."
    payload = {
        "contents": [{"parts": [{"text": prompt}]}]
    }
    
    data = json.dumps(payload).encode("utf-8")
    headers = {"Content-Type": "application/json"}

    # 4. 통신 시도 (재시도 로직 포함)
    for attempt in range(3):
        try:
            req = urllib.request.Request(url, data=data, headers=headers, method="POST")
            with urllib.request.urlopen(req) as response:
                res_body = response.read().decode("utf-8")
                res_json = json.loads(res_body) # 77라인 문법 오류 방지를 위해 완전하게 작성
                
                # 결과 추출
                try:
                    text = res_json["candidates"][0]["content"]["parts"][0]["text"]
                    return f"Selected Theme: [{today_theme}]\n\n{text}"
                except (KeyError, IndexError):
                    return f"API 응답 구조 이상: {res_body}"

        except urllib.error.HTTPError as e:
            if e.code == 429: # 사용량 제한 시 대기
                time.sleep(10)
                continue
            return f"HTTP 오류 ({e.code}): {e.reason}"
        except Exception as e:
            return f"알 수 없는 오류: {str(e)}"

    return "❌ 3회 재시도 실패"

def send_email(content):
    sender_email = os.environ.get("MY_EMAIL", "").strip()
    sender_password = os.environ.get("MY_APP_PASSWORD", "").strip()
    
    msg = MIMEMultipart()
    msg["From"] = sender_email
    msg["To"] = sender_email
    msg["Subject"] = f"[{datetime.now().strftime('%Y-%m-%d')}] 오늘의 추천 도서 (Verified)"
    msg.attach(MIMEText(content, "plain", "utf-8"))

    try:
        context = ssl.create_default_context()
        with smtplib.SMTP_SSL("smtp.gmail.com", 465, context=context) as server:
            server.login(sender_email, sender_password)
            server.sendmail(sender_email, sender_email, msg.as_string())
        print("✅ 발송 성공")
    except Exception as e:
        print(f"❌ 발송 실패: {e}")

if __name__ == "__main__":
    result = get_book_recommendation()
    send_email(result)
