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
    [성공 공식 적용]
    1. 모델명: 'gemini-1.5-flash' (latest 삭제)
    2. 엔드포인트: v1beta
    3. 라이브러리: urllib (호환성 문제 해결)
    """
    # 1. GitHub Secrets에서 새 키 가져오기
    api_key = os.environ.get("GEMINI_API_KEY", "").strip()
    if not api_key:
        return "오류: GEMINI_API_KEY가 설정되지 않았습니다."

    themes = [
        "니체의 철학을 현대적으로 해석한 책",
        "제1차 세계대전과 지정학적 변화를 다룬 역사서",
        "불교 철학의 심오함을 쉽게 풀어낸 에세이",
        "빅 히스토리(Big History) 관점의 인문학",
        "실존주의 문학의 숨겨진 명작",
        "인공지능과 인간의 미래를 다룬 통찰력 있는 책",
        "일본의 미학(Wabi-sabi)이나 장인 정신을 다룬 책",
        "의사가 쓴 죽음과 삶에 대한 에세이"
    ]
    today_theme = random.choice(themes)

    # [핵심 1] 군더더기 없는 정확한 모델명 사용
    model_name = "gemini-1.5-flash" 
    
    # [핵심 2] v1beta 엔드포인트 사용
    url = f"https://generativelanguage.googleapis.com/v1beta/models/{model_name}:generateContent?key={api_key}"

    prompt = f"""
    당신은 정형외과 의사를 위한 독서 큐레이터입니다.
    다음 주제에 맞는 책 1권을 한국어로 추천해주세요.
    
    주제: {today_theme}
    
    형식:
    1. 책 제목 / 저자
    2. 추천 이유 (3줄 요약)
    3. 인상 깊은 구절
    """

    payload = {
        "contents": [{"parts": [{"text": prompt}]}]
    }
    data = json.dumps(payload).encode("utf-8")
    headers = {"Content-Type": "application/json"}

    # [핵심 3] 재시도 로직 (429 에러 대비)
    for attempt in range(3):
        try:
            req = urllib.request.Request(url, data=data, headers=headers, method="POST")
            with urllib.request.urlopen(req) as response:
                response_body = response.read().decode("utf-8")
                response_json = json.loads(response_body)
                
                try:
                    text = response_json["candidates"][0]["content"]["parts"][0]["text"]
                    return f"Selected Theme: [{today_theme}]\n\n{text}"
                except (KeyError, IndexError):
                    return f"파싱 오류: {response_body}"

        except urllib.error.HTTPError as e:
            error_content = e.read().decode("utf-8")
            # 429 = 사용량 제한 (잠시 대기 후 재시도)
            if e.code == 429:
                print(f"⚠️ 사용량 제한(429) 감지. 5초 후 재시도... ({attempt+1}/3)")
                time.sleep(5)
                continue
            # 404 = 모델명 오류 (이 코드에선 발생 안 해야 정상)
            return f"HTTP 오류 ({e.code}): {e.reason}\n상세: {error_content}"
        except Exception as e:
            return f"알 수 없는 오류: {str(e)}"

    return "❌ 3회 재시도 실패. (새 프로젝트 키가 맞는지 확인해주세요)"

def send_email(content):
    sender_email = os.environ.get("MY_EMAIL", "").strip()
    sender_password = os.environ.get("MY_APP_PASSWORD", "").strip()
    
    # [수정] 수신자 이메일이 Secrets에 없으면 발신자와 동일하게 설정
    receiver_email = os.environ.get("RECEIVER_EMAIL", "").strip()
    if not receiver_email:
        receiver_email = sender_email

    if not sender_email or not sender_password:
        print("❌ 이메일 설정 누락")
        return

    msg = MIMEMultipart()
    msg["From"] = sender_email
    msg["To"] = receiver_email
    today = datetime.now().strftime("%Y-%m-%d")
    msg["Subject"] = f"[{today}] 오늘의 추천 도서 (Gemini 1.5)"

    body = f"""
    원장님, 성공 사례 분석 결과로 수정한 코드입니다.
    
    ==================================================
    {content}
    ==================================================
    """
    msg.attach(MIMEText(body, "plain", "utf-8"))

    try:
        context = ssl.create_default_context()
        with smtplib.SMTP_SSL("smtp.gmail.com", 465, context=context) as server:
            server.login(sender_email, sender_password)
            server.sendmail(sender_email, receiver_email, msg.as_string())
        print("✅ 이메일 발송 성공")
    except Exception as e:
        print(f"❌ 이메일 발송 실패: {e}")

if __name__ == "__main__":
    print("Gemini 1.5 Flash (검증된 설정) 호출 중...")
    result = get_book_recommendation()
    print(result)
    send_email(result)
